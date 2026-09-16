import csv
import os
import gzip
import struct
import threading

import pytest

from metatracer import annotate as annotate_mod

from metatracer.annotate import (
    HitRecord,
    MappingRow,
    add_eggnog_columns,
    chunk_reader,
    load_mapping_table,
    load_mapping_tables,
    parse_assignments_line,
    parse_gff_attributes,
    parse_hit,
    prepare_annotation_resources,
    read_eggnog_annotations,
    run_eggnog_mapper,
    validate_gff_sort_order,
    write_sorted_chunk,
)


def test_parse_assignments_line_handles_colons_and_hits():
    read_id, hits = parse_assignments_line("read:1:2:123-foo-5=1,999=2")
    assert read_id == "read:1:2"
    assert hits == ["123-foo-5=1", "999=2"]


@pytest.mark.parametrize("protein_sequence", ["", ">protein1\nMAAA\n"])
@pytest.mark.parametrize("custom_output", [False, True])
def test_skip_eggnog_preserves_deposited_annotations(tmp_path, monkeypatch, protein_sequence, custom_output):
    assignments = tmp_path / "reads.clp"
    assignments.write_text("read1:123=0\n")
    output = tmp_path / "annotated.tsv"
    proteins_out = tmp_path / "proteins" / "unique.faa" if custom_output else None
    monkeypatch.setattr(annotate_mod, "load_mapping_tables", lambda paths: ({}, {}))
    monkeypatch.setattr(annotate_mod, "prepare_annotation_resources", lambda *args: ({}, True))
    monkeypatch.setattr(annotate_mod, "build_taxid_name_map", lambda ids: {123: "Taxon"})

    def fake_merge(**kwargs):
        assert not kwargs["taxa_only"]
        with open(kwargs["out_tsv"], "w") as handle:
            handle.write("ReadID\tProtein ID\tAnnotation\nread1\tprotein1\tATPase\n")
        with open(kwargs["proteins_fasta_out"], "w") as handle:
            handle.write(protein_sequence)

    def unexpected_emapper(*args, **kwargs):
        pytest.fail("eggNOG must not run when skip_eggnog is true")

    monkeypatch.setattr(annotate_mod, "merge_sorted_chunks", fake_merge)
    monkeypatch.setattr(annotate_mod, "run_eggnog_mapper", unexpected_emapper)
    assert annotate_mod.run(
        assignments=str(assignments), map_table=[], out=str(output),
        data_dir=str(tmp_path), tmpdir=str(tmp_path / "chunks"), skip_eggnog=True,
        proteins_out=str(proteins_out) if proteins_out else None,
    ) == 0
    with output.open() as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["Annotation"] == "ATPase"
    assert row["Protein ID"] == "protein1"
    assert row["Eggnog"] == "SKIPPED"
    assert all(value == "" for key, value in row.items() if key.startswith("eggnog_"))
    saved = proteins_out or tmp_path / "annotated.tsv.proteins.faa"
    assert saved.read_text() == protein_sequence


def test_protein_output_rejects_annotation_path_collision(tmp_path):
    with pytest.raises(SystemExit, match="must be distinct"):
        annotate_mod.run("input.clp", [], str(tmp_path / "out.tsv"),
                         proteins_out=str(tmp_path / "out.tsv"))


def test_protein_output_rejects_taxa_only(tmp_path):
    with pytest.raises(SystemExit, match="cannot be used with --taxa-only"):
        annotate_mod.run("input.clp", [], str(tmp_path / "out.tsv"),
                         taxa_only=True, proteins_out=str(tmp_path / "proteins.faa"))


def test_parse_assignments_line_blank_and_no_colon():
    assert parse_assignments_line("") == ("", [])
    assert parse_assignments_line("readonly") == ("readonly", [])


def test_parse_hit_full_and_taxa_only():
    assert parse_hit("123-acc-10=2") == (123, "acc", 10, 2)
    assert parse_hit("999=3") == (999, "", 0, 3)


def test_parse_hit_invalid_raises():
    with pytest.raises(ValueError):
        parse_hit("nope")


def test_parse_gff_attributes():
    attrs = parse_gff_attributes("ID=cds1;product=Foo;Note=bar=baz")
    assert attrs["ID"] == "cds1"
    assert attrs["product"] == "Foo"
    assert attrs["Note"] == "bar=baz"


def test_write_sorted_chunk_and_reader_roundtrip(tmp_path):
    records = [
        HitRecord("r2", 2, "asm2", "contig2", "desc2", 10, 1, "k2"),
        HitRecord("r1", 1, "asm1", "contig3", "desc1", 5, 0, "k1"),
        HitRecord("r3", 3, "asm1", "contig1", "desc3", 7, 2, "k3"),
    ]
    out_path = tmp_path / "chunk.tsv"
    write_sorted_chunk(records, str(out_path))

    rows = list(chunk_reader(str(out_path)))
    order = [(r[2], r[3], r[5]) for r in rows]
    assert order == sorted(order)


def test_reference_build_manifest_does_not_need_resource_columns(tmp_path):
    manifest = tmp_path / "reference.map.tsv"
    manifest.write_text(
        "seqid\ttaxid\tassembly\theader\tdescription\n"
        "1\t101\tGCF_000001.1\tNC_1\tcontig one\n",
        encoding="utf-8",
    )
    by_seqid, by_assembly = load_mapping_table(str(manifest))
    assert by_seqid[(101, "1")].assembly == "GCF_000001.1"
    assert by_assembly["GCF_000001.1"].gff_path == ""


def test_multiple_manifests_are_combined_by_taxid_and_seqid(tmp_path):
    first = tmp_path / "first.tsv"
    second = tmp_path / "second.tsv"
    header = "seqid\ttaxid\tassembly\theader\tdescription\n"
    first.write_text(header + "1\t101\tGCF_1.1\tNC_1\tfirst\n", encoding="utf-8")
    second.write_text(header + "1\t202\tGCF_2.1\tNC_2\tsecond\n", encoding="utf-8")

    by_hit, by_assembly = load_mapping_tables([str(first), str(second)])

    assert by_hit[(101, "1")].assembly == "GCF_1.1"
    assert by_hit[(202, "1")].assembly == "GCF_2.1"
    assert set(by_assembly) == {"GCF_1.1", "GCF_2.1"}


def test_multiple_manifests_reject_ambiguous_taxid_seqid_pair(tmp_path):
    first = tmp_path / "first.tsv"
    second = tmp_path / "second.tsv"
    header = "seqid\ttaxid\tassembly\theader\tdescription\n"
    first.write_text(header + "1\t101\tGCF_1.1\tNC_1\tfirst\n", encoding="utf-8")
    second.write_text(header + "1\t101\tGCF_2.1\tNC_2\tsecond\n", encoding="utf-8")

    with pytest.raises(SystemExit, match="Ambiguous mapping across manifests"):
        load_mapping_tables([str(first), str(second)])


def test_prepare_resources_sorts_indexes_and_reports(tmp_path):
    assembly = "GCF_000001.1"
    assembly_dir = tmp_path / "ncbi_dataset" / "data" / assembly
    assembly_dir.mkdir(parents=True)
    (assembly_dir / "genomic.gff").write_text(
        "##gff-version 3\n"
        "NC_1\tRefSeq\tCDS\t20\t30\t.\t+\t0\tID=cds-b\n"
        "NC_1\tRefSeq\tCDS\t5\t10\t.\t+\t0\tID=cds-a\n",
        encoding="utf-8",
    )
    (assembly_dir / "protein.faa").write_text(">cds-a\nMKK\n", encoding="utf-8")
    mapping = MappingRow("1", 101, assembly, "NC_1", "contig", "", "")
    report = tmp_path / "resources.tsv"

    prepared, ready = prepare_annotation_resources(
        {assembly: mapping}, str(tmp_path), str(report)
    )

    assert ready
    assert prepared[assembly].gff_path.endswith(".gz")
    assert os.path.exists(prepared[assembly].gff_path + ".tbi")
    validate_gff_sort_order(prepared[assembly].gff_path)
    rows = list(csv.DictReader(report.open(), delimiter="\t"))
    assert rows[0]["status"] == "READY"
    assert rows[0]["gff_sort_status"] == "SORTED"
    assert rows[0]["gff_index_status"] == "INDEXED"


@pytest.mark.parametrize("index_case", ["valid", "internal_separator", "incompatible"])
def test_gff_separator_index_validation_and_repair(tmp_path, index_case):
    import pysam

    assembly = "GCF_1.1"
    folder = tmp_path / "ncbi_dataset" / "data" / assembly
    folder.mkdir(parents=True)
    original = folder / "genomic.gff"
    contents = "##gff-version 3\n###\n"
    contents += "NC_1\tRefSeq\tCDS\t5\t10\t.\t+\t0\tID=cds-a\n"
    if index_case == "internal_separator":
        contents += "###\n"
    contents += "NC_1\tRefSeq\tCDS\t20\t30\t.\t+\t0\tID=cds-b\n"
    original.write_text(contents)
    (folder / "protein.faa").write_text(">a\nMKK\n")
    compressed = pysam.tabix_index(str(original), preset="gff", keep_original=True)
    index = compressed + ".tbi"
    if index_case == "incompatible":
        # Simulate an index with an incompatible start-coordinate column.
        with gzip.open(index, "rb") as handle:
            data = bytearray(handle.read())
        struct.pack_into("<i", data, 16, 5)
        with pysam.BGZFile(index, "wb") as handle:
            handle.write(bytes(data))
    index_before = open(index, "rb").read()
    mapping = MappingRow("1", 101, assembly, "NC_1", "contig", "", "")
    report = tmp_path / "resources.tsv"
    prepared, ready = prepare_annotation_resources({assembly: mapping}, str(tmp_path), str(report))
    assert ready
    with report.open() as handle:
        row = next(csv.DictReader(handle, delimiter="\t"))
    assert row["gff_index_status"] == ("ALREADY_INDEXED" if index_case == "valid" else "REINDEXED")
    assert original.read_text() == contents
    assert open(index, "rb").read() == index_before
    with pysam.TabixFile(prepared[assembly].gff_path) as tb:
        assert len(list(tb.fetch("NC_1"))) == 2
    # Reuse a derivative on subsequent runs rather than rebuilding it.
    prepared, ready = prepare_annotation_resources({assembly: mapping}, str(tmp_path), str(report))
    assert ready
    with report.open() as handle:
        assert next(csv.DictReader(handle, delimiter="\t"))["gff_index_status"] == "ALREADY_INDEXED"


def test_gff_query_failure_is_not_silently_ignored():
    class BadTabix:
        def fetch(self, contig):
            raise ValueError("Failed to parse GFF")

    annotator = annotate_mod.IntervalGFFAnnotator.__new__(annotate_mod.IntervalGFFAnnotator)
    annotator._get_tabix = lambda assembly: BadTabix()
    annotator._IntervalTree = list
    annotator.by_assembly = {"asm": MappingRow("1", 101, "asm", "NC_1", "", "bad.gff.gz", "")}
    with pytest.raises(annotate_mod.GFFQueryError, match="partial CDS annotations discarded"):
        annotator._build_tree_for_contig("asm", "NC_1")


@pytest.mark.parametrize("threads", [1, 2])
def test_parallel_preparation_reports_all_assemblies_in_order(tmp_path, monkeypatch, threads):
    mappings = {}
    for assembly in ["GCF_Z.1", "GCF_A.1", "GCF_MISSING.1"]:
        folder = tmp_path / "ncbi_dataset" / "data" / assembly
        folder.mkdir(parents=True)
        mappings[assembly] = MappingRow("1", 101, assembly, "NC_1", "", "", "")
        if "MISSING" not in assembly:
            (folder / "genomic.gff").write_text(
                "NC_1\tRefSeq\tCDS\t5\t10\t.\t+\t0\tID=cds-a\n"
                "###\n"
                "NC_1\tRefSeq\tCDS\t20\t30\t.\t+\t0\tID=cds-b\n"
            )
            (folder / "protein.faa").write_text(">a\nMKK\n")
    validator = annotate_mod.validate_gff_index
    barrier = threading.Barrier(2)
    worker_names = []

    def validate(path, pysam, source_path=None):
        worker_names.append(threading.current_thread().name)
        if threads > 1:
            # Both independent assemblies must reach validation concurrently.
            barrier.wait(timeout=10)
        return validator(path, pysam, source_path)

    monkeypatch.setattr(annotate_mod, "validate_gff_index", validate)
    report = tmp_path / "resources.tsv"
    prepared, ready = prepare_annotation_resources(mappings, str(tmp_path), str(report), threads=threads)
    assert not ready
    assert list(prepared) == ["GCF_A.1", "GCF_Z.1"]
    with report.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["accession"] for row in rows] == sorted(mappings)
    assert [row["status"] for row in rows] == ["READY", "MISSING_OR_AMBIGUOUS_RESOURCE", "READY"]
    assert len(set(worker_names)) == threads


def test_preparation_rejects_invalid_thread_count(tmp_path):
    with pytest.raises(ValueError, match="threads must be at least 1"):
        prepare_annotation_resources({}, str(tmp_path), str(tmp_path / "report.tsv"), threads=0)


def test_prepare_resources_reports_missing_files(tmp_path):
    assembly = "GCF_000002.1"
    (tmp_path / "ncbi_dataset" / "data" / assembly).mkdir(parents=True)
    mapping = MappingRow("1", 102, assembly, "NC_2", "contig", "", "")
    report = tmp_path / "resources.tsv"

    prepared, ready = prepare_annotation_resources(
        {assembly: mapping}, str(tmp_path), str(report)
    )

    assert not ready
    assert prepared == {}
    rows = list(csv.DictReader(report.open(), delimiter="\t"))
    assert rows[0]["status"] == "MISSING_OR_AMBIGUOUS_RESOURCE"
    assert "GFF: no files matched" in rows[0]["message"]
    assert "protein FASTA: no files matched" in rows[0]["message"]


def test_prepare_resources_accepts_custom_patterns(tmp_path):
    assembly = "GCF_000003.1"
    resource_dir = tmp_path / "custom" / assembly
    resource_dir.mkdir(parents=True)
    (resource_dir / "features.gff").write_text(
        "##gff-version 3\nNC_3\tRefSeq\tCDS\t1\t3\t.\t+\t0\tID=x\n",
        encoding="utf-8",
    )
    (resource_dir / "proteins.faa").write_text(">x\nM\n", encoding="utf-8")
    mapping = MappingRow("1", 103, assembly, "NC_3", "contig", "", "", assembly)

    prepared, ready = prepare_annotation_resources(
        {assembly: mapping},
        str(tmp_path),
        str(tmp_path / "resources.tsv"),
        "{basepath}/custom/{accession}/*.gff",
        "{basepath}/custom/{accession}/*.faa",
    )

    assert ready
    assert prepared[assembly].protein_fa_path.endswith("proteins.faa")


def test_eggnog_annotations_join_through_protein_id(tmp_path):
    eggnog = tmp_path / "results.emapper.annotations"
    eggnog.write_text(
        "## emapper version\n"
        "#query\tseed_ortholog\tevalue\teggNOG_OGs\tDescription\n"
        "1\tCOG0001\t1e-20\tCOG1234@1|root,COG1234@2|Bacteria\tATPase\n",
        encoding="utf-8",
    )
    columns, annotations = read_eggnog_annotations(str(eggnog))
    deposited = tmp_path / "deposited.tsv"
    deposited.write_text(
        "ReadID\tProtein ID\tAnnotation\n"
        "read1\t1\tdeposited product\n"
        "read2\tNA\tNA\n",
        encoding="utf-8",
    )
    output = tmp_path / "annotated.tsv"

    add_eggnog_columns(str(deposited), str(output), columns, annotations, "SUCCESS")

    rows = list(csv.DictReader(output.open(), delimiter="\t"))
    assert rows[0]["eggnog_seed_ortholog"] == "COG0001"
    assert rows[0]["eggnog_Description"] == "ATPase"
    assert rows[0]["eggnog_OG"] == "COG1234"
    assert rows[0]["Eggnog"] == "SUCCESS"
    assert rows[1]["eggnog_seed_ortholog"] == ""
    assert rows[1]["eggnog_OG"] == ""


def test_run_eggnog_mapper_uses_unique_protein_fasta(tmp_path, monkeypatch):
    proteins = tmp_path / "unique.faa"
    proteins.write_text(">1\nMKK\n", encoding="utf-8")
    calls = []

    def fake_run(command, check):
        calls.append(command)
        (tmp_path / "metatracer.emapper.annotations").write_text(
            "#query\tDescription\n1\tATPase\n", encoding="utf-8"
        )

    monkeypatch.setattr("metatracer.annotate.shutil.which", lambda _name: "/bin/emapper.py")
    monkeypatch.setattr("metatracer.annotate.subprocess.run", fake_run)
    result = run_eggnog_mapper(
        str(proteins), str(tmp_path), "emapper.py", 4, "/db", ("--go_evidence", "all")
    )

    assert result.endswith("metatracer.emapper.annotations")
    assert calls[0][calls[0].index("-i") + 1] == str(proteins)
    assert calls[0][calls[0].index("--cpu") + 1] == "4"
    assert calls[0][calls[0].index("--data_dir") + 1] == "/db"


def test_failed_eggnog_status_preserves_deposited_annotations(tmp_path):
    deposited = tmp_path / "deposited.tsv"
    deposited.write_text(
        "ReadID\tProtein ID\tAnnotation\nread1\t1\tdeposited product\n",
        encoding="utf-8",
    )
    output = tmp_path / "annotated.tsv"
    columns = [("Description", "eggnog_Description")]

    add_eggnog_columns(str(deposited), str(output), columns, {}, "FAILED")

    rows = list(csv.DictReader(output.open(), delimiter="\t"))
    assert rows[0]["Annotation"] == "deposited product"
    assert rows[0]["Eggnog"] == "FAILED"
    assert rows[0]["eggnog_Description"] == ""
    assert rows[0]["eggnog_OG"] == ""
