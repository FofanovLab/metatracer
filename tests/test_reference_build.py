import json
import csv
from pathlib import Path

from metatracer import reference_build as rb


def test_sniff_tsv_delim(tmp_path):
    tsv = tmp_path / "report.tsv"
    tsv.write_text("a\tb\n1\t2\n", encoding="utf-8")
    assert rb._sniff_tsv_delim(tsv) == "\t"

    csv_path = tmp_path / "report.csv"
    csv_path.write_text("a,b\n1,2\n", encoding="utf-8")
    assert rb._sniff_tsv_delim(csv_path) == ","


def test_extract_from_json_obj_variants():
    asm, tax = rb._extract_from_json_obj({"assembly_accession": "GCF_1", "tax_id": "123"})
    assert asm == "GCF_1"
    assert tax == 123

    obj = {"assembly": {"assembly_accession": "GCF_2", "organism": {"taxid": "456"}}}
    asm, tax = rb._extract_from_json_obj(obj)
    assert asm == "GCF_2"
    assert tax == 456


def test_read_assembly_taxid_report_tsv(tmp_path):
    report = tmp_path / "report.tsv"
    report.write_text(
        "assembly_accession\ttax_id\nGCF_1\t123\nGCF_2\t456\n", encoding="utf-8"
    )
    mapping = rb.read_assembly_taxid_report(report)
    assert mapping == {"GCF_1": 123, "GCF_2": 456}


def test_read_assembly_taxid_report_jsonl(tmp_path):
    report = tmp_path / "report.jsonl"
    report.write_text(
        json.dumps({"assembly_accession": "GCF_1", "tax_id": "123"}) + "\n"
        + json.dumps({"assembly": {"assembly_accession": "GCF_2", "taxid": "456"}})
        + "\n",
        encoding="utf-8",
    )
    mapping = rb.read_assembly_taxid_report(report)
    assert mapping == {"GCF_1": 123, "GCF_2": 456}


def test_locate_assembly_files(tmp_path):
    asm_dir = tmp_path / "GCF_1"
    asm_dir.mkdir()
    fna = asm_dir / "foo_genomic.fna"
    gff = asm_dir / "foo_genomic.gff"
    faa = asm_dir / "foo_protein.faa"
    fna.write_text(">a\nAC\n", encoding="utf-8")
    gff.write_text("##gff-version 3\n", encoding="utf-8")
    faa.write_text(">p\nM\n", encoding="utf-8")

    genomic_fna, gff_path, protein = rb.locate_assembly_files(asm_dir)
    assert genomic_fna == str(fna)
    assert gff_path == str(gff)
    assert protein == str(faa)


def test_fasta_record_bytes_and_write(tmp_path):
    out = tmp_path / "out.fasta"
    header = ">rec1"
    seq = "ACGT"
    expected = len(header) + 1 + len(seq) + 1
    assert rb.fasta_record_bytes(header, seq, wrap=60) == expected

    with out.open("wt", encoding="utf-8", newline="\n") as handle:
        rb.write_fasta_record(handle, header, seq, wrap=60)
    assert out.read_text(encoding="utf-8") == ">rec1\nACGT\n"


def test_build_reference_plans_whole_fastas_and_writes_dual_taxonomy_manifest(
    tmp_path, monkeypatch
):
    data_dir = tmp_path / "data"
    first_dir = data_dir / "GCF_1.1"
    second_dir = data_dir / "GCF_2.1"
    first_dir.mkdir(parents=True)
    second_dir.mkdir(parents=True)
    first_fasta = first_dir / "first_genomic.fna"
    second_fasta = second_dir / "second_genomic.fna"
    first_fasta.write_text(">NC_1 first contig\n" + "A" * 600_000 + "\n")
    second_fasta.write_text(">NC_2 second contig\n" + "C" * 600_000 + "\n")
    report = tmp_path / "manifest.tsv"
    report.write_text(
        "ncbi_accession\tncbi_taxid\tgtdb_representative_code\n"
        "GCF_1.1\t101\tGTDB_100101\n"
        "GCF_2.1\t202\tGTDB_100202\n"
    )
    monkeypatch.setattr(
        rb,
        "resolve_ncbi_species_taxids",
        lambda taxids: ({taxid: taxid for taxid in taxids}, {taxid: "species" for taxid in taxids}),
    )

    lists_dir = tmp_path / "fasta-lists"
    sequence_manifest = tmp_path / "reference.tsv"
    rb.build_reference(
        data_dir=data_dir,
        report_path=report,
        out_dir=lists_dir,
        max_size_mb=1,
        map_tsv_path=sequence_manifest,
        summary_path=tmp_path / "summary.txt",
        taxonomy_map_path=tmp_path / "taxonomy.tsv",
        index_gff=False,
        force_reindex=False,
        mapping_only=False,
        taxonomy_source="ncbi",
    )

    path_lists = sorted(lists_dir.glob("*.fasta-list.txt"))
    assert len(path_lists) == 2
    assert path_lists[0].read_text().strip() == str(first_fasta.resolve())
    assert path_lists[1].read_text().strip() == str(second_fasta.resolve())
    with sequence_manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [(row["header"], row["seqid"], row["taxid"], row["alternate_taxid"], row["index"])
            for row in rows] == [
        ("NC_1", "1", "101", "100101", "0"),
        ("NC_2", "2", "202", "100202", "1"),
    ]
    assert [row["original_alternate_taxid"] for row in rows] == [
        "GTDB_100101", "GTDB_100202"
    ]
    assert "GFF_Path" not in rows[0]
    assert "Protein_Path" not in rows[0]
    assert "gff" not in rows[0]
    assert "protein_fasta" not in rows[0]
    assert not list(lists_dir.glob("*.fasta"))


def test_report_index_column_overrides_size_planning(tmp_path, monkeypatch):
    data_dir = tmp_path / "data"
    for accession, header in (("GCF_1.1", "NC_1"), ("GCF_2.1", "NC_2")):
        assembly_dir = data_dir / accession
        assembly_dir.mkdir(parents=True)
        (assembly_dir / f"{accession}_genomic.fna").write_text(
            f">{header}\nACGT\n", encoding="utf-8"
        )
    report = tmp_path / "manifest.tsv"
    report.write_text(
        "ncbi_accession\tncbi_taxid\tindex\n"
        "GCF_1.1\t101\t7\n"
        "GCF_2.1\t202\t2\n",
        encoding="utf-8",
    )
    monkeypatch.setattr(
        rb,
        "resolve_ncbi_species_taxids",
        lambda taxids: ({taxid: taxid for taxid in taxids}, {taxid: "species" for taxid in taxids}),
    )

    lists_dir = tmp_path / "lists"
    manifest = tmp_path / "reference.tsv"
    rb.build_reference(
        data_dir=data_dir,
        report_path=report,
        out_dir=lists_dir,
        max_size_mb=0,
        map_tsv_path=manifest,
        summary_path=tmp_path / "summary.txt",
        taxonomy_map_path=tmp_path / "taxonomy.tsv",
        index_gff=False,
        force_reindex=False,
        mapping_only=False,
        taxonomy_source="ncbi",
    )

    assert sorted(path.name for path in lists_dir.iterdir()) == [
        "metatracer_reference.index.2.fasta-list.txt",
        "metatracer_reference.index.7.fasta-list.txt",
    ]
    with manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert {row["accession"]: row["index"] for row in rows} == {
        "GCF_1.1": "7",
        "GCF_2.1": "2",
    }
    assert "Index max size (MB):  ignored" in (tmp_path / "summary.txt").read_text()


def test_report_index_column_requires_values_for_every_assembly(tmp_path):
    report = tmp_path / "manifest.tsv"
    report.write_text(
        "assembly_accession\ttax_id\tindex\nGCF_1.1\t101\t\n",
        encoding="utf-8",
    )
    try:
        rb.read_predefined_indices(report)
    except SystemExit as exc:
        assert "Missing index for assembly GCF_1.1" in str(exc)
    else:
        raise AssertionError("missing predefined index was accepted")


def test_taxonomy_label_conversion_preserves_groups_and_warns(caplog):
    normalized, originals = rb.normalize_taxonomy_values(
        {"GCF_1": "taxon-100", "GCF_2": "taxon-100", "GCF_3": "taxon-200"},
        "alternate_taxid",
    )
    assert normalized == {"GCF_1": 100, "GCF_2": 100, "GCF_3": 200}
    assert originals["GCF_1"] == "taxon-100"
    assert "Converted 3 non-integer alternate_taxid value(s)" in caplog.text


def test_taxonomy_label_conversion_rejects_changed_groupings():
    try:
        rb.normalize_taxonomy_values(
            {"GCF_1": "taxon-100", "GCF_2": "species-100"},
            "alternate_taxid",
        )
    except SystemExit as exc:
        assert "does not preserve groupings" in str(exc)
    else:
        raise AssertionError("colliding taxonomy labels were accepted")


def test_taxonomy_id_must_fit_unsigned_32_bits():
    try:
        rb.normalize_taxonomy_values({"GCF_1": str(2**32)}, "taxid")
    except SystemExit as exc:
        assert "unsigned 32-bit integer range" in str(exc)
    else:
        raise AssertionError("out-of-range taxonomy ID was accepted")


def test_supplied_accession_table_uses_ids_directly(tmp_path):
    table = tmp_path / "accessions.tsv"
    table.write_text(
        "accession\ttaxid\talternate_taxid\n"
        "GCF_1.1\t101\t100101\n"
        "GCF_2.1\t202\t\n",
        encoding="utf-8",
    )
    primary, alternate, original_primary, original_alternate = (
        rb.read_supplied_taxonomy_table(table)
    )
    assert primary == {"GCF_1.1": 101, "GCF_2.1": 202}
    assert alternate == {"GCF_1.1": 100101, "GCF_2.1": 202}
    assert original_primary == {"GCF_1.1": "101", "GCF_2.1": "202"}
    assert original_alternate == {"GCF_1.1": "100101", "GCF_2.1": "202"}
