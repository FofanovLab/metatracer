import json
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
