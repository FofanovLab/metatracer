import csv

import pytest

from metatracer.annotate import (
    HitRecord,
    chunk_reader,
    parse_assignments_line,
    parse_gff_attributes,
    parse_hit,
    write_sorted_chunk,
)


def test_parse_assignments_line_handles_colons_and_hits():
    read_id, hits = parse_assignments_line("read:1:2:123-foo-5=1,999=2")
    assert read_id == "read:1:2"
    assert hits == ["123-foo-5=1", "999=2"]


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
