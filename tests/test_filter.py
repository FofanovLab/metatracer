from pathlib import Path

import pytest

from metatracer.filter import (
    edit_delta_filter,
    load_taxa_file,
    main,
    parse_hit,
    parse_line,
    taxa_filter,
)


def test_parse_hit_and_line():
    hit = parse_hit("123-acc-10=2")
    assert hit.taxid == 123
    assert hit.edit == 2
    assert hit.raw == "123-acc-10=2"

    read_id, hits = parse_line("read:1:2:123-acc-10=2,999=3")
    assert read_id == "read:1:2"
    assert hits == ["123-acc-10=2", "999=3"]


def test_load_taxa_file(tmp_path):
    taxa_path = tmp_path / "taxa.txt"
    taxa_path.write_text(
        "# comment\n123\n456 extra\nbad\n\n", encoding="utf-8"
    )
    taxa = load_taxa_file(str(taxa_path))
    assert taxa == {123, 456}


def test_taxa_and_edit_delta_filters():
    hits = [parse_hit("1=0"), parse_hit("2=3"), parse_hit("3=1")]
    filtered = taxa_filter(hits, include={1, 3}, exclude={3})
    assert [h.taxid for h in filtered] == [1]

    filtered = edit_delta_filter(hits, edit_delta=1)
    assert [h.edit for h in filtered] == [0, 1]


def test_main_filters_file(tmp_path):
    input_path = tmp_path / "input.clp"
    output_path = tmp_path / "out.clp"
    include_path = tmp_path / "include.txt"
    exclude_path = tmp_path / "exclude.txt"

    input_path.write_text(
        "\n".join(
            [
                "read1:1-a-10=0,2-b-20=3,3=1",
                "read2:4=2,4=3",
                "read3:5=1",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    include_path.write_text("1\n3\n", encoding="utf-8")
    exclude_path.write_text("3\n", encoding="utf-8")

    ret = main(
        [
            "--input",
            str(input_path),
            "--out",
            str(output_path),
            "--include-taxa",
            str(include_path),
            "--exclude-taxa",
            str(exclude_path),
            "--edit-delta",
            "0",
        ]
    )
    assert ret == 0

    out_lines = output_path.read_text(encoding="utf-8").strip().splitlines()
    assert out_lines == ["read1:1-a-10=0"]
