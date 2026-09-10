import csv

import pytest
from click.testing import CliRunner

from metatracer.cli import cli
from metatracer.count import run


def read_counts(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_taxid_multihits_are_counted_as_a_separate_group(tmp_path):
    annotations = tmp_path / "sample.tsv"
    annotations.write_text(
        "sample_id\tReadID\tTaxid\n"
        "soil\tr1\t123\n"
        "soil\tr1\t1234\n"
        "soil\tr1\t1234\n"
        "soil\tr2\t123\n"
        "soil\tr3\t1234;123\n"
    )
    output = tmp_path / "counts.tsv"

    result = run([annotations], output, ["taxid"])

    assert read_counts(output) == [
        {"sample_id": "soil", "taxid": "123", "count": "1"},
        {"sample_id": "soil", "taxid": "123;1234", "count": "2"},
    ]
    assert result["count_rows"] == 2


def test_multiple_columns_count_the_joint_group(tmp_path):
    annotations = tmp_path / "sample.tsv"
    annotations.write_text(
        "ReadID\tTaxid\teggnog_OG\n"
        "r1\t123\tCOG1\n"
        "r1\t1234\tCOG2,COG1\n"
        "r2\t123\tCOG1\n"
    )
    output = tmp_path / "counts.tsv"

    run([annotations], output, ["taxid", "eggnog_OG"])

    assert read_counts(output) == [
        {
            "sample_id": "sample", "taxid": "123", "eggnog_OG": "COG1", "count": "1",
        },
        {
            "sample_id": "sample", "taxid": "123;1234",
            "eggnog_OG": "COG1;COG2", "count": "1",
        },
    ]


def test_missing_values_skip_read_and_missing_columns_fail(tmp_path):
    annotations = tmp_path / "sample.tsv"
    annotations.write_text("read_id\ttaxid\teggnog_OG\nr1\t1\t\nr2\t2\tCOG1\n")
    output = tmp_path / "counts.tsv"

    result = run([annotations], output, ["taxid", "eggnog_og"])
    assert result["skipped_reads"] == 1
    assert read_counts(output)[0]["taxid"] == "2"
    with pytest.raises(ValueError, match="missing required count columns: absent"):
        run([annotations], output, ["absent"])


def test_count_cli_accepts_repeated_columns(tmp_path):
    annotations = tmp_path / "sample.tsv"
    annotations.write_text("ReadID\tTaxid\teggnog_OG\nr1\t1\tCOG1\n")
    output = tmp_path / "counts.tsv"

    result = CliRunner().invoke(
        cli,
        ["count", "--input", str(annotations), "--output", str(output),
         "--column", "taxid", "--column", "eggnog_OG"],
    )

    assert result.exit_code == 0, result.output
    assert read_counts(output)[0]["count"] == "1"
