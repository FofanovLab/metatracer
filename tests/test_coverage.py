import csv

import pytest
from click.testing import CliRunner

from metatracer.cli import cli


def test_coverage_unions_per_genome_then_takes_max(tmp_path):
    source = tmp_path / "input.tsv"
    source.write_text(
        "read_id\tread_length\ttaxa\tGID\tposition\tedit_distance\n"
        "r1\t100\t2;2;2\t10;10;11\t0;0;0\t0;0;0\n"
        "r2\t100\t2;2\t10;11\t50;200\t0;0\n"
        "r3\t20\t3\t12\t0\t0\n"
    )
    out, genomes = tmp_path / "out.tsv", tmp_path / "genomes.tsv"
    result = CliRunner().invoke(cli, ["coverage", "--input", str(source), "--out", str(out),
                                     "--genome-out", str(genomes), "--approx-read-length", "999"])
    assert result.exit_code == 0, result.output
    rows = list(csv.DictReader(out.open(), delimiter="\t"))
    assert {key: rows[0][key] for key in ("taxid", "supporting_reads", "genomes_hit",
                "max_covered_bases", "max_coverage_gid", "max_coverage_supporting_reads")} == dict(taxid="2", supporting_reads="2", genomes_hit="2",
                          max_covered_bases="200", max_coverage_gid="11",
                          max_coverage_supporting_reads="2")
    assert rows[1]["max_covered_bases"] == "20"
    per_genome = list(csv.DictReader(genomes.open(), delimiter="\t"))
    assert [r["covered_bases"] for r in per_genome] == ["150", "200", "20"]


@pytest.mark.parametrize("length_column,row_length", [("read_length\t", "0\t"), ("", "")])
def test_missing_length_requires_fallback(tmp_path, length_column, row_length):
    source, out = tmp_path / "in.tsv", tmp_path / "out.tsv"
    source.write_text(f"read_id\t{length_column}taxa\tGID\tposition\n"
                      f"r1\t{row_length}2\t10\t5\n")
    args = ["coverage", "--input", str(source), "--out", str(out)]
    result = CliRunner().invoke(cli, args)
    assert result.exit_code != 0
    assert "--approx-read-length" in result.output
    assert not out.exists()
    assert CliRunner().invoke(cli, args + ["--approx-read-length", "75"]).exit_code == 0
    assert list(csv.DictReader(out.open(), delimiter="\t"))[0]["max_covered_bases"] == "75"


def test_mismatched_lists_fail_without_report(tmp_path):
    source, out = tmp_path / "in.tsv", tmp_path / "out.tsv"
    source.write_text("read_id\tread_length\ttaxa\tGID\tposition\n"
                      "r1\t100\t2;3\t10\t0;50\n")
    result = CliRunner().invoke(cli, ["coverage", "--input", str(source), "--out", str(out)])
    assert result.exit_code != 0
    assert "Line 2" in result.output
    assert not out.exists()


def test_edit_distributions_use_best_hit_per_read_and_normalize_length(tmp_path):
    import json

    source, out, genomes = [tmp_path / name for name in ("in.tsv", "out.tsv", "genomes.tsv")]
    source.write_text("read_id\tread_length\ttaxa\tGID\tposition\tedit_distance\n"
                      "r1\t100\t2;2;2\t10;10;11\t0;5;0\t8;4;2\n"
                      "r2\t200\t2;2\t10;11\t100;100\t0;6\n")
    result = CliRunner().invoke(cli, ["coverage", "--input", str(source), "--out", str(out),
                                     "--genome-out", str(genomes)])
    assert result.exit_code == 0, result.output
    row = list(csv.DictReader(out.open(), delimiter="\t"))[0]
    assert row["edit_reads"] == "2"
    assert float(row["edit_mean"]) == 1
    assert float(row["edit_median"]) == 1
    assert float(row["edit_p90"]) == pytest.approx(1.8)
    assert float(row["edit_zero_fraction"]) == .5
    assert float(row["edit_rate_mean"]) == .01
    assert json.loads(row["edit_histogram"]) == {"0": 1, "2": 1}
    per_genome = list(csv.DictReader(genomes.open(), delimiter="\t"))
    assert json.loads(per_genome[0]["edit_histogram"]) == {"0": 1, "4": 1}
    assert float(per_genome[1]["edit_rate_mean"]) == .025


@pytest.mark.parametrize("edits", ["-1", "abc", "1;2", ""])
def test_bad_edit_distances_fail_before_writing(tmp_path, edits):
    source, out = tmp_path / "in.tsv", tmp_path / "out.tsv"
    source.write_text("read_id\tread_length\ttaxa\tGID\tposition\tedit_distance\n"
                      f"r1\t100\t2\t10\t0\t{edits}\n")
    result = CliRunner().invoke(cli, ["coverage", "--input", str(source), "--out", str(out)])
    assert result.exit_code != 0
    assert not out.exists()
