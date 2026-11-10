import csv

import pytest
from click.testing import CliRunner

from metatracer.cli import cli
from metatracer.taxa_report_filter import EXPECTED_COLUMNS


def setup_reports(tmp_path):
    source = tmp_path / "taxa.tsv"
    with source.open("w") as stream:
        writer = csv.DictWriter(stream, fieldnames=EXPECTED_COLUMNS, delimiter="\t")
        writer.writeheader()
        for taxid, reads in [(1, 20), (2, 20), (3, 20), (4, 2)]:
            row = {key: 0 for key in EXPECTED_COLUMNS}
            row.update(taxid=taxid, total_reads=reads)
            writer.writerow(row)
    coverage = tmp_path / "coverage.tsv"
    coverage.write_text("taxid\tmax_covered_bases\n1\t99\n2\t100\n4\t200\n")
    args = ["taxa-report-filter", "--input", str(source),
            "--out", str(tmp_path / "out.tsv"),
            "--include-out", str(tmp_path / "include.txt"),
            "--exclude-out", str(tmp_path / "exclude.txt"),
            "--log", str(tmp_path / "log.txt")]
    return args, coverage


def test_coverage_cutoff_combines_with_support_and_missing_fails(tmp_path):
    args, coverage = setup_reports(tmp_path)
    result = CliRunner().invoke(cli, args + ["--coverage-report", str(coverage),
                                            "--min-covered-bases", "100",
                                            "--min-total-reads", "10"])
    assert result.exit_code == 0, result.output
    assert (tmp_path / "include.txt").read_text() == "2\n"
    assert (tmp_path / "exclude.txt").read_text() == "1\n3\n4\n"
    rows = list(csv.DictReader((tmp_path / "out.tsv").open(), delimiter="\t"))
    assert [row["taxid"] for row in rows] == ["2"]
    log = (tmp_path / "log.txt").read_text()
    assert "rows_fail_coverage: 2" in log
    assert "rows_missing_coverage: 1" in log


def test_without_coverage_preserves_behavior(tmp_path):
    args, _ = setup_reports(tmp_path)
    assert CliRunner().invoke(cli, args).exit_code == 0
    assert (tmp_path / "include.txt").read_text() == "1\n2\n3\n4\n"


def test_cutoff_requires_report(tmp_path):
    args, _ = setup_reports(tmp_path)
    result = CliRunner().invoke(cli, args + ["--min-covered-bases", "100"])
    assert result.exit_code != 0
    assert "requires --coverage-report" in result.output
    assert not (tmp_path / "out.tsv").exists()


@pytest.mark.parametrize("content", ["taxid\twrong\n1\t100\n",
    "taxid\tmax_covered_bases\n1\t100\n1\t200\n",
    "taxid\tmax_covered_bases\n1\t-1\n",
    "taxid\tmax_covered_bases\n1\tnan\n"])
def test_invalid_coverage_rejected_before_writing(tmp_path, content):
    args, coverage = setup_reports(tmp_path)
    coverage.write_text(content)
    result = CliRunner().invoke(cli, args + ["--coverage-report", str(coverage),
                                            "--min-covered-bases", "100"])
    assert result.exit_code != 0
    assert not (tmp_path / "out.tsv").exists()
