from unittest.mock import patch

from click.testing import CliRunner

from metatracer.cli import cli


def test_index_build_forwards_repeated_fastas():
    with (
        patch("metatracer.cli._which_or_die", return_value="/bin/mtsv-build"),
        patch("metatracer.cli._run") as run,
    ):
        result = CliRunner().invoke(cli, [
            "index-build", "--fasta", "a.fasta", "--fasta", "b.fasta",
            "--mapping", "manifest.tsv", "--index", "reference.index",
        ])

    assert result.exit_code == 0, result.output
    argv = run.call_args.args[1]
    assert argv.count("--fasta") == 2
    assert argv[argv.index("--mapping") + 1] == "manifest.tsv"


def test_index_build_forwards_fasta_list():
    with (
        patch("metatracer.cli._which_or_die", return_value="/bin/mtsv-build"),
        patch("metatracer.cli._run") as run,
    ):
        result = CliRunner().invoke(cli, [
            "index-build", "--fasta-list", "references.txt",
            "--mapping", "manifest.tsv", "--index", "reference.index",
        ])

    assert result.exit_code == 0, result.output
    argv = run.call_args.args[1]
    assert argv[argv.index("--fasta-list") + 1] == "references.txt"


def test_index_build_requires_a_fasta_source():
    result = CliRunner().invoke(cli, ["index-build", "--index", "reference.index"])
    assert result.exit_code == 1
    assert "at least one --fasta or a --fasta-list" in result.output
