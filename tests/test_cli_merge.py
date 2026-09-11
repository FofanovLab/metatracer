from click.testing import CliRunner

from metatracer import cli as cli_module


def test_merge_hardcodes_taxid_gi_and_hides_mode(monkeypatch):
    invocation = {}
    monkeypatch.setattr(cli_module, "_which_or_die", lambda _name: "mtsv-collapse")
    monkeypatch.setattr(
        cli_module,
        "_run",
        lambda executable, arguments: invocation.update(
            executable=executable, arguments=arguments
        ),
    )

    result = CliRunner().invoke(
        cli_module.cli,
        ["merge", "--output", "merged.clp", "input.bn"],
    )

    assert result.exit_code == 0, result.output
    assert invocation["executable"] == "mtsv-collapse"
    assert invocation["arguments"][:2] == ["--mode", "taxid-gi"]

    help_result = CliRunner().invoke(cli_module.cli, ["merge", "--help"])
    assert help_result.exit_code == 0
    assert "--mode" not in help_result.output
