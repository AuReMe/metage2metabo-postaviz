from pathlib import Path

from m2m_postaviz import cli


def test_main_non_interactive_missing_required(monkeypatch, capsys):
    called = {"build": False}

    def fake_build_dataframes(*_args, **_kwargs):
        called["build"] = True

    monkeypatch.setattr(cli, "_is_interactive_terminal", lambda: False)
    monkeypatch.setattr(cli.du, "build_dataframes", fake_build_dataframes)

    cli.main([])

    out = capsys.readouterr().out
    assert "Required arguments missing" in out
    assert called["build"] is False


def test_main_interactive_precomputed_bypass(monkeypatch, tmp_path):
    load_dir = tmp_path / "precomputed"
    load_dir.mkdir()

    captured = {"load_path": None, "ran_shiny": False}

    class DummyDataStorage:
        def __init__(self, path):
            captured["load_path"] = Path(path)

    def fake_run_shiny(_data):
        captured["ran_shiny"] = True

    monkeypatch.setattr(cli, "_is_interactive_terminal", lambda: True)
    monkeypatch.setattr(cli, "_prompt_yes_no", lambda *args, **kwargs: True)
    monkeypatch.setattr(cli, "_prompt_required_path", lambda *args, **kwargs: load_dir)
    monkeypatch.setattr(cli, "DataStorage", DummyDataStorage)
    monkeypatch.setattr(cli, "_run_shiny", fake_run_shiny)

    cli.main([])

    assert captured["load_path"] == load_dir
    assert captured["ran_shiny"] is True


def test_main_interactive_required_and_optional_skip(monkeypatch, tmp_path):
    data_dir = tmp_path / "m2m_data"
    data_dir.mkdir()

    metadata_file = tmp_path / "metadata.tsv"
    metadata_file.write_text("sample\tgroup\nS1\tA\n", encoding="utf-8")

    output_dir = tmp_path / "out"

    captured = {
        "dir_path": None,
        "metadata_path": None,
        "abundance_path": None,
        "taxonomy_path": None,
        "save_path": None,
        "metacyc_path": None,
        "ran_shiny": False,
        "data_path": None,
    }

    def fake_build_dataframes(*args):
        captured["dir_path"] = args[0]
        captured["metadata_path"] = args[1]
        captured["abundance_path"] = args[2]
        captured["taxonomy_path"] = args[3]
        captured["save_path"] = args[4]
        captured["metacyc_path"] = args[5]

    class DummyDataStorage:
        def __init__(self, path):
            captured["data_path"] = Path(path)

    def fake_run_shiny(_data):
        captured["ran_shiny"] = True

    required_values = iter([data_dir, metadata_file, output_dir])

    monkeypatch.setattr(cli, "_is_interactive_terminal", lambda: True)
    monkeypatch.setattr(cli, "_prompt_yes_no", lambda *args, **kwargs: False)
    monkeypatch.setattr(cli, "_prompt_required_path", lambda *args, **kwargs: next(required_values))
    monkeypatch.setattr(cli, "_prompt_optional_path", lambda *args, **kwargs: None)

    monkeypatch.setattr(cli.du, "build_dataframes", fake_build_dataframes)
    monkeypatch.setattr(cli, "DataStorage", DummyDataStorage)
    monkeypatch.setattr(cli, "_run_shiny", fake_run_shiny)

    cli.main([])

    assert captured["dir_path"] == data_dir
    assert captured["metadata_path"] == metadata_file
    assert captured["abundance_path"] is None
    assert captured["taxonomy_path"] is None
    assert captured["save_path"] == output_dir
    assert captured["metacyc_path"] is None
    assert captured["data_path"] == output_dir
    assert captured["ran_shiny"] is True


def test_main_with_required_args_does_not_prompt_optional(monkeypatch, tmp_path):
    data_dir = tmp_path / "m2m_data"
    data_dir.mkdir()

    metadata_file = tmp_path / "metadata.tsv"
    metadata_file.write_text("sample\tgroup\nS1\tA\n", encoding="utf-8")

    output_dir = tmp_path / "out"

    captured = {
        "dir_path": None,
        "metadata_path": None,
        "abundance_path": None,
        "taxonomy_path": None,
        "save_path": None,
        "metacyc_path": None,
        "ran_shiny": False,
    }

    def fake_build_dataframes(*args):
        captured["dir_path"] = args[0]
        captured["metadata_path"] = args[1]
        captured["abundance_path"] = args[2]
        captured["taxonomy_path"] = args[3]
        captured["save_path"] = args[4]
        captured["metacyc_path"] = args[5]

    class DummyDataStorage:
        def __init__(self, _path):
            pass

    def fake_run_shiny(_data):
        captured["ran_shiny"] = True

    def fail_prompt(*_args, **_kwargs):
        raise AssertionError("Prompt helper should not be called when required args are provided")

    monkeypatch.setattr(cli, "_is_interactive_terminal", lambda: True)
    monkeypatch.setattr(cli, "_prompt_yes_no", fail_prompt)
    monkeypatch.setattr(cli, "_prompt_required_path", fail_prompt)
    monkeypatch.setattr(cli, "_prompt_optional_path", fail_prompt)

    monkeypatch.setattr(cli.du, "build_dataframes", fake_build_dataframes)
    monkeypatch.setattr(cli, "DataStorage", DummyDataStorage)
    monkeypatch.setattr(cli, "_run_shiny", fake_run_shiny)

    cli.main([
        "-d", str(data_dir),
        "-m", str(metadata_file),
        "-o", str(output_dir),
    ])

    assert captured["dir_path"] == data_dir
    assert captured["metadata_path"] == metadata_file
    assert captured["abundance_path"] is None
    assert captured["taxonomy_path"] is None


def test_main_graceful_exit_via_exit_command(monkeypatch, capsys):
    """Test that typing 'exit' exits gracefully without errors."""
    
    def fake_prompt_yes_no(*_args, **_kwargs):
        raise cli.InteractivePromptExit("User typed exit")
    
    monkeypatch.setattr(cli, "_is_interactive_terminal", lambda: True)
    monkeypatch.setattr(cli, "_prompt_yes_no", fake_prompt_yes_no)
    
    cli.main([])
    
    out = capsys.readouterr().out
    assert "Exiting interactive mode" in out


def test_main_graceful_exit_via_keyboard_interrupt(monkeypatch, capsys):
    """Test that Ctrl+C during prompting exits gracefully."""
    
    def fake_prompt_required(*_args, **_kwargs):
        raise KeyboardInterrupt()
    
    def fake_prompt_yes_no(*_args, **_kwargs):
        return False
    
    monkeypatch.setattr(cli, "_is_interactive_terminal", lambda: True)
    monkeypatch.setattr(cli, "_prompt_yes_no", fake_prompt_yes_no)
    monkeypatch.setattr(cli, "_prompt_required_path", fake_prompt_required)
    
    cli.main([])
    
    out = capsys.readouterr().out
    assert "Exiting..." in out
