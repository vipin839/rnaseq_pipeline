"""Regression tests for F0: credentials must never reach project files, reports, manifests or logs."""
import logging

import pytest

from rnaseq_pipeline import PipelineError, logger, net, secrets, ui
from rnaseq_pipeline import config as C
from rnaseq_pipeline.project import Project

KEY = "ecdbTESTSECRET0123456789abcdef"


@pytest.fixture(autouse=True)
def registered():
    secrets.register(KEY)
    yield


def cfg_with_key():
    cfg = C.load()
    cfg["ncbi"] = {"email": "someone@example.org", "api_key": KEY}
    return cfg


def test_new_project_never_stores_credentials(tmp_path):
    p = Project.create("sec", tmp_path, base_config=cfg_with_key())
    text = p.config_path.read_text()
    assert KEY not in text and "someone@example.org" not in text


def test_old_project_is_cleaned_on_open(tmp_path):
    p = Project.create("old", tmp_path)
    C.save_yaml(cfg_with_key(), p.config_path)          # simulate a project written by an older version
    assert KEY in p.config_path.read_text()
    Project.open(p.root)
    assert KEY not in p.config_path.read_text()


def test_save_config_strips_credentials(tmp_path):
    p = Project.create("s2", tmp_path)
    p.save_config(cfg_with_key())
    assert KEY not in p.config_path.read_text()


def test_config_snapshot_is_redacted(tmp_path):
    path = C.snapshot(cfg_with_key(), tmp_path)
    text = path.read_text()
    assert KEY not in text and "someone@example.org" not in text and "redacted" in text


def test_redacted_copy_leaves_original_intact():
    cfg = cfg_with_key()
    r = secrets.redacted(cfg)
    assert KEY not in str(r) and cfg["ncbi"]["api_key"] == KEY


def test_network_error_does_not_echo_key():
    with pytest.raises(PipelineError) as e:
        net.get_text(f"http://127.0.0.1:9/x?db=sra&api_key={KEY}", retries=1, timeout=2)
    assert KEY not in str(e.value) and "api_key=***" in str(e.value)


def test_logs_and_terminal_are_scrubbed(tmp_path, capsys):
    logger.attach_project(tmp_path / "logs")
    try:
        ui.info(f"using key {KEY} for the search")
        logging.getLogger("rnaseq").error("raw record with %s", KEY)
        logger.record_command({"stage": "t", "command_str": f"curl https://x?api_key={KEY}"})
    finally:
        logger.detach_project()
    assert KEY not in capsys.readouterr().out
    for f in (tmp_path / "logs").iterdir():
        assert KEY not in f.read_text(), f.name


def test_url_parameters_scrubbed_even_if_unregistered():
    assert secrets.scrub("https://e.utils/x?api_key=abcdef123456&db=sra") == "https://e.utils/x?api_key=***&db=sra"
