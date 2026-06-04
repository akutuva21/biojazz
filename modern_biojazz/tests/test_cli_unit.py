from __future__ import annotations

import json
import os
import sys
from unittest.mock import patch

import pytest

from modern_biojazz.cli import main


def test_cli_default_run(fixtures_dir, capsys):
    seed_path = fixtures_dir / "seed_network.json"
    grounding_path = fixtures_dir / "grounding_payload.json"

    test_args = [
        "cli.py",
        "--seed",
        str(seed_path),
        "--grounding",
        str(grounding_path),
        "--generations",
        "1",
        "--population",
        "2",
        "--sim-t-end",
        "5.0",
    ]

    with patch.object(sys, "argv", test_args):
        main()

    captured = capsys.readouterr()
    output = json.loads(captured.out)

    assert "best_score" in output
    assert "history" in output
    assert "best_network" in output
    assert "grounding" in output
    assert output["grounding"] is not None


def test_cli_http_backend_missing_base_url(fixtures_dir):
    seed_path = fixtures_dir / "seed_network.json"

    test_args = [
        "cli.py",
        "--seed",
        str(seed_path),
        "--sim-backend",
        "http",
    ]

    with patch.object(sys, "argv", test_args):
        with pytest.raises(ValueError, match="--sim-base-url is required when --sim-backend=http"):
            main()


def test_cli_openai_compatible_missing_base_url(fixtures_dir):
    seed_path = fixtures_dir / "seed_network.json"

    test_args = [
        "cli.py",
        "--seed",
        str(seed_path),
        "--llm-provider",
        "openai_compatible",
    ]

    with patch.object(sys, "argv", test_args):
        with pytest.raises(ValueError, match="--llm-base-url is required when --llm-provider=openai_compatible"):
            main()


def test_cli_openai_compatible_missing_api_key(fixtures_dir, monkeypatch):
    seed_path = fixtures_dir / "seed_network.json"

    monkeypatch.delenv("OPENAI_API_KEY", raising=False)

    test_args = [
        "cli.py",
        "--seed",
        str(seed_path),
        "--llm-provider",
        "openai_compatible",
        "--llm-base-url",
        "http://localhost:8000",
        "--llm-api-key-env",
        "OPENAI_API_KEY",
    ]

    with patch.object(sys, "argv", test_args):
        with pytest.raises(ValueError, match="Environment variable OPENAI_API_KEY must be set for llm provider"):
            main()
