from __future__ import annotations

import sys
from pathlib import Path
from urllib.error import HTTPError


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from nist_ir import nist_client  # noqa: E402
from nist_ir.nist_client import DEFAULT_HEADERS  # noqa: E402


def test_default_nist_headers_include_browser_like_user_agent():
    assert "User-Agent" in DEFAULT_HEADERS
    assert "Mozilla/5.0" in DEFAULT_HEADERS["User-Agent"]
    assert "Accept" in DEFAULT_HEADERS
    assert "text/html" in DEFAULT_HEADERS["Accept"]


def test_default_fetch_text_retries_after_http_429(monkeypatch):
    calls = {"count": 0}
    waits = []

    class DummyResponse:
        def __enter__(self):
            return self

        def __exit__(self, exc_type, exc, tb):
            return False

        def read(self):
            return b"ok"

    def fake_urlopen(req, timeout=0):
        calls["count"] += 1
        if calls["count"] == 1:
            raise HTTPError(
                req.full_url,
                429,
                "Too Many Requests",
                {"Retry-After": "1"},
                None,
            )
        return DummyResponse()

    monkeypatch.setattr(nist_client, "urlopen", fake_urlopen)
    monkeypatch.setattr(nist_client.time, "sleep", lambda seconds: waits.append(seconds))

    result = nist_client.default_fetch_text("https://example.com/test", max_retries=2)
    assert result == "ok"
    assert calls["count"] == 2
    assert waits == [1]
