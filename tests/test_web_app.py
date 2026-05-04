from __future__ import annotations

import http.client
import json
import shutil
import sys
import threading
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from web_app import create_web_import_handler, parse_multipart_hess_upload, parse_multipart_upload  # noqa: E402
from http.server import ThreadingHTTPServer  # noqa: E402


def _multipart_body(filename: str, content: bytes, boundary: str = "orcaveda-test-boundary"):
    body = (
        f"--{boundary}\r\n"
        f'Content-Disposition: form-data; name="files"; filename="{filename}"\r\n'
        "Content-Type: application/octet-stream\r\n"
        "\r\n"
    ).encode("utf-8") + content + f"\r\n--{boundary}--\r\n".encode("utf-8")
    return f"multipart/form-data; boundary={boundary}", body


def _multipart_body_with_nist(filename: str, content: bytes, boundary: str = "orcaveda-test-boundary"):
    body = (
        f"--{boundary}\r\n"
        'Content-Disposition: form-data; name="include_nist_ir"\r\n'
        "\r\n"
        "1\r\n"
        f"--{boundary}\r\n"
        f'Content-Disposition: form-data; name="files"; filename="{filename}"\r\n'
        "Content-Type: application/octet-stream\r\n"
        "\r\n"
    ).encode("utf-8") + content + f"\r\n--{boundary}--\r\n".encode("utf-8")
    return f"multipart/form-data; boundary={boundary}", body


def test_parse_multipart_hess_upload():
    content_type, body = _multipart_body("H2O_freq.hess", b"$atoms\n")
    files = parse_multipart_hess_upload(content_type, body)

    assert len(files) == 1
    assert files[0].filename == "H2O_freq.hess"
    assert files[0].content == b"$atoms\n"


def test_parse_multipart_upload_reads_nist_option():
    content_type, body = _multipart_body_with_nist("H2O_freq.hess", b"$atoms\n")
    upload = parse_multipart_upload(content_type, body)

    assert upload.files[0].filename == "H2O_freq.hess"
    assert upload.fields["include_nist_ir"] == ("1",)


def test_web_import_page_embeds_interactive_viewer_after_upload():
    handler = create_web_import_handler()
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        conn = http.client.HTTPConnection("127.0.0.1", server.server_port, timeout=10)
        conn.request("GET", "/")
        response = conn.getresponse()
        html = response.read().decode("utf-8")

        assert response.status == 200
        assert "resultViewer" in html
        assert "Interactive spectrum viewer" in html
        assert "renderResult(payload)" in html
        assert "Diagnostics" in html
        assert "include_nist_ir" in html
        assert 'class="topbar"' in html
        assert 'class="artifact-strip"' in html
        assert 'class="diagnostics-toggle"' in html
        assert "height: calc(100vh - 156px)" in html
        assert "<h1>ORCAVEDA</h1>" in html
        assert "setRunningStatus()" in html
        assert "Parsing ORCA .hess data..." in html
        assert "<h2>Results</h2>" in html
        assert "<h2>Run <code>" not in html
    finally:
        server.shutdown()
        server.server_close()


def test_web_import_http_endpoint_with_fake_pipeline():
    import_root = ROOT / "outputs" / "pytest_web_app_imports"
    staging_root = ROOT / "outputs" / "pytest_web_app_staging"
    for path in (import_root, staging_root):
        if path.exists():
            shutil.rmtree(path)

    def fake_pipeline(paths, outdir):
        outdir.mkdir(parents=True, exist_ok=True)
        html_path = outdir / "fake_interactive_spectrum.html"
        json_path = outdir / "fake_spectrum_data.json"
        xlsx_path = outdir / "fake_report.xlsx"
        html_path.write_text("<html><body>viewer</body></html>", encoding="utf-8")
        json_path.write_text("{}", encoding="utf-8")
        xlsx_path.write_bytes(b"xlsx")
        manifest = {
            "xlsx_report": str(xlsx_path),
            "interactive_spectrum_html": str(html_path),
            "interactive_spectrum_data_json": str(json_path),
        }
        (outdir / "sample__integration_manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
        return {"assignment_audit": pd.DataFrame()}

    handler = create_web_import_handler(
        import_root=import_root,
        staging_root=staging_root,
        pipeline_runner=fake_pipeline,
    )
    server = ThreadingHTTPServer(("127.0.0.1", 0), handler)
    thread = threading.Thread(target=server.serve_forever, daemon=True)
    thread.start()
    try:
        conn = http.client.HTTPConnection("127.0.0.1", server.server_port, timeout=10)
        content_type, body = _multipart_body("H2O_freq.hess", b"$atoms\n")
        conn.request(
            "POST",
            "/api/hess/import?run_id=http_run",
            body=body,
            headers={"Content-Type": content_type, "Content-Length": str(len(body))},
        )
        response = conn.getresponse()
        payload = json.loads(response.read().decode("utf-8"))
        assert response.status == 200
        assert payload["run_id"] == "http_run"
        assert payload["artifact_urls"]["interactive_spectrum_html"] == "/api/runs/http_run/artifacts/interactive_spectrum_html"

        conn.request("GET", "/api/runs/http_run")
        manifest_response = conn.getresponse()
        manifest_payload = json.loads(manifest_response.read().decode("utf-8"))
        assert manifest_response.status == 200
        assert manifest_payload["input_files"] == ["H2O_freq.hess"]

        conn.request("GET", "/api/runs/http_run/artifacts/interactive_spectrum_html")
        artifact_response = conn.getresponse()
        artifact_body = artifact_response.read().decode("utf-8")
        assert artifact_response.status == 200
        assert "viewer" in artifact_body
    finally:
        server.shutdown()
        server.server_close()
