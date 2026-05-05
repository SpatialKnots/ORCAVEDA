from __future__ import annotations

from dataclasses import dataclass
from email import policy
from email.parser import BytesParser
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
import html
import json
import mimetypes
import re
import shutil
from typing import Dict, Sequence
from urllib.parse import parse_qs, urlparse
from uuid import uuid4

from web_import import (
    NistIrRunner,
    PipelineRunner,
    default_web_import_root,
    import_hess_files_for_web,
    validate_web_run_id,
)


@dataclass(frozen=True)
class UploadedFile:
    filename: str
    content: bytes


@dataclass(frozen=True)
class MultipartUpload:
    files: tuple[UploadedFile, ...]
    fields: Dict[str, tuple[str, ...]]


def parse_multipart_upload(content_type: str, body: bytes) -> MultipartUpload:
    if not content_type.lower().startswith("multipart/form-data"):
        raise ValueError("Expected multipart/form-data upload.")

    header = f"Content-Type: {content_type}\r\nMIME-Version: 1.0\r\n\r\n".encode("utf-8")
    message = BytesParser(policy=policy.default).parsebytes(header + body)
    if not message.is_multipart():
        raise ValueError("Upload body is not multipart.")

    files: list[UploadedFile] = []
    fields: Dict[str, list[str]] = {}
    for part in message.iter_parts():
        disposition = str(part.get_content_disposition() or "")
        field_name = str(part.get_param("name", header="content-disposition") or "")
        filename = str(part.get_filename() or "")
        if disposition != "form-data" or not field_name:
            continue
        if filename:
            if field_name not in {"files", "file", "hess"}:
                continue
            payload = part.get_payload(decode=True) or b""
            if not payload:
                raise ValueError(f"Uploaded .hess file is empty: {filename}")
            files.append(UploadedFile(filename=filename, content=payload))
        else:
            value = str(part.get_content() or "")
            fields.setdefault(field_name, []).append(value)

    if not files:
        raise ValueError("No .hess files were found in the upload form.")
    return MultipartUpload(
        files=tuple(files),
        fields={key: tuple(values) for key, values in fields.items()},
    )


def parse_multipart_hess_upload(content_type: str, body: bytes) -> tuple[UploadedFile, ...]:
    return parse_multipart_upload(content_type, body).files


def _truthy_form_field(fields: Dict[str, tuple[str, ...]], key: str) -> bool:
    values = fields.get(key, ())
    return any(str(value).strip().lower() in {"1", "true", "yes", "on"} for value in values)


def _safe_staging_filename(filename: str) -> str:
    name = Path(str(filename)).name.strip()
    clean = re.sub(r"[^A-Za-z0-9._+-]+", "_", name).strip("_")
    if not clean:
        raise ValueError(f"Uploaded filename cannot be made safe: {filename!r}")
    return clean


def _write_staging_files(files: Sequence[UploadedFile], staging_dir: Path) -> tuple[Path, ...]:
    staging_dir.mkdir(parents=True, exist_ok=False)
    paths: list[Path] = []
    seen: Dict[str, int] = {}
    for item in files:
        safe_name = _safe_staging_filename(item.filename)
        suffix_index = seen.get(safe_name, 0)
        seen[safe_name] = suffix_index + 1
        if suffix_index:
            path_obj = Path(safe_name)
            safe_name = f"{path_obj.stem}_{suffix_index + 1}{path_obj.suffix}"
        target = staging_dir / safe_name
        target.write_bytes(item.content)
        paths.append(target)
    return tuple(paths)


def _artifact_urls(run_id: str, artifacts: Dict[str, str]) -> Dict[str, str]:
    return {key: f"/api/runs/{run_id}/artifacts/{key}" for key in sorted(artifacts)}


def _result_payload(result) -> Dict[str, object]:
    payload = result.as_dict()
    payload["artifact_urls"] = _artifact_urls(result.run_id, result.artifacts)
    return payload


def _index_html() -> str:
    return """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>ORCAVEDA .hess Import</title>
  <style>
    body { font-family: system-ui, sans-serif; margin: 0; background: #f7f7f4; color: #202124; }
    main { width: min(1680px, calc(100vw - 32px)); margin: 16px auto; padding: 0 0 24px; }
    .topbar { display: grid; grid-template-columns: auto minmax(320px, 1fr) minmax(220px, auto); gap: 14px; align-items: center; background: #fff; border: 1px solid #d7d7d2; border-radius: 6px; padding: 10px 12px; position: sticky; top: 8px; z-index: 5; box-shadow: 0 1px 2px rgba(0,0,0,0.04); }
    h1 { font-size: 20px; margin: 0; white-space: nowrap; }
    form { display: flex; align-items: center; gap: 10px; min-width: 0; }
    .file-label { position: absolute; width: 1px; height: 1px; overflow: hidden; clip: rect(0, 0, 0, 0); }
    input[type=file] { min-width: 220px; max-width: 420px; }
    .nist-option { display: inline-flex; align-items: center; gap: 6px; white-space: nowrap; color: #374151; font-size: 14px; }
    button { background: #155e75; color: #fff; border: 0; border-radius: 4px; padding: 8px 12px; cursor: pointer; white-space: nowrap; }
    button:disabled { opacity: 0.65; cursor: wait; }
    .status { min-width: 220px; color: #374151; font-size: 14px; text-align: right; }
    .status.running { display: flex; align-items: center; justify-content: flex-end; gap: 8px; color: #155e75; }
    .status.running::before { content: ""; width: 14px; height: 14px; border: 2px solid #cfe3e8; border-top-color: #155e75; border-radius: 50%; animation: spin 0.8s linear infinite; }
    .progress-line { display: block; height: 3px; margin-top: 6px; overflow: hidden; border-radius: 999px; background: #e8ece8; }
    .progress-line::before { content: ""; display: block; width: 42%; height: 100%; border-radius: inherit; background: #155e75; animation: progressSlide 1.35s ease-in-out infinite; }
    @keyframes spin { to { transform: rotate(360deg); } }
    @keyframes progressSlide { 0% { transform: translateX(-110%); } 55% { transform: translateX(90%); } 100% { transform: translateX(260%); } }
    .error { color: #9f1239; }
    .result { margin-top: 12px; }
    .run-header { display: flex; flex-wrap: wrap; gap: 8px 16px; align-items: baseline; background: #fff; border: 1px solid #d7d7d2; border-radius: 6px 6px 0 0; padding: 10px 12px; }
    .run-header h2 { margin: 0; font-size: 18px; }
    .meta { color: #4b5563; font-size: 14px; }
    .artifact-strip { display: flex; flex-wrap: wrap; gap: 8px; align-items: center; background: #fff; border-left: 1px solid #d7d7d2; border-right: 1px solid #d7d7d2; padding: 8px 12px; }
    .artifact-strip a, .diagnostics-toggle { border: 1px solid #d7d7d2; border-radius: 999px; padding: 4px 9px; background: #f7f7f4; color: #155e75; text-decoration: none; font-size: 13px; }
    .diagnostics-toggle { cursor: pointer; font-family: inherit; }
    .diagnostics { background: #fff; border-left: 1px solid #d7d7d2; border-right: 1px solid #d7d7d2; padding: 8px 12px; overflow: auto; max-height: 120px; font-size: 12px; }
    .viewer-frame { display: block; width: 100%; height: calc(100vh - 156px); min-height: 760px; border: 1px solid #d7d7d2; border-radius: 0 0 6px 6px; background: #fff; }
    @media (max-width: 980px) {
      main { width: min(100vw - 20px, 1680px); margin-top: 10px; }
      .topbar { grid-template-columns: 1fr; position: static; }
      form { flex-wrap: wrap; }
      .status { text-align: left; }
      .viewer-frame { height: 760px; min-height: 680px; }
    }
    code { background: #eeeeea; padding: 1px 4px; border-radius: 3px; }
  </style>
</head>
<body>
<main>
  <div class="topbar">
    <h1>ORCAVEDA</h1>
    <form id="uploadForm">
      <label class="file-label" for="files">ORCA .hess files</label>
      <input id="files" name="files" type="file" accept=".hess" multiple required>
      <label class="nist-option"><input id="includeNist" name="include_nist_ir" type="checkbox"> NIST IR</label>
      <button id="submitButton" type="submit">Run Import</button>
    </form>
    <div id="status" class="status"></div>
  </div>
  <section id="result" class="result" hidden></section>
</main>
<script>
const form = document.getElementById("uploadForm");
const filesInput = document.getElementById("files");
const statusBox = document.getElementById("status");
const resultBox = document.getElementById("result");
const submitButton = document.getElementById("submitButton");

const artifactLabels = {
  interactive_spectrum_html: "Interactive spectrum viewer",
  interactive_spectrum_data_json: "Spectrum data JSON",
  run_manifest_json: "Run manifest JSON",
  xlsx_report: "XLSX report",
};

let statusTimer = null;
const runningMessages = [
  "Uploading .hess file...",
  "Parsing ORCA .hess data...",
  "Running assignment audit...",
  "Building interactive viewer...",
  "Attaching NIST IR references...",
  "Preparing results...",
];

function escapeHtml(value) {
  return String(value ?? "").replace(/[&<>"']/g, (char) => ({
    "&": "&amp;",
    "<": "&lt;",
    ">": "&gt;",
    '"': "&quot;",
    "'": "&#39;",
  }[char]));
}

function renderResult(payload) {
  const artifactUrls = payload.artifact_urls || {};
  const viewerUrl = artifactUrls.interactive_spectrum_html || "";
  const inputFiles = (payload.input_files || []).map(escapeHtml).join(", ");
  const diagnosticsItems = [`run_id:${payload.run_id || "not_reported"}`, ...(payload.diagnostics || [])];
  const diagnostics = diagnosticsItems.map((item) => `<div>${escapeHtml(item)}</div>`).join("");
  const links = Object.entries(artifactUrls)
    .map(([key, url]) => {
      const label = artifactLabels[key] || (key.startsWith("nist_reference_set_") ? `NIST reference set (${key.replace("nist_reference_set_", "")})` : key);
      return `<a href="${url}" target="_blank" rel="noopener">${escapeHtml(label)}</a>`;
    })
    .join("");
  const viewerBlock = viewerUrl
    ? `<iframe id="resultViewer" class="viewer-frame" src="${viewerUrl}" title="ORCAVEDA interactive spectrum viewer"></iframe>`
    : `<div class="error">Interactive spectrum viewer artifact was not reported.</div>`;

  resultBox.innerHTML = `
    <div class="run-header">
      <h2>Results</h2>
      <div class="meta">${inputFiles || "Input file not reported"}</div>
      ${viewerUrl ? `<a href="${viewerUrl}" target="_blank" rel="noopener">Open viewer in a new tab</a>` : ""}
    </div>
    <div class="artifact-strip">
      ${links}
      <button class="diagnostics-toggle" type="button" onclick="document.getElementById('diagnosticsPanel').hidden = !document.getElementById('diagnosticsPanel').hidden">Diagnostics</button>
    </div>
    <div id="diagnosticsPanel" class="diagnostics" hidden>${diagnostics || "No diagnostics reported."}</div>
    ${viewerBlock}
  `;
  resultBox.hidden = false;
}

function setIdleStatus(message) {
  clearInterval(statusTimer);
  statusTimer = null;
  statusBox.className = "status";
  statusBox.textContent = message || "";
}

function setRunningStatus() {
  clearInterval(statusTimer);
  let index = 0;
  const update = () => {
    const message = runningMessages[index % runningMessages.length];
    statusBox.className = "status running";
    statusBox.innerHTML = `${escapeHtml(message)}<span class="progress-line"></span>`;
    index += 1;
  };
  update();
  statusTimer = setInterval(update, 1800);
}

form.addEventListener("submit", async (event) => {
  event.preventDefault();
  resultBox.hidden = true;
  resultBox.innerHTML = "";
  setRunningStatus();
  submitButton.disabled = true;
  const data = new FormData();
  for (const file of filesInput.files) data.append("files", file);
  if (document.getElementById("includeNist").checked) data.append("include_nist_ir", "1");
  try {
    const response = await fetch("/api/hess/import", { method: "POST", body: data });
    const payload = await response.json();
    if (!response.ok) throw new Error(payload.error || "Import failed.");
    renderResult(payload);
    setIdleStatus("Import completed.");
  } catch (error) {
    clearInterval(statusTimer);
    statusTimer = null;
    statusBox.className = "status error";
    statusBox.textContent = error.message;
  } finally {
    submitButton.disabled = false;
  }
});
</script>
</body>
</html>
"""


def create_web_import_handler(
    *,
    import_root: str | Path | None = None,
    staging_root: str | Path | None = None,
    pipeline_runner: PipelineRunner | None = None,
    nist_ir_runner: NistIrRunner | None = None,
    max_upload_bytes: int = 100 * 1024 * 1024,
):
    resolved_import_root = Path(import_root) if import_root is not None else default_web_import_root()
    resolved_staging_root = Path(staging_root) if staging_root is not None else Path("outputs") / "web_upload_staging"

    class ORCAVEDAWebImportHandler(BaseHTTPRequestHandler):
        server_version = "ORCAVEDAWebImport/0.1"

        def log_message(self, format, *args):
            return

        def _send_bytes(self, status: int, content: bytes, content_type: str) -> None:
            self.send_response(status)
            self.send_header("Content-Type", content_type)
            self.send_header("Content-Length", str(len(content)))
            self.end_headers()
            self.wfile.write(content)

        def _send_json(self, status: int, payload: Dict[str, object]) -> None:
            self._send_bytes(status, json.dumps(payload, indent=2).encode("utf-8"), "application/json; charset=utf-8")

        def do_GET(self) -> None:
            parsed = urlparse(self.path)
            if parsed.path in {"/", "/import"}:
                self._send_bytes(200, _index_html().encode("utf-8"), "text/html; charset=utf-8")
                return
            if parsed.path.startswith("/api/runs/"):
                parts = [part for part in parsed.path.split("/") if part]
                if len(parts) == 3:
                    self._send_run_manifest(parts[2])
                    return
                if len(parts) == 5 and parts[3] == "artifacts":
                    self._send_artifact(parts[2], parts[4])
                    return
            self._send_json(404, {"error": "Not found."})

        def do_POST(self) -> None:
            parsed = urlparse(self.path)
            if parsed.path != "/api/hess/import":
                self._send_json(404, {"error": "Not found."})
                return

            length = int(self.headers.get("Content-Length", "0"))
            if length <= 0:
                self._send_json(400, {"error": "Upload body is empty."})
                return
            if length > max_upload_bytes:
                self._send_json(413, {"error": f"Upload exceeds {max_upload_bytes} bytes."})
                return

            try:
                run_id = validate_web_run_id(parse_qs(parsed.query).get("run_id", [uuid4().hex])[0])
            except ValueError as exc:
                self._send_json(400, {"error": str(exc)})
                return
            staging_dir = resolved_staging_root / run_id
            try:
                body = self.rfile.read(length)
                upload = parse_multipart_upload(str(self.headers.get("Content-Type", "")), body)
                staged_paths = _write_staging_files(upload.files, staging_dir)
                result = import_hess_files_for_web(
                    staged_paths,
                    import_root=resolved_import_root,
                    run_id=run_id,
                    pipeline_runner=pipeline_runner,
                    include_nist_ir=_truthy_form_field(upload.fields, "include_nist_ir"),
                    nist_ir_runner=nist_ir_runner,
                )
            except (FileExistsError, FileNotFoundError, ValueError) as exc:
                self._send_json(400, {"error": str(exc)})
                return
            finally:
                if staging_dir.exists():
                    shutil.rmtree(staging_dir)

            self._send_json(200, _result_payload(result))

        def _send_run_manifest(self, run_id: str) -> None:
            try:
                run_id = validate_web_run_id(run_id)
            except ValueError as exc:
                self._send_json(400, {"error": str(exc)})
                return
            manifest_path = resolved_import_root / run_id / "web_import_manifest.json"
            if not manifest_path.exists():
                self._send_json(404, {"error": f"Run not found: {run_id}"})
                return
            payload = json.loads(manifest_path.read_text(encoding="utf-8"))
            payload["artifact_urls"] = _artifact_urls(run_id, dict(payload.get("artifacts", {})))
            self._send_json(200, payload)

        def _send_artifact(self, run_id: str, artifact_key: str) -> None:
            try:
                run_id = validate_web_run_id(run_id)
            except ValueError as exc:
                self._send_json(400, {"error": str(exc)})
                return
            manifest_path = resolved_import_root / run_id / "web_import_manifest.json"
            if not manifest_path.exists():
                self._send_json(404, {"error": f"Run not found: {run_id}"})
                return
            payload = json.loads(manifest_path.read_text(encoding="utf-8"))
            artifacts = dict(payload.get("artifacts", {}))
            artifact_path_text = artifacts.get(artifact_key)
            if not artifact_path_text:
                self._send_json(404, {"error": f"Artifact not found: {html.escape(artifact_key)}"})
                return
            artifact_path = Path(str(artifact_path_text)).resolve()
            run_dir = (resolved_import_root / run_id).resolve()
            if run_dir not in artifact_path.parents and artifact_path != run_dir:
                self._send_json(403, {"error": "Artifact path is outside the run directory."})
                return
            if not artifact_path.exists() or not artifact_path.is_file():
                self._send_json(404, {"error": f"Artifact file is missing: {artifact_key}"})
                return
            content_type = mimetypes.guess_type(str(artifact_path))[0] or "application/octet-stream"
            self._send_bytes(200, artifact_path.read_bytes(), content_type)

    return ORCAVEDAWebImportHandler


def run_web_import_server(host: str = "127.0.0.1", port: int = 8765) -> None:
    handler = create_web_import_handler()
    server = ThreadingHTTPServer((host, port), handler)
    print(f"ORCAVEDA web import server: http://{host}:{port}/")
    server.serve_forever()


if __name__ == "__main__":
    run_web_import_server()
