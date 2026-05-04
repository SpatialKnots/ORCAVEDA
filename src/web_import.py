from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json
import re
import shutil
from typing import Callable, Dict, Sequence
from uuid import uuid4

import pandas as pd


PipelineRunner = Callable[[Sequence[Path], Path], Dict[str, pd.DataFrame]]


@dataclass(frozen=True)
class WebHessImportResult:
    run_id: str
    status: str
    input_files: tuple[str, ...]
    input_paths: tuple[str, ...]
    output_dir: str
    artifacts: Dict[str, str]
    diagnostics: tuple[str, ...]

    def as_dict(self) -> Dict[str, object]:
        return {
            "run_id": self.run_id,
            "status": self.status,
            "input_files": list(self.input_files),
            "input_paths": list(self.input_paths),
            "output_dir": self.output_dir,
            "artifacts": dict(self.artifacts),
            "diagnostics": list(self.diagnostics),
        }


def default_web_import_root() -> Path:
    return Path("outputs") / "web_imports"


def _safe_upload_name(filename: str) -> str:
    name = Path(str(filename)).name.strip()
    if not name:
        raise ValueError("Uploaded .hess filename is empty.")
    clean = re.sub(r"[^A-Za-z0-9._+-]+", "_", name).strip("_")
    if not clean:
        raise ValueError(f"Uploaded .hess filename cannot be made safe: {filename!r}")
    return clean


def _unique_destination(directory: Path, filename: str) -> Path:
    target = directory / filename
    if not target.exists():
        return target
    stem = target.stem
    suffix = target.suffix
    counter = 2
    while True:
        candidate = directory / f"{stem}_{counter}{suffix}"
        if not candidate.exists():
            return candidate
        counter += 1


def _default_pipeline_runner(paths: Sequence[Path], outdir: Path) -> Dict[str, pd.DataFrame]:
    from ORCAVEDA_patched_stage3D_v5_0 import analyze_orca_ped_like

    return analyze_orca_ped_like(paths, outdir)


def _read_integration_manifest(output_dir: Path) -> Dict[str, str]:
    manifests = sorted(output_dir.glob("*__integration_manifest.json"))
    if not manifests:
        return {}
    manifest = json.loads(manifests[0].read_text(encoding="utf-8"))
    keys = ("xlsx_report", "interactive_spectrum_html", "interactive_spectrum_data_json")
    return {key: str(manifest[key]) for key in keys if key in manifest}


def import_hess_files_for_web(
    source_paths: Sequence[str | Path],
    *,
    import_root: str | Path | None = None,
    run_id: str | None = None,
    pipeline_runner: PipelineRunner | None = None,
) -> WebHessImportResult:
    if not source_paths:
        raise ValueError("No .hess files were provided for web import.")

    resolved_run_id = str(run_id or uuid4().hex)
    root = Path(import_root) if import_root is not None else default_web_import_root()
    run_dir = root / resolved_run_id
    input_dir = run_dir / "inputs"
    output_dir = run_dir / "outputs"
    input_dir.mkdir(parents=True, exist_ok=False)
    output_dir.mkdir(parents=True, exist_ok=False)

    saved_paths: list[Path] = []
    diagnostics: list[str] = []
    for source in source_paths:
        source_path = Path(source)
        if source_path.suffix.lower() != ".hess":
            raise ValueError(f"Unsupported upload file extension for {source_path.name!r}; expected .hess.")
        if not source_path.exists():
            raise FileNotFoundError(f"Uploaded .hess source file does not exist: {source_path}")
        if not source_path.is_file():
            raise ValueError(f"Uploaded .hess source path is not a file: {source_path}")
        safe_name = _safe_upload_name(source_path.name)
        destination = _unique_destination(input_dir, safe_name)
        shutil.copy2(source_path, destination)
        saved_paths.append(destination)
        diagnostics.append(f"accepted_hess_upload:{destination.name}")

    runner = pipeline_runner or _default_pipeline_runner
    tables = runner(saved_paths, output_dir)
    artifacts = _read_integration_manifest(output_dir)
    artifacts["run_manifest_json"] = str(run_dir / "web_import_manifest.json")

    result = WebHessImportResult(
        run_id=resolved_run_id,
        status="completed",
        input_files=tuple(path.name for path in saved_paths),
        input_paths=tuple(str(path) for path in saved_paths),
        output_dir=str(output_dir),
        artifacts=artifacts,
        diagnostics=tuple(diagnostics + [f"pipeline_tables:{','.join(sorted(tables.keys()))}"]),
    )
    (run_dir / "web_import_manifest.json").write_text(json.dumps(result.as_dict(), indent=2), encoding="utf-8")
    return result
