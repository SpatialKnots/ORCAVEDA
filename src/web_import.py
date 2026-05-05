from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json
import re
import shutil
from typing import Callable, Dict, List, Sequence
from uuid import uuid4
from urllib.error import HTTPError, URLError

import pandas as pd


PipelineRunner = Callable[[Sequence[Path], Path], Dict[str, pd.DataFrame]]
NistIrRunner = Callable[[Path, Path], Sequence[Dict[str, str]]]
RUN_ID_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._+-]{0,79}$")


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


def validate_web_run_id(run_id: str) -> str:
    candidate = str(run_id or "").strip()
    if not RUN_ID_PATTERN.fullmatch(candidate):
        raise ValueError(
            "Invalid run_id; expected 1-80 characters using letters, digits, dot, underscore, plus, or hyphen."
        )
    return candidate


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


def _default_nist_ir_runner(hess_path: Path, outdir: Path) -> Sequence[Dict[str, str]]:
    from nist_ir.pipeline import nist_ir_from_hess

    return nist_ir_from_hess(hess_path, outdir)


def _read_integration_manifest(output_dir: Path) -> Dict[str, str]:
    manifests = sorted(output_dir.glob("*__integration_manifest.json"))
    if not manifests:
        return {}
    manifest = json.loads(manifests[0].read_text(encoding="utf-8"))
    keys = ("xlsx_report", "interactive_spectrum_html", "interactive_spectrum_data_json")
    return {key: str(manifest[key]) for key in keys if key in manifest}


def _find_reference_set_manifest(nist_dir: Path) -> Path | None:
    manifests = sorted(nist_dir.glob("*_reference_set.json"))
    return manifests[0] if manifests else None


def _attach_nist_references_to_viewer(
    *,
    saved_paths: Sequence[Path],
    output_dir: Path,
    artifacts: Dict[str, str],
    nist_ir_runner: NistIrRunner,
) -> List[str]:
    diagnostics: List[str] = []
    spectrum_json = artifacts.get("interactive_spectrum_data_json", "")
    spectrum_html = artifacts.get("interactive_spectrum_html", "")
    if not spectrum_json or not spectrum_html:
        diagnostics.append("nist_ir_skipped:missing_interactive_spectrum_artifacts")
        return diagnostics

    from reports import attach_nist_reference_set, safe_output_stem, write_interactive_spectrum_viewer

    spectrum_json_path = Path(spectrum_json)
    spectrum_html_path = Path(spectrum_html)
    payload = json.loads(spectrum_json_path.read_text(encoding="utf-8"))
    attached_count = 0

    for hess_path in saved_paths:
        nist_dir = output_dir / "nist_ir" / safe_output_stem(hess_path.name)
        try:
            references = nist_ir_runner(hess_path, nist_dir)
        except (LookupError, ValueError, RuntimeError, OSError, HTTPError, URLError) as exc:
            diagnostics.append(f"nist_ir_failed:{hess_path.name}:{type(exc).__name__}:{exc}")
            continue

        reference_manifest = _find_reference_set_manifest(nist_dir)
        if reference_manifest is None:
            diagnostics.append(f"nist_ir_failed:{hess_path.name}:missing_reference_set_manifest")
            continue

        artifact_key = f"nist_reference_set_{safe_output_stem(hess_path.name)}"
        artifacts[artifact_key] = str(reference_manifest)
        payload = attach_nist_reference_set(
            payload,
            reference_manifest,
            file_title=safe_output_stem(hess_path.name),
        )
        attached_count += 1
        diagnostics.append(f"nist_ir_attached:{hess_path.name}:references={len(references)}")

    if attached_count:
        write_interactive_spectrum_viewer(payload, spectrum_html_path, json_path=spectrum_json_path)
        diagnostics.append(f"nist_ir_viewer_rewritten:attached_files={attached_count}")
    else:
        diagnostics.append("nist_ir_viewer_unchanged:no_reference_sets_attached")
    return diagnostics


def import_hess_files_for_web(
    source_paths: Sequence[str | Path],
    *,
    import_root: str | Path | None = None,
    run_id: str | None = None,
    pipeline_runner: PipelineRunner | None = None,
    include_nist_ir: bool = False,
    nist_ir_runner: NistIrRunner | None = None,
) -> WebHessImportResult:
    if not source_paths:
        raise ValueError("No .hess files were provided for web import.")

    resolved_run_id = validate_web_run_id(str(run_id or uuid4().hex))
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
    if include_nist_ir:
        diagnostics.extend(
            _attach_nist_references_to_viewer(
                saved_paths=saved_paths,
                output_dir=output_dir,
                artifacts=artifacts,
                nist_ir_runner=nist_ir_runner or _default_nist_ir_runner,
            )
        )
    else:
        diagnostics.append("nist_ir_skipped:not_requested")

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
