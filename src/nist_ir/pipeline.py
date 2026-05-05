from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Callable, Dict, List
from urllib.error import HTTPError, URLError

from orca_parser import read_orca_hess

from .identifiers import hess_to_identifiers, smiles_to_identifiers
from .jcamp_parser import parse_jcamp_text
from .nist_client import (
    download_jcamp,
    extract_ir_spec_index_metadata,
    find_ir_jcamp_links,
    get_nist_page_by_id,
    get_nist_page_by_inchi,
    get_nist_page_by_inchikey,
)


def _safe_token(text: str) -> str:
    token = re.sub(r"[^A-Za-z0-9._+-]+", "_", str(text or "")).strip("_")
    return token or "unknown"


def _build_base_name(inchikey: str, nist_id: str, index: str) -> str:
    return f"{_safe_token(inchikey)}_{_safe_token(nist_id)}_IR_{_safe_token(index)}"


def _phase_tag_and_priority(description: str, metadata: Dict[str, str]) -> tuple[str, str, int]:
    text = " ".join(
        part for part in [
            str(description or "").lower(),
            str(metadata.get("STATE", "")).lower(),
        ]
        if part
    )
    if "gas" in text:
        return "gas", "gas", 100
    if "vapor" in text or "vapour" in text:
        return "vapor", "vapor", 90
    if "liquid" in text and "neat" in text:
        return "liquid_neat", "liquid (neat)", 80
    if "liquid" in text:
        return "liquid", "liquid", 70
    if "solid" in text:
        return "solid", "solid", 60
    if "solution" in text:
        return "solution", "solution", 40
    if metadata.get("STATE", "").strip():
        state = metadata["STATE"].strip().lower()
        return _safe_token(state), metadata["STATE"].strip(), 30
    return "unknown_phase", "unknown phase", 10


def nist_ir_from_smiles(
    smiles: str,
    outdir: str | Path,
    *,
    fetch_page_text: Callable[[str], str] | None = None,
    fetch_jcamp_text: Callable[[str], str] | None = None,
    nist_id: str | None = None,
) -> List[Dict[str, str]]:
    identifiers = smiles_to_identifiers(smiles)
    return nist_ir_from_identifiers(
        identifiers,
        outdir,
        fetch_page_text=fetch_page_text,
        fetch_jcamp_text=fetch_jcamp_text,
        nist_id=nist_id,
    )


def nist_ir_from_identifiers(
    identifiers: Dict[str, str],
    outdir: str | Path,
    *,
    fetch_page_text: Callable[[str], str] | None = None,
    fetch_jcamp_text: Callable[[str], str] | None = None,
    nist_id: str | None = None,
) -> List[Dict[str, str]]:
    page_attempts = []
    if nist_id:
        page_attempts.append(lambda: get_nist_page_by_id(nist_id, fetch_text=fetch_page_text))
    page_attempts.append(lambda: get_nist_page_by_inchi(identifiers["inchi"], fetch_text=fetch_page_text))
    page_attempts.append(lambda: get_nist_page_by_inchikey(identifiers["inchikey"], fetch_text=fetch_page_text))

    page = None
    matches = []
    page_errors: List[str] = []
    for attempt in page_attempts:
        try:
            candidate = attempt()
        except (HTTPError, URLError, LookupError, ValueError, RuntimeError, OSError) as exc:
            page_errors.append(f"{type(exc).__name__}:{exc}")
            continue
        candidate_matches = find_ir_jcamp_links(candidate["page_url"], candidate["html"])
        if candidate_matches:
            page = candidate
            matches = candidate_matches
            break

    if not matches:
        detail = "; ".join(page_errors) if page_errors else "no lookup errors reported"
        raise LookupError(f"No NIST IR JCAMP links found for this InChI; attempts={detail}")

    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    ir_index_meta = extract_ir_spec_index_metadata(page["html"]) if page else {}

    results: List[Dict[str, str]] = []
    for match in matches:
        jdx_text = download_jcamp(match["jcamp_url"], fetch_text=fetch_jcamp_text)
        metadata, spectrum = parse_jcamp_text(jdx_text)
        index_meta = ir_index_meta.get(str(match["index"]), {})
        phase_tag, phase_label, priority = _phase_tag_and_priority(index_meta.get("description", ""), metadata)

        base_name = f"{_build_base_name(identifiers['inchikey'], match['nist_id'], match['index'])}__{phase_tag}"
        jdx_path = outdir / f"{base_name}.jdx"
        csv_path = outdir / f"{base_name}.csv"
        meta_path = outdir / f"{base_name}_meta.json"

        jdx_path.write_text(jdx_text, encoding="utf-8")
        spectrum.to_csv(csv_path, index=False, encoding="utf-8")
        meta_path.write_text(
            json.dumps(
                {
                    "smiles": identifiers["input_smiles"],
                    "canonical_smiles": identifiers["canonical_smiles"],
                    "inchi": identifiers["inchi"],
                    "inchikey": identifiers["inchikey"],
                    "nist_page_url": page["page_url"],
                    "jcamp_url": match["jcamp_url"],
                    "nist_id": match["nist_id"],
                    "index": match["index"],
                    "phase_tag": phase_tag,
                    "phase_label": phase_label,
                    "selection_priority": priority,
                    "description": index_meta.get("description", ""),
                    "jcamp_metadata": metadata,
                },
                ensure_ascii=False,
                indent=2,
            ),
            encoding="utf-8",
        )

        results.append(
            {
                "jdx": str(jdx_path),
                "csv": str(csv_path),
                "meta_json": str(meta_path),
                "jcamp_url": match["jcamp_url"],
                "nist_id": match["nist_id"],
                "index": match["index"],
                "phase_tag": phase_tag,
                "phase_label": phase_label,
                "selection_priority": priority,
                "description": index_meta.get("description", ""),
            }
        )

    results = sorted(results, key=lambda item: (-int(item["selection_priority"]), str(item["index"])))
    manifest = {
        "smiles": identifiers["input_smiles"],
        "canonical_smiles": identifiers["canonical_smiles"],
        "inchi": identifiers["inchi"],
        "inchikey": identifiers["inchikey"],
        "nist_page_url": page["page_url"] if page else "",
        "reference_spectra": results,
        "preferred_reference": results[0] if results else None,
    }
    (outdir / f"{_safe_token(identifiers['inchikey'])}_reference_set.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )

    return results


def nist_ir_from_hess(
    hess_path: str | Path,
    outdir: str | Path,
    *,
    charge: int = 0,
    fetch_page_text: Callable[[str], str] | None = None,
    fetch_jcamp_text: Callable[[str], str] | None = None,
    nist_id: str | None = None,
) -> List[Dict[str, str]]:
    hess = read_orca_hess(hess_path)
    identifiers = hess_to_identifiers(hess, charge=charge)
    return nist_ir_from_identifiers(
        identifiers,
        outdir,
        fetch_page_text=fetch_page_text,
        fetch_jcamp_text=fetch_jcamp_text,
        nist_id=nist_id,
    )
