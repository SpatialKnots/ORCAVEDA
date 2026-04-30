from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from reports import attach_nist_reference_set, write_interactive_spectrum_viewer  # noqa: E402


def main() -> None:
    parser = argparse.ArgumentParser(description="Attach a NIST reference_set.json to an existing ORCAVEDA viewer payload.")
    parser.add_argument("spectrum_json", help="Existing ORCAVEDA __spectrum_data.json file")
    parser.add_argument("reference_set_json", help="NIST reference_set.json file")
    parser.add_argument("--file-title", help="Optional payload file title when the payload contains multiple files")
    parser.add_argument("--html-out", help="Output HTML path. Defaults next to the JSON with __with_nist suffix.")
    parser.add_argument("--json-out", help="Output JSON path. Defaults next to the JSON with __with_nist suffix.")
    args = parser.parse_args()

    spectrum_json = Path(args.spectrum_json)
    payload = json.loads(spectrum_json.read_text(encoding="utf-8"))
    payload = attach_nist_reference_set(payload, args.reference_set_json, file_title=args.file_title)

    default_json_out = spectrum_json.with_name(f"{spectrum_json.stem}__with_nist{spectrum_json.suffix}")
    json_out = Path(args.json_out) if args.json_out else default_json_out
    default_html_out = json_out.with_name(json_out.stem.replace("__spectrum_data", "__interactive_spectrum") + ".html")
    html_out = Path(args.html_out) if args.html_out else default_html_out

    write_interactive_spectrum_viewer(payload, html_out, json_path=json_out)
    print(f"JSON: {json_out}")
    print(f"HTML: {html_out}")


if __name__ == "__main__":
    main()
