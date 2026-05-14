from __future__ import annotations

import re
import time
from html import unescape
from typing import Callable, Dict, List
from urllib.parse import parse_qs, quote, urljoin, urlparse
from urllib.error import HTTPError
from urllib.request import Request, urlopen


NIST_SEARCH_URL = "https://webbook.nist.gov/cgi/cbook.cgi"
DEFAULT_TIMEOUT = 30
DEFAULT_MAX_RETRIES = 5
DEFAULT_RETRY_AFTER_SECONDS = 1
TRANSIENT_HTTP_STATUS_CODES = {429, 500, 502, 503, 504}
DEFAULT_HEADERS = {
    "User-Agent": (
        "Mozilla/5.0 (Windows NT 10.0; Win64; x64) "
        "AppleWebKit/537.36 (KHTML, like Gecko) "
        "Chrome/124.0.0.0 Safari/537.36"
    ),
    "Accept": "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
}


def _retry_after_seconds(error: HTTPError) -> int:
    header = None
    try:
        header = error.headers.get("Retry-After")
    except Exception:
        header = None
    if not header:
        return DEFAULT_RETRY_AFTER_SECONDS
    try:
        return max(DEFAULT_RETRY_AFTER_SECONDS, int(float(str(header).strip())))
    except Exception:
        return DEFAULT_RETRY_AFTER_SECONDS


def _exponential_backoff_seconds(attempt: int) -> int:
    return max(DEFAULT_RETRY_AFTER_SECONDS, 2 ** max(0, int(attempt)))


def default_fetch_text(
    url: str,
    *,
    timeout: int = DEFAULT_TIMEOUT,
    max_retries: int = DEFAULT_MAX_RETRIES,
) -> str:
    req = Request(url, headers=DEFAULT_HEADERS)
    last_error: HTTPError | None = None
    for attempt in range(max_retries):
        try:
            with urlopen(req, timeout=timeout) as resp:
                raw = resp.read()
            return raw.decode("utf-8", errors="ignore")
        except HTTPError as error:
            last_error = error
            status_code = int(getattr(error, "code", 0))
            if status_code not in TRANSIENT_HTTP_STATUS_CODES or attempt >= max_retries - 1:
                raise
            wait_retry_after = _retry_after_seconds(error) if status_code == 429 else DEFAULT_RETRY_AFTER_SECONDS
            wait_backoff = _exponential_backoff_seconds(attempt)
            time.sleep(max(wait_retry_after, wait_backoff))
    if last_error is not None:
        raise last_error
    raise RuntimeError("default_fetch_text exhausted retries without result")


def build_nist_search_url(inchi_str: str) -> str:
    encoded = quote(inchi_str, safe="")
    return f"{NIST_SEARCH_URL}?InChI={encoded}&Units=SI&Mask=80"


def build_nist_inchikey_url(inchikey: str) -> str:
    encoded = quote(inchikey, safe="")
    return f"https://webbook.nist.gov/cgi/inchi?InChIKey={encoded}&Mask=80"


def build_nist_direct_id_url(nist_id: str) -> str:
    encoded = quote(nist_id, safe="")
    return f"https://webbook.nist.gov/cgi/inchi?ID={encoded}&Mask=80"


def get_nist_page_by_inchi(
    inchi_str: str,
    *,
    fetch_text: Callable[[str], str] | None = None,
) -> Dict[str, str]:
    fetcher = fetch_text or default_fetch_text
    page_url = build_nist_search_url(inchi_str)
    html = fetcher(page_url)
    return {
        "page_url": page_url,
        "html": html,
    }


def get_nist_page_by_inchikey(
    inchikey: str,
    *,
    fetch_text: Callable[[str], str] | None = None,
) -> Dict[str, str]:
    fetcher = fetch_text or default_fetch_text
    page_url = build_nist_inchikey_url(inchikey)
    html = fetcher(page_url)
    return {
        "page_url": page_url,
        "html": html,
    }


def get_nist_page_by_id(
    nist_id: str,
    *,
    fetch_text: Callable[[str], str] | None = None,
) -> Dict[str, str]:
    fetcher = fetch_text or default_fetch_text
    page_url = build_nist_direct_id_url(nist_id)
    html = fetcher(page_url)
    return {
        "page_url": page_url,
        "html": html,
    }


def find_ir_jcamp_links(page_url: str, html: str) -> List[Dict[str, str]]:
    pattern = re.compile(r'href=["\']([^"\']*JCAMP=[^"\']+)["\']', re.IGNORECASE)
    irspec_pattern = re.compile(r'href=["\']([^"\']*Type=IR-SPEC[^"\']+)["\']', re.IGNORECASE)
    links: List[Dict[str, str]] = []
    seen = set()

    for href in pattern.findall(html or ""):
        full_url = urljoin(page_url, unescape(href))
        query = parse_qs(urlparse(full_url).query)
        if query.get("Type", [""])[0].upper() != "IR":
            continue

        jcamp_id = query.get("JCAMP", [""])[0]
        index = query.get("Index", ["0"])[0]
        key = (jcamp_id, index, full_url)
        if key in seen:
            continue
        seen.add(key)

        links.append(
            {
                "jcamp_url": full_url,
                "nist_id": jcamp_id,
                "index": index,
            }
        )

    for href in irspec_pattern.findall(html or ""):
        full_url = urljoin(page_url, unescape(href))
        query = parse_qs(urlparse(full_url).query)
        nist_id = query.get("ID", [""])[0]
        index = query.get("Index", ["0"])[0]
        if not nist_id:
            continue
        jcamp_url = f"https://webbook.nist.gov/cgi/inchi?JCAMP={quote(nist_id, safe='')}&Index={quote(index, safe='')}&Type=IR"
        key = (nist_id, index, jcamp_url)
        if key in seen:
            continue
        seen.add(key)
        links.append(
            {
                "jcamp_url": jcamp_url,
                "nist_id": nist_id,
                "index": index,
            }
        )

    return links


def extract_ir_spec_index_metadata(html: str) -> Dict[str, Dict[str, str]]:
    pattern = re.compile(
        r'href=["\'][^"\']*Type=IR-SPEC(?:&amp;|&)?Index=(\d+)[^"\']*["\'][^>]*>(.*?)</a>',
        re.IGNORECASE | re.DOTALL,
    )
    result: Dict[str, Dict[str, str]] = {}
    for index, label_html in pattern.findall(html or ""):
        label = re.sub(r"<[^>]+>", "", unescape(label_html)).strip()
        result[str(index)] = {
            "index": str(index),
            "description": " ".join(label.split()),
        }
    return result


def download_jcamp(
    jcamp_url: str,
    *,
    fetch_text: Callable[[str], str] | None = None,
) -> str:
    fetcher = fetch_text or default_fetch_text
    return fetcher(jcamp_url)
