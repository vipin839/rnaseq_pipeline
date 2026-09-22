"""Keep credentials out of project files, reports, manifests, terminal output and logs.

User-level settings such as the NCBI API key (and the email NCBI asks for) belong to the user, not to a
project: they are never written into a project, and any registered secret value is scrubbed from text
before it is printed or logged.
"""
import copy
import logging
import re

# (section, key) pairs that are credentials or personal data and must not be persisted in projects.
SECRET_KEYS = {("ncbi", "api_key")}
PERSONAL_KEYS = {("ncbi", "email")}

_values = set()
_URL_SECRET_RE = re.compile(r"(api_key|apikey|token|password)=([^&\s'\"]+)", re.I)


def register(value):
    """Remember a secret value so it is scrubbed from any text shown or logged."""
    if value and isinstance(value, str) and len(value) >= 6:
        _values.add(value)


def register_from_config(cfg):
    for section, key in SECRET_KEYS:
        register(((cfg or {}).get(section) or {}).get(key))


def scrub(text):
    """Replace registered secrets and credential-looking URL parameters with ***."""
    if not isinstance(text, str) or not text:
        return text
    for v in _values:
        if v in text:
            text = text.replace(v, "***")
    return _URL_SECRET_RE.sub(lambda m: f"{m.group(1)}=***", text)


def strip_user_settings(cfg):
    """Copy of a configuration without credentials/personal settings (what a project may store)."""
    out = copy.deepcopy(cfg or {})
    for section, key in SECRET_KEYS | PERSONAL_KEYS:
        if isinstance(out.get(section), dict):
            out[section].pop(key, None)
    return out


def redacted(cfg):
    """Copy of a configuration safe to show or archive: secrets masked, personal data marked."""
    out = copy.deepcopy(cfg or {})
    for section, key in SECRET_KEYS:
        if isinstance(out.get(section), dict) and out[section].get(key):
            out[section][key] = "*** (redacted)"
    for section, key in PERSONAL_KEYS:
        if isinstance(out.get(section), dict) and out[section].get(key):
            out[section][key] = "(set — not stored)"
    return out


class ScrubFilter(logging.Filter):
    """Logging filter that scrubs secrets from every record before it is written."""

    def filter(self, record):
        msg = record.getMessage()
        clean = scrub(msg)
        if clean != msg:
            record.msg, record.args = clean, ()
        return True
