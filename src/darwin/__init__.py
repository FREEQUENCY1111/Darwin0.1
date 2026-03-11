"""Darwin — Fast, accurate prokaryotic genome annotation."""

__app_name__ = "darwin-annotator"

try:
    from importlib.metadata import version

    __version__ = version(__app_name__)
except Exception:
    __version__ = "0.0.0-dev"
