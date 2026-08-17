"""Bundled resources (example models, icons)."""

from __future__ import annotations

import os

_HERE = os.path.dirname(os.path.abspath(__file__))


def resource_path(*parts: str) -> str:
    """Absolute path to a bundled resource (works in source and PyInstaller builds)."""
    return os.path.join(_HERE, *parts)


def example_model_path() -> str:
    """Path to the bundled E. coli core example model."""
    return resource_path("models", "e_coli_core.xml")


def app_icon_path() -> str:
    """Path to the application launcher icon (.ico, used as the window icon too)."""
    return resource_path("icons", "app_icon.ico")


def offline_universal_path() -> str:
    """Path to the bundled offline universal reaction database (Issue 4).

    Curated, cytosol-collapsed, mass-balance-filtered universal that ships with
    the app so Pathway Design works with no internet. May not exist in a checkout
    that hasn't run tools/build_offline_universal.py."""
    return resource_path("databases", "offline_universal.json")


def has_offline_universal() -> bool:
    return os.path.exists(offline_universal_path())


def escher_host_html() -> str:
    """Path to the QtWebEngine host page that renders interactive Escher maps."""
    return resource_path("escher", "escher_host.html")


def has_escher_assets() -> bool:
    """Whether the bundled Escher JS engine + host page are present."""
    return (os.path.exists(escher_host_html())
            and os.path.exists(resource_path("escher", "escher.min.js")))
