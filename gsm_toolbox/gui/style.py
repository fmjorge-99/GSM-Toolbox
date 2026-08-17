"""Centralized theme: color palette and a modern Qt stylesheet (QSS).

A clean, light, flat look with a single accent color — calm and readable for
long analysis sessions, which suits scientific desktop tools.
"""

from __future__ import annotations

# Palette ---------------------------------------------------------------
ACCENT = "#2D6CDF"        # primary actions / selection
ACCENT_DARK = "#1F52B0"
ACCENT_SOFT = "#E8F0FE"   # selection background
BG = "#F5F6F8"            # window background
SURFACE = "#FFFFFF"       # panels / cards
BORDER = "#E0E3E8"
TEXT = "#1F2733"
TEXT_MUTED = "#6B7280"
SUCCESS = "#1E8E3E"
DANGER = "#D93025"

# Network map colors
NODE_METABOLITE = "#2D6CDF"
NODE_ADDED_METABOLITE = "#34A853"   # green: metabolite added to the model (heterologous)
NODE_REACTION = "#8A93A0"
NODE_HIGHLIGHT = "#FBBC04"
EDGE = "#CBD2DA"
FLUX_FORWARD = "#E8453C"
FLUX_REVERSE = "#1A73E8"

#: Categorical series colours for plots, in fixed assignment order.
#:
#: Chosen by maximising the worst pairwise OKLab separation, not by taste. The palette
#: it replaces reused ACCENT (#2D6CDF) and FLUX_REVERSE (#1A73E8) as its first two
#: series — two blues ΔE 2.1 apart, indistinguishable on a line plot — and repeated
#: after seven series. This one is twelve long with a worst adjacent-pair ΔE of 28
#: (11.3 under simulated protanopia).
#:
#: Twelve mutually distinguishable colours is achievable for normal vision; twelve that
#: also survive red-green colour blindness is not, because protanopia and deuteranopia
#: collapse the red-green axis and leave roughly one usable dimension. Plots therefore
#: vary dash pattern and marker as well, so identity never rests on hue alone.
SERIES_COLORS = [
    "#eda100",   # amber
    "#4a3aa7",   # indigo
    "#8a8f00",   # olive
    "#8430ce",   # violet
    "#008300",   # green
    "#c026a0",   # magenta
    "#1baf7a",   # teal-green
    "#c5221f",   # red
    "#3aa8d8",   # sky
    "#d2691e",   # burnt orange
    "#0277bd",   # deep blue
    "#e87ba4",   # pink
]

#: Secondary channels, cycled against the colours so a repeat needs both to coincide.
#: Four dashes x twelve colours = 48 visually distinct series.
SERIES_DASHES = ["-", "--", "-.", ":"]
SERIES_MARKERS = ["o", "s", "^", "D", "v", "P"]


STYLESHEET = f"""
* {{
    font-family: "Segoe UI", "Inter", sans-serif;
    font-size: 13px;
    color: {TEXT};
}}

QMainWindow, QWidget {{
    background: {BG};
}}

/* ---- Docks ---- */
QDockWidget {{
    titlebar-close-icon: none;
    font-weight: 600;
}}
QDockWidget::title {{
    background: {SURFACE};
    padding: 8px 10px;
    border-bottom: 1px solid {BORDER};
}}

/* ---- Surfaces / group boxes ---- */
QGroupBox {{
    background: {SURFACE};
    border: 1px solid {BORDER};
    border-radius: 10px;
    margin-top: 14px;
    padding: 10px;
}}
QGroupBox::title {{
    subcontrol-origin: margin;
    left: 12px;
    padding: 0 4px;
    color: {TEXT_MUTED};
    font-weight: 700;
}}

/* ---- Tabs ---- */
QTabWidget::pane {{
    border: 1px solid {BORDER};
    border-radius: 10px;
    background: {SURFACE};
    top: -1px;
}}
QTabBar::tab {{
    background: transparent;
    color: {TEXT_MUTED};
    padding: 8px 16px;
    margin-right: 2px;
    border: none;
    border-top-left-radius: 8px;
    border-top-right-radius: 8px;
    font-weight: 600;
}}
QTabBar::tab:selected {{
    color: {ACCENT};
    background: {SURFACE};
    border-bottom: 2px solid {ACCENT};
}}
QTabBar::tab:hover:!selected {{
    color: {TEXT};
}}

/* ---- Tables ---- */
QTableView {{
    background: {SURFACE};
    alternate-background-color: #FAFBFC;
    gridline-color: {BORDER};
    border: 1px solid {BORDER};
    border-radius: 8px;
    selection-background-color: {ACCENT_SOFT};
    selection-color: {TEXT};
}}
QHeaderView::section {{
    background: #F0F2F5;
    color: {TEXT_MUTED};
    padding: 6px 8px;
    border: none;
    border-right: 1px solid {BORDER};
    border-bottom: 1px solid {BORDER};
    font-weight: 600;
}}
QTableView::item {{ padding: 2px 4px; }}

/* ---- Buttons ---- */
QPushButton {{
    background: {SURFACE};
    border: 1px solid {BORDER};
    border-radius: 8px;
    padding: 7px 14px;
    font-weight: 600;
}}
QPushButton:hover {{ border-color: {ACCENT}; color: {ACCENT}; }}
QPushButton:pressed {{ background: {ACCENT_SOFT}; }}
QPushButton:disabled {{ color: #B0B6BE; border-color: #ECEEF1; }}
QPushButton#primary {{
    background: {ACCENT};
    color: white;
    border: none;
}}
QPushButton#primary:hover {{ background: {ACCENT_DARK}; color: white; }}
QPushButton#primary:disabled {{ background: #AEC2EC; }}

/* ---- Inputs ---- */
QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox {{
    background: {SURFACE};
    border: 1px solid {BORDER};
    border-radius: 8px;
    padding: 6px 8px;
    selection-background-color: {ACCENT_SOFT};
}}
QLineEdit:focus, QComboBox:focus, QSpinBox:focus, QDoubleSpinBox:focus {{
    border-color: {ACCENT};
}}
QComboBox::drop-down {{
    subcontrol-origin: padding;
    subcontrol-position: top right;
    width: 22px;
    border-left: 1px solid {BORDER};
    background: #F0F2F5;
    border-top-right-radius: 8px;
    border-bottom-right-radius: 8px;
}}
QComboBox::drop-down:hover {{ background: {ACCENT_SOFT}; }}
QComboBox::down-arrow {{ image: url("__CHEVRON__"); width: 12px; height: 12px; }}

/* ---- Menus / toolbar ---- */
QMenuBar {{ background: {SURFACE}; border-bottom: 1px solid {BORDER}; }}
QMenuBar::item {{ padding: 7px 12px; background: transparent; }}
QMenuBar::item:selected {{ background: {ACCENT_SOFT}; color: {ACCENT}; }}
QMenu {{ background: {SURFACE}; border: 1px solid {BORDER}; border-radius: 8px; padding: 4px; }}
QMenu::item {{ padding: 6px 24px 6px 12px; border-radius: 6px; }}
QMenu::item:selected {{ background: {ACCENT_SOFT}; color: {ACCENT}; }}
QToolBar {{ background: {SURFACE}; border-bottom: 1px solid {BORDER}; padding: 4px; spacing: 4px; }}
QToolButton {{ padding: 6px 10px; border-radius: 8px; }}
QToolButton:hover {{ background: {ACCENT_SOFT}; }}

/* ---- Status bar ---- */
QStatusBar {{ background: {SURFACE}; border-top: 1px solid {BORDER}; color: {TEXT_MUTED}; }}

/* ---- Scrollbars ---- */
QScrollBar:vertical {{ background: transparent; width: 10px; margin: 2px; }}
QScrollBar::handle:vertical {{ background: #C4CAD2; border-radius: 5px; min-height: 30px; }}
QScrollBar::handle:vertical:hover {{ background: #A9B1BB; }}
QScrollBar:horizontal {{ background: transparent; height: 10px; margin: 2px; }}
QScrollBar::handle:horizontal {{ background: #C4CAD2; border-radius: 5px; min-width: 30px; }}
QScrollBar::add-line, QScrollBar::sub-line {{ width: 0; height: 0; }}

/* ---- Progress bar ---- */
QProgressBar {{ border: 1px solid {BORDER}; border-radius: 8px; background: {SURFACE}; text-align: center; }}
QProgressBar::chunk {{ background: {ACCENT}; border-radius: 7px; }}

QListWidget {{
    background: {SURFACE};
    border: 1px solid {BORDER};
    border-radius: 8px;
    padding: 4px;
}}
QListWidget::item {{ padding: 6px; border-radius: 6px; }}
QListWidget::item:selected {{ background: {ACCENT_SOFT}; color: {TEXT}; }}
"""


def _inject_icon_paths(qss: str) -> str:
    """Replace icon placeholders with absolute paths (QSS needs forward slashes)."""
    from .. import resources

    chevron = resources.resource_path("icons", "chevron-down.svg").replace("\\", "/")
    return qss.replace("__CHEVRON__", chevron)


STYLESHEET = _inject_icon_paths(STYLESHEET)
