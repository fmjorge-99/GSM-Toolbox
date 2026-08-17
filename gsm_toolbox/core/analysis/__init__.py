"""Analysis families for GSM ToolBox.

Each module wraps a group of constraint-based methods and returns clean result
objects (rather than leaking solver internals) so the GUI can render them
uniformly. Phase 1 ships :mod:`fba`; later phases add strain design, omics
integration, gap-filling and QC.
"""
