#!/usr/bin/env python

"""
useful functions and classes
"""

from .condition import CONDITION
from .core import Page


def all_conditions() -> list[Page]:
    """Get all pages of specified conditions.

    Returns:
        list[Page]: list of Pages.
    """
    pages: list[Page] = []
    for simul_data in CONDITION.simul_data:
        for disp_type in CONDITION.disp_type:
            for frac_up in CONDITION.frac_up:
                page = Page(simul_data=simul_data, disp_type=disp_type, frac_up=frac_up)
                pages.append(page)
    return pages
