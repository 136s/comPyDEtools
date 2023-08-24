#!/usr/bin/env python

"""
useful functions and classes
"""

from .core import AllPages


def all_conditions() -> AllPages:
    """Get all pages of specified conditions.

    Returns:
        Pages: list of Pages.
    """
    return AllPages()
