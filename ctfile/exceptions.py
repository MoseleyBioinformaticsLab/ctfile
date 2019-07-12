"""
Custom exceptions.
"""


class SpecError(Exception):
    """Specification error."""


class IsotopeSpecError(SpecError):
    """Isotope specification error."""


class ChargeSpecError(SpecError):
    """Charge specification error."""
