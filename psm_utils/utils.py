"""Various utility functions."""

from __future__ import annotations

from pyteomics.mass import nist_mass  # type: ignore[import-untyped]


def mass_to_mz(mass: float, charge: int, adduct_mass: float | None = None) -> float:
    """
    Convert mass to m/z.

    Parameters
    ----------
    mass
        Mass of the uncharged ion without adducts.
    charge
        Charge of the ion.
    adduct_mass
        Mass of the charge-carrying adduct. Defaults to the mass of a proton.

    """
    if adduct_mass is None:
        _adduct_mass = nist_mass["H"][1][0]
    else:
        _adduct_mass = float(adduct_mass)
    return (mass + charge * _adduct_mass) / charge


def mz_to_mass(mz: float, charge: int, adduct_mass: float | None = None) -> float:
    """
    Convert m/z to mass.

    Parameters
    ----------
    mz
        m/z of the charged ion and adducts.
    charge
        Charge of the ion.
    adduct_mass
        Mass of the charge-carrying adduct. Defaults to the mass of a proton.

    """
    if adduct_mass is None:
        _adduct_mass = nist_mass["H"][1][0]
    else:
        _adduct_mass = float(adduct_mass)
    return mz * charge - charge * _adduct_mass
