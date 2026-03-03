"""Tests for PSF quality metrics in utilsTestImagingPlots.

These tests verify that measure_kappa, measure_skirt_level, and
measure_epsilon return the expected values when the PSF radial
profile is a perfect Gaussian matching the beam parameters.

For a perfect Gaussian:
    - kappa should be 0 (no deviation from Gaussian within the FWHM)
    - skirt_level should be exp(-4 log 2) ~ 0.0625
    - epsilon should be 1 (clean beam = dirty beam)
"""

import unittest
import numpy as np
from phangsPipeline.utilsTestImagingPlots import (
    measure_kappa,
    measure_skirt_level,
    measure_epsilon,
    fwhm_factor,
)


def _make_gaussian_profile(bmaj, bmin, n_bins=500, max_radius_factor=3.0):
    """Create a perfect Gaussian radial profile matching beam parameters.

    Parameters
    ----------
    bmaj : float
        Beam major axis FWHM in arcseconds.
    bmin : float
        Beam minor axis FWHM in arcseconds.
    n_bins : int
        Number of radial bins.
    max_radius_factor : float
        Profile extends to this factor times the geometric-mean FWHM.

    Returns
    -------
    radii : numpy.ndarray
        Radial bin centres in arcseconds.
    profile : numpy.ndarray
        Gaussian profile values.
    """
    fwhm = np.sqrt(bmaj * bmin)
    sigma = fwhm / fwhm_factor
    max_radius = max_radius_factor * fwhm
    radii = np.linspace(0, max_radius, n_bins)
    profile = np.exp(-radii**2 / (2.0 * sigma**2))
    return radii, profile


class TestMeasureKappa(unittest.TestCase):
    """Tests for measure_kappa."""

    def test_gaussian_kappa_is_zero(self):
        """A perfect Gaussian PSF should have kappa = 0."""
        bmaj, bmin = 5.0, 3.0
        radii, profile = _make_gaussian_profile(bmaj, bmin)
        kappa = measure_kappa(radii, profile, bmaj, bmin)
        self.assertIsNotNone(kappa)
        self.assertAlmostEqual(kappa, 0.0, places=4)

    def test_gaussian_kappa_symmetric_beam(self):
        """Circular beam should also give kappa = 0."""
        bmaj, bmin = 4.0, 4.0
        radii, profile = _make_gaussian_profile(bmaj, bmin)
        kappa = measure_kappa(radii, profile, bmaj, bmin)
        self.assertIsNotNone(kappa)
        self.assertAlmostEqual(kappa, 0.0, places=4)

    def test_peaked_psf_negative_kappa(self):
        """A PSF more peaked than the Gaussian should give kappa < 0."""
        bmaj, bmin = 5.0, 3.0
        radii, profile = _make_gaussian_profile(bmaj, bmin)
        # Make core more peaked by raising to a power > 1
        peaked = profile ** 1.5
        kappa = measure_kappa(radii, peaked, bmaj, bmin)
        self.assertIsNotNone(kappa)
        self.assertLess(kappa, 0.0)

    def test_flat_psf_positive_kappa(self):
        """A PSF flatter than the Gaussian should give kappa > 0."""
        bmaj, bmin = 5.0, 3.0
        radii, profile = _make_gaussian_profile(bmaj, bmin)
        # Make core flatter by raising to a power < 1
        flat = profile ** 0.5
        kappa = measure_kappa(radii, flat, bmaj, bmin)
        self.assertIsNotNone(kappa)
        self.assertGreater(kappa, 0.0)

    def test_empty_input_returns_none(self):
        """Empty arrays should return None."""
        kappa = measure_kappa(np.array([]), np.array([]), 5.0, 3.0)
        self.assertIsNone(kappa)

    def test_all_nan_returns_none(self):
        """All-NaN profile should return None."""
        radii = np.linspace(0, 10, 100)
        profile = np.full_like(radii, np.nan)
        kappa = measure_kappa(radii, profile, 5.0, 3.0)
        self.assertIsNone(kappa)


class TestMeasureSkirtLevel(unittest.TestCase):
    """Tests for measure_skirt_level."""

    def test_gaussian_skirt_level(self):
        """A perfect Gaussian should give skirt_level = exp(-4 log 2) ~ 0.0625."""
        expected = np.exp(-4.0 * np.log(2.0))
        bmaj, bmin = 5.0, 3.0
        radii, profile = _make_gaussian_profile(bmaj, bmin)
        skirt = measure_skirt_level(radii, profile, bmaj, bmin)
        self.assertIsNotNone(skirt)
        self.assertAlmostEqual(skirt, expected, places=3)

    def test_gaussian_skirt_level_circular(self):
        """Circular beam should also give skirt ~ 0.0625."""
        expected = np.exp(-4.0 * np.log(2.0))
        bmaj, bmin = 6.0, 6.0
        radii, profile = _make_gaussian_profile(bmaj, bmin)
        skirt = measure_skirt_level(radii, profile, bmaj, bmin)
        self.assertIsNotNone(skirt)
        self.assertAlmostEqual(skirt, expected, places=3)

    def test_broader_wings_higher_skirt(self):
        """PSF with broader wings should have skirt > 0.0625."""
        expected = np.exp(-4.0 * np.log(2.0))
        bmaj, bmin = 5.0, 3.0
        radii, profile = _make_gaussian_profile(bmaj, bmin)
        # Broaden the wings: raise to power < 1
        broad = profile ** 0.7
        skirt = measure_skirt_level(radii, broad, bmaj, bmin)
        self.assertIsNotNone(skirt)
        self.assertGreater(skirt, expected)

    def test_empty_input_returns_none(self):
        """Empty arrays should return None."""
        skirt = measure_skirt_level(np.array([]), np.array([]), 5.0, 3.0)
        self.assertIsNone(skirt)


class TestMeasureEpsilon(unittest.TestCase):
    """Tests for measure_epsilon."""

    def test_gaussian_epsilon_is_one(self):
        """A perfect Gaussian PSF should give epsilon = 1."""
        bmaj, bmin = 5.0, 3.0
        radii, profile = _make_gaussian_profile(bmaj, bmin)
        eps = measure_epsilon(radii, profile, bmaj, bmin)
        self.assertIsNotNone(eps)
        self.assertAlmostEqual(eps, 1.0, places=3)

    def test_gaussian_epsilon_circular(self):
        """Circular beam should also give epsilon = 1."""
        bmaj, bmin = 4.0, 4.0
        radii, profile = _make_gaussian_profile(bmaj, bmin)
        eps = measure_epsilon(radii, profile, bmaj, bmin)
        self.assertIsNotNone(eps)
        self.assertAlmostEqual(eps, 1.0, places=3)

    def test_excess_sidelobe_epsilon_below_one(self):
        """PSF with excess flux (wider wings) should give epsilon < 1."""
        bmaj, bmin = 5.0, 3.0
        radii, profile = _make_gaussian_profile(bmaj, bmin)
        # Add extra flux in the wings
        excess = profile + 0.1 * np.exp(-radii / 5.0)
        excess /= excess[0]  # re-normalise peak to 1
        eps = measure_epsilon(radii, excess, bmaj, bmin)
        self.assertIsNotNone(eps)
        self.assertLess(eps, 1.0)

    def test_empty_input_returns_none(self):
        """Empty arrays should return None."""
        eps = measure_epsilon(np.array([]), np.array([]), 5.0, 3.0)
        self.assertIsNone(eps)

    def test_profile_with_null(self):
        """Profile that crosses zero should still compute epsilon."""
        bmaj, bmin = 5.0, 3.0
        fwhm = np.sqrt(bmaj * bmin)
        sigma = fwhm / fwhm_factor
        radii = np.linspace(0, 3.0 * fwhm, 500)
        # Gaussian core with a negative trough beyond the FWHM
        profile = np.exp(-radii**2 / (2.0 * sigma**2))
        profile[radii > 1.5 * fwhm] = -0.05
        eps = measure_epsilon(radii, profile, bmaj, bmin)
        self.assertIsNotNone(eps)
        # Should be close to 1 since the Gaussian core is identical
        self.assertAlmostEqual(eps, 1.0, delta=0.05)


if __name__ == '__main__':
    unittest.main()
