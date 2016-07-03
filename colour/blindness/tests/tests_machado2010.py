# !/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Defines unit tests for :mod:`colour.blindness.machado2010` module.
"""

from __future__ import division, unicode_literals

import numpy as np
import unittest

from colour.blindness import (
    CVD_MATRICES_MACHADO_2010,
    cvd_matrix_Machado2010,
    anomalous_trichromacy_cmfs_Machado2010,
    anomalous_trichromacy_matrix_Machado2010)
from colour.characterisation import DISPLAYS_RGB_PRIMARIES
from colour.colorimetry import LMS_CMFS
from colour.utilities import ignore_numpy_errors

__author__ = 'Colour Developers'
__copyright__ = 'Copyright (C) 2013 - 2015 - Colour Developers'
__license__ = 'New BSD License - http://opensource.org/licenses/BSD-3-Clause'
__maintainer__ = 'Colour Developers'
__email__ = 'colour-science@googlegroups.com'
__status__ = 'Production'

__all__ = ['TestAnomalousTrichromacyCmfsMachado2010',
           'TestAnomalousTrichromacyMatrixMachado2010',
           'TestCvdMatrixMachado2010']


class TestAnomalousTrichromacyCmfsMachado2010(unittest.TestCase):
    """
    Defines :func:`colour.blindness.machado2010.\
anomalous_trichromacy_cmfs_Machado2010` definition unit tests methods.
    """

    def test_anomalous_trichromacy_cmfs_Machado2010(self):
        """
        Tests :func:`colour.blindness.machado2010.\
    anomalous_trichromacy_cmfs_Machado2010` definition.
        """

        cmfs = LMS_CMFS.get('Smith & Pokorny 1975 Normal Trichromats')
        np.testing.assert_almost_equal(
            anomalous_trichromacy_cmfs_Machado2010(
                cmfs,
                np.array([0, 0, 0]))[450],
            cmfs[450],
            decimal=7)

        np.testing.assert_almost_equal(
            anomalous_trichromacy_cmfs_Machado2010(
                cmfs,
                np.array([1, 0, 0]))[450],
            np.array([0.03631700, 0.06350000, 0.91000000]),
            decimal=7)

        np.testing.assert_almost_equal(
            anomalous_trichromacy_cmfs_Machado2010(
                cmfs,
                np.array([0, 1, 0]))[450],
            np.array([0.03430000, 0.06178404, 0.91000000]),
            decimal=7)

        np.testing.assert_almost_equal(
            anomalous_trichromacy_cmfs_Machado2010(
                cmfs,
                np.array([0, 0, 1]))[450],
            np.array([0.03430000, 0.06350000, 0.92270240]),
            decimal=7)

        np.testing.assert_almost_equal(
            anomalous_trichromacy_cmfs_Machado2010(
                cmfs,
                np.array([10, 0, 0]))[450],
            np.array([0.05447001, 0.06350000, 0.91000000]),
            decimal=7)

        np.testing.assert_almost_equal(
            anomalous_trichromacy_cmfs_Machado2010(
                cmfs,
                np.array([0, 10, 0]))[450],
            np.array([0.03430000, 0.04634036, 0.91000000]),
            decimal=7)

        np.testing.assert_almost_equal(
            anomalous_trichromacy_cmfs_Machado2010(
                cmfs,
                np.array([0, 0, 10]))[450],
            np.array([0.03430000, 0.06350000, 1.00000000]),
            decimal=7)


class TestAnomalousTrichromacyMatrixMachado2010(unittest.TestCase):
    """
    Defines :func:`colour.blindness.machado2010.\
anomalous_trichromacy_matrix_Machado2010` definition unit tests methods.
    """

    def test_anomalous_trichromacy_matrix_Machado2010(self):
        """
        Tests :func:`colour.blindness.machado2010.\
anomalous_trichromacy_matrix_Machado2010` definition.
        """

        cmfs = LMS_CMFS.get('Smith & Pokorny 1975 Normal Trichromats')
        primaries = DISPLAYS_RGB_PRIMARIES['Typical CRT Brainard 1997']
        np.testing.assert_almost_equal(
            anomalous_trichromacy_matrix_Machado2010(
                cmfs,
                primaries,
                np.array([0, 0, 0])),
            np.identity(3),
            decimal=7)

        np.testing.assert_allclose(
            anomalous_trichromacy_matrix_Machado2010(
                cmfs,
                primaries,
                np.array([2, 0, 0])),
            CVD_MATRICES_MACHADO_2010.get('Protanomaly').get(0.1),
            rtol=0.0001,
            atol=0.0001)

        np.testing.assert_allclose(
            anomalous_trichromacy_matrix_Machado2010(
                cmfs,
                primaries,
                np.array([20, 0, 0])),
            CVD_MATRICES_MACHADO_2010.get('Protanomaly').get(1.0),
            rtol=0.0001,
            atol=0.0001)

        np.testing.assert_allclose(
            anomalous_trichromacy_matrix_Machado2010(
                cmfs,
                primaries,
                np.array([0, 2, 0])),
            CVD_MATRICES_MACHADO_2010.get('Deuteranomaly').get(0.1),
            rtol=0.0001,
            atol=0.0001)

        np.testing.assert_allclose(
            anomalous_trichromacy_matrix_Machado2010(
                cmfs,
                primaries,
                np.array([0, 20, 0])),
            CVD_MATRICES_MACHADO_2010.get('Deuteranomaly').get(1.0),
            rtol=0.0001,
            atol=0.0001)

        # TODO: Check inconsistency with ground truth values, d_S domain seems
        # to be [5, 59] instead of [0, 20].
        np.testing.assert_allclose(
            anomalous_trichromacy_matrix_Machado2010(
                cmfs,
                primaries,
                np.array([0, 0, 5])),
            CVD_MATRICES_MACHADO_2010.get('Tritanomaly').get(0.1),
            rtol=0.0001,
            atol=0.0001)

        np.testing.assert_allclose(
            anomalous_trichromacy_matrix_Machado2010(
                cmfs,
                primaries,
                np.array([0, 0, 59])),
            CVD_MATRICES_MACHADO_2010.get('Tritanomaly').get(1.0),
            rtol=0.001,
            atol=0.001)


class TestCvdMatrixMachado2010(unittest.TestCase):
    """
    Defines :func:`colour.blindness.machado2010.cvd_matrix_Machado2010`
    definition unit tests methods.
    """

    def test_cvd_matrix_Machado2010(self):
        """
        Tests :func:`colour.blindness.machado2010.cvd_matrix_Machado2010`
        definition.
        """

        np.testing.assert_almost_equal(
            cvd_matrix_Machado2010('Protanomaly', 0.0),
            np.array([[1, 0, 0],
                      [0, 1, 0],
                      [0, 0, 1]]),
            decimal=7)

        np.testing.assert_almost_equal(
            cvd_matrix_Machado2010('Deuteranomaly', 0.1),
            np.array([[0.86643500, 0.17770400, -0.04413900],
                      [0.04956700, 0.93906300, 0.01137000],
                      [-0.00345300, 0.00723300, 0.99622000]]),
            decimal=7)

        np.testing.assert_almost_equal(
            cvd_matrix_Machado2010('Tritanomaly', 1.0),
            np.array([[1.25552800, -0.07674900, -0.17877900],
                      [-0.07841100, 0.93080900, 0.14760200],
                      [0.00473300, 0.69136700, 0.30390000]]),
            decimal=7)

        np.testing.assert_almost_equal(
            cvd_matrix_Machado2010('Tritanomaly', 0.55),
            np.array([[1.06088700, -0.01504350, -0.04584350],
                      [-0.01895750, 0.96774750, 0.05121150],
                      [0.00317700, 0.27513700, 0.72168600]]),
            decimal=7)

    @ignore_numpy_errors
    def test_nan_cvd_matrix_Machado2010(self):
        """
        Tests :func:`colour.blindness.machado2010.cvd_matrix_Machado2010`
        definition nan support.
        """

        for case in [-1.0, 0.0, 1.0, -np.inf, np.inf, np.nan]:
            cvd_matrix_Machado2010('Tritanomaly', case)


if __name__ == '__main__':
    unittest.main()
