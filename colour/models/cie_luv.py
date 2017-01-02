#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
CIE L*u*v* Colourspace
======================

Defines the *CIE L\*u\*v\** colourspace transformations:

-   :func:`XYZ_to_Luv`
-   :func:`Luv_to_XYZ`
-   :func:`Luv_to_uv`
-   :func:`Luv_uv_to_xy`
-   :func:`Luv_to_LCHuv`
-   :func:`LCHuv_to_Luv`

See Also
--------
`CIE Luv Colourspace Jupyter Notebook
<http://nbviewer.jupyter.org/github/colour-science/colour-notebooks/\
blob/master/notebooks/models/cie_luv.ipynb>`_

References
----------
.. [1]  CIE TC 1-48. (2004). CIE 1976 uniform colour spaces. In CIE 015:2004
        Colorimetry, 3rd Edition (p. 24). ISBN:978-3-901-90633-6
.. [2]  CIE TC 1-48. (2004). CIE 1976 uniform chromaticity scale diagram (UCS
        diagram). In CIE 015:2004 Colorimetry, 3rd Edition (p. 24).
        ISBN:978-3-901-90633-6
"""

from __future__ import division, unicode_literals

import numpy as np

from colour.algebra import cartesian_to_polar, polar_to_cartesian
from colour.colorimetry import ILLUMINANTS
from colour.constants import CIE_E, CIE_K
from colour.models import xy_to_xyY, xyY_to_XYZ
from colour.utilities import tsplit, tstack

__author__ = 'Colour Developers'
__copyright__ = 'Copyright (C) 2013-2016 - Colour Developers'
__license__ = 'New BSD License - http://opensource.org/licenses/BSD-3-Clause'
__maintainer__ = 'Colour Developers'
__email__ = 'colour-science@googlegroups.com'
__status__ = 'Production'

__all__ = ['XYZ_to_Luv',
           'Luv_to_XYZ',
           'Luv_to_uv',
           'Luv_uv_to_xy',
           'Luv_to_LCHuv',
           'LCHuv_to_Luv']


def XYZ_to_Luv(XYZ,
               illuminant=ILLUMINANTS.get(
                   'CIE 1931 2 Degree Standard Observer').get('D50')):
    """
    Converts from *CIE XYZ* tristimulus values to *CIE L\*u\*v\** colourspace.

    Parameters
    ----------
    XYZ : array_like
        metadata : {'type': 'CIE XYZ', 'symbol': 'XYZ', 'extent': (0, 1)}
        *CIE XYZ* tristimulus values.
    illuminant : array_like, optional
        metadata : {'type': 'CIE xyY', 'symbol': 'xyY', 'extent': (0, 1)}
        Reference *illuminant* *xy* chromaticity coordinates or *CIE xyY*
        colourspace array.

    Returns
    -------
    ndarray
        metadata : {'type': 'CIE Luv', 'symbol': 'L*u*v*', 'extent':
        ((0, 100), (-100, 100), (-100, 100))}
        *CIE L\*u\*v\** colourspace array.

    Notes
    -----
    metadata : {'method_name': 'CIE 1976', 'method_strict_name': 'CIE 1976'}

    Examples
    --------
    >>> XYZ = np.array([0.07049534, 0.10080000, 0.09558313])
    >>> XYZ_to_Luv(XYZ)  # doctest: +ELLIPSIS
    array([ 37.9856291..., -28.8021959...,  -1.3580070...])
    """

    X, Y, Z = tsplit(XYZ)
    X_r, Y_r, Z_r = tsplit(xyY_to_XYZ(xy_to_xyY(illuminant)))

    y_r = Y / Y_r

    L = np.where(y_r > CIE_E, 116 * y_r ** (1 / 3) - 16, CIE_K * y_r)

    u = (13 * L * ((4 * X / (X + 15 * Y + 3 * Z)) -
                   (4 * X_r / (X_r + 15 * Y_r + 3 * Z_r))))
    v = (13 * L * ((9 * Y / (X + 15 * Y + 3 * Z)) -
                   (9 * Y_r / (X_r + 15 * Y_r + 3 * Z_r))))

    Luv = tstack((L, u, v))

    return Luv


def Luv_to_XYZ(Luv,
               illuminant=ILLUMINANTS.get(
                   'CIE 1931 2 Degree Standard Observer').get('D50')):
    """
    Converts from *CIE L\*u\*v\** colourspace to *CIE XYZ* tristimulus values.

    Parameters
    ----------
    Luv : array_like
        metadata : {'type': 'CIE Luv', 'symbol': 'L*u*v*', 'extent':
        ((0, 100), (-100, 100), (-100, 100))}
        *CIE L\*u\*v\** colourspace array.
    illuminant : array_like, optional
        metadata : {'type': 'CIE xyY', 'symbol': 'xyY', 'extent': (0, 1)}
        Reference *illuminant* *xy* chromaticity coordinates or *CIE xyY*
        colourspace array.

    Returns
    -------
    ndarray
        metadata : {'type': 'CIE XYZ', 'symbol': 'XYZ', 'extent': (0, 1)}
        *CIE XYZ* tristimulus values.

    Notes
    -----
    metadata : {'method_name': 'CIE 1976', 'method_strict_name': 'CIE 1976'}

    Examples
    --------
    >>> Luv = np.array([37.9856291 , -28.80219593,  -1.35800706])
    >>> Luv_to_XYZ(Luv)  # doctest: +ELLIPSIS
    array([ 0.0704953...,  0.1008    ,  0.0955831...])
    """

    L, u, v = tsplit(Luv)
    X_r, Y_r, Z_r = tsplit(xyY_to_XYZ(xy_to_xyY(illuminant)))

    Y = np.where(L > CIE_E * CIE_K, ((L + 16) / 116) ** 3, L / CIE_K)

    a = 1 / 3 * ((52 * L / (u + 13 * L *
                            (4 * X_r / (X_r + 15 * Y_r + 3 * Z_r)))) - 1)
    b = -5 * Y
    c = -1 / 3.0
    d = Y * (39 * L / (v + 13 * L *
                       (9 * Y_r / (X_r + 15 * Y_r + 3 * Z_r))) - 5)

    X = (d - b) / (a - c)
    Z = X * a + b

    XYZ = tstack((X, Y, Z))

    return XYZ


def Luv_to_uv(Luv,
              illuminant=ILLUMINANTS.get(
                  'CIE 1931 2 Degree Standard Observer').get('D50')):
    """
    Returns the :math:`u\'v\'` chromaticity coordinates from given
    *CIE L\*u\*v\** colourspace array.

    Parameters
    ----------
    Luv : array_like
        metadata : {'type': 'CIE Luv', 'symbol': 'L*u*v*', 'extent':
        ((0, 100), (-100, 100), (-100, 100))}
        *CIE L\*u\*v\** colourspace array.
    illuminant : array_like, optional
        metadata : {'type': 'CIE xyY', 'symbol': 'xyY', 'extent': (0, 1)}
        Reference *illuminant* *xy* chromaticity coordinates or *CIE xyY*
        colourspace array.

    Returns
    -------
    ndarray
        metadata : {'type': "CIE u'v'", 'symbol': "u'v'", 'extent': (0, 1)}
        :math:`u\'v\'` chromaticity coordinates.

    Notes
    -----
    metadata : {'method_name': 'CIE 1976', 'method_strict_name': 'CIE 1976'}

    Examples
    --------
    >>> Luv = np.array([37.9856291 , -28.80219593,  -1.35800706])
    >>> Luv_to_uv(Luv)  # doctest: +ELLIPSIS
    array([ 0.1508531...,  0.4853297...])
    """

    X, Y, Z = tsplit(Luv_to_XYZ(Luv, illuminant))

    uv = tstack((4 * X / (X + 15 * Y + 3 * Z),
                 9 * Y / (X + 15 * Y + 3 * Z)))

    return uv


def Luv_uv_to_xy(uv):
    """
    Returns the *xy* chromaticity coordinates from given *CIE L\*u\*v\**
    colourspace :math:`u\'v\'` chromaticity coordinates.

    Parameters
    ----------
    uv : array_like
        metadata : {'type': "CIE u'v'", 'symbol': "u'v'", 'extent': (0, 1)}
        *CIE L\*u\*v\* u\'v\'* chromaticity coordinates.

    Returns
    -------
    ndarray
        metadata : {'type': 'CIE xy', 'symbol': 'xy', 'extent': (0, 1)}
        *xy* chromaticity coordinates.

    Notes
    -----
    metadata : {'method_name': 'CIE 1976', 'method_strict_name': 'CIE 1976'}

    References
    ----------
    .. [3]  Wikipedia. (n.d.). The reverse transformation. Retrieved from
            http://en.wikipedia.org/wiki/CIELUV#The_reverse_transformation

    Examples
    --------
    >>> uv = np.array([0.150853098829857, 0.485329708543180])
    >>> Luv_uv_to_xy(uv)  # doctest: +ELLIPSIS
    array([ 0.2641477...,  0.3777000...])
    """

    u, v = tsplit(uv)

    xy = tstack((9 * u / (6 * u - 16 * v + 12),
                 4 * v / (6 * u - 16 * v + 12)))

    return xy


def Luv_to_LCHuv(Luv):
    """
    Converts from *CIE L\*u\*v\** colourspace to *CIE LCH(uv)* colourspace.

    Parameters
    ----------
    Luv : array_like
        metadata : {'type': 'CIE Luv', 'symbol': 'L*u*v*', 'extent':
        ((0, 100), (-100, 100), (-100, 100))}
        *CIE L\*u\*v\** colourspace array.

    Returns
    -------
    ndarray
        metadata : {'type': 'CIE LCHuv', 'symbol': 'LCH(uv)', 'extent':
        ((0, 100), (0, 360), (0, 360))}
        *CIE LCH(uv)* colourspace array.

    Notes
    -----
    metadata : {'method_name': 'CIE 1976', 'method_strict_name': 'CIE 1976'}

    Examples
    --------
    >>> Luv = np.array([37.9856291 , -28.80219593,  -1.35800706])
    >>> Luv_to_LCHuv(Luv)  # doctest: +ELLIPSIS
    array([  37.9856291...,   28.8341927...,  182.6994640...])
    """

    L, u, v = tsplit(Luv)

    C, H = tsplit(cartesian_to_polar(tstack((u, v))))

    LCHuv = tstack((L, C, np.degrees(H) % 360))

    return LCHuv


def LCHuv_to_Luv(LCHuv):
    """
    Converts from *CIE LCH(uv)* colourspace to *CIE L\*u\*v\** colourspace.

    Parameters
    ----------
    LCHuv : array_like
        metadata : {'type': 'CIE LCHuv', 'symbol': 'LCH(uv)', 'extent':
        ((0, 100), (0, 360), (0, 360))}
        *CIE LCH(uv)* colourspace array.

    Returns
    -------
    ndarray
        metadata : {'type': 'CIE Luv', 'symbol': 'L*u*v*', 'extent':
        ((0, 100), (-100, 100), (-100, 100))}
        *CIE L\*u\*v\** colourspace array.

    Notes
    -----
    metadata : {'method_name': 'CIE 1976', 'method_strict_name': 'CIE 1976'}

    Examples
    --------
    >>> LCHuv = np.array([37.98562910, 28.83419279, 182.69946404])
    >>> LCHuv_to_Luv(LCHuv)  # doctest: +ELLIPSIS
    array([ 37.9856291..., -28.8021959...,  -1.3580070...])
    """

    L, C, H = tsplit(LCHuv)

    u, v = tsplit(polar_to_cartesian(tstack((C, np.radians(H)))))

    Luv = tstack((L, u, v))

    return Luv
