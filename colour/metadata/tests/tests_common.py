#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Defines unit tests for :mod:`colour.metadata.common` module.
"""

from __future__ import division, unicode_literals

import unittest

import colour
from colour.metadata.common import (
    parse_parameters_field_metadata,
    parse_returns_field_metadata,
    parse_notes_field_metadata,
    set_metadata,
    filter_metadata_registry)
from colour.metadata.common import NotesMetadata, ParameterMetadata

__author__ = 'Colour Developers'
__copyright__ = 'Copyright (C) 2013-2016 - Colour Developers'
__license__ = 'New BSD License - http://opensource.org/licenses/BSD-3-Clause'
__maintainer__ = 'Colour Developers'
__email__ = 'colour-science@googlegroups.com'
__status__ = 'Production'

__all__ = ['TestParseParametersFieldMetadata',
           'TestParseReturnsFieldMetadata',
           'TestParseNotesFieldMetadata',
           'TestSetMetadata',
           'TestFilterMetadataRegistry']


class TestParseParametersFieldMetadata(unittest.TestCase):
    """
    Defines :func:`colour.metadata.common.parse_parameters_field_metadata`
    definition units tests methods.
    """

    def test_parse_parameters_field_metadata(self):
        """
        Tests :func:`colour.metadata.common.parse_parameters_field_metadata`
        definition.
        """

        field = (['Lstar', u'numeric or array_like'],
                 "metadata : {'type': 'Lightness', 'symbol': 'L^\\star', \
'extent': 100} *Lightness* :math:`L^\\star`")
        self.assertTupleEqual(
            tuple(parse_parameters_field_metadata(field)),
            ('Lightness', 'L^\\star', 100))


class TestParseReturnsFieldMetadata(unittest.TestCase):
    """
    Defines :func:`colour.metadata.common.parse_returns_field_metadata`
    definition units tests methods.
    """

    def test_parse_returns_field_metadata(self):
        """
        Tests :func:`colour.metadata.common.parse_returns_field_metadata`
        definition.
        """

        field = (['Lstar', u'numeric or array_like'],
                 "metadata : {'type': 'Lightness', 'symbol': 'L^\\star', \
'extent': 100} *Lightness* :math:`L^\\star`")
        self.assertTupleEqual(
            tuple(parse_returns_field_metadata(field)),
            ('Lightness', 'L^\\star', 100))


class TestParseNotesFieldMetadata(unittest.TestCase):
    """
    Defines :func:`colour.metadata.common.parse_notes_field_metadata`
    definition units tests methods.
    """

    def test_parse_notes_field_metadata(self):
        """
        Tests :func:`colour.metadata.common.parse_notes_field_metadata`
        definition.
        """

        field = ("metadata : {'method_name': 'Wyszecki 1963', \
'method_strict_name':", "'Wyszecki (1963)'}")
        self.assertTupleEqual(
            tuple(parse_notes_field_metadata(field)),
            ('Wyszecki 1963', 'Wyszecki (1963)'))


class TestSetMetadata(unittest.TestCase):
    """
    Defines :func:`colour.metadata.common.set_metadata` definition units
    tests methods.
    """

    def test_set_metadata(self):
        """
        Tests :func:`colour.metadata.common.set_metadata` definition.
        """

        def fn_a(argument_1):
            """
            Summary of docstring.

            Description of docstring.

            Parameters
            ----------
            argument_1 : object
                metadata : {'type': 'type', 'symbol': 'symbol',
                    'extent': 'extent'}
                Description of `argument_1`.

            Returns
            -------
            object
                metadata : {'type': 'type', 'symbol': 'symbol',
                    'extent': 'extent'}
                Description of `object`.

            Notes
            -----
            metadata : {'method_name': 'method_name',
                'method_strict_name': 'method_strict_name'}
            """

            return argument_1

        set_metadata(fn_a)

        self.assertTrue(hasattr(fn_a, '__metadata__'))
        self.assertDictEqual(
            dict(fn_a.__metadata__),
            {'returns': [ParameterMetadata(
                type='type',
                symbol='symbol',
                extent='extent')],
                'notes': [NotesMetadata(
                    method_name='method_name',
                    method_strict_name='method_strict_name')],
                'parameters': [ParameterMetadata(
                    type='type',
                    symbol='symbol',
                    extent='extent')]})


class TestFilterMetadataRegistry(unittest.TestCase):
    """
    Defines :func:`colour.metadata.common.filter_metadata_registry`
    definition units tests methods.
    """

    def test_filter_metadata_registry(self):
        """
        Tests :func:`colour.metadata.common.filter_metadata_registry`
        definition.
        """

        self.assertSetEqual(
            set(filter_metadata_registry('Luminance',
                                         categories='parameters',
                                         attributes='type')),
            set([colour.lightness_Glasser1958,
                 colour.lightness_Wyszecki1963,
                 colour.lightness_CIE1976]))

        self.assertSetEqual(
            set(filter_metadata_registry('Luminance',
                                         categories='parameters',
                                         attributes='type',
                                         any_parameter=True)),
            set([colour.lightness_Glasser1958,
                 colour.lightness_Wyszecki1963,
                 colour.lightness_CIE1976,
                 colour.luminance_CIE1976]))

        self.assertSetEqual(
            set(filter_metadata_registry('Luminance',
                                         categories='returns',
                                         attributes='type')),
            set([colour.luminance_ASTMD153508,
                 colour.luminance_Newhall1943,
                 colour.luminance_CIE1976]))

        self.assertSetEqual(
            set(filter_metadata_registry('Luminance',
                                         categories=('parameters', 'returns'),
                                         attributes='type')),
            set([colour.lightness_Glasser1958,
                 colour.lightness_Wyszecki1963,
                 colour.lightness_CIE1976,
                 colour.luminance_ASTMD153508,
                 colour.luminance_Newhall1943,
                 colour.luminance_CIE1976]))

        self.assertSetEqual(
            set(filter_metadata_registry('Luminance',
                                         categories=('parameters', 'returns'),
                                         attributes='type',
                                         any_parameter=True)),
            set([colour.lightness_Glasser1958,
                 colour.lightness_Wyszecki1963,
                 colour.lightness_CIE1976,
                 colour.luminance_CIE1976,
                 colour.luminance_ASTMD153508,
                 colour.luminance_Newhall1943,
                 colour.luminance_CIE1976]))

        self.assertSetEqual(
            set(filter_metadata_registry('CIE 1976',
                                         categories='notes',
                                         attributes='method_name')),
            set([colour.Lab_to_LCHab,
                 colour.Lab_to_XYZ,
                 colour.LCHab_to_Lab,
                 colour.LCHuv_to_Luv,
                 colour.lightness_CIE1976,
                 colour.luminance_CIE1976,
                 colour.Luv_to_LCHuv,
                 colour.Luv_to_uv,
                 colour.Luv_to_XYZ,
                 colour.Luv_uv_to_xy,
                 colour.XYZ_to_Lab,
                 colour.XYZ_to_Luv]))


if __name__ == '__main__':
    unittest.main()
