# -*- coding: utf-8 -*-

"""Tests for population."""

from tests.cases import TemporaryCacheClassMixin


class TestPopulate(TemporaryCacheClassMixin):
    """Test population of the database."""

    def test_contents(self):
        """Test the contents."""
        self.assertEqual(14, self.manager.count_proteins())
        self.assertEqual(10, self.manager.count_interactions())
