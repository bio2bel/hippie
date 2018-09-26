# -*- coding: utf-8 -*-

"""Test cases for Bio2BEL HIPPIE."""

import os

from bio2bel.testing import AbstractTemporaryCacheClassMixin
from bio2bel_hippie import Manager

HERE = os.path.abspath(os.path.dirname(__file__))
TEST_HIPPIE_URL = os.path.join(HERE, 'hippie_test.txt')
TEST_UNIPROT_URL = os.path.join(HERE, 'uniprot_test.txt')


class TemporaryCacheClassMixin(AbstractTemporaryCacheClassMixin):
    """A temporary cache that contains HIPPIE."""

    Manager = Manager

    @classmethod
    def populate(cls):
        """Populate the test HIPPIE database."""
        cls.manager.populate(url=TEST_HIPPIE_URL)
