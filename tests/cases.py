# -*- coding: utf-8 -*-

r"""Test cases for Bio2BEL HIPPIE.

To generate the uniprot test data, run:

>>> import pandas as pd
>>> import bio2bel_uniprot
>>> hdf = pd.read_csv('hippie_test.txt', sep='\t', header=None)
>>> proteins= set(hdf[0]) | set(hdf[2])
>>> udf = bio2bel_uniprot.get_mappings_df()
>>> udf[udf[1].isin(proteins)].to_csv('uniprot_test.tsv', sep='\t', index=None, header=None)
"""

import os

from bio2bel.testing import AbstractTemporaryCacheClassMixin
from bio2bel_hippie import Manager

HERE = os.path.abspath(os.path.dirname(__file__))
TEST_HIPPIE_URL = os.path.join(HERE, 'hippie_test.txt')
TEST_UNIPROT_URL = os.path.join(HERE, 'uniprot_test.tsv')


class TemporaryCacheClassMixin(AbstractTemporaryCacheClassMixin):
    """A temporary cache that contains HIPPIE."""

    Manager = Manager

    @classmethod
    def populate(cls):
        """Populate the test HIPPIE database."""
        cls.manager.populate(url=TEST_HIPPIE_URL, uniprot_url=TEST_UNIPROT_URL)
