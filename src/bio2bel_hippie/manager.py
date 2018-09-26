# -*- coding: utf-8 -*-

"""Bio2BEL HIPPIE."""

import itertools as itt
import logging
import time
from typing import Mapping, Optional

from tqdm import tqdm

from bio2bel import AbstractManager
from bio2bel.manager.bel_manager import BELManagerMixin
from bio2bel.manager.flask_manager import FlaskMixin
from bio2bel_uniprot import get_slim_mappings_df
from pybel import BELGraph
from .constants import MODULE
from .models import Base, Interaction, Protein
from .parser import get_df

__all__ = [
    'Manager',
]

logger = logging.getLogger(__name__)


class Manager(AbstractManager, BELManagerMixin, FlaskMixin):
    """Manager for HIPPIE."""

    module_name = MODULE
    _base = Base
    flask_admin_models = [Protein, Interaction]

    def is_populated(self) -> bool:
        """Check if the database is populated."""
        return 0 < self.count_proteins()

    def count_proteins(self) -> int:
        """Count the number of proteins in the database."""
        return self._count_model(Protein)

    def count_interactions(self) -> int:
        """Count the number of interactions in the database."""
        return self._count_model(Interaction)

    def summarize(self) -> Mapping[str, int]:
        """Summarize the contents of the database."""
        return dict(
            proteins=self.count_proteins(),
            interactions=self.count_interactions(),
        )

    def populate(self, url: Optional[str] = None, uniprot_url: Optional[str] = None) -> None:
        """Populate the database."""
        print('get uniprot mappings')
        up_mappings_df = get_slim_mappings_df(url=uniprot_url)
        up_id_to_up_acc = {}
        up_id_to_tax_id = {}
        for idx, (up_acc, up_id, tax_id) in up_mappings_df[['UniProtKB-AC', 'UniProtKB-ID', 'NCBI-Taxon']].iterrows():
            up_id_to_up_acc[up_id] = up_acc
            up_id_to_tax_id[up_id] = tax_id

        print('get hippie')
        df = get_df(url=url)

        entrez_protein = {
            protein.entrez_id: Protein
            for protein in self.session.query(Protein)
        }

        i = itt.chain(
            df[['source_uniprot_id', 'source_entrez_id']].iterrows(),
            df[['target_uniprot_id', 'target_entrez_id']].iterrows(),
        )

        for idx, (uniprot_id, entrez_id) in tqdm(i, total=(2 * len(df.index)), desc='proteins'):
            protein = entrez_protein.get(entrez_id)
            if protein is None:
                entrez_protein[entrez_id] = Protein(
                    entrez_id=entrez_id,
                    uniprot_id=uniprot_id,
                    uniprot_accession=up_id_to_up_acc.get(uniprot_id),
                    taxonomy_id=up_id_to_tax_id.get(uniprot_id),
                )

        logger.info('committing protein models')
        time_commit_proteins_start = time.time()
        self.session.add_all(list(entrez_protein.values()))
        self.session.commit()
        logger.info('committed protein models in %.2f seconds', time.time() - time_commit_proteins_start)

        _columns = ['source_entrez_id', 'target_entrez_id', 'confidence']
        for idx, (source_entrez_id, target_entrez_id, confidence) in tqdm(df[_columns].iterrows(), total=len(df.index)):
            interaction = Interaction(
                source=entrez_protein.get(source_entrez_id),
                target=entrez_protein.get(target_entrez_id),
                confidence=confidence,
            )
            self.session.add(interaction)

        logger.info('committing interaction models')
        time_commit_interactions_start = time.time()
        self.session.commit()
        logger.info('committed interaction models in %.2f seconds', time.time() - time_commit_interactions_start)

    def to_bel(self) -> BELGraph:
        """Convert to a BEL graph."""
        bel_graph = BELGraph(
            name='HIPPIE',
            version='2.1',
        )

        # TODO should rely on entrez package to deal with semantics

        for interaction in tqdm(self._get_query(Interaction), total=self.count_interactions(), desc='interactions'):
            interaction.add_to_bel_graph(bel_graph)

        return bel_graph
