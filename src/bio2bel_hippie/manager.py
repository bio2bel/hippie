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
from pybel import BELGraph
from .constants import MODULE
from .models import Base, Interaction, Protein
from .parser import get_df_preprocessed

__all__ = [
    'Manager',
]

logger = logging.getLogger(__name__)


class Manager(AbstractManager, BELManagerMixin, FlaskMixin):
    """Protein-protein physical interactions."""

    module_name = MODULE
    edge_model = Interaction
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
        logger.info('Getting HIPPIE data')
        df = get_df_preprocessed(url=url)

        protein_tuples = set(map(tuple, itt.chain(
            df[['source_uniprot_id', 'source_uniprot_entry_name', 'source_entrez_id', 'source_tax_id', 'source_hgnc_id',
                'source_hgnc_symbol']].values,
            df[['target_uniprot_id', 'target_uniprot_entry_name', 'target_entrez_id', 'target_tax_id', 'target_hgnc_id',
                'target_hgnc_symbol']].values,
        )))

        uniprot_id_to_protein = {}
        it = tqdm(protein_tuples, desc='HIPPIE: making protein models')
        for uniprot_id, uniprot_entry_name, entrez_id, tax_id, hgnc_id, symbol in it:
            uniprot_id_to_protein[uniprot_id] = Protein(
                entrez_id=entrez_id,
                uniprot_entry_name=uniprot_entry_name,
                uniprot_id=uniprot_id,
                taxonomy_id=tax_id,
                symbol=symbol,
                hgnc_id=hgnc_id,
            )
        logger.info('Made entrez to protein dict')

        logger.info('committing protein models')
        time_commit_proteins_start = time.time()
        self.session.add_all(list(uniprot_id_to_protein.values()))
        self.session.commit()
        logger.info('committed protein models in %.2f seconds', time.time() - time_commit_proteins_start)

        _columns = ['source_uniprot_id', 'target_uniprot_id', 'confidence']
        it = tqdm(df[_columns].values, total=len(df.index), desc='HIPPIE: making PPI models')
        for source_uniprot_id, target_uniprot_id, confidence in it:
            interaction = Interaction(
                source=uniprot_id_to_protein[source_uniprot_id],
                target=uniprot_id_to_protein[target_uniprot_id],
                confidence=confidence,
            )
            self.session.add(interaction)

        logger.info('committing interaction models')
        time_commit_interactions_start = time.time()
        self.session.commit()
        logger.info('committed interaction models in %.2f seconds', time.time() - time_commit_interactions_start)

    def to_bel(self, namespace: Optional[str] = None) -> BELGraph:
        """Convert to a BEL graph."""
        bel_graph = BELGraph(
            name='HIPPIE',
            version='2.1',
        )

        for interaction in tqdm(self._get_query(Interaction), total=self.count_interactions(), desc='interactions'):
            interaction.add_to_bel_graph(bel_graph, namespace=namespace)

        return bel_graph
