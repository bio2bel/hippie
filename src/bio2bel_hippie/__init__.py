# -*- coding: utf-8 -*-

"""Bio2BEL HIPPIE."""

import logging
import os
from typing import Mapping, Optional

import itertools as itt
import time
from sqlalchemy import Column, Float, ForeignKey, Integer, String
from sqlalchemy.ext.declarative import DeclarativeMeta, declarative_base
from sqlalchemy.orm import backref, relationship
from tqdm import tqdm

import pybel.dsl
from bio2bel import AbstractManager, get_data_dir
from bio2bel.downloading import make_df_getter, make_downloader
from bio2bel.manager.bel_manager import BELManagerMixin
from bio2bel.manager.flask_manager import FlaskMixin
from pybel import BELGraph

logger = logging.getLogger(__name__)

MODULE = 'hippie'
PROTEIN_TABLE_NAME = f'{MODULE}_protein'
INTERACTION_TABLE_NAME = f'{MODULE}_interaction'

DATA_DIR = get_data_dir(MODULE)
URL = 'http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/hippie_current.txt'
PATH = os.path.join(DATA_DIR, 'hippie_current.txt')

HEADER = [
    'source_uniprot_id',
    'source_entrez_id',
    'target_uniprot_id',
    'target_entrez_id',
    'confidence',
    'metadata',
]

# Downloaders
download_database = make_downloader(URL, PATH)
get_df = make_df_getter(URL, PATH, sep='\t', names=HEADER)

# SQLAlchemy stuff
Base: DeclarativeMeta = declarative_base()


class Protein(Base):
    """Represents a protein."""

    __tablename__ = PROTEIN_TABLE_NAME
    id = Column(Integer, primary_key=True)

    entrez_id = Column(String, nullable=False, index=True, unique=True)
    uniprot_id = Column(String, nullable=True)

    def as_pybel(self) -> pybel.dsl.protein:
        """Serializes this protein as a protein."""
        return pybel.dsl.protein(
            namespace='ncbigene',
            name=self.entrez_id,
        )

    def __repr__(self):
        return f'<Protein entrez_id={self.entrez_id}, uniprot_id={self.uniprot_id}>'


class Interaction(Base):
    """Represents a protein-protein interaction."""

    __tablename__ = INTERACTION_TABLE_NAME
    id = Column(Integer, primary_key=True)

    source_id = Column(Integer, ForeignKey(f'{Protein.__tablename__}.id'), nullable=False)
    source = relationship(Protein, foreign_keys=[source_id], backref=backref('out_edges', lazy='dynamic'))

    target_id = Column(Integer, ForeignKey(f'{Protein.__tablename__}.id'), nullable=False)
    target = relationship(Protein, foreign_keys=[target_id], backref=backref('in_edges', lazy='dynamic'))

    confidence = Column(Float)

    # TODO parse experiments
    # experiments = ...

    # TODO parse sources
    # sources = ...

    def add_to_bel_graph(self, graph: BELGraph) -> None:
        """Add this interaction to a BEL graph as a complex abundance."""
        node = pybel.dsl.ComplexAbundance([
            self.source.as_pybel(),
            self.target.as_pybel(),
        ])
        graph.add_node_from_data(node)


class Manager(AbstractManager, BELManagerMixin, FlaskMixin):
    """Manager for HIPPIE."""

    module_name = MODULE
    _base = Base
    flask_admin_models = [Protein, Interaction]

    def is_populated(self) -> bool:
        return 0 < self.count_proteins()

    def count_proteins(self) -> int:
        return self._count_model(Protein)

    def count_interactions(self) -> int:
        return self._count_model(Interaction)

    def summarize(self) -> Mapping[str, int]:
        return dict(
            proteins=self.count_proteins(),
            interactions=self.count_interactions(),
        )

    def populate(self, url: Optional[str] = None):
        """Populate the database."""
        df = get_df(url=url)

        entrez_protein = {
            protein.entrez_id: Protein
            for protein in self.session.query(Protein)
        }

        i = itt.chain(
            df[['source_uniprot_id', 'source_entrez_id']].iterrows(),
            df[['target_uniprot_id', 'target_entrez_id']].iterrows()
        )

        for idx, (uniprot_id, entrez_id) in tqdm(i, total=(2 * len(df.index)), desc='proteins'):
            protein = entrez_protein.get(entrez_id)
            if protein is None:
                entrez_protein[entrez_id] = Protein(
                    entrez_id=entrez_id,
                    uniprot_id=uniprot_id,
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

    def to_bel(self) -> BELGraph():
        """Convert to a BEL graph."""
        bel_graph = BELGraph(
            name='HIPPIE',
            version='2.1',
        )

        # TODO should rely on entrez package to deal with semantics

        for interaction in tqdm(self._get_query(Interaction), total=self.count_interactions(), desc='interactions'):
            interaction.add_to_bel_graph(bel_graph)

        return bel_graph


main = Manager.get_cli()

if __name__ == '__main__':
    main()
