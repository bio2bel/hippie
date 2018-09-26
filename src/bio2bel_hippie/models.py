# -*- coding: utf-8 -*-

"""SQLAlchemy models for Bio2BEL HIPPIE."""

import logging
from typing import Optional

from sqlalchemy import Column, Float, ForeignKey, Integer, String
from sqlalchemy.ext.declarative import DeclarativeMeta, declarative_base
from sqlalchemy.orm import backref, relationship

import pybel.dsl
from pybel import BELGraph
from .constants import MODULE

logger = logging.getLogger(__name__)

PROTEIN_TABLE_NAME = f'{MODULE}_protein'
INTERACTION_TABLE_NAME = f'{MODULE}_interaction'

# SQLAlchemy stuff
Base: DeclarativeMeta = declarative_base()


class Protein(Base):
    """Represents a protein."""

    __tablename__ = PROTEIN_TABLE_NAME
    id = Column(Integer, primary_key=True)

    entrez_id = Column(String, nullable=False, index=True, unique=True)
    uniprot_id = Column(String, nullable=True)
    uniprot_accession = Column(String, nullable=True)
    taxonomy_id = Column(String)

    def as_pybel(self) -> Optional[pybel.dsl.protein]:
        """Serialize as a PyBEL protein."""
        if not self.uniprot_id:
            logger.debug('no uniprot id for %d %s', self.id, self.entrez_id)
            return
        return pybel.dsl.protein(
            namespace='uniprot',
            name=self.uniprot_id,
        )

    def __repr__(self):  # noqa: D105
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

    def add_to_bel_graph(self, graph: BELGraph) -> bool:
        """Add this interaction to a BEL graph as a complex abundance."""
        source_dsl = self.source.as_pybel()
        target_dsl = self.target.as_pybel()
        if source_dsl is None or target_dsl is None:
            return False

        node = pybel.dsl.ComplexAbundance([source_dsl, target_dsl])
        graph.add_node_from_data(node)
        return True
