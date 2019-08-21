# -*- coding: utf-8 -*-

"""Parsers for HIPPIE."""

import pandas as pd

import bio2bel_hgnc
from bio2bel.downloading import make_df_getter, make_downloader
from bio2bel_hippie.constants import HEADER, PATH, URL
from bio2bel_uniprot import get_slim_mappings_df

__all__ = [
    'download_database',
    'get_df',
]

download_database = make_downloader(URL, PATH)

get_df = make_df_getter(
    URL,
    PATH,
    sep='\t',
    names=HEADER,
    dtype={
        'source_entrez_id': str,
        'target_entrez_id': str,
    },
)


def get_df_preprocessed(url=None, uniprot_url=None) -> pd.DataFrame:
    """Get the HIPPIE dataframe enriched with UniProt and HGNC information."""
    up_mappings_df = get_slim_mappings_df(url=uniprot_url)
    up_entry_name_to_id = {}
    up_entry_name_to_tax_id = {}
    # acc and id are wrong for uniprot. don't worry.
    for up_id, up_entry_name, tax_id in up_mappings_df[['UniProtKB-ID', 'UniProtKB-AC', 'NCBI-Taxon']].values:
        up_entry_name_to_id[up_entry_name] = up_id
        up_entry_name_to_tax_id[up_entry_name] = str(tax_id)

    df = get_df(url=url)

    df['source_uniprot_id'] = df['source_uniprot_entry_name'].map(up_entry_name_to_id.get)
    df['source_tax_id'] = df['source_uniprot_entry_name'].map(up_entry_name_to_tax_id.get)

    df['target_uniprot_id'] = df['target_uniprot_entry_name'].map(up_entry_name_to_id.get)
    df['target_tax_id'] = df['target_uniprot_entry_name'].map(up_entry_name_to_tax_id.get)

    hgnc_manager = bio2bel_hgnc.Manager()
    if not hgnc_manager.is_populated():
        raise RuntimeError('HGNC is not populated. Run `bio2bel_hgnc populate`.')

    entrez_id_to_hgnc_symbol = hgnc_manager.build_entrez_id_to_hgnc_symbol_mapping()
    entrez_id_to_hgnc_id = hgnc_manager.build_entrez_id_to_hgnc_id_mapping()

    df['source_hgnc_id'] = df['source_entrez_id'].map(entrez_id_to_hgnc_symbol.get)
    df['source_hgnc_symbol'] = df['source_entrez_id'].map(entrez_id_to_hgnc_id.get)
    df['target_hgnc_id'] = df['target_entrez_id'].map(entrez_id_to_hgnc_symbol.get)
    df['target_hgnc_symbol'] = df['target_entrez_id'].map(entrez_id_to_hgnc_id.get)

    return df
