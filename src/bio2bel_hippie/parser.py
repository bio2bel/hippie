# -*- coding: utf-8 -*-

"""Parsers for HIPPIE."""

from bio2bel.downloading import make_df_getter, make_downloader
from .constants import HEADER, PATH, URL

__all__ = [
    'download_database',
    'get_df',
]

download_database = make_downloader(URL, PATH)

get_df = make_df_getter(URL, PATH, sep='\t', names=HEADER)
