# -*- coding: utf-8 -*-

"""Constants for Bio2BEL HIPPIE."""

import os

from bio2bel import get_data_dir

MODULE = 'hippie'
DATA_DIR = get_data_dir(MODULE)

URL = 'http://cbdm-01.zdv.uni-mainz.de/~mschaefer/hippie/hippie_current.txt'
PATH = os.path.join(DATA_DIR, 'hippie_current.txt')
HEADER = [
    'source_uniprot_entry_name',
    'source_entrez_id',
    'target_uniprot_entry_name',
    'target_entrez_id',
    'confidence',
    'metadata',
]
