# -*- coding: utf-8 -*-

"""Constants for Bio2BEL HIPPIE."""

import os

from bio2bel import get_data_dir

MODULE = 'hippie'
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
