import os
import shutil
import tempfile
import re
import requests
from typing import Optional, List, Dict

import numpy

# noinspection PyPackageRequirements
from Bio.PDB import PDBList
from furret.chains import ChainGroups
from furret.utilities import Citation


class Structure:
    def __init__(self, uniprot: str,
                 sequence: Optional[str] = None,
                 the_chains: Optional[ChainGroups] = None) -> None:
        the_string = uniprot.lower()
        if the_string and not the_string.isalnum():
            raise ValueError('Uniprot code should consist only of letters and digits')
        self.uniprot: str = the_string
        the_string = sequence
        if the_string:
            the_string = the_string.lower()
            if not the_string.isalpha():
                raise ValueError('sequence should consist only of letters')
        self.sequence: Optional[str] = the_string
        self.chains: Optional[ChainGroups] = the_chains
        self.method: Optional[str] = None
        self.downloaded: bool = False
        self.resolution: float = numpy.nan
        self.citations: List[Citation] = []

    def set_chains(self, the_chain_group: ChainGroups) -> None:
        self.chains: ChainGroups = the_chain_group

    def add_citation(self, citation: Citation) -> None:
        self.citations.append(citation)


class PDB(Structure):
    pdb_lister = PDBList()

    # fucking stupid people think that you may need an empty random directory to pollute your hierarchy
    shutil.rmtree(pdb_lister.obsolete_pdb, ignore_errors=True)
    # prevent a bunch of unwandted prints
    pdb_lister._verbose = False

    tmp_dir = os.path.join(tempfile.gettempdir(), '.{}'.format(hash(os.times())))
    os.makedirs(tmp_dir)
    if not os.path.exists(tmp_dir):
        raise PermissionError(f'''Unable to create {tmp_dir} temporary directory''')

    def __init__(self, uniprot: str,
                 sequence: Optional[str] = None,
                 the_chains: Optional[ChainGroups] = None,
                 code: Optional[str] = None) -> None:
        super().__init__(uniprot, sequence=sequence, the_chains=the_chains)
        the_string = code
        the_string = the_string.lower()
        if the_string and not the_string.isalnum():
            raise ValueError('PDB code should consist only of letters and digits')
        self.code: str = the_string
        if self.chains:
            self.coverage = 100 * len(self.chains.merged()) / len(self.sequence)
        else:
            self.coverage = 0.0

    # @timeout_retries(60, 5)
    def retrieve_file(self, directory: Optional[str] = None) -> bool:
        if not self.code:
            return False
        if not directory:
            directory = os.getcwd()
        else:
            os.makedirs(directory, exist_ok=True)
        if not os.access(directory, os.W_OK):
            raise PermissionError(f'''Can't write in {directory}''')
        file_retrieved = os.path.join(PDB.tmp_dir, 'pdb' + self.code + '.ent')
        file_name = os.path.join(directory, self.code + '.pdb')
        attempt = 0
        while attempt < 10:
            # noinspection PyBroadException
            try:
                PDB.pdb_lister.retrieve_pdb_file(self.code, file_format='pdb', pdir=PDB.tmp_dir)
                break
            except:
                attempt += 1
                if attempt == 10:
                    raise ConnectionAbortedError
        if os.path.exists(file_retrieved):
            shutil.move(file_retrieved, file_name)
            self.downloaded = True
        else:
            return False


# class PDBsm(PDB):
#     suffix = '_sm'
#
#     def __init__(self, code='', uniprot='',sequence='', ):
#         super(PDBsm, self).__init__(code, uniprot, sequence)


class Model(Structure):  # single structure within a uniprot accession
    def __init__(self, uniprot: str, sequence: Optional[str] = None,
                 data: Optional[Dict[str, any]] = None):
        assert data, "Need data to initialize Model"
        super().__init__(uniprot, sequence=sequence)

        self.similarity: float = data['similarity'] if 'similarty' in data else 0.0
        self.oligo_state: str = data['oligo-state'] if 'oligo-state' in data else ''
        self.coverage: float = data['coverage'] if 'coverage' in data else 0.0
        self.alignment: str = data['alignment'] if 'alignment' in data else ''
        # the following corrects wrong coding of newlines in uniprot
        self.alignment = re.sub(r'\\n', '\n', self.alignment) if self.alignment else ''
        self.gmqe: float = data['gmqe'] if 'gmqe' in data else 0.0
        self.qmean_norm: float = data['qmean_norm'] if 'qmean_norm' in data else 0.0
        self.request: str = data['coordinates'] if 'coordinates' in data else ''
        self.identity: str = data['identity'] if 'identity' in data else 0.0
        self.template: str = data['template'] if 'template' in data else ''
        self.qmean: float = data['qmean'] if 'qmean' in data else 0.0
        from_res = data['from'] if 'from' in data else None
        to_res = data['to'] if 'to' in data else None
        if from_res and to_res:
            chain = self.template.split('.')
            chain = chain[2] if len(chain) == 3 else 'SM'
            self.set_chains(ChainGroups([f'{chain}={from_res}-{to_res}']))

    def retrieve_file(self, directory=None):  #
        if not self.uniprot:
            return None
        if not directory:
            directory = os.getcwd()
        else:
            os.makedirs(directory, exist_ok=True)
        if not os.access(directory, os.W_OK):
            raise PermissionError(f'''Can't write in {directory}''')
        template = '_' + self.template if self.template else ''

        attempt = 0
        response = None
        while attempt < 10:
            # noinspection PyBroadException
            try:
                response = requests.get(self.request)
                if response:
                    break
            except:
                attempt += 1
                if attempt == 10:
                    raise ConnectionAbortedError
        file_name = os.path.join(directory, self.uniprot + template + '.pdb')
        with open(file_name, 'w') as the_model:
            the_model.write(response.text)

        self.downloaded = True
