from datetime import datetime
import requests
import xmltodict

import pickle
import sys
import shutil
from multiprocessing import Pool
from typing import List
from requests import ReadTimeout, ConnectTimeout, HTTPError, Timeout, ConnectionError

from PyQt5.QtWidgets import QApplication, QStatusBar


import furret.config as config
from furret.utilities import format_filename, validate_string, seq2fasta, Obj
from furret.tables import *
from furret.meme import *
from typing import Dict


# noinspection DuplicatedCode
class Query:
    def __init__(self, the_query: str, status: QStatusBar) -> None:

        self.query = the_query
        self.time = datetime.now().isoformat()
        name = format_filename(self.time + '_' + the_query)
        self.querydir = os.path.join(config.working_directory, name)
        try:
            os.makedirs(self.querydir)
        except OSError:
            print(f"Query directory '{self.querydir}' exists!")
            print("This should't have happened")
            sys.exit(1)
        self.xmlfile = os.path.join(self.querydir, 'uniprot.xml')
        self.dump = os.path.join(self.querydir, 'query.pickle')
        self.tbldir = os.path.join(self.querydir, 'Tables')
        # os.makedirs(self.tbldir)
        self.tbldump = os.path.join(self.querydir, 'tables.pickle')
        self.seqdir = os.path.join(self.querydir, 'Sequences')
        os.makedirs(self.seqdir)
        self.structdir = os.path.join(self.querydir, 'Structures')
        # os.makedirs(self.structdir)
        self.prepdir = os.path.join(self.querydir, 'Prepared')
        self.imported = os.path.join(self.querydir, 'Imported')
        self.famdir = os.path.join(self.querydir, 'Families')
        # os.makedirs(self.famdir)
        self.famstrdir = os.path.join(self.querydir, 'Families_Structures')
        # os.makedirs(self.famstrdir)
        self.motivedir = os.path.join(self.querydir, 'Motives')
        # os.makedirs(self.motivedir)
        QApplication.processEvents()

        status.showMessage(f'Quering {self.query}, please be patient')
        QApplication.processEvents()
        response = requests.get(f'https://www.uniprot.org/uniprot/?query={self.query}&format=xml')
        text = response.text
        status.showMessage(f'Saving XML for {self.query}')
        QApplication.processEvents()
        with open(self.xmlfile, 'wt') as xmlfile:
            xmlfile.write(text)
        status.showMessage(f'Parsing XML')
        QApplication.processEvents()
        force_list = ('entry', 'accession', 'reference', 'dbReference', 'property', 'keyword', 'scope', 'name')
        data = xmltodict.parse(text, force_list=force_list)
        entries = data['uniprot']['entry']
        self.proteins: Dict[str, Protein] = {}
        for entry in entries:
            uniprot = entry['accession'][0]
            status.showMessage(f'Reading {uniprot}')
            QApplication.processEvents()
            self.proteins[uniprot] = Protein(entry)
            sequence_file = os.path.join(self.seqdir, uniprot+'.fasta')
            with open(sequence_file, encoding="ascii", mode='wt') as fasta:
                fasta.write(seq2fasta(self.proteins[uniprot].sequence, uniprot))
        self.tables = None
        self.process_tables(status)
        status.showMessage(f'Saving Query')
        self.save()
        status.showMessage(f'Done.')

    def save(self):
        pickle.dump(self, open(self.dump, 'wb'))
        with open(os.path.join(self.querydir, 'query.txt'), 'wt') as textfile:
            textfile.write(self.query + '\n')
            textfile.write(self.time + '\n')

    def process_tables(self, status):
        the_tables = Obj()
        os.makedirs(self.tbldir, exist_ok=True)
        status.showMessage(f'Generating keywords table')
        QApplication.processEvents()
        the_tables.keywords = process_keywords(self.proteins, self.tbldir)
        status.showMessage(f'Generating GO table')
        QApplication.processEvents()
        the_tables.go = process_go(self.proteins, self.tbldir)
        status.showMessage(f'Generating databases table')
        QApplication.processEvents()
        the_tables.db = process_links(self.proteins, self.tbldir)
        status.showMessage(f'Generating sequences table')
        QApplication.processEvents()
        the_tables.sequences = process_sequences(self.proteins, self.tbldir)
        status.showMessage(f'Generating PDB table')
        QApplication.processEvents()
        the_tables.pdb = process_pdb(self.proteins, self.tbldir)
        status.showMessage(f'Generating Swiss Models table')
        QApplication.processEvents()
        the_tables.sm = process_sm(self.proteins, self.tbldir)
        status.showMessage(f'Generating citations scopes table')
        QApplication.processEvents()
        the_tables.cit_scopes = process_cit_scopes(self.proteins, self.tbldir)
        status.showMessage(f'Generating citations table')
        QApplication.processEvents()
        the_tables.citations = process_citations(self.proteins, self.tbldir)
        status.showMessage(f'Generating comments table')
        QApplication.processEvents()
        the_tables.comments = process_comments(self.proteins, self.tbldir)
        status.showMessage(f'Generating organisms table')
        QApplication.processEvents()
        the_tables.organisms = process_organisms(self.proteins, self.tbldir)
        status.showMessage(f'Generating family equivalence table')
        QApplication.processEvents()
        generate_families_equivalence_table(the_tables, self.tbldir)
        # pickle.dump(the_tables, open(self.tbldump, 'wb'))
        self.tables = the_tables
        self.save()

    def download_structures(self, status: QStatusBar) -> None:

        # @timeout_retries(60, 5)
        def retrieve_best_sm():
            attempt = 1
            response = None
            while attempt < 10:
                try:
                    response = requests.get(
                        'https://swissmodel.expasy.org/repository/uniprot/' + accession.upper() + '.pdb',
                        timeout=30)
                    break
                except (ConnectTimeout, HTTPError, ReadTimeout, Timeout, ConnectionError):
                    attempt += 1
                    if attempt == 10:
                        raise ConnectionAbortedError
            if response.ok:
                model_name = os.path.join(the_dir, f"{accession}_SM.pdb")
                with open(model_name, 'w') as text_file:
                    text_file.write(response.text)

        total = len(self.proteins)
        for count, (accession, the_protein) in enumerate(self.proteins.items()):
            the_dir = os.path.join(self.structdir, accession)
            if the_protein.experimental_structures:
                for pdb in the_protein.experimental_structures:
                    if not pdb.downloaded:
                        status.showMessage(f'{count} of {total} Downloading structures for {accession}: {pdb.code}')
                        QApplication.processEvents()
                        os.makedirs(the_dir, exist_ok=True)
                        pdb.retrieve_file(the_dir)
                        pdb.downloaded = True
            elif the_protein.models and the_protein.models[0].downloaded is False:
                status.showMessage(f'{count} of {total} Downloading Model for {accession}')
                QApplication.processEvents()
                os.makedirs(the_dir, exist_ok=True)
                retrieve_best_sm()
                for model in the_protein.models:
                    model.downloaded = True
        self.save()
        status.showMessage(f'Done.')

    def gen_fam_seq(self, status: QStatusBar) -> None:
        db_list = self.tables.db['Database'].unique()
        for db in db_list:
            df = self.tables.db.loc[self.tables.db['Database'] == db]
            values = df['Value'].unique()
            names = [validate_string(n) for n in values]
            if len(set(names)) != len(values):
                raise ValueError('Ambiguous validated value in {db}')
            for i, (value, name) in enumerate(zip(values, names)):
                subdirectory = os.path.join(self.famdir, db, name)
                if not os.path.exists(subdirectory):
                    os.makedirs(subdirectory, exist_ok=True)
                if i % 1 == 0:
                    status.showMessage(f'Processing {name}')
                    QApplication.processEvents()
                text = ""
                hits = df.loc[df['Value'] == value]
                codes = hits['Uniprot'].unique()
                proteins = self.tables.sequences.loc[self.tables.sequences['Uniprot'].isin(codes)]
                for _, protein in proteins.iterrows():
                    text += seq2fasta(protein['Sequence'], protein['Uniprot'])
                filename = f'{db}_{name}.fasta'
                filename = os.path.join(subdirectory, filename)
                with open(filename, "wt") as handle:
                    handle.write(text)

        status.showMessage(f'Done.')

    def gen_meme(self, status: QStatusBar) -> None:

        db_list = self.tables.db['Database'].unique()
        jobs: List[MemeJob] = []
        for db in db_list:
            df = self.tables.db.loc[self.tables.db['Database'] == db]
            values = df['Value'].unique()
            family_names = [validate_string(n) for n in values]
            if len(set(family_names)) != len(values):
                raise ValueError('Ambiguous validated value in {db}')
            for i, (value, family_name) in enumerate(zip(values, family_names)):
                subdirectory = os.path.join(self.famdir, db, family_name)
                if not os.path.exists(subdirectory):
                    os.makedirs(subdirectory, exist_ok=True)
                # if i % 1 == 0:
                    # status.showMessage(f'Processing {family_name}')
                    # QApplication.processEvents()
                text = ""
                hits = df.loc[df['Value'] == value]
                codes = hits['Uniprot'].unique()
                proteins = self.tables.sequences.loc[self.tables.sequences['Uniprot'].isin(codes)]
                count = 0
                for _, protein in proteins.iterrows():
                    count += 1
                    text += seq2fasta(protein['Sequence'], protein['Uniprot'])
                filename = f'{db}_{family_name}.fasta'
                filename = os.path.join(subdirectory, filename)
                with open(filename, "wt") as handle:
                    handle.write(text)
                if count > 1:

                    os.makedirs(self.motivedir, exist_ok=True)
                    motivedir = os.path.join(self.motivedir, db, family_name)
                    os.makedirs(motivedir, exist_ok=True)
                    htmlfile = os.path.join(motivedir, 'meme.html')
                    txtfile = os.path.join(motivedir, 'meme.txt')
                    if os.path.isfile(txtfile) and os.path.isfile(htmlfile):
                        continue
                    jobs.append(MemeJob(filename, motivedir, config.meme_executable))
        status.showMessage(f'Processing Meme Motifs')
        QApplication.processEvents()
        p = Pool()
        _ = p.map(process_meme, jobs)
        status.showMessage(f'Done.')

    def gen_fam_struct(self, status: QStatusBar) -> None:
        # struct dir is where structures are, directory is where to put results
        db_list = self.tables.db['Database'].unique()
        for db in db_list:
            df = self.tables.db.loc[self.tables.db['Database'] == db]
            values = df['Value'].unique()
            names = [validate_string(n) for n in values]
            if len(set(names)) != len(values):
                raise ValueError('Ambiguous validated value in {db}')
            for i, (value, name) in enumerate(zip(values, names)):
                subdirectory = os.path.join(self.famstrdir, db, name)
                if not os.path.exists(subdirectory):
                    os.makedirs(subdirectory, exist_ok=True)
                if i % 10 == 0:
                    status.showMessage(f'Processing {name}')
                    QApplication.processEvents()
                text = ""
                hits = df.loc[df['Value'] == value]
                codes = hits['Uniprot'].unique()
                proteins = self.tables.sequences.loc[self.tables.sequences['Uniprot'].isin(codes)]
                for _, protein in proteins.iterrows():
                    if protein["Fragment"] is not '':
                        continue
                    text += seq2fasta(protein['Sequence'], protein['Uniprot'])
                    path = os.path.join(self.prepdir, protein['Uniprot'])
                    if os.path.exists(path):
                        ll = os.listdir(path)
                        for n in ll:
                            if n[-4:] in ['.PDB', '.pdb']:
                                the_name = os.path.join(path, n)
                                shutil.copy2(the_name, subdirectory)
                filename = f'{db}_{name}_nofragments.fasta'
                filename = os.path.join(subdirectory, filename)
                with open(filename, "wt") as handle:
                    handle.write(text)
        status.showMessage(f'Done.')

    def gen_summary(self, status: QStatusBar) -> None:

        with open(os.path.join(self.querydir, 'summary.txt'), 'wt') as out:

            n_pdbs = 0
            n_models = 0

            for _, p in self.proteins.items():
                if p.experimental_structures:
                    n_pdbs += 1
                elif p.models:
                    n_models += 1

            out.write(f'Query: <<{self.query}>>\n\n')
            out.write(f'{len(self.proteins)} proteins were found on UniProt database.\n')
            out.write(f'{n_pdbs} have at least 1 experimental structure on PDB database\n')
            out.write(f'{n_models} have a 3d model.\n')
        status.showMessage('Done.')
        QApplication.processEvents()
