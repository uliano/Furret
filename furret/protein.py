import requests
import pandas
import numpy
from io import StringIO
import json
from furret.utilities import Link, Citation, Comment, Go, Keyword
from furret.chains import ChainGroups
from furret.structure import PDB, Model
from typing import List


class Protein:  # we look for swiss models ONLY IF no PDB is found
    # consider caching the response
    response = requests.get('https://www.uniprot.org/keywords/?query=*&format=tab&force=true&compress=no')
    if response.ok:
        text = response.text
        kw_table = pandas.read_csv(StringIO(text), sep='\t')
        del text
    else:
        # logger.error('Unable to access keyword table, skipping keyword processing.')
        kw_table = None
    del response

    def __init__(self, entry, database_list=('Gene3D', 'InterPro', 'Pfam', 'SUPFAM', 'PROSITE', 'PRINTS',
                                             'SMART', 'TIGRFAMs', 'CDD', 'PANTHER', 'PIRSF')):
        def to_float(s: str) -> str:
            try:
                result = float(s)
            except ValueError:
                result = numpy.nan
            return result

        def get_elements_of_type(the_list, the_type=None) -> List[str]:
            assert the_type
            elements = []
            for the_element in the_list:
                if the_element['@type'] != the_type:
                    continue
                elements.append(the_element)
            return elements

        def get_value_for_label_of_type(the_list, the_label=None, the_type=None) -> str:
            assert the_label and the_type
            elements = get_elements_of_type(the_list, the_type=the_type)
            if len(elements) == 0:
                raise KeyError(f'''@type {the_type} not found in {self.accession}''')
            elif len(elements) == 1:
                return elements[0][the_label]
            else:
                raise KeyError(f'''More than one @type {the_type} found in {self.accession}''')

        def get_comments_list(the_list) -> List[Comment]:
            comments = []
            for the_element in the_list:
                try:
                    the_type = the_element['@type']
                except TypeError:
                    continue
                text = None
                try:
                    text = the_element['text']['#text']
                except (TypeError, KeyError) as _:
                    try:
                        text = the_element['text']
                    except KeyError:
                        pass
                if text:
                    the_comment = Comment(the_type, text)
                    comments.append(the_comment)
            return comments

        def get_citation_list(the_list) -> List[Citation]:
            citations = []
            for the_element in the_list:
                citation = the_element['citation']
                title = citation['title'] if 'title' in citation else None
                scope = the_element['scope'] if 'scope' in the_element else None
                try:
                    pubmed = get_value_for_label_of_type(citation['dbReference'], the_label='@id', the_type='PubMed')
                except KeyError:
                    pubmed = None
                try:
                    doi = get_value_for_label_of_type(citation['dbReference'], the_label='@id', the_type='DOI')
                except KeyError:
                    doi = None
                if not pubmed and not doi:
                    continue
                the_citation = Citation(title=title, scope=scope, pubmed=pubmed, doi=doi)
                citations.append(the_citation)
            return citations

        def get_go(go_entries) -> List[Go]:
            go_list = []
            for the_element in go_entries:
                go_id = the_element['@id']
                go = get_value_for_label_of_type(the_element['property'],
                                                 the_type='term', the_label='@value').split(':')
                go_type = go[0]
                go_value = go[1]
                the_go = Go(go_id, go_type, go_value)
                go_list.append(the_go)
            return go_list

        def get_keywords(keyword_entries) -> List[Keyword]:
            keyword_list = []
            for keyword in keyword_entries:
                key_id = keyword['@id']
                row = self.kw_table.loc[self.kw_table['Keyword ID'] == key_id]
                category = row.iloc[0]['Category']
                value = keyword['#text']
                the_keyword = Keyword(key_id, value, category)
                keyword_list.append(the_keyword)
            return keyword_list

        def get_uniprot_pdb(pdb_entries) -> List[PDB]:
            pdb_list = []
            for the_element in pdb_entries:
                pdb_id = the_element['@id']
                properties = the_element['property']
                method = ''
                resolution = numpy.nan
                the_chains = None
                for the_property in properties:
                    property_type = the_property['@type']
                    property_value = the_property['@value']
                    if property_type == 'method':
                        method = the_property['@value']
                    elif property_type == 'resolution':
                        resolution = to_float(property_value)
                    elif property_type == 'chains':
                        chain_groups = property_value.split(',')
                        the_chains = ChainGroups(chain_groups)
                    else:
                        assert False, f'Unknown PDB type {property_type} in PDB' + \
                                      f' code {pdb_id} in accession{self.accession}'
                the_pdb = PDB(self.accession, code=pdb_id, sequence=self.sequence, the_chains=the_chains)
                the_pdb.resolution = resolution
                the_pdb.method = method
                pdb_list.append(the_pdb)
            return pdb_list

        def get_swiss_models() -> List[Model]:
            answer = requests.get(f'https://swissmodel.expasy.org/repository/uniprot/{self.accession}.json')
            j = json.load(StringIO(answer.text))
            structures = j['result']['structures']
            if len(structures) == 0:
                return []
            swiss_model_list = []
            for struct in structures:
                swiss_model_list.append(Model(uniprot=self.accession, sequence=self.sequence, data=struct))
            return swiss_model_list

        def get_best_coverage(structures):
            if len(structures) == 0:
                return
            cov = structures[0].coverage
            idx = 0
            for i, the_structure in enumerate(structures):
                if the_structure.coverage > cov:
                    cov = the_structure.coverage
                    idx = i
            return cov, structures[idx].code

        self.accession: str = entry['accession'][0]
        try:
            fullname = entry['protein']['recommendedName']['fullName']
        except KeyError:
            try:
                fullname = entry['protein']['submittedName']['fullName']
            except KeyError:
                fullname = ''
        # noinspection PyTypeChecker
        self.name = fullname if type(fullname) is str else fullname['#text']
        self.organism: str = get_value_for_label_of_type(entry['organism']['name'], the_label='#text', the_type='scientific')
        self.citations = get_citation_list(entry['reference'])
        try:
            self.comments = get_comments_list(entry['comment'])
        except KeyError:
            self.comments = []
        dbreferences = entry.get('dbReference', [])
        self.go = get_go(get_elements_of_type(dbreferences, the_type='GO'))
        self.protein_existence: str = entry['proteinExistence']['@type']
        keywords = entry.get('keyword', [])
        self.keywords = get_keywords(keywords)
        self.sequence: str = entry['sequence']['#text'].replace('\n', '')
        self.precursor: str = entry['sequence'].get('@precursor', '')
        self.fragment: str = entry['sequence'].get('@fragment', '')
        self.experimental_structures = get_uniprot_pdb(get_elements_of_type(dbreferences, the_type='PDB'))
        self.models = get_swiss_models() if not self.experimental_structures else []
        self.coverage = get_best_coverage(self.experimental_structures)
        self.links = []
        for db in database_list:
            for element in get_elements_of_type(dbreferences, the_type=db):
                self.links.append(Link(database=db, data=element))

    def best_model(self):
        best = None
        gmqe_max = -1.0
        for model in self.models:
            if model.gmqe > gmqe_max:
                best = model
                gmqe_max = model.gmqe
        return best
