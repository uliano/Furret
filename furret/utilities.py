import multiprocessing.pool
from dataclasses import dataclass
import functools
# noinspection PyPackageRequirements
from Bio.Entrez import efetch, read
import xmltodict
import furret.config as config
import pandas

from typing import Optional


def fetch_abstract(pmid):
    handle = efetch(db='pubmed', id=pmid, retmode='xml')
    xml_data = read(handle)
    try:
        article = xml_data['PubmedArticle'][0]['MedlineCitation']['Article']
        abstract = article['Abstract']['AbstractText'][0]
        title = article['ArticleTitle']
        return abstract, title
    except (IndexError, KeyError) as _:
        return None


def validate_string(s):
    valid_chars = '-_.() abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789'
    return ''.join(c for c in s if c in valid_chars)


def seq2fasta(seq: str, label: str) -> str:

    def group_by(s, size):
        res = ""
        length = len(s)
        steps = list(range(0, length, size))
        last, *steps = steps
        for step in steps:
            res += s[last:step] + '\n'
            last, *steps = steps
        res += s[last:length] + '\n'
        return res

    return f'>{label}\n' + group_by(seq, 60)


def timeout_retries(max_timout, max_retries):
    def timeout_decorator(function):
        @functools.wraps(function)
        def wrapper(*args, **kwargs):
            pool = multiprocessing.pool.ThreadPool(processes=1)
            async_results = pool.apply_async(function, args, kwargs)
            for _ in range(max_retries):
                try:
                    return async_results.get(max_timout)
                except TimeoutError:
                    continue
            else:
                raise TimeoutError(f'''No answer after {max_retries} retries  when calling {function}''')
        return wrapper
    return timeout_decorator


def find_uniprot(codes, protein_dict):
    if len(codes) == 0:
        codes = [codes]
    result = []
    for p in protein_dict.values():
        if p.accession in codes:
            result.append(p)
    return result


class Obj:
    pass


@dataclass
class Gos:
    molecular_function: pandas.DataFrame
    biological_process: pandas.DataFrame
    cellular_component: pandas.DataFrame


class Link:
    def __init__(self, database='', data=None):
        assert data and database
        self.database: str = database
        self.id: str = data['@id']
        self.name: str = data['property'][0]['@value']


@dataclass
class PubMed:
    id: str

    def retrieve_abstract(self):
        # noinspection PyUnresolvedReferences
        answer = requests.get('https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=' +
                              f'{self.id}&tool=my_tool&email={config.entrez_email}&retmode=xml')
        pm = xmltodict.parse(answer.text)
        abstract = pm['PubmedArticleSet']['PubmedArticle']['MedlineCitation']['Article']['Abstract']['AbstractText']
        return abstract


@dataclass
class Citation:
    pubmed: Optional[PubMed] = ''
    doi: Optional[str] = ''
    title: Optional[str] = ''
    scope: Optional[str] = ''


@dataclass
class Comment:
    type: str
    text: str


@dataclass
class Go:
    id: str
    type: str
    value: str


@dataclass
class Keyword:
    id: str
    value: str
    category: str


def format_filename(filename):
    return "".join([c for c in filename if c.isalpha() or c.isdigit() or c == ' ']).rstrip()
