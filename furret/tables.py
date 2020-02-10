import os
import pandas
from furret.utilities import Gos
from furret.protein import Protein
from typing import Dict, Optional


def process_keywords(protein_dict: Dict[str, Protein], output_dir: str, output_file: Optional[str] = 'keywords.xlsx')\
        -> pandas.DataFrame:

    keywords = []
    for i, p in protein_dict.items():
        for k in p.keywords:
            element = (p.accession, k.id, k.value, k.category)
            keywords.append(element)

    cols = ("Uniprot", "ID", "Keyword", "Category")

    df = pandas.DataFrame(keywords, columns=cols)
    df_count = df['Keyword'].value_counts()
    names = list(df['Category'].unique())
    counts = []

    for name in names:
        dd = df.loc[df['Category'] == name]
        counts.append(dd['Keyword'].value_counts())

    writer = pandas.ExcelWriter(os.path.join(output_dir, output_file), engine='xlsxwriter')
    df.to_excel(writer, sheet_name='Keywords (All)', index=False)
    df_count.to_excel(writer, sheet_name='Keywords (All Counts)')

    for name, count in zip(names, counts):
        count.to_excel(writer, sheet_name=name)

    writer.save()

    print(f'''File {output_file} written''')
    return df


def process_go(protein_dict, output_dir, output_file='go.xlsx'):
    molecular_function = []
    cellular_component = []
    biological_process = []

    for i, p in protein_dict.items():
        for go in p.go:
            element = (p.accession, go.id, go.value)
            if go.type == 'F':
                molecular_function.append(element)
            elif go.type == 'C':
                cellular_component.append(element)
            elif go.type == 'P':
                biological_process.append(element)

    cols = ("Uniprot", "ID", "GO")

    mf = pandas.DataFrame(molecular_function, columns=cols)
    cc = pandas.DataFrame(cellular_component, columns=cols)
    bp = pandas.DataFrame(biological_process, columns=cols)

    mf_counts = mf['GO'].value_counts()
    cc_counts = cc['GO'].value_counts()
    bp_counts = bp['GO'].value_counts()

    writer = pandas.ExcelWriter(os.path.join(output_dir, output_file), engine='xlsxwriter')

    mf.to_excel(writer, sheet_name='Molecular Function (All)', index=False)
    mf_counts.to_excel(writer, sheet_name='Molecular Function (Counts)')
    bp.to_excel(writer, sheet_name='Biological Process (All)', index=False)
    bp_counts.to_excel(writer, sheet_name='Biological Process (Counts)')
    cc.to_excel(writer, sheet_name='Cellular Component (All)', index=False)
    cc_counts.to_excel(writer, sheet_name='Cellular Component (Counts)')

    writer.save()
    print('Done with GO')
    return Gos(mf, bp, cc)


def process_links(protein_dict, output_dir, output_file='databases.xlsx'):

    links = []
    for i, p in protein_dict.items():
        for link in p.links:
            element = (p.accession, link.id, link.name, link.database)
            links.append(element)
    cols = ("Uniprot", "ID", "Value", "Database")
    df = pandas.DataFrame(links, columns=cols)
    names = list(df['Database'].unique())
    df_count = df['Value'].value_counts()
    counts = []
    dbs = []

    for name in names:
        dd = df.loc[df['Database'] == name]
        count = dd['Value'].value_counts()
        ids1 = [df.loc[df['Value'] == aa] for aa in count.index.values]
        ids = [x.iloc[0, 1] for x in ids1]
        s = pandas.Series(ids, index=count.index.values)
        the_count = pandas.concat([s, count], axis=1)
        the_count.columns = ['ID', 'Number']

        counts.append(the_count)

        dbs.append(dd)

    writer = pandas.ExcelWriter(os.path.join(output_dir, output_file), engine='xlsxwriter')
    df.to_excel(writer, sheet_name='Databases (All)', index=False)
    df_count.to_excel(writer, sheet_name='Databases (All Counts)')

    for name, count, db in zip(names, counts, dbs):
        count.to_excel(writer, sheet_name=name)

    writer.save()
    return df


def process_sequences(protein_dict, output_dir, output_file='sequences.xlsx'):
    sequences = []
    for p in protein_dict.values():
        the_structure = ''
        if p.experimental_structures:
            the_structure = 'PDB'
        if p.models:
            the_structure = 'Swiss Model'
        coverage, pdb = p.coverage if p.coverage else (0.0, '')
        element = (p.accession, p.fragment, p. precursor, the_structure, len(p.sequence), coverage, pdb, p.sequence)
        sequences.append(element)
    cols = ("Uniprot", "Fragment", "Precursor", "Structure", "Length", "Crystal Coverage", "Best PDB", "Sequence")

    se = pandas.DataFrame(sequences, columns=cols)
    se['Length_by_10'] = 10 * (pandas.to_numeric(se['Length'], errors='coerce') // 10 + 1)
    len_counts = se['Length_by_10'].value_counts().sort_index()
    se = se.drop(['Length_by_10'], axis=1)
    se['Length'] = se['Length'].astype(float)
    pr_counts = se['Precursor'].value_counts()
    fr_counts = se['Fragment'].value_counts()

    writer = pandas.ExcelWriter(os.path.join(output_dir, output_file), engine='xlsxwriter')

    se.to_excel(writer, sheet_name='Sequences (All)', index=False)
    pr_counts.to_excel(writer, sheet_name='Precursor (Counts)')
    fr_counts.to_excel(writer, sheet_name='Fragment (Counts)')
    len_counts.to_excel(writer, sheet_name='Lenght by 10 res')

    writer.save()

    return se


def process_pdb(protein_dict, output_dir, output_file='pdb.xlsx'):
    structures = []
    for p in protein_dict.values():
        for s in p.experimental_structures:
            structures.append((p.accession, s.code, s.method, s.resolution, s.coverage))
    cols = ('Uniprot', 'PDB', 'Method', 'Resolution', "Coverage")
    # noinspection DuplicatedCode
    dframe = pandas.DataFrame(structures, columns=cols)
    uniprot_counts = dframe['Uniprot'].value_counts()
    pdb_counts = dframe['PDB'].value_counts()
    method_counts = dframe['Method'].value_counts()

    writer = pandas.ExcelWriter(os.path.join(output_dir, output_file), engine='xlsxwriter')

    dframe.to_excel(writer, sheet_name='PDB (All)', index=False)
    uniprot_counts.to_excel(writer, sheet_name='Uniprot (Counts)')
    pdb_counts.to_excel(writer, sheet_name='PDB (Counts)')
    method_counts.to_excel(writer, sheet_name='Method (Counts)')

    writer.save()

    return dframe


def process_sm(protein_dict, output_dir, output_file='swiss_models.xlsx'):
    models = []
    for p in protein_dict.values():
        for m in p.models:
            models.append((p.accession, m.template, m.identity, m.oligo_state,
                           m.coverage, m.qmean, m.qmean_norm, m.gmqe))
    cols = ('Uniprot', 'Template', 'Identity', 'Oligo', 'Coverage', 'Qmean',
            'Qmean Norm', 'GMQE')

    df = pandas.DataFrame(models, columns=cols)
    uniprot_counts = df['Uniprot'].value_counts()
    templete_counts = df['Template'].value_counts()
    method_counts = df['Method'].value_counts()

    writer = pandas.ExcelWriter(os.path.join(output_dir, output_file), engine='xlsxwriter')

    df.to_excel(writer, sheet_name='Swiss Models (All)', index=False)
    uniprot_counts.to_excel(writer, sheet_name='Uniprot (Counts)')
    templete_counts.to_excel(writer, sheet_name='Template (Counts)')
    method_counts.to_excel(writer, sheet_name='Method (Counts)')

    writer.save()

    return df


def process_cit_scopes(protein_dict, output_dir, output_file='citation_scopes.xlsx'):
    scopes = []
    for p in protein_dict.values():
        for c in p.citations:
            scopes.append((p.accession, c.pubmed, c.doi, c.scope, c.title))
    cols = ('Uniprot', 'Pubmed', 'DOI', 'Scope', 'Title')

    df = pandas.DataFrame(scopes, columns=cols)
    uniprot_counts = df['Uniprot'].value_counts()
    templete_counts = df['Scope'].value_counts()

    writer = pandas.ExcelWriter(os.path.join(output_dir, output_file), engine='xlsxwriter')

    df.to_excel(writer, sheet_name='Scopes (All)', index=False)
    uniprot_counts.to_excel(writer, sheet_name='Uniprot (Counts)')
    templete_counts.to_excel(writer, sheet_name='Scope (Counts)')

    writer.save()

    return df


def process_citations(protein_dict, output_dir, output_file='citations.xlsx'):
    citations = []
    for p in protein_dict.values():
        for c in p.citations:
            citations.append((p.accession, c.pubmed, c.doi, c.title))
    cols = ('Uniprot', 'Pubmed', 'DOI', 'Title')

    df = pandas.DataFrame(citations, columns=cols)
    uniprot_counts = df['Uniprot'].value_counts()

    writer = pandas.ExcelWriter(os.path.join(output_dir, output_file), engine='xlsxwriter')

    df.to_excel(writer, sheet_name='Citations (All)', index=False)
    uniprot_counts.to_excel(writer, sheet_name='Uniprot (Counts)')

    writer.save()

    return df


def process_comments(protein_dict, output_dir, output_file='comments.xlsx'):
    comments = []
    for p in protein_dict.values():
        for c in p.comments:
            comments.append((p.accession, c.type, c.text))
    cols = ('Uniprot', 'Type', 'Text')

    df = pandas.DataFrame(comments, columns=cols)
    uniprot_counts = df['Uniprot'].value_counts()
    type_counts = df['Type'].value_counts()

    writer = pandas.ExcelWriter(os.path.join(output_dir, output_file), engine='xlsxwriter')

    df.to_excel(writer, sheet_name='Comments (All)', index=False)
    uniprot_counts.to_excel(writer, sheet_name='Uniprot (Counts)')
    type_counts.to_excel(writer, sheet_name='Type (Counts)')
    writer.save()

    return df


def process_organisms(protein_dict, output_dir, output_file='organisms.xlsx'):
    organisms = []
    for p in protein_dict.values():
        organisms.append((p.accession, p.organism))
    cols = ('Uniprot', 'Organism')
    df = pandas.DataFrame(organisms, columns=cols)
    orgnism_counts = df['Organism'].value_counts()

    writer = pandas.ExcelWriter(os.path.join(output_dir, output_file), engine='xlsxwriter')

    df.to_excel(writer, sheet_name='Organisms (All)', index=False)
    orgnism_counts.to_excel(writer, sheet_name='Organism (Counts)')

    writer.save()

    return df


def generate_families_equivalence_table(the_tables, table_dir, output_file='Pfam_identities.xlsx'):
    db_list = the_tables.db['Database'].unique()
    triples = []
    for db in db_list:
        df = the_tables.db.loc[the_tables.db['Database'] == db]
        # df contiene solo le entry di db

        values = df['Value'].unique()
        for family in values:
            hits = df.loc[df['Value'] == family]
            codes = hits['Uniprot'].unique()
            uniprot_set = set()
            # codice uniprot del database db che hanno valore family

            proteins = the_tables.sequences.loc[the_tables.sequences['Uniprot'].isin(codes)]
            for _, the_protein in proteins.iterrows():
                if the_protein["Fragment"] is not '':
                    continue
                uniprot_set.add(the_protein['Uniprot'])

            triples.append((db, family, uniprot_set))

    pfam = [(db, family, uniprot_set) for (db, family, uniprot_set) in triples if db == 'Pfam']
    ignore = []
    pfams = dict()
    for pf in pfam:
        if pf[1] in ignore:  # we already added this
            continue
        pfams[pf[1]] = []
        for fam in triples:
            if pf[0] == fam[0] and pf[1] == fam[1]:  # they are the same, skip!
                continue
            if pf[2] == fam[2]:  # we have an hit
                if fam[0] == pf[0]:  # ouch the hit is in pfam, add to ignore!
                    ignore.append(fam[1])
                pfams[pf[1]].append((fam[0], fam[1]))

    writer = pandas.ExcelWriter(os.path.join(table_dir, output_file), engine='xlsxwriter')
    for key, fams in pfams.items():
        cols = ("Database", "Family")
        df = pandas.DataFrame(fams, columns=cols)
        df.to_excel(writer, sheet_name=key, index=False)
    writer.save()
