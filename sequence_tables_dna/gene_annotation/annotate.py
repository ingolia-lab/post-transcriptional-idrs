import pandas as pd

# http://sgd-archive.yeastgenome.org/curation/literature/go_terms.tab
# http://sgd-archive.yeastgenome.org/curation/literature/gene_association.sgd.gaf.gz
# http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab

# Fig 2a
go_terms = [ 'GO:0000288', 'GO:0061157', 'GO:0043488', 'GO:0061014', 'GO:0017148' ]

class GOTerms:
    def __init__(self, go_terms_path='go_terms.tab'):
        self.terms = pd.read_csv(go_terms_path, sep='\t', index_col=False, comment='!',
            names = ['GOID', 'GO_Term', 'GO_Aspect', 'GO_Term_Definition'])

    def term(self, goid):
        return self.terms.loc[ self.terms['GOID'] == goid, 'GO_Term' ]

# class GOAssoc:
#     def __init__(self, gaf_path='gene_association.sgd.gaf'):
#         self.gaf = pd.read_csv(gaf_path, sep='\t', index_col=False, comment='!',
#             names = ['DB', 'DB_Object_ID', 'DB_Object_Symbol', 'NOT',
#                 'GO_ID', 'DB_Reference', 'Evidence', 'With',
#                 'Aspect', 'Name', 'Synonym', 'Object_Type',
#                 'taxon', 'Date', 'Assigned_by'])

#     def gene_assocs(self, gene):
#         assocs = self.gaf.loc[ self.gaf['DB_Object_Symbol'] == gene, 'GO_ID' ].unique()
#         return assocs

class PantherFile:
    def __init__(self, goid):
        gonum = goid.replace('GO:', '')
        file = f'pantherGeneList {gonum}.txt'
        self.gene_list = pd.read_csv(file, sep='\t', index_col=False,
            header=None, skip_blank_lines=True)
        self.genes = self.gene_list[1].str.upper().unique()
        print(f'Read {len(self.genes)} genes in {goid} from {file}')

    def contains_gene(self, gene):
        return gene in self.genes

class SGD:
    def __init__(self, sgd_path='SGD_features.tab'):
        self.sgd = pd.read_csv(sgd_path, sep='\t', index_col=False,
                          names=["sgdid", "type", "qual", "name", "gene", "alias",
                                      "parent", "sgdid2", "chrom", "start", "end",
                                      "strand", "genpos", "cver", "sver", "desc"])
        self.sgd = self.sgd[ self.sgd['name'].notna() ].copy()
        self.gene_to_yorf = pd.Series(data=self.sgd['name'].array, index=self.sgd['gene'])
        self.yorf_to_gene = pd.Series(data=self.sgd['gene'].array, index=self.sgd['name'])
        self.yorf_to_desc = pd.Series(data=self.sgd['desc'].array, index=self.sgd['name'])

    def gene_name(self, yorf):
        gene = self.yorf_to_gene.get(yorf, default=yorf)
        return (yorf if pd.isna(gene) else gene)

    def yorf(self, gene):
        if gene in self.sgd['name'].values:
            return gene
        else:
            yorf = self.gene_to_yorf.get(gene)
            if yorf is None:
                raise RuntimeError(f'Unknown gene {gene}')
            return yorf

sgd = SGD()
go = GOTerms()

orfs = pd.read_csv('repressive_orfs.csv', index_col=False)
orfs.rename(columns={'ORF': 'protein'}, inplace=True)
orfs['gene'] = orfs['protein'].str.replace(r'p$', '', regex=True).str.upper()
orfs['yorf'] = orfs['gene'].map(sgd.yorf)

for goid in go_terms:
    term = go.term(goid)
    print(f'{goid}\t{term}')
    pf = PantherFile(goid)
    orfs[goid] = orfs['gene'].map(pf.contains_gene)

orfs['desc'] = orfs['yorf'].map(sgd.yorf_to_desc)

orfs.to_csv('repressive_orfs_annot.csv', index=False)

fragments = pd.read_csv('repressive_fragments.csv', index_col=False)
fragments.rename(columns={'yorf': 'protein'}, inplace=True)
fragments['gene'] = fragments['protein'].str.replace(r'p$', '', regex=True).str.upper()
fragments['yorf'] = fragments['gene'].map(sgd.yorf)

for goid in go_terms:
    term = go.term(goid)
    print(f'{goid}\t{term}')
    pf = PantherFile(goid)
    fragments[goid] = fragments['gene'].map(pf.contains_gene)

fragments['desc'] = fragments['yorf'].map(sgd.yorf_to_desc)

fragments.to_csv('repressive_fragments_annot.csv', index=False)
