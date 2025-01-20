#!/bin/bash

set -e
set -x

# curl -O 'http://sgd-archive.yeastgenome.org/curation/literature/gene_association.README'
# curl -O 'http://sgd-archive.yeastgenome.org/curation/literature/gene_association.sgd.gaf.gz'
# gunzip gene_association.sgd.gaf.gz

curl -O 'http://sgd-archive.yeastgenome.org/curation/literature/go_terms.README'
curl -O 'http://sgd-archive.yeastgenome.org/curation/literature/go_terms.tab'

curl -O 'http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/SGD_features.README'
curl -O 'http://sgd-archive.yeastgenome.org/curation/chromosomal_feature/SGD_features.tab'
