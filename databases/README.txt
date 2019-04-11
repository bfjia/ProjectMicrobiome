//qiime2 needs to be on the same version the db is built with.

//Build Naive Bayes DB:
1. Download silva trainned classifier (make sure to go to the latest doc verion): https://docs.qiime2.org/2018.2/data-resources/ 
2. tar -xvf $_
3. done.

//To use:
qiime feature-classifier classify-sklearn --i-classifier "$classifierPath" --i-reads "$REPSEQ".qza --o-classification "$TAXONOMYSKLEARN" --p-n-jobs "$THREADS" --verbose



//Build VSearch DB
1. Download silva qiime2 compatible release from https://www.arb-silva.de/download/archive/qiime
2. tar -xvf $_ 
2. Make sequence artefact: qiime tools import --type FeatureData[Sequence] --input-path $SILVA_DIR/rep_set/rep_set_16S_only/99/99_otus_16S.fasta --output-path $databaseDir/silva.seq.qza
3. Make taxonomy artefact: qiime tools import --type FeatureData[Taxonomy] --input-path $SILVA_DIR/taxonomy/16S_only/99/consensus_taxonomy_all_levels.txt --output-path $databaseDir/silva.taxonomy.qza
4. done.

//To use
qiime feature-classifier classify-consensus-vsearch --i-query "$REPSEQ".qza --i-reference-reads "$databaseDir/silva.seq.qza" --i-reference-taxonomy "$databaseDir/silva.taxonomy.qza" --p-perc-identity 0.99 --o-classification "$TAXONOMYVSEARCH".qza --p-threads "$THREADS" --verbose
