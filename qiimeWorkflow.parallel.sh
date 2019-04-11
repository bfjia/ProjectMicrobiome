#!/bin/bash -l
#
#SBATCH --workdir=/home/jjjjia/scratch/beluga_newprimer/
#SBATCH --account=rrg-fiona-ad
#SBATCH --job-name=beluga_workflow
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=3500M
#SBATCH --time=71:30:00
#SBATCH --mail-user=<bja20@sfu.ca>
#SBATCH --mail-type=ALL

source activate qiime2-2018.2

#qiime2 workflow
#parameters:  $1=outputName, $2=threads, $3=trim-left-forward(22), $4=trim-left-reverse(22), 
#parameters: $5=trunc-len-forward(200), $6=trunc-len-reverse(180), $7=depthFrequency(5000), $8=diversitysamplingdepth(40000)
#things to note:
#$1=MANIFEST(CSV: sample-id,absolute-filepath,direction), must be in unix format
#$9=metadata(TSV: sampleid,x,y,z), must be in utf8 and unix format

#LETS ASSIGN SOME VARIABLES
RESULTNAME="${1:-16sAnalysis}"
THREADS="${2:-48}"
TRIMLF="${3:-20}"
TRIMLR="${4:-230}"
TRIMRF="${5:-20}"
TRIMRR="${6:-200}"
DEPTH="${7:-1}"
DIVERSITYDEPTH="${8:-10000}"

echo $RESULTNAME
echo $THREADS
echo $TRIMLF
echo $TRIMLR
echo $TRIMRF
echo $TRIMRR
echo $DEPTH
echo $DIVERSITYDEPTH

#sbatch qiimeWorkflow.sh beluga_newprimer 48 25 25 220 200 1 40000

mkdir $RESULTNAME
cd $RESULTNAME

#This command imports the FASTQ files into a QIIME artifact
#For 20 million reads, this took 2hrs. 300MB RAM @ 32CPU
echo "importing sequences"

manifestPath="../csv/manifest.csv"
#paired end sequence
#qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path "$manifestPath" --output-path "$RESULTNAME" --source-format PairedEndFastqManifestPhred33
#forward only
qiime tools import --type 'SampleData[SequencesWithQuality]' --input-path "$manifestPath" --output-path "$RESULTNAME" --source-format SingleEndFastqManifestPhred33

#Using DADA2 to analyze quality scores of 10000 random samples

READQZA="$RESULTNAME".qza
qiime demux summarize --p-n 10000 --i-data "$READQZA" --o-visualization "$RESULTNAME".qzv

echo "Denoising sequences"
#Denoising with DADA2. Using quality score visualizations, you can choose trunc-len-f and trunc-len-r (note: sequences < trunc-len in length are discarded!)
# The drop-off for the forward reads was not so bad, but there is a significant drop-off in quality for the reverse reads, so let's trim 10bp
#For 20 million reads, this took 3.5hours, 24GB RAM @ 32CPU

FEATURETABLE="$RESULTNAME".featuretable
REPSEQ="$RESULTNAME".repseq

#paired end denoise
#qiime dada2 denoise-paired --p-n-threads "$THREADS" --verbose --i-demultiplexed-seqs "$READQZA" --o-table "$FEATURETABLE" --o-representative-sequences "$REPSEQ" --p-trim-left-f "$TRIMLF" --p-trim-left-r "$TRIMLR" --p-trunc-len-f "$TRIMRF" --p-trunc-len-r "$TRIMRR" --verbose

#forward only denoise
qiime dada2 denoise-single --p-n-threads "$THREADS" --verbose --i-demultiplexed-seqs "$READQZA" --o-table "$FEATURETABLE" --o-representative-sequences "$REPSEQ" --p-trim-left "$TRIMLF" --p-trunc-len "$TRIMRF" --verbose

#This visualization shows us the sequences/sample spread
qiime feature-table summarize --i-table "$FEATURETABLE".qza --o-visualization "$FEATURETABLE".qzv

#Filter out sequences with few samples, use the number of reads from the lowest sample from the $resultname.qzv file
FILTEREDFEATURETABLE="$FEATURETABLE".filtered."$DEPTH"
qiime feature-table filter-samples --i-table "$FEATURETABLE".qza --p-min-frequency "$DEPTH" --o-filtered-table "$FILTEREDFEATURETABLE"

qiime feature-table summarize --i-table "$FILTEREDFEATURETABLE".qza --o-visualization $FILTEREDFEATURETABLE.qzv

echo "assigning taxonomy"
sklearnSilvaPath="../databases/silva119.sklearn.qza"
vsearchSilvaPath="../databases/silva132.vsearch.qza"
vsearchSilvaTaxonomyPath="../databases/silva132taxonomy.vsearch.qza"
TAXONOMYSKLEARN="$RESULTNAME".taxonomy.sklearn
TAXONOMYVSEARCH="$RESULTNAME".taxonomy.vsearch


#Classify against it with Naive Bayes
#for 20million reads, this took 30minutes, 28GB ram @ 32CPU
qiime feature-classifier classify-sklearn --i-classifier "$sklearnSilvaPath" --i-reads "$REPSEQ".qza --o-classification "$TAXONOMYSKLEARN" --p-n-jobs "$THREADS" --verbose

#Taxa bar plots
qiime taxa barplot --i-table "$FEATURETABLE".qza --i-taxonomy "$TAXONOMYSKLEARN".qza --m-metadata-file "$metadataPath" --o-visualization "$TAXONOMYSKLEARN".qzv
qiime taxa barplot --i-table "$FILTEREDFEATURETABLE".qza --i-taxonomy "$TAXONOMYSKLEARN".qza --m-metadata-file "$metadataPath" --o-visualization "$TAXONOMYSKLEARN".qzv

source deactivate #because we are running two different versions of qiime, need to change environment

#run vsearch for taxonomy classification
source activate qiime2-2018.8

qiime feature-classifier classify-consensus-vsearch --i-query "$REPSEQ".qza --i-reference-reads "$vsearchSilvaPath" --i-reference-taxonomy "$vsearchSilvaTaxonomyPath" --p-perc-identity 0.99 --o-classification "$TAXONOMYVSEARCH".qza --p-threads "$THREADS" --verbose

#Taxa bar plots
metadataPath="../csv/metadata.tsv"

qiime taxa barplot --i-table "$FEATURETABLE".qza --i-taxonomy "$TAXONOMYVSEARCH".qza --m-metadata-file "$metadataPath" --o-visualization "$TAXONOMYVSEARCH".qzv
qiime taxa barplot --i-table "$FILTEREDFEATURETABLE".qza --i-taxonomy "$TAXONOMYVSEARCH".qza --m-metadata-file "$metadataPath" --o-visualization "$TAXONOMYVSEARCH".qzv

source deactivate

echo "constructing tree"
#Steps for generating a phylogenetic tree

ROOTEDTREE="$RESULTNAME".rootedTree

qiime alignment mafft --i-sequences "$REPSEQ".qza --o-alignment "$RESULTNAME".alignedrepseq
qiime alignment mask --i-alignment "$RESULTNAME".alignedrepseq.qza --o-masked-alignment "$RESULTNAME".maskedalignedrepseq
qiime phylogeny fasttree --i-alignment "$RESULTNAME".maskedalignedrepseq.qza --o-tree "$RESULTNAME".unrooted_tree
qiime phylogeny midpoint-root --i-tree "$RESULTNAME".unrooted_tree.qza --o-rooted-tree "$ROOTEDTREE"

echo "calculating diversity"
#Generate alpha/beta diversity measures at 41000 sequences/sample
#Also generates PCoA plots automatically
#non filtered first
ALPHADIV="$RESULTNAME".AlphaDiversity."$DIVERSITYDEPTH"
qiime diversity core-metrics-phylogenetic --i-phylogeny "$ROOTEDTREE.qza" --i-table "$FEATURETABLE".qza --p-sampling-depth "$DIVERSITYDEPTH" --m-metadata-file "$metadataPath" --output-dir "$ALPHADIV"
#Test for between-group differences
qiime diversity alpha-group-significance --i-alpha-diversity "$ALPHADIV"/faith_pd_vector.qza --m-metadata-file "$metadataPath" --o-visualization "$ALPHADIV".alpha_PD_significance.qzv
qiime diversity alpha-group-significance --i-alpha-diversity "$ALPHADIV"/shannon_vector.qza --m-metadata-file "$metadataPath" --o-visualization "$ALPHADIV".alpha_shannon_significance.qzv
#qiime diversity beta-group-significance --i-distance-matrix "$ALPHADIV"/bray_curtis_distance_matrix.qza --m-metadata-file "$METADATA" --m-metadata-column bee_type --o-visualization beta_bray_beetype_significance
#Alpha rarefaction curves show taxon accumulation as a function of sequence depth
qiime diversity alpha-rarefaction --i-phylogeny "$ROOTEDTREE".qza --i-table "$FEATURETABLE".qza --p-max-depth "$DIVERSITYDEPTH" --o-visualization "$ALPHADIV".alpha_rarefaction.qzv --m-metadata-file "$metadataPath"

#then filtered
ALPHADIV="$RESULTNAME".AlphaDiversity."$DIVERSITYDEPTH"."$DEPTH"
qiime diversity core-metrics-phylogenetic --i-phylogeny "$ROOTEDTREE.qza" --i-table "$FILTEREDFEATURETABLE".qza --p-sampling-depth "$DIVERSITYDEPTH" --m-metadata-file "$metadataPath" --output-dir "$ALPHADIV"

#Test for between-group differences
qiime diversity alpha-group-significance --i-alpha-diversity "$ALPHADIV"/faith_pd_vector.qza --m-metadata-file "$metadataPath" --o-visualization "$ALPHADIV".alpha_PD_significance.qzv
qiime diversity alpha-group-significance --i-alpha-diversity "$ALPHADIV"/shannon_vector.qza --m-metadata-file "$metadataPath" --o-visualization "$ALPHADIV".alpha_shannon_significance.qzv
#qiime diversity beta-group-significance --i-distance-matrix "$ALPHADIV"/bray_curtis_distance_matrix.qza --m-metadata-file "$METADATA" --m-metadata-column bee_type --o-visualization beta_bray_beetype_significance
#Alpha rarefaction curves show taxon accumulation as a function of sequence depth
qiime diversity alpha-rarefaction --i-phylogeny "$ROOTEDTREE".qza --i-table "$FILTEREDFEATURETABLE".qza --p-max-depth "$DIVERSITYDEPTH" --o-visualization "$ALPHADIV".alpha_rarefaction.qzv --m-metadata-file "$metadataPath"



echo "finishing up"
mkdir visualizations
mv *.qzv visualizations
mv "$ALPHADIV"/*.qzv visualizations
tar -zcvf visualizations.tar.gz visualizations/*

#get a tsv file to export to STAMP
#start with non filtered
mkdir stamptables.sklearn
cd stamptables.sklearn
LEVEL=3
for LEVEL in 3 4 5 6 7;
do
	STAMPNAME="$RESULTNAME".level"$LEVEL".sklearn.taxontable
	qiime taxa collapse --i-table "$FEATURETABLE".qza --i-taxonomy "$TAXONOMYSKLEARN".qza --p-level "$LEVEL" --o-collapsed-table "$STAMPNAME";
	qiime tools export --output-dir "$STAMPNAME".biom "$STAMPNAME".qza	;
	biom convert -i "$STAMPNAME".biom/feature-table.biom -o "$STAMPNAME".tsv --to-tsv;
	tail -n +2 "$STAMPNAME".tsv
done;
cd ../

#filtered
mkdir stamptables."$DEPTH".sklearn
cd stamptables."$DEPTH".sklearn
LEVEL=3
for LEVEL in 3 4 5 6 7;
do
	STAMPNAME="$RESULTNAME"."$DEPTH".level"$LEVEL".sklearn.taxontable
	qiime taxa collapse --i-table "$FILTEREDFEATURETABLE".qza --i-taxonomy "$TAXONOMYSKLEARN".qza --p-level "$LEVEL" --o-collapsed-table "$STAMPNAME";
	qiime tools export --output-dir "$STAMPNAME".biom "$STAMPNAME".qza	;
	biom convert -i "$STAMPNAME".biom/feature-table.biom -o "$STAMPNAME".tsv --to-tsv;
	tail -n +2 "$STAMPNAME".tsv
done;
cd../

source deactivate
source activate qiime2-2018.8

mkdir stamptables.vsearch
cd stamptables.vsearch
LEVEL=3
for LEVEL in 3 4 5 6 7;
do
	STAMPNAME="$RESULTNAME".level"$LEVEL".vsearch.taxontable
	qiime taxa collapse --i-table "$FEATURETABLE".qza --i-taxonomy "$TAXONOMYVSEARCH".qza --p-level "$LEVEL" --o-collapsed-table "$STAMPNAME";
	qiime tools export --output-path "$STAMPNAME".biom "$STAMPNAME".qza	;
	biom convert -i "$STAMPNAME".biom/feature-table.biom -o "$STAMPNAME".tsv --to-tsv;
	tail -n +2 "$STAMPNAME".tsv	
done;
cd ../

#vsearch
#filtered
mkdir stamptables."$DEPTH".vsearch
cd stamptables."$DEPTH".vsearch
LEVEL=3
for LEVEL in 3 4 5 6 7;
do
	STAMPNAME="$RESULTNAME"."$DEPTH".level"$LEVEL".vsearch.taxontable
	qiime taxa collapse --i-table "$FILTEREDFEATURETABLE".qza --i-taxonomy "$TAXONOMYVSEARCH".qza --p-level "$LEVEL" --o-collapsed-table "$STAMPNAME";
	qiime tools export --output-path "$STAMPNAME".biom "$STAMPNAME".qza	;
	biom convert -i "$STAMPNAME".biom/feature-table.biom -o "$STAMPNAME".tsv --to-tsv;
	tail -n +2 "$STAMPNAME".tsv	
done;
cd ../

source deactivate

#copy phyloseq related objects 
#required inputs:
#1. metadata in utf8
#2. feature table
#3. taxonomy
#4. rooted tree

source load r3.5.1


echo "done"
