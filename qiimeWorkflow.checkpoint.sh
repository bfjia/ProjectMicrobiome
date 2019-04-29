#!/bin/bash -le
#
#SBATCH --workdir=/home/jjjjia/scratch/
#SBATCH --account=rrg-fiona-ad
#SBATCH --job-name=qiime2_workflow
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --mem-per-cpu=3500M
#SBATCH --time=71:50:00
#SBATCH --mail-user=<bja20@sfu.ca>
#SBATCH --mail-type=ALL

#qiime2 workflow
#parameters:  $1=outputName, $2=threads, $3=trim-left-forward(22), $4=trim-left-reverse(22), 
#parameters: $5=trunc-len-forward(200), $6=trunc-len-reverse(180), $7=depthFrequency(5000), 
#parameters: $8=diversitysamplingdepth(40000), $9=areReadsForwardOnly(false)
#things to note:
#$1=MANIFEST(CSV: sample-id,absolute-filepath,direction), must be in unix format
#$9=metadata(TSV: sampleid,x,y,z), must be in utf8 and unix format
#run command: sbatch qiimeWorkflow.sh frog_complete 48 25 25 220 200 1 40000 false

#LETS ASSIGN SOME VARIABLES
RESULTNAME="$1"
THREADS="$2"
TRIMLF="$3"
TRIMLR="$4"
TRIMRF="$5"
TRIMRR="$6"
DEPTH="$7"
DIVERSITYDEPTH="$8"
FORWARDONLY="$9"

echo $RESULTNAME
echo $THREADS
echo $TRIMLF
echo $TRIMLR
echo $TRIMRF
echo $TRIMRR
echo $DEPTH
echo $DIVERSITYDEPTH
echo ${FORWARDONLY,,} #requires bash4

#checkpoints
checkPointFile=./checkpoint
restartedFromCheckpoint="false"
#database
dbDir=~/project/jjjjia/databases/qiimeWorkflowDB
#csv
manifestPath="../csv/manifest.csv"
metadataPath="../csv/metadata.tsv"
#import
READQZA="$RESULTNAME".qza
#denoise
FEATURETABLE="$RESULTNAME".featuretable
REPSEQ="$RESULTNAME".repseq
DADASTAT="$RESULTNAME".dadaStat
FILTEREDFEATURETABLE="$FEATURETABLE".filtered."$DEPTH"
#sklearn
sklearnSilvaPath="../databases/silva119.sklearn.qza"
TAXONOMYSKLEARN="$RESULTNAME".taxonomy.sklearn
FILTEREDTAXONOMYSKLEARN="$RESULTNAME".taxonomy.sklearn.filtered."$DEPTH"
#vsearch
vsearchSilvaPath="../databases/silva132.vsearch.qza"
vsearchSilvaTaxonomyPath="../databases/silva132taxonomy.vsearch.qza"
TAXONOMYVSEARCH="$RESULTNAME".taxonomy.vsearch
FILTEREDTAXONOMYVSEARCH="$RESULTNAME".taxonomy.vsearch.filtered."$DEPTH"
#tree
ROOTEDTREE="$RESULTNAME".rootedTree
#diversity
ALPHADIV="$RESULTNAME".AlphaDiversity."$DIVERSITYDEPTH"
FILTEREDALPHADIV="$RESULTNAME".AlphaDiversity."$DIVERSITYDEPTH"."$DEPTH"



#function as a hack for labels and goto statements.
#labels start with '#' and end with ':' to avoid syntax errors. e.g. "#LabelName:"
function jumpTo ()
{
    label=$1    
    cmd=$(sed -n "/#$label:/{:a;n;p;ba};" $0 | grep -v ':$')
    #echo "$cmd"
    eval "$cmd"
    exit
}

#checkpoint logics

if [ -f $checkPointFile ]; then
    echo "checkpoint file found..."
    step=`cat $checkPointFile`;
    echo $step
    if [ $step = "finish" ]; then
        echo "previous analysis completed without error, remove result folder to restart analysis. exiting"
        exit 0
    else
        echo "attempting to restart workflow from $step"
        restartedFromCheckpoint="true"
        jumpTo $step
    fi
else
    #no checkpoint, start fresh
    jumpTo begin
fi

#begin:
echo "begin" > $checkPointFile

#restarted from last checkpoint, needs to remove files
if [ $restartedFromCheckpoint = "true" ]; then
    echo "removing residual files from previous run"
    rm -rf ./databases
    rm -rf $RESULTNAME
fi

source activate qiime2-2018.2

#lets move the database to current directory
date
cp -R $dbDir ./databases

mkdir $RESULTNAME
cd $RESULTNAME

#unix-ify all the csv
dos2unix ../csv/*
#assign csv path variables


#import:
echo "import" > $checkPointFile


#restarted from last checkpoint, needs to remove files
if [ $restartedFromCheckpoint = "true"]; then
    echo "removing residual files from previous run"
    rm -rf "$RESULTNAME".qza
    rm -rf "$RESULTNAME".qzv
fi

#paired end sequence
#This command imports the FASTQ files into a QIIME artifact
#For 20 million reads, this took 2hrs. 300MB RAM @ 32CPU
echo "importing sequences"
date

if [ $FORWARDONLY = "false"]; then
    qiime tools import --type 'SampleData[PairedEndSequencesWithQuality]' --input-path "$manifestPath" --output-path "$RESULTNAME" --source-format PairedEndFastqManifestPhred33
else
#forward only
    qiime tools import --type 'SampleData[SequencesWithQuality]' --input-path "$manifestPath" --output-path "$RESULTNAME" --source-format SingleEndFastqManifestPhred33
fi

#Using DADA2 to analyze quality scores of 10000 random samples
qiime demux summarize --p-n 10000 --i-data "$READQZA" --o-visualization "$RESULTNAME".qzv


#denoise:
echo "denoise" > $checkPointFile 


#restarted from last checkpoint, needs to remove files
if [ $restartedFromCheckpoint = "true"]; then
    echo "removing residual files from previous run"
    rm -rf "$FEATURETABLE".qza
    rm -rf "$FEATURETABLE".qzv
    rm -rf "$FILTEREDFEATURETABLE".qza
    rm -rf "$FILTEREDFEATURETABLE".qzv
    rm -rf "$REPSEQ".qza
fi

echo "Denoising sequences"
date
#Denoising with DADA2. Using quality score visualizations, you can choose trunc-len-f and trunc-len-r (note: sequences < trunc-len in length are discarded!)
# The drop-off for the forward reads was not so bad, but there is a significant drop-off in quality for the reverse reads, so let's trim 10bp
#For 20 million reads, this took 3.5hours, 24GB RAM @ 32CPU



#paired end denoise
if [ $FORWARDONLY = "false"]; then
    qiime dada2 denoise-paired --p-n-threads "$THREADS" --verbose --i-demultiplexed-seqs "$READQZA" --o-table "$FEATURETABLE" --o-representative-sequences "$REPSEQ" --p-trim-left-f "$TRIMLF" --p-trim-left-r "$TRIMLR" --p-trunc-len-f "$TRIMRF" --p-trunc-len-r "$TRIMRR"
#--o-denoising-stats "$DADASTAT"
else
#forward only denoise
    qiime dada2 denoise-single --p-n-threads "$THREADS" --verbose --i-demultiplexed-seqs "$READQZA" --o-table "$FEATURETABLE" --o-representative-sequences "$REPSEQ" --p-trim-left "$TRIMLF" --p-trunc-len "$TRIMRF"
#--o-denoising-stats "$DADASTAT"
fi

#This visualization shows us the sequences/sample spread
qiime feature-table summarize --i-table "$FEATURETABLE".qza --o-visualization "$FEATURETABLE".qzv

#Filter out sequences with few samples, use the number of reads from the lowest sample from the $resultname.qzv file
qiime feature-table filter-samples --i-table "$FEATURETABLE".qza --p-min-frequency "$DEPTH" --o-filtered-table "$FILTEREDFEATURETABLE"

qiime feature-table summarize --i-table "$FILTEREDFEATURETABLE".qza --o-visualization $FILTEREDFEATURETABLE.qzv

#sklearn:
echo "sklearn" > $checkPointFile
#restarted from last checkpoint, needs to remove files
if [ $restartedFromCheckpoint = "true"]; then
    echo "removing residual files from previous run"
    rm -rf "$TAXONOMYSKLEARN".qza
    rm -rf "$TAXONOMYSKLEARN".qzv
    rm -rf "$FILTEREDTAXONOMYSKLEARN".qzv
fi

echo "assigning taxonomy"
#QIIME group has a 515f 806r 99% pre-clustered GreenGenes database

#Classify against it with Naive Bayes
#for 20million reads, this took 30minutes, 28GB ram @ 32CPU

date
qiime feature-classifier classify-sklearn --i-classifier "$sklearnSilvaPath" --i-reads "$REPSEQ".qza --o-classification "$TAXONOMYSKLEARN" --p-n-jobs "$THREADS" --verbose

qiime taxa barplot --i-table "$FEATURETABLE".qza --i-taxonomy "$TAXONOMYSKLEARN".qza --m-metadata-file "$metadataPath" --o-visualization "$TAXONOMYSKLEARN".qzv
qiime taxa barplot --i-table "$FILTEREDFEATURETABLE".qza --i-taxonomy "$TAXONOMYSKLEARN".qza --m-metadata-file "$metadataPath" --o-visualization "$FILTEREDTAXONOMYSKLEARN".qzv

#because we are running two different versions of qiime, need to change environment
source deactivate

#vsearch:
echo "vsearch" > $checkPointFile
#restarted from last checkpoint, needs to remove files
if [ $restartedFromCheckpoint = "true"]; then
    echo "removing residual files from previous run"
    rm -rf "$TAXONOMYVSEARCH".qza
    rm -rf "$TAXONOMYVSEARCH".qzv
    rm -rf "$FILTEREDTAXONOMYVSEARCH".qzv
fi

source activate qiime2-2018.8

date
qiime feature-classifier classify-consensus-vsearch --i-query "$REPSEQ".qza --i-reference-reads "$vsearchSilvaPath" --i-reference-taxonomy "$vsearchSilvaTaxonomyPath" --p-perc-identity 0.99 --o-classification "$TAXONOMYVSEARCH".qza --p-threads "$THREADS" --verbose


#Taxa bar plots
date

qiime taxa barplot --i-table "$FEATURETABLE".qza --i-taxonomy "$TAXONOMYVSEARCH".qza --m-metadata-file "$metadataPath" --o-visualization "$TAXONOMYVSEARCH".qzv
qiime taxa barplot --i-table "$FILTEREDFEATURETABLE".qza --i-taxonomy "$TAXONOMYVSEARCH".qza --m-metadata-file "$metadataPath" --o-visualization "$FILTEREDTAXONOMYVSEARCH".qzv

source deactivate

#tree:
echo "tree" > $checkPointFile
#restarted from last checkpoint, needs to remove files
if [ $restartedFromCheckpoint = "true"]; then
    echo "removing residual files from previous run"
    rm -rf "$RESULTNAME".alignedrepseq.qza
    rm -rf "$RESULTNAME".maskedalignedrepseq.qza
    rm -rf "$RESULTNAME".unrooted_tree.qza
    rm -rf "$ROOTEDTREE".qza

fi

echo "constructing tree"
date

source activate qiime2-2018.2
#Steps for generating a phylogenetic tree

qiime alignment mafft --i-sequences "$REPSEQ".qza --o-alignment "$RESULTNAME".alignedrepseq
qiime alignment mask --i-alignment "$RESULTNAME".alignedrepseq.qza --o-masked-alignment "$RESULTNAME".maskedalignedrepseq
qiime phylogeny fasttree --i-alignment "$RESULTNAME".maskedalignedrepseq.qza --o-tree "$RESULTNAME".unrooted_tree
qiime phylogeny midpoint-root --i-tree "$RESULTNAME".unrooted_tree.qza --o-rooted-tree "$ROOTEDTREE"

#diversity:
echo "diversity" > $checkPointFile
#restarted from last checkpoint, needs to remove files
if [ $restartedFromCheckpoint = "true"]; then
    echo "removing residual files from previous run"
    rm -rf "$ALPHADIV"
    rm -rf "$FILTEREDALPHADIV"
fi

echo "calculating diversity"
date
#Generate alpha/beta diversity measures at 41000 sequences/sample
#Also generates PCoA plots automatically
#non filtered first


qiime diversity core-metrics-phylogenetic --i-phylogeny "$ROOTEDTREE.qza" --i-table "$FEATURETABLE".qza --p-sampling-depth "$DIVERSITYDEPTH" --m-metadata-file "$metadataPath" --output-dir "$ALPHADIV"
#Test for between-group differences
qiime diversity alpha-group-significance --i-alpha-diversity "$ALPHADIV"/faith_pd_vector.qza --m-metadata-file "$metadataPath" --o-visualization "$ALPHADIV"/"$ALPHADIV".alpha_PD_significance.qzv
qiime diversity alpha-group-significance --i-alpha-diversity "$ALPHADIV"/shannon_vector.qza --m-metadata-file "$metadataPath" --o-visualization "$ALPHADIV"/"$ALPHADIV".alpha_shannon_significance.qzv
#qiime diversity beta-group-significance --i-distance-matrix "$ALPHADIV"/bray_curtis_distance_matrix.qza --m-metadata-file "$METADATA" --m-metadata-column bee_type --o-visualization beta_bray_beetype_significance
#Alpha rarefaction curves show taxon accumulation as a function of sequence depth
qiime diversity alpha-rarefaction --i-phylogeny "$ROOTEDTREE".qza --i-table "$FEATURETABLE".qza --p-max-depth "$DIVERSITYDEPTH" --o-visualization "$ALPHADIV"/"$ALPHADIV".alpha_rarefaction.qzv --m-metadata-file "$metadataPath"

#then filtered

qiime diversity core-metrics-phylogenetic --i-phylogeny "$ROOTEDTREE.qza" --i-table "$FILTEREDFEATURETABLE".qza --p-sampling-depth "$DIVERSITYDEPTH" --m-metadata-file "$metadataPath" --output-dir "$FILTEREDALPHADIV"
#Test for between-group differences
qiime diversity alpha-group-significance --i-alpha-diversity "$FILTEREDALPHADIV"/faith_pd_vector.qza --m-metadata-file "$metadataPath" --o-visualization "$FILTEREDALPHADIV"/"$FILTEREDALPHADIV".alpha_PD_significance.qzv
qiime diversity alpha-group-significance --i-alpha-diversity "$FILTEREDALPHADIV"/shannon_vector.qza --m-metadata-file "$metadataPath" --o-visualization "$FILTEREDALPHADIV"/"$FILTEREDALPHADIV".alpha_shannon_significance.qzv
#qiime diversity beta-group-significance --i-distance-matrix "$FILTEREDALPHADIV"/bray_curtis_distance_matrix.qza --m-metadata-file "$METADATA" --m-metadata-column bee_type --o-visualization beta_bray_beetype_significance
#Alpha rarefaction curves show taxon accumulation as a function of sequence depth
qiime diversity alpha-rarefaction --i-phylogeny "$ROOTEDTREE".qza --i-table "$FILTEREDFEATURETABLE".qza --p-max-depth "$DIVERSITYDEPTH" --o-visualization "$FILTEREDALPHADIV"/"$FILTEREDALPHADIV".alpha_rarefaction.qzv --m-metadata-file "$metadataPath"

#stamp:
echo "stamp" > $checkPointFile
#restarted from last checkpoint, needs to remove files
if [ $restartedFromCheckpoint = "true"]; then
    echo "removing residual files from previous run"
    rm -rf stamptables.sklearn
    rm -rf stamptables."$DEPTH".vsearch
    rm -rf stamptables.sklearn
    rm -rf stamptables."$DEPTH".vsearch
fi

#get a tsv file to export to STAMP
#start with non filtered
mkdir stamptables.sklearn
cd stamptables.sklearn
LEVEL=3
for LEVEL in 3 4 5 6 7;
do
	STAMPNAME="$RESULTNAME".level"$LEVEL".sklearn.taxontable
	qiime taxa collapse --i-table ../"$FEATURETABLE".qza --i-taxonomy ../"$TAXONOMYSKLEARN".qza --p-level "$LEVEL" --o-collapsed-table "$STAMPNAME";
	qiime tools export --output-dir "$STAMPNAME".biom "$STAMPNAME".qza;
	biom convert -i "$STAMPNAME".biom/feature-table.biom -o "$STAMPNAME".tsv --to-tsv;
	tail -n +2 "$STAMPNAME".tsv > "$STAMPNAME".tsv
done;
cd ../

#filtered
mkdir stamptables."$DEPTH".sklearn
cd stamptables."$DEPTH".sklearn
LEVEL=3
for LEVEL in 3 4 5 6 7;
do
	STAMPNAME="$RESULTNAME"."$DEPTH".level"$LEVEL".sklearn.taxontable
	qiime taxa collapse --i-table ../"$FILTEREDFEATURETABLE".qza --i-taxonomy ../"$TAXONOMYSKLEARN".qza --p-level "$LEVEL" --o-collapsed-table "$STAMPNAME";
	qiime tools export --output-dir "$STAMPNAME".biom "$STAMPNAME".qza	;
	biom convert -i "$STAMPNAME".biom/feature-table.biom -o "$STAMPNAME".tsv --to-tsv;
	tail -n +2 "$STAMPNAME".tsv > "$STAMPNAME".tsv
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
	qiime taxa collapse --i-table ../"$FEATURETABLE".qza --i-taxonomy ../"$TAXONOMYVSEARCH".qza --p-level "$LEVEL" --o-collapsed-table "$STAMPNAME";
	qiime tools export --output-path "$STAMPNAME".biom "$STAMPNAME".qza	;
	biom convert -i "$STAMPNAME".biom/feature-table.biom -o "$STAMPNAME".tsv --to-tsv;
	tail -n +2 "$STAMPNAME".tsv	> "$STAMPNAME".tsv
done;
cd ../

#filtered
mkdir stamptables."$DEPTH".vsearch
cd stamptables."$DEPTH".vsearch
LEVEL=3
for LEVEL in 3 4 5 6 7;
do
	STAMPNAME="$RESULTNAME"."$DEPTH".level"$LEVEL".vsearch.taxontable
	qiime taxa collapse --i-table ../"$FILTEREDFEATURETABLE".qza --i-taxonomy ../"$TAXONOMYVSEARCH".qza --p-level "$LEVEL" --o-collapsed-table "$STAMPNAME";
	qiime tools export --output-path "$STAMPNAME".biom "$STAMPNAME".qza	;
	biom convert -i "$STAMPNAME".biom/feature-table.biom -o "$STAMPNAME".tsv --to-tsv;
	tail -n +2 "$STAMPNAME".tsv	> "$STAMPNAME".tsv
done;
cd ../
date
source deactivate

#copy phyloseq related objects 
#required inputs:
#1. metadata in utf8
#2. feature table
#3. taxonomy
#4. rooted tree

#source load r3.5.1

#clean:
echo "clean" > $checkPointFile

#mv the useful files to scp folder and zip it
echo "finishing up"

if [ $restartedFromCheckpoint = "true"]; then
    echo "removing residual files from previous run"
    rm -rf scp.tar.gz
fi

date
mkdir scp
mv *.qzv scp
mv "$ALPHADIV"/*.qzv scp
mv "$FILTEREDALPHADIV"/*.qzv scp
mv *.featuretable*.qza scp
mv *.taxonomy*.qza scp
mv *.rootedTree*.qza scp
mv *.repseq*.qza scp
mv stamptables* scp

tar -zcvf scp.tar.gz scp/*

echo "finish" > $checkPointFile
echo "done"
