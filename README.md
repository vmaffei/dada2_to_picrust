
# DADA2 -> picrust via ASR

## Readme
Inputs:
http://greengenes.secondgenome.com/downloads/database/13_5
- greengenes 13.5 fasta (gg_13_5.fasta, see link)
- IMG ids to gg ids mapping file (gg_13_5_img.txt, see link)
- picrust tutorial files (found in picrust/tutorials/picrust_starting_files.zip)

Outputs:
- fasta containing 1) ASR-ready 13.5 16S sequences, 2) denoised DADA2 sequences
- new 16S_13_5_precalculated.tab
- new ko_13_5_precalculated.tab
- PICRUSt-normalized .biom and -metagenome predictions .biom

Points of improvement:
- 16S sequence alignment (chosen aligner / options, are better available?)
- Add KEGG BRITE hierarchy pathway information (cannot use categorize_by_function.py without this; currently waiting to hear back from KEGG about their API download policy)
- Calculate new genome prediction confidence intervals (NSTI's)

Multiple lines in part 3 of this workflow were taken from the official picrust genome prediction tutorial https://picrust.github.io/picrust/tutorials/genome_prediction.html


## Part 1: Start with R
1. Make study sequence db; subset `gg_13_5.fasta` for sequences with matching IMG kegg data (subsetting to the essentials greatly speeds up run-time)
2. Output sample count sequences for predict_metagenomes.py later

**Note:** Make sure to modify file locations and directories as you like them; this pipeline as written assumes the input files and folders (listed above) are in a folder named `genome_prediction` which is within your current directory.

```R
# Dependencies
library(ShortRead)
library(biom) # use Joey's biom latest dev version; library(devtools); install_github("joey711/biom")
library(dada2)
# 1) Make study db
# read in gg 13_5 db
ref <- readFasta("genome_prediction/gg_13_5.fasta")
# read in IMG to ko
img_ko <- read.table(file = "genome_prediction/picrust_starting_files/IMG_ko_counts.tab", sep = "\t", row.names = 1, header = TRUE, comment.char = '')
# read in gg id to IMG
id_img <- read.table(file = "genome_prediction/gg_13_5_img.txt", sep = "\t", header = TRUE, comment.char = '')
# intersect gg id and ko list
img_ko_sub <- img_ko[rownames(img_ko) %in% id_img$img_genome_id, ]
id_img_sub <- id_img[id_img$img_genome_id %in% rownames(img_ko_sub), ]
# pull db ids
ids <- id(ref)
# subset db for gg w/ ko
ref_sub <- ref[as.character(ids) %in% id_img_sub$X.gg_id, ]
# build study db
ids_db <- as.character(id(ref_sub))
seqs_db <- as.character(sread(ref_sub))
# grab study seqs
load(file = "seqtab.nochim.robj")
seqs_study <- colnames(seqtab.nochim)
ids_study <- paste("study", 1:ncol(seqtab.nochim), sep = "_")
# merge db and study seqs
ids_out <- c(ids_db, ids_study)
seqs_out <- c(seqs_db, seqs_study)
fasta <- ShortRead(sread = DNAStringSet(seqs_out), id = BStringSet(ids_out))
writeFasta(fasta, file = "genome_prediction/gg_13_5_study_db.fasta")
# (optional) output gg id-to-IMG mapping info; modify this to filter later for predicted traits (predict_traits.py -l )
write.table(as.character(id(fasta)), file = "genome_prediction/traits_sample_filter.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
# 2) output seq variant count data as biom;
# note: use Joey's biom latest dev version; library(devtools); install_github("joey711/biom")
seqtab_biom <- t(seqtab.nochim)
rownames(seqtab_biom) <- ids_study
biom_object <- biom::make_biom(data = seqtab_biom)
biom::write_biom(biom_object, biom_file = "genome_prediction/sample_counts.biom")
```
## Part 2: Align seqs and build tree
**Note:** I chose pynast to align and fasttree to build the tree since they are quick, but any alignment / tree building combo is fair game for optimization
```sh
# align w/ pynast using QIIME scripts; the included options lessen alignment restrictions to prevent alignment failure
align_seqs.py -e 90 -p 0.1 -i ./genome_prediction/gg_13_5_study_db.fasta -o ./genome_prediction/gg_13_5_study_db.fasta.aligned
# build tree with fasttree; options are taken from greengenes 13_5 readme notes
FastTree -nt -gamma -fastest -no2nd -spr 4 ./genome_prediction/gg_13_5_study_db.fasta.aligned > ./genome_prediction/study_tree.tree
```
## Part 3: Create new precalculated files
**Adapted from** https://picrust.github.io/picrust/tutorials/genome_prediction.html

```sh
# format 16S copy number data
format_tree_and_trait_table.py -t ./genome_prediction/study_tree.tree -i ./genome_prediction/picrust_starting_files/IMG_16S_counts.tab -m ./genome_prediction/gg_13_5_img.txt -o ./genome_prediction/format/16S/
# format kegg IMG data
format_tree_and_trait_table.py -t ./genome_prediction/study_tree.tree -i ./genome_prediction/picrust_starting_files/IMG_ko_counts.tab -m ./genome_prediction/gg_13_5_img.txt -o ./genome_prediction/format/KEGG/
# perform ancestral state reconstruction
ancestral_state_reconstruction.py -i ./genome_prediction/format/16S/trait_table.tab -t ./genome_prediction/format/16S/pruned_tree.newick -o ./genome_prediction/asr/16S_asr_counts.tab
ancestral_state_reconstruction.py -i ./genome_prediction/format/KEGG/trait_table.tab -t ./genome_prediction/format/KEGG/pruned_tree.newick -o ./genome_prediction/asr/KEGG_asr_counts.tab
# predict traits
predict_traits.py -i ./genome_prediction/format/16S/trait_table.tab -t ./genome_prediction/format/16S/reference_tree.newick -r ./genome_prediction/asr/16S_asr_counts.tab -o ./genome_prediction/predict_traits/16S_precalculated.tab
predict_traits.py -i ./genome_prediction/format/KEGG/trait_table.tab -t ./genome_prediction/format/KEGG/reference_tree.newick -r ./genome_prediction/asr/KEGG_asr_counts.tab -o ./genome_prediction/predict_traits/ko_precalculated.tab
# copy/move precalculated file to replace PICRUSt precalculated data files
cp ./genome_prediction/predict_traits/16S_precalculated.tab ./genome_prediction/predict_traits/16S_13_5_precalculated.tab
cp ./genome_prediction/predict_traits/ko_precalculated.tab ./genome_prediction/predict_traits/ko_13_5_precalculated.tab
# gzip files
gzip ./genome_prediction/predict_traits/16S_13_5_precalculated.tab
gzip ./genome_prediction/predict_traits/ko_13_5_precalculated.tab
```
**Important note from last step:** ensure the final gzipped files are named exactly as above so they may successfully drop in and replace the official picrust precalculated files. On that note be sure to backup your original picrust precalculated files in a safe/convenient place. If you replace them on accident, no worries just re-download https://picrust.github.io/picrust/picrust_precalculated_files.html.
### From here go ahead and move or copy the new .tab.gz files to your PICRUSt precalculated files directory.
```sh
cp ./genome_prediction/predict_traits/16S_13_5_precalculated.tab.gz picrust/data/
cp ./genome_prediction/predict_traits/ko_13_5_precalculated.tab.gz picrust/data/
```
## Part 4: Run PICRUSt (finally!)
```sh
# run PICRUSt
normalize_by_copy_number.py -i ./genome_prediction/sample_counts.biom -o ./picrust/norm_counts.biom
predict_metagenomes.py -i ./picrust/norm_counts.biom -o ./picrust/meta_counts.biom
```
**Final note:** the last 2 key missing pieces are:
1. KEGG BRITE pathway information for each KEGG ID
2. Prediction's confidence intervals

Happy to discuss any part of this / incorporate suggested modifications.
