
# DADA2 -> PICRUSt via ASR (v2.1)

## Readme
Updates from (v2)
- Added KEGG metadata for categorize_by_function.py
Updates from (v1):
- Validation against shotgun metagenome data available
	- DADA2-PICRUSt performs as well as VSEARCH-PICRUSt (see figures below)
- NSTI values now calculated
- No longer need to move new 16S\_ or ko\_ precaculated files into PICRUSt data dir
- Parallel options listed for Part 2
- Alignment filtering added
- Predictions made from the full official ko_13_5 and 16S_13_5 precalculated files, which maximize prediction accuracy (for the time being)

Inputs:
http://greengenes.secondgenome.com/downloads/database/13_5
- greengenes 13.5 fasta (gg_13_5.fasta, see link)
- PICRUSt precalculated files (both ko_13_5_ and 16S_13_5_precalculated.tab)

Outputs:
- fasta containing 1) ASR-ready 13.5 16S sequences, 2) denoised DADA2 sequences
- new 16S_13_5_precalculated.tab
- new ko_13_5_precalculated.tab
- PICRUSt-normalized .biom and -metagenome predictions .biom

Points of improvement:
- Add KEGG BRITE hierarchy pathway information (cannot use categorize_by_function.py without this)
- Add phyloseq import_genome function?

Multiple lines in part 3 of this workflow were taken from the official picrust genome prediction tutorial https://picrust.github.io/picrust/tutorials/genome_prediction.html

## Part 0: Generate gg_16S_counts.tab and gg_ko_counts.tab files
**Note:** This step needs to be performed only once. Once these files are on-hand, they may be reused in subsequent projects. These files are used in Part 3 only.
```sh
# build gg_16S_counts.tab and gg_ko_counts.tab
gunzip -c 16S_13_5_precalculated.tab.gz | sed 's/\#OTU_IDs/taxon_oid/g' > gg_16S_counts.tab
gunzip -c ko_13_5_precalculated.tab.gz | sed 's/\#OTU_IDs/GenomeID/g' > gg_ko_counts.tab
# remove empty lines in gg_ko_counts.tab
sed -i '/^\s*$/d' gg_ko_counts.tab
# save KEGG Pathways and Description metadata
gunzip -c ko_13_5_precalculated.tab.gz | tail -n 2 > kegg_meta
```

## Part 1: Start with R
1. Add study sequences to `gg_13_5.fasta`
2. Output sample count sequences for predict_metagenomes.py later

**Note:** Make sure to modify file locations and directories as you like them; this pipeline as written assumes the input files and folders (listed above) are in a folder named `genome_prediction` which is within your current directory. Also, double check that vsearch is in your path. If not, just swap "vsearch" below with an absolute path to the vsearch binary.

```R
# Dependencies: ShortRead & biom
library(ShortRead)
library(biom) # note: use Joey's biom latest dev version; library(devtools); install_github("joey711/biom");
# 1) Make study db
# grab study seqs
load(file = "genome_prediction/seqtab.nochim.robj")
seqs_study <- colnames(seqtab.nochim)
ids_study <- paste("study", 1:ncol(seqtab.nochim), sep = "_")
# merge db and study seqs
db_out <- data.frame(ids=ids_study,seqs=seqs_study,count=colSums(seqtab.nochim))
fasta <- ShortRead(sread = DNAStringSet(db_out$seqs), id = BStringSet(db_out$ids))
# write study fasta for filtering
writeFasta(fasta, file = "genome_prediction/gg_13_5_study_db.fasta.pre")
# filter sequences that diverge from gg_13_5 by 97%
# depending on how well greengenes covers your study sequences, consider reducing 97% to 70 or 50%
system('vsearch --usearch_global genome_prediction/gg_13_5_study_db.fasta.pre --db genome_prediction/gg_13_5.fasta --matched genome_prediction/gg_13_5_study_db.fasta --id 0.97')
id_filtered <- as.character(id(readFasta("genome_prediction/gg_13_5_study_db.fasta")))
db_out_filt <- db_out[db_out$ids%in%id_filtered,]
seqtab_biom <- t(seqtab.nochim)
# 2) output seq variant count data as biom;
# subset seqtab and output sample count biom
seqtab_biom <- seqtab_biom[rownames(seqtab_biom)%in%db_out_filt$seqs,]
rownames(seqtab_biom) <- db_out_filt[db_out_filt$seqs%in%rownames(seqtab_biom),"ids"]
biom_object <- biom::make_biom(data = seqtab_biom)
biom::write_biom(biom_object, biom_file = "sample_counts.biom")
# create final study db
system('cat genome_prediction/gg_13_5.fasta >> genome_prediction/gg_13_5_study_db.fasta')
```
## Part 2: Align seqs and build tree
**Note:** I chose pynast to align and fasttree to build the tree since they are quick, but any alignment / tree building combo is fair game for optimization
```sh
# align w/ pynast using QIIME scripts; the included options lessen alignment restrictions to prevent alignment failure
# minimum sequence length set by -e
# alignment runtime greatly reduced by parallelization: parallel_align_seqs_pynast.py -O and # of cores
parallel_align_seqs_pynast.py -e 90 -p 0.1 -i ./genome_prediction/gg_13_5_study_db.fasta -o ./genome_prediction/gg_13_5_study_db.fasta.aligned -O 45
# filter alignment with default settings; consider lane filtering by entropy using -e and a low entropy value of ~0.01-0.02
# note: FastTree and/or PICRUSt produce weird errors (segfaults) if -e filters too many lanes
filter_alignment.py -i ./genome_prediction/gg_13_5_study_db.fasta.aligned/gg_13_5_study_db_aligned.fasta -o ./genome_prediction/gg_13_5_study_db.fasta.aligned.filtered/
# build tree with fasttree; options are taken from greengenes 13_5 readme notes
# tree building runtime greatly reduced by parallelization: use FastTreeMP w/ same options instead of FastTree
FastTreeMP -nt -gamma -fastest -no2nd -spr 4 ./genome_prediction/gg_13_5_study_db.fasta.aligned.filtered/gg_13_5_study_db_aligned_pfiltered.fasta > ./genome_prediction/study_tree.tree
```
## Part 3: Create new precalculated files
**Adapted from** https://picrust.github.io/picrust/tutorials/genome_prediction.html

```sh
# format 16S copy number data
format_tree_and_trait_table.py -t ./genome_prediction/study_tree.tree -i gg_16S_counts.tab -o ./genome_prediction/format/16S/
# format kegg IMG data
format_tree_and_trait_table.py -t ./genome_prediction/study_tree.tree -i gg_ko_counts.tab -o ./genome_prediction/format/KEGG/
# perform ancestral state reconstruction
ancestral_state_reconstruction.py -i ./genome_prediction/format/16S/trait_table.tab -t ./genome_prediction/format/16S/pruned_tree.newick -o ./genome_prediction/asr/16S_asr_counts.tab -c ./genome_prediction/asr/asr_ci_16S.tab
ancestral_state_reconstruction.py -i ./genome_prediction/format/KEGG/trait_table.tab -t ./genome_prediction/format/KEGG/pruned_tree.newick -o ./genome_prediction/asr/KEGG_asr_counts.tab -c ./genome_prediction/asr/asr_ci_KEGG.tab
# collect study sequence ids for predict_traits.py -g (greatly reduces runtime)
grep "study_[0-9]*" ./genome_prediction/gg_13_5_study_db.fasta.aligned.filtered/gg_13_5_study_db_aligned_pfiltered.fasta -o | tr "\n" "," > study_ids
# predict traits
predict_traits.py -i ./genome_prediction/format/16S/trait_table.tab -t ./genome_prediction/format/16S/reference_tree.newick -r ./genome_prediction/asr/16S_asr_counts.tab -o ./genome_prediction/predict_traits/16S_precalculated.tab -a -c ./genome_prediction/asr/asr_ci_16S.tab -g "$(< study_ids)"
predict_traits.py -i ./genome_prediction/format/KEGG/trait_table.tab -t ./genome_prediction/format/KEGG/reference_tree.newick -r ./genome_prediction/asr/KEGG_asr_counts.tab -o ./genome_prediction/predict_traits/ko_precalculated.tab -a -c ./genome_prediction/asr/asr_ci_KEGG.tab -g "$(< study_ids)"
# add KEGG metadata
cat kegg_meta >> ./genome_prediction/predict_traits/ko_precalculated.tab
```
## Part 4: Run PICRUSt (finally!)
```sh
# run PICRUSt
normalize_by_copy_number.py -i sample_counts.biom -o norm_counts.biom -c ./genome_prediction/predict_traits/16S_precalculated.tab
predict_metagenomes.py -i norm_counts.biom -o meta_counts_asr.biom -c ./genome_prediction/predict_traits/ko_precalculated.tab
# optional: agglomerate counts by KEGG pathway level
categorize_by_function.py -i meta_counts_asr.biom -o cat_2_counts_asr.biom -c KEGG_Pathways -l 2
```

Happy to discuss any part of this / incorporate suggested modifications.

## DADA2 -> PICRUSt validation

### Approach & Workflow
We tested the accuracy of two experimental pipelines against that of the original PICRUSt pipeline. We used two paired-16S-shotgun metagenome datasets available on NCBI SRA (see links below to Jovel et. al. 2016 and HMP). Both the Jovel dataset and the HMP dataset downloaded for this test were generated by Illumina sequencers. The experimental pipelines are labeled "DADA2...ASR" and "DADA2...kmer." The original pipeline is labeled "VSEARCH...pick_closed." In DADA2...ASR we applied the genome prediction tutorial on DADA2-clustered reads (see above). In DADA2...kmer we used RDP to assign a best-hit to individual sequences in the greengenes database and filtered for hits with > 80% bootstrapped confidence. This effectively assigns a greengenes OTU ID to DADA2-clustered reads for subsequent PICRUSt analysis. We produced predicted KEGG gene counts from all three pipelines using the "ace_pic," default ASR method and an entropy filter of 0.02 (-e 0.02 in Part 2). We compared predicted KEGG counts to HUMAnN2 translated search KEGG gene count profiles from paired shotgun metagenomic samples. We pre-processed the shotgun reads with KneadData and combined forward and reverse reads before running HUMAnN2. We also used the legacy KEGG protein database as provided by the Huttenhower lab for compatibility with the current version of PICRUSt. We found that both experimental pipelines performed essentially as well as the original pipeline using Spearman rho as a measure of accuracy between normalized 16S and shotgun metagenome KEGG gene counts (Figure 1 & 2). These fit pretty closely with those published in figure 3 of the original PICRUSt paper (http://www.nature.com/nbt/journal/v31/n9/full/nbt.2676.html).
### Jovel et. al. 2016; https://www.ncbi.nlm.nih.gov/pubmed/27148170
![Alt text](https://github.com/vmaffei/dada2_to_picrust/blob/master/jovel_asr_kmer_closed.png?raw=true "Jovel et. al. 2016 Results")
### HMP; http://www.hmpdacc.org/
![Alt text](https://github.com/vmaffei/dada2_to_picrust/blob/master/hmp_asr_kmer_closed_all.png?raw=true "HMP Results")
