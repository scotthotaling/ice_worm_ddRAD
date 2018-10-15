######################################
Iceworm ddRAD analysis notebok
Written by Scott Hotaling
Scott.Hotaling@wsu.edu
Last updated = 15 October 2018
######################################

## 1. Demultiplexing the data
	> Unpack the data:
		tar -xzvf Iceworm_RawData_201829.tgz
	> Two libraries:
		1 = ATCACG
		2 = AGTCAA		
	> Move each library to its own folder 
	> Prepare demultiplex map file for each library (using process_radtags)

/// Demultiplex with "process_radtags" [script = demultiplex_ice_worm.sh]
	# -r = paired-end reads are interleaved in single files
	# -c = “clean data” – removes any read with an uncalled base
	# -q = reads with quality score below 10 (phred33) in a sliding window (15% length of the read) 
	# are discarded 

process_radtags \
        -f ~/1_ATCACG_L006_R1_001.fastq \
        -o ~/Lib1_Demultiplex \
        -b ~/L1-Iceworm-Demultiplex.txt \
        --renz_1 nlaIII --renz_2 ecoRI -r -c -q --retain_header &

process_radtags \
        -f ~/2_AGTCAA_L006_R1_001.fastq \
        -o ~/Lib2_Demultiplex \
        -b ~/L2-Iceworm-Demultiplex.txt \
        --renz_1 nlaIII --renz_2 ecoRI -r -c -q --retain_header

// L1-Iceworm-Demultiplex.txt
ATTAT	MS1
CACCA	MS2
CCTCG	MS3
CTTGA	MS4
GGATA	MS5
GGTGT	MS6
TATGT	MS7
TGCTT	MS8
AACTGG	MS9
ACAACT	MS10
ATAGAT	MS11
CAGATA	MS12
GAAGTG	MS13
GGCTTA	MS14
TCTTGG	MS15
TCACTG	MS16
ACCAGGA	MS17
CCACTCA	MS18
CCGAACA	MS19
CTAAGCA	MS20
CTCGCGG	MS21
GCGTCCT	MS22
GGAACGA	MS23
TAGCCAA	MS24
ACTGCGAT	MS25
ATGAGCAA	MS26
GCCTACCT	MS27
TAGCGGAT	MS28
TGACGCCA	MS29
ACGGTACT	MS30
AAGACGCT	MS31
TCAGAGAT	MS32

// L2-Iceworm-Demultiplex.txt
ATTAT	MS33
CACCA	MS34
CCTCG	MS35
CTTGA	MS37
GGATA	MS38
GGTGT	MS39
TATGT	MS40
TGCTT	MS41
AACTGG	MS42
ACAACT	MS43
ATAGAT	MS44
CAGATA	MS45
GAAGTG	MS46
GGCTTA	MS47
TCTTGG	MS48
TCACTG	MS49
ACCAGGA	MS50
CCACTCA	MS52
CCGAACA	MS53
CTAAGCA	MS54
CTCGCGG	MS55
GCGTCCT	MS56
GGAACGA	MS57
TAGCCAA	MS58
ACTGCGAT	MS59
ATGAGCAA	MS61
GCCTACCT	MS62
TAGCGGAT	MS63
TGACGCCA	MS64
ACGGTACT	EA1
AAGACGCT	MO1
TCAGAGAT	MO3

> Run it: sbatch -p Long -N 1 -n 16 demultiplex_ice_worm.sh
	
> Combine demultiplexed libraries into one folder for downstream analysis		
		>> mkdir Demultiplexed
		>> [in Lib1_Demultiplex] cp *.fq ../Demultiplexed/
		>> [in Lib2_Demultiplex] cp *.fq ../Demultiplexed/


## 2. Running STACKS pipeline to call SNPs ("denovo_map")

/// Demulti-plex with "process_radtags" (script = denovo_map.sh)

denovo_map.pl -m 10 -n 2 -T 16 -b 1 -S \
-D "de novo map ice worms" \
-e "path to stacks" \
-o ~/b1_output \
-s ~/MS1.fq -s ~/MS2.fq -s ~/MS3.fq -s ~/MS4.fq -s ~/MS5.fq -s ~/MS8.fq -s ~/MS9.fq -s ~/MS10.fq \
-s ~/MS11.fq -s ~/MS12.fq -s ~/MS13.fq -s ~/MS14.fq -s ~/MS15.fq -s ~/MS16.fq -s ~/MS17.fq \
-s ~/MS18.fq -s ~/MS19.fq -s ~/MS20.fq -s ~/MS21.fq -s ~/MS22.fq -s ~/MS23.fq -s ~/MS24.fq \
-s ~/MS25.fq -s ~/MS26.fq -s ~/MS27.fq -s ~/MS28.fq -s ~/MS29.fq -s ~/MS30.fq -s ~/MS31.fq \
-s ~/MS32.fq -s ~/MS33.fq -s ~/MS34.fq -s ~/MS35.fq -s ~/MS37.fq -s ~/MS38.fq -s ~/MS39.fq \
-s ~/MS40.fq -s ~/MS41.fq -s ~/MS42.fq -s ~/MS43.fq -s ~/MS44.fq -s ~/MS45.fq -s ~/MS46.fq \
-s ~/MS47.fq -s ~/MS48.fq -s ~/MS49.fq -s ~/MS50.fq -s ~/MS52.fq -s ~/MS53.fq -s ~/MS54.fq \
-s ~/MS55.fq -s ~/MS56.fq -s ~/MS57.fq -s ~/MS58.fq -s ~/MS59.fq -s ~/MS61.fq -s ~/MS62.fq \
-s ~/MS63.fq -s ~/MS64.fq

> sbatch -p Long -N 1 -n 16 denovo_map.sh
	

## 3. Running STACKS "populations" to output SNPs.

	> Example map file. [ice_worm_map-1.txt]
		
MS1	1
MS2	4
MS3	5
MS4	6
MS5	3
MS8	1
MS9	4
MS10	4
MS11	4
MS12	4
MS13	7
MS14	7
MS15	7
MS16	8
MS17	5
MS18	5
MS19	9
MS20	9
MS21	2
MS22	2
MS23	4
MS24	4
MS25	4
MS26	4
MS27	4
MS28	4
MS29	4
MS30	7
MS31	7
MS32	8
MS33	8
MS34	8
MS35	8
MS37	5
MS38	6
MS39	6
MS40	6
MS41	6
MS42	9
MS43	9
MS44	9
MS45	2
MS46	3
MS47	3
MS48	3
MS49	1
MS50	1
MS52	4
MS53	4
MS54	4
MS55	7
MS56	7
MS57	7
MS58	8
MS59	5
MS61	5
MS62	6
MS63	6
MS64	6
	
	> Run populations.
	
~/stacks-1.46/populations -P b1_output -M iceworm_map_1.txt -O <output path> -b 1 -p 5 --fasta --vcf --plink --fstats
		
		
## 4. Additional filtering and data manipulation in PLINK
	> Flags of interest
		--geno = filters SNPs for freq of missing data
		--mind = filters  individuals for freq of missing data
		--hardy = filters SNPs for HWE
		--maf = filters SNPs for MAF
		--recode = outputs new PED/MAP files with filters applied
	
	> Filtering for HWE, MAF, and missing data:
		./plink --file 1.5.plink --missing --hardy --maf 0.025 --geno 0.5 --recode --noweb


## 5. DAPC
	> Make input files via PLINK:
		./plink --file 1.5.1 --recodeA --noweb
		Re-order/pull IDs in Text Wrangler+Excel
		
// Script for DAPC analyses

library(adegenet)
library(ape)
rm(list=ls())

### set the working directory
setwd("path to input files")

### read in the data
	snp_data <- read.PLINK("1.5.1-oneSNP.raw")
	snp_data

### assign geographic groups to each sample
geog_grps <- c(1,1,1,1,2,2,2,2,3,3,3,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,8,8,8,8,8,8,9,9,9,9,9)
snp_geogr_grps <- as.factor(geog_grps)

### select number of clusters (answer the prompt with K for the lowest BIC)
clust_snp_data <- find.clusters(snp_data, n.pca=50, max.n.clust=10)

### assign clusters
snp_clusters <- as.factor(clust_snp_data$grp)
snp_clusters

### choosing PCs for DAPC
no_pcs <- dapc(snp_data, n.pc = 18, n.da = 5)
a_score <- optim.a.score(no_pcs)

### final DAPC - # PC from above, n.da = number of groups-1
final.dapc.ice <- dapc(snp_data, clust_snp_data$grp, n.pca = 4, n.da = 5)
scatter(final.dapc.ice, scree.pca = TRUE)

### Visualize as populations as the six assigned clusters
compoplot(final.dapc.ice, only.grp=NULL, subset=NULL, new.pred=NULL, col=c("red", "blue", "forest green", "orange", "purple", "gray"), lab=NULL, legend=TRUE, txt.leg=NULL, ncol=4, posi=NULL, cex.axis=.75, cleg=.5, bg=transp("white"))	
	
	
	
## 6. Admixture
	> Make a folder with Admixture executable in it

	> Re-code PED to Admixture format with PLINK
		./plink --file 1.5.1 --recode12 --out 1.5.1.Admix --noweb
			> Copy the output to the Admixture folder with the executable
			
	> Run one replicate each of K=1-12 with current time as the random seed	
		for K in 1 2 3 4 5 6 7 8 9 10 11 12; do ./admixture -s time 1.5.1.oneSNP.Admix.ped $K --cv >> run??.log; done
			
	> Repeat 25 times and re-name the output each time so it does not get overwritten
	
mv 1.5.1.oneSNP.Admix.1.P run25.Admix.k1.P \
mv 1.5.1.oneSNP.Admix.2.P run25.Admix.k2.P \
mv 1.5.1.oneSNP.Admix.3.P run25.Admix.k3.P \
mv 1.5.1.oneSNP.Admix.4.P run25.Admix.k4.P \
mv 1.5.1.oneSNP.Admix.5.P run25.Admix.k5.P \
mv 1.5.1.oneSNP.Admix.6.P run25.Admix.k6.P \
mv 1.5.1.oneSNP.Admix.7.P run25.Admix.k7.P \
mv 1.5.1.oneSNP.Admix.8.P run25.Admix.k8.P \
mv 1.5.1.oneSNP.Admix.9.P run25.Admix.k9.P \
mv 1.5.1.oneSNP.Admix.10.P run25.Admix.k10.P \
mv 1.5.1.oneSNP.Admix.11.P run25.Admix.k11.P \
mv 1.5.1.oneSNP.Admix.12.P run25.Admix.k12.P \
mv 1.5.1.oneSNP.Admix.1.Q run25.Admix.k1.Q \
mv 1.5.1.oneSNP.Admix.2.Q run25.Admix.k2.Q \
mv 1.5.1.oneSNP.Admix.3.Q run25.Admix.k3.Q \
mv 1.5.1.oneSNP.Admix.4.Q run25.Admix.k4.Q \
mv 1.5.1.oneSNP.Admix.5.Q run25.Admix.k5.Q \
mv 1.5.1.oneSNP.Admix.6.Q run25.Admix.k6.Q \
mv 1.5.1.oneSNP.Admix.7.Q run25.Admix.k7.Q \
mv 1.5.1.oneSNP.Admix.8.Q run25.Admix.k8.Q \
mv 1.5.1.oneSNP.Admix.9.Q run25.Admix.k9.Q \
mv 1.5.1.oneSNP.Admix.10.Q run25.Admix.k10.Q \
mv 1.5.1.oneSNP.Admix.11.Q run25.Admix.k11.Q \
mv 1.5.1.oneSNP.Admix.12.Q run25.Admix.k12.Q

	> Grab CV errors across runs so comparisons among runs can be made (I used Excel)

		e.g. grep # "CV error" run25.log		
		
	> Plot results in R
	
		setwd(#working directory)
		tbl=read.table(#run of interest.Q)
		barplot(t(as.matrix(tbl)), col=rainbow(# of clusters, K), xlab="Individual #", ylab="Ancestry", border=NA)	 


## 7. SVDQuartets
	> Get input data
		populations -P b2_output -M SVD_map.txt -O output -b 3 -W whitelist --write_random_snp --phylip_var

	> Convert Phylip to Nexus
		# Example input file = example-data-snakes.nex (obtained from data_snakes.phy -- Laura Kubatko data)
		
	> Site: https://hcv.lanl.gov/content/sequence/FORMAT_CONVERSION/form.html
			> Choose file.
			> Input = Phylip relaxed - sequential
			> Output = Nexus sequential
			> Line width = as wide as possible

	> Open *.nex files, edit to match # example-data-snakes.nex
			> Delete/modify header part
			> Add charpartition (list of SNPs)
	
	> Run SVDQuarters via PAUP*.
		> Open PAUP: 
			./paup4a159_osx
			
		> Execute data: 
			exe ice_worms.nex
			
		> List options
			svdq ?
		
		> Run the program
			svdquartets nquartets=500000 bootstrap=standard loci=loci nthreads=4 treeFile=ice_worms;
						Time = 6.75 hours.
						All samples grouped by individual
					# nquartets = number of quartets
						[500000 = Exhaustive! 100%.]
					# bootstrap = standard
					# loci = charpartition name (for loci)
					### nreps = number of reps (default = 100)
					nthreads = 4 [number of cpu threads (default = 2)]
				
			> save consensus tree: savetrees file=1.5.1-oneSNP-6048_1.tre;
			> save consensus tree: savetrees file=1.5.1-oneSNP-6048_1.bootstrap.tre savebootp=NodeLabels;

		> Visualized/edited in FigTree.
		
	
## 8. GenoDive - Testing IBD with Mantel tests
		> Mapped populations in Google Earth
		> Calculate straight-line Euclidean distances between sites.
		> Make a geographic distance matrix that matches your FST matrix.
	
	Two matrices:
		> ice_worm_fst.txt
		> ice_worm_geography.txt
		
	In GenoDive GUI:
		1. Open GenoDive
		2. Load matrices
		3. Analysis > Mantel test
			> Correlate with geographic matrix
			> Mantels r
			> 9999 permutations
	
	Performed Mantel tests for 1) all, 2) eastern only, 3) western only, and 4) western only with no Vancouver Island


## 9. fineRADstructure
		# Useful: https://github.com/millanek/fineRADstructure
		# Also: http://cichlid.gurdon.cam.ac.uk/fineRADstructure.html
		# And: http://evomics.org/learning/population-and-speciation-genomics/2018-population-and-speciation-genomics/fineradstructure-activity/

		> Make input file following http://evomics.org/learning/population-and-speciation-genomics/2018-population-and-speciation-genomics/fineradstructure-activity/
			ice_worm_fineRAD.data
			
		> Calculate co-ancestry matrix
			./RADpainter paint ice_worm_fineRAD.data

		> Run fineSTRUCTURE MCMC clusering algorithm.
		
		> Run fineSTRUCTURE tree-building algorithm.
		
		> Visualize results 
			See script below.
			
		
##################################################################
## A simple R example for plotting fineRADstructure output
## Author: Milan Malinsky (millanek@gmail.com), adapted from a Finestructure R Example by Daniel Lawson (dan.lawson@bristol.ac.uk) and using his library of R functions
## Date: 04/04/2016
## Notes:
##    These functions are provided for help working with fineSTRUCTURE output files
## but are not a fully fledged R package for a reason: they are not robust
## and may be expected to work only in some specific cases - often they may require 
## at least minor modifications! USE WITH CAUTION!
## SEE FinestrictureLibrary.R FOR DETAILS OF THE FUNCTIONS
##
## Licence: GPL V3
## 
##    This program is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.

##    This program is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.

##    You should have received a copy of the GNU General Public License
##    along with this program.  If not, see <http://www.gnu.org/licenses/>.


### 1) EDIT THE FOLLOWING THREE LINES TO PROVIDE PATHS TO THE fineRADstructure OUTPUT 
setwd("") ## The directory where the files are located
chunkfile<-"" ## RADpainter output file
mcmcfile<-"" ## finestructure mcmc file
treefile<-".xml" ## finestructure tree file
### 2) EDIT THIS PATH TO WHERE YOU WANT THE PLOTS:
plotsFolder <- "~/fineRAD_Plots"
### 3) SET VALUES FOR THESE VARIABLES: "analysisName" will be included in output plots
analysisName <- "ice_worms_fineRAD";  maxIndv <- 10000; maxPop<-10000


### 4) EDIT THE PATH TO YOUR COPY of FinestructureLibrary.R
source("", chdir = TRUE) # read in the R functions, which also calls the needed packages

### 5) EXECUTE THE CODE ABOVE AND THE REST OF THE CODE BELOW
## make some colours
some.colors<-MakeColorYRP() # these are yellow-red-purple
some.colorsEnd<-MakeColorYRP(final=c(0.2,0.2,0.2)) # as above, but with a dark grey final for capped values
###### READ IN THE CHUNKCOUNT FILE
dataraw<-as.matrix(read.table(chunkfile,row.names=1,header=T,skip=1)) # read in the pairwise coincidence 
###### READ IN THE MCMC FILES
mcmcxml<-xmlTreeParse(mcmcfile) ## read into xml format
mcmcdata<-as.data.frame.myres(mcmcxml) ## convert this into a data frame
###### READ IN THE TREE FILES
treexml<-xmlTreeParse(treefile) ## read the tree as xml format
ttree<-extractTree(treexml) ## extract the tree into ape's phylo format

## Reduce the amount of significant digits printed in the posteror assignment probabilities (numbers shown in the tree):
ttree$node.label[ttree$node.label!=""] <-format(as.numeric(ttree$node.label[ttree$node.label!=""]),digits=2)
 # convert to dendrogram format
tdend<-myapetodend(ttree,factor=1)
## Now we work on the MAP state
mapstate<-extractValue(treexml,"Pop") # map state as a finestructure clustering
mapstatelist<-popAsList(mapstate) # .. and as a list of individuals in populations
popnames<-lapply(mapstatelist,NameSummary) # population names IN A REVERSIBLE FORMAT (I.E LOSSLESS)
## NOTE: if your population labels don't correspond to the format we used (NAME<number>) YOU MAY HAVE TROUBLE HERE. YOU MAY NEED TO RENAME THEM INTO THIS FORM AND DEFINE YOUR POPULATION NAMES IN popnamesplot BELOW
popnamesplot<-lapply(mapstatelist,NameMoreSummary) # a nicer summary of the populations
names(popnames)<-popnamesplot # for nicety only
names(popnamesplot)<-popnamesplot # for nicety only
popdend<-makemydend(tdend,mapstatelist) # use NameSummary to make popdend
popdend<-fixMidpointMembers(popdend) # needed for obscure dendrogram reasons
popdendclear<-makemydend(tdend,mapstatelist,"NameMoreSummary")# use NameMoreSummary to make popdend
popdendclear<-fixMidpointMembers(popdendclear) # needed for obscure dendrogram reasons

	
########################
## Plot 1: COANCESTRY MATRIX
fullorder<-labels(tdend) # the order according to the tree
datamatrix<-dataraw[fullorder,fullorder] # reorder the data matrix

tmpmat<-datamatrix 
tmpmat[tmpmat>maxIndv]<-maxIndv #  # cap the heatmap
pdf(file=paste(plotsFolder,analysisName,"-SimpleCoancestry.pdf",sep=""),height=25,width=25)
plotFinestructure(tmpmat,dimnames(tmpmat)[[1]],dend=tdend,cols=some.colorsEnd,cex.axis=1.1,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=1.2))
dev.off()

########################
## Plot 2: POPULATIONS AND COANCESTRY AVERAGES
popmeanmatrix<-getPopMeanMatrix(datamatrix,mapstatelist)

tmpmat<-popmeanmatrix
tmpmat[tmpmat>maxPop]<-maxPop # cap the heatmap
pdf(file=paste(plotsFolder,analysisName,"-PopAveragedCoancestry.pdf",sep=""),height=20,width=20)
plotFinestructure(tmpmat,dimnames(tmpmat)[[1]],dend=tdend,cols=some.colorsEnd,cex.axis=1.1,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=1.2))
dev.off()


########################
## Plot 3: POPULATIONS AND COANCESTRY AVERAGES WITH PERHAPS MORE INFORMATIVE LABELS
mappopcorrectorder<-NameExpand(labels(popdend))
mappopsizes<-sapply(mappopcorrectorder,length)
labellocs<-PopCenters(mappopsizes)
xcrt=0
ycrt=45

pdf(file=paste(plotsFolder,analysisName,"-PopAveragedCoancestry2.pdf",sep=""),height=25,width=25)
plotFinestructure(tmpmat,dimnames(tmpmat)[[1]],labelsx=labels(popdendclear),labelsatx=labellocs,xcrt=xcrt,cols=some.colorsEnd,ycrt=ycrt,dend=tdend,cex.axis=1.1,edgePar=list(p.lwd=0,t.srt=90,t.off=-0.1,t.cex=1.2),hmmar=c(3,0,0,1))
dev.off()
