[![DOI](https://zenodo.org/badge/508619527.svg)](https://zenodo.org/badge/latestdoi/508619527)


# Aneuploidy detection method used for detecting 'simple' and 'complex' aneuploidies
Script uses simply alignment and a hard threshold for deviation in order to detect changes in coverage <br/>
The hard threshold is calculated using the ploidy ; 0.7*(1/n) ; with 0.7 set as a conservative measure of deviation <br/>
For example, for a triploid 0.7*1/3=0.21, meaning a deviation (+-) of 0.21 times the median genome-wide coverage would be considered as a region exhibiting a change in copy number <br/>

The script needs a few inputs:
  1. A genome assembly
  2. Paired illumina reads (you can use other read types, just will need to modify the alignment steps)
  3. A tsv file with the sample name in the first column and the ploidy value in the second column e.g. '2' for diploid

The script requires only three installed tools in path:
  1. bwa (alignment with mem)
  2. samtools (alignment output pre-processing)
  3. bedtools (calculating coverage, formatting bins etc)

As the user there are only a few variables that are required as input depending on your sample
  1. The path to the reference genome = genome
  2. The path to the tsv file with sample names and ploidy status
  3. A path to all the illumina files (put all of them in a single folder and name them as '*sample*_1.fq.gz' and '*sample*_2.fq.gz' to simplify the process)
  4. The size of the window used for binning (the region median-averaged to smooth the coverage distribution) = window
      NOTE: Defined later, any region considered to contain a copy number change must be bigger than this size
  5. The size of the slide for each window (the distance shifted along the genome each step in order to re-calculate the coverage median) = slide
  6. A prefix which I usually just use the reference sample name etc
  7. A name of a mitochondrial genome if present within the reference to be removed
  8. A file containing the centromere positions

Additionally, half a window in size will be removed from the ends of all contigs/scaffolds/chromosomes within the genome provided  <br/>
Personally I am/was generally working on yeast, and developed this using *Saccharomyces cerevisiae*. In this respect, for the binning and slide, I used 30kb windows with a 10kb slide because subtelomeric regions in the reference genome are ~30kb. Therefore I found that the trimming of the edges using a half-window size generally worked well in order to heavily smooth out subtelomeric regions which often containing higher coverage. <br/>

### First steps here are to set variables, generate the alignment and calculate the coverage, both raw and binned.

    ###SET THESE VARIABLES
    ##path to reference genome
    reference="path/to/genome.fa"
    ##prefix name
    prefix="genome"
    ##path to folder containing all the illumina files ending with ${sample}_1.fq.gz or ${sample}_2.fq.gz
    illumina="path/to/illumina_folder"
    ##path to tsv file containing the sample names in column 1 and ploody status in column 2
    list="path/to/list.tsv"
    ##the window size for median-average binning (in bp)
    window=30000
    ##distance for the window to move before recalculating the window (in bp)
    slide=10000
    ##name of mitochondrial contig/scaffold/genome
    mito="chrMT"
    ##File containing positions of centromeres in your reference, with 1-3 columns in order chromosome,start-position,end-position
    centromeres="centromeres.tsv"
    
    ##some preprocessing
    ##index reference genome
    bwa index ${reference}
    ##same as path to reference genoe but remove the fasta suffix (removes both .fasta and .fa incase)
    reference2=$( echo $reference | sed 's/.fa$//g ; s/.fasta$//g'  )
    ##generate a folder to place the bam files
    mkdir ${prefix}.illumina_alignment
    ##generate a folder to place the coverage files
    mkdir ${prefix}.illumina_alignment.cov
    
    ##now take all strains and map; them transforming the data directly into a sorted bam file
    ##read in the strain list and sample just the strain column and run the loop per strain for mapping
    cat $list | cut -f1 | while read strain
    do
      bwa mem ${reference} ${illumina}/${strain}_1.fq.gz ${illumina}/${strain}_2.fq.gz | samtools sort -o ${prefix}.illumina_alignment/${strain}.bwamem_${prefix}.sorted.bam -
    done
    ##calculate the coverage per bam file using bedtools genomecov, using the -ibam option
    bedtools genomecov -d -ibam ${prefix}.illumina_alignment/${strain}.bwamem_${prefix}.sorted.bam | gzip > ${prefix}.illumina_alignment.cov/${strain}.bwamem_${prefix}.cov.tsv.gz
    ##calculate genome wide median
    median=$( zcat ${prefix}.illumina_alignment.cov/${strain}.bwamem_${prefix}.cov.tsv.gz | awk -v mito="$mito" '{if($1 != mito && $3 != "0") print $3}' | sort -n | awk '{ a[i++]=$1} END{x=int((i+1)/2); if(x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1];}' )
    ##generate another file with the genome wide median coverage normalised coverage in the last column
    zcat ${prefix}.illumina_alignment.cov/${strain}.bwamem_${prefix}.cov.tsv.gz | awk -v median="$median" '{print $0"\t"$3/median}' | gzip > ${prefix}.illumina_alignment.cov/${strain}.bwamem_${prefix}.cov_medianNORM.tsv.gz


    ##reformatting window and slide numbers for directory generation in kb
    window2=$( echo $window | awk '{print $0/1000}' )
    slide2=$( echo $slide | awk '{print $0/1000}' )
    
    ##generate a folder to place the binned coverage files
    mkdir ${prefix}.illumina_alignment.cov.${window2}kbwindow_${slide2}kbsliding

    
    ##make a window file 
    ##first need a reference bed file, can use the fai index generated by samtools
    samtools faidx ${reference}
    cat ${reference}.fai | cut -f1-2 > ${reference2}.bed
    ##now split the reference bed file up into the above prescribed bins
    bedtools makewindows -w ${window} -s ${slide} -g ${reference2}.bed > ${reference2}.${window2}kbwindow_${slide2}kbslide.bed
    ##get the coverage file (slightly modify by giving a range for the single basepair coverage value) then use bedtools map to overlap with the reference derived window file to create median-averaged bins
    zcat ${prefix}.illumina_alignment.cov/${strain}.bwamem_${prefix}.cov.tsv.gz | awk '{print $1"\t"$2"\t"$2"\t"$3}' |\
    bedtools map -b - -a ${reference2}.${window2}kbwindow_${slide2}kbslide.bed -c 4 -o median | awk -v median="$median" '{print $0"\t"$4/median}' > ${prefix}.illumina_alignment.cov.${window2}kbwindow_${slide2}kbsliding/${strain}.bwamem_${prefix}.cov_medianNORM.${window2}kbwindow_${slide2}kbsliding.tsv
    done
    
    
    
### Next steps are now to actually analyse the coverage for changes that would be considered, based on the ploidy, deviating sufficiently to be considered a change in copy number
    
    ##As described above the hard threshold for this is 0.7*1/n, where n is the ploidy
    ##remove any files potentially produced by an earlier round
    rm raw_output*
    ##add header to coming files
    echo "strain,threshold,chromosome,start,end,size,region_cov,median_cov,change_cov,change,aneuploidy,centromere" | sed 's/,/\t/g' > raw_output.tsv
    echo "strain,chromosome,centromere_bin_size,sum_bin_size,chr_size_unadjusted,proportion_unadjusted,size_diff_unadjusted,chr_size_adjusted,proportion_adjusted,size_diff_adjusted,aneuploidy" | sed 's/,/\t/g' > raw_output.sumcen.tsv
    echo "strain,chromosome,largest_bin_size,sum_bin_size,chr_size_unadjusted,proportion_unadjusted,size_diff_unadjusted,chr_size_adjusted,proportion_adjusted,size_diff_adjusted,aneuploidy" | sed 's/,/\t/g' > raw_output.sumall.tsv


    cat $list | awk -F "\t" '{print $1}' | while read strain
    do
	    echo $strain
	    ##get the ploidy value for each strain then define by which how much the coverage should change for a loss or gain of a portion
	    ploidy=$( cat $list | awk -F "\t" -v strain="$strain" '{if($1 == strain) {print $2}}' | awk '{print ((1/$0)*0.7)}' )
	    ##get the median coverage
	    median=$( zcat ${prefix}.illumina_alignment.cov/${strain}.bwamem_${prefix}.cov.tsv.gz | awk -v mito="$mito" '{if($1 != mito && $3 > 0) {print $3}}' | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'  )
	    ##get the windows and remove the mitochondiral genome
	    cat ${prefix}.illumina_alignment.cov.${window2}kbwindow_${slide2}kbsliding/${strain}.bwamem_${prefix}.cov_medianNORM.${window2}kbwindow_${slide2}kbsliding.tsv |\
	    grep -v ${mito} | awk '{print $1}' | sort -u | while read chromosome
	    do
		    ##get the centromere position
		    cenpos=$( cat ${centromeres} | awk -v chr="$chr" '{if($1 == chr) {print ($2+$3)/2}}' )
		    ##get the chromosome size
		    chrsize=$( grep $chromosome ${reference2}.bed | awk '{print $2}' )
		    ##initial filtering of the alignment
		    ##-get those regions which deviate from the ploidy sufficiently to be considered a region of duplication
		    ##-get those regions which are larger than the window size
		    ##-remove those regions with low coverage, here we use less than 10 percent of the coverage
		    ##-remove the first and last half window size in each chromosome
		    ##-remove any alignments smaller than the bin size; can happen due to wierd binning; although shouldnt be an issue now
		    cat ${prefix}.illumina_alignment.cov.${window2}kbwindow_${slide2}kbsliding/${strain}.bwamem_${prefix}.cov_medianNORM.${window2}kbwindow_${slide2}kbsliding.tsv |\
			  awk -v slide="$slide" -v mito="$mito" '{if($1 ~ mito ) {} else {print $0}}' |\
			  awk -v window="$window" -v chr="$chromosome" -v median="$median" -v ploidy="$ploidy" -v chrsize="$chrsize" -v slide="$slide" '{if($1 == chr && $4 >= median*(1+ploidy) && $4 != "." && ($3-$2) >= window && $4 > (0.1*median) && $2 >= (window/2) && $3 <= (chrsize-(window/2)) || $1 == chr && $4 <= median*(1-ploidy) && $4 != "." && ($3-$2) >= window && $4 > (0.1*median) && $2 >= (window/2) && $3 <= (chrsize-(window/2)) ) {print $0"\t"median"\t"$4-median"\t"$4/median}}' |\
			  awk '{if($7 > 1) {print $0"\tincrease"} else {print $0"\tdecrease"}}' > temp.txt
			
		    ##split the dataset into those with an increase and those with a decrease
		    ##concatenate the deviating windows
		    ##to conatenate windows need to overlap

		    cat temp.txt | grep 'increase' |\
			  awk -v median="$median" -v slide="$slide"  '{if($1==chr && $2<=end+slide && $8==diff) {end=$3; print $0"\t"count} else if(chr=="") {count=count+1; chr=$1; end=$3; diff=$8; print $0"\t"count} else if($1!=chr || $2>end+slide || $8!=diff) {count=count+1; print $0"\t"count; chr=$1; cov=$4; start=$2 ;end=$3; diff=$8}}' > temp2.txt
		    cat temp2.txt | awk '{print $9}' | sort -u | while read group
		    do
			    cat temp2.txt |\
				  awk -v group="$group" '{if($9 == group) {print $0}}' |\
				  awk -v median="$median" -v strain="$strain" -v ploidy="$ploidy" '{if(count=="") {chr=$1; start=$2; end=$3; cov=cov+$4; count=count+1} else if(count != "") {chr=$1; end=$3; cov=cov+$4; count=count+1}} END{print strain"\t"ploidy"\t"chr"\t"start"\t"end"\t"end-start"\t"cov/count"\t"median"\t"(cov/count)/median}' |\
				  awk -v cenpos="$cenpos" '{if($4 <= cenpos && $5 >= cenpos && int(($9-1)/$2) > 0) {print $0"\tincrease\t+"int(($9-1)/$2)$3"\tcentromere"} else if($4 <= cenpos && $5 >= cenpos && int(($9-1)/$2) < 0) {print $0"\tincrease\t"int(($9-1)/$2)$3"\tcentromere"} else if(int(($9-1)/$2) > 0) {print $0"\tincrease\t+"int(($9-1)/$2)$3} else if(int(($9-1)/$2) < 0) {print $0"\tincrease\t"int(($9-1)/$2)$3}}' | awk '{gsub(/chr/,"*",$11);print}' | sed 's/ /\t/g' |\
				  awk -v window="$window" '{if($6 > window) print $0}'
		    done >> raw_output.tsv

		    cat temp.txt | grep 'decrease' |\
			  awk -v median="$median" -v slide="$slide" '{if($1==chr && $2<=end+slide && $8==diff) {end=$3; print $0"\t"count} else if(chr=="") {count=count+1; chr=$1; end=$3; diff=$8; print $0"\t"count} else if($1!=chr || $2>end+slide || $8!=diff) {count=count+1; print $0"\t"count; chr=$1; cov=$4; start=$2 ;end=$3; diff=$8}}' > temp2.txt
		    cat temp2.txt | awk '{print $9}' | sort -u | while read group
		    do
			    cat temp2.txt |\
				  awk -v group="$group" '{if($9 == group) {print $0}}' |\
				  awk -v median="$median" -v strain="$strain" -v ploidy="$ploidy" '{if(count=="") {chr=$1; start=$2; end=$3; cov=cov+$4; count=count+1} else if(count != "") {chr=$1; end=$3; cov=cov+$4; count=count+1}} END{print strain"\t"ploidy"\t"chr"\t"start"\t"end"\t"end-start"\t"cov/count"\t"median"\t"(cov/count)/median}' |\
				  awk -v cenpos="$cenpos" '{if($4 <= cenpos && $5 >= cenpos && int(($9-1)/$2) > 0) {print $0"\tdecrease\t+"int(($9-1)/$2)$3"\tcentromere"} else if($4 <= cenpos && $5 >= cenpos && int(($9-1)/$2) < 0) {print $0"\tdecrease\t"int(($9-1)/$2)$3"\tcentromere"} else if(int(($9-1)/$2) > 0) {print $0"\tdecrease\t+"int(($9-1)/$2)$3} else if(int(($9-1)/$2) < 0) {print $0"\tdecrease\t"int(($9-1)/$2)$3}}' | awk '{gsub(/chr/,"*",$11);print}' | sed 's/ /\t/g' |\
				  awk -v window="$window" '{if($6 > window) print $0}'
		    done >> raw_output.tsv
		
		    rm temp.txt
		    rm temp2.txt
      
	    done

      ##removed the gff-gene count step from here
      
	    ##start summing up the entire chromosome event
	    ##no filtering here but keep track of initial region around centromere with increased coverage
	    awk -v strain="$strain" '{if($1 == strain) print}' raw_output.tsv | grep 'centromere' | awk '{print $3}' | while read chr
	    do
		    ##get the other name of the chromosome
		    chrsize=$( grep $chr ${reference2}.bed | awk '{print $2}' )
		    ##get initial size of concatenated region around centromere
		    initialsize=$( awk -v strain="$strain" -v chr="$chr" '{if($1 == strain && $3==chr && $13 == "centromere") print $6}' raw_output.tsv )
		    ##get size of chromosome minus the 10kb taken from each end
		    sizetotal=$( awk -v chr="$chr" -v slide="$slide" -v window="$window" '{if($1 == chr) {print $2-(2*(window/2))}}' ${reference2}.bed  )
		    ##get the size of the above chromosome minus also the regions which had less than 10 percent of the coverage
		    sizefinal=$( zcat ${prefix}.illumina_alignment.cov/${strain}.bwamem_${prefix}.cov.tsv.gz | awk -v window="$window" -v chrsize="$chrsize" '{if($2 >= (window/2) && $2 <= (chrsize-(window/2))) {print}}' | awk -v median="$median" -v chr="$chr" -v size="$sizetotal" '{if($1 == chr && $3 <= (.1*median) ) count=count+1} END{print size-count}' )
		    ##get the direction of the change, +1 or +2 or -1 etc
		    change=$( awk -v strain="$strain" -v chr="$chr" '{if($1 == strain && $3 == chr) print}' raw_output.tsv |  grep 'centromere' | awk '{print $11}'  )
		    ##concantenate all the regions in the chromosome of interest that had the same change of direction and the number of genes
		    cat raw_output.tsv | awk -v strain="$strain" -v chr="$chr" -v sizet="$sizetotal" -v size="$sizefinal" -v change="$change" -v initialsize="$initialsize" '{if($1 == strain && $3 == chr && $11==change) {sum=sum+$6; genes=genes+$12}} END{print strain"\t"chr"\t"initialsize"\t"sum"\t"sizet"\t"sum/sizet"\t"sizet-sum"\t"size"\t"sum/size"\t"size-sum"\t"change"\t"genes}' >> raw_output.sumcen.tsv
	    done  

	    ##same thing as above but now use all chromosomes, not just those covering a centromere
	    ##helps with detecting complex aneuploidies
	    ##instead of originally a centromere region being the seed, will use the largest overlapping region
	    awk -v strain="$strain" '{if($1 == strain) print}' raw_output.tsv | awk '{print $3}' | sort -u | while read chr
	    do
		    chrsize=$( grep $chr ${reference2}.bed | awk '{print $2}' )
		    ##get initial size of concatenated region around centromere
		    initialsize=$( awk -v strain="$strain" -v chr="$chr" '{if($1 == strain && $3==chr && $6 > size) {size=$6}} END{print size}' raw_output.tsv )
		    ##get size of chromosome minus the 10kb taken from each end
		    sizetotal=$( awk -v chr="$chr" -v slide="$slide" -v window="$window" '{if($1 == chr) {print $2-(2*(window/2))}}' ${reference2}.bed  )
		    ##get the size of the above chromosome minus also the regions which had less than 10 percent of the coverage
		    sizefinal=$( zcat ${prefix}.illumina_alignment.cov/${strain}.bwamem_${prefix}.cov.tsv.gz | awk -v window="$window" -v chrsize="$chrsize" '{if($2 >= (window/2) && $2 <= (chrsize-(window/2))) {print}}' | awk -v median="$median" -v chr="$chr" -v size="$sizetotal" '{if($1 == chr && $3 <= (.1*median) ) count=count+1} END{print size-count}' )
		    ##get the direction of the change, +1 or +2 or -1 etc
		    change=$( awk -v strain="$strain" -v chr="$chr" -v initialsize="$initialsize" '{if($1 == strain && $3 == chr && $6 == initialsize) print}' raw_output.tsv | head -n1 | awk '{print $11}'  )
		    #concantenate all the regions in the chromosome of interest that had the same change of direction
		    cat raw_output.tsv | awk -v strain="$strain" -v chr="$chr" -v sizet="$sizetotal" -v size="$sizefinal" -v change="$change" -v initialsize="$initialsize" '{if($1 == strain && $3 == chr && $11==change) {sum=sum+$6; genes=genes+$12}} END{print strain"\t"chr"\t"initialsize"\t"sum"\t"sizet"\t"sum/sizet"\t"sizet-sum"\t"size"\t"sum/size"\t"size-sum"\t"change"\t"genes}' >> raw_output.sumall.tsv
	    done

	    ##same thing as above but now use only chromosomal regions withoout centromere
	    ##helps with detecting complex aneuploidies
	    ##instead of originally a centromere region being the seed, will use the largest overlapping region
	    awk -v strain="$strain" '{if($1 == strain) print}' raw_output.tsv | grep -v centromere | awk '{print $11}' | sort -u | while read event
	    do
		    ##get the chromosome involved
		    chr=$( awk -v strain="$strain" -v event="$event" '{if($1 == strain && $11 == event) print}' raw_output.tsv | head -n1 | awk '{print $3}'  )
		    chrsize=$( grep $chr ${reference2}.bed | awk '{print $2}' )
		    ##get initial size of concatenated region around centromere
		    initialsize=$( awk -v strain="$strain" -v event="$event" '{if($1 == strain && $11==event && $6 > size) {size=$6}} END{print size}' raw_output.tsv )
		    ##get size of chromosome minus the 10kb taken from each end
		    sizetotal=$( awk -v chr="$chr" -v slide="$slide" -v window="$window" '{if($1 == chr) {print $2-(2*(window/2))}}' ${reference2}.bed  )
		    ##get the size of the above chromosome minus also the regions which had less than 10 percent of the coverage
		    sizefinal=$( zcat ${prefix}.illumina_alignment.cov/${strain}.bwamem_${prefix}.cov.tsv.gz | awk -v window="$window" -v chrsize="$chrsize" '{if($2 >= (window/2) && $2 <= (chrsize-(window/2))) {print}}' | awk -v median="$median" -v chr="$chr" -v size="$sizetotal" '{if($1 == chr && $3 <= (.1*median) ) count=count+1} END{print size-count}' )
		    ##get the other name of the chromosome
		    ##concantenate all the regions in the chromosome of interest that had the same change of direction
		    cat raw_output.tsv | awk -v strain="$strain" -v chr="$chr" -v sizet="$sizetotal" -v size="$sizefinal" -v event="$event" -v initialsize="$initialsize" '{if($1 == strain && $3 == chr && $11==event) {sum=sum+$6}} END{print strain"\t"chr"\t"initialsize"\t"sum"\t"sizet"\t"sum/sizet"\t"sizet-sum"\t"size"\t"sum/size"\t"size-sum"\t"event}' >> raw_output.sumnocen.tsv

	    done

    done
    rm raw_output.temp.tsv

    ##need to re process the list with no centromeres to remove some that made it through
    cat raw_output.sumcen.tsv raw_output.sumnocen.tsv | awk '{print $1"\t"$11}' | sort | uniq -c | awk '{if($1 == "1") print $2"\t"$3}' | while read line
    do
      strain=$( echo $line | awk '{print $1}' )
      event=$( echo $line | awk '{print $2}' )
      cat raw_output.sumnocen.tsv | awk -v strain="$strain" -v event="$event" '{if($1 == strain && $11 == event) {print $0}}' >> raw_output.sumnocen.temp.tsv
    done
    echo "strain,chromosome,largest_bin_size,sum_bin_size,chr_size_unadjusted,proportion_unadjusted,size_diff_unadjusted,chr_size_adjusted,proportion_adjusted,size_diff_adjusted,aneuploidy" | sed 's/,/\t/g' > raw_output.sumnocen.tsv
    cat raw_output.sumnocen.temp.tsv >> raw_output.sumnocen.tsv
    rm raw_output.sumnocen.temp.tsv

### Now that we have our regions of copy-number changes (with and without coverage a centromere) we can generate our candidate simple and complex aneuploidy list
    
    ##first we simply remove any centromere-related (CR) event smaller than 50kb
    ##next we label all of the remaining CRs as complex if the difference between the size of the summed CR region and the chromosome size (minus the 30kb removed from the ends) is greater than 100kb. If it is smaller than 100kb, the CR is labelled as simple. 
    echo "strain,chromosome,largest_bin_size,sum_bin_size,chr_size_unadjusted,proportion_unadjusted,size_diff_unadjusted,chr_size_adjusted,proportion_adjusted,size_diff_adjusted,aneuploidy,aneuploidy_type" | sed 's/,/\t/g' > raw_output.sumcen.candidates.tsv
	cat raw_output.sumcen.tsv | tail -n+2 | awk '{if($7 > 100000 && $3 > 50000) {print $0"\tcomplex"} else if($3 > 50000) {print $0"\tsimple"}}' >> raw_output.sumcen.candidates.tsv

### These aneuploidies, considered complex or simple, now need to be manually currated in order to remove false positives. To do this we can plot images of the relative coverage and scrutinise them manually to indicate whether they are truely simple, truely complex, or even aneuploid.
    
    ##in order to do this the relative coverage for each aneuploidy chromosome will be calculated and then plotted.
    ##this will allow me to manually view the alignment for each of these aneuploidies and remove any that do not conform due to strange alignment properties such as the smiley effect or large fluctuations
    ##getting relative coverage per chromosome
    ##calculate the median coverage and use this to calculate the relative coverage across the whole genome
    ##printing relative coverage info next to details on the aneuploidy etc
    echo "chromosome;start;end;coverage;coverage_rel;strain;proportion_chromosome;loss_chromosome;ploidy;aneuploidy_type" | sed 's/;/\t/g' > raw_output.sumcen.candidates.aneuploidy_relativecov.tsv
    cat raw_output.sumcen.candidates.tsv | tail -n+2 | cut -f2 | sort -u | while read chr
    do
      cat raw_output.sumcen.candidates.tsv | awk -v chr="$chr" '{if($2 == chr) print $1"\t"$6"\t"$7"\t"$12}' | while read line
      do
        strain=$( echo $line | awk '{print $1}')
        proportion=$( echo $line | awk '{print $2}'  )
        loss=$( echo $line | awk '{print $3}'  )
        candidacy=$( echo $line | awk '{print $4}'  )
        ploidy=$( cat $list | awk -F "\t" -v strain="$strain" '{if($1 == strain) print $3}' )
        median=$( zcat ${prefix}.illumina_alignment.cov/${strain}.bwamem_${prefix}.cov.tsv.gz | awk -v mito="$mito" '{if($1 == mito && $3 > 0) {print $3}}' | sort -n | awk ' { a[i++]=$1; } END { x=int((i+1)/2); if (x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1]; }'  )
        cat ${prefix}.illumina_alignment.cov.${window2}kbwindow_${slide2}kbsliding/${strain}.bwamem_${prefix}.cov_medianNORM.${window2}kbwindow_${slide2}kbsliding.tsv | awk -v median="$median" -v strain="$strain" -v prop="$proportion" -v loss="$loss" -v chr="$chr" -v ploidy="$ploidy" -v candidacy="$candidacy" '{if($1 == chr ) {print $0"\t"$4/median"\t"strain"\t"prop"\t"loss"\t"ploidy"\t"candidacy}}' | sed 's/ /\t/g' 
      done >> raw_output.sumcen.candidates.aneuploidy_relativecov.tsv
    done 

    
NOT FINISHED IMPORTING STEPS ONTO GITHUB, CURRENTLY DOING FILE3, START STEP FOR GENERATING RELATIVE COVERAGE PLOTS
