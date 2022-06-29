# aneuploidy_detection
Script uses simply alignment and a hard threshold for deviation in order to detect changes in coverage


The script needs a few inputs:
  1. A genome assembly
  2. Paired illumina reads (you can use other read types, just will need to modify the alignment steps)
  3. A tsv file with the sample name in the first column and the ploidy value in the second column e.g. '2' for diploid

The script requires only three installed tools in path:
  1. bwa (alignment)
  2. samtools (alignment output pre-processing)
  3. bedtools (calculating coverage, formatting bins etc)

As the user there are only a few variables that are required as input depending on your sample
  1. The path to the reference genome = genome
  2. The path to the tsv file with sample names and ploidy status
  3. A path to all the illumina files (put all of them in a single folder and name them as '*sample*_1.fq.gz' and '*sample*_2.fq.gz' to simplify the process)
  4. The size of the window used for binning (the region median-averaged to smooth the coverage distribution) = window
  5. The size of the slide for each window (the distance shifted along the genome each step in order to re-calculate the coverage median) = slide
  6. 



    reference=~/Documents/projects/Gcandidum/reference/CLIB918.fa
    ##same as above but without the file extension
    reference2=~/Documents/projects/Gcandidum/reference/CLIB918
    ##reference name
    refname="CLIB918"
    bwa index ${reference}
    ## now take all 114 strains and map them transformaing the data directly into a bam file
    ls ../illumina_*/*.gz | sed 's/_[12].fq.gz//g' | sort -u | while read strain1
    do
    strain2=$( echo $strain1 | awk -F "/" '{print $3}' )
    bwa mem ${reference} ${strain1}_1.fq.gz ${strain1}_2.fq.gz | samtools sort -o illumina_vs_${refname}/${strain2}.bwamem_${refname}.sorted.bam -
    done
    mkdir illumina_vs_${refname}_cov 
    ls illumina_vs_${refname}/ | awk -F ".bwamem" '{print $1}' | sort -u | while read strain
    do
    bedtools genomecov -d -ibam illumina_vs_${refname}/${strain}.bwamem_${refname}.sorted_RG_markdup.bam | gzip > illumina_vs_${refname}_cov/${strain}.bwamem_${refname}.cov.tsv.gz
    median=$( zcat illumina_vs_${refname}_cov/${strain}.bwamem_${refname}.cov.tsv.gz | awk '{if($3 != "0") print $3}' | sort -n | awk '{ a[i++]=$1} END{x=int((i+1)/2); if(x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1];}' )
    zcat illumina_vs_${refname}_cov/${strain}.bwamem_${refname}.cov.tsv.gz | awk -v median="$median" '{print $0"\t"$3/median}' | gzip > illumina_vs_${refname}_cov/${strain}.bwamem_${refname}.cov_medianNORM.tsv.gz
    done
    window=30000
    ##distance for the window to move before recalculating the window
    slide=10000
    ##reformatting window and slide numbers for directory generation in kb
    window2=$( echo $window | awk '{print $0/1000}' )
    slide2=$( echo $slide | awk '{print $0/1000}' )
    ##make a window file 
    ##first need a reference bed file, can use the fai index
    cat ${reference}.fai | cut -f1-2 > ${refname}.bed
    ##now split it up into the above prescribed bins
    bedtools makewindows -w ${window} -s ${slide} -g ${refname}.bed > ${refname}.${window2}kbwindow_${slide2}kbslide.bed
    ##now the binning and finding the median
    ##plus we can immediately add the normalised median coverage as before
    mkdir illumina_vs_${refname}_cov_${window2}kbwindow_${slide2}kbsliding
    ls illumina_vs_${refname}_cov/ | awk -F "." '{print $1}' | sort -u | while read strain
    do
    median=$( zcat illumina_vs_${refname}_cov/${strain}.bwamem_${refname}.cov.tsv.gz | awk '{if($3 != "0") print $3}' | sort -n | awk '{ a[i++]=$1} END{x=int((i+1)/2); if(x < (i+1)/2) print (a[x-1]+a[x])/2; else print a[x-1];}' )
    ##get the coverage file then use bedtools map to overlap with the reference derived window file created a single time above
    zcat illumina_vs_${refname}_cov/${strain}.bwamem_${refname}.cov.tsv.gz | awk '{print $1"\t"$2"\t"$2"\t"$3}' |\
    bedtools map -b - -a ${refname}.${window2}kbwindow_${slide2}kbslide.bed -c 4 -o median | awk -v median="$median" '{print $0"\t"$4/median}' > illumina_vs_${refname}_cov_${window2}kbwindow_${slide2}kbsliding/${strain}.bwamem_${refname}.cov_medianNORM.${window2}kbwindow_${slide2}kbsliding.tsv
    done
