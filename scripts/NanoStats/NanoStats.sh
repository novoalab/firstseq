## FROM SAMTOOLS MPILEUP TO STATS ##
## Eva Maria Novoa, 2019 ##

if [[ $1 == "-h" ]]; then
echo "Usage: $0 <bam> <ref>"
echo "Expects mpileup in the same folder, in the form of <file.bam>"
fi

# INPUT
i=`echo $1 | sed 's/\.bam//'`
echo "# Analyzing $1"

ref=$2


# SCRIPT
#Create mpileup file
samtools mpileup -d 5000000 -Q 0 -f $ref -s $i.bam > $i.mpileup
####
###############################
## PART 1: QUALITIES
###############################

echo "# Convert qualities to ordinal numbers"
python2 scripts/1.base_quality_chr2ord.keep_all.py $i.mpileup > $i.mpileup.digit

echo "# Get mean and median quality scores (produces output with header)"
python2 scripts/2.get_mean_median_qualities.py $i.mpileup.digit > $i.mean_median_quals

###############################
## PART 2: COVERAGE AND REF_NUC
###############################

echo "# Get general file info (coverage, ref_nuc)"
echo -ne "chr\tpos\tref_nuc\tcoverage\n" > $i.info
cut -f 1-4 $i.mpileup >> $i.info

###############################
## PART 3: RT-STOP
###############################

echo "# Get number of rt-stop, insertion and deletions (produces file with header)"
echo "0" > $i.mpileup.digit.rtstop.tmp
cat $i.mpileup.digit | awk {'print $5'} | awk -F "^" '{print NF-1}' >> $i.mpileup.digit.rtstop.tmp 
cat $i.mpileup.digit | awk {'print $5'} | awk -F "-" '{print NF-1}'  > $i.mpileup.digit.del
cat $i.mpileup.digit | awk {'print $5'} | awk -F "+" '{print NF-1}'  > $i.mpileup.digit.ins
sed '$ d' $i.mpileup.digit.rtstop.tmp > $i.mpileup.digit.rtstop # remove last line as we pushed everything 1
echo -ne "rtstop\tins\tdel\n" > $i.3moreinfo
paste $i.mpileup.digit.rtstop $i.mpileup.digit.del $i.mpileup.digit.ins >> $i.3moreinfo

###############################
## PART 4: MISMATCHES
###############################

echo "# Get mismatches"
perl scripts/3.pileup2base.pl  $i.mpileup 0 $i.mismatches.tmp 
# Previous script generates both mismatches for FWD and REV strands separately, so they will be merged afterwards
cut -f 4-12 $i.mismatches.tmp | tail -n +2 > $i.mismatches.tmp2
cut -f 1-3 $i.mismatches.tmp | tail -n +2 > $i.mismatches.tmp3
echo -ne "chr\tpos\tref_nuc\tA\tT\tC\tG\n" > $i.mismatches
 awk 'BEGIN{print "sum" > "kk"}NR>0{print $1+$5"\t"$2+$6"\t"$3+$7"\t"$4+$8}' $i.mismatches.tmp2 > $i.mismatches.tmp4
paste $i.mismatches.tmp3 $i.mismatches.tmp4 >> $i.mismatches
rm $i.mismatches.tmp $i.mismatches.tmp2 $i.mismatches.tmp3 $i.mismatches.tmp4 kk
 
###############################
## PART 5: INDELS
###############################

echo "# Identity of insertions and deletions"
perl scripts/4.pileup2baseindel.pl -i $i.mpileup
cut -f 12-13 sample1.txt > $i.indel_info

echo "# Merge mismatches and indel info"
paste $i.mismatches $i.indel_info > $i.mismatch_and_indels


###############################
## PART 6: CLEAN UP
###############################

# Merge all files
paste $i.info $i.mean_median_quals $i.3moreinfo $i.mismatches | awk '{$10=$11=$12=""; print $0}' > $i.STATS      



# Remove temporary files
rm $i.info $i.mean_median_quals $i.3moreinfo $i.mpileup.digit.rtstop $i.mpileup.digit.del $i.mpileup.digit.ins $i.mpileup.digit.rtstop.tmp $i.mpileup.digit $i.mismatch_and_indels $i.mismatches $i.indel_info $i.mpileup 

# Script completed!
echo "Script completed! Have a nice day! :-)"
