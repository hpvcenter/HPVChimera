if_fasta_file() {
 file=$1
 extension="${file##*.}"
 if [ ! -f "$file" ] || [ $extension != "fasta" ]; then
     echo "This is not a .fasta file or the $file does exist."
 fi

}

extract_id_indivial_files() {

cat $fasta_file |sed $'s/[^[:print:]\t]//g'  |awk '{if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fa")}print $0 > filename}'
echo -e "Fasta file well divided"

}

divide_in_three(){
echo "-------- Divide in three -----------"

cd $root_folder
cd workspace
ls | grep -v fasta| grep -v output | grep -v div| while read filename;do
    extension="${filename%.*}"
    is_candidate=$(cat $result_file |grep -w $filename| wc -l) 
    if [ "$is_candidate" -eq "0" ]; then
      echo -e "$filename is a candidate"
      long=$(cat $filename | tail -n +2|sed "/\w+/g" |wc -m)
      sequence=$(cat $filename | tail -n +2)
      L1=$(expr $long / 3)
      L11=$(expr $L1 - 1)
      L2=$(expr $long \* 2)
      L2=$(expr $L2 / 3)
      echo $filename --- $long "----" $L1 "----" $L2
      var1="${sequence:0:$L1}"
      var2="${sequence:$L1:$L2-$L1}"
      var3="${sequence:$L2}"
      mkdir -p div
      echo -e ">"$extension"_1">div/"$extension"_1.fasta
      echo -e $var1>>div/"$extension"_1.fasta
      echo -e ">"$extension"_2">div/"$extension"_2.fasta
      echo -e $var2>>div/"$extension"_2.fasta
      echo -e ">"$extension"_3">div/"$extension"_3.fasta
      echo -e $var3>>div/"$extension"_3.fasta
   fi
done

}

perc_identity_90()
{
 ls | grep fasta | while read filename;do
   mkdir -p output
   blastn -db $hpv_database -query "$filename" -perc_identity 90 -out output/"$filename.out"

 done
}

blast_85(){
 echo -e "------ blast 85 ---------"
  ls | grep -v fasta| while read filename;do
    extension="${filename%.*}"  
    mkdir -p output
    echo $filename blastn -perc_identity 85  -outfmt
    blastn -db $hpv_database -query "$filename" -perc_identity 85  -outfmt "6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore" -out output/"$extension.out"
  done
}

blast_85_div()
{
 echo -e "------ blast 85 div ---------"
 cd div
 ls|sort | while read filename;do
    extension="${filename%.*}"
    mkdir -p output
    echo -e "blastn for $filename"
    blastn -db /tmp/chimera/HPV_COMPLETE_GENOME_REFERENCES.fasta -query $filename -perc_identity 85  -outfmt "6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore" -out output/"$extension""_concoverageb.out"
    done 
}


get_only_first_line_unique()
{
  echo -e "------ get only unique line -----"
  cd output
  ls | grep -v fasta| while read filename;do
     mkdir -p unique
     extension="${filename%.*}"
     echo $extension
     if [ -s "$filename" ];then
       bit_score_first=$(cat $filename | head -1 | awk -F " " '{print $13}')
       count=1
       echo $count>/tmp/count_txt.log
       echo "--------------------------" 
       echo "Getting first value in order to compare BitScore"
       cat $filename|tail -n +2 | while read line;do
         bit_score=$(echo $line | awk -F " " '{print $13}')
         echo $bit_score  $bit_score_first
         if [ "$bit_score" == "$bit_score_first" ]; then
           count=$((count+1)) 
           echo $count>/tmp/count_txt.log
         else 
           count=1
           break
         fi
         
       done
     fi  
     count=$(cat /tmp/count_txt.log)
     rm -rf /tmp/count_txt.log
     echo -e "Number of Evalue found" $count
     cat $filename | head -$count >unique/"$extension"_unique.out
     echo -e "Unique head for $filename"
  done
}

check_coverage()
{ 
  echo -e "----- Check Similarity -----"
  cd unique
  ls | grep -v fasta| while read filename;do
    extension="${filename%.*}"
    extension="${extension%_*}"
    if [ -s "$filename" ] 
    then
        value=$(cat $filename | awk -F " " {'print $4'})
	echo "----------- $extension has some data. Value: $value"
        echo "$value" "-lt" "$value_to_discriminate"  $filename
        if [ "$value" -lt "$value_to_discriminate" ]; then
          
          echo "$extension.fa -> Chimera: $value"
          echo "$extension.fa -> Chimera: The top hit alignment does not cover >60% of contigs sequence ">>$result_file
        else
          echo "$extension.fa -> NO Chimera: $value"
        fi
   
    else
	echo "$extension.fa -> Chimera: The top hit alignment does not cover >60% of contigs sequence"
        echo "$extension.fa -> Chimera: The top hit alignment does not cover >60% of contigs sequence">>$result_file
        
    fi
  done
}

last_check()
{
  echo -e "--------- Last check ----------"
  cd unique
  list=($(ls | awk -F "_" '{print $1}'|sort -u | uniq))

  for i in "${list[@]}";do
    cd $root_folder/workspace/div/output/unique
    echo "---------"
    new_file=$(echo $i"_")
    echo $new_file
    HPV=$(cat $new_file* | awk -F " " '{print $2}' | sort -u | wc -l)
    HPV_type=$(cat $new_file* | awk -F " " '{print $2}' | sort -u)
    number_or_repetition=$(cat $new_file* | awk -F " " '{print $2}'|sort| uniq -c |sort -nr|head -1|awk -F " " '{print $1}')
    probable_hpv=$(cat $new_file* | awk -F " " '{print $2}'|sort| uniq -c |sort -nr|head -1|awk -F " " '{print $2}')
    if [ "$HPV" -ne "1" ]; then
       if [ $number_or_repetition -eq "3" ]; then
          cd $root_folder/workspace/output/unique
          HPV_type_original=$(cat $new_file*| awk -F " " '{print $2}')
          if [ "$probable_hpv" != "$HPV_type_original" ]; then
            echo "$i.fa -> Chimera: Contigs segments dont have same TopHit as the total contig">>$result_file
            echo "$i.fa -> Chimera-------$probable_hpv----$HPV_type_original"
          else
            cd $root_folder/workspace/div/output/unique
            four_column=$(cat $new_file* |grep $probable_hpv   |awk -F " " '{print $4}')
            three_column=$(cat $new_file* | grep $probable_hpv | awk -F " " '{print $3}' | sort)
            three_column_1=$(echo $three_column  | awk -F " " '{print $1}')
            three_column_2=$(echo $three_column  | awk -F " " '{print $2}')
            three_column_3=$(echo $three_column  |  awk -F " " '{print $3}')

            four_column_1=$(echo $four_column  | awk -F " " '{print $1}')
            four_column_2=$(echo $four_column  | awk -F " " '{print $2}')
            four_column_3=$(echo $four_column  | awk -F " " '{print $3}')
            echo $i "----"$four_column_1 "------------ " $four_column_2 "------------ " $four_column_3 "------------ "
            if [ "$four_column_1" -lt "70" ] || [ "$four_column_2" -lt "70" ] || [ "$four_column_3" -lt "70" ]; then
             echo "$i.fa -> Chimera: At least one segments sequence shows <70% coverage">>$result_file
             echo "$i.fa -> Chimera-------$probable_hpv----$HPV_type_original"
           else
             echo $three_column
             cal=$(echo "$three_column_3 $three_column_1" | awk '{print $1-$2}')
             echo $cal
             if [ 1 -eq "$(echo "${three_column_1} < ${val_identity_ninety}" | bc)" ] && [ 1 -eq "$(echo "${cal} > ${val_identity_five}" | bc)" ]; then
               echo "$i.fa -> Chimera: Contigs segment show <90% identity and >5% difference with another segment">>$result_file
               echo "$i.fa -> Chimera-------$probable_hpv----$HPV_type_original"

             else
               echo "$i.fa -> $HPV_type_original, no chimera detected">>$result_file
               echo "$i.fa -> $HPV_type_original, no chimera detected"
             fi
           fi
        fi
      else
         echo "$i.fa -> Chimera: Contigs segments dont have same TopHit">>$result_file
         echo "$i.fa -> Chimera"
       fi
    elif [ "$HPV" -eq "1" ]; then
      cd $root_folder/workspace/output/unique
      HPV_type_original=$(cat $new_file*| awk -F " " '{print $2}')
      if [ "$HPV_type" != "$HPV_type_original" ]; then
         echo "$i.fa -> Chimera: Contigs segments dont have same TopHit as the total contig">>$result_file
         echo "$i.fa -> Chimera-------$HPV_type----$HPV_type_original"
      else
         cd $root_folder/workspace/div/output/unique
         four_column=$(cat $new_file* | awk -F " " '{print $4}')
         three_column=$(cat $new_file* | awk -F " " '{print $3}' | sort)
         three_column_1=$(echo $three_column | awk -F " " '{print $1}')
         three_column_2=$(echo $three_column | awk -F " " '{print $2}')
         three_column_3=$(echo $three_column | awk -F " " '{print $3}')
        
         four_column_1=$(echo $four_column | awk -F " " '{print $1}')
         four_column_2=$(echo $four_column | awk -F " " '{print $2}')
         four_column_3=$(echo $four_column | awk -F " " '{print $3}')
         echo $i "----"$four_column_1 "------------ " $four_column_2 "------------ " $four_column_3 "------------ "
         if [ "$four_column_1" -lt "70" ] || [ "$four_column_2" -lt "70" ] || [ "$four_column_3" -lt "70" ]; then
           echo "$i.fa -> Chimera: At least one segments sequence shows <70% coverage">>$result_file
           echo "$i.fa -> Chimera-------$HPV_type----$HPV_type_original"
         else
           echo $three_column
           cal=$(echo "$three_column_3 $three_column_1" | awk '{print $1-$2}')
           echo $cal
           if [ 1 -eq "$(echo "${three_column_1} < ${val_identity_ninety}" | bc)" ] && [ 1 -eq "$(echo "${cal} > ${val_identity_five}" | bc)" ]; then
             echo "$i.fa -> Chimera: Contigs segment show <90% identity and >5% difference with another segment">>$result_file
             echo "$i.fa -> Chimera-------$HPV_type----$HPV_type_original"
           
           else 
             echo "$i.fa -> $HPV_type_original, no chimera detected">>$result_file
             echo "$i.fa -> $HPV_type_original, no chimera detected"
           fi
         fi
         
         
      fi
      

    fi
  done


}
help()
{ 
  echo "***********************************"
  echo -e "Chimera Check"
  echo -e "Three arguemtns: "
  echo -e "1. Assembled contigs to analyze in fasta file"
  echo -e "2. Folder to work"
  echo -e "3. HPV database fasta file."
  echo -e "''''''''''''''''''"

}
header()
{
  echo "***********************************"
  echo -e "HPV Chimera Check"
  echo "************************************"

}
end(){

echo "HPV Chimera Check Finished"
echo "-------------------------"
echo "Check result in: $folder_to_start/chimera/HPVchimera_results.txt"

}

if [ "$#" -ne 3 ]; then
    echo "Illegal number of parameters"
    help
    exit 0
fi



if [ "$1" == "-h" ] ; then
    echo "Usage: `basename $0` [-h]"
    help
    exit 0
fi

fasta_file=$1
folder_to_start=$2
hpv_database=$3
mkdir -p $folder_to_start/chimera
val_identity_ninety=90.000
val_identity_five=5.000
root_folder=$folder_to_start/chimera
value_to_discriminate=60
result_file=$folder_to_start/chimera/HPVchimera_results.txt
rm -rf $result_file
touch $result_file
if_fasta_file $fasta_file
if_fasta_file $hpv_database
cp $fasta_file $root_folder
cd $folder_to_start/chimera/
rm -rf workspace
mkdir -p workspace
cp $fasta_file workspace
cd workspace
header
extract_id_indivial_files
blast_85
get_only_first_line_unique
check_coverage
divide_in_three
blast_85_div
cd output
get_only_first_line_unique
last_check
end
