# HPVChimera

This script aims to check putative human papillomavirus (HPV) contigs for chimera presence.

Requirements for the script:
Fasta file with contigs´ sequences
HPVdatabase
Nucleotide-Nucleotide BLAST 2.10.1+

You will need to make the database for BLAST:
makeblastdb -in HPV_COMPLETE_GENOMES.fasta -parse_seqids -blastdb_version 5  -title "HPVcomplete" -dbtype nucl

This script requires three parameters as input.
1. Fasta file with your contigs´s sequences
2. Workspace folder
3. HPVdatabase (HPV_COMPLETE_GENOMES.fasta)

Example

./HPVChimera.sh contigssequences.fasta /$workspacefolder HPVdatabase.fasta

After being executed, you will see the message "Check result in: /$workspacefolder/chimera/result.txt"

Summary of steps:

The putative HPV contig´s sequences are compared to a database of known HPV sequences with BLAST.
The database (HPVdatabase) comprises HPV sequences from all oficially established HPV types present in the International HPV Reference Center database (hpvcenter.se) together with all sequences from non-oficially established HPV types found in the PaVe database (https://pave.niaid.nih.gov/).
Steps that are performed:
1. Check if the contig´s sequence shows at least 85% sequence identity to any of the HPV types present in the database. We decided to set the cut-off at 85% identity, to give a slight margin for the 90% homology within the HPV L1 gene required for 2 sequences to be considered the same HPV type. 
2. Check if the contig´s sequence shows at least 60% aligning coverage to the top hit HPV type.
3. Contigs´s sequences are divided in 3 equal segments. All three segments´s sequences from every contig are compared to the same database of known HPV sequences and chimeras are reported if: a) the segments do not share the same top hit when being compared to the HPV database,  b) the top hit obtained from the three segments is not the same as the top hit obtained from the corresponding (total) contig when being compared to the HPV database, c) the alignment coverage is <70% for any of the 3 segments and, d) if at least one of the segments has less than 90% similarity, and at least one of the segments has more than 90% similarity, and if the difference between these segments in terms of similarity to corresponding overlapping parts is more than 5% (for example if segment 1 was 88% similar and segment 2 was 94% similar).



Example of output:(HPVchimera_results.txt)

Name of 1st contig:     HPV85, no chimera detected.
Name of 2nd contig:     Chimera: Contig sequence does not show >85% sequence identity to any HPV type given.
Name of 3rd contig:     Chimera: The top hit alignment does not cover >60% of contigs sequence.
Name of 4th contig:     Chimera: Contigs segments dont have same TopHit.
Name of 5th contig:     Chimera: Contigs segments dont have same TopHit as the total contig.
Name of 6th contig:     Chimera: At least one segments sequence shows <75% coverage.
Name of 7th contig:     Chimera: Contigs segment show <90% identity and >5% difference with another segment.
