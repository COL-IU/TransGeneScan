#!/bin/bash
# $1 = reference; $2 = reads; $3 = ptt; $4 = TGSHome; $5 = n; $6 = k; $7 = t

if [ ! $# == 6 ]; then
        echo "Usage: pileline.sh <refFile(abs path)> <pairedReadsPrefix> <TGSHome> <n> <k> <t>"
        exit
else

	bwa index $1
	samtools faidx $1
	bwa aln -t $6 -n $4 -k $5 $1 $2_1.fastq > alignment.1.sai
	bwa aln -t $6 -n $4 -k $5 $1 $2_2.fastq > alignment.2.sai
	bwa sampe $1 alignment.1.sai alignment.2.sai $2_1.fastq $2_2.fastq > alignment.sam

	samtools view -f 83 -q 1 -bhS alignment.sam > alignment.p.1.bam
	samtools view -f 163 -q 1 -bhS alignment.sam > alignment.p.2.bam
	samtools merge alignment.p.bam alignment.p.1.bam alignment.p.2.bam
	samtools sort alignment.p.bam alignment.p.sorted
	samtools mpileup -f $1 alignment.p.sorted.bam > alignment.p.pileup
	awk '{ total += $4 } END { print total/NR }' alignment.p.pileup

	samtools view -f 99 -q 1 -bhS alignment.sam > alignment.n.1.bam
	samtools view -f 147 -q 1 -bhS alignment.sam > alignment.n.2.bam
	samtools merge alignment.n.bam alignment.n.1.bam alignment.n.2.bam 
	samtools sort alignment.n.bam alignment.n.sorted
	samtools mpileup -f $1 alignment.n.sorted.bam > alignment.n.pileup
	awk '{ total += $4 } END { print total/NR }' alignment.n.pileup

	${3}scripts/printBases.pl alignment.p.pileup > alignment.p.printBases
	${3}scripts/printBases.pl alignment.n.pileup > alignment.n.printBases 

	${3}scripts/transConsensus.pl alignment.p.printBases alignment.p.consensus
	${3}scripts/transConsensus.pl alignment.n.printBases alignment.n.consensus

	${3}scripts/getTranscripts.py alignment.p.consensus fwdTranscripts $1
	${3}scripts/getNegTranscripts.py alignment.n.consensus revTranscripts $1

	cat fwdTranscripts.fasta revTranscripts.fasta > transcripts.fasta

	rm alignment*
	rm fwdTranscripts.fasta revTranscripts.fasta
	rm ${1}.*

fi
