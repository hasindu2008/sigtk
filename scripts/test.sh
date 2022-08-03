#!/bin/bash

# terminate script
die() {
	echo "$1" >&2
	echo
	exit 1
}

if [ "$1" = 'mem' ]; then
    mem=1
else
    mem=0
fi

ex() {
    if [ $mem -eq 1 ]; then
        valgrind --leak-check=full --error-exitcode=1 "$@"
    else
        "$@"
    fi
}

echo "DNA sref"
ex  ./sigtk sref test/nCoV-2019.reference.fasta > a.tsv || die "Running the tool failed"
diff -q test/sref_dna.exp a.tsv || die "diff failed"

echo "RNA sref"
ex ./sigtk sref test/rnasequin_sequences_2.4.fa --rna > a.tsv || die "Running the tool failed"
diff -q test/sref_rna.exp a.tsv  || die "diff failed"

echo "DNA seg"
ex ./sigtk seg test/sp1_dna.blow5 > test/seg_dna.txt || die "Running the tool failed"
diff test/seg_dna.txt test/seg_dna.exp

echo "RNA seg"
ex ./sigtk seg test/sequin_rna.blow5 > test/seg_rna.txt || die "Running the tool failed"
diff -q test/seg_rna.txt test/seg_rna.exp  || die "diff failed"

echo "DNA jnn"
ex ./sigtk jnn test/sp1_dna.blow5 > test/jnn_dna.txt || die "Running the tool failed"
diff test/jnn_dna.txt test/jnn_dna.exp

echo "RNA jnn"
ex ./sigtk jnn test/sequin_rna.blow5 > test/jnn_rna.txt || die "Running the tool failed"
diff -q test/jnn_rna.txt test/jnn_rna.exp  || die "diff failed"

echo "Tests passed"