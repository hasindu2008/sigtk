#!/bin/bash

# terminate script
die() {
	echo "$1" >&2
	echo
	exit 1
}

echo "DNA sref"
./sigtk sref test/nCoV-2019.reference.fasta > a.tsv || die "Running the tool failed"
diff -q test/sref_dna.exp a.tsv || die "diff failed"

echo "RNA sref"
./sigtk sref test/rnasequin_sequences_2.4.fa --rna > a.tsv || die "Running the tool failed"
diff -q test/sref_rna.exp a.tsv  || die "diff failed"

echo "Tests passed"