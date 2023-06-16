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


handle_tests() {
	numfailed=$(wc -l < diff.txt)
	numcases=$(wc -l < ${ORIG})
	numres=$(wc -l < ${RES})
	echo "$numfailed of $numcases test cases deviated."
	missing=$(echo "$numcases-$numres" | bc)
	echo "$missing entries in the truthset are missing in the testset"
	failp=$(echo "$numfailed*100/$numcases" | bc)
	[ "$failp" -gt ${THRESH} ] && die "${1}: Validation failed"
	echo "Validation passed"
}

execute_test() {

	ORIG=$1
	RES=$2
	THRESH=$3
	diff -y --suppress-common-lines ${ORIG} ${RES} > diff.txt || handle_tests $testdir

}

echo "DNA sref"
ex  ./sigtk sref test/nCoV-2019.reference.fasta > a.tsv || die "Running the tool failed"
execute_test test/sref_dna.exp a.tsv 5 || die "diff failed"

echo "RNA sref"
ex ./sigtk sref test/rnasequin_sequences_2.4.fa --rna > a.tsv || die "Running the tool failed"
execute_test test/sref_rna.exp a.tsv 5 || die "diff failed"

echo "DNA prefix"
ex ./sigtk prefix test/sp1_dna.blow5 > test/prefix_dna.txt || die "Running the tool failed"
execute_test test/prefix_dna.txt test/prefix_dna.exp 5

echo "RNA prefix"
ex ./sigtk prefix test/sequin_rna.blow5 > test/prefix_rna.txt || die "Running the tool failed"
execute_test test/prefix_rna.txt test/prefix_rna.exp 20 || die "diff failed"

echo "DNA jnn"
ex ./sigtk jnn test/sp1_dna.blow5 > test/jnn_dna.txt || die "Running the tool failed"
#execute_test test/jnn_dna.txt test/jnn_dna.exp 5 || die "diff failed"

echo "RNA jnn"
ex ./sigtk jnn test/sequin_rna.blow5 > test/jnn_rna.txt || die "Running the tool failed"
execute_test test/jnn_rna.txt test/jnn_rna.exp 20 || die "diff failed"

echo "DNA event"
ex ./sigtk event test/sp1_dna.blow5 05d90f17-f4a6-4349-924c-3ffd3457a99d > test/event_dna.txt || die "Running the tool failed"
execute_test test/event_dna.txt test/event_dna.exp 5 || die "diff failed"

echo "RNA event"
ex ./sigtk event -c test/sequin_rna.blow5 > test/event_rna.txt || die "Running the tool failed"
execute_test test/event_rna.txt test/event_rna.exp 5 || die "diff failed"

echo "Tests passed"