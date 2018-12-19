#! /bin/sh

# Run tantan with various inputs and check the outputs.

cd $(dirname $0)

# Make sure we use this version of tantan:
PATH=../src:$PATH

countLowercaseLetters () {
    grep -v '^>' "$@" | tr -cd a-z | wc -c | tr -d ' '
}

{
    tantan hg19_chrM.fa
    echo
    tantan -p -xx titin_human.fa
    echo
    tantan -p -f2 titin_human.fa
    echo
    tantan -c hg19_chrM.fa | countLowercaseLetters
    echo
    tantan -w50 -e.005 -d1 hg19_chrM.fa | countLowercaseLetters
    echo
    tantan -m atMask.mat -r.01 hg19_chrM.fa | countLowercaseLetters
    echo
    tantan -p -a11 -b2 titin_human.fa | countLowercaseLetters
    echo
    tantan panda.fastq
    echo
    tantan -i2 -j3 hg19_chrM.fa | countLowercaseLetters
    echo
    tantan -f4 panda.fastq
    echo
    tantan -f4 -b12 hard.fa
    echo
    tantan -f4 -n1 hg19_chrM.fa
} 2>&1 |
diff -u tantan_test.out -
