#! /bin/sh

# Run tantan with various inputs and check the outputs.

cd $(dirname $0)
PATH=$PATH:../src

countLowercaseLetters () {
    grep -v '^>' "$@" | tr -cd a-z | wc -c
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
} |
diff tantan_test.out -
