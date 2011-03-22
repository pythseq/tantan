tantan
======

Introduction
------------

tantan is a tool for masking simple regions (low complexity and
short-period tandem repeats) in biological sequences.

The aim of tantan is to prevent false predictions when searching for
homologous regions between two sequences.  Simple repeats often align
strongly to each other, causing false homology predictions.

Setup
-----

Using the command line, go into tantan's "src" directory and type
"make".  This assumes you have a C++ compiler.

Usage
-----

* Suppose you have some nucleotide sequences (DNA or RNA) in a
  FASTA-format file called "ntseqs.fa".  You can identify simple
  repeats like this:

    tantan ntseqs.fa > masked.fa

  This will create a new FASTA file called "masked.fa" that replaces
  all masked regions with lowercase letters.  (tantan also works on
  FASTQ-format, though it does not use the quality data.)

* To mask proteins effectively, tantan needs to use different
  algorithm parameters than for nucleotides.  You have to tell it when
  you have proteins, using the "-p" option:

    tantan -p aaseqs.fa > masked.fa

  If you omit "-p" and the sequences look proteinaceous, tantan will
  print a warning message.

* By default, tantan indicates repetitive regions with lowercase
  letters.  You can make it replace repetitive letters with (say) "N"
  by using the "-x" option:

    tantan -x N ntseqs.fa > masked.fa

* By default, tantan does not preserve lowercase letters in the input
  sequences.  You can tell it to preserve them by using the "-c"
  option.  So the output will have the union of the lowercase in the
  input and the lowercase assigned by tantan:

    tantan -c ntseqs.fa > masked.fa

* tantan's masking rate is usually OK, but you can alter it by
  changing the "-r" parameter from its default value of 0.005.  Higher
  values increase the amount of masking, and lower values decrease it.
  This increases the masking rate:

    tantan -r 0.02 ntseqs.fa > masked.fa

* Finally, to mask extremely AT-rich DNA, you should change tantan's
  scoring matrix.  The "test" directory contains a matrix "atMask.mat"
  that works well for DNA with ~80% A+T, such as Plasmodium and
  Dictyostelium genomes.  We recommend masking such DNA like this:

    tantan -m atMask.mat -r 0.01 atrich.fa > masked.fa

The preceding examples cover all of tantan's options that you should
ever need.

Recommendations for homology search
-----------------------------------

1) Mask *both* (sets of) sequences.

2) If for some reason you wish to mask only one (set of) sequence(s),
   increase "-r" to 0.02 (0.05 for AT-rich DNA).

3) For DNA-versus-protein alignment, increase "-r" for the proteins to
   0.02.  If for some reason you mask only one (set of) sequence(s),
   make sure it's the proteins.

4) If you merely wish to find homologs, "hard masking" works well.
   Hard masking means replacing masked letters with a dummy letter,
   typically N for DNA and X for proteins.  Alternatively, some
   alignment tools have options to treat lowercase as hard-masked.

5) If you also wish to align the homologs, "soft masking" is
   attractive.  Soft-masking means that lowercase is treated as masked
   in earlier stages of the alignment algorithm but not later stages.
   However, it is important that the masking be lifted only *after*
   homology has been decided.  Most alignment tools cannot do this
   properly; FASTA and LAST can.

For more information, please read the tantan publication (see below).

Options
-------

-p  interpret the sequences as proteins
-x  letter to use for masking, instead of lowercase
-c  preserve uppercase/lowercase in non-masked regions
-m  file for letter pair scores (scoring matrix)
-r  probability of a repeat starting per position
-e  probability of a repeat ending per position
-w  maximum tandem repeat period to consider
-d  probability decay per period (period-(i+1) / period-i)
-a  gap existence cost
-b  gap extension cost
-s  minimum repeat probability for masking
-f  output type: 0=masked sequence, 1=repeat probabilities,
                 2=repeat counts, 3=BED
-h, --help  show help message, then exit
--version   show version information, then exit

Advanced issues
---------------

When tantan masks tandem repeats, it tends to leave the first
(left-most) repeat unit unmasked.  This sometimes allows us to find
homologs we would otherwise miss:

  TGCAAGCTA TTAGGCTTAGGTCAGTGC ttaagcttaggtcagtgc AACATA
  ||| ||| | |||||||||||||||||| ||| |||||||||||||| ||| ||
  TGCTAGCAA TTAGGCTTAGGTCAGTGC ttaggcttaggtcagtgc AACGTA

However, there is a danger of non-equivalent repeat units being
unmasked.  This happens especially if we mask DNA on one strand but
align it on the other strand:

                     TGCAAGCTA TTAGGCTTAGGTCAGTGC ttaagcttaggtcagtgc AACATA
                               ||||||||||||||||||
  TGCTAGCAA ttaggcttaggtcagtgc TTAGGCTTAGGTCAGTGC AACGTA

(My thanks to Junko Tsuji and Paul Horton for finding these issues.)

Miscellaneous
-------------

tantan is distributed under the GNU General Public License, either
version 3 of the License, or (at your option) any later version.  For
details, see COPYING.txt.

If you use tantan in your research, please cite:
"A new repeat-masking method enables specific detection of homologous
sequences", MC Frith, Nucleic Acids Research (in press).

tantan's website is: http://www.cbrc.jp/tantan/

If you have any questions, comments, or problems concerning tantan,
please email: tantan (ATmark) cbrc (dot) jp.
