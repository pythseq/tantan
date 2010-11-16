tantan
======

Introduction
------------

tantan is a tool for masking simple regions (low complexity and
short-period tandem repeats) in biological sequences.

Setup
-----

Using the command line, go into tantan's "src" directory and type
"make".  This assumes you have a C++ compiler.

Usage
-----

* Suppose you have some nucleotide sequences (DNA or RNA) in a
  FASTA-format file called "ntseqs.fa".  You can mask them like this:

    tantan ntseqs.fa > masked.fa

  This will put the masked sequences in a new file called "masked.fa".

* To mask proteins, you need to use the "-p" option:

    tantan -p aaseqs.fa > masked.fa

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
ever need.  The algorithm has additional parameters that are listed in
the next section for completeness.

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
-f  output type: 0=masked sequence, 1=repeat probabilities, 2=repeat counts
-h, --help  show help message, then exit
--version   show version information, then exit

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
