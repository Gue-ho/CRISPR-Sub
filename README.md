# CRISPR-Sub

CRISPR-Cas9 induces DNA cleavages at desired target sites in a guide RNA-dependent manner; DNA editing occurs through the resulting activity of DNA repair processes including non-homologous end joining (NHEJ), which is dominant in mammalian cells. NHEJ repair frequently causes small insertions and deletions (indels) near DNA cleavage sites but only rarely causes nucleotide substitutions. High-throughput sequencing is the primary means of assessing indel and substitution frequencies in bulk populations of cells in the gene editing field. However, it is difficult to detect bona fide substitutions, which are embedded among experimentally-induced substitution errors, in high-throughput sequencing data. Here, we developed a novel analysis method, named CRISPR-Sub, to statistically detect Cas9-mediated substitutions in high-throughput sequencing data by comparing Mock- and CRISPR-treated samples.


CIRSPR-Sub needs python3 and below module in python3:
    
    xlsxwriter, scipy, numpy



# Usage

CRISPR-Sub can run with:

    python3 nt_substitution.py {reference sequence} {target sequence} {fastqjoin file 1} {fastqjoin file 2} {output file}

test code:

    python3 nt_substitution.py aggggataccaccgatctctgtgatctgcgactgttttctctgtctgtgcaggtccacagtatggcattgcccgtgaagatgtggtcctgaatcgtattcttggggaaggcttttttggggaggtctatgaaggtgtctacacaaatcatgtgagttctaggatcttcccttacactcctcttccacatgtctgtagggtgagacagagctcgaa GGTCCTGAATCGTATTCTTGggg test.fastqjoin test_con.fastqjoin output

# License
-------
CRISPR-Sub is licensed under the new BSD licence.

Copyright (c) 2019, Gue-ho Hwang and Sangsu Bae
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

* Neither the name of the Hanyang University nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNERS OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
