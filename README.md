Purpose of VcfHunter
====================

Scaffremodler regroup a severals programs which principal aims is to detect use paired reads to 
link different genomic regions. These main programs are accompanied with a severals others 
programs that can be used in complement to improve scaffold assemblies or to detect large 
structural variations between a reference sequence and a re-sequenced genome.
<br>

Installation
------------

All proposed tools described here are written in python and work on linux system
To install the tools:
1- unzip the folder
2- go to the bin directory and open the loca_programs.conf file
3- set the path to each programs required (listed below)
To run fully all programs are needed:
1- bowtie can be found at http://bowtie-bio.sourceforge.net/index.shtml
2- bowtie2 can be found at http://bowtie-bio.sourceforge.net/bowtie2/index.shtml
3- SortSam.jar, MarkDuplicates.jar and FilterSamReads.jar belong to Picard Tools and can be
found at http://broadinstitute.github.io/picard/
4- bwa can be found at http://bio-bwa.sourceforge.net
5- samtools can be found at http://www.htslib.org
6- circos-0.67-7 is required and can be found at http://circos.ca/software/download/circos/
perl, python and java are required. Biopython is also required.
7- bamgrepreads can be found at https://code.google.com/p/variationtoolkit/wiki/BamGrepReads

<br>

How to cite
-----------
Please cite either:

Martin G, Baurens F-C, Droc G, Rouard M, Cenci A, Kilian A, Hastie A, Doležel J, Aury J-M, Alberti A, et al. 2016. **Improvement of the banana “Musa acuminata” reference sequence using NGS data and semi-automated bioinformatics methods.** *BMC Genomics 17:1–12.* https://bmcgenomics.biomedcentral.com/articles/10.1186/s12864-016-2579-4

Guillaume Martin, Françoise Carreel, Olivier Coriton, Catherine Hervouet, Céline Cardi, Paco Derouault, Danièle Roques, Frédéric Salmon, Mathieu Rouard, Julie Sardos, Karine Labadie, Franc-Christophe Baurens, Angélique D’Hont. 2017. **Evolution of the banana genome (Musa acuminata) is impacted by large chromosomal translocations.** *Molecular Biology and Evolution* https://academic.oup.com/mbe/article/34/9/2140/3852084

<br>

Author of scripts
-----------

Guillaume Martin (CIRAD)
Paco Derouault

<br>

License
-----------
Licencied under GPLv3

<br>
