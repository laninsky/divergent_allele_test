# divergent_allele_test v1.0.1
Are the combination of alleles in your population the most divergent you could get out of all the alleles in the total meta-population? 

#How it works#

This R script calculates average pairwise distance for all unique combinations of alleles in your metapopulation (for the number of alleles observed in your actual population). It then tells you whether your observed alleles are as divergent as could be expected given the metapopulation alleles. E.g. if your population has three alleles, are those three alleles the most diverse combination you could "get" from the metapopulation?

More mathematically, it compares the divergence of observed unique alleles compared to the divergence of all combinations of (n choose k) alleles, where n = total number of alleles in metapopulation, and k = number of observed alleles in each population.

#Things to note#

1) This program "eats" arlequin files that have a haplotype block at the top of the file. Check out the example data if you are unsure on this format.

2) You need to make sure the library stringr is loaded in to your R environment

3) The program will look for and consolidate any identical haplotypes

4) It counts indels as a state for the calculation of pairwise distance (so make sure sequences of different lengths are 'padded' by missing data, not dashes).

5) If there are haplotypes/alleles that do not reciprocally match (i.e. a haplotype/allele with missing data can match to more than one haplotype/allele, and the haplotypes/alleles it matches to do not match to each other... I'm just going to use the word 'haplotype' from here on out...), the program will identify the haplotype with the largest amount of missing data. In most cases this will be the culprit, but a better idea is to make sure all of your haplotypes/alleles are unique before running this program.

6) Pariwise distance will differ slightly from arlequin if there is missing data, as divergent_allele_test calculates the proportional difference between haplotypes only based on non-missing data, whereas arlequin averages the differences over the entire length of alignment to calculate proportional differences.

7) Pairwise distance will also differ slightly between the programs if there are indels, because Arlequin only counts indels as part of the haplotype alignment length if they are variable within the population in question, and divergent_allele_test counts them as part of the length if they are variable across all haplotype definitions (over all of your populations).

8) The program can handle spaces and tabs (or a mixture of these) as delimiters in your arlequin file. However, make sure that closing braces i.e. "}" do not occur on the same line as your data. This might give you a weird non-specific error.

9) For comparisons of divergence, the program uses only 10 decimal points due to floating point issues.

If arlequin can accept the file, and divergent_allele_test can't, then email alana dot alexander at ku dot edu. If you are having trouble with arlequin, make sure your hyphens haven't been made into em-dashes by word.

How to run it

The easiest way for people less familiar with R, is to paste the entire function into R. You can then call the function by:

divergent_allele_test(working_dir,file_name)

where:

working_dir == pathway to the folder with your arlequin file e.g. "C:/blahblahblah"

file_name == the name of your arlequin file (in haplotype list format) e.g. "data.txt"

Demo arlequin file "ATL_by_region_394.arp" located in example folder

#Suggested citation#

This code was first published in:
TBD

If you could cite the pub, and the progam as below, that would be lovely:

Alexander, A. 2015. divergent_allele_test v1.0.0. Available from https://github.com/laninsky/divergent_allele_test

This pipeline also wouldn't be possible without:

R: R Core Team. 2015. R: A language and environment for statistical computing. URL http://www.R-project.org/. R Foundation for Statistical Computing, Vienna, Austria. https://www.r-project.org/

Stringr:  Hadley Wickham (2012). stringr: Make it easier to work with strings..
  R package version 0.6.2. http://CRAN.R-project.org/package=stringr (for up-to-date citation information run citation("stringr" in R)

#Version history#

v1.0.1: changed language in the readme and R code to reflect we are calculating average pairwise distance over unique alleles (rather than saying nucleotide diversity seeing as this is not a population metric where we expect haplotypes to be identical).

v1.0.0: ready to rock 'n' roll


