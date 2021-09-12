# delta_probe
analyzes DNA basepair sequences to find an approximate 'best-fit' probe to help with strain identification.

Takes a formatted .txt holding full DNA sequences of two different types and uses hashing to analyze and pick the best fitting probe of K characters.
Currently this is setup to find a 100 character probe that identifies delta variant strains from regular covid-19 strains.  The best fit probe is determined based on
an error rate which is calculated by determining the probe's false positive and false negative rates.

Future plans:
Currently this can only check for exact sequence matching.  Future plans include allowing up to any 1 character mismatch between the probe and sequence to still be considered a 'match'.
project will be updated to be more efficient by adding a polynomial rolling hash function instead of using substr()
