# delta_probe
analyzes DNA basepair sequences to find an approximate 'best-fit' probe to help with strain identification.

Takes a formatted .txt holding full DNA sequences of and uses a hash map to analyze and pick the best fitting probe of K characters.  For every length K probe, a one-character mismatch is allowed, meaning it is still considered a match if the probe has at most 1 mismatch along its entire length.

Currently this project is setup to find a 100 character probe that identifies delta variant strains from regular covid-19 strains.  Technically this should work with any text string as long as the input file is structured and labeled the same way as this project's input file.  The best fit probe is determined based on
an error rate which is calculated by determining the probe's false positive and false negative rates.

This project makes use of a polynomial rolling hash to speed up the comparisons.  However, running time will increase dramatically with larger and larger probe sizes.  
