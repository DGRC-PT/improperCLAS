# improperCLAS


Simple tool for identification and classification of improper pairs from SAM files.

The script starts by sampling 10,000 proper-pairs to calculate median and standard
deviation (SD) of the insert size (IS) (optionally, it can be indicated manually by the user).
Then, it categorizes the improper pairs as (i) interchromosomal (interchr) or chimeric, if the reads of the same pair are mapped in
different chromosomes; (ii) inversion (inv) improper-pairs, if both reads of the same
pair are mapped in the same orientation; (iii) deletion (del) for inward facing and (iv)
duplication (dup) for outward facing read pairs with an IS larger than the median
IS+3*IS SD, previously calculated.

For more information see + [David et al., 2020](https://pubmed.ncbi.nlm.nih.gov/32030560/) 


## Dependencies:
+ python3
+ python sys & statistics

## Usage:


In the command line:
<pre><code> python3 improperCLAS.py [SAMFILE] [IMPROPER PAIRS DELETIONS OUTPUT] [IMPROPER PAIRS DUPLICATIONS OUTPUT] [IMPROPER PAIRS INVERSIONS OUTPUT] [IMPROPER PAIRS TRANSLOCATIONS OUTPUT] [INSERT SIZE (optinal)] [INSERT SIZE STANDAR DEVIATION (optinal)]

</code></pre>


## License:

GPLv2


## Found a bug?

Or maybe just wanto to drop some feedback? Just open an issue on github