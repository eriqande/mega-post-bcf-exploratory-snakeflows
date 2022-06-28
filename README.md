mega-post-bcf-exploratory-snakeflows
================

I am just starting to put this together. The basic idea is that when
exploring NGS data that comes to you in a BCF file, there are lots of
ways you might want to tweak it in the process. For example:

-   Use different BCF files
-   Use different subsets of individuals
-   Use different filtering criteria
-   Thin the number of markers to a different level
-   Thin the number of markers to the same level, but do different
    randomly sampled replicates of that.

On top of that, once you start analyzing those data you will want to try
different things:

-   Different values of K with ngsAdmix
-   Different replicate runs of ngsAdmix
-   Different parameter settings of other packages
-   Some you will want to break over chromosomes to parallelize
-   Others, you can’t break over chromosomes, etc.

So, these are all just the things I am thinking about as I start
noodling around with a few things here.

I am thinking about a wildcarding/directory structure that looks
something like this for the basic BCF maneuvers:

``` python
"bcf_{bcf_file}/samp-sub_{sample_subset}/filt_{filter_crit}/tf_{thinning_factor}/tf-seed_{thin_factor_seed}/"
```
