# SSR
Sparse low-rank separated representation models for learning from data

Commands to issue for compiling and running the codes (synthetic/real datasets).

(1) SYNTHETIC DATASETS

 (1.1) SSR-CD algorithm:
      make clean;
      make synth_ssrcd;
      ./synth_ssrcd

 (1.2) Block-SSR algorithm:
      make clean;
      make synth_bssr;
      ./synth_bssr

 (1.3) Standard ALS:
      make clean;
      make synth_als;
      ./synth_als

(2) Machine learning benchmark DATASETS

 (2.1) SSR-CD with splits of the data (used for comparisons with Yang et al. 2015):
      make clean;
      make real_ssrcd;
      ./real_ssrcd
       
 (2.2) SSR-CD without splits (simplified version):
      make clean;
      make real_single_ssrcd;
      ./real_single_ssrcd
      
      
  The regression datasets with test-train splits can be downloaded from https://people.orie.cornell.edu/andrew/code
