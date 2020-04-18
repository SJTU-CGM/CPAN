# CPAN: A pan-genome based method for cancer genomics study

## Introduction
Cancer genomics studies are essentially relied on the human reference genome. However, tens of million base-pairs are still missing in the current human reference genome sequence. Pan-genome analysis has been proved its powerful in plant and human et al. We expand the pan-genome analysis to study cancer genomes. 

The CPAN homepage is http://cgm.sjtu.edu.cn/cpan/CPAN.html.
## Installation
- Requirements
    + R 3.1 or later (https://www.r-project.org/)
        R is utilized for visualization and statistical tests in CPAN. Please install R first and make sure R and Rscript are under your PATH.
        Several R packages are needed: ggplot2, reshape2 et al.
    + perl 
- Installation procedures
    + Download the CPAN from github:
        ```
        git clone git@github.com:SJTU-CGM/CPAN.git
        ```
    + Comile:
        ```
        make
        ```
        You will find executable files `cpan` et al. in `bin/` directory.
    + Add `bin/` to PATH and add `lib/` to LD_LIBRARY_PATH by adding the following text to `~/.bash_profile`:
        ```
        export PATH=$PATH:/path/to/CPAN/bin/:
        export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/CPAN/lib/:
        export PERL5LIB=$PERL5LIB:/path/to/CPAN/lib/:
        ```
    + And run:
        ```
        source ~/bash_profile.
        ```
    + 
    + 