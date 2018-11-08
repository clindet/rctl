
<!-- README.md is generated from README.Rmd. Please edit that file -->

## ngsjs

[ngsjs](https://github.com/JhuangLab/ngsjs) is a command line scripts
library to facilitate the construction, distribution and executiton of
reproducible next-generation sequencing (NGS) data analysis workflows
including [snakemake](https://snakemake.readthedocs.io/en/stable/),
[nextflow](https://www.nextflow.io/),
[bpipe](https://github.com/ssadedin/bpipe) and
[WDL](https://github.com/openwdl/wdl).

This is a experimental project exploring the simple way to construct,
distribute and execute various next-generation sequencing (NGS) data
analysis workflows.

Using the [Common workflow language (CWL)](https://www.commonwl.org/) to
construct bioinformatics data analysis workflows is attractive and
valuable for improving the performance and reproducible of biological
data analysis works.

Now, the difficulties that need to be solved are list:

  - Lack of integration and unify of massive data analysis workflows.
  - Lack of the unified distribution platform for various data analysis
    workflows (e.g. snakemake, nextflow, Galaxy, etc.).
  - Reuse of workflows language codes (e.g. commands, input and output
    information) on other programming platform are still slighly
    complicated.
  - The readability and reusable will also decreased when massive Python
    and R codes mixed with the workflows language codes.

Here, we proposed that using [node](https://nodejs.org/en/) to
distribute the bioinformatics data analysis workflows. The creation,
update and upload of a node package are very simple. Well-tested and
high-performance distribution tools of node packages, such as
[npm](https://www.npmjs.com/) and [yarn](https://www.yarnpkg.com), are
providing the service for more than 831,195 node packages.

Besides, we esbalish a framwork to integrate various data analysis
workflows and command line scripts. The functional files used in this
framework as shown below:

  - rbashful: Built in the ngsjs package (or user custom script) to
    dynamically render and execute the env.toml commands and cli.yaml
    workflows.
  - cli.yaml: Workflow controlor (e.g. bashful supported style).
  - env.toml: Store the fields and values of input and output
    parameters; the function commands (Python, R or others) indexed by
    unique keys.
  - others: Other workflow files and scripts (such as snakemake,
    nextflow workflows, scripts). The workflow files and scripts is
    better to upload as the node package, and the deployment functions
    for required extra files.

**Scripts:**

  - rsession: Get ouput of `sessionInfo()` and
    `sessioninfo::session\_info()`

<!-- end list -->

``` bash
rsession -f 2
#> ─ Session info ───────────────────────────────────────────────────────────────
#>  setting  value                       
#>  version  R version 3.5.1 (2018-07-02)
#>  os       macOS  10.14.1              
#>  system   x86_64, darwin15.6.0        
#>  ui       X11                         
#>  language (EN)                        
#>  collate  en_US.UTF-8                 
#>  ctype    en_US.UTF-8                 
#>  tz       Asia/Shanghai               
#>  date     2018-11-09                  
#> 
#> ─ Packages ───────────────────────────────────────────────────────────────────
#>  package     * version date       lib source        
#>  assertthat    0.2.0   2017-04-11 [1] CRAN (R 3.5.0)
#>  cli           1.0.1   2018-09-25 [1] CRAN (R 3.5.0)
#>  crayon        1.3.4   2017-09-16 [1] CRAN (R 3.5.0)
#>  getopt        1.20.2  2018-02-16 [1] CRAN (R 3.5.0)
#>  optparse    * 1.6.0   2018-06-17 [1] CRAN (R 3.5.0)
#>  pacman      * 0.5.0   2018-10-22 [1] CRAN (R 3.5.0)
#>  rstudioapi    0.8     2018-10-02 [1] CRAN (R 3.5.0)
#>  sessioninfo   1.1.1   2018-11-05 [1] CRAN (R 3.5.1)
#>  withr         2.1.2   2018-03-15 [1] CRAN (R 3.5.0)
#> 
#> [1] /Library/Frameworks/R.framework/Versions/3.5/Resources/library
```

  - rinstall: Install R packages and
    [BioInstaller](https://github.com/JhuangLab/BioInstaller) resources
    using `install.packages()` and R packages `devtools`, `BiocManager`
    and `BioInstaller`

<!-- end list -->

``` bash
rinstall -h
#> Usage: /usr/local/bin/rinstall [options] [params]
#> Examples:
#> /usr/local/bin/rinstall -p ini
#> /usr/local/bin/rinstall ini,yaml
#> /usr/local/bin/rinstall -f 2 JhuangLab/ngstk
#> /usr/local/bin/rinstall -f 2 -e "force = TRUE, ref = 'develop'" JhuangLab/ngstk
#> /usr/local/bin/rinstall -f 3 ggtree; /usr/local/bin/rinstall rinstall -f 4 ggtree
#> /usr/local/bin/rinstall -f 5 -e "show.all.names=T"
#> /usr/local/bin/rinstall -f 5 -e "show.all.versions=T" db_annovar_avsnp
#> /usr/local/bin/rinstall -f 5 -e "download.dir='/tmp/avsnp', extra.list=list(buildver='hg19')" db_annovar_avsnp
#> 
#> 
#> Options:
#>  -v, --verbose
#>      Print extra output [default FALSE]
#> 
#>  -f FUNC, --func=FUNC
#>      Index or name of used function [e.g. install.packages (1), devtools::install_github (2), BiocManager::install (3), pacman::p_load (4), BioInstaller::install.bioinfo (5)].
#> 
#>  -p PKGS, --pkgs=PKGS
#>      Package or item names [e.g. ggplot2,stringr or JhuangLab/BioInstaller (mode is devtools::install_github)].
#> 
#>  -e EXTRA, --extra=EXTRA
#>       Extra parameters [e.g. ref='develop'].
#> 
#>  -d, --doc
#>      Print functions document
#> 
#>  -h, --help
#>      Show this help message and exit
```

  - rbashful: Using the GO program
    [bashful](https://github.com/wagoodman/bashful), yaml and toml and R
    scripts to stitch together commands and bash snippits and run them
    with a bit of
style.

![](https://raw.githubusercontent.com/wagoodman/bashful/master/doc/demo.gif)

``` bash
rbashful -h
#> Usage: /usr/local/bin/rbashful [options]
#> 
#> 
#> Options:
#>  -v, --verbose
#>      Print extra output [FALSE]
#> 
#>  -c CLI-YAML, --cli-yaml=CLI-YAML
#>      bashful used YAML file [cli.yaml]
#> 
#>  -t ENV-TOML, --env-toml=ENV-TOML
#>      TOML file stores environment variables [env.toml]
#> 
#>  -e EXTRA-LIST, --extra-list=EXTRA-LIST
#>      Need to replaced environment variables
#> 
#>  -p, --parse-cli-yaml
#>      Replace cli config keys [FALSE]
#> 
#>  -o OUTPUT-CLI-YAML, --output-cli-yaml=OUTPUT-CLI-YAML
#>      Output file of parsed cli YAML file [*.parsed.yaml]
#> 
#>  -n CMD-NAME, --cmd-name=CMD-NAME
#>      Run CMDs section using name [NULL]
#> 
#>  --auto-create-dir
#>      Auto create dir in env.toml output section [FALSE]
#> 
#>  -h, --help
#>      Show this help message and exit
```

## Requirements

  - [node](https://nodejs.org/en/)
  - [R](https://cran.r-project.org/)
  - [GO](https://golang.org/)

### R packages

  - optparse
  - configr
  - devtools
  - BiocManager
  - sessioninfo

## Installation

You need to install the [node](https://nodejs.org/en/),
[R](https://cran.r-project.org/) and [GO](https://golang.org/) to run
all functions of [ngsjs](https://github.com/JhuangLab/ngsjs).

``` bash
npm install -g ngsjs
# or
yarn global add ngsjs

# If you not to globaly install ngsjs 
# You need to set the PATH
echo "export NGSJS_ROOT=/path/node_modules/nodejs" >> ~/.bashrc
echo "export PATH=$PATH:${NGSJS_ROOT}/bin" >> ~/.bashrc

# Current dir is /path
npm install ngsjs
# or
yarn add ngsjs
```

## Usage

Before using the `ngsjs` command line tools, you should run the
`resolve_r_deps` command to install the extra R packages used in `ngsjs`
scripts.

``` bash
# install the extra R packages used in `ngsjs` scripts
resolve_r_deps
```

### rsession

``` bash
# Print commandline help of rsession
rsession -h

# Print rsession R document (Just like ?sessionInfo in R client)
rsession -d

# Print R sessionInfo()
# The followed three lines are equivalent.
rsession
rsession -f 1
rsession -f sessionInfo

# The followed two lines are equivalent.
rsession -f 2 -e 'include_base=TRUE'
rsession -f sessioninfo::session_info -e 'include_base=TRUE'
```

### rinstall

``` bash
# Print commandline help of rinstall
rinstall -h

# Print rinstall R document (Just like ?sessionInfo in R client)
rinstall -d

# Install CRAN R package yaml (install.package)
rinstall yaml
rinstall -f 1 yaml

# Install R package ngstk from GitHub JhuangLab/ngstk (devtools::install_github)
rinstall -f 2 JhuangLab/ngstk

# Install R package ngstk from GitHub JhuangLab/ngstk (install.package)
# devtools::install_github with force = TRUE, ref = 'develop'
rinstall -f 2 -e "force = TRUE, ref = 'develop'" JhuangLab/ngstk

# Install Bioconductor package ggtree (BiocManager)
# BiocManager::install('ggtree')
rinstall -f 3 ggtree

# Install R packages (pacman)
# pacman::p_load(ggtree)
rinstall -f 4 ggtree

# Install and download BioInstaller resources

# Show all BioInstaller default resources name
rinstall -f BioInstaller::install.bioinfo -e "show.all.names=T"
rinstall -f 5 -e "show.all.names=T"

# Show ANNOVAR refgene and avsnp versions
rinstall -f BioInstaller::install.bioinfo -e "show.all.versions=T" db_annovar_refgene
rinstall -f 5 -e "show.all.versions=T" db_annovar_avsnp

# Show ANNOVAR hg19 refgene and avsnp
rinstall -f 5 -e "download.dir='/tmp/refgene', extra.list=list(buildver='hg19')" db_annovar_refgene
rinstall -f 5 -e "download.dir='/tmp/avsnp', extra.list=list(buildver='hg19')" db_annovar_avsnp
```

### rbashful

View a rbashful demo
[here](https://github.com/JhuangLab/ngsjs/tests/rbashful/rnaseq_splicing).

``` r
source_dir <- "~/Documents/repositories/ljf/github/ngsjs/tests/rbashful/rnaseq_splicing"

# View the cli.yaml
cat(paste0(readLines(sprintf("%s/cli.yaml", source_dir)), 
           collapse = "\n"), sep = "\n")
#> config:
#>     log-path: '{{logfn}}'
#> x-reference-data:
#>     all-apps: &ids
#>         - '{{id}}'
#> tasks:
#>     - name: Download genomes (hg38)
#>       tags: download_hg38
#>       parallel-tasks: 
#>         - cmd: rbashful -e "genome_version='94'" -n download_hg38_reffa
#> 
#>     - name: Download genomes (hg19)
#>       tags: download_hg19
#>       parallel-tasks: 
#>         - cmd: rbashful -e "genome_version='75'" -n download_hg19_reffa
#> 
#>     - name: STAR Indexing 
#>       tags: index
#>       parallel-tasks: 
#>         - cmd: rbashful -e "genome_version='75'" -n hg19_star_2_5_3a_rerffa_index
#>         - cmd: rbashful -e "genome_version='94'" -n hg38_star_2_5_3a_rerffa_index
#> 
#>     - name: Clean STAR Indexing 
#>       tags: clean
#>       parallel-tasks: 
#>         - cmd: rbashful -e "genome_version='75'" -n clean_star_2_5_3a_rerffa_index
#>         - cmd: rbashful -e "genome_version='94'" -n clean_star_2_5_3a_rerffa_index
#> 
#>     - name: STAR alignment hg19 and hg38
#>       tags: star_alignment
#>       parallel-tasks:
#>         - cmd: rbashful -e 'genome_version="75", id="<replace>"' & rbashful -e 'genome_version="94", id="<replace>"'
#>           for-each: *ids

# View the env.toml
cat(paste0(readLines(sprintf("%s/env.toml", source_dir)), 
           collapse = "\n"), sep = "\n")
#> title = "Environment variables for genome index (STAR 2.5.3a)"
#> [input]
#> star_bin="/opt/bin/aligner/STAR/STAR-2.5.3a/bin/Linux_x86_64/STAR"
#> hg19_dir="/u4/jhuangdata/reference/ensembl/75"
#> hg38_dir="/u4/jhuangdata/reference/ensembl/94"
#> hg19="!!glue {config$input$hg19_dir}/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa"
#> hg38="!!glue {config$input$hg38_dir}/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
#> hg19_gtf="!!glue {config$input$hg19_dir}/Homo_sapiens.GRCh37.75.gtf"
#> hg38_gtf="!!glue {config$input$hg38_dir}/Homo_sapiens.GRCh38.94.gtf"
#> 
#> [output]
#> star_index_out_dir="/u4/jhuangdata/reference/ensembl/{{genome_version}}/star_index"
#> 
#> [cmds]
#> 
#> # CMDs to download hg19 and hg38 genome
#> download_hg19_reffa = "!!glue cd {config$input$hg19_dir} && wget ftp://ftp.ensembl.org/pub/release-75/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa.gz ftp://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz"
#> download_hg38_reffa = "!!glue cd {config$input$hg38_dir} && wget ftp://ftp.ensembl.org/pub/release-94/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz ftp://ftp.ensembl.org/pub/release-94/gtf/homo_sapiens/Homo_sapiens.GRCh38.94.gtf.gz"
#> 
#> # CMDs to generate hg19 and hg38 genome index
#> hg19_star_2_5_3a_rerffa_index = "!!glue {config$input$star_bin} --runMode genomeGenerate --genomeDir {config$output$star_index_out_dir} --genomeFastaFiles {config$input$hg19} --sjdbGTFfile {config$input$hg19_gtf} --sjdbOverhang 100 --runThreadN 30"
#> hg38_star_2_5_3a_rerffa_index = "!!glue {config$input$star_bin} --runMode genomeGenerate --genomeDir {config$output$star_index_out_dir} --genomeFastaFiles {config$input$hg38} --sjdbGTFfile {config$input$hg38_gtf} --sjdbOverhang 100 --runThreadN 30"
#> 
#> # CMD to clean genome index
#> clean_star_2_5_3a_rerffa_index = "!!glue rm -r {config$output$star_index_out_dir}/*"

# View the submit.sh
cat(paste0(readLines(sprintf("%s/submit.sh", source_dir)), 
           collapse = "\n"), sep = "\n")
#> #!/usr/bin/env bash
#> 
#> id=12
#> workdir=~/Documents/repositories/ljf/github/ngsjs/tests/rbashful/rnaseq_splicing
#> # Run the default cmd in env.toml
#> rbashful -c ${workdir}/cli.yaml \
#>     --env-toml ${workdir}/env.toml \
#>     --cmd-name default \
#>     -v
#> 
#> # Parse the bashful YAML file (Output cli.parsed.yaml) 
#> # and then run the default cmd in env.toml
#> rbashful -c ${workdir}/cli.yaml \
#>     --env-toml ${workdir}/env.toml \
#>     --extra-list "id=${id},logfn=\"${id}.log\"" \
#>     --parse-cli-yaml \
#>     -v
#> 
#> # Running the parsed yaml (given tags cmds)
#> rbashful -c ${workdir}/cli.parsed.yaml \
#>     --env-toml ${workdir}/env.toml \
#>     --cmd-name start_tags \
#>     --extra-list "tags='star_alignment'" -v
```

## Maintainer

[Jianfeng Li](https://github.com/Miachol)
