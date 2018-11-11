
<!-- README.md is generated from README.Rmd. Please edit that file -->

<p align="center">
  <a href="https://github/JhuangLab/ngsjs">
    <img
      alt="ngsjs"
      src="doc/images/ngsjs-logo.svg"
      width="400"
    />
  </a>
</p>

<p align="center">
   <a href="https://circleci.com/gh/JhuangLab/ngsjs/tree/dev"><img src="https://img.shields.io/circleci/project/github/JhuangLab/ngsjs/dev.svg" alt="Build Status"></a>
  <a href="https://npmcharts.com/compare/ngsjs?minimal=true"><img src="https://img.shields.io/npm/dm/ngsjs.svg" alt="Downloads"></a>
  <a href="https://www.npmjs.com/package/ngsjs"><img src="https://img.shields.io/npm/v/ngsjs.svg" alt="Version"></a>
  <a href="https://www.npmjs.com/package/ngsjs"><img src="https://img.shields.io/npm/l/ngsjs.svg" alt="License"></a>
</p>

[ngsjs](https://github.com/JhuangLab/ngsjs) is a command line scripts
library to facilitate the construction, distribution and execution of
reproducible next-generation sequencing (NGS) data analysis workflows
including [snakemake](https://snakemake.readthedocs.io/en/stable/),
[nextflow](https://www.nextflow.io/),
[bpipe](https://github.com/ssadedin/bpipe) and
[WDL](https://github.com/openwdl/wdl).

# ngsjs

This is an experimental project exploring a simple way to construct,
distribute and execute various next-generation sequencing (NGS) data
analysis workflows.

Using the [Common workflow language (CWL)](https://www.commonwl.org/) to
construct bioinformatics data analysis workflows is attractive and
valuable for improving the performance and reproducible of biological
data analysis works.

Now, there are several difficulties for developing bioinformatics
workflow that needs to be solved:

  - Lack of integration and unify of massive data analysis workflows.
  - Lack of the unified distribution platform for various data analysis
    workflows (e.g. snakemake, nextflow, Galaxy, etc.).
  - Reuse of workflows language codes (e.g. commands, input and output
    information) on other programming platform are still complicated.
  - The readability and reusable will also be decreased when massive
    Python and R codes mixed with the workflows language codes.

Here, we proposed that using [node](https://nodejs.org/en/) to
distribute the bioinformatics data analysis workflows. The creation,
update and upload of a node package are very simple. Well-tested and
high-performance distribution tools of node packages, such as
[npm](https://www.npmjs.com/) and [yarn](https://www.yarnpkg.com), are
providing the service for more than 831,195 node packages.

Besides, we establish a framework to integrate various data analysis
workflows and command line scripts. The functional files used in this
framework as shown below:

  - rbashful: Built in the ngsjs package (or user custom script) to
    dynamically render and execute the env.toml commands and cli.yaml
    workflows.
  - cli.yaml: Workflow controller (e.g. bashful supported style).
  - env.toml: Store the fields and values of input and output
    parameters; the function commands (Python, R or others) indexed by
    unique keys.
  - others: Other workflow files and scripts (such as snakemake,
    nextflow workflows, scripts). We recommend users to upload the
    workflow/script files as the node package, and provide the download
    functions for required extra files.

**Scripts:**

  - rdeps: Install `ngsjs` required R packages
  - rsession: Get output of `sessionInfo()` and
    `sessioninfo::session\_info()`
  - rinstall: Install R packages and
    [BioInstaller](https://github.com/JhuangLab/BioInstaller) resources
    using `install.packages()` and R packages `devtools`, `BiocManager`
    and `BioInstaller`
  - rbashful: Using the GO program
    [bashful](https://github.com/wagoodman/bashful), yaml and toml and R
    scripts to stitch together commands and bash snippits and run them
    with a bit of style
  - rconfig: Using the R package `configr` to parse and generate json,
    ini, yaml, and toml format configuration files

## Requirements

  - [node](https://nodejs.org/en/)
  - [R](https://cran.r-project.org/)
  - [GO](https://golang.org/)

### R packages

  - optparse
  - configr
  - devtools
  - pacman
  - BiocManager
  - sessioninfo

## Installation

You need to install the [node](https://nodejs.org/en/),
[R](https://cran.r-project.org/) and [GO](https://golang.org/) for
running all [ngsjs](https://github.com/JhuangLab/ngsjs) executable
files.

``` bash
npm install -g ngsjs
# or
curl --compressed -o- -L https://yarnpkg.com/install.sh | bash
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

Before try your `ngsjs` command line tools, you need run the `rdeps`
getting all the extra R packages required by `ngsjs`.

``` bash
# install the extra R packages used in `ngsjs` scripts
rdeps
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

### rconfig

``` r
# Use configr::read.config parsing json format file
# Reture the list object output
rconfig package.json
rconfig -i package.json

# Use configr::read.config parsing json format file with the custom R function
rconfig -i test.json -r 'function(x){x[["a"]] + x[["b"]]}'
rconfig -i test.json -r 'function(x){x[["a"]]}'
rconfig -i test.json -r 'function(x){x[["b"]]}'

# Use configr::write.config parsing json format file
rconfig -f "configr::write.config" test.json -e "config.dat=list(a=1, b=2), write.type='json'"
rconfig -f 2 test.json -e "config.dat=list(a=1, b=2), write.type='json'"
```

### rbashful

[bashful](https://github.com/wagoodman/bashful) is a GO program and used
by `rbashful`, so you need to install it before use the
`rbashful`.

**Ubuntu/Debian**

``` bash
wget https://github.com/wagoodman/bashful/releases/download/v0.0.10/bashful_0.0.10_linux_amd64.deb
sudo apt install ./bashful_0.0.10_linux_amd64.deb
```

**RHEL/Centos**

``` bash
wget https://github.com/wagoodman/bashful/releases/download/v0.0.10/bashful_0.0.10_linux_amd64.rpm
rpm -i bashful_0.0.10_linux_amd64.rpm
```

**Mac**

``` bash
brew tap wagoodman/bashful
brew install bashful
```

or download a Darwin build from the releases page.

**Go
tools**

``` bash
go get github.com/wagoodman/bashful
```

![](https://raw.githubusercontent.com/wagoodman/bashful/master/doc/demo.gif)

View a `rbashful` demo
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

## Snippets of ngsjs scripts output

### rdeps

``` bash
rdeps
#> INFO [2018-11-09 22:01:54] All basic dependences (R packages) were resolved.
#> INFO [2018-11-09 22:01:54] optparse, configr, devtools, BiocManager, sessioninfo, BioInstaller, glue, futile.logger
```

### rinstall

``` bash
rinstall
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
#> Description:
#> rinstall is an R-based tool to install or download R packages and other resources supported by R package BioInstaller.
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

### rbashful

``` bash
rbashful
#> INFO [2018-11-09 22:01:55] No commands were ran.
#> Usage: /usr/local/bin/rbashful [options] [params]
#> Examples:
#> /usr/local/bin/rbashful -c ${workdir}/cli.yaml --env-toml ${workdir}/env.toml --cmd-name default -v
#> Description:
#> rbashful is an extend bashful tool for style bash commands.
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

### rconfig

``` bash
rconfig
#> Usage: /usr/local/bin/rconfig [options] [params]
#> /usr/local/bin/rconfig package.json
#> /usr/local/bin/rconfig -i package.json
#> /usr/local/bin/rconfig -f 2 test.json -e "config.dat=list(a=1, b=2), write.type='json'"
#> /usr/local/bin/rconfig -f "configr::write.config" test.json -e "config.dat=list(a=1, b=2), write.type='json'"
#> /usr/local/bin/rconfig -i test.json -r 'function(x){x[["a"]] + x[["b"]]}'
#> /usr/local/bin/rconfig -i test.json -r 'function(x){x[["a"]]}'
#> /usr/local/bin/rconfig -i test.json -r 'function(x){x[["b"]]}'
#> Description:
#> rconfig is an R-based tool to parse and generate configuration file.
#> 
#> Options:
#>  -v, --verbose
#>      Print extra output [default FALSE]
#> 
#>  -f FUNC, --func=FUNC
#>      Index or name of used function [e.g. configr::read.config (1), configr::write.config].
#> 
#>  -i CONFIG, --config=CONFIG
#>      Input or output configuationo file.
#> 
#>  -e EXTRA, --extra=EXTRA
#>      Extra parameters [e.g. ref='develop'].
#> 
#>  -r RFUNC, --rfunc=RFUNC
#>      R function to process the parsed configuation file [e.g. function(x){return(x[[1]])} ].
#> 
#>  -d, --doc
#>      Print functions document
#> 
#>  -h, --help
#>      Show this help message and exit
```

### rsession

``` bash
rsession -f 2 -e 'include_base=TRUE'
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
#>  ! package     * version date       lib source        
#>    assertthat    0.2.0   2017-04-11 [1] CRAN (R 3.5.0)
#>    base        * 3.5.1   2018-07-05 [?] local         
#>    cli           1.0.1   2018-09-25 [1] CRAN (R 3.5.0)
#>  P compiler      3.5.1   2018-07-05 [1] local         
#>    crayon        1.3.4   2017-09-16 [1] CRAN (R 3.5.0)
#>  P datasets    * 3.5.1   2018-07-05 [1] local         
#>    getopt        1.20.2  2018-02-16 [1] CRAN (R 3.5.0)
#>  P graphics    * 3.5.1   2018-07-05 [1] local         
#>  P grDevices   * 3.5.1   2018-07-05 [1] local         
#>  P methods     * 3.5.1   2018-07-05 [1] local         
#>    optparse    * 1.6.0   2018-06-17 [1] CRAN (R 3.5.0)
#>    pacman      * 0.5.0   2018-10-22 [1] CRAN (R 3.5.0)
#>    rstudioapi    0.8     2018-10-02 [1] CRAN (R 3.5.0)
#>    sessioninfo   1.1.1   2018-11-05 [1] CRAN (R 3.5.1)
#>  P stats       * 3.5.1   2018-07-05 [1] local         
#>  P utils       * 3.5.1   2018-07-05 [1] local         
#>    withr         2.1.2   2018-03-15 [1] CRAN (R 3.5.0)
#> 
#> [1] /Library/Frameworks/R.framework/Versions/3.5/Resources/library
#> 
#>  P ── Loaded and on-disk path mismatch.
```

## How to contribute?

Please fork the [GitHub ngsjs
repository](https://github.com/JhuangLab/ngsjs), modify it, and submit a
pull request to us.

## Maintainer

[Jianfeng Li](https://github.com/Miachol)

## License

MIT
