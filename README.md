
<!-- README.md is generated from README.Rmd. Please edit that file -->

<p align="center">
  <a href="https://github/ngsjs/ngsjs">
    <img
      alt="ngsjs"
      src="doc/images/ngsjs-logo.svg"
      width="400"
    />
  </a>
</p>

<p align="center">
  <a href="https://www.npmjs.com/package/ngsjs"><img src="https://img.shields.io/badge/lifecycle-experimental-orange.svg" alt="Life cycle: experimental">
  <a href="https://circleci.com/gh/ngsjs/ngsjs/tree/master"><img src="https://img.shields.io/circleci/project/github/ngsjs/ngsjs/master.svg" alt="Build Status"></a>
  <a href="https://npmcharts.com/compare/ngsjs?minimal=true"><img src="https://img.shields.io/npm/dm/ngsjs.svg" alt="Downloads"></a>
  <a href="https://www.npmjs.com/package/ngsjs"><img src="https://img.shields.io/npm/v/ngsjs.svg" alt="Version"></a>
  <a href="https://www.npmjs.com/package/ngsjs"><img src="https://img.shields.io/npm/l/ngsjs.svg" alt="License"></a>
</p>

[ngsjs](https://github.com/ngsjs/ngsjs) is a set of command line tools,
NGS data analysis workflows \[[WDL](https://github.com/openwdl/wdl),
[Nextflow](https://www.nextflow.io/),
[snakemake](https://snakemake.readthedocs.io/en/stable/), and
[bpipe](https://github.com/ssadedin/bpipe)\], and R shiny plugins/R
markdown document for exploring next-generation sequencing data.

# ngsjs

Now, there are several difficulties for next-generation sequencing (NGS)
data analysis projects that needs to be solved:

  - Standardized project management, directory structured，recording and
    checking of raw data and analysis result, standardized logging for
    input, output and commands
  - Construction and redeployment of computing environment including all
    required tools, databases and other files.
  - Lack of integration and unify of massive data analysis workflows.
  - Lack of the unified distribution platform for various data analysis
    workflows (e.g. snakemake, nextflow, Galaxy, etc.).
  - Reuse of workflows language codes (e.g. commands, input and output
    information) on other programming platform are still complicated.
  - The readability and reusable will also be decreased when massive
    Python and R codes mixed with the workflows language codes.

This is an experimental project to providing a set of tools for the
exploring next-generation sequencing (NGS) data. We aim to integrate and
develop command line tools, NGS data analysis workflows
\[[WDL](https://github.com/openwdl/wdl),
[Nextflow](https://www.nextflow.io/),
[snakemake](https://snakemake.readthedocs.io/en/stable/), and
[bpipe](https://github.com/ssadedin/bpipe)\], and R shiny plugins/R
markdown document.

<p align="center">
  <img 
      alt="Best practice of reproducible NGS data analysis projects"
      src="https://raw.githubusercontent.com/Miachol/ftp/master/files/images/ngsjs/reproducible_NGS_data_analysis_projects_best_practice.jpg"
  />
</p>

We proposed that using [node](https://nodejs.org/en/) to distribute the
bioinformatics data analysis required workflows (e.g [Common workflow
language (CWL)](https://www.commonwl.org/)) and user created command
line scripts in data analysis process. The creation, update and upload
of a node package are very simple. Well-tested and high-performance
distribution tools of node packages, such as
[npm](https://www.npmjs.com/) and [yarn](https://www.yarnpkg.com), are
providing the service for more than 831,195 node packages.

**Command line scripts supported
now:**

| tool         | function                                                                                                                                                                         |
| ------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| rdeps        | Getting all `ngsjs` required R packages                                                                                                                                          |
| rsession     | `sessionInfo()` and `sessioninfo::session\_info()`                                                                                                                               |
| rinstall     | Install R packages and [BioInstaller](https://github.com/ngsjs/BioInstaller) resources using `install.packages()` and R packages `devtools`, `BiocManager` and `BioInstaller`    |
| rbashful     | Using the GO program [bashful](https://github.com/wagoodman/bashful), yaml and toml and R scripts to stitch together commands and bash snippits and run them with a bit of style |
| rconfig      | Parsing and generating json, ini, yaml, and toml format configuration files                                                                                                      |
| rclrs        | Generating colors for visulization using a theme key                                                                                                                             |
| rmv          | Formating the file names.                                                                                                                                                        |
| ranystr      | Generating any counts and any length random strings (e.g. Ies1y7fpgMVjsAyBAtTT)                                                                                                  |
| rtime\_stamp | Generating time stamp (e.g. 2018\_11\_15\_22\_43\_25\_, 2018/11/15/, 2018/11/15/22/).                                                                                            |
| rdownload    | Parallel download URLs with logs                                                                                                                                                 |
| rbin         | Collecting R packages inst/bin files to a directory, e.g. PATH                                                                                                                   |

We are collecting the CWL language created workflows and publish on the
[npm](https://www.npmjs.com/):

  - [ngsjs-wkfl-wdl](https://github.com/ngsjs/ngsjs-wkfl-wdl)
  - [ngsjs-wkfl-nextflow](https://github.com/ngsjs/ngsjs-wkfl-nextflow)
  - [ngsjs-wkfl-snakemake](https://github.com/ngsjs/ngsjs-wkfl-snakemake)
  - [ngsjs-wkfl-bpipe](https://github.com/ngsjs/ngsjs-wkfl-bpipe)

Besides, we are developing a framework to integrate various data
analysis workflows and command line scripts:

  - rbashful: A ngsjs command line tool to dynamically render env.toml
    and cli.yaml for a unified downstream analysis environment shared
    between all integrated tools, workflows, scripts.
  - cli.yaml: Process controller with the
    [bashful](https://github.com/wagoodman/bashful) style.
  - env.toml: Store the fields and values of input and output
    parameters; the core command line commands indexed by unique keys.
  - others

## Requirements

  - [node](https://nodejs.org/en/)
  - [R](https://cran.r-project.org/)
  - [GO](https://golang.org/)

### R packages

  - optparse
  - configr
  - stringi
  - futile.logger
  - glue
  - ngstk
  - BioInstaller
  - devtools
  - pacman
  - BiocManager
  - sessioninfo
  - future

## Installation

You need to install the [node](https://nodejs.org/en/),
[R](https://cran.r-project.org/) and [GO](https://golang.org/) for
running all [ngsjs](https://github.com/ngsjs/ngsjs) executable files.

``` bash
# Use conda to manage the env
conda install go nodejs \
echo 'export NODE_PATH="/path/miniconda2/lib/node_modules/"\n' >> ~/.bashrc \
&& npm install -g npm \
&& npm install -g yarn \
&& echo 'export PATH=$PATH:~/.yarn/bin/\n' >> ~/.bashrc

# Other see https://nodejs.org/en/download/package-manager/
# For: Ubuntu
apt update
apt install -y npm golang

# For MacOS
brew install node go
```

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

Before try your `ngsjs` command line tools, you need run the `rdeps`
getting all the extra R packages required by `ngsjs`.

``` bash
# install the extra R packages used in `ngsjs` scripts
rdeps
#> INFO [2018-11-17 22:50:04] All basic dependences (R packages) were resolved.
#> INFO [2018-11-17 22:50:04] devtools, BiocManager, sessioninfo, glue, futile.logger, stringi, future, configr, ngstk, BioInstaller, ngsjs
```

Then you can use the `ngsjs` to run all sub-commands.

``` bash
ngsjs -h
#> Description:
#>  Main interface of ngsjs package tools. View the ngsjs homepage https://github.com/ngsjs/ngsjs for more detail.
#>  Now, total 11 subcmds are supported: rbashful rconfig rdeps rinstall rsession rclrs rmv rtime_stamp ranystr rdownload rbin.
#> 
#> Usage: 
#>  /usr/local/bin/ngsjs [subcmds] [options]
#>  /usr/local/bin/ngsjs rconfig -h
#> 
#> Options:
#>  -l, --list-all-subcmds
#>                Print all supported subcmds of ngsjs.
#>  -h, --help
#>                Show this help message and exit
#> 
#> Commands:
#>  rbashful      An R-based tool to extend GO bashful tool for style bash commands
#>  rconfig       An R-based tool to parse and generate configuration file
#>  rdeps         An R-based tool to install ngsjs command line tools required R packages.
#>  rinstall      An R-based tool to install or download R packages and other resources supported by R package BioInstaller.
#>  rsession      An R-based tool to show R environment using sessionInfo() and sessioninfo::session_info()
#>  rclrs         An R-based tool to generate colors for visulization using a theme key
#>  rmv           An R-based tool to format file names
#>  rtime_stamp   An R-based tool to generate time stamp
#>  ranystr       An R-based tool to generate any counts and any urls random strings
#>  rdownload     An R-based tool to concurrently download urls with logging.
#>  rbin          An R-based tool to collect R packages bin files.
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
[here](https://github.com/ngsjs/ngsjs/test/rbashful/rnaseq_splicing).

``` r
source_dir <- "/Users/ljf/Documents/repositories/ljf/github/ngsjs/examples/rbashful/rnaseq_splicing/02_leafcutter_majiq"

# View the cli.yaml
cat(paste0(readLines(sprintf("%s/cli.yaml", source_dir)), 
           collapse = "\n"), sep = "\n")
#> config:
#>     log-path: "{{submit_log}}"
#> tasks:
#>     - name: "RNA-seq splicing analysis step 01: parse the task YAML and link bams"
#>       parallel-tasks: 
#>       - cmd: ./scripts/01_link_bams.R --env-parsed-yaml "{{env_parsed_yaml}}" --task-name "majiq_splicing_01_link_bams"

# View the env.toml
cat(paste0(readLines(sprintf("%s/env.toml", source_dir)), 
           collapse = "\n"), sep = "\n")
#> title = "Environment variables for running splicing analysis (from bam start)"
#> 
#> global_vars = ["project_dir", "workdir", "submit_yaml", "log_dir", 
#>                "env_parsed_yaml", "env_toml", "rbashful", "submit_log"]
#> 
#> # Global var (Replace all '{{var}}' in cli.yaml).
#> project_dir = "/home/ljf/projects/eqtf" 
#> workdir = "!!glue {config$project_dir}/analysis/rnaseq_splicing"
#> submit_yaml = "!!glue {config$workdir}/cli.yaml"
#> log_dir = "!!glue {config$workdir}/logs" 
#> env_toml = "!!glue {config$workdir}/env.toml"
#> env_parsed_yaml = "!!glue {config$workdir}/env.parsed.yaml"
#> rbashful = "rbashful"
#> 
#> # Sample specific var
#> id = "{{id}}"
#> [input]
#> samples_info_fn = "{{workdir}}/samples"
#> leafcutter_root = "/opt/bin/caller/leafcutter"
#> bam2_junc_script = "!!glue {config$input$leafcutter_root}/scripts/bam2junc.sh"
#> rawbamdir = "{{project_dir}}/analysis/rnaseq/output/bam/merge"
#> rawbamfile = "!!glue {config$input$rawbamdir}/{config$id}{config$input$bam_suffix}"
#> bamdir = "{{project_dir}}/analysis/rnaseq_splicing/bam_tmp"
#> bam_suffix = ".bam_AddGroup.bam_MarkDup.bam_SplitNtrim.bam_IndelRealigner.bam_PrintReads.bam"
#> bamfile = "!!glue {config$input$bamdir}/{config$id}{config$input$bam_suffix}"
#> genome = "hg19"
#> genome_path = "/u4/jhuangdata/reference/ucsc/hg19"
#> hg19_gtf = "!!glue {config$input$genome_path}/Homo_sapiens.GRCh37.75.gtf"
#> hg19_gff = "!!glue {config$input$genome_path}/Homo_sapiens.GRCh37.75.gff"
#> 
#> [output]
#> leafcutter_out_dir = "!!glue {config$project_dir}}/analysis/rnaseq_splicing/leafcutter/bamfile_juncs/"
#> majiq_build_out_dir = "!!glue {config$project_dir}/analysis/rnaseq_splicing/output/majiq_build"
#> leafcutter_bamfile_junc = "!!glue {config$output$leafcutter_out_dir}/{config$id}.junc"
#> 
#> [cmds]
#> bam2junc = "!!glue sh {config$input$bam2_junc_script} {config$input$bamfile} {config$output$leafcutter_bamfile_junc}"
#> generate_majiq_conf = """
#> !!glue rconfig majiq.ini -e {{dqm}}extra.list=list(bamdir='{config$input$bamdir}', \
#>          genome = '{config$input$genome}', genome_path='{config$input$genome_path}'){{dqm}} \
#>          -r {{dqm}}write.config(x, 'majiq.parsed.ini'){{dqm}}\
#> """
#> majiq_splicing_01_link_bams = """
#> !!glue rawbam={config$input$rawbamdir}/{config$id}{config$input$bam_suffix}; \
#> if [ -f $rawbam ]
#> then
#>   ln -s $rawbam {config$input$bamdir}
#> fi
#> """
#> 
#> majiq_builder_step = """
#> !!glue majiq build {config$input$hg19_gff} -conf majiq.parsed.ini --nthreads 10 --output {config$output$majiq_build_out_dir}
#> """

# View the submit.sh
cat(paste0(readLines(sprintf("%s/submit", source_dir)), 
           collapse = "\n"), sep = "\n")
#> #! /usr/bin/env Rscript
#> 
#> # yarn global add ngsjs
#> # system("rdeps")
#> pkgs <- c("stringi", "configr")
#> pacman::p_load(pkgs, character.only = TRUE)
#> 
#> env_toml = "/u7/home/ljf/projects/eqtf/analysis/rnaseq_splicing/env.toml"
#> 
#> config <- read.config(env_toml, extra.list = list('dqm'='"'), 
#>                       global.vars.field = NULL, rcmd.parse = TRUE)
#> rm(env_toml)
#> for (i in 1:1) config <- parse.extra(config, glue.parse = TRUE, 
#>                                      global.vars.field = NULL)
#> for (i in 1:1) config <- parse.extra(config, glue.parse = TRUE)
#> config$submit_log <- paste0(file.path(config[["workdir"]], "logs", format(Sys.time(), "%Y_%m_%d_%H_%M_%S_")),
#>                      stri_rand_strings(1, 20), ".log")
#> attach(config)
#> x <- write.config(config, env_parsed_yaml, write.type = "yaml", indent = 4)
#> if (!x) {stop(sprintf("Generating %s failed.", env_parsed_yaml))}
#> 
#> extra_list <- paste0(global_vars, "='", unname(sapply(config[global_vars],
#>                      function(x){return(x)[1]})), "'", collapse = ", ")
#> print(extra_list)
#> cmd <- sprintf('./rbashful -c %s -e "%s" -p >> %s', submit_yaml, extra_list, submit_log)
#> if (!dir.exists(dirname(submit_log))) dir.create(dirname(submit_log))
#> cat(cmd, file = submit_log, append = TRUE)
#> message(cmd, sep = "\n")
#> system(cmd)
#> 
#> cmd <- "bashful run cli.parsed.yaml"
#> message(cmd, sep = "\n")
#> system(cmd)
```

``` bash
rbashful -h
#> Description:
#>  rbashful is an extend bashful tool for style bash commands.
#> 
#> Usage: 
#>  /usr/local/bin/rbashful [options] [params]
#>  /usr/local/bin/rbashful -c ${workdir}/cli.yaml --env-toml ${workdir}/env.toml --cmd-name default -v
#> 
#> Options:
#>  -v, --verbose
#>                Print extra output [FALSE]
#>  -c CLI-YAML, --cli-yaml=CLI-YAML
#>                bashful used YAML file [cli.yaml]
#>  -t ENV-TOML, --env-toml=ENV-TOML
#>                TOML file stores environment variables [env.toml]
#>  -e EXTRA-LIST, --extra-list=EXTRA-LIST
#>                Need to replaced environment variables
#>  -p, --parse-cli-yaml
#>                Replace cli config keys [FALSE]
#>  -o OUTPUT-CLI-YAML, --output-cli-yaml=OUTPUT-CLI-YAML
#>                Output file of parsed cli YAML file [*.parsed.yaml]
#>  -n CMD-NAME, --cmd-name=CMD-NAME
#>                Run CMDs section using name [NULL]
#>  --auto-create-dir
#>                Auto create dir in env.toml output section [FALSE]
#>  -h, --help
#>                Show this help message and exit
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

``` bash
rsession -f 2 -e 'include_base=TRUE'

rsession -h
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
#>  date     2018-11-17                  
#> 
#> ─ Packages ───────────────────────────────────────────────────────────────────
#>  ! package     * version date       lib source                          
#>    assertthat    0.2.0   2017-04-11 [1] CRAN (R 3.5.0)                  
#>    base        * 3.5.1   2018-07-05 [?] local                           
#>    cli           1.0.1   2018-09-25 [1] CRAN (R 3.5.0)                  
#>    codetools     0.2-15  2016-10-05 [1] CRAN (R 3.5.1)                  
#>  P compiler      3.5.1   2018-07-05 [1] local                           
#>    configr       0.3.4.1 2018-11-14 [1] Github (Miachol/configr@0df7b68)
#>    crayon        1.3.4   2017-09-16 [1] CRAN (R 3.5.0)                  
#>    data.table    1.11.8  2018-09-30 [1] CRAN (R 3.5.0)                  
#>  P datasets    * 3.5.1   2018-07-05 [1] local                           
#>    digest        0.6.18  2018-10-10 [1] CRAN (R 3.5.0)                  
#>    future        1.10.0  2018-10-17 [1] CRAN (R 3.5.0)                  
#>    getopt        1.20.2  2018-02-16 [1] CRAN (R 3.5.0)                  
#>    globals       0.12.4  2018-10-11 [1] CRAN (R 3.5.0)                  
#>    glue          1.3.0   2018-07-17 [1] CRAN (R 3.5.0)                  
#>  P graphics    * 3.5.1   2018-07-05 [1] local                           
#>  P grDevices   * 3.5.1   2018-07-05 [1] local                           
#>    ini           0.3.1   2018-05-20 [1] CRAN (R 3.5.0)                  
#>    jsonlite      1.5     2017-06-01 [1] CRAN (R 3.5.0)                  
#>    listenv       0.7.0   2018-01-21 [1] CRAN (R 3.5.0)                  
#>    magrittr      1.5     2014-11-22 [1] CRAN (R 3.5.0)                  
#>  P methods     * 3.5.1   2018-07-05 [1] local                           
#>    ngstk       * 0.2.2.4 2018-11-17 [1] Github (JhuangLab/ngstk@2548f83)
#>    optparse      1.6.0   2018-06-17 [1] CRAN (R 3.5.0)                  
#>    pacman      * 0.5.0   2018-10-22 [1] CRAN (R 3.5.0)                  
#>  P parallel      3.5.1   2018-07-05 [1] local                           
#>    Rcpp          1.0.0   2018-11-07 [1] CRAN (R 3.5.0)                  
#>    RcppTOML      0.1.5   2018-10-31 [1] CRAN (R 3.5.0)                  
#>    rstudioapi    0.8     2018-10-02 [1] CRAN (R 3.5.0)                  
#>    sessioninfo   1.1.1   2018-11-05 [1] CRAN (R 3.5.1)                  
#>  P stats       * 3.5.1   2018-07-05 [1] local                           
#>    stringi       1.2.4   2018-07-20 [1] CRAN (R 3.5.0)                  
#>    stringr       1.3.1   2018-05-10 [1] CRAN (R 3.5.0)                  
#>  P tools         3.5.1   2018-07-05 [1] local                           
#>  P utils       * 3.5.1   2018-07-05 [1] local                           
#>    withr         2.1.2   2018-03-15 [1] CRAN (R 3.5.0)                  
#>    yaml          2.2.0   2018-07-25 [1] CRAN (R 3.5.0)                  
#> 
#> [1] /Library/Frameworks/R.framework/Versions/3.5/Resources/library
#> 
#>  P ── Loaded and on-disk path mismatch.
#> Description:
#>  rsession is an R-based tool to show R environment using sessionInfo() and sessioninfo::session_info().
#> 
#> Usage: 
#>  /usr/local/bin/rsession [options] [params]
#>  /usr/local/bin/rsession
#>  /usr/local/bin/rsession -f 1
#>  /usr/local/bin/rsession -f 2 -e 'include_base=TRUE'
#>  /usr/local/bin/rsession -d
#> 
#> Options:
#>  -v, --verbose
#>                Print extra output [default FALSE]
#>  -f FUNC, --func=FUNC
#>                Function name [e.g. sessionInfo (1), sessioninfo::session_info (2)]
#>  -e EXTRA, --extra=EXTRA
#>                Extra parameters [e.g. include_base=TRUE]
#>  -d, --doc
#>                Print functions document
#>  -h, --help
#>                Show this help message and exit
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

# Install R package ngstk from GitHub ngsjs/ngstk (devtools::install_github)
rinstall -f 2 ngsjs/ngstk

# Install R package ngstk from GitHub ngsjs/ngstk (install.package)
# devtools::install_github with force = TRUE, ref = 'develop'
rinstall -f 2 -e "force = TRUE, ref = 'develop'" ngsjs/ngstk

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

``` bash
rinstall -h
#> Description:
#>  rinstall is an R-based tool to install or download R packages and other resources supported by R package BioInstaller.
#> 
#> Usage: 
#>  /usr/local/bin/rinstall [options] [params]
#>  Examples:
#>  /usr/local/bin/rinstall -p ini
#>  /usr/local/bin/rinstall ini,yaml
#>  /usr/local/bin/rinstall -f 2 JhuangLab/ngstk
#>  /usr/local/bin/rinstall -f 2 -e "force = TRUE, ref = 'develop'" JhuangLab/ngstk
#>  /usr/local/bin/rinstall -f 3 ggtree; /usr/local/bin/rinstall rinstall -f 4 ggtree
#>  /usr/local/bin/rinstall -f 5 -e "show.all.names=T"
#>  /usr/local/bin/rinstall -f 5 -e "show.all.versions=T" db_annovar_avsnp
#>  /usr/local/bin/rinstall -f 5 -e "download.dir='/tmp/avsnp', extra.list=list(buildver='hg19')" db_annovar_avsnp
#> 
#> Options:
#>  -v, --verbose
#>                Print extra output [default FALSE]
#>  -f FUNC, --func=FUNC
#>                Index or name of used function [e.g. install.packages (1), devtools::install_github (2), BiocManager::install (3), pacman::p_load (4), BioInstaller::install.bioinfo (5)].
#>  -p PKGS, --pkgs=PKGS
#>                Package or item names [e.g. ggplot2,stringr or JhuangLab/BioInstaller (mode is devtools::install_github)].
#>  -e EXTRA, --extra=EXTRA
#>                 Extra parameters [e.g. ref='develop'].
#>  -d, --doc
#>                Print functions document
#>  -h, --help
#>                Show this help message and exit
```

### rdownload

``` bash
rdownload "https://img.shields.io/npm/dm/ngsjs.svg,https://img.shields.io/npm/v/ngsjs.svg,https://img.shields.io/npm/l/ngsjs.svg"

rdownload "https://img.shields.io/npm/dm/ngsjs.svg,https://img.shields.io/npm/v/ngsjs.svg,https://img.shields.io/npm/l/ngsjs.svg" --destfiles "/tmp/ngsjs1.svg,ngsjs2.svg,ngsjs3.svg"

rdownload --urls "https://img.shields.io/npm/dm/ngsjs.svg , https://img.shields.io/npm/v/ngsjs.svg,https://img.shields.io/npm/l/ngsjs.svg" \
          --destfiles "ngsjs1.svg,ngsjs2.svg,ngsjs3.svg" --max-cores 1
```

``` bash
rdownload -h
#> Description:
#>  rdownload is an R-based tool to concurrently download urls with logging.
#> 
#> Usage: 
#>  /usr/local/bin/rdownload [options] [params]
#>  /usr/local/bin/rdownload "https://img.shields.io/npm/dm/ngsjs.svg,https://img.shields.io/npm/l/ngsjs.svg"
#>  /usr/local/bin/rdownload "https://img.shields.io/npm/dm/ngsjs.svg,https://img.shields.io/npm/l/ngsjs.svg" \ 
#>          --destfiles "ngsjs1.svg,ngsjs2.svg"
#>  /usr/local/bin/rdownload --urls "https://img.shields.io/npm/dm/ngsjs.svg,https://img.shields.io/npm/l/ngsjs.svg" \ 
#>          --destfiles "ngsjs1.svg,ngsjs2.svg"
#>  /usr/local/bin/rdownload --urls "https://img.shields.io/npm/dm/ngsjs.svg,https://img.shields.io/npm/l/ngsjs.svg" \ 
#>          --destfiles "ngsjs1.svg,ngsjs2.svg" --max-cores 1
#> 
#> Options:
#>  -v, --verbose
#>                Print extra output [default FALSE]
#>  -f FUNC, --func=FUNC
#>                Index or name of used function [e.g. ngstk::par_download (1).
#>  -u URLS, --urls=URLS
#>                URLs of of a resource to be downloaded (multiple files split by ',' or ';'
#>  --destfiles=DESTFILES
#>                Filenames of downloaded files, default use the basename(urls)
#>  --max-cores=MAX-CORES
#>                Define the maxium used cores [future::availableCores()]
#>  -r RFUNC, --rfunc=RFUNC
#>                R function (input param 'x') to process the returned colors or urls [e.g. function(x){return(x[[2]])} or x[[2]]]
#>  -e EXTRA, --extra=EXTRA
#>                Extra parameters [...]
#>  -d, --doc
#>                Print functions document
#>  -h, --help
#>                Show this help message and exit
```

### rconfig

``` r
# Use configr::read.config parsing json format file
# Reture the list object output
rconfig package.json
rconfig -c package.json

# Use configr::read.config parsing json format file with the custom R function
rconfig -c test.json -r 'function(x){x[["a"]] + x[["b"]]}'
rconfig -c test.json -r 'function(x){x[["a"]]}'
rconfig -c test.json -r 'function(x){x[["b"]]}'
rconfig -c test.json -r 'x[["b"]]'

# Use configr::write.config parsing json format file
rconfig -f "configr::write.config" test.json -e "config.dat=list(a=1, b=2), write.type='json'"
rconfig -f 2 test.json -e "config.dat=list(a=1, b=2), write.type='json'"
```

``` bash
rconfig -f "configr::fetch.config" "https://raw.githubusercontent.com/Miachol/configr/master/inst/extdata/config.global.toml"

rconfig -h
#> trying URL 'https://raw.githubusercontent.com/Miachol/configr/master/inst/extdata/config.global.toml'
#> Content type 'text/plain; charset=utf-8' length 303 bytes
#> ==================================================
#> downloaded 303 bytes
#> 
#> List of 7
#>  $ global_vars: chr [1:4] "gvar_1" "gvar_2" "gvar_3" "gvar_5"
#>  $ gvar_1     : chr "G1"
#>  $ gvar_2     : chr "G2"
#>  $ gvar_3     : chr "G3"
#>  $ gvar_5     : chr "G5"
#>  $ subsection :List of 4
#>   ..$ value_1: chr "G1/value_1"
#>   ..$ value_2: chr "G2/value_2"
#>   ..$ value_3: chr "G3/value_3"
#>   ..$ value_5: chr "G5/value_5"
#>  $ title      : chr "Demo of global vars of configuration files"
#> Description:
#>  rconfig is an R-based tool to parse and generate configuration file.
#> 
#> Usage: 
#>  /usr/local/bin/rconfig [options] [params]
#>  /usr/local/bin/rconfig package.json
#>  /usr/local/bin/rconfig -c package.json
#>  /usr/local/bin/rconfig -f 2 test.json -e "config.dat=list(a=1, b=2), write.type='json'"
#>  /usr/local/bin/rconfig -f "configr::write.config" test.json -e "config.dat=list(a=1, b=2), write.type='json'"
#>  /usr/local/bin/rconfig -i test.json -r 'function(x){x[["a"]] + x[["b"]]}'
#>  /usr/local/bin/rconfig -i test.json -r 'function(x){x[["a"]]}'
#>  /usr/local/bin/rconfig -i test.json -r 'function(x){x[["b"]]}'
#> 
#> Options:
#>  -v, --verbose
#>                Print extra output [default FALSE]
#>  -f FUNC, --func=FUNC
#>                Index or name of used function [e.g. configr::read.config (1), configr::fetch.config (2), configr::write.config (3)].
#>  -c CFG, --cfg=CFG
#>                Input or output configuationo file.
#>  -e EXTRA, --extra=EXTRA
#>                Extra parameters [e.g. extra.list=list(key='value')].
#>  -r RFUNC, --rfunc=RFUNC
#>                R function (input param 'x') to process the parsed configuation file [e.g. function(x){return(x[[1]])} or x[[1]] ].
#>  -d, --doc
#>                Print functions document
#>  -h, --help
#>                Show this help message and exit
```

### rclrs

``` bash
# Show default and red/blue theme colors
rclrs default
rclrs -t default
rclrs -t red_blue

# Show default theme colors (extract the first element)
rclrs -t default -r 'x[1]'

# Show all supported theme
rclrs --show-all-themes
#> #0073c3
#> #efc000
#> #696969
#> #ce534c
#> #7ba6db
#> #035892
#> #052135
#> #666633
#> #660000
#> #990000
#> #0073c3
#> #efc000
#> #696969
#> #ce534c
#> #7ba6db
#> #035892
#> #052135
#> #666633
#> #660000
#> #990000
#> #c20b01
#> #196abd
#> #0073c3
#> List of 35
#>  $ Title                       : chr "ngstk theme configuration file (colors)"
#>  $ adobe_color_cc_1            :List of 1
#>   ..$ colors: chr [1:10] "#FFE350" "#E8740C" "#FF0000" "#9C0CE8" ...
#>  $ ball_subtype_colors         :List of 4
#>   ..$ colors : chr [1:10] "#0071bc" "#82a6d7" "#003b64" "#e8c122" ...
#>   ..$ doi    : chr "10.1200/JCO.2016.70.7836"
#>   ..$ figure : chr "figure1"
#>   ..$ journal: chr "Journal of Clinical Oncology"
#>  $ default                     :List of 1
#>   ..$ colors: chr [1:10] "#0073c3" "#efc000" "#696969" "#ce534c" ...
#>  $ ggsci_aaas_default          :List of 1
#>   ..$ colors: chr [1:10] "#3B4992" "#EE0000" "#008B45" "#631879" ...
#>  $ ggsci_d3_category10         :List of 1
#>   ..$ colors: chr [1:10] "#1F77B4" "#FF7F0E" "#2CA02C" "#D62728" ...
#>  $ ggsci_d3_category20         :List of 1
#>   ..$ colors: chr [1:20] "#1F77B4" "#FF7F0E" "#2CA02C" "#D62728" ...
#>  $ ggsci_d3_category20b        :List of 1
#>   ..$ colors: chr [1:20] "#393B79" "#637939" "#8C6D31" "#843C39" ...
#>  $ ggsci_d3_category20c        :List of 1
#>   ..$ colors: chr [1:20] "#3182BD" "#E6550D" "#31A354" "#756BB1" ...
#>  $ ggsci_futurama_planetexpress:List of 1
#>   ..$ colors: chr [1:12] "#FF6F00" "#C71000" "#008EA0" "#8A4198" ...
#>  $ ggsci_gsea_default          :List of 1
#>   ..$ colors: chr [1:12] "#4500AD" "#2700D1" "#6B58EF" "#8888FF" ...
#>  $ ggsci_igv_alternating       :List of 1
#>   ..$ colors: chr [1:2] "#5773CC" "#FFB900"
#>  $ ggsci_igv_default           :List of 1
#>   ..$ colors: chr [1:51] "#5050FF" "#CE3D32" "#749B58" "#F0E685" ...
#>  $ ggsci_jama_defalut          :List of 1
#>   ..$ colors: chr [1:7] "#374E55" "#DF8F44" "#00A1D5" "#B24745" ...
#>  $ ggsci_jco_default           :List of 1
#>   ..$ colors: chr [1:10] "#0073C2" "#EFC000" "#868686" "#CD534C" ...
#>  $ ggsci_lancet_lanonc         :List of 1
#>   ..$ colors: chr [1:9] "#00468B" "#ED0000" "#42B540" "#0099B4" ...
#>  $ ggsci_locuszoom             :List of 1
#>   ..$ colors: chr [1:7] "#D43F3A" "#EEA236" "#5CB85C" "#46B8DA" ...
#>  $ ggsci_nejm_default          :List of 1
#>   ..$ colors: chr [1:8] "#BC3C29" "#0072B5" "#E18727" "#20854E" ...
#>  $ ggsci_npg_nrc               :List of 1
#>   ..$ colors: chr [1:10] "#E64B35" "#4DBBD5" "#00A087" "#3C5488" ...
#>  $ ggsci_rickandmorty_schwifty :List of 1
#>   ..$ colors: chr [1:12] "#FAFD7C" "#82491E" "#24325F" "#B7E4F9" ...
#>  $ ggsci_simpsons_springfield  :List of 1
#>   ..$ colors: chr [1:16] "#FED439" "#709AE1" "#8A9197" "#D2AF81" ...
#>  $ ggsci_startrek_uniform      :List of 1
#>   ..$ colors: chr [1:7] "#CC0C00" "#5C88DA" "#84BD00" "#FFCD00" ...
#>  $ ggsci_tron_legacy           :List of 1
#>   ..$ colors: chr [1:7] "#FF410D" "#6EE2FF" "#F7C530" "#95CC5E" ...
#>  $ ggsci_uchicago_dark         :List of 1
#>   ..$ colors: chr [1:9] "#800000" "#767676" "#CC8214" "#616530" ...
#>  $ ggsci_uchicago_default      :List of 1
#>   ..$ colors: chr [1:9] "#800000" "#767676" "#FFA319" "#8A9045" ...
#>  $ ggsci_uchicago_light        :List of 1
#>   ..$ colors: chr [1:9] "#800000" "#D6D6CE" "#FFB547" "#ADB17D" ...
#>  $ ggsci_ucscgb_default        :List of 1
#>   ..$ colors: chr [1:26] "#FF0000" "#FF9900" "#FFCC00" "#00FF00" ...
#>  $ nature_brest_signatures     :List of 4
#>   ..$ colors : chr [1:12] "#3d1572" "#7d4594" "#e84286" "#f7c0ba" ...
#>   ..$ doi    : chr "10.1038/nature17676"
#>   ..$ figure : chr "figure3b"
#>   ..$ hournal: chr "Nature"
#>  $ ng_mutations                :List of 4
#>   ..$ colors : chr [1:9] "#609ec2" "#b56248" "#d0cb6c" "#9cb46f" ...
#>   ..$ doi    : chr "10.1038/ng.3900"
#>   ..$ figure : chr "figure4a"
#>   ..$ journal: chr "Nature Gentics"
#>  $ nm_lines                    :List of 4
#>   ..$ colors : chr [1:6] "#3b7ab5" "#e7211c" "#ff831d" "#2ee0d1" ...
#>   ..$ doi    : chr "10.1038/nmeth.4083"
#>   ..$ figure : chr "figure2"
#>   ..$ journal: chr "Nature Methods"
#>  $ proteinpaint_chromHMM_state :List of 2
#>   ..$ colors : chr [1:16] "#c0222c" "#f12424" "#ff00c7" "#d192fb" ...
#>   ..$ journal: chr "Nature Gentics"
#>  $ proteinpaint_domains        :List of 2
#>   ..$ colors : chr [1:8] "#a6d854" "#8dd3c7" "#fb8072" "#80b1d3" ...
#>   ..$ journal: chr "Nature Gentics"
#>  $ proteinpaint_mutations      :List of 2
#>   ..$ colors : chr [1:10] "#3987cc" "#ff7f0e" "#db3d3d" "#6633ff" ...
#>   ..$ journal: chr "Nature Gentics"
#>  $ proteinpaint_significance   :List of 2
#>   ..$ colors : chr [1:7] "#aaaaaa" "#e99002" "#5bc0de" "#f04124" ...
#>   ..$ journal: chr "Nature Gentics"
#>  $ red_blue                    :List of 1
#>   ..$ colors: chr [1:2] "#c20b01" "#196abd"
```

``` bash
rclrs -h
#> Description:
#>  rclrs is an R-based tool to generate colors for visulization using a theme key.
#> 
#> Usage: 
#>  /usr/local/bin/rclrs [options] [params]
#>  /usr/local/bin/rclrs default
#>  /usr/local/bin/rclrs -t default
#>  /usr/local/bin/rclrs -t default -r 'x[1]'
#>  /usr/local/bin/rclrs -t red_blue
#>  /usr/local/bin/rclrs --show-all-themes
#> 
#> Options:
#>  -v, --verbose
#>                Print extra output [default FALSE]
#>  -f FUNC, --func=FUNC
#>                Index or name of used function [e.g. ngstk::set_colors (1).
#>  -t THEME, --theme=THEME
#>                Input the theme name and return the colors.
#>  -s, --show-all-themes
#>                Show all included themes.
#>  -r RFUNC, --rfunc=RFUNC
#>                R function (input param 'x') to process the returned colors or theme [e.g. function(x){return(x[[2]])} or x[[2]]].
#>  -e EXTRA, --extra=EXTRA
#>                Extra parameters [e.g. theme_config_file = 'your_color_cfg.toml'].
#>  -d, --doc
#>                Print functions document
#>  -h, --help
#>                Show this help message and exit
```

### rmv

``` bash
# do.rename is used to preview the new filenames
# 
rmv "`ls`" -e "do.rename = F, prefix = 'prefix', suffix = 'suffix'"
rmv "`ls`" -e "do.rename = F, replace = list(old =c('-', '__'), new = c('_', '_'))"
rmv "`ls`" -e "do.rename = F, toupper = TRUE"
rmv "`ls`" -e "do.rename = F, tolower = TRUE"

rmv "`ls`" -e "do.rename=T, replace=list(old='new', new='old')"
```

``` bash
rmv -h
#> Description:
#>  rmv is an R-based tool to format file names.
#> 
#> Usage: 
#>  /usr/local/bin/rmv [options] [params]
#>  /usr/local/bin/rmv "`ls`" -e "do.rename = FALSE, prefix = 'prefix', suffix = 'suffix'"
#>  /usr/local/bin/rmv "`ls`" -e "do.rename = FALSE, replace = list(old =c('-', '__'), new = c('_', '_'))"
#>  /usr/local/bin/rmv "`ls`" -e "do.rename = FALSE, toupper = TRUE"
#>  /usr/local/bin/rmv "`ls`" -e "do.rename = FALSE, tolower = TRUE"
#>  /usr/local/bin/rmv -e "files_dir = '.', pattern = '.*.txt', do.rename=F, replace=list(old='old', new='new')"
#>  /usr/local/bin/rmv "`ls`" -e "do.rename=T, replace=list(old='old', new='new')"
#> 
#> Options:
#>  -v, --verbose
#>                Print extra output [default FALSE]
#>  -f FUNC, --func=FUNC
#>                Index or name of used function [e.g. ngstk::format_filenames (1).
#>  -l OLD-FILES, --old-files=OLD-FILES
#>                Input the files need to be renamed (string will be split by '\n', ',' and ';').
#>  -r RFUNC, --rfunc=RFUNC
#>                R function (input param 'x') to process the returned colors or old-files [e.g. function(x){return(x[[2]])} or x[[2]]].
#>  -e EXTRA, --extra=EXTRA
#>                Extra parameters [e.g. toupper = TRUE)].
#>  -d, --doc
#>                Print functions document
#>  -h, --help
#>                Show this help message and exit
```

### rtime\_stamp

``` bash
rtime_stamp

rtime_stamp -r 'x[[1]]'

rtime_stamp -r 'x[[1]][1]'

rtime_stamp -t '%Y_%d'

rtime_stamp -e "extra_flag=c('*')"

rtime_stamp -h
#> [[1]]
#> [1] "2018_11_17_22_50_12_" "2018_11_17_22_50_"    "2018_11_17_22_"      
#> [4] "2018_11_17_"          "2018_11_"             "2018_"               
#> 
#> [[2]]
#> [1] "2018-11-17-22-50-12-" "2018-11-17-22-50-"    "2018-11-17-22-"      
#> [4] "2018-11-17-"          "2018-11-"             "2018-"               
#> 
#> [[3]]
#> [1] "2018/11/17/22/50/12/" "2018/11/17/22/50/"    "2018/11/17/22/"      
#> [4] "2018/11/17/"          "2018/11/"             "2018/"               
#> 
#> [[4]]
#> [1] "2018@11@17@22@50@12@" "2018@11@17@22@50@"    "2018@11@17@22@"      
#> [4] "2018@11@17@"          "2018@11@"             "2018@"               
#> 
#> 2018_11_17_22_50_13_
#> 2018_11_17_22_50_
#> 2018_11_17_22_
#> 2018_11_17_
#> 2018_11_
#> 2018_
#> 2018_11_17_22_50_13_
#> [[1]]
#> [1] "2018_17"
#> 
#> [[2]]
#> [1] "2018-17"
#> 
#> [[3]]
#> [1] "2018/17"
#> 
#> [[4]]
#> [1] "2018@17"
#> 
#> [[1]]
#> [1] "2018_11_17_22_50_14_" "2018_11_17_22_50_"    "2018_11_17_22_"      
#> [4] "2018_11_17_"          "2018_11_"             "2018_"               
#> 
#> [[2]]
#> [1] "2018*11*17*22*50*14*" "2018*11*17*22*50*"    "2018*11*17*22*"      
#> [4] "2018*11*17*"          "2018*11*"             "2018*"               
#> 
#> Description:
#>  rtime_stamp is an R-based tool to generate time stamp.
#> 
#> Usage: 
#>  /usr/local/bin/rtime_stamp [options] [params]
#>  /usr/local/bin/rtime_stamp
#>  /usr/local/bin/rtime_stamp -r 'x[[1]]'
#>  /usr/local/bin/rtime_stamp -r 'x[[1]][1]'
#>  /usr/local/bin/rtime_stamp -t '%Y_%d'
#>  /usr/local/bin/rtime_stamp -e "extra_flat=c('-')"
#> 
#> Options:
#>  -v, --verbose
#>                Print extra output [default FALSE]
#>  -f FUNC, --func=FUNC
#>                Index or name of used function [e.g. ngstk::time_stamp (1).
#>  -t TEMPLATE, --template=TEMPLATE
#>                Input the template and return the time stamps.
#>  -r RFUNC, --rfunc=RFUNC
#>                R function (input param 'x') to process the returned colors or template [e.g. function(x){return(x[[2]])} or x[[2]]].
#>  -e EXTRA, --extra=EXTRA
#>                Extra parameters [e.g. extra_flag = c('-', '/')].
#>  -d, --doc
#>                Print functions document
#>  -h, --help
#>                Show this help message and exit
```

### ranystr

``` bash
./bin/ranystr

./bin/ranystr -l 30

./bin/ranystr -l 20 -n 3
#> UnHysWEDkaqDsWE0xW7Z
#> Vw61Zwck8WzsquXlmPX134Qmec7ehp
#> FftocIUgArDT5o8K77Fq
#> kNclR3FdmxUGxALDqd1E
#> HsOXKVMZjDcpGpLnV066
```

``` bash
ranystr -h
#> Description:
#>  ranystr is an R-based tool to generate any counts and any length random strings.
#> 
#> Usage: 
#>  /usr/local/bin/ranystr [options] [params]
#>  /usr/local/bin/ranystr
#>  /usr/local/bin/ranystr -l 30
#>  /usr/local/bin/ranystr -l 20 -n 3
#> 
#> Options:
#>  -v, --verbose
#>                Print extra output [default FALSE]
#>  -f FUNC, --func=FUNC
#>                Index or name of used function [e.g. stringi::stri_rand_strings (1).
#>  -n NUM, --num=NUM
#>                Counts of random strings [1].
#>  -l LENGTH, --length=LENGTH
#>                Length of one of random string [20].
#>  -p PATTERN, --pattern=PATTERN
#>                Character vector specifying character classes to draw elements ['[A-Za-z0-9]'].
#>  -r RFUNC, --rfunc=RFUNC
#>                R function (input param 'x') to process the returned colors or length [e.g. function(x){return(x[[2]])} or x[[2]]].
#>  -e EXTRA, --extra=EXTRA
#>                Extra parameters [...].
#>  -d, --doc
#>                Print functions document
#>  -h, --help
#>                Show this help message and exit
```

``` bash
# Collect system.files("extdata", "bin", package = "ngstk")
# multiple packages (i.e. ngstk,configr) 
rbin ngstk

rbin --destdir /tmp/path ngstk

rbin -h
#> Copying ngstk bin/ demo_bin.sh to /Users/ljf/.ngstk/bin
#> Please set /Users/ljf/.ngstk/bin in your PATH to use the bin files.
#> Linux/Mac OS X: echo 'export PATH=$PATH:/Users/ljf/.ngstk/bin\n' >> ~/.bashrc
#> R users: echo 'Sys.setenv(PATH="/usr/local/Cellar/hugo/0.48/bin:/Users/ljf/Bioinfo/miniconda3/bin:/usr/local/Cellar/gnu-sed/4.4/bin:/Library/Frameworks/R.framework/Versions/3.5/Resources/library/ngstk/extdata/tools/rbash:/Users/ljf/Bioinfo/spack/bin:/Users/ljf/Bioinfo/miniconda3/bin:/usr/local/Cellar/gnu-sed/4.4/bin:/Library/Frameworks/R.framework/Versions/3.5/Resources/library/ngstk/extdata/tools/rbash:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/Library/TeX/texbin:/usr/local/go/bin:/opt/X11/bin:/Users/ljf/.ngstk/bin")\n' >> ~/.Rprofile
#> ngstk.demo_bin.sh 
#>              TRUE 
#> Copying ngstk bin/ demo_bin.sh to /private/tmp/path
#> Please set /private/tmp/path in your PATH to use the bin files.
#> Linux/Mac OS X: echo 'export PATH=$PATH:/private/tmp/path\n' >> ~/.bashrc
#> R users: echo 'Sys.setenv(PATH="/usr/local/Cellar/hugo/0.48/bin:/Users/ljf/Bioinfo/miniconda3/bin:/usr/local/Cellar/gnu-sed/4.4/bin:/Library/Frameworks/R.framework/Versions/3.5/Resources/library/ngstk/extdata/tools/rbash:/Users/ljf/Bioinfo/spack/bin:/Users/ljf/Bioinfo/miniconda3/bin:/usr/local/Cellar/gnu-sed/4.4/bin:/Library/Frameworks/R.framework/Versions/3.5/Resources/library/ngstk/extdata/tools/rbash:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/Library/TeX/texbin:/usr/local/go/bin:/opt/X11/bin:/private/tmp/path")\n' >> ~/.Rprofile
#> ngstk.demo_bin.sh 
#>              TRUE 
#> Description:
#>  rbin is an R-based tool to collect R packages bin files.
#> 
#> Usage: 
#>  /usr/local/bin/rbin [options] [params]
#> 
#> Options:
#>  -v, --verbose
#>                Print extra output [default FALSE]
#>  -f FUNC, --func=FUNC
#>                Index or name of used function [e.g. ngstk::bin (1)].
#>  -c PKGS, --pkgs=PKGS
#>                Package names used to copy inst/bin directory files to PATH.
#>  --destdir=DESTDIR
#>                Destination directory to store the inst/bin files [~/.ngstk/bin].
#>  -e EXTRA, --extra=EXTRA
#>                Extra parameters [e.g. extra.list=list(key='value')].
#>  -r RFUNC, --rfunc=RFUNC
#>                R function (input param 'x') to process the parsed configuation file [e.g. function(x){return(x[[1]])} or x[[1]] ].
#>  -d, --doc
#>                Print functions document
#>  -h, --help
#>                Show this help message and exit
```

## How to contribute?

Please fork the [GitHub ngsjs
repository](https://github.com/ngsjs/ngsjs), modify it, and submit a
pull request to us.

## Maintainer

[Jianfeng Li](https://github.com/Miachol)

## License

MIT
