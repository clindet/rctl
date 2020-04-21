
<!-- README.md is generated from README.Rmd. Please edit that file -->

<p align="center">

<a href="https://www.npmjs.com/package/rctl"><img src="https://img.shields.io/badge/lifecycle-experimental-orange.svg" alt="Life cycle: experimental">
<a href="https://npmcharts.com/compare/rctl?minimal=true"><img src="https://img.shields.io/npm/dm/rctl.svg" alt="Downloads"></a>
<a href="https://www.npmjs.com/package/rctl"><img src="https://img.shields.io/npm/v/rctl.svg" alt="Version"></a>
<a href="https://www.npmjs.com/package/rctl"><img src="https://img.shields.io/npm/l/rctl.svg" alt="License"></a>

</p>

# rctl

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

This project is part of [openanno](https://github.com/openanno), and aim
to integrate and develop command line tools based on R and JavaScript
ecosystem.

<p align="center">

<img 
      alt="Best practice of reproducible NGS data analysis projects"
      src="https://raw.githubusercontent.com/Miachol/ftp/master/files/images/rctl/reproducible_NGS_data_analysis_projects_best_practice.jpg"
  />

</p>

We proposed that using [node](https://nodejs.org/en/) to distribute the
bioinformatics data analysis required workflows (e.g [Common workflow
language (CWL)](https://www.commonwl.org/)) or the command line scripts
that created by users. The creation, update and upload of a node package
are very simple. Well-tested and high-performance distribution tools of
node packages, such as [npm](https://www.npmjs.com/) and
[yarn](https://www.yarnpkg.com), are providing the service for more than
1,264,413 packages.

**Command line scripts supported now:**

| tool         | function                                                                                                                                                                     |
| ------------ | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| rdeps        | Getting all `rctl` required R packages                                                                                                                                       |
| rsession     | `sessionInfo()` and `sessioninfo::session\_info()`                                                                                                                           |
| rinstall     | Install R packages and [BioInstaller](https://github.com/rctl/BioInstaller) resources using `install.packages()` and R packages `devtools`, `BiocManager` and `BioInstaller` |
| rconfig      | Parsing and generating json, ini, yaml, and toml format configuration files                                                                                                  |
| rbashful     | run bashful                                                                                                                                                                  |
| rclrs        | Generating colors for visulization using a theme key                                                                                                                         |
| rmv          | Formating the file names.                                                                                                                                                    |
| ranystr      | Generating any counts and any length random strings (e.g. Ies1y7fpgMVjsAyBAtTT)                                                                                              |
| rtime\_stamp | Generating time stamp (e.g. 2018\_11\_15\_22\_43\_25\_, 2018/11/15/, 2018/11/15/22/).                                                                                        |
| rdownload    | Parallel download URLs with logs                                                                                                                                             |
| rbin         | Collecting R packages inst/bin files to a directory, e.g. PATH                                                                                                               |

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

``` bash
# Use conda to manage the env
conda install go nodejs \
echo 'export NODE_PATH="/path/miniconda2/lib/node_modules/"\n' >> ~/.bashrc \
&& npm install -g yarn \
&& echo 'export PATH=$PATH:~/.yarn/bin/\n' >> ~/.bashrc

# Other see https://nodejs.org/en/download/package-manager/
# For: Ubuntu
apt update
apt install -y npm

# For MacOS
brew install node
```

``` bash
npm install -g rctl
# or
yarn global add rctl

# If you not to globaly install rctl 
# You need to set the PATH
echo "export rctl_ROOT=/path/node_modules/nodejs" >> ~/.bashrc
echo "export PATH=$PATH:${rctl_ROOT}/bin" >> ~/.bashrc

# Current dir is /path
npm install rctl
# or
yarn add rctl
```

## Usage

Before try your `rctl` command line tools, you need run the `rdeps`
getting all the extra R packages required by `rctl`.

``` bash
# install the extra R packages used in `rctl` scripts
rdeps
```

Then you can use the `rctl` to run all sub-commands.

``` bash
rctl -h
#> Description:
#>  Main interface of rctl package tools. View the rctl homepage https://github.com/openanno/rctl for more detail.
#>  Now, total 11 subcmds are supported: rbashful rconfig rdeps rinstall rsession rclrs rmv rtime_stamp ranystr rdownload rbin.
#> 
#> Usage: 
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rctl [subcmds] [options]
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rctl rconfig -h
#> 
#> Options:
#>  -l, --list-all-subcmds
#>                Print all supported subcmds of rctl.
#>  -h, --help
#>                Show this help message and exit
#> 
#> Commands:
#>  rbashful      An R-based tool to extend GO bashful tool for style bash commands
#>  rconfig       An R-based tool to parse and generate configuration file
#>  rdeps         An R-based tool to install rctl command line tools required R packages.
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
by `rbashful`, so you need to install it before use the `rbashful`.

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

**Go tools**

``` bash
go get github.com/wagoodman/bashful
```

![](https://raw.githubusercontent.com/wagoodman/bashful/master/doc/demo.gif)

View a `rbashful` demo
[here](https://github.com/openanno/rctl/test/rbashful/rnaseq_splicing).

``` r
source_dir <- "/Users/apple/Documents/repositories/rctl/examples/rbashful/rnaseq_splicing/02_leafcutter_majiq"

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
#> # yarn global add rctl
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
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rbashful [options] [params]
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rbashful -c ${workdir}/cli.yaml --env-toml ${workdir}/env.toml --cmd-name default -v
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
#>  version  R version 3.6.3 (2020-02-29)
#>  os       macOS Catalina 10.15.4      
#>  system   x86_64, darwin15.6.0        
#>  ui       X11                         
#>  language (EN)                        
#>  collate  en_US.UTF-8                 
#>  ctype    en_US.UTF-8                 
#>  tz       Asia/Shanghai               
#>  date     2020-04-21                  
#> 
#> ─ Packages ───────────────────────────────────────────────────────────────────
#>  ! package     * version date       lib source                          
#>    assertthat    0.2.1   2019-03-21 [1] CRAN (R 3.6.0)                  
#>  R base        * 3.6.3   <NA>       [?] <NA>                            
#>    cli           2.0.2   2020-02-28 [1] CRAN (R 3.6.3)                  
#>  R codetools     0.2-16  <NA>       [2] <NA>                            
#>  R compiler      3.6.3   <NA>       [2] <NA>                            
#>    configr       0.3.4.1 2019-11-12 [1] Github (Miachol/configr@0df7b68)
#>    crayon        1.3.4   2017-09-16 [1] CRAN (R 3.6.0)                  
#>    data.table    1.12.8  2019-12-09 [1] CRAN (R 3.6.0)                  
#>  R datasets    * 3.6.3   <NA>       [2] <NA>                            
#>    digest        0.6.25  2020-02-23 [1] CRAN (R 3.6.3)                  
#>    fansi         0.4.1   2020-01-08 [1] CRAN (R 3.6.0)                  
#>    future        1.16.0  2020-01-16 [1] CRAN (R 3.6.3)                  
#>    getopt        1.20.3  2019-03-22 [1] CRAN (R 3.6.0)                  
#>    globals       0.12.5  2019-12-07 [1] CRAN (R 3.6.0)                  
#>    glue          1.3.2   2020-03-12 [1] CRAN (R 3.6.3)                  
#>  R graphics    * 3.6.3   <NA>       [2] <NA>                            
#>  R grDevices   * 3.6.3   <NA>       [2] <NA>                            
#>    ini           0.3.1   2018-05-20 [1] CRAN (R 3.6.0)                  
#>    jsonlite      1.6.1   2020-02-02 [1] CRAN (R 3.6.3)                  
#>    listenv       0.8.0   2019-12-05 [1] CRAN (R 3.6.0)                  
#>    magrittr      1.5     2014-11-22 [1] CRAN (R 3.6.0)                  
#>  R methods     * 3.6.3   <NA>       [2] <NA>                            
#>    ngstk       * 0.2.3   2019-11-12 [1] Github (JhuangLab/ngstk@e5953e3)
#>    optparse      1.6.4   2019-09-16 [1] CRAN (R 3.6.0)                  
#>    pacman      * 0.5.1   2019-03-11 [1] CRAN (R 3.6.0)                  
#>  R parallel      3.6.3   <NA>       [2] <NA>                            
#>    Rcpp          1.0.4   2020-03-17 [1] CRAN (R 3.6.3)                  
#>    RcppTOML      0.1.6   2019-06-25 [1] CRAN (R 3.6.0)                  
#>    rstudioapi    0.11    2020-02-07 [1] CRAN (R 3.6.3)                  
#>    sessioninfo   1.1.1   2018-11-05 [1] CRAN (R 3.6.0)                  
#>  R stats       * 3.6.3   <NA>       [2] <NA>                            
#>    stringi       1.4.6   2020-02-17 [1] CRAN (R 3.6.3)                  
#>    stringr       1.4.0   2019-02-10 [1] CRAN (R 3.6.0)                  
#>  R tools         3.6.3   <NA>       [2] <NA>                            
#>  R utils       * 3.6.3   <NA>       [2] <NA>                            
#>    withr         2.1.2   2018-03-15 [1] CRAN (R 3.6.0)                  
#>    yaml          2.2.1   2020-02-01 [1] CRAN (R 3.6.3)                  
#> 
#> [1] /Users/apple/.R/library/3.6
#> [2] /Library/Frameworks/R.framework/Versions/3.6/Resources/library
#> 
#>  R ── Package was removed from disk.
#> Description:
#>  rsession is an R-based tool to show R environment using sessionInfo() and sessioninfo::session_info().
#> 
#> Usage: 
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rsession [options] [params]
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rsession
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rsession -f 1
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rsession -f 2 -e 'include_base=TRUE'
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rsession -d
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

# Install R package ngstk from GitHub rctl/ngstk (devtools::install_github)
rinstall -f 2 rctl/ngstk

# Install R package ngstk from GitHub rctl/ngstk (install.package)
# devtools::install_github with force = TRUE, ref = 'develop'
rinstall -f 2 -e "force = TRUE, ref = 'develop'" rctl/ngstk

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
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rinstall [options] [params]
#>  Examples:
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rinstall -p ini
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rinstall ini,yaml
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rinstall -f 2 JhuangLab/ngstk
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rinstall -f 2 -e "force = TRUE, ref = 'develop'" JhuangLab/ngstk
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rinstall -f 3 ggtree; /Users/apple/.config/yarn/global/node_modules/.bin/rinstall rinstall -f 4 ggtree
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rinstall -f 5 -e "show.all.names=T"
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rinstall -f 5 -e "show.all.versions=T" db_annovar_avsnp
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rinstall -f 5 -e "download.dir='/tmp/avsnp', extra.list=list(buildver='hg19')" db_annovar_avsnp
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
rdownload "https://img.shields.io/npm/dm/rctl.svg,https://img.shields.io/npm/v/rctl.svg,https://img.shields.io/npm/l/rctl.svg"

rdownload "https://img.shields.io/npm/dm/rctl.svg,https://img.shields.io/npm/v/rctl.svg,https://img.shields.io/npm/l/rctl.svg" --destfiles "/tmp/rctl1.svg,rctl2.svg,rctl3.svg"

rdownload --urls "https://img.shields.io/npm/dm/rctl.svg , https://img.shields.io/npm/v/rctl.svg,https://img.shields.io/npm/l/rctl.svg" \
          --destfiles "rctl1.svg,rctl2.svg,rctl3.svg" --max-cores 1
```

``` bash
rdownload -h
#> Description:
#>  rdownload is an R-based tool to concurrently download urls with logging.
#> 
#> Usage: 
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rdownload [options] [params]
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rdownload "https://img.shields.io/npm/dm/rctl.svg,https://img.shields.io/npm/l/rctl.svg"
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rdownload "https://img.shields.io/npm/dm/rctl.svg,https://img.shields.io/npm/l/rctl.svg" \ 
#>          --destfiles "rctl1.svg,rctl2.svg"
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rdownload --urls "https://img.shields.io/npm/dm/rctl.svg,https://img.shields.io/npm/l/rctl.svg" \ 
#>          --destfiles "rctl1.svg,rctl2.svg"
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rdownload --urls "https://img.shields.io/npm/dm/rctl.svg,https://img.shields.io/npm/l/rctl.svg" \ 
#>          --destfiles "rctl1.svg,rctl2.svg" --max-cores 1
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
#> Error in download.file(link, config.download.tmp) : 
#>   cannot open URL 'https://raw.githubusercontent.com/Miachol/configr/master/inst/extdata/config.global.toml'
#> Calls: eval -> eval -> do.call -> <Anonymous> -> download.file
#> In addition: Warning message:
#> In download.file(link, config.download.tmp) :
#>   URL 'https://raw.githubusercontent.com/Miachol/configr/master/inst/extdata/config.global.toml': status was 'Couldn't connect to server'
#> Execution halted
#> Description:
#>  rconfig is an R-based tool to parse and generate configuration file.
#> 
#> Usage: 
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rconfig [options] [params]
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rconfig package.json
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rconfig -c package.json
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rconfig -f 2 test.json -e "config.dat=list(a=1, b=2), write.type='json'"
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rconfig -f "configr::write.config" test.json -e "config.dat=list(a=1, b=2), write.type='json'"
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rconfig -i test.json -r 'function(x){x[["a"]] + x[["b"]]}'
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rconfig -i test.json -r 'function(x){x[["a"]]}'
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rconfig -i test.json -r 'function(x){x[["b"]]}'
#> 
#> Options:
#>  -v, --verbose
#>                Print extra output [default FALSE]
#>  -f FUNC, --func=FUNC
#>                Index or name of used function [e.g. configr::read.config (1), configr::fetch.config (2), configr::write.config (3), configr::get.config.type (4)].
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
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rclrs [options] [params]
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rclrs default
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rclrs -t default
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rclrs -t default -r 'x[1]'
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rclrs -t red_blue
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rclrs --show-all-themes
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
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rmv [options] [params]
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rmv "`ls`" -e "do.rename = FALSE, prefix = 'prefix', suffix = 'suffix'"
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rmv "`ls`" -e "do.rename = FALSE, replace = list(old =c('-', '__'), new = c('_', '_'))"
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rmv "`ls`" -e "do.rename = FALSE, toupper = TRUE"
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rmv "`ls`" -e "do.rename = FALSE, tolower = TRUE"
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rmv -e "files_dir = '.', pattern = '.*.txt', do.rename=F, replace=list(old='old', new='new')"
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rmv "`ls`" -e "do.rename=T, replace=list(old='old', new='new')"
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
#> [1] "2020_04_21_20_06_47_" "2020_04_21_20_06_"    "2020_04_21_20_"      
#> [4] "2020_04_21_"          "2020_04_"             "2020_"               
#> 
#> [[2]]
#> [1] "2020-04-21-20-06-47-" "2020-04-21-20-06-"    "2020-04-21-20-"      
#> [4] "2020-04-21-"          "2020-04-"             "2020-"               
#> 
#> [[3]]
#> [1] "2020/04/21/20/06/47/" "2020/04/21/20/06/"    "2020/04/21/20/"      
#> [4] "2020/04/21/"          "2020/04/"             "2020/"               
#> 
#> [[4]]
#> [1] "2020@04@21@20@06@47@" "2020@04@21@20@06@"    "2020@04@21@20@"      
#> [4] "2020@04@21@"          "2020@04@"             "2020@"               
#> 
#> 2020_04_21_20_06_49_
#> 2020_04_21_20_06_
#> 2020_04_21_20_
#> 2020_04_21_
#> 2020_04_
#> 2020_
#> 2020_04_21_20_06_51_
#> [[1]]
#> [1] "2020_21"
#> 
#> [[2]]
#> [1] "2020-21"
#> 
#> [[3]]
#> [1] "2020/21"
#> 
#> [[4]]
#> [1] "2020@21"
#> 
#> [[1]]
#> [1] "2020_04_21_20_06_53_" "2020_04_21_20_06_"    "2020_04_21_20_"      
#> [4] "2020_04_21_"          "2020_04_"             "2020_"               
#> 
#> [[2]]
#> [1] "2020*04*21*20*06*53*" "2020*04*21*20*06*"    "2020*04*21*20*"      
#> [4] "2020*04*21*"          "2020*04*"             "2020*"               
#> 
#> Description:
#>  rtime_stamp is an R-based tool to generate time stamp.
#> 
#> Usage: 
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rtime_stamp [options] [params]
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rtime_stamp
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rtime_stamp -r 'x[[1]]'
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rtime_stamp -r 'x[[1]][1]'
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rtime_stamp -t '%Y_%d'
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rtime_stamp -e "extra_flat=c('-')"
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
#> B6LkGM4RYuXnkwu0ECwE
#> G9q20T6T7pWuNp8YRjK6zOgbUAWyGh
#> CNiG1JomKCZySyc6MZJZ
#> lwpwBCEMnZAZiZbo9rB5
#> WWfgOA50HVPCQs63jteL
```

``` bash
ranystr -h
#> Description:
#>  ranystr is an R-based tool to generate any counts and any length random strings.
#> 
#> Usage: 
#>  /Users/apple/.config/yarn/global/node_modules/.bin/ranystr [options] [params]
#>  /Users/apple/.config/yarn/global/node_modules/.bin/ranystr
#>  /Users/apple/.config/yarn/global/node_modules/.bin/ranystr -l 30
#>  /Users/apple/.config/yarn/global/node_modules/.bin/ranystr -l 20 -n 3
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
#> Copying ngstk bin/ demo_bin.sh to /Users/apple/.ngstk/bin
#> Please set /Users/apple/.ngstk/bin in your PATH to use the bin files.
#> Linux/Mac OS X: echo 'export PATH=$PATH:/Users/apple/.ngstk/bin\n' >> ~/.bashrc
#> R users: echo 'Sys.setenv(PATH="/usr/local/Cellar/spack/bin:/Users/apple/miniconda3/bin:/Users/apple/.yarn/bin:/Users/apple/.config/yarn/global/node_modules/.bin:/Users/apple/Documents/gopath/bin:/usr/local/Cellar/gnu-sed/4.7/bin:/usr/local/opt/gnu-sed/libexec/gnubin:/Users/apple/.go/bin:/usr/local/bin:/usr/local/opt/coreutils/libexec/gnubin:/usr/local/Cellar/spack/bin:/usr/local/Cellar/modules/4.3.0/bin:/Users/apple/miniconda3/bin:/Users/apple/.yarn/bin:/Users/apple/.config/yarn/global/node_modules/.bin:/Users/apple/repositories/local/mac/go/bin:/usr/local/Cellar/gnu-sed/4.7/bin:/usr/local/opt/gnu-sed/libexec/gnubin:/Users/apple/.go/bin:/usr/local/bin:/usr/local/opt/coreutils/libexec/gnubin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/Library/TeX/texbin:/usr/local/go/bin:/opt/X11/bin:/Users/apple/.ngstk/bin")\n' >> ~/.Rprofile
#> ngstk.demo_bin.sh 
#>              TRUE 
#> Copying ngstk bin/ demo_bin.sh to /private/tmp/path
#> Please set /private/tmp/path in your PATH to use the bin files.
#> Linux/Mac OS X: echo 'export PATH=$PATH:/private/tmp/path\n' >> ~/.bashrc
#> R users: echo 'Sys.setenv(PATH="/usr/local/Cellar/spack/bin:/Users/apple/miniconda3/bin:/Users/apple/.yarn/bin:/Users/apple/.config/yarn/global/node_modules/.bin:/Users/apple/Documents/gopath/bin:/usr/local/Cellar/gnu-sed/4.7/bin:/usr/local/opt/gnu-sed/libexec/gnubin:/Users/apple/.go/bin:/usr/local/bin:/usr/local/opt/coreutils/libexec/gnubin:/usr/local/Cellar/spack/bin:/usr/local/Cellar/modules/4.3.0/bin:/Users/apple/miniconda3/bin:/Users/apple/.yarn/bin:/Users/apple/.config/yarn/global/node_modules/.bin:/Users/apple/repositories/local/mac/go/bin:/usr/local/Cellar/gnu-sed/4.7/bin:/usr/local/opt/gnu-sed/libexec/gnubin:/Users/apple/.go/bin:/usr/local/bin:/usr/local/opt/coreutils/libexec/gnubin:/usr/local/bin:/usr/bin:/bin:/usr/sbin:/sbin:/Library/TeX/texbin:/usr/local/go/bin:/opt/X11/bin:/private/tmp/path")\n' >> ~/.Rprofile
#> ngstk.demo_bin.sh 
#>              TRUE 
#> Description:
#>  rbin is an R-based tool to collect R packages bin files.
#> 
#> Usage: 
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rbin [options] [params]
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rbin ngstk
#>  /Users/apple/.config/yarn/global/node_modules/.bin/rbin --destdir /usr/local/bin/ ngstk
#> 
#> Options:
#>  -v, --verbose
#>                Print extra output [default FALSE]
#>  -f FUNC, --func=FUNC
#>                Index or name of used function [e.g. ngstk::bin (1)].
#>  -c PKGS, --pkgs=PKGS
#>                Package names used to copy inst/bin directory files to PATH.
#>  -o DESTDIR, --destdir=DESTDIR
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

Please fork the [GitHub rctl
repository](https://github.com/openanno/rctl), modify it, and submit a
pull request to us.

## Maintainer

[Jianfeng Li](https://github.com/Miachol)

## License

MIT
