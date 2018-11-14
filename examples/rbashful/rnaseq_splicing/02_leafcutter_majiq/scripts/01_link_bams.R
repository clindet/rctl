#!/usr/bin/env Rscript

suppressMessages(if(!require('pacman')) install.packages('pacman'))
pkgs <- c("optparse", "configr", "stringr")
suppressMessages(pacman::p_load(pkgs, character.only = TRUE))
option_list <- list(
  make_option(c("--env-parsed-yaml"), default = "env.parsed.yaml", help = "Parsed ENV YAML file."),
  make_option(c("--task-name"), default = "majiq_splicing_01_link_bams", help = "Task name.")
)

usage <- "%prog --env-parsed-yaml env.parsed.yaml"
description <- sprintf("Description:\n%s",
                       "Script to create tempdir and store the temp bams for analysis."
                      )

opt_parser_obj <- OptionParser(option_list = option_list,
                        usage = usage, description = description)
opt <- parse_args(opt_parser_obj)

# Section to run majiq_splicing analysis

## set required tools
config <- read.config(opt[["env-parsed-yaml"]])
attach(config)

## set output dir and other parameters
task_name <- opt[["task-name"]]

rawbamdir <- input[['rawbamdir']]
logfn <- str_replace(submit_log, ".log$", sprintf('.%s.log', task_name))
logdir <- dirname(logfn)
if (!dir.exists(logdir)) dir.create(logdir)

samples_info_fn <- input[['samples_info_fn']]
ids <- readLines(samples_info_fn)
sapply(ids, function(x) {
  config_input_tmp <- parse.extra(config$input, extra.list = list(id=as.character(x)), global.vars.field = NULL)
  rawbamfile <- normalizePath(config_input_tmp$rawbamfile, mustWork =  FALSE)
  fn <- normalizePath(config_input_tmp$bamfile, mustWork = FALSE)
  if (!dir.exists(dirname(fn))) dir.create(dirname(fn), recursive = TRUE)
  if (file.exists(rawbamfile) && !file.exists(fn)) file.link(rawbamfile, fn) 
})
warnings()
