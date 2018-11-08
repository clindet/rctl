#!/usr/bin/env bash

id=12
workdir=~/Documents/repositories/ljf/github/ngsjs/tests/rbashful/rnaseq_splicing
# Run the default cmd in env.toml
rbashful -c ${workdir}/cli.yaml \
    --env-toml ${workdir}/env.toml \
    --cmd-name default \
    -v

# Parse the bashful YAML file (Output cli.parsed.yaml) 
# and then run the default cmd in env.toml
rbashful -c ${workdir}/cli.yaml \
    --env-toml ${workdir}/env.toml \
    --extra-list "id=${id},logfn=\"${id}.log\"" \
    --parse-cli-yaml \
    -v

# Running the parsed yaml (given tags cmds)
rbashful -c ${workdir}/cli.parsed.yaml \
    --env-toml ${workdir}/env.toml \
    --cmd-name start_tags \
    --extra-list "tags='star_alignment'" -v
