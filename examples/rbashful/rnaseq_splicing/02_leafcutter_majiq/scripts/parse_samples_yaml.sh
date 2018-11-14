#!/usr/bin/env bash

id=${1}
workdir=${2}
rawbamdir=${3}
extra_list="workdir=\"${workdir}\", logfn=\"${workdir}/logs/${id}.log\", rawbamdir = \"${rawbamdir}\""

rbashful -c ${workdir}/cli.yaml \
    --extra-list "${extra_list}" \
    --env-toml ${workdir}/env.toml \
    -p \
    -o ${workdir}/cli_tmp/${id}.parsed.yaml \
    -v
