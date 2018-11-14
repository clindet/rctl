#!/usr/bin/env bash

workdir=${1}
extra_list="workdir=\"${workdir}\", logfn=\"${workdir}/logs/majiq_builder_step.log\""

rbashful -c ${workdir}/cli_tmp/${id}.parsed.yaml \
    --extra-list "${extra_list}" \
    --env-toml ${workdir}/env.toml \
    -n "majiq_builder_step" -v
