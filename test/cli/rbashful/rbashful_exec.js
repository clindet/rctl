#!/usr/bin/env node

var spawnSync = require('child_process').spawnSync;

function rbashful(cmdStr = "rbashful") {
    var rbashful = spawnSync(cmdStr);
    return(rbashful)
}

module.exports = rbashful;
