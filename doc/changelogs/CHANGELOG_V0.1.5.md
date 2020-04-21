# rctl 0.1.5 ChangeLog

## 2018-11-16, Version 0.1.5

**New scripts:**  

`./bin/rtime_stamp`
`./bin/ranystr
`./bin/rdownload`

**New features:**

- support rmv only input the `extra` params
- add R package `./src/rctl` to resolve the R dependence and reduce the redundant R codes
- add `configr::fetch.config` in `./bin/rconfig`
- add a complete set of RNA splcing analysis pipeline
- './bin/rtime_stamp' support to generate time stamp
- `./bin/ranstr` support to generate random strings
- `./bin/rdownload` support to parallel download urls
- `rctl --list-all-subcmds` return all sub-commands split by space

**Test related:**

- add `./scripts/test/start_cli_test` to test the command line tools default output

**Minor bug fix:**

- rename `./bin/ngstk` subcomd `rcolors` to `rclrs`
- add `rmv` in `package.json` test
