# ngsjs 0.1.5 ChangeLog

## 2018-11-12, Version 0.1.5 (Current)

**New scripts:**  

`./bin/rtime_stamp`
`./bin/ranystr
`./bin/rdownload`

**New features:**

- support rmv only input the `extra` params
- add R package `./src/ngsjs` to resolve the R dependence and reduce the redundant R codes
- add `configr::fetch.config` in `./bin/rconfig`
- add a complete set of RNA splcing analysis pipeline
- './bin/rtime_stamp' support to generate time stamp
- `./bin/ranstr` support to generate random strings
- `./bin/rdownload` support to parallel download urls

**Minor bug fix:

- rename `./bin/ngstk` subcomd `rcolors` to `rclrs`
- add `rmv` in `package.json` test
