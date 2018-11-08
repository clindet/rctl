# ngsjs 0.1.2 ChangeLog

## 2018-11-08, Version 0.1.2 (Current)

### Notable changes

**New scripts:**  

- `bin/rdeps`
- `bin/rsession`
- `bin/rbashful`
- `bin/rinstall`

**New features:**

- Support automatically install [ngsjs](https://github.com/JhuangLab/ngsjs) scripts required R packages
- Support `sessionInfo()` and `sessioninfo::session_info()` in `bin/rsession`
- Support [bashful](https://github.com/wagoodman/bashful) to style and parallelize the bash commands using `bin/rbashful`
- Support `devtools`, `BiocManager`, `BioInstaller` and `pacman` R packages (pkgs, repo, name) in `bin/rinstall`
