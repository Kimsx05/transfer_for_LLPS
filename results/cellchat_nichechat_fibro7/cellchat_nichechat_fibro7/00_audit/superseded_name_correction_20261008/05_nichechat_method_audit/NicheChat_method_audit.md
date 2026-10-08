# NicheChat implementation audit

Audit date: 2026-10-07 (Asia/Shanghai). Status: **AUDIT_ONLY**. This audit is completed before any NicheChat inference.

No verifiable biological cell-communication implementation named NicheChat was found in the audited environment. This is a scoped search result, not a claim that no such method exists anywhere.

## Search evidence

- Project and historical outputs: `/home/jinshengxi/LLPS` and `/home/jinshengxi/bladder_metabolism`; R/Python/shell scripts, Markdown, README and package DESCRIPTION files. No NicheChat script, repository reference, method documentation, or prior result was identified. The search excludes this newly generated result directory. See `local_project_search.txt` and `../00_audit/project_file_inventory.txt`.
- Installed R packages: all project libraries and default libraries enumerated in `../00_audit/installed_R_packages.tsv`; no NicheChat package.
- Installed Python distributions: 68 distributions in `../00_audit/installed_python_packages.json`; no NicheChat package.
- Broader local package/repository filenames: `/home/jinshengxi`, `/usr/local/lib`, `/opt`, excluding Codex/cache files; see `../00_audit/method_file_search.txt`. Some system-owned directories under `/opt` were unreadable; they are not assumed to contain a method.
- Public queries: `"NicheChat" cell communication github`, `"NicheChat" single cell ligand receptor`, `site:github.com "NicheChat"`.
- [GitHub repository search](https://api.github.com/search/repositories?q=NicheChat) returned two repositories (raw response in `github_repository_search.json`). [niyoseris/nichechat](https://github.com/niyoseris/nichechat) documents a Flask/Socket.IO social chat application, unrelated to biology. [Superunknown1103/NicheChat](https://github.com/Superunknown1103/NicheChat) is described in repository metadata as an app connecting people with common interests. Neither is an eligible cell-communication method.
- The user has been asked for the intended NicheChat repository, paper or local source path; no method source was provided at the time this audit was written.

## Requested method properties

| Property | Verified result |
|---|---|
| Biological package / repository name | NOT IDENTIFIED |
| Biological source URL | NOT IDENTIFIED |
| Language | UNKNOWN |
| Version / commit | UNKNOWN |
| Input format | UNKNOWN |
| Spatial information required | UNKNOWN |
| Communication / niche score definition | UNKNOWN |
| scRNA support | UNKNOWN |
| Spatial support | UNKNOWN |
| Joint scRNA + spatial support | UNKNOWN |
| Ligand–receptor output | UNKNOWN |
| Sender / receiver scoring | UNKNOWN |
| Per-sample inference support | UNKNOWN |
| Compatibility with current scRNA data | NOT ESTABLISHED |
| Method parameters / significance cutoff | NOT ESTABLISHED |

## Decision and missing information

NicheChat inference is stopped at the audit stage. Required to continue: the exact biological repository/package or paper, verifiable source version, documented score definition, input requirements, official inference API, and significance rules. Installing a similarly named social chat package would not resolve this ambiguity.

NicheNet, CellChat, CellPhoneDB, LIANA, COMMOT and other communication methods are not used as replacements. NicheChat sample coverage is 0; attempted inference samples 0; inference failures 0 because inference was not attempted. Figures C/D and cross-method support classifications cannot be produced until two verified methods have run. Unavailable NicheChat results will be recorded as NA/NOT_RUN, never as negative support or CELLCHAT_ONLY.

CellChat can continue independently once its dependency installation is approved. If a valid NicheChat source is later provided, append a separate validated-method audit before running it; do not reinterpret this initial audit as method validation.
