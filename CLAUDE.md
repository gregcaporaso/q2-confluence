# CLAUDE.md

Repository-local instructions for the q2-confluence plugin (registered name `confluence`, module `q2_confluence`); the global `~/.claude/CLAUDE.md` also applies.
The plugin is a nested repository inside the thermo-hec project, whose `CLAUDE.md` and memory hold the project context; git commands here use `git -C q2-confluence` from the project directory.

## Conventions set in review (developer, 2026-10-13)

 * Only registered actions are public functions; every other function, including usage-example functions, fixture factories, nested helpers and test helpers, starts with an underscore.
 * Every function body opens with `# TODO: this is AI-generated code; human review pending`, and every registered action's description ends with `EXPERIMENTAL: this is AI-generated code; human review pending.`; the developer deletes both as functions pass review. New functions and actions get the same markers.
 * An output name never repeats an input or parameter name: q2cli classifies signature entries by name and renders such an output as a second `--i-` option. Follow the stock convention of descriptive output names (`relative_frequency_table`, `collapsed_tree`, `integrated_metadata`); `tests/test_plugin_setup.py` enforces it. Render `qiime confluence <action> --help` for every action whose signature changes.
 * Metadata and TSV column names use hyphens (`original-sample-id`, `unassigned-share`); action, parameter and output identifiers stay Python identifiers. `sample-id` is a reserved metadata ID header and cannot be a column name (query-reads uses `sample`).
 * Merged sample IDs are `<id>-<key>` with `original-sample-id` and `profile-type` columns (renamed from `source` on 2026-10-16); pairing relies on identical original IDs across sources.
 * Provisional actions are registered with `deprecated=True` and a description naming their intended home (q2-feature-table, q2-phylogeny, MOSHPIT).
 * Stock actions over custom code; one concern per function and per commit; no rerun-only modes. Every action keeps one usage example, and `tests/test_examples.py` executes them all.
 * The workflow runs on proportions: no counts conversion and no rarefaction anywhere (the `frequency` action was dropped on 2026-10-13); `relative-frequency` renormalizes after feature filters, and the paired report computes its distances once on the relative table.
 * Commit messages name @gregcaporaso as the reviewer who requested the change and say what Claude implemented, with the `Co-Authored-By` trailer; never push.

## Environment and checks

 * Named conda env `q2-confluence-dev` (`/opt/homebrew/Caskroom/miniforge/base/envs/q2-confluence-dev`, plugin installed editable). Use it without activation: `PATH=$D/bin:$PATH CONDA_PREFIX=$D MPLCONFIGDIR=$TMPDIR/mpl NUMBA_CACHE_DIR=$TMPDIR/numba`.
 * The q2cli cache under that env is read-only in the sandbox: for `qiime ... --help` and `qiime dev refresh-cache`, unset `CONDA_PREFIX` and set `HOME=$TMPDIR/home`; q2cli then caches under that HOME.
 * Before every commit: `flake8`, `q2lint`, `python -m pytest -q` (about five minutes; run only the affected test module while iterating). Run the full suite in the background and do not edit files while it runs.
 * `sed -i ''` on macOS has no `\b`; use Python for word-boundary edits. Tokenize-based renames shift continuation-line indents (flake8 E128); check them.
 * q2cli does not create missing output directories; a page or script must `mkdir` them.
 * The docs build (`make html`) needs a local port bind and network access the sandbox blocks; the developer builds it.

## Ownership

 * `docs/` prose is the developer's once they write it; until then `docs/tutorial.md` and `docs/data-preparation.md` are Claude's drafts and may be updated for signature changes. `docs/intro.md` is a placeholder. The tutorial's commands are describe-usage blocks (rendered for every interface by q2doc); validate them by rendering through q2cli's CLIUsage driver, and keep the shell alternatives in the tip boxes in step with them.
 * Code, tests, examples, action help text and `README.md` are Claude's to keep accurate in the same commit as the change they describe.
 * The thermo-hec project's `slurm/tutorial.sbatch` renders this repository's tutorial page to CLI commands with q2cli's usage driver and runs them on bacon; a change to the page or to an action's CLI is tested there (re-ship the plugin as a `git archive` tarball at the commit under test).
