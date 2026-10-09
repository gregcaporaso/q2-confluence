# la-confluence

A [QIIME 2](https://qiime2.org) plugin by Greg Caporaso (jgc@caplab.dev) that integrates 16S rRNA amplicon and shotgun metagenome taxonomic profiles on one GTDB release.
The registered plugin name is `confluence` (`qiime confluence --help`); the Python package is `la_confluence`.

Status: in development, not yet released.
Actions marked deprecated at registration are provisional and name the plugin they are expected to move to.

## Development installation

Create the environment (QIIME 2 2026.7 `qiime2` distribution, sylph, and the documentation toolchain), then install the plugin into it:

```shell
conda env create -n la-confluence-dev --file environment-files/la-confluence-qiime2-dev.yml
conda activate la-confluence-dev
make dev
qiime dev refresh-cache
qiime confluence --help
```

On Apple Silicon, prefix the `conda env create` with `CONDA_SUBDIR=osx-64` and run `conda config --env --set subdir osx-64` after activating.

Run the tests and the linters with:

```shell
make test
make lint
```

## Documentation

The documentation is a MyST book built with q2doc, in the layout of the [gut-to-soil tutorial](https://github.com/qiime2/gut-to-soil-tutorial), and is published on ReadTheDocs from `.readthedocs.yaml`.
Build it locally, in the environment above with the plugin installed, with:

```shell
q2doc refresh-cache
make html
make serve
```

## About

`la-confluence` was [created from the QIIME 2 plugin template](https://develop.qiime2.org/en/stable/plugins/tutorials/create-from-template.html) on 2026-10-09.
It is a community plugin: see the [QIIME 2 Library](https://library.qiime2.org) for support options and the [QIIME 2 Forum](https://forum.qiime2.org) for help.
