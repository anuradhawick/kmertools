# Contributing to kmertools project

We love to have your contributions to the kmertools project, whether it's:
* Reporting a bug
* Submitting a fix
* Proposing new features

## Clone and install kmertools onto your machine

First, make sure you have [git](https://github.com/git-guides/install-git) and [rust](https://www.rust-lang.org/tools/install) installed on your machine.

On GitHub, [fork](https://docs.github.com/en/get-started/quickstart/fork-a-repo) the kmertools repository and [clone](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) it to your machine.

```shell
# clone repository to your local machine
git clone https://github.com/anuradhawick/kmertools.git
```

Move to the kmertools directory 

```shell
cd kmertools
```

Now build kmertools using cargo. Make sure to have [`rust`](https://docs.conda.io/en/latest/) installed.

```shell
cargo build
```

## Test kmertools build

Use the following command to test the build. All tests should pass.

```shell
cargo test
```

## Coding Style

We use [Clippy](https://doc.rust-lang.org/clippy/) to lint code in kmertools.

Before committing, make sure to run Clippy as follows.

```shell
cargo clippy
```

## Report bugs using GitHub's issues

We use GitHub issues to track public bugs. Report a bug by opening a new issue in GitHub [issues](https://github.com/anuradhawick/kmertools/issues). You will get to select between templates for bug reports and feature requests. If none of these templates matches what you want to report, you can use the custom issue template.

## Committing code

Once you have finished coding and all the tests pass, commit your code and make a pull request. 

```bash
# Add changed/added files
git add <file name>

# Commit changes
git commit -m "<commit message>"

# Push changes
git push
```

Make sure to follow the commit style of [c3dev](https://github.com/cogent3/c3dev/wiki#style-for-commit-messages). Relevant prefixes are replicated below for convenience.

| **Commit Prefix** | **For**                                       |
|-------------------|-----------------------------------------------|
| DEV:              | development tool or utility                   |
| DOC:              | documentation                                 |
| TST:              | addition or modification of tests             |
| REL:              | related to a release                          |
| MAINT:            | maintenance commit (refactoring, typos, etc.) |
| FIX:              | fix for bugs                                  |
| GIT:              | git related                                   |
| REV:              | revert an earlier commit                      |


Your contribution will be reviewed before accepting it. 

## License

By contributing, you agree that your contributions will be licensed under the GPL-3.0 license.

## References

This document was adapted from the open-source contribution guidelines for [Transcriptase](https://github.com/briandk/transcriptase-atom/blob/master/CONTRIBUTING.md) and [c3dev](https://github.com/cogent3/c3dev/wiki/How-to-Contribute-Code).
