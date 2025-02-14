# Contributing to spant development

Contributions to spant are most welcome and generally take the form of bug reports/feature requests, documentation and suggestions for code changes via a pull request.

## Bug reports/feature requests

Please report bugs and make feature requests via <https://github.com/martin3141/spant/issues/>.

When filing an issue, the most important thing is to include a minimal reproducible example so that we can quickly verify the problem, and then figure out how to fix it.

## Documentation

If you would like to help new users with a particular aspect of MRS data processing or analysis please consider contributing a short document in R markdown format to be added to the package. For examples please see the doc folder in the main GitHub repository.

Function reference documentation is derived from roxygen2 code comments and are best contributed via a pull request, as described below.

## Pull requests

To contribute a change to spant, follow these steps:

1. Create a branch in git and make your changes.
1. Push branch to github and issue pull request (PR).
1. Discuss the pull request.
1. Iterate until either we accept the PR or decide that it's not a good fit for spant.

If you're not familiar with git or github, please start by reading <http://r-pkgs.had.co.nz/git.html>. Please use the following style <http://adv-r.had.co.nz/Style.html>.

Please ensure all checks pass with `devtools::check(args = c('--as-cran'))` and consider adding one or more unit tests to the repository tests directory to confirm expected behavior with `devtools::test()`. CI checks are set to run following pushes to the master branch using GitHub Actions, see <https://orchid00.github.io/actions_sandbox> for details.