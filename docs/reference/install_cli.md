# Install the spant command-line interface scripts to a system path.

This should be run following each new install of spant to ensure
consistency. Typical command line usage : sudo Rscript -e
"spant::install_cli()" Note sudo resets the HOME environment variable by
default, which can cause "there is no package called 'spant'" errors if
spant is only installed in your personal library. Use "sudo -E Rscript
-e \\spant::install_cli()\\" to preserve HOME and avoid this issue.

## Usage

``` r
install_cli(path = NULL)
```

## Arguments

- path:

  optional path to install the scripts. Defaults to : "/usr/local/bin".
