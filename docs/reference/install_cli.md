# Install the spant command-line interface scripts to a system path.

This should be run following each new install of spant to ensure
consistency. Typical command line usage : sudo Rscript -e
"spant::install_cli()"

## Usage

``` r
install_cli(path = NULL)
```

## Arguments

- path:

  optional path to install the scripts. Defaults to : "/usr/local/bin".
