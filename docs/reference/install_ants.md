# Install ANTs / ANTsX from : https://github.com/ANTsX/ANTs/releases

Install ANTs / ANTsX from : https://github.com/ANTsX/ANTs/releases

## Usage

``` r
install_ants(platform, version = "2.6.5")
```

## Arguments

- platform:

  see the releases page for supported platforms. Common platforms
  include : "macos-14-ARM64-clang", "almalinux9-X64-gcc",
  "centos7-X64-gcc" or "ubuntu-24.04-X64-gcc". Note, It may be necessary
  to increase the timeout with "options(timeout = 1000)" on slower
  connections.

- version:

  ANTs version to install, defaults to "2.6.5".
