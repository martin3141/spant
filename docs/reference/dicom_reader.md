# A very simple DICOM reader.

Note this reader is very basic and does not use a DICOM dictionary or
try to convert the data to the correct datatype. For a more robust and
sophisticated reader use the oro.dicom package.

## Usage

``` r
dicom_reader(
  input,
  tags = list(sop_class_uid = "0008,0016"),
  endian = "little",
  debug = FALSE
)
```

## Arguments

- input:

  either a file path or raw binary object.

- tags:

  a named list of tags to be extracted from the file. eg tags \<-
  list(spec_data = "7FE1,1010", pat_name = "0010,0010")

- endian:

  can be "little" or "big".

- debug:

  print out some debugging information, can be "little" or "big".

## Value

a list with the same structure as the input, but with tag codes replaced
with the corresponding data in a raw format.
