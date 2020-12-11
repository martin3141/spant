# stolen from oro.dicom to convert unsigned int tags to character arrays
rawToHex <- function(bytes) toupper(paste(rev(bytes), collapse = ""))

# A very simple DICOM reader. tags is a named list of the tag data to be
# extracted from the file eg tags <- list(spec_data = "7FE1,1010",
#                                         pat_name = "0010,0010")
# A named list is returned containing the raw binary data is returned.
dicom_reader <- function(fname, tags) {
  
  # read the full file as raw
  fsize <- file.info(fname)$size
  fraw  <- readBin(fname, "raw", n = as.integer(fsize), endian = "little")
  
  # skip the 128 byte preamble and read the next 4 bytes (should be "DICM")
  dicom_prefix <- rawToChar(fraw[129:132])
  
  if (dicom_prefix != "DICM") stop("File is not DICOM.")
  
  # strip the preamble
  data_set <- fraw[133:fsize]
  
  # prepare the output structure
  out <- tags
  out[] <- list(NULL)
  
  short_vrs <- c("OB","OW","SQ","UN","UT")
  long_vrs  <- c("AE","AS","AT","CS","DA","DS","DT","FL","FD","IS","LO",
                 "LT","OF","PN","SH","SL","SS","ST","TM","UI","UL","US")
  
  # current position in the data
  pos <- 0
  
  while (pos < (fsize - 133 + 1)) {
    
    # read the tag ID
    group   <- rawToHex(data_set[pos + 1:2])
    element <- rawToHex(data_set[pos + 3:4])
    
    # tentatively read in the VR - only used for explicit VRs
    vr <- rawToChar(data_set[pos + 5:6])
    if (vr %in% short_vrs) {
      # explicit VR with two bytes of zero padding
      length_raw <- data_set[pos + 9:12]
      length <- read_uint32(length_raw)
      # move pos to start of the value
      pos <- pos + 12
    } else if (vr %in% long_vrs) {
      # explicit VR without zero padding
      length_raw <- data_set[pos + 7:8]
      length <- readBin(length_raw, "integer", size = 2, signed = FALSE)
      pos <- pos + 8
    } else {
      # looks like implicit VR
      length_raw <- data_set[pos + 5:8]
      length <- read_uint32(length_raw)
      pos <- pos + 8
    }
    
    tag_str <- paste(group, element, sep = ",") 
    
    # check if current tag is the in list of those to be returned
    search_res <- grep(tag_str, tags) 
    if (length(search_res) == 1) {
      start_byte <- pos + 1
      end_byte   <- pos + length
      out[[search_res]] <- data_set[start_byte:end_byte]
    }
    
    pos <- pos + length
    
    # cat(paste(group, element, start_byte, end_byte, length, "\n"))
  }
  return(out)
}