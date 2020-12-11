# stolen from oro.dicom to convert unsigned int tags to character arrays
rawToHex <- function(bytes) toupper(paste(rev(bytes), collapse = ""))

read_raw <- function(fname, endian = "little") {
  # read the full file as raw
  fsize <- file.info(fname)$size
  return(readBin(fname, "raw", n = as.integer(fsize), endian = endian))
}

# extracted from the file eg tags <- list(spec_data = "7FE1,1010",
#                                         pat_name = "0010,0010")
# A named list is returned containing the raw binary data is returned.
#' @export
dicom_reader <- function(input, tags = list(sop_class_uid = "0008,0016"),
                         endian = "little", debug = FALSE) {
  
  # assume as file path if input is character, otherwise assume raw
  if (is.character(input)) {
    fsize <- file.info(input)$size
    fraw <- readBin(input, "raw", n = as.integer(fsize), endian = endian)
  } else {
    fraw <- input 
  }
  
  # skip the 128 byte preamble and read the next 4 bytes (should be "DICM")
  dicom_prefix <- rawToChar(fraw[129:132])
  
  if (dicom_prefix != "DICM") {
    stop("Could not find DICM prefix - input is probally not DICOM.")
  }
  
  # strip the preamble
  fraw <- fraw[133:length(fraw)]
  
  # prepare the output structure
  out <- tags
  out[] <- list(NULL)
  
  short_vrs <- c("OB","OW","SQ","UN","UT")
  long_vrs  <- c("AE","AS","AT","CS","DA","DS","DT","FL","FD","IS","LO",
                 "LT","OF","PN","SH","SL","SS","ST","TM","UI","UL","US")
  
  # current position in the data
  pos <- 0
  
  if (debug) debug_table <- NULL
  
  while (pos < length(fraw)) {
    
    # read the tag ID
    group   <- rawToHex(fraw[pos + 1:2])
    element <- rawToHex(fraw[pos + 3:4])
    tag_str <- paste(group, element, sep = ",") 
    
    # tentatively read in the VR - only used for explicit VRs
    vr <- rawToChar(fraw[pos + 5:6])
    
    # 5600,0020 is a hack that seems to be needed for Siemens SpectroscopyData
    if (vr %in% short_vrs | tag_str == "5600,0020") {
      # explicit VR with two bytes of zero padding
      length_raw <- fraw[pos + 9:12]
      length <- read_uint32(length_raw)
      # move pos to start of the value
      pos <- pos + 12
    } else if (vr %in% long_vrs) {
      # explicit VR without zero padding
      length_raw <- fraw[pos + 7:8]
      length <- readBin(length_raw, "integer", size = 2, signed = FALSE)
      pos <- pos + 8
    } else {
      # looks like implicit VR
      length_raw <- fraw[pos + 5:8]
      length <- read_uint32(length_raw)
      pos <- pos + 8
    }
    
    if (rawToHex(length_raw) %in% c("FFFFFFFF", "00000A00")) length <- 0 
      
    if (length == 0) {
      start_byte <- NA
      end_byte   <- NA
    } else {
      start_byte <- pos + 1
      end_byte   <- pos + length
    }
  
    if (debug) {
      tag <- paste("(", tag_str, ")", sep = "")
      x <- data.frame(tag, start_byte, end_byte, length)
      debug_table <- rbind(debug_table, x)
    }
    
    # check if current tag is the in list of those to be returned
    search_res <- grep(tag_str, tags) 
    if (length(search_res) == 1) {
      out[[search_res]] <- fraw[start_byte:end_byte]
    }
    
    pos <- pos + length
    
  }
  
  if (debug) print(debug_table)
  
  return(out)
}