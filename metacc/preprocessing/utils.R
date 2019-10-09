convert_chr_format <- function(chr, to) {
  # Function to convert the chr from short to long format and viceversa
  # to: "short" or "long"
  chr <- as.character(chr)
  stopifnot(to %in% c("short","long"))
  short_alphabet <- c(1:19,"X","Y","MT")
  long_alphabet <- paste("chr",short_alphabet,sep="")
  if (to == "short") {
    if (all(chr %in% short_alphabet)) { 
      return(chr) 
    } else {
      stopifnot(all(chr %in% long_alphabet))
      names(short_alphabet) <- long_alphabet
      return(unname(short_alphabet[chr]))
    }
  }
  if (to == "long") {
    if (all(chr %in% long_alphabet)) { 
      return(chr) 
    } else {
      stopifnot(all(chr %in% short_alphabet))
      names(long_alphabet) <- short_alphabet
      return(unname(long_alphabet[chr]))
    }
  }
}