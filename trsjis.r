files <- list.files()
target <- grep("^code", files)
for(idx in target){
    sf <- paste("s_", files[idx], sep = "")
    dat <- readLines(files[idx])
    for(i in 1:length(dat)){
        dat[i] <- iconv(dat[i], to = "shift-jis")
    }
    writeLines(dat,sf)
}
