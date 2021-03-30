### FUNCTIONS

# round values while preserving their sum
# from: https://stackoverflow.com/questions/32544646/round-vector-of-numerics-to-integer-while-preserving-their-sum
smart_round <- function(x) {
  y <- floor(x)
  indices <- tail(order(x-y), round(sum(x)) - sum(y))
  y[indices] <- y[indices] + 1
  y
}

# bray-curtis similarity
bray_curtis <- function(p,q) sum(pmin(p,q))

# jaccard similarity

# sequence similarity
seq_sim <- function(seq1,seq2) {
  
  seq1 <- strsplit(seq1,"")[[1]]
  seq2 <- strsplit(seq2,"")[[1]]
  
  if (length(seq1) == length(seq2)) return(sum(seq1==seq2)/length(seq1))
  else {
    
    if (length(seq1)>length(seq2)) {
      long_seq <- seq1
      short_seq <- seq2
    } else {
      long_seq <- seq2
      short_seq <- seq1
    }
    
    length_diff <- length(long_seq) - length(short_seq)
    sim <- rep(NA,length_diff+1)
    for (i in 1:(length_diff+1)) {
      long_seq_i <- long_seq[i:(i+length(short_seq)-1)]
      sim[i] <- seq_sim(short_seq,long_seq_i)
    }
    
    return(max(sim))
    
  }
  
}

# sum rows of OTU table
sumOTUs <- function(otus_to_sum) {
  out <- cbind(X = otus_to_sum$X[1], ### FIXME: something more elegant could be done with the sequences?
               data.frame(t(colSums(otus_to_sum[,2:ncol(otus_to_sum)]))))
  return(out)
}

# merge taxa based on similarity
mergeTaxa <- function(otus,taxa,similarity_threshold=0.995) {
  
  # collapse taxa with sequence similarity over threshold
  i <- 1
  while (i<nrow(otus)) {
    
    # what entries of the taxonomy table match the i-th row? (same taxonomy down to the genus and sequence similarity over 99%)
    sim <- sapply(otus$X[(i+1):nrow(otus)],FUN=function(x) seq_sim(x,otus$X[i]))
    n <- which(sim > similarity_threshold &
               taxa$Kingdom[(i+1):nrow(taxa)] == taxa$Kingdom[i] &
               taxa$Phylum[(i+1):nrow(taxa)] == taxa$Phylum[i] &
               taxa$Class[(i+1):nrow(taxa)] == taxa$Class[i] &
               taxa$Order[(i+1):nrow(taxa)] == taxa$Order[i] &
               taxa$Family[(i+1):nrow(taxa)] == taxa$Family[i] &
               taxa$Genus[(i+1):nrow(taxa)] == taxa$Genus[i])
    n <- n + i
    
    # if there are matches, resolve them: sum abundances in OTU table and remove redundant entries
    if(length(n)>0) {
      taxa <- taxa[-n,]
      otus[i,] <- sumOTUs(otus[c(i,n),])
      otus <- otus[-n,]
    }
    
    # advance to the next row
    i <- i+1
    
    # print progress
    print(paste(100*i/nrow(otus),'%',collapse=" "))
    
  }
  
  return(list(otus=otus,taxa=taxa))

}





