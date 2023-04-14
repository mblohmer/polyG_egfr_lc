# loading needed libraries
library(tidyverse) 
library(data.table)

dir_list <-  list.files("../data/batch2/",
                        pattern= "_R$")

normal_samples_vector <- c("CC001_5", "T790M_III_1_N", "T790M_III_1_N2", "CC002_11", "CC005_11", "WGS00C_N20", 
                      "WGS00E_N22", "CC006_3", "WGS00D_N8")

# create a table to which others can be joined 
markerlengths <- NULL

# samples for the batch2 were not renamed, thus not all normal start with "N"
## this function identifies the correct normal samples 

select_normal <- function(subject) {
  if (subject=="CC001") {normal_sample <- "CC001_5"
  } else if (subject=="T790M") {normal_sample <- "T790M_III_1_N"
  } else if (subject=="CC002") {normal_sample <- "CC002_11"
  } else if (subject=="CC005") {normal_sample <- "CC005_11"
  } else if (subject=="WGS00C") {normal_sample <- "WGS00C_N20"
  } else if (subject=="WGS00E") {normal_sample <- "WGS00E_N22"
  } else if (subject=="CC006") {normal_sample <- "CC006_3"
  } else if (subject=="WGS00D") {normal_sample <- "WGS00D_N8"}
  normal_sample
}

for (dir in dir_list) {
  
  subject <-  str_split(dir, "_") %>% purrr::map(1) %>% unlist   

  #select normal
  normal <- select_normal(subject)
  normal <- paste0("^", normal)
  
  ## path to poly-G raw data directory (marker length files)
  marker_dir <- paste0("../data/batch2/", dir, "/repre_repli_data/")
  
  ## load marker lengths and get the average length of each marker in each sample
  get_marker_lengths <- function(subject,marker_dir) {
    message(subject)
    
    f <- function(file) {
      x <- fread(file)
      x$length <- 1:nrow(x)
      x <- melt(x, id.var='length')
      toavg <- function(x) {
        avg <- sum(x$length * x$value) / sum(x$value)
        list(avg=avg)
      }
      x <- x[,toavg(.SD),by=variable]
      marker <- tail(strsplit(file,'[/]')[[1]],1)
      marker <- strsplit(marker,'_')[[1]][1]
      x$marker <- marker
      x
    }
    marker_files <- dir(marker_dir, full.names=T)
    l <- lapply(marker_files, f)
    l <- rbindlist(l)
    l <- l[,c(1,3,2),with=F]
    names(l) <- c('sample','marker','length')
    #l$sample <- gsub(subject,'',l$sample)
    l$sample <- as.character(l$sample)
    #l$sample <- sapply(l$sample, function(s) strsplit(s,'_')[[1]][1])
    l$subject <- subject
    l
  }
  markers <- get_marker_lengths(subject,marker_dir)
  
  root_sample <- select_normal(subject)

  markers <- markers %>% 
    group_by(marker) %>% 
    filter(any(str_detect(sample, normal))) %>% 
    mutate(length = length - length[str_detect(sample, normal)][1]) %>% 
    as.data.table
  
  markers$subject <- subject

  markerlengths <- bind_rows(markerlengths, markers) 
  
}  

markerlengths <- markerlengths %>% 
  mutate(sample = str_sub(sample, end = -3))

# Calculating correlations ---------------------------------------------------

# find samples
samples <- markerlengths$sample %>% unique

combos_wide <- combn(samples, m= 2) %>% as.data.frame()

## make combo table longer
combos_long <- data.frame(a=as.character(combos_wide[1,]), b=as.character(combos_wide[2,]))

## only include combinations within the same sample
combos <- combos_long %>% 
  separate(a, c("sample_a", NA), sep="_", remove = FALSE) %>% 
  separate(b, c("sample_b", NA), sep="_", remove = FALSE) %>% 
  filter(sample_a==sample_b, !(a %in% normal_samples_vector),
         !(b %in% normal_samples_vector)) 

# function for cor
get_cor_for_combination <- function(i, combos, markerlengths) {
  
  sample_a <- combos$a[i]
  sample_b <- combos$b[i]
  patient <- combos$sample_a[i]
  
  markerlengths_a  <-  markerlengths %>% 
    filter(sample == sample_a)
  
  markerlengths_b  <-  markerlengths %>% 
    filter(sample == sample_b)
  
  # find markers that are in both samples
  common_markers <- intersect(markerlengths_a$marker, markerlengths_b$marker)
  
  length_a  <-  markerlengths_a %>% 
    filter(marker %in% common_markers) %>% 
    arrange(marker) %>% 
    pull(length)
  
  length_b  <-  markerlengths_b %>% 
    filter(marker %in% common_markers) %>% 
    arrange(marker) %>% 
    pull(length)
  
  cor <- cor(length_a, length_b)
  
  list(a=sample_a, b=sample_b, cor=cor, patient=patient)
}

# calculating the correlations 
cor_tbl <- lapply(1:nrow(combos), get_cor_for_combination, combos, markerlengths)
cor_tbl <- bind_rows(cor_tbl)

write_tsv(cor_tbl, "../results/batch2_cor.tsv")
