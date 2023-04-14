library(data.table)
library(tidyverse)



# Batch1 Markerlengths ----------------------------------------------------
dir_list <-  list.files("../data/batch1/",
                        pattern= "_R$")


# create a table to which others can be joined 
markerlengths <- NULL

for (dir in dir_list) {
    
    subject <-  str_split(dir, "_") %>% purrr::map(1) %>% unlist   
    
    ## path to poly-G raw data directory (marker length files)
    marker_dir <- paste0("../data/batch1/", dir, "/repre_repli_data/")
    
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
        l$sample <- as.character(l$sample)
        l$subject <- subject
        l
    }
    
    markers <- get_marker_lengths(subject,marker_dir) %>% 
      mutate(sample = str_sub(sample, end = -3))
    
    markers$sample <- gsub(subject,'',markers$sample)
    
    # subtracting the length of the normal sample all other samples, per marker and subject
    markers <- markers %>% 
        group_by(subject, marker) %>% 
        mutate(length=length-length[str_detect(sample, "^N")][1]) %>% 
        filter(!str_detect(sample, "^N")) %>% 
        ungroup %>% 
        as.data.table
    
    markerlengths_int  <- dcast(sample ~ marker, value.var='length', data=markers)
    markerlengths_int$subject <- subject
    markerlengths <- bind_rows(markerlengths, markerlengths_int)
    
}  

# adding subject name to the sample
markerlengths$sample <- paste0(markerlengths$sample, "_", markerlengths$subject)
markerlengths_batch1 <- markerlengths


# Batch2 Markerlengths ----------------------------------------------------

dir_list <-  list.files("../data/batch2/",
                        pattern= "_R$")
normal_samples <- c("CC001_5|T790M_III_1_N|CC002_11|CC005_11|WGS00C_N20|WGS00E_N22|CC006_3|WGS00D_N8")

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
    l$sample <- as.character(l$sample)
    l$subject <- subject
    l
  }
  markers <- get_marker_lengths(subject,marker_dir) %>% 
    mutate(sample = str_sub(sample, end = -3))
  
  root_sample <- select_normal(subject)
  
  markers <- markers %>% 
    group_by(marker) %>% 
    filter(any(str_detect(sample, normal))) %>% 
    mutate(length = length - length[str_detect(sample, normal)][1]) %>% 
    filter(!str_detect(sample, normal)) %>% 
    as.data.table
  
  markerlengths_int  <- dcast(sample ~ marker, value.var='length', data=markers)
  markerlengths_int$subject <- subject
  markerlengths <- bind_rows(markerlengths, markerlengths_int)
  
}  

markerlengths_batch2 <- markerlengths

# reorder names to make sample name first for easier processing
split_names <- str_split(markerlengths_batch2$sample, "_")
markerlengths_batch2$sample <- paste(map(split_names, 2), map(split_names, 1), sep = "_")

# Joinig batch1 and batch2 markerlengths ----------------------------------

joined_lengths <- bind_rows(markerlengths_batch1, markerlengths_batch2)

## get all pairwise (unique) combinations only for primary samples, and then calculate the variance in the differences between their mean lengths
samples <- unique(joined_lengths$sample)
combos_wide <- combn(samples, m= 2) %>% as.data.frame()

## make combo table longer
combos_long <- data.frame(s1=as.character(combos_wide[1,]), s2=as.character(combos_wide[2,]))

## remove all combintions where both primaries come from the same sample
x <- combos_long %>% 
    separate(s1, c(NA, "sample_a"), sep="_", remove = FALSE) %>% 
    separate(s2, c(NA, "sample_b"), sep="_", remove = FALSE) %>% 
    mutate(sample_a=ifelse(sample_a=="WGS004", "WGS00D", sample_a),
           sample_b=ifelse(sample_b=="WGS004", "WGS00D", sample_b),
           sample_a=ifelse(sample_a=="WGS005", "WGS00E", sample_a),
           sample_b=ifelse(sample_b=="WGS005", "WGS00E", sample_b)) %>% 
    filter(sample_a!=sample_b,!str_detect(s1, "N"), !str_detect(s2, "N")) %>% 
    dplyr::select(s1, s2) 

## moving sample names to the rownams 
m <- joined_lengths %>% 
    dplyr::select(-subject) %>% 
    column_to_rownames("sample")


    ## get tumor-normal length (for each tumor) and calculate their pearson correlation
    get_correlation <- function(i, x, m) {
        s1 <- x$s1[i]
        s2 <- x$s2[i]
        mean_lengths1 <- as.numeric(m[s1,])
        mean_lengths2 <- as.numeric(m[s2,])

        # count on how markers observation the correlation is based 
        one_minus_two <- mean_lengths1 - mean_lengths2
        marker_number <- sum(!is.na(one_minus_two))

        if (marker_number>0) {
            r <- cor(mean_lengths1, mean_lengths2, method='pearson', use="complete.obs")
        list(a=s1,b=s2,cor=r,markers=marker_number)
        }
    }

# applying correlation
    result_list <- lapply(1:nrow(x), get_correlation, x, m)
    result <- rbindlist(result_list)
    
    result_cor_15_markers <- result[markers>15] 
 
write_tsv(result_cor_15_markers, "../results/interpatient_cor.tsv")

