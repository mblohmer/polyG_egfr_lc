
# loading needed libraries
library(tidyverse) 
library(data.table)
library(readxl)

dir_list <-  list.files("data/batch1/",
                        pattern= "_R$")

# create a table to which others can be joined 
markerlengths <- NULL

for (dir in dir_list) {
  
  subject <-  str_split(dir, "_") %>% purrr::map(1) %>% unlist   

  #select normal
  #normal <- select_normal(subject)
  #normal <- paste0("^", normal)
  
  ## path to poly-G raw data directory (marker length files)
  marker_dir <- paste0("data/batch1/", dir, "/repre_repli_data/")
  
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
    l$sample <- gsub(subject,'',l$sample)
    l$sample <- as.character(l$sample)
    l$sample <- sapply(l$sample, function(s) strsplit(s,'_')[[1]][1])
    l$subject <- subject
    l
  }
  markers <- get_marker_lengths(subject,marker_dir)

  normal <- str_subset(unique(markers$sample), "N")

  markers <- markers %>% 
    group_by(marker) %>% 
    filter(any(str_detect(sample, normal))) %>% 
    mutate(length = length - length[str_detect(sample, normal)][1]) %>% 
    as.data.table
  
  markers$subject <- subject

  markerlengths <- bind_rows(markerlengths, markers) 
  
}  

# subjects where renamed during data generation and their name needs to be changed back 
## to the original subject id 
rename_sample <- Vectorize(function(sample) {
  
  if (sample=="WGS002") return("WGS00B") 
  else if (sample=="WGS004") return("WGS00D")
  else if (sample=="WGS005") return("WGS00E")
  else if (sample=="CC003") return("CC003")
  
  
})

markerlengths <- markerlengths %>% 
  mutate(subject=rename_sample(subject), sample=paste0(subject, sample)) 
  

# Continue with the rest of the script ------------------------------------

# find number of markers
n_markers <- length(unique(markerlengths$marker))

# Calculating correlations ---------------------------------------------------

# find samples
samples <- markerlengths$sample %>% unique

combos_wide <- combn(samples, m= 2) %>% as.data.frame()

## make combo table longer
combos_long <- data.frame(a=as.character(combos_wide[1,]), b=as.character(combos_wide[2,]))

## only include combinations within the same sample

combos <- combos_long %>% 
  mutate(sample_a = str_sub(a, end=-3), sample_b = str_sub(b, end=-3)) %>% 
  filter(sample_a==sample_b, !str_detect(a, "N"), !str_detect(b, "N")) 

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
cor_tbl <- rbindlist(cor_tbl)

# calculating bootstrap values

# function to repeat 1000 times for bootstrapping:
sampled_cor <- function(markers_nested, sample_a, sample_b, marker_n) {
  
  markers_sampled <-  markers_nested %>% 
    slice_sample(n=marker_n, replace = TRUE) %>% 
    tidyr::unnest(cols = c(data))
  
  length_a  <-  markers_sampled %>% 
    filter(sample == sample_a) %>% 
    arrange(marker) %>% 
    pull(length)
  
  length_b  <-  markers_sampled %>% 
    filter(sample == sample_b) %>% 
    arrange(marker) %>% 
    pull(length)
  
  cor <- cor(length_a, length_b)
  
  list(a=sample_a, b=sample_b, cor=cor)

}

get_cor_bootstrap <- function(i, combos, markerlengths) {
  
  sample_a <- combos$a[i]
  sample_b <- combos$b[i]

  markerlengths_a  <-  markerlengths %>% 
    filter(sample == sample_a)
  
  markerlengths_b  <-  markerlengths %>% 
    filter(sample == sample_b)
  
  # find markers that are in both samples
  common_markers <- intersect(markerlengths_a$marker, markerlengths_b$marker)
  marker_n <- length(common_markers)
  
  markers_nested  <-  markerlengths_a %>% 
    filter(marker %in% common_markers) %>% 
    bind_rows(markerlengths_b %>% 
                  filter(marker %in% common_markers)) %>% 
    group_by(marker) %>% 
    nest() %>% 
    ungroup
  
  replicates_1000 <- replicate(1000, sampled_cor(markers_nested, sample_a, sample_b, marker_n), simplify = FALSE)
  bind_rows(replicates_1000)

}

bootstrap_cor <- lapply(1:nrow(combos), get_cor_bootstrap, combos, markerlengths)
bootstrap_cor_1000 <- bind_rows(bootstrap_cor)

ci <- bootstrap_cor_1000 %>%
  group_by(a, b) %>%
  summarise(quantile_0.025=quantile(cor, 0.025), quantile_0.975=quantile(cor, 0.975))

# changing names to published names
new_ids <- read_excel("data/naxerova_newnames.xlsx") %>%
  select(5:7) %>%
  rename(nax_id=1, pat_id=2, samp_id=3) %>% 
  mutate(
    nax_id = case_when(
      str_detect(nax_id, "WGS004") ~ str_replace(nax_id, "WGS004", "WGS00D"),
      str_detect(nax_id, "WGS002") ~ str_replace(nax_id, "WGS002", "WGS00B"),
      str_detect(nax_id, "WGS005") ~ str_replace(nax_id, "WGS005", "WGS00E"),
      TRUE ~ nax_id
    ),
    pat_id=str_replace(pat_id, " ", "_"), 
    new_id=paste0(pat_id,"_", samp_id)) %>% 
  select(nax_id, new_id)


left_join(cor_tbl, ci) %>%
  left_join(new_ids, by = c("a" = "nax_id")) %>%
  mutate(a = new_id) %>%
  select(-new_id) %>%
  left_join(new_ids, by = c("b" = "nax_id")) %>%
  mutate(b = new_id) %>%
  select(-new_id) %>%
  mutate(patient=str_extract(a, "Patient_[:digit:]+")) %>% 
  relocate(patient) %>%
  write_tsv("results/batch1_cor.tsv")
