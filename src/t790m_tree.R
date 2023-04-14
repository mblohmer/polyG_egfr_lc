# loading needed libraries
library(tidyverse) 
library(ggtree)
library(tidytree)
library(ape)
library(phangorn)
library(Rphylip)
library(data.table)
library(pheatmap)
library(ggpubr)

# Finding the polyG mean lengths ------------------------------------

dir <-  list.files("../data/batch2/",
                        pattern= "T790M_0.1_20220607_R")

subject <-  str_split(dir, "_") %>% purrr::map(1) %>% unlist   

  #select normal
  normal_sample <- c("T790M_III_1_N2")
  normal <- paste0("^", normal_sample)
  
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
  markers <- get_marker_lengths(subject,marker_dir)

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
  
  root_sample <- select_normal(subject)

  markerlengths <- markers %>% 
    group_by(marker) %>% 
    filter(any(str_detect(sample, normal))) %>% 
    mutate(length = length - length[str_detect(sample, normal)][1]) %>% 
    as.data.table
  
# Continue with the rest of the script ------------------------------------

markerlengths <- markerlengths %>% 
  mutate(sample = str_sub(sample, end = -3))

# find number of markers
n_markers <- length(unique(markerlengths$marker))

# remove any markers that didn't amplify in all samples
filtered_markerlengths <- markerlengths %>% 
  group_by(subject, marker) %>% 
  mutate(subject_markers = length(unique(sample))) %>% 
  group_by(subject) %>% 
  mutate(sample_markers=max(subject_markers)) %>% 
  filter(sample_markers==subject_markers) %>% 
  select(-c(subject_markers, sample_markers)) %>% 
  add_count(sample, name="markers") %>% 
  ungroup() %>% 
  as.data.table()


# Plotting Mean Length Heatmap --------------------------------------------

# rotating mean length table for heatmap plotting
filtered_matrix <- filtered_markerlengths %>%
  mutate(sample=str_remove_all(sample, "T790M_")) %>% 
  filter(sample!="III_1_N") %>% 
  ungroup() %>% 
  select(-subject, -markers)  %>% 
  pivot_wider(names_from = marker, values_from=length)  %>% 
  column_to_rownames("sample")  %>% 
  as.matrix()

# setting minmax to center the heatmap color at 0
minmax <- max(abs(filtered_matrix), na.rm = TRUE)

# creating heatmap
hm <- pheatmap(filtered_matrix,
                         clustering_distance_rows="manhattan",
               clustering_distance_cols="manhattan",
                         fontsize = 6,
                         color = colorRampPalette(c("seagreen", "white", "purple3"))(42),
                         na_col = "yellow",
                         breaks=seq(-minmax, minmax, length.out=43), 
                        legend_breaks = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6,  minmax), 
                        legend_labels = c(-0.6, -0.4, -0.2, 0, 0.2, 0.4, 0.6, "Mean\nlength\n"),
                        treeheight_col=15,
                        treeheight_row=15)
# saving heatmap
ggsave(paste0("../plots/T790M_heatmap.pdf"), hm,
       width = 5, height = 4
)


# Mean length plots -------------------------------------------------------

filtered_markerlengths %>%
  mutate(sample=str_remove_all(sample, "T790M_")) %>% 
  filter(sample=="707"|sample=="710") %>% 
  select(-subject, -markers) %>% 
  pivot_wider(names_from = sample, values_from = length) %>%
  ggplot(aes(`710`, `707`)) +
  geom_smooth(method = "lm", size = 3, alpha = 0.2) +
  geom_point(size = 6) +
  stat_cor(aes(label = ..r.label..)) +
  labs(x = "Mean length sample 710", y = "Mean length sample 707") +
  theme_classic() 
ggsave("../plots/T790M_low_cor.pdf", height=4, width=4)

filtered_markerlengths %>%
  mutate(sample=str_remove_all(sample, "T790M_")) %>% 
  filter(sample=="716"|sample=="710") %>% 
  select(-subject, -markers) %>% 
  pivot_wider(names_from = sample, values_from = length) %>%
  ggplot(aes(`710`, `716`)) +
  geom_smooth(method = "lm", size = 3, alpha = 0.2) +
  geom_point(size = 6) +
  stat_cor(aes(label = ..r.label..)) +
  labs(x = "Mean length sample 710", y = "Mean length sample 716") +
  theme_classic() 
ggsave("../plots/T790M_high_cor.pdf", height=4, width=4)

# Calculating L1 ---------------------------------------------------

# find samples
samples <- filtered_markerlengths$sample %>% unique

combos_wide <- combn(samples, m= 2) %>% as.data.frame()

## make combo table longer
combos <- data.frame(a=as.character(combos_wide[1,]), b=as.character(combos_wide[2,]))

# setting up function to calculate L1
get_l1_for_combination <- function(i, combos, markerlengths) {

  sample_a <- combos$a[i]
  sample_b <- combos$b[i]
  
  marker_diff  <-  markerlengths %>% 
    filter(sample %in% c(sample_a, sample_b))  %>% 
    group_by(marker) %>% 
    summarize(length_diff=abs(length[sample==sample_a]-length[sample==sample_b]), 
              .groups="drop")
  
  
  n_markers  <-  n_distinct(marker_diff$marker)

  l1  <- marker_diff %>% 
    summarize(L1=sum(length_diff,  na.rm = TRUE)/n_markers) %>% 
    pull(L1) 
  
  list(a=sample_a, b=sample_b, l1=l1, marker=n_markers)
}

# caluclting l1 for all sample combinations

l1_table <- lapply(1:nrow(combos), get_l1_for_combination, combos, filtered_markerlengths)
l1_table <- rbindlist(l1_table)


# Tree for T790M with colored tips ----------------------------------------
subject <- "T790M"

# colors for T790M
colors  <- c("darkred", "red", "pink", "darkgreen", "darkgoldenrod1",
             "darkgreen", "darkgreen", "darkgreen", "darkgreen", "darkgreen","darkgreen", "darkgreen")
names(colors)  <- c("T790M_707", "T790M_709", "T790M_708", "T790M_710", "T790M_713",
                    "T790M_714A", "T790M_714B", "T790M_714C", "T790M_715A", "T790M_715B","T790M_716", "T790M_717")


# Tree with only 2004 samples  --------------------------------------------



#convert variance to cell divisions:
filtered_l1_table <- l1_table %>% 
  select(-marker) %>% 
  # selecting only the 2004 samples
  filter(str_detect(a, "T790M_III_1_N2|707|708|709|713|710"), str_detect(b, "T790M_III_1_N2|707|708|709|713|710"))

# adding variance of one sample with itself (0)
zero_tbl <- tibble("a"=unique(filtered_l1_table$a), "b"=unique(filtered_l1_table$a), "l1"=0) %>% 
  bind_rows(tibble("a"=unique(filtered_l1_table$b), "b"=unique(filtered_l1_table$b), "l1"=0)) %>% 
  distinct

dist_mat_int <- filtered_l1_table %>% 
  bind_rows(zero_tbl) %>% 
  pivot_wider(names_from=a, values_from=l1)

# reordering the matrix to be in the correct order for the tree
dist_mat <- dist_mat_int[match(colnames(dist_mat_int), dist_mat_int$b) %>% na.omit, ] %>% 
  column_to_rownames("b")

# building an nj tree
tree <- nj(as.dist(dist_mat))

# root tree

root_sample <-  "T790M_III_1_N2"
tree <- phytools::reroot(tree, node.number=which(tree$tip.label==root_sample))

tree_tbl <- as_tibble(tree) 

# plotting the tree

p <- ggtree(tree, root.position=0) +
  geom_tiplab(aes(color=label), size=5, show.legend = FALSE) +
  scale_colour_manual(values=colors) +
  theme_tree2() +
  xlim(0,0.16) +
  geom_strip('T790M_708', 'T790M_710', barsize=2, color='darkred', 
             label="2004 \nresection", offset.text=0.002, offset = 0.02, extend=0.2) +
  xlab("Manhattan Distance from zygote") 

suppressMessages(p <- p + labs(title="T790M poly-G Tree", x="Manhattan Distance from zygote")) 
p
ggsave(paste0("../plots/Fig1e_2004_samples.pdf"),p, width=11,height=5)


# Tree with 710 clade only ------------------------------------------------


#convert variance to cell divisions:
filtered_l1_table <- l1_table %>% 
  select(-marker) %>% 
  # selecting only the 2004 samples
  filter(!str_detect(a, "T790M_III_1_N$|707|708|709|713"), !str_detect(b, "T790M_III_1_N$|707|708|709|713"))

# adding variance of one sample with itself (0)
zero_tbl <- tibble("a"=unique(filtered_l1_table$a), "b"=unique(filtered_l1_table$a), "l1"=0) %>% 
  bind_rows(tibble("a"=unique(filtered_l1_table$b), "b"=unique(filtered_l1_table$b), "l1"=0)) %>% 
  distinct

dist_mat_int <- filtered_l1_table %>% 
  bind_rows(zero_tbl) %>% 
  pivot_wider(names_from=a, values_from=l1)

# reordering the matrix to be in the correct order for the tree
dist_mat <- dist_mat_int[match(colnames(dist_mat_int), dist_mat_int$b) %>% na.omit, ] %>% 
  column_to_rownames("b")

# building an nj tree
tree <- nj(as.dist(dist_mat))

# root tree

root_sample <-  "T790M_III_1_N2"
tree <- phytools::reroot(tree, node.number=which(tree$tip.label==root_sample))

tree_tbl <- as_tibble(tree) 

# plotting the tree
detach("package:ggpubr", unload=TRUE)
p <- ggtree(tree %>% rotate(16), root.position=0) +
  geom_tiplab(aes(label=label), size=5, show.legend = FALSE) +
  scale_colour_manual(values=colors) +
  theme_tree2() +
  xlim(0,0.16) +
  xlab("Manhattan Distance from zygote") +
  geom_strip('T790M_716', 'T790M_714C', barsize=2, color='darkgreen', 
             label="2014 \nresection", offset.text=0.002, offset = 0.018, extend=0.2) +
  geom_strip('T790M_710', 
             'T790M_710', barsize=2, color='darkred', 
             label="2004 \nresection", offset.text=0.002, offset = 0.018, extend=0.2)

suppressMessages(p <- p + labs(title="T790M poly-G Tree", x="Manhattan Distance from zygote")) 
p
ggsave(paste0("../plots/Fig1e_710_clade.pdf"),p, width=11,height=5)

