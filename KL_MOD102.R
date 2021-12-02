koala_gene <- "AGCCTTAAATAACGACCTTC"

n_nucleo <- function(nucleo, sequence) {
  result <- gregexpr(nucleo, sequence)
  if (length(result[[1]]) == 1 && result[[1]] == -1) {
    result <- 0
  } else {
    result <- length(result[[1]])
  }
  result
}
n_nucleo("A", koala_gene)

# Q1
nucs <- c("A", "T", "C", "G")
nuc_freq <- sapply(nucs, FUN=n_nucleo, sequence = koala_gene)
nuc_freq
# A T C G 
# 7 5 6 2 


# Q2
gc_content <- function(sequence) {
  total <- nchar(koala_gene)
  gc <- n_nucleo("C", sequence) + n_nucleo("G", sequence)
  gc/total
}
gc_content("AAG")
# 0.05
gc_content("TGCCCATGGG")
# 0.35


# Q3
marsupial_gene <- "AGCCTTTCAACACGACCTTC"

n_mutations <- function(my_gene, reference_gene) {
  split_my_gene <- strsplit(my_gene, "")[[1]]
  split_ref_gene <- strsplit(reference_gene, "")[[1]]
  sum(split_my_gene != split_ref_gene)
}
n_mutations(koala_gene, marsupial_gene)
# 4


# Q4 Code modified from: https://r4ds.had.co.nz/functions.html#functions-are-for-humans-and-computers
f1 <- function(string, prefix) {
  substr(string, 1, nchar(prefix)) == prefix
}
f2 <- function(x) {
  if (length(x) <= 1) return(NULL)
  x[-length(x)]
}
f3 <- function(x) {
  rep("A", length.out = length(x))
}

# Function f1 checks if a string has a certain prefix. 
# Function f2 returns the input without the last element. 
# Function f3 returns a list of "A" that is the same length as a list x.
# To improve function f3, I would turn "A" into a general parameter that you can input into the function.
# Better function names would be:
#   f1 -> check_prefix 
#   f2 -> delete_last_element
#   f3 -> repeated_vector


# Q5
# I would revise the name for f2 to be rm_last_element to make it a bit shorter than "delete_last_element"




library(tidyverse)
library(palmerpenguins)
ggplot(drop_na(penguins, flipper_length_mm, bill_length_mm), 
       aes(flipper_length_mm, bill_length_mm, color = species)) +
  geom_point(aes(shape = species), size = 2, alpha = 0.75) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  labs(x = "Flipper length (mm)",
       y = "Bill length (mm)",
       color = "Penguin species",
       shape = "Penguin species") +
  scale_color_manual(values = c(Adelie = "#E4811E",
                                Chinstrap = "#B05CC5",
                                Gentoo = "#417175")) +
  theme_minimal() +
  theme(legend.box.background = element_rect(color = NA, fill = "white"),
        legend.justification = c(1, 0),
        legend.position = c(0.95, 0.05))

ith_bill <- function(i, bills) {
  sort(bills, decreasing = TRUE)[i]
}


# Q6
chinstrap <- penguins %>% filter(species=='Chinstrap')
gentoo <- penguins %>% filter(species=='Gentoo')
adelie <- penguins %>% filter(species=='Adelie')
longestBill_ind = 3

ith_bill(longestBill_ind, chinstrap$bill_length_mm)
# 54.2
ith_bill(longestBill_ind, gentoo$bill_length_mm)
# 55.1
ith_bill(longestBill_ind, adelie$bill_length_mm)
# 45.6




ith_bill_impure <- function(i) {
  bills <- penguins$bill_length_mm
  sort(bills, decreasing = TRUE)[i]
}

# Q7
# It's not a pure function because the output also depends on the 
# bill_length_mm column from the penguins data frame.
# Therefore, ith_bill_impure

ith_bill_impure(3)
# 55.9
penguins$bill_length_mm[100:300] <- NA # deleting some data
ith_bill_impure(3)
# 53.5 <- same input, but now the output is different!



# Q8
# sqrt() - pure 
# rnorm() - impure; output is stochastic
# lm() - pure
# read.csv() - impure; if you make changes to the input csv file, then the output will be different
# write.csv() - impure; creates/modifies external file


rm(penguins)
library(palmerpenguins)

# Q9
penguin_morpho <- select(penguins, bill_length_mm, flipper_length_mm, species)
penguin_morpho <- drop_na(penguin_morpho)
penguin_clust <- kmeans(select(penguin_morpho, -species), 
                        centers = n_distinct(penguins$species))
penguin_morpho$cluster <- penguin_clust$cluster
cluster_species <- penguin_morpho %>% 
  group_by(species) %>% 
  summarize(cluster = names(which.max(table(cluster)))) %>% 
  transmute(clustered_species = species, 
            cluster = as.numeric(cluster))
penguin_morpho <- left_join(penguin_morpho, cluster_species, by = "cluster")
n_misclass <- sum(penguin_morpho$species != penguin_morpho$clustered_species)
n_total <- nrow(penguin_morpho)
n_misclass / n_total
# 0.1608187



# Rewritten:
prepare_data <- function(penguins) {
  penguin_morpho <- penguins %>% select(bill_length_mm, flipper_length_mm, species) %>% drop_na()
  penguin_morpho
}

find_clusters <- function(penguin_morpho) {
  penguin_clust <- kmeans(select(penguin_morpho, -species), 
                          centers = n_distinct(penguin_morpho$species))
  penguin_morpho$cluster <- penguin_clust$cluster
  cluster_species = penguin_morpho %>% 
    group_by(species) %>% 
    summarize(cluster = names(which.max(table(cluster)))) %>% 
    transmute(clustered_species = species, 
              cluster = as.numeric(cluster))
  penguin_morpho <- left_join(penguin_morpho, cluster_species, by = "cluster")
  penguin_morpho
}

assess_clusters <- function(penguin_morpho) {
  n_misclass <- sum(penguin_morpho$species != penguin_morpho$clustered_species)
  n_total <- nrow(penguin_morpho)
  n_misclass / n_total
}

penguin_morpho <- prepare_data(penguins) %>% find_clusters()
misclass <- assess_clusters(penguin_morpho)
misclass
# 0.1608187 (same answer)





# body_condition() quantifies the body condition of a penguin as its residual
# (normalized) mass relative to its structural size (bill length, bill depth, 
# and flipper length). Use lm() to fit a regression and resid() to get the 
# residuals. 
body_condition <- function(mass, bill_length, bill_depth, flipper_length) {
  mass_normalized <- (mass - mean(mass)) / sd(mass)
  mod <- lm(mass_normalized ~ bill_length + bill_depth + flipper_length)
  resid(mod)
}

penguins_condition <- penguins %>% 
  drop_na() %>% 
  group_by(species) %>% 
  mutate(condition = body_condition(body_mass_g, bill_length_mm, bill_depth_mm, flipper_length_mm)) %>% 
  ungroup()

# condition_trend() calculates the slope of body condition over time. Use lm() 
# to fit a regression and coef() to get the slope.
condition_trend <- function(body_condition, year) {
  mod <- lm(body_condition ~ year)
  coef(mod)[[2]]
}

penguins_trend <- penguins_condition %>% 
  group_by(species) %>% 
  summarize(trend = condition_trend(condition, year)) %>% 
  mutate(trend_lbl = sprintf("slope=%0.2f", trend))

ggplot(penguins_condition, aes(factor(year), condition)) + 
  geom_violin() + 
  geom_smooth(aes(group = species), method = "lm", se = FALSE) + 
  geom_text(aes(label = trend_lbl), penguins_trend, x = 2, y = 0.05) +
  facet_grid(rows = vars(species)) +
  theme_minimal()


