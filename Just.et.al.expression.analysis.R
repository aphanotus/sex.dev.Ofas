# Analysis script
# for Just, Laslo, Lee, Yarnell, Zhang, and Angelini
# Distinct developmental mechanisms influence sexual dimorphisms in the milkweed bug Oncopeltus fasciatus

# Tissue- and isoform-specific gene expression analysis 

# # Clear the R workspace
# rm(list=ls())

# Load packages
library(tidyverse)
library(viridis)

# Importing the data
exp.adult <- read.csv("qpcr.tissue.teneral.adults.csv")
exp.juv <- read.csv("qpcr.tissue.juveniles.csv")

exp.adult$tissue <- sub('ventral thorax','legs',exp.adult$tissue)
exp.adult$tissue <- sub('genitalia','genital',exp.adult$tissue)

############
# Normalize expression to the reference gene
############

normalize.maintaining.zeros <- function(target, ref) {
  target <- as.numeric(target)
  ref <- as.numeric(ref)
  if (all(target[!is.na(target)]==0)) { return(target) }
  else {
    target[which(is.na(ref) | ref==0)] <- NA
    nonzeros <- which(target>0)
    n <- target - ref
    n[-nonzeros] <- 0
    n[which(is.na(target))] <- NA
    return(n)
  }
}
mean.omit.zero <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) { return(NA) } 
  mean(x[which(x!=0)], na.rm = TRUE)
}
min.omit.zero <- function(x) {
  x <- as.numeric(x)
  if (all(is.na(x))) { return(NA) } 
  min(x[which(x!=0)], na.rm = TRUE)
}


norm.adult <- exp.adult %>% 
  mutate(instar = "adult", day = "d1") %>% 
  select(sample, tissue, instar, day, sex)

{
  norm.adult$dsx1 <- normalize.maintaining.zeros(exp.adult$dsx1_3t4, exp.adult$ef1a_3t4)
  norm.adult$dsx2 <- normalize.maintaining.zeros(exp.adult$dsx2_1t2, exp.adult$ef1a_1t2)
  norm.adult$dsx2b <- normalize.maintaining.zeros(exp.adult$dsx2b_5t6, exp.adult$ef1a_5t6)
  norm.adult$dsx3 <- normalize.maintaining.zeros(exp.adult$dsx3_7t8, exp.adult$ef1a_7t8)
  norm.adult$dsx4 <- normalize.maintaining.zeros(exp.adult$dsx4_9t10, exp.adult$ef1a_9t10)
  norm.adult$dsx6 <- normalize.maintaining.zeros(exp.adult$dsx6_11t12, exp.adult$ef1a_11t12)
  norm.adult$dmrt99B1 <- normalize.maintaining.zeros(exp.adult$dmrt99B1_3t4, exp.adult$ef1a_3t4)
  norm.adult$dmrt99B2 <- normalize.maintaining.zeros(exp.adult$dmrt99B2_13t14, exp.adult$ef1a_13t14)
  norm.adult$dmrt99B3 <- normalize.maintaining.zeros(exp.adult$dmrt99B3_15t16, exp.adult$ef1a_15t16)
  norm.adult$dmrt99B5 <- normalize.maintaining.zeros(exp.adult$dmrt99B5_17t18, exp.adult$ef1a_17t18)
  norm.adult$dmrt93B1 <- normalize.maintaining.zeros(exp.adult$dmrt93B1_3t4, exp.adult$ef1a_3t4)
  norm.adult$dmrt93B3 <- normalize.maintaining.zeros(exp.adult$dmrt93B3_19t20, exp.adult$ef1a_19t20)
  norm.adult$dmrt93B4a <- normalize.maintaining.zeros(exp.adult$dmrt93B4a_21t22, exp.adult$ef1a_21t22)
  norm.adult$dmrt93B4b <- normalize.maintaining.zeros(exp.adult$dmrt93B4b_23t24, exp.adult$ef1a_23t24)
  norm.adult$dmrt93B5 <- normalize.maintaining.zeros(exp.adult$dmrt93B5_25t26, exp.adult$ef1a_25t26)
  
  # intersex and fruitless were assayed multiple times
  # The code below takes the average after normalization
  tmp.fru1 <- normalize.maintaining.zeros(exp.adult$fru_1t2, exp.adult$ef1a_1t2)
  tmp.fru2 <- normalize.maintaining.zeros(exp.adult$fru_5t6, exp.adult$ef1a_5t6)
  tmp.ix1 <- normalize.maintaining.zeros(exp.adult$ix_3t4, exp.adult$ef1a_3t4)
  tmp.ix2 <- normalize.maintaining.zeros(exp.adult$ix_7t8, exp.adult$ef1a_7t8)

  na.i <- which(is.na(tmp.fru1) & is.na(tmp.fru2))
  tmp.fru1[which(is.na(tmp.fru1))] <- 0
  tmp.fru2[which(is.na(tmp.fru2))] <- 0
  norm.adult$fru <- (tmp.fru1+tmp.fru2)/2
  norm.adult$fru[na.i] <- NA
  
  na.i <- which(is.na(tmp.ix1) & is.na(tmp.ix2))
  tmp.ix1[which(is.na(tmp.ix1))] <- 0
  tmp.ix2[which(is.na(tmp.ix2))] <- 0
  norm.adult$ix <- (tmp.ix1+tmp.ix2)/2
  norm.adult$ix[na.i] <- NA

  rm(tmp.ix1, tmp.ix2, tmp.fru1, tmp.fru2, na.i)
}

norm.juv <- exp.juv %>% 
  mutate(sample = sample_name) %>% 
  select(sample, tissue, instar, day, sex) 

{
  norm.juv$dsx1 <- normalize.maintaining.zeros(exp.juv$dsx1, exp.juv$ef1a)
  norm.juv$dsx6 <- normalize.maintaining.zeros(exp.juv$dsx6, exp.juv$ef1a)
  norm.juv$fru <- normalize.maintaining.zeros(exp.juv$fru, exp.juv$ef1a)
}

# Reorder the columns to match
norm.adult <- norm.adult %>% 
  select(
    sample, tissue, instar, day, sex, dsx1, dsx6, fru, 
    dmrt99B1, dmrt99B2, dmrt99B3, dmrt99B5, dmrt93B1, dmrt93B3, dmrt93B4a, dmrt93B4b, dmrt93B5, 
    dsx2, dsx2b, dsx3, dsx4, ix
  )

# Join the datasets
norm.q <- bind_rows(norm.juv, norm.adult)

############
# Standardize gene expression
############

# Find the global minimum "floor" value 
# as a fraction of the mean value, which will be defined as 1
x <- unlist(c(exp.juv[,8:10],exp.adult[,4:22]))
plot(x)
abline(v=length(unlist(exp.juv[,8:10]))+0.5)
abline(h=mean.omit.zero(x), lty = 2, col = "darkred")

relative.min <- min.omit.zero(x) / mean.omit.zero(x)

# Find the min and mean non-zero, normalized values
x <- unlist(norm.q[,-c(1:5)])
plot(x)

n.min <- min.omit.zero(x)
n.mean <- mean.omit.zero(x)

zeros <- norm.q
zeros[,1:5] <- FALSE
zeros[,-c(1:5)] <- (norm.q[,-c(1:5)] == 0)
zeros <- as.matrix(zeros)

std.q <- norm.q
std.q[,-c(1:5)] <- (norm.q[,-c(1:5)] - n.min*(1+relative.min)) / (n.mean - n.min*(1+relative.min))
std.q[zeros] <- 0
plot(unlist(c(std.q[,-c(1:5)])))

write_csv(std.q,'standardized.tissue.specific.expression.csv')

############
# Plots
############

colors.indicating.sex <- c(
  juvenile = "#808076",  # "concrete"
  female   = "#ff9933",  # "neon carrot" or "deep saffron"
  male     = "#666699"   # "scampi" or "mostly desaturated dark blue"
)
# scales::show_col(colors.indicating.sex)

q.long <- std.q %>% 
  pivot_longer(cols = 6:(dim(std.q)[2]), names_to = "target") %>% 
  filter(!(target %in% c("abdA","abdB"))) %>% 
  mutate(tissue = gsub("_"," ",tissue))

q.long$tissue <- sub('dorsal thorax','thorax',q.long$tissue)
q.long$tissue <- factor(
  sub(" ","\n",q.long$tissue), 
  levels = c(
    "sternites","viscera","posterior\nabdomen",
    'head','thorax','legs','ovaries','ovipositor','testes','genital\ncapsule'
  )
)

q.long$target <- factor(
  q.long$target, 
  levels = c(
    "dsx1", "dsx2", "dsx2b", "dsx3", "dsx4", "dsx6", "ix", "fru",
    "dmrt99B1", "dmrt99B2", "dmrt99B3", "dmrt99B5",
    "dmrt93B1", "dmrt93B3", "dmrt93B4a", "dmrt93B4b", "dmrt93B5"       
  )
)

expression.box.plots <- q.long %>% 
  filter(!grepl("dmrt",target)) %>% 
  filter(!(instar == "L5" & sex == " j")) %>%
  mutate(target = sub("dsx","dsx exon",target)) %>% 
  mutate(sex = sub(" j","juvenile",sex)) %>% 
  mutate(sex = sub("^f$","female",sex)) %>% 
  mutate(sex = sub("^m$","male",sex)) %>% 
  mutate(instar = sub("L"," instar ",instar)) %>% 
  mutate(day    = sub("d","day ",day)) %>% 
  mutate(stage = paste0(instar,"\n",day)) %>% 
  ggplot(aes(sex, value, color=sex)) +
  theme_bw() + 
  theme(
    legend.position="none",
    panel.grid.minor = element_blank()
  ) +
  facet_grid(target ~ stage*tissue, scales = 'free_x', space = 'free_x') +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.25, size=2, alpha = 0.75) +
  scale_color_manual(values = colors.indicating.sex) +
  labs(y = "relative expression")
# expression.box.plots

ggsave("plots/expression.box.plots.png", expression.box.plots, height = 12, width = 14, scale = 1)
ggsave("plots/expression.box.plots.pdf", expression.box.plots, height = 12, width = 14, scale = 1)

# For DMRTs
expression.dmrt.box.plots <- q.long %>% 
  filter(grepl("dmrt",target)) %>% 
  mutate(sex = sub(" j","juvenile",sex)) %>% 
  mutate(sex = sub("^f$","female",sex)) %>% 
  mutate(sex = sub("^m$","male",sex)) %>% 
  mutate(instar = sub("L"," instar ",instar)) %>% 
  mutate(day    = sub("d","day ",day)) %>% 
  mutate(stage = paste0(instar,"\n",day)) %>% 
  filter(!is.na(value)) %>% 
  ggplot(aes(sex, value, color=sex)) +
  theme_bw() + 
  theme(
    legend.position="none",
    panel.grid.minor = element_blank()
  ) +
  facet_grid(target ~ stage*tissue, scales = 'free_x', space = 'free_x') +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.25, size=2, alpha = 0.75) +
  scale_color_manual(values = colors.indicating.sex) +
  labs(y = "relative expression")

ggsave("plots/expression.dmrt.box.plots.png", expression.dmrt.box.plots, height = 9, width = 6.5, scale = 1)
ggsave("plots/expression.dmrt.box.plots.pdf", expression.dmrt.box.plots, height = 9, width = 6.5, scale = 1)


# Heat map
q.means <- q.long %>% 
  filter(!(instar == "L5" & sex == " j")) %>%
  mutate(instar = sub("L"," instar ",instar)) %>% 
  mutate(day    = sub("d","day ",day)) %>% 
  mutate(stage = paste0(instar,", ",day)) %>% 
  group_by(target, tissue, stage, sex) %>% 
  summarise(mean.exp = mean(value, na.rm=TRUE)) %>% 
  filter(!grepl('dmrt',target)) %>% 
  droplevels() 

{
  facet.labels <- c(` j` = "indeterminate", f = "female", m = "male")
  target.labels <- sub('dsx','dsx exon',rev(levels(q.means$target)))
  target.labels <- sub('ix','intersex',target.labels)
  target.labels <- sub('fru','fruitless',target.labels)
}

exp.heatmap <- ggplot(q.means, aes(x=tissue, y=target, z=mean.exp)) +
  theme_minimal() + 
  theme(
    strip.text.x = element_text(size=8),
    axis.title.x = element_blank(),
    axis.text.x  = element_text(size=8, face="bold", angle=50, vjust=1, hjust=0),
    axis.title.y = element_blank(),
    axis.text.y = element_text(size=8, face="italic"),
    axis.ticks = element_blank(),
    legend.title = element_text(size=8),
    panel.grid = element_blank()
  ) +
  facet_grid(.~stage*sex, scales="free_x", space="free_x", labeller=labeller(sex = facet.labels)) +
  geom_tile(aes(fill = mean.exp)) +
  scale_fill_viridis('relative\nexpression', option = 'plasma', direction = 1) +
  scale_x_discrete(position = "top") +
  scale_y_discrete(name='qPCR target\n',
                   limits = rev(levels(q.means$target)), 
                   labels = target.labels ) 

exp.heatmap

ggsave('plots/expression.heatmap.png', exp.heatmap, scale=1.25, width=9, height=3.125)
ggsave('plots/expression.heatmap.pdf', exp.heatmap, scale=1.25, width=9, height=3.125)

# Distribution of fourth instar gene expression

fourth.instar.expression.histograms <- q.long %>% 
  filter(instar == "L4" & !is.na(value)) %>% 
  mutate(target = sub('dsx','dsx exon',target)) %>% 
  group_by(target, tissue) %>%
  ggplot(aes(x=value)) +
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  facet_grid(target ~ tissue) +
  geom_histogram(binwidth = 0.25) +
  labs(x="relative expression")

fourth.instar.expression.histograms
  
ggsave('fourth.instar.expression.histograms.pdf', fourth.instar.expression.histograms, 
       scale=1, width=5, height=4.5)


############
# Stats
############

targets.of.interest <- c("dsx1","dsx2","dsx2b","dsx3","dsx4","dsx6","ix","fru")

q.stats <- q.long %>% 
  filter(instar != "L4", sex != " j", !is.na(sex)) %>% 
  filter(target %in% targets.of.interest, !is.na(value)) %>% 
  mutate(stage = paste0(instar,"-",day)) %>% 
  mutate(tissue = sub("\\\n","_",tissue)) %>% 
  mutate(tissue = sub("testes","gonad",tissue)) %>% 
  mutate(tissue = sub("ovaries","gonad",tissue)) %>% 
  mutate(tissue = sub("genital_capsule","terminalia",tissue)) %>% 
  mutate(tissue = sub("ovipositor","terminalia",tissue)) %>% 
  mutate(stage.tissue = paste0(stage,"-",tissue))

{
  sink(file = "expression.wrst.csv")
  cat("target,stage,tissue,f,m,W,p\n")
  unlist(lapply(targets.of.interest, function(target.i) {
    unlist(lapply(unique(q.stats$stage.tissue), function(stage.tissue.i) {
      dat.i <- q.stats %>% filter(target == target.i, stage.tissue == stage.tissue.i)
      if (dim(dat.i)[1]>4) {
        cat(paste0(target.i,",",dat.i$stage[1],",",dat.i$tissue[1],",",paste0(c(by(dat.i$sex,dat.i$sex,length)), collapse = ","),","))
        wrst.i <- wilcox.test(value ~ sex, data = dat.i)
        cat(wrst.i$statistic,",",wrst.i$p.value,"\n")
      }
    }))
  }))
  sink()
}

expression.wrst <- read.csv("expression.wrst.csv")

expression.wrst$padj <- expression.wrst %>% 
  group_by(target) %>% 
  summarise(padj = p.adjust(p, method = "fdr", n = n())) %>% 
  pull(padj)

write.csv(expression.wrst, "expression.wrst.csv", quote = FALSE, row.names = FALSE)


