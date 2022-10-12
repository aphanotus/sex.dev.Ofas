# Analysis script
# for Just, Laslo, Lee, Yarnell, Zhang, and Angelini
# Distinct developmental mechanisms influence sexual dimorphisms in the milkweed bug Oncopeltus fasciatus

# Analysis of phenotypes

# # Clear the R workspace
# rm(list=ls())

# Load packages
library(tidyverse)

# Define plot style elements
nice.theme <- theme_bw() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size=11, face="bold")
  )

colors.indicating.sex <- c(
  female = "#ff9933",  # "neon carrot" or "deep saffron"
  male   = "#666699"   # "scampi" or "mostly desaturated dark blue"
)

sigstar <- function (x) {symnum(x, cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "** ", "*  ", "   "), na = '   ') }

# Load data on genitalia length
gen <- read.csv("genitalia.length.data.csv")
gen <- gen %>% filter(stage != "L5")
gen$normalized.genitalia.length <- gen$genitalia.length / gen$femur.length

# Re-order the dsRNA treatment levels
gen$treatment <- factor(gen$treatment,
                        levels=unique(gen$treatment)[c(8,1,7,6,5,3,4,2)] )

# Find mean control values
average.genitalia <- c(
  female = filter(gen, treatment=='AmpR', sex=="female") %>% pull(normalized.genitalia.length) %>% mean(na.rm=TRUE),
  male   = filter(gen, treatment=='AmpR', sex=="male") %>% pull(normalized.genitalia.length) %>% mean(na.rm=TRUE)
)

# Generate the plot
genitalia.plot <- gen %>%
  ggplot(aes(x=treatment, y=normalized.genitalia.length)) +
  nice.theme +
  facet_grid(sex~.) +
  geom_hline(yintercept = average.genitalia["female"], size=0.65, color=colors.indicating.sex["female"]) +
  geom_hline(yintercept = average.genitalia["male"], size=0.65, color=colors.indicating.sex["male"]) +
  geom_violin(color='gray35', fill='gray90', size=0.5, alpha = 0.5) +
  geom_jitter(position = position_jitter(width = .15), alpha = 0.85 ) +
  scale_color_manual(values=c('black','gray50')) +
  xlab("dsRNA treatment") +
  ylab("normalized genitalia length")

genitalia.plot
ggsave("genitalia.plot.pdf", genitalia.plot,
       width = 7, height = 5.5, scale = 1)

# Statistical tests

# For controls, are males and females different?
with( filter(gen, treatment=='AmpR'),
      boxplot(normalized.genitalia.length ~ sex)
)
# Yes, they are non-ovarlapping
with( filter(gen, treatment=='AmpR'),
      wilcox.test(normalized.genitalia.length ~ sex)
)
# W = 2756, p-value < 2.2e-16

f <- gen %>% filter(sex == "female")
with(f, kruskal.test(normalized.genitalia.length, treatment))
# Kruskal-Wallis chi-squared = 113.41, df = 7, p-value < 2.2e-16
AmpR.f <- f %>% filter(treatment == "AmpR") %>% pull(normalized.genitalia.length)
ix.f <- f %>% filter(treatment == "ix") %>% pull(normalized.genitalia.length)
fru.f <- f %>% filter(treatment == "fru") %>% pull(normalized.genitalia.length)
dsx.f <- f %>% filter(treatment == "dsx") %>% pull(normalized.genitalia.length)
dmrt99B.f <- f %>% filter(treatment == "dmrt99B") %>% pull(normalized.genitalia.length)
dmrt93B.f <- f %>% filter(treatment == "dmrt93B") %>% pull(normalized.genitalia.length)
triple.f <- f %>% filter(treatment == "dmrtTriple") %>% pull(normalized.genitalia.length)
wilcox.test(AmpR.f, ix.f, alternative = "greater")      # W = 2093, p-value = 1.098e-12 ***
wilcox.test(AmpR.f, fru.f, alternative = "greater")     # W = 703, p-value = 3.225e-06 ***
wilcox.test(AmpR.f, dsx.f, alternative = "greater")     # W = 1835, p-value = 6.533e-05 ***
wilcox.test(AmpR.f, dmrt99B.f, alternative = "greater")    # W = 271, p-value = 0.9985
wilcox.test(AmpR.f, dmrt93B.f, alternative = "greater")    # W = 501, p-value = 0.8911
wilcox.test(AmpR.f, triple.f, alternative = "greater")  # W = 2268, p-value = 0.0008763 **
# FDR correction
x <- c(1.098e-12, 3.225e-06, 6.533e-05, 0.0008763)
sigstar(p.adjust(x, method = "fdr", n=6)) # *** *** *** **

m <- gen %>% filter(sex == "male")
with(m, kruskal.test(normalized.genitalia.length, treatment))
# Kruskal-Wallis chi-squared = 100.99, df = 7, p-value < 2.2e-16
AmpR.m <- m %>% filter(treatment == "AmpR") %>% pull(normalized.genitalia.length)
ix.m <- m %>% filter(treatment == "ix") %>% pull(normalized.genitalia.length)
fru.m <- m %>% filter(treatment == "fru") %>% pull(normalized.genitalia.length)
dsx.m <- m %>% filter(treatment == "dsx") %>% pull(normalized.genitalia.length)
dmrt99B.m <- m %>% filter(treatment == "dmrt99B") %>% pull(normalized.genitalia.length)
dmrt93B.m <- m %>% filter(treatment == "dmrt93B") %>% pull(normalized.genitalia.length)
triple.m <- m %>% filter(treatment == "dmrtTriple") %>% pull(normalized.genitalia.length)
wilcox.test(AmpR.m, ix.m, alternative = "greater")      # W = 1486, p-value = 7.93e-11 ***
wilcox.test(AmpR.m, fru.m, alternative = "greater")     # W = 695, p-value = 2.323e-06 ***
wilcox.test(AmpR.m, dsx.m, alternative = "greater")     # W = 1301, p-value = 3.432e-08 ***
wilcox.test(AmpR.m, dmrt99B.m, alternative = "greater")    # W = 499, p-value = 0.8076
wilcox.test(AmpR.m, dmrt93B.m, alternative = "greater")    # W = 347, p-value = 0.9487
wilcox.test(AmpR.m, triple.m, alternative = "greater")  # W = 1955, p-value = 3.336e-08
# FDR correction
x <- c(7.93e-11, 2.323e-06, 3.432e-08, 3.336e-08)
sigstar(p.adjust(x, method = "fdr", n=6)) # *** *** *** ***


# Load data on sternite curvature
curv <- read.csv("sternite.curvature.csv")

# Re-order the dsRNA treatment levels
curv$treatment <- factor(curv$treatment,
                         levels=unique(curv$treatment)[c(8,1,7,6,2,4,5,3)] )

# Find mean control values
average.curvature <- c(
  female = filter(curv, treatment=='AmpR', sex=="female") %>% pull(curvature) %>% mean(na.rm=TRUE),
  male   = filter(curv, treatment=='AmpR', sex=="male") %>% pull(curvature) %>% mean(na.rm=TRUE)
)

# Generate the plot
curvature.plot <- curv %>%
  ggplot(aes(x=treatment, y=curvature)) +
  nice.theme +
  facet_grid(sex~.) +
  geom_hline(yintercept = average.curvature["female"], size=0.65, color=colors.indicating.sex["female"]) +
  geom_hline(yintercept = average.curvature["male"], size=0.65, color=colors.indicating.sex["male"]) +
  geom_violin(color='gray35', fill='gray90', size=0.5, alpha = 0.5) +
  geom_jitter(position = position_jitter(width = .15), alpha = 0.85 ) +
  scale_color_manual(values=c('black','gray50')) +
  xlab("dsRNA treatment") +
  ylab("A4 sternite curvature")

curvature.plot
ggsave("curvature.plot.pdf", curvature.plot,
       width = 7, height = 5.5, scale = 1)

# Statistical tests

# For controls, are males and females different?
with( filter(curv, treatment=='AmpR'),
      boxplot(curvature ~ sex)
)
with( filter(curv, treatment=='AmpR'),
      wilcox.test(curvature ~ sex)
)
# W = 2133, p-value < 2.2e-16

f <- curv %>% filter(sex == "female")
with(f, kruskal.test(curvature, treatment))
# Kruskal-Wallis chi-squared = 43.157, df = 7, p-value = 3.111e-07
AmpR.f <- f %>% filter(treatment == "AmpR") %>% pull(curvature)
ix.f <- f %>% filter(treatment == "ix") %>% pull(curvature)
fru.f <- f %>% filter(treatment == "fru") %>% pull(curvature)
dsx.f <- f %>% filter(treatment == "dsx") %>% pull(curvature)
dmrt99B.f <- f %>% filter(treatment == "dmrt99B") %>% pull(curvature)
dmrt93B.f <- f %>% filter(treatment == "dmrt93B") %>% pull(curvature)
triple.f <- f %>% filter(treatment == "dmrtTriple") %>% pull(curvature)
wilcox.test(AmpR.f, ix.f, alternative = "greater")      # W = 537, p-value = 0.9326
wilcox.test(AmpR.f, fru.f, alternative = "greater")     # W = 639, p-value = 0.01072   *
wilcox.test(AmpR.f, dsx.f, alternative = "greater")     # W = 2692, p-value = 0.004728 *
wilcox.test(AmpR.f, dmrt99B.f, alternative = "greater")    # W = 474, p-value = 0.2512
wilcox.test(AmpR.f, dmrt93B.f, alternative = "greater")    # W = 546, p-value = 0.1588
wilcox.test(AmpR.f, triple.f, alternative = "greater")  # W = 1611, p-value = 0.04582 ns
# FDR correction
x <- c(0.01072, 0.004728, 0.04582)
sigstar(p.adjust(x, method = "fdr", n=6)) # *   *   ns

m <- curv %>% filter(sex == "male")
with(m, kruskal.test(curvature, treatment))
# Kruskal-Wallis chi-squared = 31.906, df = 7, p-value = 4.228e-05
AmpR.m <- m %>% filter(treatment == "AmpR") %>% pull(curvature)
ix.m <- m %>% filter(treatment == "ix") %>% pull(curvature)
fru.m <- m %>% filter(treatment == "fru") %>% pull(curvature)
dsx.m <- m %>% filter(treatment == "dsx") %>% pull(curvature)
dmrt99B.m <- m %>% filter(treatment == "dmrt99B") %>% pull(curvature)
dmrt93B.m <- m %>% filter(treatment == "dmrt93B") %>% pull(curvature)
triple.m <- m %>% filter(treatment == "dmrtTriple") %>% pull(curvature)
wilcox.test(AmpR.m, ix.m, alternative = "less")      # W = 267, p-value = 0.06832
wilcox.test(AmpR.m, fru.m, alternative = "less")     # W = 472, p-value = 0.8305
wilcox.test(AmpR.m, dsx.m, alternative = "less")     # W = 446, p-value = 1.117e-07 ***
wilcox.test(AmpR.m, dmrt99B.m, alternative = "less")    # W = 500, p-value = 0.3649
wilcox.test(AmpR.m, dmrt93B.m, alternative = "less")    # W = 362, p-value = 0.2503
wilcox.test(AmpR.m, triple.m, alternative = "less")  # W = 665, p-value = 0.004148 *
# FDR correction
x <- c(1.117e-07, 0.004148)
sigstar(p.adjust(x, method = "fdr", n=6)) # *** *


# Load data on abdominal pigmentation
mel <- read.csv("ventral.pigmentation.csv")
mel <- mel %>% filter(stage != "L5")

# Re-order the dsRNA treatment levels
mel$treatment <- factor(mel$treatment,
                        levels=unique(mel$treatment) )

# Find mean control values
average.melanism <- c(
  female = filter(mel, treatment=='AmpR', sex=="female") %>% pull(A54.melanic.ratio) %>% mean(na.rm=TRUE),
  male   = filter(mel, treatment=='AmpR', sex=="male") %>% pull(A54.melanic.ratio) %>% mean(na.rm=TRUE)
)

# Generate the plot
melanism.plot <- mel %>%
  ggplot(aes(x=treatment, y=A54.melanic.ratio)) +
  nice.theme +
  facet_grid(sex~.) +
  geom_hline(yintercept = average.melanism["female"], size=0.65, color=colors.indicating.sex["female"]) +
  geom_hline(yintercept = average.melanism["male"], size=0.65, color=colors.indicating.sex["male"]) +
  geom_violin(color='gray35', fill='gray90', size=0.5, alpha = 0.5) +
  geom_jitter(position = position_jitter(width = .15), alpha = 0.85 ) +
  scale_color_manual(values=c('black','gray50')) +
  xlab("dsRNA treatment") +
  ylab("ratio of abdominal melanism\n(melanic area of A5 vs. A4)")

melanism.plot
ggsave("genitalia.plot.pdf", genitalia.plot,
       width = 7, height = 5.5, scale = 1)

# Statistical tests

# For controls, are males and females different?
with( filter(mel, treatment=='AmpR'),
      boxplot(A54.melanic.ratio ~ sex)
)
# Yes, they are non-ovarlapping

f <- mel %>% filter(sex == "female")
with(f, kruskal.test(A54.melanic.ratio, treatment))
# Kruskal-Wallis chi-squared = 49.277, df = 7, p-value = 2.003e-08
AmpR.f <- f %>% filter(treatment == "AmpR") %>% pull(A54.melanic.ratio)
ix.f <- f %>% filter(treatment == "ix") %>% pull(A54.melanic.ratio)
fru.f <- f %>% filter(treatment == "fru") %>% pull(A54.melanic.ratio)
dsx.f <- f %>% filter(treatment == "dsx") %>% pull(A54.melanic.ratio)
dmrt99B.f <- f %>% filter(treatment == "dmrt99B") %>% pull(A54.melanic.ratio)
dmrt93B.f <- f %>% filter(treatment == "dmrt93B") %>% pull(A54.melanic.ratio)
triple.f <- f %>% filter(treatment == "dmrtTriple") %>% pull(A54.melanic.ratio)
wilcox.test(AmpR.f, ix.f, alternative = "less")      # W = 775, p-value = 0.8926
wilcox.test(AmpR.f, fru.f, alternative = "less")     # W = 318, p-value = 0.4093
wilcox.test(AmpR.f, dsx.f, alternative = "less")     # W = 575, p-value = 0.0002185 ***
wilcox.test(AmpR.f, dmrt99B.f, alternative = "less")    # W = 212, p-value = 0.1125
wilcox.test(AmpR.f, dmrt93B.f, alternative = "less")    # W = 373, p-value = 0.3967
wilcox.test(AmpR.f, triple.f, alternative = "less")  # W = 679, p-value = 7.828e-05 ***
# FDR correction
x <- c(0.0002185, 7.828e-05)
sigstar(p.adjust(x, method = "fdr", n=6)) # *** ***

m <- mel %>% filter(sex == "male")
with(m, kruskal.test(A54.melanic.ratio, treatment))
# Kruskal-Wallis chi-squared = 37.832, df = 7, p-value = 3.262e-06
AmpR.m <- m %>% filter(treatment == "AmpR") %>% pull(A54.melanic.ratio)
ix.m <- m %>% filter(treatment == "ix") %>% pull(A54.melanic.ratio)
fru.m <- m %>% filter(treatment == "fru") %>% pull(A54.melanic.ratio)
dsx.m <- m %>% filter(treatment == "dsx") %>% pull(A54.melanic.ratio)
dmrt99B.m <- m %>% filter(treatment == "dmrt99B") %>% pull(A54.melanic.ratio)
dmrt93B.m <- m %>% filter(treatment == "dmrt93B") %>% pull(A54.melanic.ratio)
triple.m <- m %>% filter(treatment == "dmrtTriple") %>% pull(A54.melanic.ratio)
wilcox.test(AmpR.m, ix.m, alternative = "greater")      # W = 384, p-value = 0.01732 *
wilcox.test(AmpR.m, fru.m, alternative = "greater")     # W = 214, p-value = 0.7549
wilcox.test(AmpR.m, dsx.m, alternative = "greater")     # W = 330, p-value = 7.617e-07 ***
wilcox.test(AmpR.m, dmrt99B.m, alternative = "greater")    # W = 419, p-value = 0.05999
wilcox.test(AmpR.m, dmrt93B.m, alternative = "greater")    # W = 320, p-value = 0.2133
wilcox.test(AmpR.m, triple.m, alternative = "greater")  # W = 989, p-value = 0.0009465 **
# FDR correction
x <- c(0.01732, 7.617e-07, 0.0009465)
sigstar(p.adjust(x, method = "fdr", n=6)) # *   *** **

# Among dsx RNAi specimens, is there still any sexual dimorphism?
dsx.results <- mel %>% filter(treatment == "dsx")
with(dsx.results, boxplot(A54.melanic.ratio ~ sex, ylab="ratio of abdominal mexlanism\n(melanic area of A5 vs. A4)"))
abline(h=average.melanism["female"], col=colors.indicating.sex["female"], lwd=3)
abline(h=average.melanism["male"], col=colors.indicating.sex["male"], lwd=3)
with(dsx.results, wilcox.test(A54.melanic.ratio ~ sex))
# W = 199, p-value = 0.2469


# Correlation between genitalia length and sternite curvature?
matching.ids <- which(curv$specimen.id %in% gen$specimen.id)
curv$normalized.genitalia.length <- NA
for (i in 1:length(matching.ids)) {
  j <- which(gen$specimen.id == curv$specimen.id[matching.ids[i]])
  if (length(j) == 1) {
    curv$normalized.genitalia.length[matching.ids[i]] <- gen$normalized.genitalia.length[j]
  }
}

curvature.vs.genitalia.plot <- curv %>%
  filter(!is.na(normalized.genitalia.length)) %>%
  ggplot(aes(x=normalized.genitalia.length, y=curvature, color=sex)) +
  nice.theme +
  theme(legend.justification=c(0.5,0.5), legend.position=c(5/6,1/6)) +
  facet_wrap(.~treatment) +
  geom_smooth(method=lm, fill = "gray85") +
  geom_point(alpha = 0.85) +
  scale_color_manual(values=colors.indicating.sex) +
  ylab("A4 sternite curvature") +
  xlab("normalized genitalia length")

curvature.vs.genitalia.plot
ggsave("curvature.vs.genitalia.plot.by.dsRNA.pdf", curvature.vs.genitalia.plot,
       width = 7, height = 7, scale = 1)


# Correlation between sternite curvature and abdominal melanism?
matching.ids <- which(mel$specimen.id %in% curv$specimen.id)
mel$curvature <- NA
for (i in 1:length(matching.ids)) {
  j <- which(curv$specimen.id == mel$specimen.id[matching.ids[i]])
  if (length(j) == 1) {
    mel$curvature[matching.ids[i]] <- curv$curvature[j]
  }
}

# Rough plots
with(mel, plot(curvature, A54.melanic.ratio, col = as.factor(sex),
                    xlab = "A4 sternite curvature", ylab = "ratio of abdominal melanism"))
x <- which(mel$sex=="female")
abline(with(mel[x,], lm(A54.melanic.ratio ~ curvature)), col = 1)
with(mel[x,], cor.test(curvature, A54.melanic.ratio, method = "spearman"))
# S = 2324224, p-value = 0.01464, rho = -0.1612643
x <- which(mel$sex=="male")
abline(with(mel[x,], lm(A54.melanic.ratio ~ curvature)), col = 2)
with(mel[x,], cor.test(curvature, A54.melanic.ratio, method = "spearman"))
# S = 729266, p-value = 0.01379, rho = -0.1981011

# Just for no-dsRNA controls
x <- which(mel$sex=="female" & mel$treatment=="none")
with(mel[x,], cor.test(curvature, A54.melanic.ratio, method = "spearman"))
# S = 2184, p-value = 0.7194, rho = -0.07905138
x <- which(mel$sex=="male" & mel$treatment=="none")
with(mel[x,], cor.test(curvature, A54.melanic.ratio, method = "spearman"))
# S = 1970, p-value = 0.2195, rho = -0.2792208

# AmpR dsRNA controls
x <- which(mel$sex=="female" & mel$treatment=="AmpR")
with(mel[x,], cor.test(curvature, A54.melanic.ratio, method = "spearman"))
# S = 6578, p-value = 0.07376, rho = -0.3262097
x <- which(mel$sex=="male" & mel$treatment=="AmpR")
with(mel[x,], cor.test(curvature, A54.melanic.ratio, method = "spearman"))
# S = 6990, p-value = 0.02303, rho = -0.4092742

melanism.vs.curvature.plot <- mel %>%
  filter(!is.na(curvature)) %>%
  ggplot(aes(x=curvature, y=A54.melanic.ratio, color=sex)) +
  nice.theme +
  theme(legend.justification=c(0.5,0.5), legend.position=c(5/6,1/6)) +
  facet_wrap(.~treatment) +
  geom_smooth(method=lm, fill = "gray85") +
  geom_point(alpha = 0.85) +
  scale_color_manual(values=colors.indicating.sex) +
  xlab("A4 sternite curvature") +
  ylab("ratio of abdominal melanism\n(melanic area of A5 vs. A4)")

melanism.vs.curvature.plot
ggsave("melanism.vs.curvature.plot.pdf", melanism.vs.curvature.plot,
       width = 7, height = 7, scale = 1)
