library(aod)
library(R.matlab)
library(pwr)
library(rsq)
library(lme4)

setwd('/Users/david/OneDrive/Documents/Yale/NTB_lab/batters_problem/vma_continuous')

### DEFINE PLOT FUNCTION ###
plot_bars <- function(dat_set){
  # dat_set <- subset(dat_set, nocue == 0)
  for (i_block in c(1,2,3,4,5)){
    dat_block <- subset(dat_set, block == (i_block+1))
    dat_L <- subset(dat_block, hpt == 0)
    dat_H <- subset(dat_block, hpt == 1)
    
    barplot(c(100+mean(dat_L$metric[dat_L$predictive == 1]),
              100+mean(dat_L$metric[dat_L$predictive == 0])), ylim=c(-10, 10))
    
    barplot(c(100+mean(dat_H$metric[dat_H$predictive == 1]),
              100+mean(dat_H$metric[dat_H$predictive == 0])), ylim=c(-10, 10))

  }
} 

plot_timecourse <- function(dat_cue, dat_nocue){
  # dat_cue <- subset(dat_set, nocue == 0)
  # dat_nocue <- subset(dat_set, nocue == 1)
  for (i_block in c(1,2,3,4,5)){
    block_cue_L_P <- subset(dat_cue, block == (i_block+1) & hpt == 0 & predictive == 1)
    block_cue_L_N <- subset(dat_cue, block == (i_block+1) & hpt == 0 & predictive == 0)
    
    block_cue_H_P <- subset(dat_cue, block == (i_block+1) & hpt == 1 & predictive == 1)
    block_cue_H_N <- subset(dat_cue, block == (i_block+1) & hpt == 1 & predictive == 0)
    
    # block_nocue_L_P <- subset(dat_nocue, block == i_block & hpt == 0 & predictive == 1)
    # block_nocue_L_N <- subset(dat_nocue, block == i_block & hpt == 0 & predictive == 0)
    # 
    # block_nocue_H_P <- subset(dat_nocue, block == i_block & hpt == 1 & predictive == 1)
    # block_nocue_H_N <- subset(dat_nocue, block == i_block & hpt == 1 & predictive == 0)
    
    if (i_block == 1){
      plot(i_block, mean(block_cue_H_P$metric), pch = 0, col = "blue", xlim = c(0, 6), ylim = c(-100, -90))
    }else{
      par(new = TRUE)
      plot(i_block, mean(block_cue_H_P$metric), pch = 0, col = "blue", xlim = c(0, 6), ylim = c(-100, -90))
    }
    
    par(new = TRUE)
    plot(i_block, mean(block_cue_H_N$metric), pch = 0, col = "red", xlim = c(0, 6), ylim = c(-100, -90))
    par(new = TRUE)
    plot(i_block, mean(block_cue_L_P$metric), pch = 9, col = "blue", xlim = c(0, 6), ylim = c(-100, -90))
    par(new = TRUE)
    plot(i_block, mean(block_cue_L_N$metric), pch = 9, col = "red", xlim = c(0, 6), ylim = c(-100, -90))
  }
  
  # Now for no cue cases:
  block_nocue_L_1 <- subset(dat_nocue, block == 1 & hpt == 0)
  block_nocue_H_1 <- subset(dat_nocue, block == 1 & hpt == 1)
  
  block_nocue_L_n <- subset(dat_nocue, block > 1 & hpt == 0)
  block_nocue_H_n <- subset(dat_nocue, block > 1 & hpt == 1)
  
  par(new = TRUE)
  plot(0, mean(block_nocue_L_1$metric), pch = 9, col = "black", xlim = c(0, 6), ylim = c(-100, -90))
  par(new = TRUE)
  plot(0, mean(block_nocue_H_1$metric), pch = 0, col = "black", xlim = c(0, 6), ylim = c(-100, -90))
  par(new = TRUE)
  plot(6, mean(block_nocue_L_n$metric), pch = 9, col = "black", xlim = c(0, 6), ylim = c(-100, -90))
  par(new = TRUE)
  plot(6, mean(block_nocue_H_n$metric), pch =0, col = "black", xlim = c(0, 6), ylim = c(-100, -90))
}

data_1 <- read.csv('H_metric_all_v1.txt', header = FALSE, na.strings = 'NaN')
data_1 <- na.omit(data_1)
nsub1 <- max(unique(data_1$V8)) #V8 is subject number.

data_2 <- read.csv('H_metric_all_v3.txt', header = FALSE, na.strings = 'NaN')
data_2 <- na.omit(data_2)
data_2$V8 <- data_2$V8 + nsub1
nsub1_2 <- max(unique(data_2$V8))

data_3 <- read.csv('H_metric_all_v4.txt', header = FALSE, na.strings = 'NaN')
data_3 <- na.omit(data_3)
data_3$V8 <- data_3$V8 + nsub1_2

data_comb <- rbind(data_1, rbind(data_2, data_3))

head(data_comb)
# V1=correct or not, V2=trial type, V3=PT, V4=subject

dat_ <- data.frame(
  metric = data_comb$V1,
  predictive = as.factor(data_comb$V2),
  pt = as.numeric(data_comb$V3),
  hpt = as.factor(data_comb$V4),
  nocue = as.factor(data_comb$V5),
  trial = data_comb$V6,
  block = data_comb$V7,
  subject = data_comb$V8)

dat_ <- na.omit(dat_)
dat0 <- subset(dat_, nocue == 1)
dat <- subset(dat_, nocue == 0)
# dat <- dat_

# gmod1 <- lmer(metric ~ hpt + predictive + block | subject, data = dat, REML = FALSE)
# gmod1_ <- lmer(metric ~ hpt + predictive | subject, data = dat, REML = FALSE)
gmod1 <- lm(metric ~ hpt + predictive + block, family = gaussian, data = dat)
gmod_all <- lm(metric ~ hpt*predictive*block, family = gaussian, data = dat)

summary(gmod1)
summary(gmod_all)

#### break down by H vs L-PT

# HPT:
dat_LPT <- subset(dat, hpt == 0)

gm_lpt <- lm(metric ~  predictive + block, family = gaussian, data = dat_LPT)
gm_all_lpt <- lm(metric ~ predictive*block, family = gaussian, data = dat_LPT)

summary(gm_lpt)
summary(gm_all_lpt)

dat_HPT <- subset(dat, hpt == 1)

gm_hpt <- lm(metric ~  predictive + block, family = gaussian, data = dat_HPT)
gm_all_hpt <- lm(metric ~ predictive*block, family = gaussian, data = dat_HPT)

summary(gm_hpt)
summary(gm_all_hpt)

plot_bars(dat)
plot_timecourse(dat, dat0)
# plot_timecourse(dat0)

# calculate power of GLM test:
m1 <- lm(metric ~ block, data = dat)
m2 <- lm(metric ~ predictive, data = dat)
m3 <- lm(metric ~ hpt, data = dat)


pwr.f2.test(u = 3, v = 80 - 3 - 1, f2 = .02/(1 - .02), sig.level = .05)
pwr.f2.test(u = 3, f2 = .02/(1 - .02), sig.level = .05, power = 0.8)

##      Multiple regression power calculation 
## 
##               u = 2 (degrees of freedom = number of independent variables)
##               v = 49.88971 (degrees of freedom = n - u - 1)
##              f2 = 0.4285714 (effect size (r^2/(1-r^2)))
##       sig.level = 0.001 (usually .05)
##           power = 0.8 (desired power)
