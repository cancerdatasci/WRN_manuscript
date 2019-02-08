library(magrittr)
library(plyr)
library(tidyverse)
library(lme4)

df <- read_csv('~/CPDS/WRN_manuscript/data/KM12_in_vivo_combined.csv') %>%
  reshape2::melt(id.vars = c('Tumor', 'Dox', 'Day'), value.name = 'Volume', variable.name = 'Mouse') %>% 
  mutate(Mouse = paste0(Tumor, '.', Mouse, '.', Dox)) 

used_samples <- df %>% group_by(Mouse) %>% summarise(n = sum(!is.na(Volume))) %>% filter(n >= 5) %>% .[['Mouse']]
df %<>% filter(Mouse %in% used_samples)

all_slopes <- ddply(df, .(Mouse, Tumor, Dox), function(aa) {
  lm(Volume ~ Day, aa) %>% 
    broom::tidy()
}) 
group_slopes <- all_slopes %>% 
  filter(term == 'Day') %>% 
  group_by(Tumor, Dox) %>% 
  summarise(mean_slope = mean(estimate))

#####compare DOX vs No-DOX for shWRN
m1 <- lmer(Volume ~ Day + (0 + Day|Mouse), df %>% filter(Tumor == 'shWRN1'), REML = FALSE)
m2 <- lmer(Volume ~ Day + Day:Dox + (0 + Day|Mouse), df %>% filter(Tumor == 'shWRN1'), REML = FALSE)
anova(m1, m2)

###compare DOX vs No-DOX for shWRN-C911
# mod <- lmer(Volume ~ Day + Day:Dox + (1|Mouse), df %>% filter(Tumor == 'shWRN1-C911'), REML = TRUE)
# summary(mod)
m1 <- lmer(Volume ~ Day + (0 + Day|Mouse), df %>% filter(Tumor == 'shWRN1-C911'), REML = FALSE)
m2 <- lmer(Volume ~ Day + Day:Dox + (0 + Day|Mouse), df %>% filter(Tumor == 'shWRN1-C911'), REML = FALSE)
anova(m1, m2)

###Compare shWRN vs C911 with Dox
# mod <- lmer(Volume ~ Day + Day:Tumor + (1|Mouse), df %>% filter(Dox == 'Y'), REML = TRUE)
# summary(mod)
m1 <- lmer(Volume ~ Day + (0 + Day|Mouse), df %>% filter(Dox == 'Y'), REML = FALSE)
m2 <- lmer(Volume ~ Day + Day:Tumor + (0 + Day|Mouse), df %>% filter(Dox == 'Y'), REML = FALSE)
anova(m1, m2)


#compare DOX vs No-DOX for shWRN
mod <- lm(Volume ~ Day + Day:Dox, df %>% filter(Tumor == 'shWRN1'))
summary(mod)

