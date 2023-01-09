#CC model and vaccines 


#setwd('/work/Exam-stuff')

#install.packages("pacman")
pacman::p_load(R2jags, parallel, polspline, tidyverse, glue, patchwork, bayestestR)


########### INFO ABOUT THE GAME ##########
groupSize <- 4
ntrials <- 10
pi <- 1.6 # multiplication factor in game
ntokens <- 20
vals <- seq(0,ntokens,1) #possible values to contribute - from 0 to 20 tokens


########## LOADING DATA AND MERGING ###########
df <- read.csv("Exam-stuff/Data/pgg_vacc_data.csv", header = T)


########## DATA PREPARATION ########## 
ngroups <- length(unique(df$groupid))

## group average contribution on each trial 
Gga_np <- df %>% 
  group_by(groupid, period) %>% 
  summarise(mean_send = mean(senderscontribution)) 

## group average contribution on each trial for each subject
# this fixes problem with some groups not having 4 subject on each trial 
Ggas_np <- Gga_np %>% 
  slice(rep(row_number(), 4)) %>% 
  arrange(groupid, period)

## Nation for every group
Nation_df <- df %>% 
  group_by(groupid) %>% 
  summarise(nation = mean(nation))

## calculating winnings per group with a function

calculate_winnings_df <- function(df, pi, Nation_df){
  # sum of contributions for each group on each trial times pi
  win_gs <- df %>% 
    group_by(period, groupid) %>% 
    summarise(sum_con_trial = sum(senderscontribution)) %>% 
    mutate(win_trial = sum_con_trial*pi)
  
  # sum of winnings for each group
  winnings_m <- win_gs %>% 
    group_by(groupid) %>% 
    summarise(win = sum(win_trial))
  
  # merging with nation to later find winnings at nation-level
  winnings <- merge(winnings_m, Nation_df)
  
  return(winnings)
  
}

#Group winnings over all trials per group
winnings <- calculate_winnings_df(df, pi, Nation_df)


### making arrays for JAGS
Ggas <- array(data = Ggas_np$mean_send, 
              dim = c(groupSize, ntrials, ngroups))
Nation <- array(data = unlist(Nation_df$nation))

########## CC MODEL with vaccination rate ###########
nnations <- length(unique(df$nation))
ngroups <- length(unique(df$groupid))

Ga <- Ggas

df <- arrange(df, groupid, period)
c <- array(data = df$senderscontribution, dim = c(groupSize, ntrials, ngroups))

data <- list("groupSize", "ngroups", "ntrials", "nnations","c","Ga","X","Nation","invSigma") 
params <- c("beta0_alpha","betaX_alpha","beta0_rho","betaX_rho","beta0_omega","betaX_omega") 


#Run JAGS

#standardise percent vaccinated people
df$std_perc_vacc <- (df$percent_vacc - mean(df$percent_vacc)) / sd(df$percent_vacc)

vacc_df <- df %>% 
  group_by(groupid) %>% 
  summarise(vacc = mean(std_perc_vacc))

vacc <- vacc_df$vacc 
X = array(vacc)

invSigma <- solve(t(X)%*%X) # required for JZS priors

CC.samples <- jags.parallel(data, inits=NULL, params,
                               model.file ="Exam-stuff/CC_corr1.txt",
                               n.chains=32, n.iter=15000, n.burnin=5000, n.thin=1, n.cluster=4)




# --------------------- plotting ---------------------------- #

cc_samples <- data.frame(CC.samples$BUGSoutput$sims.list$betaX_alpha,
                            CC.samples$BUGSoutput$sims.list$betaX_omega,
                            CC.samples$BUGSoutput$sims.list$betaX_rho)
names(cc_samples) <- c("Initial belief", "Belief learning rate", 
                          "Cooperation preference")


cc_samples <- cc_samples_output

a1 <-  cc_samples %>% 
  ggplot()+
  geom_density(aes(x = `Initial belief`))+
  geom_point(aes(x = mean(`Initial belief`), y=0), 
             col ="black")+
  geom_segment(aes(x = ci(`Initial belief`, ci = 0.95)$CI_low, y = 0, 
                   xend = ci(`Initial belief`, ci = 0.95)$CI_high, yend = 0), 
               col = 'black')+
  geom_vline(xintercept = 0.0, linetype = 'dotted', color = 'black', size = .4)+
  theme_bw()+
  scale_x_continuous(limits = c(-1, 1))+
  scale_y_continuous(limits = c(0,5))


a1


a2 <-  cc_samples %>% 
  ggplot()+
  geom_density(aes(x = `Belief learning rate`))+
  geom_point(aes(x = mean(`Belief learning rate`), y=0), 
             col ="green")+
  geom_segment(aes(x = ci(`Belief learning rate`, ci = 0.95)$CI_low, y = 0, 
                   xend = ci(`Belief learning rate`, ci = 0.95)$CI_high, yend = 0), 
               col = 'green')+
  geom_vline(xintercept = 0.0, linetype = 'dotted', color = 'black', size = .4)+
  theme_bw()+
  scale_x_continuous(limits = c(-1, 1.5))+
  scale_y_continuous(limits = c(0,5))

a2


a3 <-  cc_samples %>% 
  ggplot()+
  geom_density(aes(x = `Cooperation preference`))+
  geom_point(aes(x = mean(`Cooperation preference`), y=0), 
             col ="black")+
  geom_segment(aes(x = ci(`Cooperation preference`, ci = 0.95)$CI_low, y = 0, 
                   xend = ci(`Cooperation preference`, ci = 0.95)$CI_high, yend = 0), 
               col = 'black')+
  geom_vline(xintercept = 0.0, linetype = 'dotted', color = 'black', size = .4)+
  theme_bw()+
  scale_x_continuous(limits = c(-1, 1))+
  scale_y_continuous(limits = c(0,5))


a3


a_all <- a1 + a2 +a3
cc_plot <- a_all + plot_annotation(
  title = 'Posterior distributions for Conditional Cooperation parameters')

cc_plot


all_plots <- w1 /cc_plot
all_plots <- all_plots + plot_annotation(tag_levels = 'A')

all_plots

ggsave("Exam-stuff/cc_post_dist.png", cc_plot)

ggsave("Exam-stuff/all plots_post_dist.png", all_plots)



