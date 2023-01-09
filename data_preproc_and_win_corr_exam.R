seed_id = 1982
set.seed(seed_id)

#install.packages("pacman")
pacman::p_load(R2jags, parallel, polspline, ggplot2, glue, tidyverse)

# defining a function for calculating the maximum of the posterior density (not exactly the same as the mode)
MPD <- function(x) {
  density(x)$x[which(density(x)$y==max(density(x)$y))]
}


##### Preprocessing #####

raw <- read.csv("216377/Module5/data/HerrmannThoeniGaechterDATA.csv", skip = 3) # Public goods game data
vaccData <- read.csv("Exam-stuff/Data/clean_data.csv")


### removing punishment condition - not relevant
raw <- raw %>% 
  filter(p == "N-experiment")

### setting up cities with their corresponding nation id + getting their data 
city = c("Melbourne", "Minsk", "Chengdu", "Copenhagen", "Bonn", "Athens", 
         "Seoul", "Samara", "Zurich", "St. Gallen", "Istanbul", "Nottingham", 
         "Boston", "Dnipropetrovs'k", "Muscat", "Riyadh" )
nation = c(1,2,3,4,5,6,7,8,9,9,10,11,12,13, 14, 15)

# doing this separately so we avoid the repeating nation (#9) in these lists
# bc then they become super useful later (hint: covariates at nation level)
id = seq(1,15,1)
#trust = c(.40, .42, .55, .67, .38, .24, .27, .24, .37, .16, .29, .27, .36)

### first round of merge
data <- data.frame(id) %>% 
  rename(nation = id)
cities = data.frame(city, nation)
data = merge(data, cities)


### second round
# merge automatically removes countries without data
# if you want to retain add all = T after on = 'city'
df <- merge(raw, data , on = "city") %>% 
  # select every third row
  filter(row_number() %% 3 == 1)


#adding vaccination data 
vacci <- vaccData %>% 
  select(percent_vacc, location)


#Creating nation col that matches with the numbers, matching the cities above ^
vacci$nation <- c(1, 2, 3, 4, 11, 5, 6, 14, 8, 15, 7, 9, 10, 13, 12)


#Merging vaccine data with PGG data 

df <- merge(df, vacci, on = "nation")

#write_csv(df, "Exam-stuff/Data/pgg_vacc_data.csv")


######### You can start from here: #############

df <- read.csv("Exam-stuff/Data/pgg_vacc_data.csv", header = T)



#### setting up variables ####

groupSize <- 4
ntrials <- 10
pi <- 1.6 # how much the group project was multiplied by
ntokens <- 20
vals <- seq(0,ntokens,1) #vals = values
#vals <- seq(1,21,1) #possible values to contribute - from 0 to 20 tokens



## group average contribution on each trial 
ngroups = length(unique(df$groupid))

Gga_np <- df %>% #ga = group average - np = no punishment
  group_by(groupid, period) %>% 
  summarise(mean_send = mean(senderscontribution))

# make it an array 
Gga <- array(data = Gga_np$mean_send, dim = c(ntrials, ngroups))



## looking at problematic groups with weird data from collection. 
df_gga <- merge(Gga_np, df, all = TRUE)
count_probs <- df_gga %>%
  group_by(groupid) %>%
  count(subjectid) %>%
  filter(n != 10) %>%
  mutate(diff = abs(n - 10))

problem_groups <- unique(count_probs$groupid)
length(problem_groups) #51 problematic groups. 
sum(count_probs$diff) #530 ???? 

## group average contribution on each trial for each subject
# this fixes problem with some groups not having 4 subject on each trial 
Ggas_np <- Gga_np %>% 
  slice(rep(row_number(), 4)) %>% 
  arrange(groupid, period)



## Nation for every group - giving each group a column for which nation they are from. 
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
Nation <- array(data = unlist(Nation_df$nation))
# also standardizing winnings
Y = array((winnings$win - mean(winnings$win)) / sd(winnings$win))



########### JAGS - WIN CORRELATION #######

#Number of groups and nations
ngroups = length(unique(df$groupid))
nnations = 15


#Things for JAGS
data <- list("ngroups", "Y", "nnations","X","Nation","invSigma") #data put into jags
params <- c("beta0","betaX") #parameters to track in jags


#### Run jags correlation model

#standardise percent vaccinated people
df$std_perc_vacc <- (df$percent_vacc - mean(df$percent_vacc)) / sd(df$percent_vacc)

vacc_df <- df %>% 
  group_by(groupid) %>% 
  summarise(vacc = mean(std_perc_vacc))

vacc <- vacc_df$vacc 
X = array(vacc)

invSigma <- solve(t(X)%*%X) # required for JZS priors

Y = array((winnings$win - mean(winnings$win)) / sd(winnings$win)) #standardising winnings


#Run jags model 
win.samples <- jags.parallel(data, inits=NULL, params,
                                model.file ="Exam-stuff/win_corr.txt",
                                n.chains=3, n.iter=15000, n.burnin=5000, n.thin=1, n.cluster=4)



 cor.test(vacc, winnings$win, method=c("pearson"))

# ----------------------Plotting ---------------------- #

samples_vacc <- data.frame(win.samples$BUGSoutput$sims.list$betaX)
names(samples_vacc) <- "Vacc"
summary_vacc <- data.frame(win.samples$BUGSoutput$summary)

w1 <- samples_vacc %>% 
  ggplot(aes(x = Vacc))+
  geom_density()+
  geom_point(aes(x=summary_vacc[2,1], y=0), 
             col ="black")+
  geom_segment(aes(x = summary_vacc[2,3], y = 0, 
                   xend = summary_vacc[2,7], yend = 0), 
               col = 'black')+
  geom_vline(xintercept = 0.0, linetype = 'dotted', color = 'black', size = .4)+
  theme_bw()+
  scale_x_continuous()+
  scale_y_continuous()+
  labs(x = "% of people vaccinated, standardized",
       title = 'Correlation between vaccination rates and group winnings')


w1

ggsave("win_vacc_cor.png", w1)



y<-rpois(10000,2)
plot(ecdf(X))

# Correlation tests 

cor.test(vacc, winnings$win, method=c("spearman")) #p = 0.02
cor.test(vacc, winnings$win, method=c("kendal")) #p = .03
cor.test(vacc, winnings$win, method=c("pearson")) #p = .04


