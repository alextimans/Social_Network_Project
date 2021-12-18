#### Analysis of paper: Schöne et al. 2021 ####
# Alexander Timans, MSc Statistics, Nov-Dec 2021
# for GESS Social Media Networks course, ETH Zurich

#### Package imports + wd ####
setwd("~/Desktop/GESS/analysis_env") # fit as needed

library(readr); library(lattice); library(tidyverse); library(plyr)
library(dplyr); library(lme4); library(lmerTest); library(MuMIn);
library(pbkrtest); library(r2glmm); library(performance); library(scales);
library(bestNormalize); library(ordinal); library(glmmTMB); library(vader);
library(VGAM); library(psych); library(effects); library(nlme);



#### Functions ####
transform_cols <- function(data) { 
  # obtain missing transformations of data for reproduction of analysis
  
  data_tf <- data %>% mutate(
    retweets = round(exp(retweets_log) - 1),
    retweets_exp1 = round(exp(retweets_log - 1)), # incorrect formula used by author
    retweets_reciprocal = 1 - reciprocallink(retweets + 1),
    retweets_reciprocal_exp1 = 1 - reciprocallink(retweets_exp1 + 1), # incorrect
    followers_log = log(followers + 1),
    sen_pos_og = sen_pos + 1,
    sen_neg_og = -(sen_neg + 1)
    )
  
  data_tf <- data_tf %>% relocate(retweets, retweets_reciprocal, # re-arrange cols
                                  retweets_exp1, retweets_reciprocal_exp1,
                                  .before = retweets_log) %>% 
    relocate(followers_log, .after = followers) %>%
    relocate(sen_pos_og, sen_neg_og, .after = sen_neg)

  return(data_tf)
}

# Political affiliation assignment
# Nf=#ppl followed, Nc=#Conserv followed, Nl=#Liberals followed
# eps1 = diff ratio Nc & Nl; eps2: diff ratio (Nc + Nl)/Nf
polit_aff <- function(Nf, Nc, Nl, eps1=1.2, eps2=0.2){
  user <- ""
  if(Nc > Nl){
    r = Nc/Nl
    if(r > eps1) user="Con" else user="Neu"
  } else if(Nl > Nc){
    r = Nl/Nc
    if(r > eps1) user="Lib" else user="Neu"
  } else if(Nc == Nl) user="Neu" 
  
  r2 = (Nc + Nl)/Nf
  if(r2 <= eps2) user="Neu"
  return(user)
}



#### Data import + transform ####
d.ferguson <- read_csv("original_data/ferguson_anonymise_cleaned.csv",
                       col_types = cols(...1 = col_skip()))
d.ferguson <- transform_cols(d.ferguson)
#df_f <- d.ferguson

d.hillary <- read_csv("original_data/hillary_anonymise_cleaned.csv",
                      col_types = cols(...1 = col_skip()))
d.hillary <- transform_cols(d.hillary)
#df_h <- d.hillary

# with binary affiliation (on import encoded as char)
d.trump <- read_csv("original_data/trump_anonymise_cleaned.csv",
                    col_types = cols(...1 = col_skip()))
d.trump <- transform_cols(d.trump) # get extra cols
#df_t <- d.trump

d.samesex <- read_csv("original_data/same_sex_anonymise_cleaned.csv",
                      col_types = cols(...1 = col_skip()))
d.samesex <- transform_cols(d.samesex) # get extra cols
#df_g <- d.samesex

# text analysis; trump data has only 4 cols (missing user info + sen pos/neg)
d.trump.text <- read_csv("original_data/trump_tweets_for_textanalsis.csv",
                         col_types = cols(...1 = col_skip()))
d.samesex.text <- read_csv("original_data/same_sex_tweets_for_textanalsis.csv",
                           col_types = cols(...1 = col_skip(), X = col_skip(),
                                            user_id = col_character()))
# column valence = sen_pos + sen_neg (see d.samesex.text) 



#### Data exploration ####
# NA's per column
sapply(d.trump, function(y) sum(length(which(is.na(y)))))

# Distribution of retweets/skewness
hist(d.trump$retweets)
hist(d.trump$retweets[d.trump$retweets <= 5])
hist(d.trump$retweets_log)
hist(d.trump$retweets_reciprocal) #freq=F
lines(density(d.trump$retweets_reciprocal))

sum(d.trump$retweets == 0)/nrow(d.trump)
sum(d.trump$retweets <= 5)/nrow(d.trump)
summary(d.trump$retweets)
sort(unique(d.trump$retweets))

hist(d.trump$retweets_reciprocal[d.trump$retweets_reciprocal <= 0.2])
sum(d.trump$retweets_reciprocal <= 0.2)/nrow(d.trump)
sort(unique(d.trump$retweets_reciprocal))[1:20]

skew(d.trump$retweets)
skew(d.trump$retweets_log)
skew(d.trump$retweets_reciprocal)

# Binary affiliation subsampling
unique(d.trump$binary_affiliation)
length(which(d.trump$binary_affiliation == 'Mid'))
nrow(d.trump[(d.trump$retweets >= 1 & d.trump$sen_neg >= 2), ]) #11932, not 30830 ?
count(d.trump$binary_affiliation); count(d.samesex$binary_affiliation)

# Distribution of followers
hist(d.trump$followers)
hist(d.trump$followers_log)
hist(scale(d.trump$followers_log))

boxplot(d.trump$followers, horizontal = T)
summary(d.trump$followers)
sum(d.trump$followers <= 800)/nrow(d.trump)
quantile(d.trump$followers, probs=c(0.8, 0.9, 0.95, 0.99, 1))
d.trump.filt <- d.trump[d.trump$followers <= 800,]
boxplot(d.trump.filt$followers, horizontal = T)



#### Data plots ####
# various exploratory plots
boxplot(retweets ~ sen_neg, d.trump)
boxplot(retweets_reciprocal ~ sen_neg, d.trump)
ggplot(d.trump, aes(x = sen_neg, y = retweets, group = sen_pos, col = sen_pos)) + 
  geom_point() + stat_summary(fun = mean, geom = "line")

d.trump.filt <- d.trump[(d.trump$retweets <= 5),] #& (d.trump$sen_pos != 4),]
boxplot(retweets ~ sen_neg, d.trump.filt)

with(d.trump, interaction.plot(x.factor = sen_neg,
                               trace.factor = sen_pos,
                               response = retweets))
with(d.trump, interaction.plot(x.factor = sen_neg, 
                               trace.factor = sen_pos,
                               response = retweets_reciprocal))
stripchart(retweets ~ sen_neg, vertical=TRUE, pch=1, data=d.trump)

# Effects plot (code as run by authors in analysis_manuscript.Rmd called MLM graphs)
eff <- as.data.frame(Effect(c("sen_pos","sen_neg"), mod.trump))
eff <- eff %>% filter(sen_pos == 0 | sen_neg == 0)
eff$intensity <- (eff$sen_pos + eff$sen_neg)
eff <- rbind(eff, eff[1,]) # add first row again
eff$valence[1:5] <- "sen_pos"; eff$valence[6:10] <- "sen_neg"
View(eff)

plot <- ggplot(eff, aes(x=intensity, y=fit, colour=valence)) + 
  theme(plot.title = element_text(hjust=0.5)) + 
  geom_line(size=1) + 
  labs(x="Emotional Intensity", y="Number of Retweets 1-(reciprocal+1)", 
       colour="") + 
  scale_color_manual(labels=c("Negative Language", "Positive Language"), 
                     values=c("red", "green3")) + 
  ggtitle("Election win") + 
  geom_ribbon(alpha=0.1, aes(ymin=lower, ymax=upper, colour=valence), 
              linetype=0,) + ylim(0.12, 0.22)
plot



#### Model analysis ####

# COMMENTS
# Working with correct retweet and reciprocal values
# Main focus on Study 1 (Hillary/Trump data)
# Findings should mostly also apply to Study 2
# Fitting with lmer(...) calls lmerTest::lmer



#### Original models as fitted by authors ####
# Trump
mod.trump <- lmer(retweets_reciprocal ~ sen_pos + sen_neg 
                 + sen_pos:sen_neg + scale(followers_log) + (1|user_id), 
                 d.trump)
summary(mod.trump)
r.squaredGLMM(mod.trump)

# Hillary
mod.hillary <- lmer(retweets_reciprocal ~ sen_pos + sen_neg 
                   + sen_pos:sen_neg + scale(followers_log) + (1|user_id), 
                   d.hillary)
summary(mod.hillary)
r.squaredGLMM(mod.hillary)



#### Models with retweets on original scale
mod.trump.retweet <- lmer(retweets ~ sen_pos + sen_neg 
                          + sen_pos:sen_neg + scale(followers_log) + (1|user_id), 
                          d.trump)
summary(mod.trump.retweet)
r.squaredGLMM(mod.trump.retweet)

mod.hillary.retweet <- lmer(retweets ~ sen_pos + sen_neg 
                            + sen_pos:sen_neg + scale(followers_log) + (1|user_id), 
                            d.hillary)
summary(mod.hillary.retweet)
r.squaredGLMM(mod.hillary.retweet)



#### Test SentiStrength score transformation ####
# Models by authors: Transformation of sen_pos, sen_neg to [0, 4]
# Comparing to original sen_pos in [1, 5] & sen_neg in [-5, -1]
mod.trump.og.scale <- lmer(retweets_reciprocal ~ sen_pos_og + sen_neg_og 
                  + sen_pos_og:sen_neg_og + scale(followers_log) + (1|user_id), 
                  d.trump)
summary(mod.trump.og.scale)
r.squaredGLMM(mod.trump.og.scale)

mod.hillary.og.scale <- lmer(retweets_reciprocal ~ sen_pos_og + sen_neg_og 
              + sen_pos_og:sen_neg_og + scale(followers_log) + (1|user_id), 
              d.hillary)
summary(mod.hillary.og.scale)
r.squaredGLMM(mod.hillary.og.scale)



#### Residual analysis ####
# Tukey-Anscombe plot
plot(mod.trump, main="Tukey-Anscombe")
# plot(resid(mod.trump) ~ fitted(mod.trump), main="Tukey-Anscombe")
# abline(h=0, lty=2)

# Scale-location plot
xyplot(sqrt(abs(resid(mod.trump))) ~ fitted(mod.trump), 
       type=c("p", "smooth", "g"), col.line="black")

# QQ plot
qqnorm(resid(mod.trump), main="Residuals")
qqline(resid(mod.trump), lty=2)
qqnorm(ranef(mod.trump)$user_id[, 1], main = "Random effects of user_id")
qqline(ranef(mod.trump)$user_id[, 1], lty=2)
# qqmath(mod.trump)

# Histograms vs. normal distribution
ran <- ranef(mod.trump); res <- resid(mod.trump)
#sd correspond to estimated std dev. from mod.trump for respective effect
hist(ran$user_id[,1], freq=F, breaks=60, main='Histogram random effect user_id vs. normal distr.', 
     xlab='user_id random intercept')
curve(dnorm(x, mean=0, sd=0.1506), col=2, lwd=2, add=TRUE)
hist(res, freq=F, breaks=60, main='Histogram residuals vs. normal distr.', 
     xlab='residuals')
curve(dnorm(x, mean=0, sd=0.2111), col=2, lwd=2, add=TRUE)

# observe data skewed to right and piled up too much in the center
skew(ran$user_id[,1]); skew(res) # positive, so skewed right
kurtosi(ran$user_id[,1]); kurtosi(res) # light tails

# Residuals vs predictors 
xyplot(resid(mod.trump) ~ sen_neg, data=d.trump, jitter.x=TRUE,
       abline=0, type=c("p", "a", "g"))
xyplot(resid(mod.trump) ~ sen_pos, data=d.trump, jitter.x=TRUE,
       abline=0, type=c("p", "a", "g"))



#### Model tests ####
# Term deletions and model comparisons
m4 <- lmer(retweets_reciprocal ~ sen_pos + sen_neg
           + scale(followers_log) + (1|user_id), d.trump)
m3 <- lmer(retweets_reciprocal ~ sen_neg + scale(followers_log) + (1|user_id), d.trump)
m2 <- lmer(retweets_reciprocal ~ sen_pos + scale(followers_log) + (1|user_id), d.trump)
m1 <- lmer(retweets_reciprocal ~ scale(followers_log) + (1|user_id), d.trump)
m0 <- lmer(retweets_reciprocal ~ 1 + (1|user_id), d.trump)

drop1(m4, test="LRT") # suggests that sen_neg cannot be dropped

anova(m0, m1, m2, m3, m4, mod.trump, type="3", ddf="Kenward-Roger") 
# suggests full model, m4, m2 do not add info
# best model by AIC: m3; best model by BIC: m1

MuMIn::model.sel(m0, m1, m2, m3, m4, mod.trump) #best model by AICc: m1

performance::compare_performance(m0, m1, m2, m3, m4, mod.trump) # long run time
# best model by AIC: m1; by BIC: m1 
# R^2 score is the same across models except m0 - why? because all include followers

# step-wise procedure
dat.step <- na.omit(d.trump[, -16]) # remove single NA disregarding col 16
mod.step <- lmer(retweets_reciprocal ~ sen_pos + sen_neg + sen_pos:sen_neg 
               + scale(followers_log) + (1|user_id), dat.step)
step(mod.step)
# suggests to drop interaction and sen_pos, but not sen_neg and not random eff
ranova(mod.step) # suggests not to drop random eff

# Additional marginal/pseudo R^2 score
performance::r2(mod.trump) # identical results to r.squaredGLMM



#### Impact of followers as covariate ####
# As seen in model comparisons, followers contains largest explanatory power by far
# Huge bump in model selection criterias when including it, after that barely
summary(mod.trump) # original model for comparison
r.squaredGLMM(mod.trump) # original model for comparison
anova(mod.trump) # followers obtains largest sum of squares part

mod.trump.nofollow <- lmer(retweets_reciprocal ~ sen_pos + sen_neg 
                  + sen_pos:sen_neg + (1|user_id), d.trump)
summary(mod.trump.nofollow)
r.squaredGLMM(mod.trump.nofollow) # much smaller than original model


#### Encoding of sen_pos & sen_neg as factors ####
lvls <- c("0", "1", "2", "3", "4")
# simple factors
# for ordered set ordered=TRUE
d.trump.fact <- d.trump
d.trump.fact$sen_neg <- factor(d.trump.fact$sen_neg, ordered=FALSE, levels=lvls)
d.trump.fact$sen_pos <- factor(d.trump.fact$sen_pos, ordered=FALSE, levels=lvls)
str(d.trump.fact$sen_neg); levels(d.trump.fact$sen_neg)

# Trump
mod.trump.fact <- lmer(retweets_reciprocal ~ sen_pos + sen_neg 
                  + sen_pos:sen_neg + scale(followers_log) + (1|user_id), 
                  d.trump.fact)
summary(mod.trump.fact)
r.squaredGLMM(mod.trump.fact)
anova(mod.trump.fact)

# Hillary
d.hillary.fact <- d.hillary
d.hillary.fact$sen_neg <- factor(d.hillary.fact$sen_neg, ordered=FALSE, levels=lvls)
d.hillary.fact$sen_pos <- factor(d.hillary.fact$sen_pos, ordered=FALSE, levels=lvls)
str(d.hillary.fact$sen_neg); levels(d.hillary.fact$sen_neg)

mod.hillary.fact <- lmer(retweets_reciprocal ~ sen_pos + sen_neg 
                         + sen_pos:sen_neg + scale(followers_log) + (1|user_id), 
                         d.hillary.fact)
summary(mod.hillary.fact)
r.squaredGLMM(mod.hillary.fact)
# summary returns polynomial contrasts for ordered factors
# capture trends in the data -> non-linear trends?



#### Different models ####
# Testing different models with sen_neg, sen_pos 
# as (simple) factors where suitable, i.e. data=d.trump.fact

# Non-parametric regression, not really adequate for inference statements
mod.trump.loess <- loess(retweets ~ sen_pos + sen_neg, data=d.trump, degree=2) # runtime ~30m


# Linear Mixed Model with polynomials up to degree 4 for main effects
mod.trump.poly <- lmer(retweets_reciprocal ~ poly(sen_pos, 4) + poly(sen_neg, 4)
                       + sen_pos:sen_neg + scale(followers_log) + (1|user_id), 
                       d.trump)
summary(mod.trump.poly)
r.squaredGLMM(mod.trump.poly)
anova(mod.trump.poly)


# Linear Mixed Model with "best" normalization transform for response
(best_norm <- bestNormalize(d.trump$retweets, warn=TRUE, quiet=FALSE))
retweets_norm <- best_norm$x.t
mod.trump.norm <- lmer(retweets_norm ~ sen_pos + sen_neg 
                  + sen_pos:sen_neg + scale(followers_log) + (1|user_id), 
                  d.trump.fact)
summary(mod.trump.norm)
r.squaredGLMM(mod.trump.norm)
# remove outliers based on TA residual plot
out <- which(resid(mod.trump.norm) >= 10) | which(fitted(mod.trump.norm) >= 15)
out <- as.integer(names(out)) # 102 outliers
d.trump.fact.out <- d.trump.fact[-out, ]; retweets_norm_out <- retweets_norm[-out]
# rerun without outliers
mod.trump.norm.out <- lmer(retweets_norm_out ~ sen_pos + sen_neg 
                       + sen_pos:sen_neg + scale(followers_log) + (1|user_id), 
                       d.trump.fact.out)
summary(mod.trump.norm.out)
# residual analysis looks better, although still strong heteroscedasticity


# Generalized Linear Mixed Model with zero inflated negative binomial
mod.trump.zi <- glmmTMB(retweets ~ sen_pos + sen_neg 
                        + sen_pos:sen_neg + scale(followers_log) + (1|user_id), 
                        data=d.trump.fact,
                        family=nbinom2,
                        ziformula=~1)
summary(mod.trump.zi)
performance::r2(mod.trump.zi) # runtime few mins
# remove outliers based on TA residual plot
out <- which(resid(mod.trump.zi) >= 500) | which(fitted(mod.trump.zi) >= 2000)
out <- as.integer(names(out)) # 173 outliers
d.trump.fact.out <- d.trump.fact[-out, ]
# rerun without outliers
mod.trump.zi.out <- glmmTMB(retweets ~ sen_pos + sen_neg 
                            + sen_pos:sen_neg + scale(followers_log) + (1|user_id), 
                            data=d.trump.fact.out,
                            family=nbinom2,
                            ziformula=~1)
summary(mod.trump.zi.out)
performance::r2(mod.trump.zi.out) # runtime few mins
# residual analysis looks better (plots are readable) without outliers


# Non-linear mixed model - currently not runable due to syntax issues
# for "best" normalization approach with removed outliers
out <- which(retweets_norm > 10) # ~400 outliers
retweets_norm <- retweets_norm[-c(out, 278188)]
hist(retweets_norm, freq=FALSE, xlim=c(-3, 3), main="Histogram retweets_norm vs. std. normal")
curve(dnorm(x, mean=0, sd=1), col=2, lwd=2, add=TRUE)

d.trump.fact.nlme <- d.trump.fact[-out, ]
d.trump.fact.nlme <- na.omit(d.trump.fact.nlme[, -16]) # remove single NA
mod.trump.nlme <- nlme(retweets_norm ~ SSasymp(sen_pos, sen_neg, sen_pos:sen_neg,
                                              scale(followers_log), user_id), 
                      data=d.trump.fact.nlme,
                      fixed=sen_pos + sen_neg + sen_pos:sen_neg + scale(followers_log) ~ 1,
                      random=user_id ~ 1)
summary(mod.trump.nlme)
r.squaredGLMM(mod.trump.nlme)



#### Vader sentiment scores comparison ####
d.text <- na.omit(d.trump.text) # remove 1 NA
samp.int <- sample.int(nrow(d.text), 10000) # sample 10000 values
# applying vader sentiment analysis scoring
vader.score <- vader_df(d.text$tweet_body[samp.int]) # contains 42 NA not sure why
vader.score.lemma <- vader_df(d.text$tweet_body_lemmatized[samp.int])
cor(vader.score$compound, vader.score.lemma$compound)

valence.score <- ceiling(vader.score$compound * 4) # rescale to [0, 4] and round to int
# since vader scores are normalized to [0, 1] we are optimistic about importance
valence.score.lemma <- ceiling(vader.score.lemma$compound * 4)
valence.senti <- d.text$valence[samp.int]
cor(valence.score, valence.senti)
cor(valence.score.lemma, valence.senti) # 0.549, quite low

# table with both sentiment scores from SentiStrength and VADER for given text samples
text.sentiments <- data.frame(tweet=d.text$tweet_body[samp.int],
                              tweet_lemma=d.text$tweet_body_lemmatized[samp.int],
                              vader=valence.score.lemma,
                              sentistrength=valence.senti)
# since the data containing tweet information and the data containing tweet texts
# is separated, we cannot refit a model on retweet predictions using VADER scores
# and compare this to SentiStrength, which would be interesting



#### Some presentation figures/tables ####
d.samp <- d.trump[sample(nrow(d.trump), 20),
                  c("tweet_id", "time", "retweets", "user_id", "followers", "sen_pos_og", "sen_neg_og")]
d.samp.text <- d.trump.text[sample(nrow(d.trump.text), 10), 1:2]
barchart(d.trump$binary_affiliation, horizontal=FALSE, main="Political affiliation, Trump tweets subset")



#### References ####
# (excl. papers, see presentation for full list) 
# https://data.library.virginia.edu/understanding-ordered-factors-in-a-linear-model/
# http://bbolker.github.io/mixedmodels-misc/glmmFAQ.html#model-diagnostics
# https://stat.ethz.ch/~meier/teaching/anova/random-and-mixed-effects-models.html
# https://doi.org/10.1111/2041-210X.13434
# https://ethz.ch/content/dam/ethz/special-interest/math/statistics/sfs/Education/Advanced%20Studies%20in%20Applied%20Statistics/course-material-1921/Regression/MixedModels_Lab.pdf
# https://bookdown.org/animestina/phd_july_19/testing-the-assumptions.html
# https://www.tidytextmining.com/topicmodeling.html
# https://cran.r-project.org/web/packages/glmmTMB/vignettes/glmmTMB.pdf
# https://cran.r-project.org/web/packages/bestNormalize/vignettes/bestNormalize.html

