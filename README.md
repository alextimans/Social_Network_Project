# Social Network Project
## To Do

The progress board can be accessed [here](https://github.com/SamuelsGitHub/Social_Network_Project/projects/1). The ToDo's have been inspired by the `To Do` section in `analysis_paper.md`.

## Detailed analysis of paper
(taken from `analysis_paper.md`)

### 1. Dataset

- `df_trump.csv`: 28.42% users wrote multiple tweets, so 
278188 - 82970 * (1 - 0.2842) = 218798.1 tweets written by 23580 users (avg 9.28 tweets/user which tweet multiple times)
- `df_hillary.csv`: 29.88% users wrote multiple tweets, so
346730 - 161194 * (1 - 0.2988) = 233700.8 tweets written by 45811.33 users (avg 5.1 tweets/user which tweet multiple times)
- `df_ferguson.csv`: 39.47% users wrote multiple tweets, so
552911 - 181179 * (1 - 0.3947) = 443243.4 tweets written by 71511.35 users (avg 6.2 tweets/user which tweet multiple times)
- `df_samesex.csv`: 31.89% users wrote multiple tweets, so
980018 - 389111 * (1 - 0.3189) = 714994.5 tweets written by 124087.5 users (avg 5.76 tweets/user which tweet multiple times)

**Conclusion:** Dataset not as diverse as it seems, although still somewhat OK

### 2. SentiStrength

From [SentiStrength paper (Thelwall, Buckley et al., 2010)](https://github.com/SamuelsGitHub/Social_Network_Project/blob/main/paper_sentistrength.pdf):

- "Overall, it seems that SentiStrength is not good at identifying negative emotion"; "results for negative sentiment are not significant"
- Accuracy for positive emotion 60.6%, for negative emotion 72.8% (2% above baseline, slightly outperformed by SVM
- Rule- and lexicon-based algorithm (not ML)
- "Gold standard", i.e. test labels: mean of sentiment scores as assigned by three female operators

From SentiStrength follow-up (Someone et al., resource where ???):

- SentiStrength is outperformed on Twitter data by "Best ML algorithm" for positive (60.7% vs. 56.3%) and negative (64.3% vs. 61.7%) emotion accuracy and also has lower correlations with human-coded scores (~2% difference)

From DAN2 paper (Someone et al., 2016, resource where ???):

- DAN2 algorithm based on Dynamic ANNs with feature engineering achieves 85.5% accuracy for a 5-class sentiment classification on Twitter data, outperforming SVMs (the standard industry ML algo for this) by 7%.

From SentiBench/Ribero paper [SentiBench paper (Ribeiro et al., 2016)](https://github.com/SamuelsGitHub/Social_Network_Project/blob/main/paper_sentibench.pdf):

- Paper was referenced by Jonas as justification for using SentiStrength
- For a 3-class sentiment classification on the Twitter part from a "benchmark" data set, SentiStrength is clearly outperformed by VADER; Macro-F1 46.77% vs. 58.12%, Acc 57.83% vs. 60.21%
- For 2-class sentiment classification on the same data, SentiStrength shows terrible coverage in comparison to VADER, i.e. most tweets are not classified; Coverage 34.65% vs. 94.4%, Acc 96.97% vs. 99.04%;
- For general sentiment analysis methods on 3 classes, VADER ranks #1 and SentiStrength ranks #15
- For social network-specific sentiment analysis on 3 classes, VADER ranks #3 and SentiStrength ranks #11
- From the paper: "Although SentiStrength presented good Macro-F1 values, its coverage is usually low as this method tends to classify a high number of instances as neutral". Thus when Jonas filters out low valences, he indeed removes drastically large amounts of data.
- Jonas even used and compared VADER to SentiStrength, and even though they (only) have a correlation of 52-62% in the scoring, he decided to stick to SentiStrength. Not sure why, given the above arguments.
- Jonas mentions importance of accurate sentiment analysis tool in "Limitations", but then quotes wrong number about how good SentiStrength is referring to the SentiBench/Ribeiro paper.

**Conclusions:**

- Considering Jonas’ research is from 2021, he could have used more recent and powerful approaches to obtain sentiment intensity scores
- Sentiment analysis libraries for Python: NLTK, SpaCy, TextBlob, CoreNLP, gensim (topic analysis), Flair

### 3. Preprocessing + valence analysis

- SentiStrength returns scores in [1,5] for positive emotion and [-5,-1] for negative emotion. Jonas decided to scale both to [0, 4], i.e. subtract 1 and take absolute values for the retweet predictions.
- His argument: easier to interpret regression model intercepts. But then again he uses retweets_reciprocal which makes direct interpretation of coefficients virtually impossible.
- Does this scaling change the effect of the features?

- Valence analysis is essentially comparing the means of sen_pos and sen_neg and quantifying a significant difference. The important part in the LMER output are the coefficients and specifically their signs. Intercept = mean of sen_neg.
- Why not simply do a standard paired t.test instead? `t.test(x=df$sen_neg, y=df$sen_pos, alternative="greater", paired=TRUE)` gives the desired result and is much more straightforward.
- `t.test` yields different result for `d.trump`?

### 4. Retweet prediction 

- Why is user’s #followers included, if Jonas explicitly mentions "users with more followers generally have more retweets regardless of the emotional content of their tweets"
- The variable has (relatively) high explanatory power and contains most of the explained variance of his final fitted models (see e.g. drop1 RSS)
- Is it included just to bump up the R^2 values to make it more believable?
- Features sen_pos and sen_neg were not encoded as factors?

#### Models and transformations tested by Jonas: 
- `lmer + retweets_reciprocal` (FINAL MODEL)
- `lmer + retweets_log`
- `glmer + family=poisson(link="log") + retweets`
- `glmer.nb (neg binomial) + retweets`
- `GAMM + family=gaussian + retweets_reciprocal + s(sen_pos/neg, k=4)`

#### Potential models/transformations (relating to skewed data):
- Simple normalization of data
- Square root transformation
- Box-Cox transformation
- Transformation of covariates (see point below)
- Multinomial/categorical distribution
- Pareto distribution
- Ordinal regression (e.g. Cumulative Link Mixed Model CLMM)
- Non-parametric/non-linear models (Comp Stats-style)

Personal opinion:
- Seems to me Jonas decided to use Linear Mixed Models without questioning
- Potentially misspecified model and faulty p-value interpretation
- TA plot looks weird (see reference: https://stats.stackexchange.com/questions/260081/residual-pattern-in-mixed-model) -> Covariate transform?
- QQ plot looks like systematic deviation
- Are Linear Mixed Model assumptions violated ?
- The fact that the target variable has to be so strongly transformed to make it work may hint that the true relationship is in fact non-linear.


#### Incorrect backtransform (MISTAKE)
- Initially, the model outputs for retweet prediction were not accurately reproducible, values being slightly off. Then I found out why.
- Jonas did a mistake when computing retweets from the backtransform of `retweets_log`: he used `exp(retweets_log-1)` instead of `exp(retweets_log)-1`
- when using `retweets = round(exp(retweets_log - 1))` and running the models with `retweets_reciprocal` based on these incorrect values, I obtain his results.
- Not a huge difference in the final interpretation/significance, but coefficient values etc. definitely change
- Since he made this mistake, the estimated effects used for the plots are also wrong.

### 5. Political affiliation
Comments are mainly based on `d.trump`.

- The random sample of 30830 tweets represents only 11% of d.trump, it could be an unfortunate sample, it is not evident to generalise the observed affiliation distribution to whole dataset.
- For `d.samesex`, the sample is 4816 tweets and is supposedly generalisable to dataset of size 980 018.

- Jonas writes "able to estimate affiliation for 85.9% users tested" and "tested 30830 tweets"; but then goes on to fully separate ALL the 30830 tweets into Dem/Rep, not 85.9% of them?

- Binary assignment to Rep/Dem was made purely based on `max(#Dem follows, #Rep follows)`, and if equal then Dem ("Mid" labels).
- This is a very naive way of assigning affiliation: consider someone following 51 Rep and 50 Dem from the list -> Rep; consider someone following 1000 people of which 1 is Dem and 2 are Rep on the list -> Rep;
- More sensible to work with thresholds, e.g. min. ratio of recognized accounts/total followed accounts, and min. difference ratio between Dem and Rep accounts to state a valid affiliation.

- List of accounts labelled Dem/Rep from Bail et al. is unbalanced (more Dems) and may not be extensive enough to properly base conclusions
- Furthermore, their public figure list attributes accounts a liberal/conservative score on a spectrum, NOT a partisan Dem/Rep label.
- Thus when saying that most tweets were made by "conservatives", this does not exclude Dem conservatives and vice versa.
- Claims for intragroup vs. intergroup effects cannot be made because we cannot properly disentangle political views from party membership.
- Jonas mentions the fact that intragroup interactions may have not been properly captured in "Limitations" - I agree.

- When filtering data for the stated subset (30830 tweets) with correct/incorrect retweet values I only obtain 11932 tweets. 
- When filtering for "Rep" one can see that the stated conditions (min. 1 retweet, min 2 intensity score) are not fulfilled, explaining why the subsample is larger than when I filter precisely. Why?

- Rerunning the linear model on the "Rep" subset isn’t evidence.
- Jonas writes "repeated analysis on the subsample of users identified as conservatives … analysis led to same results as more inclusive analysis". 
- This is not true? Subsample analysis returns everything significant (incl. interaction and sen_pos) -> can one even evaluate main effects then?

—> Using the above points, it is a far reach to claim having intragroup effects. No additional statistical test for this were done.

### 6. Text analysis (content of tweets)
Comments mainly based on d.trump.

- Issues with assignment of data to election victory/defeat based on very limited amount of hashtags
- E.g. multiple negative tweets (order by valence) from Hillary supporters using #MAGA or #PresidentTrump, i.e. again strong entanglement of intergroup/intragroup effects.

- One can use examples from d.samesex.text and d.trump.text to show that SentiStrength is doing a somewhat poor job. E.g. look at content of samesex tweets ranked as both sen_pos=5 and sen_neg=-5.
- These tweets then get a valence (sum of both) of 0, classifying them as "neutral" even though they clearly are NOT emotionally neutral.
- Jonas’ valence approach is a very naive way of assigning a binary label, it relies too much on the inaccurate scoring of SentiStrength.

- By later filtering out tweets classified as "neutral" he essentially discards a big chunk of data that is very expressive i.e. has high emotional intensity both ways.
- Instead he only selects those that go either way, but as we can see from examples these are somewhat arbitrary and not well classified by SentiStrength.
- Filtering out "neutral" for samesex unigrams decreases data from 446133 to 148007 samples (loses over 2/3 of data!)

- When you look at samesex tweets with valence=-4 (thus classified as very negative) it is many people being super happy and "crying tears of joy … LoveWins".
- When you look at "negative" top viral tweets for samesex, most of them are super positive.
- When looking at ap_top_terms for samesex bigrams for the two topics, they are mostly all positive -> more like one topic.

- Binarization of data into pos/neg doesn’t make use of [1,5] scores
- Maybe Jonas should have used a sentiment analysis tool that does a simpler binary classification with higher accuracy rather than the granular classification and then aggregating these values which can produce inconsistent assignments.
- Loss of much info by going from: 
2 factors, 10 lvls (sen_neg, sen_pos) -> 1 factor, 9 lvls (valence in [-4,4])
-> 1 factor, 3 lvls (pos, neg, neutral) -> 1 factor, 2 lvls (filter neutral)

- Is looking at the most unique unigrams/bigrams really representative of the topic content? 
- Isn’t that kind of cherrypicking, where you select the most extreme and potentially rare cases? Shouldn’t one just look at the most frequent grams/highest beta score?
- In that case, both topics would be revealed to pretty much be intermixed, without a clear distinction and with much overlap of grams, even though the gamma values supposedly claim a clear distinction.
- High gamma values claiming "means that words appearing in a category have higher prob to only be used in one of the topics" are disproved by the overlap of ap_top_terms?
- There is something strange going on with the LDA and only having 2 topics when giving it data with a binary labeled topic column... I get the feeling that the results are somewhat predetermined by the inputs and topic size.
- There is very high accordance of ap_top_terms with the gay_dtm values with highest counts, i.e. words that have high counts are also identified as words representing the topic (makes sense) 
-> Again, unique words are cherrypicking? Or is this common practice in topic modelling?
