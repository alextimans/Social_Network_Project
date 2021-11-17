# Social Network Project
## Progress board

The progress board can be accessed [here](https://github.com/SamuelsGitHub/Social_Network_Project/projects/1).

## Detailed analysis of paper

### Dataset

- `df_trump.csv`: 28.42% users wrote multiple tweets, so 
278188 - 82970 * (1 - 0.2842) = 218798.1 tweets written by 23580 users (avg 9.28 tweets/user which tweet multiple times)
- `df_hillary.csv`: 29.88% users wrote multiple tweets, so
346730 - 161194 * (1 - 0.2988) = 233700.8 tweets written by 45811.33 users (avg 5.1 tweets/user which tweet multiple times)
- `df_ferguson.csv`: 39.47% users wrote multiple tweets, so
552911 - 181179 * (1 - 0.3947) = 443243.4 tweets written by 71511.35 users (avg 6.2 tweets/user which tweet multiple times)
- `df_samesex.csv`: 31.89% users wrote multiple tweets, so
980018 - 389111 * (1 - 0.3189) = 714994.5 tweets written by 124087.5 users (avg 5.76 tweets/user which tweet multiple times)

**Conclusion:** Dataset not as diverse as it seems, although still somewhat OK

### SentiStrength

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

## Rest coming soon