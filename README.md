## Social Network Project

This repository contains our analysis of the paper "Negativity Spreads More than Positivity on Twitter after both Positive and Negative Political Situations" by Schöne et al. (2021) in the context of the lecture "Applied Network Science: Social Media Networks" at ETH Zurich, winter semester 2021.

─> Read `DOCUMENTATION.pdf` for more information.

### Repository structure
This repository has the following structure:
```
├── README.md
├── DOCUMENTATION.pdf
├── cited_papers/
├── our_work/
└── paper_material/
```

### Cited papers
Includes the main paper we reviewed by Schöne et al. (2021).
```
cited_papers/
├── Bail et al., 2018.pdf
├── Hutto et al., 2014.pdf
├── Ribeiro et al., 2016.pdf
├── Schöne et al. 2021 supp.pdf
├── Schöne et al., 2021.pdf
├── Sentiment Analysis ML.pdf
├── Sentiment Analysis RF.pdf
├── Thelwall et al., 2010.pdf
├── Thelwall, 2014.pdf
└── Zimbra et al., 2016.pd
```

### Contributions
```
our_work/
├── presentation/
│   └── Social_Networks_Presentation.pdf
└── scripts/
    ├── Social_Networks_Analysis1.R
    ├── Social_Networks_Analysis2.ipynb
    └── original_data -> ../../paper_material/data/
```
Contains our work:
- `Social_Networks_Presentation.pdf`: our final presentation
- `Social_Networks_Analysis1.R`: R Script for statistical analysis written by Alexander Timans
- `Social_Networks_Analysis2.ipynb`: Jupyter Notebook for exploratory data analysis written by Samuel Anzalone

### Additional directory
```
paper_material/
├── data/
│   ├── ferguson_anonymise_cleaned.csv
│   └── ...
└── r_scripts/
    ├── analysis_manuscript.Rmd
    └── ...
```
Contains the supplementary material from the paper of Schöne et al. (2021):
- `data/`: contains the tweet datasets
- `r_scripts/`: contains R Scripts
