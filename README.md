This repository contains our analysis of the paper "Negativity Spreads More than Positivity on Twitter after both Positive and Negative Political Situations" by Schöne et al. (2021) done in the context of the lecture "Applied Network Science: Social Media Networks".

**Read `DOCUMENTATION.pdf` for more information.**

# Repository Description

The repository has the following structure

```
.
├── README.md
├── DOCUMENTATION.pdf
├── cited_papers/
├── our_work/
└── paper_material/
```

List of files on the top level directory:
- `README.md`: this text file
- `DOCUMENTATION.pdf`: documentation of our work

## `cited_papers` directory

```
cited_papers/
├── Bail et al., 2018.pdf
├── Hutto et al., 2014.pdf
├── Ribeiro et al., 2016.pdf
├── Schoöne et al. 2021 supp.pdf
├── Schöne et al., 2021.pdf
├── Sentiment Analysis ML.pdf
├── Sentiment Analysis RF.pdf
├── Thelwall et al., 2010.pdf
├── Thelwall, 2014.pdf
└── Zimbra et al., 2016.pd
```

Contains cited papers used in our work and in the paper from Schöne et al.


## `our_work` directory

```
our_work/
├── analysis_paper.md
├── presentation/
│   └── Social_Networks_Presentation.pdf
└── scripts/
    ├── Social_Networks_Analysis1.R
    ├── Social_Networks_Analysis2.ipynb
    └── original_data -> ../../paper_material/data/
```

Contains our work:
- `analysis_paper.md`: analysis of the paper from Schöne et al.
- `presentation/Social_Networks_Presentation.pdf`: our final presentation
- `scripts/Social_Networks_Analysis1.R`: R Script for statistical analysis written by Alex Timans, see `DOCUMENTATION.pdf` for more information
- `scripts/data_analysis.ipynb`: Jupyter Notebook for exploratory data analysis written by Samuel Anzalone, see `DOCUMENTATION.pdf` for more information

## `paper_material` directory

```
paper_material/
├── data/
│   ├── ferguson_anonymise_cleaned.csv
│   └── ...
└── r_scripts/
    ├── analysis_manuscript.Rmd
    └── ...
```

Contains the supplementary material the paper from Schöne et al.:
- `data/`: contains the tweet datasets
- `r_scripts/`: contains the R Scripts used