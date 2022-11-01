# UPSCb-common
This repository is a library of common scripts and templates created to avoid reinventing the wheel then and again.

## Structure
There are three important directories:

```
|- pipeline
|- src
  |- bash
  |- R
|- templates
  |- R
```

1. `pipeline` is a directory that contains scripts written (mostly) in bash to be copied into your project and used through a SLURM queueing system.
2. `src` (source) contains a number of subdirectories sorted by programming language holding a number of utilities
3. `templates` contains a number of subdirectories sorted by programming language holding a number of templates to be re-used. Most notably it contains examples of data analysis in R (_e.g._ BiologicalQA.R, DifferentialExpression.R, _etc._) and some bash examples to submit jobs to the queue.

# Onboarding
This instructions are specific to the users at the Umeå Plant Science Centre, you are welcome to get inspired by them though.

## Access
1. Follows the instructions [there](https://youtu.be/hYtIKIIwRss) and our discussions in Slack.

## Setup your project
We use `git` and `SLURM` to ensure reproducible research, here are gists on how to enable that in your projects:

1. [Git setup](https://gist.github.com/nicolasDelhomme/46a1053d277510b95692318bd1732b6d)
2. [SLURM usage](https://gist.github.com/nicolasDelhomme/6fbff1e4db3c7ee4b3bb4f710667fd0d)
3. A [video](https://youtu.be/3XMHTixiszE) of a tech seminar on both (duration ~1h)
4. A description of our server structure
5. A video on how to [download data from Novogene](https://youtu.be/A6JcORYs9L0)

## Rules of conduct

# Visulization examples
1. [Gene Ontology](https://gist.github.com/amnzr/7d859ae127c30e13fef3198c20287da2) 

# Templates
Before **re-inventing the wheel**, check the templates directory! A number of useful templates are available there:

1. `R/empty.R` to initiate an R script with the Rmd header and session info blocks.
2. `R/BiologicalQA.R` to do the initial Exploratory Data Analysis (EDA) of your RNA-Seq data
3. `R/DifferentialExpression.R` to perform a Differential Expression (DE) of your data (follows the previous one)
4. `bash/runTemplate.sh` to prototype script to be run on an HPC (High Perf. Compute) using SLURM (a job manager)

# Troubleshooting

## Resolving Issues by yourself
1. [Rstudio common issues](https://gist.github.com/nicolasDelhomme/5bde1e878b2eaa3def1cced06076b7db)




## Contact us
For UPSC members, ask us to be added to our Slack channel as well as mailing list. These are the two channels we use to communicate about server updates and downtime (as well as other technical issues), but also those we use to discuss projects, provide support, _etc._
