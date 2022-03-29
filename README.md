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
