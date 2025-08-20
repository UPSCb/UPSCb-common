# R-Studio

As mentioned in the [Onboarding](2.%20Onboarding.md) tab, we host our R-Studio on the server using Posit. It is hosted on `https://aspseq.upsc.se:8080/`. 

>!!! Reminder
R-Studio should be used only for data exploration and visualization. ==**DO NOT LOG ON TO R-Studio**== to

> - Run data analysis pipeline
> - Run or submit `SLURM` jobs

### Running memory intensive R scripts

Sometimes you actually need lots of memory to successfully run an R script. No problem, just *__do not run__* it on AspSeq! 

In this case you can run the R script using `SLURM`. Here is how:

1. Save your R script to an R file, and make sure that it exports its results to a file.

2. Write an `.sh` file like the following one, which contains the instructions to execute your R script e.g 'my_script.R':

```bash
#!/bin/bash
#SBATCH --account=YOUR_PROJECT_CODE
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --mem=10G

R=/mnt/picea/home/singularity/R-4.4.3.sif 
Rlib=/mnt/picea/home/rstudio/Modules/apps/compilers/R/4.4.3/lib/R/library:/usr/local/lib/R/library

export R_LIBS=/mnt/picea/home/rstudio/Modules/apps/compilers/R/4.4.3/lib/R/library

apptainer exec -e -B /mnt:/mnt -B $Rlib $R Rscript \
--vanilla my_script.R
```

### Installing packages

Most commonly used packages are already installed on the server. However, there maybe instances where you need to install your own R-Packages. In most cases this is as easy as logging on to `aspseq` and running `>install.packages("package1")`. 

In some cases, where you need to run the R-script through the `SLURM` queue, the process has a few more steps. When run via `SLURM` R is run thorugh a container. So it is important that the packages you install are installed in the library for the right version of R. 


### R scripts

Several R scripts for generating a variety of plots or doing different analysis are available in the UPSCb-common repository and should be available in your home directory once you have cloned it as shown in the [Onboarding](2.%20Onboarding.md) tab. 

Other important scrpits are BiologicalQA.R/Rmd, and DifferentialExpression_WithGOenrichment.Rmd. These are available in the `Templates` directory of the repository and more information about them is in the [Templates](templates.md) tab. 
 

## Resolving Issues by yourself

While we are always happy to help, here are some [Rstudio common issues](https://gist.github.com/nicolasDelhomme/5bde1e878b2eaa3def1cced06076b7db) and how to try and resolve them. 

