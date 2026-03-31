# Bash templates

This directory lists bash boiler plate code (_i.e._ templates) that can be reused to avoid reinventing the wheel.

## Rationale

Really, just to avoid reinventing the wheel, speed up implementation and benefit from a common structure and others' input and knowledge. That means contributions are more than welcome :-)

## Usage

Most of the template will have very clear comments on what needs to be changed and tokens to be replaced to adjust the template to your need. Some of these templates are fairly old (check the last update date, it is also mentioned in this doc).

### Cleanup

File:

* `cleanup.sh`

This is not a template per-se, but something that could be turned into one. The reason it is here is that it contains code that can be used to clear-up project at scale, leveraging SLURM Job Arrays for example.

### Job Array

Files:

* `jobArrayTemplate.sh`
* `jobArrayTask.sh`

This is a couple files to implement a simple job array. The jobArrayTemplate is the script that will submit the job array, while the jobArrayTask is the unit task to be executed. The way it is implemented now makes it rely on getting a list of the input files to set the size of the array (in the JobArrayTemplate) and select the right input for the task (in the jobArrayTask). A simpler way is commented out in the templates, where if you name your input files as `prefix.task_number` where task_number goes from 0 to the total number of tasks to run X minus 1, you can bypass listing input files.

### Script template

File:

* `template.sh`

This is a script used as a wrapper around the command you want to run. It enables you to build a list of commands to run and by default perform a dry-run. It handles basic sanity checks too. Make sure to adjust the variable at the top `ARGs` to the number of arguments your script expects (not options, really just the size of `$@`, _i.e._ `$#`, once options have been processed).

### Seidr

Files:

* submitSeidrAggregate.sh
* submitSeidrBackbone.sh
* submitSeidRoc.sh

These files are old template meant to help running some of the steps from `seidr`. These were set to run on the UPSCb bioinfo cluster. They can easily be reused as there is a seidr singularity container available on both the UPSCb clusters and HPC2N.

### Semaphore

Files:

* process.sh
* semaphore.sh

This is a couple scripts that allows setting up a semaphore to run jobs in parallel in a mostly automated way. This is relevant for tasks where one cannot use a queueing system, such as data transfer between HPCs, ENA data submission, or the like.

### Singularity runner

File:

* runTemplate.sh

This is an old template, ancestor of the `template.sh` file [mentioned above](#script-template). It only has some extra handling for Singularity/Apptainer that can easily be ported to the newer template.
