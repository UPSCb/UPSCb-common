#!/usr/bin/env nextflow 

expr = file(params.expr)
genes = file(params.genes)
targets = params.targets
if (targets != '' ) {
  targets = file(params.targets)
  targets = '-t ' + targets
}

if (params.mi) {

  process mi {
    errorStrategy params.error_strategy
    cpus params.mi_settings.cores
    memory params.mi_settings.memory
    publishDir params.out + '/networks/mi'
    clusterOptions '-n ' + params.mi_settings.tasks + ' -A ' + params.slurm_account
    time params.mi_settings.ptime

    input:
    file expr
    file genes
    val targets

    output:
    file 'mi_network.tsv' into mi_network_raw, mi_for_aracne, mi_for_clr
    
    script:
    if (params.executor == 'slurm')
    {
      """
      ${params.mi_settings.preamble}
      srun \
      mi -m RAW -o mi_network.tsv \
         -b ${params.mi_settings.bins} -B $(expr $(expr ${params.ngenes} '/' ${ params.mi_settings.tasks}) '+' 1) -s ${params.mi_settings.spline} \
         -i ${expr} -g ${genes} -O ${params.mi_settings.cores} ${targets}
      ${params.mi_settings.epilog}
      """
    }
    else
    {
      """
      ${params.mi_settings.preamble}
      mpirun -np ${params.mi_settings.tasks} \
      mi -m RAW -o mi_network.tsv \
         -b ${params.mi_settings.bins} -s ${params.mi_settings.spline} \
         -i ${expr} -g ${genes} -O ${params.mi_settings.cores} ${targets}
      ${params.mi_settings.epilog}
      """
    }
  }

  process mi_import {
    errorStrategy params.error_strategy
    validExitStatus 0,3
    cpus params.mi_settings.importcores
    memory params.mi_settings.importmem
    publishDir params.out + '/networks/mi'
    clusterOptions '-n ' + params.mi_settings.tasks + ' -A ' + params.slurm_account
    time params.mi_settings.itime

    input:
    file genes
    file mi_net from mi_network_raw
    output:
    file 'mi_network.sf' into mi_sf

    script:
    if (targets == '')
    {
      """
      ${params.mi_settings.preamble}
      seidr import -F lm -u -r -z -n ${params.mi_settings.importname} \
                   -i ${mi_net} -g ${genes} -o mi_network.sf \
                   -O ${params.mi_settings.importcores}
      ${params.mi_settings.epilog}
      """
    }
    else
    {
      """
      ${params.mi_settings.preamble}
      seidr import -F el -u -r -z -n ${params.mi_settings.importname} \
                   -i ${mi_net} -g ${genes} -o mi_network.sf \
                   -O ${params.mi_settings.importcores}
      ${params.mi_settings.epilog}
      """
    }
  }
} else { 
  Channel.empty().set {mi_sf} 
}

if (params.clr) {

  process clr {
    errorStrategy params.error_strategy
    cpus params.clr_settings.cores
    memory params.clr_settings.memory
    publishDir params.out + '/networks/clr'
    clusterOptions '-n ' + params.clr_settings.tasks + ' -A ' + params.slurm_account
    time params.clr_settings.ptime

    input:
    file expr
    file genes
    file mi_net from mi_for_clr
    val targets

    output:
    file 'clr_network.tsv' into clr_network_raw
    
    script:
    if (params.executor == 'slurm')
    {
      """
      ${params.clr_settings.preamble}
      srun \
      mi -m CLR -o clr_network.tsv \
         -M ${mi_net} -B $(expr $(expr ${params.ngenes} '/' ${ params.clr_settings.tasks}) '+' 1) \
         -b ${params.clr_settings.bins} -s ${params.clr_settings.spline} \
         -i ${expr} -g ${genes} -O ${params.clr_settings.cores} ${targets}
      ${params.clr_settings.epilog}
      """
    }
    else
    {
      """
      ${params.clr_settings.preamble}
      mpirun -np ${params.clr_settings.tasks} \
      mi -m CLR -o clr_network.tsv \
         -M ${mi_net} \
         -b ${params.clr_settings.bins} -s ${params.clr_settings.spline} \
         -i ${expr} -g ${genes} -O ${params.clr_settings.cores} ${targets}
      ${params.clr_settings.epilog}
      """
    }
  }

  process clr_import {
    errorStrategy params.error_strategy
    validExitStatus 0,3
    cpus params.clr_settings.importcores
    memory params.clr_settings.importmem
    publishDir params.out + '/networks/clr'
    clusterOptions '-n ' + params.clr_settings.tasks + ' -A ' + params.slurm_account
    time params.clr_settings.itime

    input:
    file genes
    file clr_net from clr_network_raw
    output:
    file 'clr_network.sf' into clr_sf

    script:
    if (targets == '')
    {
      """
      ${params.clr_settings.preamble}
      seidr import -F lm -u -r -z -n ${params.clr_settings.importname} \
                   -i ${clr_net} -g ${genes} -o clr_network.sf \
                   -O ${params.clr_settings.importcores}
      ${params.clr_settings.epilog}
      """
    }
    else
    {
      """
      ${params.clr_settings.preamble}
      seidr import -F el -u -r -z -n ${params.clr_settings.importname} \
                   -i ${clr_net} -g ${genes} -o clr_network.sf \
                   -O ${params.clr_settings.importcores}
      ${params.clr_settings.epilog}
      """
    }
  }
} else { 
  Channel.empty().set {clr_sf} 
}

if (params.aracne) {

  process aracne {
    errorStrategy params.error_strategy
    cpus params.aracne_settings.cores
    memory params.aracne_settings.memory
    publishDir params.out + '/networks/aracne'
    clusterOptions '-n ' + params.aracne_settings.tasks + ' -A ' + params.slurm_account
    time params.aracne_settings.ptime

    input:
    file expr
    file genes
    file mi_net from mi_for_aracne
    val targets

    output:
    file 'aracne_network.tsv' into aracne_network_raw

    script:
    if (params.executor == 'slurm')
    {
      """
      ${params.aracne_settings.preamble}
      srun \
      mi -m ARACNE -o aracne_network.tsv \
         -M ${mi_net} -B ${params.ngenes} -B -B $(expr $(expr ${params.ngenes} '/' ${ params.aracne_settings.tasks}) '+' 1)  \
         -b ${params.aracne_settings.bins} -s ${params.aracne_settings.spline} \
         -i ${expr} -g ${genes} -O ${params.aracne_settings.cores} ${targets}
      ${params.aracne_settings.epilog}
      """
    }
    else
    {
      """
      ${params.aracne_settings.preamble}
      mpirun -n ${params.aracne_settings.tasks} \
      mi -m ARACNE -o aracne_network.tsv \
         -M ${mi_net} \
         -b ${params.aracne_settings.bins} -s ${params.aracne_settings.spline} \
         -i ${expr} -g ${genes} -O ${params.aracne_settings.cores} ${targets}
      ${params.aracne_settings.epilog}
      """
   }
  }

  process aracne_import {
    errorStrategy 'finish'
    validExitStatus 0,3
    cpus params.aracne_settings.importcores
    memory params.aracne_settings.importmem
    publishDir params.out + '/networks/aracne'
    clusterOptions '-n ' + params.aracne_settings.tasks + ' -A ' + params.slurm_account
    time params.aracne_settings.itime

    input:
    file genes
    file aracne_net from aracne_network_raw

    output:
    file 'aracne_network.sf' into aracne_sf

    script:
    if (targets == '')
    {
      """
      ${params.aracne_settings.preamble}
      seidr import -F lm -u -r -z -n ${params.aracne_settings.importname} \
                   -i ${aracne_net} -g ${genes} -o aracne_network.sf \
                   -O ${params.aracne_settings.importcores}
      ${params.aracne_settings.epilog}
      """
    }
    else
    {
      """
      ${params.aracne_settings.preamble}
      seidr import -F el -u -r -z -n ${params.aracne_settings.importname} \
                   -i ${aracne_net} -g ${genes} -o aracne_network.sf \
                   -O ${params.aracne_settings.importcores}
      ${params.aracne_settings.epilog}
      """
    }
  }
} else { 
  Channel.empty().set {aracne_sf} 
}

if (params.anova) {

  process anova {
    errorStrategy params.error_strategy
    cpus params.anova_settings.cores
    memory params.anova_settings.memory
    publishDir params.out + '/networks/anova'
    clusterOptions '-n ' + params.anova_settings.tasks + ' -A ' + params.slurm_account
    time params.anova_settings.ptime

    input:
    file expr
    file genes
    file params.anova_settings.meta_file
    val targets

    output:
    file 'anova_network.tsv' into anova_network_raw

    """
    ${params.anova_settings.preamble}
    export OMP_NUM_THREADS=1
    anoverence -i ${expr} -g ${genes} -e ${params.anova_settings.meta_file} \
               -w ${params.anova_settings.weight} ${targets}
    ${params.anova_settings.epilog}
    """
  }

  process anova_import {
    errorStrategy params.error_strategy
    validExitStatus 0,3
    cpus params.anova_settings.importcores
    memory params.anova_settings.importmem
    publishDir params.out + '/networks/anova'
    clusterOptions '-n 1 -A ' + params.slurm_account
    time params.anova_settings.itime

    input:
    file genes
    file anova_net from anova_network_raw

    output:
    file 'anova_network.sf' into anova_sf

    script:
    if (targets == '')
    {
      """
      ${params.anova_settings.preamble}
      seidr import -F lm -u -r -z -n ${params.anova_settings.importname} \
                   -i ${aracne_net} -g ${genes} -o anova_network.sf \
                   -O ${params.anova_settings.importcores}
      ${params.anova_settings.epilog}
      """
    }
    else
    {
      """
      ${params.anova_settings.preamble}
      seidr import -F el -u -r -z -n ${params.anova_settings.importname} \
                   -i ${aracne_net} -g ${genes} -o anova_network.sf \
                   -O ${params.anova_settings.importcores}
      ${params.anova_settings.epilog}
      """
    }
  }
} else { 
  Channel.empty().set {anova_sf} 
}

if (params.pearson) {

  process pearson {
    errorStrategy params.error_strategy
    cpus params.pearson_settings.cores
    memory params.pearson_settings.memory
    publishDir params.out + '/networks/pearson'
    clusterOptions '-n ' + params.pearson_settings.tasks + ' -A ' + params.slurm_account
    time params.pearson_settings.ptime

    input:
    file expr
    file genes
    val targets

    output:
    file 'pearson_network.tsv' into pearson_network_raw

    """
    ${params.pearson_settings.preamble}
    export OMP_NUM_THREADS=1
    correlation -m pearson -i ${expr} -g ${genes} -o pearson_network.tsv \
                ${params.pearson_settings.scale} \
                ${params.pearson_settings.absolute} ${targets}
    ${params.pearson_settings.epilog}
    """
  }

  process pearson_import {
    errorStrategy params.error_strategy
    validExitStatus 0,3
    cpus params.pearson_settings.importcores
    memory params.pearson_settings.importmem
    publishDir params.out + '/networks/pearson'
    clusterOptions '-n 1 -A ' + params.slurm_account
    time params.pearson_settings.itime

    input:
    file genes
    file pearson_net from pearson_network_raw

    output:
    file 'pearson_network.sf' into pearson_sf

    script:
    if (targets == '')
    {
      """
      ${params.pearson_settings.preamble}
      seidr import -F lm -A -u -r -z -n ${params.pearson_settings.importname} \
                   -i ${pearson_net} -g ${genes} -o pearson_network.sf \
                   -O ${params.pearson_settings.importcores}
      ${params.pearson_settings.epilog}
      """
    }
    else
    {
      """
      ${params.pearson_settings.preamble}
      seidr import -F el -A -u -r -z -n ${params.pearson_settings.importname} \
                   -i ${pearson_net} -g ${genes} -o pearson_network.sf \
                   -O ${params.pearson_settings.importcores}
      ${params.pearson_settings.epilog}
      """
    }
  }
} else { 
  Channel.empty().set {pearson_sf} 
}

if (params.spearman) {

  process spearman {
    errorStrategy params.error_strategy
    cpus params.spearman_settings.cores
    memory params.spearman_settings.memory
    publishDir params.out + '/networks/spearman'
    clusterOptions '-n ' + params.spearman_settings.tasks + ' -A ' + params.slurm_account
    time params.spearman_settings.ptime

    input:
    file expr
    file genes
    val targets

    output:
    file 'spearman_network.tsv' into spearman_network_raw

    """
    ${params.spearman_settings.preamble}
    export OMP_NUM_THREADS=1
    correlation -m spearman -i ${expr} -g ${genes} -o spearman_network.tsv \
                ${params.spearman_settings.scale} \
                ${params.spearman_settings.absolute} ${targets}
    ${params.spearman_settings.epilog}
    """
  }

  process spearman_import {
    errorStrategy params.error_strategy
    validExitStatus 0,3
    cpus params.spearman_settings.importcores
    memory params.spearman_settings.importmem
    publishDir params.out + '/networks/spearman'
    clusterOptions '-n 1 -A ' + params.slurm_account
    time params.spearman_settings.itime

    input:
    file genes
    file spearman_net from spearman_network_raw

    output:
    file 'spearman_network.sf' into spearman_sf

    script:
    if (targets == '')
    {
      """
      ${params.spearman_settings.preamble}
      seidr import -F lm -A -u -r -z -n ${params.spearman_settings.importname} \
                   -i ${spearman_net} -g ${genes} -o spearman_network.sf \
                   -O ${params.spearman_settings.importcores}
      ${params.spearman_settings.epilog}
      """
    }
    else
    {
      """
      ${params.spearman_settings.preamble}
      seidr import -F el -A -u -r -z -n ${params.spearman_settings.importname} \
                   -i ${spearman_net} -g ${genes} -o spearman_network.sf \
                   -O ${params.spearman_settings.importcores}
      ${params.spearman_settings.epilog}
      """
    }
  }
} else { 
  Channel.empty().set {spearman_sf} 
}

if (params.elnet) {

  process elnet {
    errorStrategy params.error_strategy
    cpus params.elnet_settings.cores
    memory params.elnet_settings.memory
    publishDir params.out + '/networks/elnet'
    clusterOptions '-n ' + params.elnet_settings.tasks + ' -A ' + params.slurm_account
    time params.elnet_settings.ptime

    input:
    file expr
    file genes
    val targets

    output:
    file 'elnet_network.tsv' into elnet_network_raw


    script:
    if (params.executor == 'slurm')
    {
      """
      ${params.elnet_settings.preamble}
      srun  \
      el-ensemble -o elnet_network.tsv \
         -l ${params.elnet_settings.min_lambda} -a ${params.elnet_settings.alpha} \
         -n ${params.elnet_settings.nlambda} \
         -X ${params.elnet_settings.max_experiment_size} \
         -x ${params.elnet_settings.min_experiment_size} \
         -P ${params.elnet_settings.max_predictor_size} \
         -p ${params.elnet_settings.min_predictor_size} \
         -e ${params.elnet_settings.ensemble} \
         -i ${expr} -g ${genes} \
         -O ${params.elnet_settings.cores} ${targets}
      ${params.elnet_settings.epilog}
      """
    }
    else
    {
      """
      ${params.elnet_settings.preamble}
      mpirun -np ${params.elnet_settings.tasks}  \
      el-ensemble -o elnet_network.tsv \
         -l ${params.elnet_settings.min_lambda} -a ${params.elnet_settings.alpha} \
         -n ${params.elnet_settings.nlambda} \
         -X ${params.elnet_settings.max_experiment_size} \
         -x ${params.elnet_settings.min_experiment_size} \
         -P ${params.elnet_settings.max_predictor_size} \
         -p ${params.elnet_settings.min_predictor_size} \
         -e ${params.elnet_settings.ensemble} \
         -i ${expr} -g ${genes} \
         -O ${params.elnet_settings.cores} ${targets}
      ${params.elnet_settings.epilog}
      """
   }
  }

  process elnet_import {
    errorStrategy params.error_strategy
    validExitStatus 0,3
    cpus params.elnet_settings.importcores
    memory params.elnet_settings.importmem
    publishDir params.out + '/networks/elnet'
    clusterOptions '-n 1 -A ' + params.slurm_account
    time params.elnet_settings.itime

    input:
    file genes
    file elnet_net from elnet_network_raw
    output:
    file 'elnet_network.sf' into elnet_sf

    script:
    if (targets == '')
    {
      """
      ${params.elnet_settings.preamble}
      seidr import -F m -r -z -n ${params.elnet_settings.importname} \
                   -i ${elnet_net} -g ${genes} -o elnet_network.sf \
                   -O ${params.elnet_settings.importcores}
      ${params.elnet_settings.epilog}
      """
    }
    else
    {
      """
      ${params.elnet_settings.preamble}
      seidr import -F el -r -z -n ${params.elnet_settings.importname} \
                   -i ${elnet_net} -g ${genes} -o elnet_network.sf \
                   -O ${params.elnet_settings.importcores}
      ${params.elnet_settings.epilog}
      """
    }
  }
} else { 
  Channel.empty().set {elnet_sf} 
}

if (params.svm) {

  process svm {
    errorStrategy params.error_strategy
    cpus params.svm_settings.cores
    memory params.svm_settings.memory
    publishDir params.out + '/networks/svm'
    clusterOptions '-n ' + params.svm_settings.tasks + ' -A ' + params.slurm_account
    time params.svm_settings.ptime

    input:
    file expr
    file genes
    val targets

    output:
    file 'svm_network.tsv' into svm_network_raw

    script:
    if (params.executor == 'slurm')
    {
      """
      ${params.svm_settings.preamble}
      srun \
      svm-ensemble -o svm_network.tsv \
         -X ${params.svm_settings.max_experiment_size} \
         -x ${params.svm_settings.min_experiment_size} \
         -P ${params.svm_settings.max_predictor_size} \
         -p ${params.svm_settings.min_predictor_size} \
         -e ${params.svm_settings.ensemble} \
         -i ${expr} -g ${genes} \
         -O ${params.svm_settings.cores} ${targets}
      ${params.svm_settings.epilog}
      """
    }
    else
    {
      """
      ${params.svm_settings.preamble}
      mpirun -np ${params.svm_settings.tasks} \
      svm-ensemble -o svm_network.tsv \
         -X ${params.svm_settings.max_experiment_size} \
         -x ${params.svm_settings.min_experiment_size} \
         -P ${params.svm_settings.max_predictor_size} \
         -p ${params.svm_settings.min_predictor_size} \
         -e ${params.svm_settings.ensemble} \
         -i ${expr} -g ${genes} \
         -O ${params.svm_settings.cores} ${targets}
      ${params.svm_settings.epilog}
      """
    }
  }

  process svm_import {
    errorStrategy params.error_strategy
    validExitStatus 0,3
    cpus params.svm_settings.importcores
    memory params.svm_settings.importmem
    publishDir params.out + '/networks/svm'
    clusterOptions '-n 1 -A ' + params.slurm_account
    time params.svm_settings.itime

    input:
    file genes
    file svm_net from svm_network_raw

    output:
    file 'svm_network.sf' into svm_sf

    script:
    if (targets == '')
    {
      """
      ${params.svm_settings.preamble}
      seidr import -F m -r -z -n ${params.svm_settings.importname} \
                   -i ${svm_net} -g ${genes} -o svm_network.sf \
                   -O ${params.svm_settings.importcores}
      ${params.svm_settings.epilog}
      """
    }
    else
    {
      """
      ${params.svm_settings.preamble}
      seidr import -F el -r -z -n ${params.svm_settings.importname} \
                   -i ${svm_net} -g ${genes} -o svm_network.sf \
                   -O ${params.svm_settings.importcores}
      ${params.svm_settings.epilog}
      """
    }
  }
} else { 
  Channel.empty().set {svm_sf} 
}

if (params.llr) {

  process llr {
    errorStrategy params.error_strategy
    cpus params.llr_settings.cores
    memory params.llr_settings.memory
    publishDir params.out + '/networks/llr'
    clusterOptions '-n ' + params.llr_settings.tasks + ' -A ' + params.slurm_account
    time params.llr_settings.ptime

    input:
    file expr
    file genes
    val targets

    output:
    file 'llr_network.tsv' into llr_network_raw

    script:
    if (params.executor == 'slurm')
    {
      """
      ${params.llr_settings.preamble}
      srun \
      llr-ensemble -o llr_network.tsv \
      -B $(expr $(expr ${params.ngenes} '/' ${ params.llr_settings.tasks}) '+' 1)  \
     #    -X ${params.llr_settings.max_experiment_size} \
     #    -x ${params.llr_settings.min_experiment_size} \
     #    -P ${params.llr_settings.max_predictor_size} \
     #    -p ${params.llr_settings.min_predictor_size} \
         -e ${params.llr_settings.ensemble} \
         -i ${expr} -g ${genes} \
         -O ${params.llr_settings.cores} ${targets}
      ${params.llr_settings.epilog}
      """
    }
    else
    {
      """
      ${params.llr_settings.preamble}
      mpirun -np ${params.llr_settings.tasks} \
      llr-ensemble -o llr_network.tsv \
         -i ${expr} -g ${genes} \
         -O ${params.llr_settings.cores} ${targets}
      ${params.llr_settings.epilog}
      """
    }
  }

  process llr_import {
    errorStrategy params.error_strategy
    validExitStatus 0,3
    cpus params.llr_settings.importcores
    memory params.llr_settings.importmem
    publishDir params.out + '/networks/llr'
    clusterOptions '-n 1 -A ' + params.slurm_account
    time params.llr_settings.itime

    input:
    file genes
    file llr_net from llr_network_raw

    output:
    file 'llr_network.sf' into llr_sf

    script:
    if (targets == '')
    {
      """
      ${params.llr_settings.preamble}
      seidr import -F m -r -z -n ${params.llr_settings.importname} \
                   -i ${llr_net} -g ${genes} -o llr_network.sf \
                   -O ${params.llr_settings.importcores}
      ${params.llr_settings.epilog}
      """
    }
    else
    {
      """
      ${params.llr_settings.preamble}
      seidr import -F el -r -z -n ${params.llr_settings.importname} \
                   -i ${llr_net} -g ${genes} -o llr_network.sf \
                   -O ${params.llr_settings.importcores}
      ${params.llr_settings.epilog}
      """
    }
  }
} else { 
  Channel.empty().set {llr_sf} 
}

if (params.pcor) {

  process pcor {
    errorStrategy params.error_strategy
    cpus params.pcor_settings.cores
    memory params.pcor_settings.memory
    publishDir params.out + '/networks/pcor'
    clusterOptions '-n ' + params.pcor_settings.tasks + ' -A ' + params.slurm_account
    time params.pcor_settings.ptime

    input:
    file expr
    file genes
    val targets

    output:
    file 'pcor_network.tsv' into pcor_network_raw

    """
    pcor -i ${expr} -g ${genes} -o pcor_network.tsv \
         ${params.pcor_settings.absolute} ${targets}
    """
  }

  process pcor_import {
    errorStrategy params.error_strategy
    validExitStatus 0,3
    cpus params.pcor_settings.importcores
    memory params.pcor_settings.importmem
    publishDir params.out + '/networks/pcor'
    clusterOptions '-n 1 -A ' + params.slurm_account
    time params.pcor_settings.itime

    input:
    file genes
    file pcor_net from pcor_network_raw

    output:
    file 'pcor_network.sf' into pcor_sf

    script:
    if (targets == '')
    {
      """
      ${params.pcor_settings.preamble}
      seidr import -F lm -A -u -r -z -n ${params.pcor_settings.importname} \
                   -i ${pcor_net} -g ${genes} -o pcor_network.sf \
                   -O ${params.pcor_settings.importcores}
      ${params.pcor_settings.epilog}
      """
    }
    else
    {
      """
      ${params.pcor_settings.preamble}
      seidr import -F el -A -u -r -z -n ${params.pcor_settings.importname} \
                   -i ${pcor_net} -g ${genes} -o pcor_network.sf \
                   -O ${params.pcor_settings.importcores}
      ${params.pcor_settings.epilog}
      """
    }
  }
} else { 
  Channel.empty().set {pcor_sf} 
}

if (params.narromi) {

  process narromi {
    errorStrategy params.error_strategy
    cpus params.narromi_settings.cores
    memory params.narromi_settings.memory
    publishDir params.out + '/networks/narromi'
    clusterOptions '-n ' + params.narromi_settings.tasks + ' -A ' + params.slurm_account
    time params.narromi_settings.ptime

    input:
    file expr
    file genes
    val targets

    output:
    file 'narromi_network.tsv' into narromi_network_raw

    script:
    if (params.executor == 'slurm')
    {
      """
      ${params.narromi_settings.preamble}
      srun \
      narromi -o narromi_network.tsv \
      -B $(expr $(expr ${params.ngenes} '/' ${ params.narromi_settings.tasks}) '+' 1)
         -a ${params.narromi_settings.alpha} \
         -m ${params.narromi_settings.method} \
         -i ${expr} -g ${genes} \
         -O ${params.narromi_settings.tasks} ${targets}
      ${params.narromi_settings.epilog}
      """
    }
    else
    {
      """
      ${params.narromi_settings.preamble}
      mpirun -np ${params.narromi_settings.tasks} \
      narromi -o narromi_network.tsv \
         -a ${params.narromi_settings.alpha} \
         -m ${params.narromi_settings.method} \
         -i ${expr} -g ${genes} \
         -O ${params.narromi_settings.tasks} ${targets}
      ${params.narromi_settings.epilog}
      """
    }
  }

  process narromi_import {
    errorStrategy params.error_strategy
    validExitStatus 0,3
    cpus params.narromi_settings.importcores
    memory params.narromi_settings.importmem
    publishDir params.out + '/networks/narromi'
    clusterOptions '-n 1 -A ' + params.slurm_account
    time params.narromi_settings.itime

    input:
    file genes
    file narromi_net from narromi_network_raw

    output:
    file 'narromi_network.sf' into narromi_sf

    script:
    if (targets == '')
    {
      """
      ${params.narromi_settings.preamble}
      seidr import -F m -r -z -n ${params.narromi_settings.importname} \
                   -i ${narromi_net} -g ${genes} -o narromi_network.sf \
                   -O ${params.narromi_settings.importcores}
      ${params.narromi_settings.epilog}
      """
    }
    else
    {
      """
      ${params.narromi_settings.preamble}
      seidr import -F el -r -z -n ${params.narromi_settings.importname} \
                   -i ${narromi_net} -g ${genes} -o narromi_network.sf \
                   -O ${params.narromi_settings.importcores}
      ${params.narromi_settings.epilog}
      """
    }
  }
} else { 
  Channel.empty().set {narromi_sf} 
}

if (params.tigress) {

  process tigress {
    errorStrategy params.error_strategy
    cpus params.tigress_settings.cores
    memory params.tigress_settings.memory
    publishDir params.out + '/networks/tigress'
    clusterOptions '-n ' + params.tigress_settings.tasks + ' -A ' + params.slurm_account
    time params.tigress_settings.ptime

    input:
    file expr
    file genes
    val targets

    output:
    file 'tigress_network.tsv' into tigress_network_raw

    script:
    if (params.executor == 'slurm')
    {
      """
      ${params.tigress_settings.preamble}
      srun \
      tigress -o tigress_network.tsv \
      -B $(expr $(expr ${params.ngenes} '/' ${ params.genie3_settings.tasks}) '+' 1) \
         -n ${params.tigress_settings.nlambda} \
         -l ${params.tigress_settings.min_lambda} \
         ${params.tigress_settings.scale} \
         -i ${expr} -g ${genes} \
         -O ${params.tigress_settings.cores} ${targets}
      ${params.tigress_settings.epilog}
      """
    }
    else
    {
      """
      ${params.tigress_settings.preamble}
      mpirun -np ${params.tigress_settings.tasks} \
      tigress -o tigress_network.tsv \
         -n ${params.tigress_settings.nlambda} \
         -l ${params.tigress_settings.min_lambda} \
         ${params.tigress_settings.scale} \
         -i ${expr} -g ${genes} \
         -O ${params.tigress_settings.cores} ${targets}
      ${params.tigress_settings.epilog}
      """ 
    }
  }

  process tigress_import {
    errorStrategy params.error_strategy
    validExitStatus 0,3
    cpus params.tigress_settings.importcores
    memory params.tigress_settings.importmem
    publishDir params.out + '/networks/tigress'
    clusterOptions '-n 1 -A ' + params.slurm_account
    time params.tigress_settings.itime

    input:
    file genes
    file tigress_net from tigress_network_raw

    output:
    file 'tigress_network.sf' into tigress_sf

    script:
    if (targets == '')
    {
      """
      ${params.tigress_settings.preamble}
      seidr import -F m -r -z -n ${params.tigress_settings.importname} \
                   -i ${tigress_net} -g ${genes} -o tigress_network.sf \
                   -O ${params.tigress_settings.importcores}
      ${params.tigress_settings.epilog}
      """
    }
    else
    {
      """
      ${params.tigress_settings.preamble}
      seidr import -F el -r -z -n ${params.tigress_settings.importname} \
                   -i ${tigress_net} -g ${genes} -o tigress_network.sf \
                   -O ${params.tigress_settings.importcores}
      ${params.tigress_settings.epilog}
      """
    }
  }
} else { 
  Channel.empty().set {tigress_sf} 
}

if (params.genie3) {

  process genie3 {
    errorStrategy params.error_strategy
    cpus params.genie3_settings.cores
    memory params.genie3_settings.memory
    publishDir params.out + '/networks/genie3'
    clusterOptions '-n ' + params.genie3_settings.tasks + ' -A ' + params.slurm_account
    time params.genie3_settings.ptime

    input:
    file expr
    file genes
    val targets

    output:
    file 'genie3_network.tsv' into genie3_network_raw

    script:
    if (params.executor == 'slurm')
    {
      """
      ${params.genie3_settings.preamble}
      srun \
      genie3 -o genie3_network.tsv \
      -B $(expr $(expr ${params.ngenes} '/' ${ params.genie3_settings.tasks}) '+' 1) \
         -p ${params.genie3_settings.min_prop} \
         -a ${params.genie3_settings.alpha} \
         -N ${params.genie3_settings.min_node_size} \
         # -m ${params.genie3_settings.mtry} \
         -n ${params.genie3_settings.ntree} \
         ${params.genie3_settings.scale} \
         -i ${expr} -g ${genes} \
         -O ${params.genie3_settings.cores} ${targets}
      ${params.genie3_settings.epilog}
      """
    }
    else
    {
      """
      ${params.genie3_settings.preamble}
      mpirun -np ${params.genie3_settings.tasks} \
      genie3 -o genie3_network.tsv \
         -p ${params.genie3_settings.min_prop} \
         -a ${params.genie3_settings.alpha} \
         -N ${params.genie3_settings.min_node_size} \
         -m ${params.genie3_settings.mtry} \
         -n ${params.genie3_settings.ntree} \
         ${params.genie3_settings.scale} \
         -i ${expr} -g ${genes} \
         -O ${params.genie3_settings.cores} ${targets}
      ${params.genie3_settings.epilog}
      """
    }
  }

  process genie3_import {
    errorStrategy params.error_strategy
    validExitStatus 0,3
    cpus params.genie3_settings.importcores
    memory params.genie3_settings.importmem
    publishDir params.out + '/networks/genie3'
    clusterOptions '-n 1 -A ' + params.slurm_account
    time params.genie3_settings.itime

    input:
    file genes
    file genie3_net from genie3_network_raw

    output:
    file 'genie3_network.sf' into genie3_sf

    script:
    if (targets == '')
    {
      """
      ${params.genie3_settings.preamble}
      seidr import -F m -r -z -n ${params.genie3_settings.importname} \
                   -i ${genie3_net} -g ${genes} -o genie3_network.sf \
                   -O ${params.genie3_settings.importcores}
      ${params.genie3_settings.epilog}
      """
    }
    else
    {
      """
      ${params.genie3_settings.preamble}
      seidr import -F el -r -z -n ${params.genie3_settings.importname} \
                   -i ${genie3_net} -g ${genes} -o genie3_network.sf \
                   -O ${params.genie3_settings.importcores}
      ${params.genie3_settings.epilog}
      """
    }
  }
} else { 
  Channel.empty().set {genie3_sf} 
}

if (params.plsnet) {

  process plsnet {
    errorStrategy params.error_strategy
    cpus params.plsnet_settings.cores
    memory params.plsnet_settings.memory
    publishDir params.out + '/networks/plsnet'
    clusterOptions '-n ' + params.plsnet_settings.tasks + ' -A ' + params.slurm_account
    time params.plsnet_settings.ptime

    input:
    file expr
    file genes
    val targets

    output:
    file 'plsnet_network.tsv' into plsnet_network_raw

    script:
    if (params.executor == 'slurm')
    {
      """
      ${params.plsnet_settings.preamble}
      srun \
      plsnet -o plsnet_network.tsv \
      -B $(expr $(expr ${params.ngenes} '/' ${ params.genie3_settings.tasks}) '+' 1) \
         -c ${params.plsnet_settings.components} \
       #  -p ${params.plsnet_settings.predictors} \
         -e ${params.plsnet_settings.ensemble} \
         ${params.plsnet_settings.scale} \
         -i ${expr} -g ${genes} \
         -O ${params.plsnet_settings.cores} ${targets}
      ${params.plsnet_settings.epilog}
      """
    }
    else
    {
      """
      ${params.plsnet_settings.preamble}
      mpirun -np ${params.plsnet_settings.tasks} \
      plsnet -o plsnet_network.tsv \
         -c ${params.plsnet_settings.components} \
         -p ${params.plsnet_settings.predictors} \
         -e ${params.plsnet_settings.ensemble} \
         ${params.plsnet_settings.scale} \
         -i ${expr} -g ${genes} \
         -O ${params.plsnet_settings.cores} ${targets}
      ${params.plsnet_settings.epilog}
      """
    }
  }

  process plsnet_import {
    errorStrategy params.error_strategy
    validExitStatus 0,3
    cpus params.plsnet_settings.importcores
    memory params.plsnet_settings.importmem
    publishDir params.out + '/networks/plsnet'
    clusterOptions '-n 1 -A ' + params.slurm_account
    time params.plsnet_settings.itime
    
    input:
    file genes
    file plsnet_net from plsnet_network_raw

    output:
    file 'plsnet_network.sf' into plsnet_sf

    script:
    if (targets == '')
    {
      """
      ${params.plsnet_settings.preamble}
      seidr import -F m -r -z -n ${params.plsnet_settings.importname} \
                   -i ${plsnet_net} -g ${genes} -o plsnet_network.sf \
                   -O ${params.plsnet_settings.importcores}
      ${params.plsnet_settings.epilog}
      """
    }
    else
    {
      """
      ${params.plsnet_settings.preamble}
      seidr import -F el -r -z -n ${params.plsnet_settings.importname} \
                   -i ${plsnet_net} -g ${genes} -o plsnet_network.sf \
                   -O ${params.plsnet_settings.importcores}
      ${params.plsnet_settings.epilog}
      """
    }
  }
} else { 
  Channel.empty().set {plsnet_sf} 
}

if (params.tomsimilarity) {

  process tomsimilarity {
    errorStrategy params.error_strategy
    cpus params.tomsimilarity_settings.cores
    memory params.tomsimilarity_settings.memory
    publishDir params.out + '/networks/tomsimilarity'
    clusterOptions '-n ' + params.tomsimilarity_settings.tasks + ' -A ' + params.slurm_account
    time params.tomsimilarity_settings.ptime

    input:
    file expr
    file genes
    val targets

    output:
    file 'tomsimilarity_network.tsv' into tomsimilarity_network_raw

    script:
    if (params.executor == 'slurm')
    {
      """
      ${params.tomsimilarity_settings.preamble}
      srun \
      tomsimilarity -o tomsimilarity_network.tsv \
         -m ${params.tomsimilarity_settings.method} \
         -b ${params.tomsimilarity_settings.sft} \
         -M ${params.tomsimilarity_settings.max_power} \
         --sft-cutoff ${params.tomsimilarity_settings.sftcutoff} \
         -T ${params.tomsimilarity_settings.tomtype} \
         ${params.tomsimilarity_settings.scale} \
         ${params.tomsimilarity_settings.absolute} \
         -i ${expr} -g ${genes} \
          ${targets}
      ${params.tomsimilarity_settings.epilog}
      """
    }
    else
    {
      """
      ${params.tomsimilarity_settings.preamble}
      mpirun -np ${params.tomsimilarity_settings.tasks} \
      tomsimilarity -o tomsimilarity_network.tsv \
         -m ${params.tomsimilarity_settings.method} \
         -b ${params.tomsimilarity_settings.sft} \
         -M ${params.tomsimilarity_settings.max_power} \
         --sft-cutoff ${params.tomsimilarity_settings.sftcutoff} \
         -T ${params.tomsimilarity_settings.tomtype} \
         ${params.tomsimilarity_settings.scale} \
         ${params.tomsimilarity_settings.absolute} \
         -i ${expr} -g ${genes} \
          ${targets}
      ${params.tomsimilarity_settings.epilog}
      """
    }
  }

  process tomsimilarity_import {
    errorStrategy params.error_strategy
    validExitStatus 0,3
    cpus params.tomsimilarity_settings.importcores
    memory params.tomsimilarity_settings.importmem
    publishDir params.out + '/networks/tomsimilarity'
    clusterOptions '-n 1 -A ' + params.slurm_account
    time params.tomsimilarity_settings.itime
    
    input:
    file genes
    file tomsimilarity_net from tomsimilarity_network_raw

    output:
    file 'tomsimilarity_network.sf' into tomsimilarity_sf

    script:
    if (targets == '')
    {
      """
      ${params.tomsimilarity_settings.preamble}
      seidr import -F lm -A -z -r -u -n ${params.tomsimilarity_settings.importname} \
                   -i ${tomsimilarity_net} -g ${genes} -o tomsimilarity_network.sf \
                   -O ${params.tomsimilarity_settings.importcores}
      ${params.tomsimilarity_settings.epilog}
      """
    }
    else
    {
      """
      ${params.tomsimilarity_settings.preamble}
      seidr import -F el -A -z -r -u -n ${params.tomsimilarity_settings.importname} \
                   -i ${tomsimilarity_net} -g ${genes} -o tomsimilarity_network.sf \
                   -O ${params.tomsimilarity_settings.importcores}
      ${params.tomsimilarity_settings.epilog}
      """
    }
  }
} else { 
  Channel.empty().set {tomsimilarity_sf} 
}

if (params.aggregate) {

  process aggregate {

    errorStrategy params.error_strategy
    validExitStatus 0,3
    cpus params.aggregate_settings.cores
    memory params.aggregate_settings.importmem
    publishDir params.out + '/aggregated'
    clusterOptions '-n 1 -A ' + params.slurm_account
    time params.aggregate_settings.itime

    input:

    file mi_file from mi_sf.collect().ifEmpty([])
    file clr_file from clr_sf.collect().ifEmpty([])
    file aracne_file from aracne_sf.collect().ifEmpty([])
    file anova_file from anova_sf.collect().ifEmpty([])
    file tomsimilarity_file from tomsimilarity_sf.collect().ifEmpty([])
    file pearson_file from pearson_sf.collect().ifEmpty([])
    file elnet_file from elnet_sf.collect().ifEmpty([])
    file spearman_file from spearman_sf.collect().ifEmpty([])
    file svm_file from svm_sf.collect().ifEmpty([])
    file llr_file from llr_sf.collect().ifEmpty([])
    file pcor_file from pcor_sf.collect().ifEmpty([])
    file narromi_file from narromi_sf.collect().ifEmpty([])
    file tigress_file from tigress_sf.collect().ifEmpty([])
    file genie3_file from genie3_sf.collect().ifEmpty([])
    file plsnet_file from plsnet_sf.collect().ifEmpty([])

    output:
    file 'aggregated.sf' into aggregated_sf

    script:
      """
      ${params.aggregate_settings.preamble}
      srun \
      seidr aggregate -o aggregated.sf \
         -m ${params.aggregate_settings.method} \
         ${mi_file} \
         ${clr_file} \
         ${aracne_file} \
         ${anova_file} \
         ${tomsimilarity_file} \
         ${pearson_file} \
         ${elnet_file} \
         ${spearman_file} \
         ${svm_file} \
         ${llr_file} \
         ${pcor_file} \
         ${narromi_file} \
         ${tigress_file} \
         ${genie3_file} \
         ${plsnet_file} \
         ${params.aggregate_settings.keep} \
         -O ${params.aggregate_settings.cores}
      ${params.aggregate_settings.epilog}
      """
  }
}

process HardThreshold {
  errorStrategy params.error_strategy
  validExitStatus 0,3
  publishDir params.out + '/aggregated'
  clusterOptions '-n 16 -t 12:00:00 --mem 96G -A ' + params.slurm_account

  input:
  file aggregated_file from aggregated_sf.collect()

  output:
  file 'thresholded.tsv' into thresholded_tsv

  script:
    """
    srun \
    seidr threshold -f --in-file ${aggregated_file} \
    -m 0.1 -M 0.9 -O 16 -n 10000 -o . thresholded.tsv
    """
  }
