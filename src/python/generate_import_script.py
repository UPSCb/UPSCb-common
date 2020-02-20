#!/usr/bin/env python3
import argparse
import re
import shutil
import os

"""
Bunch of defaults for some algorithms
"""
aracne = {
  'r': True,
  'z': True,
  'u': True,
  'a': False,
  'F': 'lm',
  'n': 'ARACNE',
  'o': 'aracne.sf'
}

clr = {
  'r': True,
  'z': True,
  'u': True,
  'a': False,
  'F': 'lm',
  'n': 'CLR',
  'o': 'clr.sf'
}

elnet = {
  'r': True,
  'z': True,
  'u': False,
  'a': False,
  'F': 'm',
  'n': 'ELNET',
  'o': 'elnet.sf'
}

genie3 = {
  'r': True,
  'z': True,
  'u': False,
  'a': False,
  'F': 'm',
  'n': 'GENIE3',
  'o': 'genie3.sf'
}

llr = {
  'r': True,
  'z': True,
  'u': False,
  'a': False,
  'F': 'm',
  'n': 'LLR',
  'o': 'llr.sf'
}

mi = {
  'r': True,
  'z': True,
  'u': True,
  'a': False,
  'F': 'lm',
  'n': 'MI',
  'o': 'mi.sf'
}

narromi = {
  'r': True,
  'z': True,
  'u': False,
  'a': False,
  'F': 'm',
  'n': 'NARROMI',
  'o': 'narromi.sf'
}

svm = {
  'r': True,
  'z': True,
  'u': False,
  'a': False,
  'F': 'm',
  'n': 'SVM',
  'o': 'svm.sf'
}

pcor = {
  'r': True,
  'z': True,
  'u': True,
  'a': True,
  'F': 'lm',
  'n': 'PCOR',
  'o': 'pcor.sf'
}

pearson = {
  'r': True,
  'z': True,
  'u': True,
  'a': True,
  'F': 'lm',
  'n': 'PEARSON',
  'o': 'pearson.sf'
}

spearman = {
  'r': True,
  'z': True,
  'u': True,
  'a': True,
  'F': 'lm',
  'n': 'SPEARMAN',
  'o': 'spearman.sf'
}

plsnet = {
  'r': True,
  'z': True,
  'u': False,
  'a': False,
  'F': 'm',
  'n': 'PLSNET',
  'o': 'plsnet.sf'
}

tigress = {
  'r': True,
  'z': True,
  'u': False,
  'a': False,
  'F': 'm',
  'n': 'TIGRESS',
  'o': 'tigress.sf'
}

class inputCtl(object):
  """docstring for inputCtl"""
  def __init__(self, args):
    super(inputCtl, self).__init__()
    self.args = args

    if not os.path.exists(args.input):
      raise FileNotFoundError(args.input)

    if not os.path.exists(args.genes):
      raise FileNotFoundError(args.genes)

    if not shutil.which('seidr') and not args.exec:
      raise RuntimeError('No executable for seidr found. Specify one with ' +
                         'the --exec option or add it to your PATH.')

    args.input = os.path.abspath(args.input)
    args.genes = os.path.abspath(args.genes)

  def estimateRAM(self):
    ng = 0
    for line in open(self.args.genes, 'r'):
      for field in line.split('\t'):
        ng = ng + 1

    return int(((ng ** 2) / 2 * 48) / (1024 ** 2) + 4096)


class cmdGenerator(object):
  """docstring for pearsonGenerator"""
  def __init__(self, args, defaults):
    super(cmdGenerator, self).__init__()
    self.cmd = []
    if args.reverse:
      self.cmd.append('-r')
    elif (not args.no_reverse) and (defaults['r']):
      self.cmd.append('-r')
    if args.drop_zero:
      self.cmd.append('-z')
    elif (not args.no_drop_zero) and (defaults['z']):
      self.cmd.append('-z')
    if args.undirected:
      self.cmd.append('-u')
    elif (not args.no_undirected) and (defaults['u']):
      self.cmd.append('-u')
    if args.absolute:
      self.cmd.append('-A')
    elif (not args.no_absolute) and (defaults['a']):
      self.cmd.append('-A')
    if args.format:
      self.cmd.append('-F')
      self.cmd.append(args.format)
    else:
      self.cmd.append('-F')
      self.cmd.append(defaults['F'])
    if args.name:
      self.cmd.append('-n')
      self.cmd.append(args.name)
    else:
      self.cmd.append('-n')
      self.cmd.append(defaults['n'])
    if args.output:
      self.cmd.append('-o')
      self.cmd.append(args.output)
    else:
      self.cmd.append('-o')
      self.cmd.append(os.path.join(os.path.dirname(args.input), defaults['o']))

  def get(self):
    return self.cmd



class scriptGenerator(object):
  """
  Class that generates a bash script for SLURM submission of seidr import
  based of file names.
  """
  def __init__(self, args):
    super(scriptGenerator, self).__init__()

    ictl = inputCtl(args)

    self.cmd = [args.exec, 'import']

    print(args.shebang)

#    print('#SBATCH --mem=', ictl.estimateRAM(), 'M', sep='')
    print('#SBATCH -n 1')
    print('#SBATCH -c', args.ncpus)

    if args.slurmctl[0]:
      for line in args.slurmctl[0]:
        print(line)

    print('\n##############################################################')
    print('########## PRE IMPORT COMMANDS ###############################')
    print('##############################################################\n')

    if args.pre:
      for line in args.pre:
        print(' '.join(line))

    print('\n##############################################################')
    print('########## ALWAYS SET OMP_NUM_THREADS ########################')
    print('##############################################################\n')

    print('export OMP_NUM_THREADS=1')

    print('\n##############################################################')
    print('########## MAIN IMPORT CMD ###################################')
    print('##############################################################\n')
 
    if args.algorithm == 'auto':
      if re.search('aracne', args.input, flags=re.IGNORECASE):
        gen = cmdGenerator(args, aracne)
        for arg in gen.get():
          self.cmd.append(arg)
      elif re.search('clr', args.input, flags=re.IGNORECASE):
        gen = cmdGenerator(args, clr)
        for arg in gen.get():
          self.cmd.append(arg)
      elif re.search('elnet', args.input, flags=re.IGNORECASE):
        gen = cmdGenerator(args, elnet)
        for arg in gen.get():
          self.cmd.append(arg)
      elif re.search('genie3', args.input, flags=re.IGNORECASE):
        gen = cmdGenerator(args, genie3)
        for arg in gen.get():
          self.cmd.append(arg)
      elif re.search('llr', args.input, flags=re.IGNORECASE):
        gen = cmdGenerator(args, llr)
        for arg in gen.get():
          self.cmd.append(arg)          
      elif re.search('narromi', args.input, flags=re.IGNORECASE):
        gen = cmdGenerator(args, narromi)
        for arg in gen.get():
          self.cmd.append(arg)
      elif re.search('mi', args.input, flags=re.IGNORECASE):
        gen = cmdGenerator(args, mi)
        for arg in gen.get():
          self.cmd.append(arg)    
      elif re.search('svm', args.input, flags=re.IGNORECASE):
        gen = cmdGenerator(args, svm)
        for arg in gen.get():
          self.cmd.append(arg)
      elif re.search('pearson', args.input, flags=re.IGNORECASE):
        gen = cmdGenerator(args, pearson)
        for arg in gen.get():
          self.cmd.append(arg)
      elif re.search('spearman', args.input, flags=re.IGNORECASE):
        gen = cmdGenerator(args, spearman)
        for arg in gen.get():
          self.cmd.append(arg)
      elif re.search('pcor', args.input, flags=re.IGNORECASE):
        gen = cmdGenerator(args, pcor)
        for arg in gen.get():
          self.cmd.append(arg)
      elif re.search('plsnet', args.input, flags=re.IGNORECASE):
        gen = cmdGenerator(args, plsnet)
        for arg in gen.get():
          self.cmd.append(arg)
      elif re.search('tigress', args.input, flags=re.IGNORECASE):
        gen = cmdGenerator(args, tigress)
        for arg in gen.get():
          self.cmd.append(arg)
      else:
        raise ValueError("Couldn't autodetect algorithm. Rerun with --algorithm.")
    else:
      ct = ['aracne', 'clr', 'elnet', 'genie3', 'llr', 'mi', 'narromi', 'svm',
            'pearson', 'spearman', 'pcor', 'plsnet', 'tigress']
      if not args.algorithm in ct:
        raise ValueError("Only the following arguments are supported with --algorithm: " +
                         ",".join(ct))
      if args.algorithm == 'aracne':
        gen = cmdGenerator(args, aracne)
        for arg in gen.get():
          self.cmd.append(arg)
      if args.algorithm == 'clr':
        gen = cmdGenerator(args, clr)
        for arg in gen.get():
          self.cmd.append(arg)
      if args.algorithm == 'elnet':
        gen = cmdGenerator(args, elnet)
        for arg in gen.get():
          self.cmd.append(arg)
      if args.algorithm == 'genie3':
        gen = cmdGenerator(args, genie3)
        for arg in gen.get():
          self.cmd.append(arg)
      if args.algorithm == 'llr':
        gen = cmdGenerator(args, llr)
        for arg in gen.get():
          self.cmd.append(arg)
      if args.algorithm == 'mi':
        gen = cmdGenerator(args, mi)
        for arg in gen.get():
          self.cmd.append(arg)    
      if args.algorithm == 'narromi':
        gen = cmdGenerator(args, narromi)
        for arg in gen.get():
          self.cmd.append(arg)
      if args.algorithm == 'svm':
        gen = cmdGenerator(args, svm)
        for arg in gen.get():
          self.cmd.append(arg)
      if args.algorithm == 'pearson':
        gen = cmdGenerator(args, pearson)
        for arg in gen.get():
          self.cmd.append(arg)
      if args.algorithm == 'spearman':
        gen = cmdGenerator(args, spearman)
        for arg in gen.get():
          self.cmd.append(arg)
      if args.algorithm == 'pcor':
        gen = cmdGenerator(args, pcor)
        for arg in gen.get():
          self.cmd.append(arg)
      if args.algorithm == 'plsnet':
        gen = cmdGenerator(args, plsnet)
        for arg in gen.get():
          self.cmd.append(arg)
      if args.algorithm == 'tigress':
        gen = cmdGenerator(args, tigress)
        for arg in gen.get():
          self.cmd.append(arg)

    self.cmd.append('-i')
    self.cmd.append(args.input)

    self.cmd.append('-g')
    self.cmd.append(args.genes)

    self.cmd.append('-O')
    self.cmd.append(str(args.ncpus))

    print(' '.join(self.cmd))    

if __name__ == "__main__":
  parser = argparse.ArgumentParser()

  groupR = parser.add_mutually_exclusive_group()
  groupZ = parser.add_mutually_exclusive_group()
  groupU = parser.add_mutually_exclusive_group()
  groupA = parser.add_mutually_exclusive_group()
  
  groupR.add_argument("-r", "--reverse", action="store_true", 
                      help="Force 'reverse' argument")
  groupZ.add_argument("-z", "--drop-zero", action="store_true", 
                      help="Force 'drop-zero' argument")
  groupU.add_argument("-u", "--undirected", action="store_true", 
                      help="Force 'undirected' argument")
  groupA.add_argument("-a", "--absolute", action="store_true", 
                      help="Force 'absolute' argument")

  groupR.add_argument("-R", "--no-reverse", action="store_true", 
                      help="Force disable 'reverse' argument")
  groupZ.add_argument("-Z", "--no-drop-zero", action="store_true", 
                      help="Force disable 'drop-zero' argument")
  groupU.add_argument("-U", "--no-undirected", action="store_true", 
                      help="Force disable 'undirected' argument")
  groupA.add_argument("-A", "--no-absolute", action="store_true", 
                      help="Force disable 'absolute' argument")

  parser.add_argument("-n", "--name", type=str, 
                      help="Override name argument")
  parser.add_argument("-f", "--format", type=str, 
                      help="Override format argument")

  parser.add_argument('-t', '--exec', type=str, help='Seidr executable', 
                      default=shutil.which('seidr'))

  parser.add_argument('-m', '--algorithm', type=str, help='Algorithm set to use', 
                      default='auto')

  parser.add_argument('-c', '--ncpus', type=int, help='Number of sorting threads', 
                      default=1)

  parser.add_argument('-p', '--pre', type=str, 
                      help='Lines to add before the import command',
                      nargs='*', action='append')

  parser.add_argument('-s', '--shebang', type=str, help='Shebang string',
                      default='#!/bin/bash -l')

  parser.add_argument("-i", "--input", type=str, help="Input file", 
                      required=True)
  parser.add_argument("-g", "--genes", type=str, help="Input gene list", 
                      required=True)

  parser.add_argument('-o', '--output', type=str, 
                      help='Override output file name')

  parser.add_argument("slurmctl", nargs="*", type=str, action='append',
                      help="Additional SLURM directives (do not set memory or CPU)")

  args = parser.parse_args()

  scriptGenerator(args)
