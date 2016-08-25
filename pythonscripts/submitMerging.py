#!/usr/bin/env python
import os,sys
import optparse
import commands
import time
import re

#command line configuration
usage = 'usage: %prog [options]'
parser = optparse.OptionParser(usage)
parser.add_option('-q', '--queue'      ,dest='queue'  ,help='batch queue'          ,default='1nd')
parser.add_option('-i', '--inputF'      ,dest='inputF' ,help='input file list'     ,default='events_merge_test.list', type='string')
parser.add_option('-m', '--nfiles'      ,dest='nfiles'   ,help='number of files per job'       ,default=500,    type=int)
parser.add_option('-j', '--jobs'       ,dest='jobs'   ,help='number of jobs'       ,default=1,    type=int)
parser.add_option(      '--proxy'      ,dest='proxy'  ,help='proxy to be used'     ,default=None, type='string')
parser.add_option('-o', '--output'     ,dest='output' ,help='output directory'     ,default='/afs/cern.ch/work/l/lvermunt/public/Merge')

(opt, args) = parser.parse_args()

#prepare working directory
cmsswBase=os.environ['CMSSW_BASE']
workBase=os.getcwd()
jobsBase='%s/FARM%s'%(workBase,time.time())
os.system('mkdir -p %s'%jobsBase)

#init a new proxy if none has been passed
if opt.proxy is None:
    print 'Initiating a new proxy'
    os.system('voms-proxy-init --voms cms --valid 72:00')
    os.system('cp /tmp/x509up_u`id -u \`whoami\`` %s/proxyforprod' % workBase)
    print 'Production proxy is now available in %s/proxyforprod (valid for 72h)' % workBase

#loop over the required number of jobs
inputfile='%s'% (opt.inputF)

f = open(inputfile,'r')
num_lines = sum(1 for line in f)
print 'num_lines %d' % num_lines
f.close()

ff = open(inputfile,'r')
newfile='new_%s' % (opt.inputF)
newf = open(newfile,'w')
for line in ff:
    newf.write(line)


njobs = num_lines/opt.nfiles + 1
print 'njobs %d' % njobs

for n in xrange(0,njobs):
    
    outdir = '%s' % (opt.output)
                       
    #create bash script to be submitted
    scriptFile = open('%s/runJob_%d.sh'%(jobsBase,n), 'w')
    scriptFile.write('#!/bin/bash\n')
    scriptFile.write('export X509_USER_PROXY=%s/proxyforprod\n' % workBase)
    scriptFile.write('cd %s/src\n'%cmsswBase)
    scriptFile.write('eval `scram r -sh`\n')
    scriptFile.write('cd -\n')

    scriptFile.write('cp /afs/cern.ch/user/m/mverweij/public/tools/runMergingForest.C %s \n' % (jobsBase))
    scriptFile.write('cp /afs/cern.ch/user/m/mverweij/public/tools/mergeFileMerger.C %s \n' % (jobsBase))
    scriptFile.write('cp /afs/cern.ch/user/m/mverweij/public/tools/runMergingForest.C . \n')
    scriptFile.write('cp /afs/cern.ch/user/m/mverweij/public/tools/mergeFileMerger.C . \n')
    scriptFile.write('cp %s/%s .\n' % (workBase,newfile))
    scriptFile.write('root -b -q "runMergingForest.C(\\"%s\\",%d,%d,\\"HiForest_%d.root\\")" \n' % (newfile,n*opt.nfiles,opt.nfiles,n))

    scriptFile.write('cmsMkdir %s\n' % outdir)
    scriptFile.write('ls\n')
    #scriptFile.write('cmsStage -f HiForest_%d.root %s/HiForest_%d.root\n' % (n,outdir,n) )
    scriptFile.write('cp -f HiForest_%d.root %s/HiForest_%d.root\n' % (n,outdir,n) )
    scriptFile.write('rm HiForest_*root\n')
    scriptFile.close()
    
    #preare to run it
    os.system('chmod u+rwx %s/runJob_%d.sh' % (jobsBase,n))
    
    #submit it to the batch or run it locally
    if opt.queue=='':
        print 'Job #%d will run locally' % n
        os.system('%s/runJob_%d.sh' % (jobsBase,n) )
    else:
        print 'Job #%d will run remotely' % n
        os.system("bsub -q %s -R \"swp>1000 && pool>30000\" -J MERGE_%d \'%s/runJob_%d.sh\'" % (opt.queue,n,jobsBase,n) )
       
