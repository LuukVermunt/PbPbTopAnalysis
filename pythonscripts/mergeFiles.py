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
parser.add_option('-j', '--jobs'       ,dest='jobs'   ,help='number of jobs'       ,default=1,    type=int)
parser.add_option('-i', '--inputF'      ,dest='inputF' ,help='input file list'     ,default='events_merge.list', type='string')
parser.add_option(      '--proxy'      ,dest='proxy'  ,help='proxy to be used'     ,default=1, type='string')
parser.add_option('-o', '--output'     ,dest='output' ,help='output directory'     ,default='/afs/cern.ch/work/l/lvermunt/public/Data110816')

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
newfile='new_%s' % (opt.inputF)
f = open(inputfile,'r')
newf = open(newfile,'w')

jobCounter=0
lineCounter=0

scriptFile = open('%s/runJob_mergeForest_%d.sh'%(jobsBase,jobCounter), 'w')
scriptFile.write('#!/bin/bash\n')
scriptFile.write('export X509_USER_PROXY=%s/proxyforprod\n' % workBase)
scriptFile.write('cd %s/src\n'%cmsswBase)
scriptFile.write('eval `scram r -sh`\n')
scriptFile.write('cd -\n')
#scriptFile.write('hadd -f %s/HiForest_%d.root '%(opt.output, jobCounter))          

for line in f:
    if not line.find('root://eoscms//eos/cms/') :
        if not line.startswith('#') :
            if lineCounter < 50 :
                #print line
                newf.write(line.replace('root://eoscms//eos/cms/', '# root://eoscms//eos/cms/'))
                if lineCounter < 2 :
                    if lineCounter == 0:
                        scriptFile.write('hadd -f %s/HiForest_%d.root '%(opt.output, jobCounter))
                        line = line.replace('\n','')
                        scriptFile.write(line),
                        scriptFile.write(' ')
                        lineCounter+=1
                    else :
                        scriptFile.write(line)
                        lineCounter+=1
                elif lineCounter%2 != 0 :
                    scriptFile.write('hadd -f %s/HiForest_%d.root '%(opt.output, jobCounter))
                    scriptFile.write('%s/HiForest_temp_%d.root '%(opt.output, jobCounter))
                    scriptFile.write(line)
                    scriptFile.write('rm %s/HiForest_temp_%d.root \n'%(opt.output, jobCounter))
                    lineCounter+=1
                else :
                    scriptFile.write('hadd -f %s/HiForest_temp_%d.root '%(opt.output, jobCounter))
                    scriptFile.write('%s/HiForest_%d.root '%(opt.output, jobCounter))
                    scriptFile.write(line)
                    scriptFile.write('rm %s/HiForest_%d.root \n'%(opt.output, jobCounter))
                    lineCounter+=1
            else :
                newf.write(line)


scriptFile.close()

#preare to run it
os.system('chmod u+rwx %s/runJob_mergeForest_%d.sh' % (jobsBase,jobCounter))

#submit it to the batch or run it locally
#if opt.queue=='':
#    print 'Job #%d will run locally' % jobCounter
#    os.system('%s/runJob_mergeForest_%d.sh' % (jobsBase,jobCounter) )
#else:
#    print 'Job #%d will run remotely' % jobCounter
#    os.system("bsub -q %s -R \"swp>1000 && pool>30000\" -J FTT1lX_%d \'%s/runJob_mergeForest_%d.sh\'" % (opt.queue,jobCounter,jobsBase,jobCounter) )

#jobCounter+=1