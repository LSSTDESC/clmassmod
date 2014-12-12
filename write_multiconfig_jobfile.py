#!/usr/bin/env python
#####################
# This program creates a job file that multiconfig_nfwfit uses to run multiple fits
#  with different configurations
#####################


import sys, os, json, argparse, glob, stat

####################

## add preconfiged functions for slac, condor running

####################

def setupCondor(configs, jobdir, jobname, simdir, simfiles, input_extensions): 
    
    if not os.path.exists(jobdir):
        os.mkdir(jobdir)

    configfiles = ['{0}/{1}/config.sh'.format(simdir, config) for config in configs]

    for i, simfile in enumerate(simfiles):

        inputfiles = ['{0}{1}'.format(simfile, x) for x in input_extensions]
        
        jobparams = createJobParams(simfile,
                                    configfiles,
                                    inputfiles = inputfiles,
                                    workbase = './',
                                    stripCatExt = False)
        writeJobfile(jobparams, '{0}/{1}.{2}.job'.format(jobdir, jobname, i))

    condorfile = '''executable = /vol/braid1/vol1/dapple/mxxl/oldmxxl/nfwfit_condorwrapper.sh
universe = vanilla
Error = {jobdir}/{jobname}.$(Process).stderr
Output = {jobdir}/{jobname}.$(Process).stdout
Log = {jobdir}/{jobname}.$(Process).batch.log
Arguments = {jobdir}/{jobname}.$(Process).job
queue {njobs}
'''.format(jobdir = jobdir, jobname = jobname, njobs = len(simfiles))

    with open('{0}/{1}.submit'.format(jobdir, jobname), 'w') as output:
        output.write(condorfile)



####################

def setupCondor_MXXL(configs, jobdir, jobname, simdir):

    simdirbase = '/vol/braid1/vol1/dapple/mxxl'
    
    input_extensions = '.convergence_map .shear_1_map .shear_2_map'.split()


    for i, snap in enumerate('snap54 snap41'.split()):

        simdir = '{0}/mxxl{1}'.format(simdirbase, snap)

        simfiles = [os.path.splitext(x)[0] for x in glob.glob('{0}/halo_*convergence_map'.format(simdir))]

        setupCondor(configs, jobdir, '{0}.mxxl.{1}'.format(jobname, snap), simdir, simfiles, input_extensions)


###################

def setupCondor_BK11(configs, jobdir, jobname):

    simdirbase = '/vol/braid1/vol1/dapple/mxxl'
    
    input_extensions = ['']

    for i, snap in enumerate('snap124 snap141'.split()):

        simdir = '{0}/{1}/intlength400'.format(simdirbase, snap)

        simfiles = glob.glob('{0}/haloid*.fit'.format(simdir))

        setupCondor(configs, jobdir, '{0}.bk11.{1}'.format(jobname, snap), simdir, simfiles, input_extensions)


###################

def setupCondor_BCC(configs, jobdir, jobname):

    simdir = '/vol/braid1/vol1/dapple/mxxl/bcc'
    
    input_extensions = ['']

    simfiles = glob.glob('{0}/cluster_*.hdf5'.format(simdir))

    setupCondor(configs, jobdir, '{0}.bcc'.format(jobname), simdir, simfiles, input_extensions)


###################


def setupLSF(configs, jobdir, jobname, simdir, simfiles):

    configfiles = ['{0}/{1}/config.sh'.format(simdir, config) for config in configs]

    for simfile in simfiles:

        jobparams = createJobParams(catalogname = simfile,
                                    confignames = configfiles,
                                    inputfiles = [simfile],
                                    workbase = '/scratch/')

        jobid = '{0}.{1}'.format(jobname, jobparams['outbasename'])

        jobfile = '{0}/{1}.job'.format(jobdir, jobid)
        writeJobfile(jobparams, jobfile)

        logfile = '{0}/{1}.log'.format(jobdir, jobid)

        lsffile = '{0}/p300.{1}'.format(jobdir, jobid)
        lsfcommand = '''#!/bin/bash
bsub -q long -oo {logfile} ./multiconfig_nfwfit.py {jobfile}

'''.format(logfile = logfile, jobfile = jobfile)

        with open(lsffile, 'w') as output:
            output.write(lsfcommand)

        os.chmod(lsffile, stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR)
        


        

####################

def setupLSF_BCC(configs, jobdir, jobname):

    simdir = '/nfs/slac/g/ki/ki02/dapple/bcc_clusters/recentered'

    simfiles = glob.glob('{0}/cluster_*.hdf5'.format(simdir))

    setupLSF(configs, jobdir, '{0}.bcc'.format(jobname), simdir, simfiles)

###################

def setupLSF_BK11(configs, jobdir, jobname):

    for snap in 'snap124 snap141'.split():

        simdir = '/nfs/slac/g/ki/ki02/dapple/beckersims/{0}/intlength400'.format(snap)

        simfiles = glob.glob('{0}/haloid*.fit'.format(simdir))

        setupLSF(configs, jobdir, '{0}.bk11.{1}'.format(jobname, snap), simdir, simfiles)

    


####################

def createJobParams(catalogname, confignames, inputfiles, outputExt = '.out', workbase = None, stripCatExt = True):


    catalogbase = os.path.basename(catalogname)
    if stripCatExt:
        catalogbase, ext = os.path.splitext(catalogbase)

    workdir = None


    jobparams = dict(inputfiles = inputfiles, outputExt = outputExt, workbase = workbase, 
                    catname = catalogname, outbasename = catalogbase, 
                    configurations = confignames)

    return jobparams


####################



def writeJobfile(jobparams, jobfile):

    with open(jobfile, 'w') as output:

        json.dump(jobparams, output)

###################


#def main(argv = sys.argv):
#
#    parser = argparse.ArgumentParser()
#
#    parser.add_argument('inputname')
