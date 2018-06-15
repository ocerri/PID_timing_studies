#!/usr/bin/env python
import os, sys, subprocess
import argparse
import commands
import time
#____________________________________________________________________________________________________________
### processing the external os commands
def processCmd(cmd, quite = 0):
    status, output = commands.getstatusoutput(cmd)
    if (status !=0 and not quite):
        print 'Error in processing command:\n   ['+cmd+']'
        print 'Output:\n   ['+output+'] \n'
    return output


#_____________________________________________________________________________________________________________
def main():
    parser = argparse.ArgumentParser()

    parser.add_argument ('--pt_min', help='ptmin',  default='0.')
    parser.add_argument ('--pt_max', help='ptmax',  default='-1')
    parser.add_argument ('--pycard', help='rel. path to example/Pythia8 directory', default='pp2HardQCD.pythia')
    parser.add_argument ('--dcard', help='delphes card (rel. path to card dir)', default='CMS_PhaseII_trackNtime_lightTree.tcl')
    parser.add_argument ('--outdir', help='output directory ', default='/afs/cern.ch/user/o/ocerri/cernbox/PID_timing_studies/_root/jobs/pp2HardQCD_PU100')
    parser.add_argument ('--force_production', action='store_true', default=False, help='Proceed even if the directory is already existing')
    parser.add_argument ('--njobs', help='number of jobs', default=2)
    parser.add_argument ('--nev', help='number of events per job', default=100)
    parser.add_argument ('--queue', help='lsf queue', default='8nh')
    parser.add_argument ('--memory', help='min virtual memory', default='4000')
    parser.add_argument ('--disk', help='min disk space', default='2000')
    parser.add_argument ('--stopmass', help='STOPMAss in GeV', default='1000')
    parser.add_argument ('--st_seed', help='starting seed', default=1, type=int)

    args = parser.parse_args()

    pt_min        = args.pt_min
    pt_max        = args.pt_max
    pycard        = args.pycard
    dcard         = args.dcard
    outdir        = args.outdir
    njobs         = int(args.njobs)
    nev           = args.nev
    queue         = args.queue
    mem           = args.memory
    disk          = args.disk
    stopmass      = args.stopmass

    if not os.path.exists(outdir):
        os.makedirs(outdir)
        os.makedirs(outdir+'/std/')
        os.makedirs(outdir+'/cfg/')
    elif not args.force_production:
        sys.exit('Output dir: "'+outdir+'" exists.')

    os.system('chmod +x submitDelphesJob.sh')
    print '[Submitting jobs]'
    jobCount=0
    for job in xrange(njobs):

        print 'Submitting job '+str(job+1)+' / '+str(njobs)
        job += args.st_seed
        basename = os.path.basename(outdir) + '_'+str(job)
        outputFile = outdir+'/'+basename
        seed=str(job)

        cmd = 'bsub -oo '+outdir+'/std/'+basename +'.out -eo '+outdir+'/std/'+basename +'.err -q '+queue
        cmd += ' -R "rusage[mem={}:pool={}]"'.format(mem,disk)
        cmd += ' -J '+basename
        cmd += ' /afs/cern.ch/user/o/ocerri/cernbox/PID_timing_studies/Generation/submission_scripts/submitDelphesJob.sh'
        cmd += ' '+dcard+' '+pycard+' '+pt_min+' '+pt_max+' '+nev+' '+seed+' '+outputFile+' '+str(stopmass)

        print cmd

        # submitting jobs

        output = processCmd(cmd)
        kk = 0
        while (kk<5 and ('error' in output)):
            kk += 1
            time.sleep(1.0);
            output = processCmd(cmd)
            if ('error' not in output):
                print 'Submitted after retry - job '+str(jobCount+1)
            else:
                print output

        jobCount += 1
        time.sleep(0.01);

#_______________________________________________________________________________________
if __name__ == "__main__":
    main()
