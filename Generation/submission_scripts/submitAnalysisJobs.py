#!/usr/bin/env python
import os, sys, subprocess
from optparse import OptionParser

def main():
    parser = OptionParser()

    parser.add_option ('-i','--input', help='input directory containing delphes trees (most likely on eos)',
                       dest='input',
                       default='/eos/experiment/fcc/hh/generation/DelphesStandalone/Test/')

    parser.add_option ('-o','--output', help='output directory (default will have same name as input dir)',
                       dest='output',
                       default='histos/Test')

    parser.add_option ('-q','--queue', help='lsf queue',
                       dest='queue',
                       default='1nh')

    parser.add_option("-c","--collect", help="collects jobs for given output name",
		      dest="collect", action="store_true", 
		      default=False)

    (options, args) = parser.parse_args()
    input_dir     = options.input
    output_dir    = options.output
    queue         = options.queue
    collect       = options.collect

    if collect:
       print 'Collecting jobs for process: '+output_dir
       hadd_dir = output_dir +'/out/'
       basename = os.path.basename(output_dir)
       outfile = hadd_dir + basename + '.root'
       hadd_files = hadd_dir + basename + '_*.root'
       cmd ='hadd -f '+ outfile + ' ' + hadd_files
       os.system(cmd)
       cmd = 'rm -rf '+hadd_files+' '+output_dir+'/std'
       os.system(cmd)
       sys.exit('Collection of jobs done.')

    # first create output dir
    if not os.path.exists(output_dir):
       os.makedirs(output_dir)
       os.makedirs(output_dir+'/std/')
       os.makedirs(output_dir+'/out/')
    else:
       sys.exit('Output dir: "'+output_dir+'" exists.')

    # find list of input files
    cmd = "find {} -type f -name '*.root'".format(options.input)
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
    lst = process.communicate()[0]
    list_of_files = lst.splitlines()

    # just send one job per ntup file
    njobs = len(list_of_files)

    for job in xrange(njobs):
       
       basename = os.path.basename(output_dir) + '_'+str(job)
       currentDir = os.getcwd()
       inputFile = list_of_files[job]
       outputFile = output_dir+'/out/'+basename+'.root'
       cmd = 'bsub -o '+output_dir+'/std/'+basename +'.out -e '+output_dir+'/std/'+basename +'.err -q '+queue
       cmd +=' -J '+basename+' "submitAnalysisJob.sh '+inputFile+' '+outputFile+'"'
       
       print cmd
       # submitting jobs
       os.system(cmd)


#_______________________________________________________________________________________
if __name__ == "__main__":
    main()
