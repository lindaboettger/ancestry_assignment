# REQUIRES `use Python-2.7`

print("I am in the python script")

# 1. have python submit the job to qstat
# 2. have python check every X seconds to see if the job is still in the queue
# 3. when the job is no longer in the queue, check for output
# 4a. if output, then end the python 
# 4b. if no output, then increment the failed counter and resubmit same job
# 5. if fail counter reaches some number, then write a file that says 'this failed X times, aborting' and quit.

import subprocess
import sys
import time
import os.path
import os

print(sys.argv)

win = sys.argv[1]
query = sys.argv[2]
out = sys.argv[3]
chr_num = sys.argv[4]
script = sys.argv[5]

sleep_sec = 60
failures = 0
max_failures = 5

err_file = '/seq/vgb/linda/mutts/err/mutt_errSM{chr_num}'.format(chr_num=chr_num)
out_file = '/seq/vgb/linda/mutts/out/mutt_outSM{chr_num}'.format(chr_num=chr_num)
if os.path.isfile(err_file):
    os.remove(err_file)
if os.path.isfile(out_file):
    os.remove(out_file)

command_string = '''
    qsub -S /bin/bash -v chr={chr_num},win={win},query={query},out={out}
       -o /seq/vgb/linda/mutts/out/mutt_outSM{chr_num}
       -e /seq/vgb/linda/mutts/err/mutt_errSM{chr_num}
       -N SM_chr_{chr_num}
       -p -10
       -wd /seq/vgb/linda/mutts/err
       {script}
'''.format(win=win, query=query, out=out, chr_num=chr_num, script=script)
print(command_string)

# look for this job name in the qstat output
job_name = 'SM_chr_'+str(chr_num)

# check for this file (full path)
outfile_name = '{out}_{chr_num}.cfg'.format(
        out=out, chr_num=chr_num)

# this will be our error file
err_file = outfile_name + '.failed'

# start first attempt
qsub_output = subprocess.check_output(command_string.split())

while True:

    time.sleep(60)
    
    qstat_output = subprocess.check_output('qstat -r'.split())

    # check if job is running/queued (in qstat output) still,
    # if not, it either failed or finished
    if job_name not in qstat_output:
        # check if it finished, outfile_name will exist
        if os.path.isfile(outfile_name):
            # job finished, so...
            quit()
        else:
            # add one to the failures counter
            failures += 1
            if failures < max_failures:
                # re-run job
                qsub_output = subprocess.check_output(command_string.split())
            else:
                # if we've had max_failures failures,
                # then write an error message file and die
                with open(err_file, 'w') as err_filehandle:
                    err_filehandle.write(('Chr {chr_num} failed. '+
                                '{outfile_name} does not exist.').format(
                            chr_num=chr_num, outfile_name=outfile_name))
                quit()










