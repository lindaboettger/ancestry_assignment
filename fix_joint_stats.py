import os
import glob
import re


pct_files = sys.argv[1]

file_list = glob.glob(pct_files)

for fn in file_list:

    new_fn = fn.replace('txt','fixed.txt')

    fh_out = open(new_fn, 'w')
    fh_in = open(fn, 'r')

    for line in fh_in:
        line = line.strip("\n")
        matchobj = re.search(
                '(^\d+\s)([a-zA-Z \-]+)(\s\d+.*$)',
                line)

        if not matchobj:
            print('ERR - Did not match: '+line)
            continue

        groups = list(matchobj.groups())
        groups[1] = groups[1].replace(' ','_')
        fh_out.write(''.join(groups)+"\n")
    fh_in.close()
    fh_out.close()

