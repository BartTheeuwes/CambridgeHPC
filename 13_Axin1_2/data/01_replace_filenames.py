import glob
import os
import re

#' Works only in the current folder
filenames = glob.glob('*fq.gz')

def renamer(filename, dry_run = True, change_extension = None):
    pattern = '(SLX-\d{5})\.(SI[a-zA-Z0-9]{4})\.(.*)\.(s_\d*)\.([ri]_\d)\.(.*)'
    slx = re.sub(pattern, '\\1', filename)

    ind = re.sub(pattern, '\\2', filename)

    lane = re.sub(pattern, '\\4', filename)
    lane = re.sub('(s_)(.*)', '\\2', lane)
    lane = 'L'+lane.zfill(3)

    read = re.sub(pattern, '\\5', filename)
    read = re.sub('r_', 'R', read)
    read = re.sub('i_', 'I', read)

    if change_extension is None:
        extension = re.sub(pattern, '\\6', filename)
    else: extension = change_extension
   
    renamed = f'{ind}_S1_{lane}_{read}_001.{extension}'
    print(f'Renaming: {filename} to {renamed}')
   
    if dry_run == False:
        os.rename(filename, renamed)
      #     return(renamed)
      #     print([slx, ind, lane, read, extension])

for i in filenames:
         renamer(i, dry_run = False, change_extension = 'fastq.gz')
