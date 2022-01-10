"""
File: pipe_by_increment.py
--------------------------

"""
import os
import subprocess
import sys

if os.path.exists(sys.argv[1]):
    meta_merge_snps_dir = sys.argv[1]
    merge_snps_dirs = sorted([dir for dir in os.listdir(meta_merge_snps_dir) if not dir.startswith('.')])

    job_dir = os.path.join(os.getcwd(), 'postprocess_pipe_snps_jobs')
    if not os.path.exists(job_dir):
        os.mkdir(job_dir)

    for snps_dir in merge_snps_dirs:
        sample = snps_dir[:7]
        job_name = sample + '_pipe_snps'
        job_filename = sample + '_pipe_snps.job'
        job_file_path = os.path.join(job_dir, job_filename)     
        snps_dir_path = os.path.join(meta_merge_snps_dir, snps_dir)

        job_file = open(job_file_path, 'w')
        job_file.write('#!/bin/bash\n')
        job_file.write('#SBATCH --nodes=2\n')
        job_file.write('#SBATCH --mem=5G\n')
        job_file.write('#SBATCH --job-name=%s\n' % job_name)
        job_file.write('#SBATCH --error=sherlock_out/%s_postprocess.err\n' % job_name)
        job_file.write('#SBATCH --output=sherlock_out/%s_postprocess.out\n' % job_name)
        job_file.write('#SBATCH --time=1-18:00:00\n')
        job_file.write('#SBATCH --qos=normal\n\n')

        job_file.write('ml py-numpy/1.19.2_py36\n')
        job_file.write('ml py-pandas/1.0.3_py36\n\n')

        job_file.write('python3 $HOME/postprocess_pipeline.py ' + snps_dir_path)
        job_file.close()
        
        subprocess.call(['sbatch', '-c', '4', job_file_path])

else:
    sys.exit('Must enter proper file path')