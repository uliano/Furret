import os
import subprocess
from dataclasses import dataclass


@dataclass
class MemeJob:
    fasta_file: str
    destination_directory: str
    meme_executable: str


def process_meme(job: MemeJob) -> bool:
    cwd = os.getcwd()
    os.chdir(job.destination_directory)
    cmd = f'{job.meme_executable} "{job.fasta_file}" -protein -oc . -nostatus -time 18000 -mod zoops -nmotifs 100 ' \
          f'-minw 6 -maxw 50 -objfun classic -markov_order 0 '
    print(f'Processing Meme for {os.path.basename(job.fasta_file)}')
    subprocess.call(cmd, shell=True)
    os.chdir(cwd)
    return True
