import sys
import datetime
import os


def print_log(msg, outdir=None):
    ''' Adds a timestamp to a string and prints it to log and stdout'''
    now = datetime.datetime.now().strftime("%d/%m/%Y %H:%M:%S")
    if outdir is not None:
        with open(os.path.join(outdir, "log.txt"), "a") as outfile:
            print(f"{now}\t{msg}", file=outfile)
    print(f"{now}\t{msg}")
    sys.stdout.flush()


def print_stdout_stderr(stdout, stderr, outdir):
    """ Prints stdout and stderr of a subprocess to log and stdout"""
    for msg in (stdout, stderr):
        if len(msg) > 0:
            print_log(msg.decode(), outdir)
            
            
def create_dir(outdir):
    if os.path.isdir(outdir) is False:
        os.mkdir(outdir)