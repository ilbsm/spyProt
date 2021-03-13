# source   jcollado http://stackoverflow.com/questions/1191374/subprocess-with-timeout?rq=1
from os import SEEK_END
import io
import subprocess
from subprocess import TimeoutExpired
import threading
# source https://stackoverflow.com/questions/4789837/how-to-terminate-a-python-subprocess-launched-with-shell-true
import psutil


def kill(proc_pid):
    process = psutil.Process(proc_pid)
    for proc in process.children(recursive=True):
        proc.kill()
    process.kill()


class Command(object):
    '''Run Command in the Operating System

       cmd: string - command to run
       '''

    def __init__(self, cmd):
        self.cmd = cmd
        self.process = None
        self.o = None
        self.e = None
        self.rcode = -100
        self.logfile = None
        self.logerrfile = None

    def run(self, timeout, logfile=None):
        '''Run the command provided in constructor

           Parameters
           ==========
           timeout: int - (seconds) - how long to wait before killing the command
           logfile: string - path to file to which to write the std out from running command

           '''
        self.logfile = logfile

        def target():
            try:
                if logfile:
                    logfile_handle = open(logfile, 'w+')
                    self.logerrfile = logfile + ".err"
                    logerrfile_handle = open(self.logerrfile, 'w+')
                    self.process = subprocess.Popen(self.cmd, shell=True, stdout=logfile_handle,
                                                    stderr=logerrfile_handle)
                else:
                    self.process = subprocess.Popen(self.cmd, shell=True, stdout=subprocess.PIPE,
                                                    stderr=subprocess.PIPE)
                self.o, self.e = self.process.communicate(timeout=timeout)
                if logfile:
                    logfile_handle.flush()
                    logfile_handle.close()
                    logerrfile_handle.flush()
                    logerrfile_handle.close()
            except TimeoutExpired:
                print("timeout occured for process pid: " + str(self.process.pid))
                self.process.terminate()
                self.rcode = -400
            except OSError as e:
                print("OS Error: " + str(e))
                self.rcode = -100
                return

        thread = threading.Thread(
            target=target)  # thread, not process - to easily comuniacate with PIPE and return rcode.
        thread.start()

        thread.join(timeout)
        if thread.is_alive():
            # self.process.terminate() # is stoping main process not those run from it - then thread.join() is waiting for shell (and other children) to end!
            # os.killpg(os.getpgid(self.process.pid), signal.SIGKILL) # Send the signal to all the process groups - including this one which is wating for next thread.join()...
            kill(
                self.process.pid)  # using psutil, but could be written without it - searching in processes for children.
            thread.join()
            raise Exception('Processing binary file has been terminated by timeout. Error? Loop?')
        if self.process:
            self.rcode = self.process.returncode

    def getReturnCode(self):
        return self.rcode

    def getOut(self):
        if self.o:
            return self.o.decode('utf8')
        elif self.logfile:
            with io.open(self.logfile, 'r', encoding="utf8", errors='ignore') as f:
                return "\n".join(line.strip() for line in f.read().splitlines())
        else:
            return None

    def getError(self):
        if self.e:
            return self.e.decode('utf8')
        elif self.logerrfile:
            with io.open(self.logerrfile, 'r', encoding="utf8", errors='ignore') as f:
                print("\n".join(line.strip() for line in f.read().splitlines()))
        else:
            return None

    def getTrimmedOut(self):
        if self.logfile:
            with io.open(self.logfile, 'r', encoding="utf8", errors='ignore') as f:
                try:
                    f.seek(-2000, SEEK_END)
                except:
                    pass
                return "..." + "\n".join(line.strip() for line in f.read().splitlines())
        if self.o is None:
            return None
        out = self.getOut()
        if out and len(out) > 500:
            return "..." + out[-500:]
        else:
            return out

    def getCmd(self):
        return self.cmd


if __name__ == "__main__":
    command = Command("sleep 2; echo 'Process finished'; sleep  2; echo 'gsgsg'; echo 'sdgsg'; echo 'sagwrgsg'")
    command.run(timeout=15, logfile="/tmp/log")
    print("Got out: " + command.getTrimmedOut())
    command.run(timeout=3, logfile="/tmp/log.2")
    print(command.getOut())
    print('BYE')
