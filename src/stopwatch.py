import time
import tools




class tStopWatch:
  def __init__(self):
    self.Elapsed = 0.
    self.LastStart = None

  def start(self):
    assert self.LastStart is None
    self.LastStart = time.time()
    return self

  def stop(self):
    assert self.LastStart is not None
    self.Elapsed += time.time() - self.LastStart
    self.LastStart = None
    return self

  def elapsed(self):
    if self.LastStart:
      return time.time() - self.LastStart + self.Elapsed
    else:
      return self.Elapsed


class tJob:
  def __init__(self, name):
    self.Name = name
    self.StopWatch = tStopWatch().start()
    if self.isVisible():
      print "%s..." % name

  def done(self):
    elapsed = self.StopWatch.elapsed()
    JOB_TIMES[self.Name] += elapsed
    if self.isVisible():
      print " " * (len(self.Name) + 2), elapsed, "seconds"
  
  def isVisible(self):
    if PRINT_JOBS.get():
      return not self.Name in HIDDEN_JOBS
    else:
      return self.Name in VISIBLE_JOBS






HIDDEN_JOBS = []
VISIBLE_JOBS = []
JOB_TIMES = tools.tDictionaryWithDefault(lambda x: 0)
PRINT_JOBS = tools.tReference(True)
