import time




class tStopWatch:
  def __init__( self ):
    self.Elapsed = 0.
    self.LastStart = None

  def start( self ):
    assert self.LastStart is None
    self.LastStart = time.time()
    return self

  def stop( self ):
    assert self.LastStart is not None
    self.Elapsed += time.time() - self.LastStart
    self.LastStart = None
    return self

  def elapsed( self ):
    if self.LastStart:
      return time.time() - self.LastStart + self.Elapsed
    else:
      return self.Elapsed


class tJob:
  def __init__( self, name ):
    print "%s..." % name
    self.Name = name
    self.StopWatch = tStopWatch().start()

  def done( self ):
    elapsed = self.StopWatch.elapsed()
    print " " * ( len( self.Name ) + 2), elapsed, "seconds"

