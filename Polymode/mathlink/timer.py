from __future__ import division
import logging, datetime, time

DefaultName = 1

class RunTime:
	def __init__(self, time=0, label=""):
		self.time = time
		self.count = 1
		self.label = label

	def average(self):
		return self.time/self.count
	
	def __iadd__(self, t):
		self.time += t
		self.count += 1
		return self

	def __cmp__(self, y):
		return cmp(self.time, y.time)

	def __str__(self):
		return "%s \t %.4fs \t %g \t %.6fs" % \
			(self.label, self.time, self.count, self.average())

class timer:
	def __init__(self):
		self.start_time = self.timenow()
		self.starttimes = {}
		self.laptimes = {}
		self.times = {}
		self.names = []

	def __popname__(self,name):
		if name==DefaultName:
			name = self.names.pop()
		else:
			try:
				self.names.remove(name)
			except:
				pass
		return name
		
	def __pushname__(self,name):
		self.names.append(name)

	def __lastname__(self,name):
		if not name:
			name = self.names[-1]
		return name
	
	def __str__(self):
	 	sorted_keys = self.times.keys(); sorted_keys.sort()
		info_str = " ".join( [ "%s: %.4fs," % (name,self.get_seconds(name)) 
							for name in sorted_keys ] )
		return info_str
	
	def report(self):
		"Print information on the timed events"
	 	if len(self.times)<1: return

	 	ts = sorted(self.times.values(), reverse=1)
		maxlen = max([len(str(t.label)) for t in ts])
		
		print "Name%s Time \t\tCalls \tTime per call" % (" "*(maxlen-4))
		
		for t in ts:
			gap = " "*(maxlen-len(str(t.label)))
			print "%s%s %.4fs \t%d \t%.6fs" % \
				(t.label, gap, t.time, t.count, t.average())
	
	def timenow(self):
		now = time.time()
		return now

	def reset(self, name=None):	
		if name is None:		#Clear all
			self.times.clear()
		else:
			name = self.__popname__(name)
			if self.times.has_key(name):
				self.times.pop(name)

	def start(self, name=DefaultName):
		self.__pushname__(name)
		self.starttimes[name] = self.timenow()

	def stop(self, name=DefaultName):
		name = self.__popname__(name)
		now = self.timenow()
		end_time = now-self.starttimes.pop(name)
		self.starttimes[name] = now
		
		try:
			self.times[name] += end_time
		except KeyError:
			self.times[name] = RunTime(end_time, name)
		
		return end_time

	def lap(self,name=DefaultName):
		name = self.__lastname__(name)
		now = self.timenow()
		self.laptimes[name] = now
		return now - self.starttimes[name]

	def lapdelta(self,name=DefaultName):
		name = self.__lastname__(name)
		now = self.timenow()
		try:
			laptime = now - self.laptimes[name]
		except KeyError:
			laptime = now - self.starttimes[name]
		self.laptimes[name] = now
		return laptime

	def get_seconds(self,name=DefaultName):
		name = self.__lastname__(name)
		return self.times[name].time

	def get_time(self,name=DefaultName,show_micro=False):
		name = self.__lastname__(name)
		time = self.times[name].time
		
		second = int(round(time)%60)
		minute = int((time//60)%60)
		hour = int(time//3600)

		if show_micro:
			microsecond = int(round( 1e6*(time%1) ))
			second = int( (time//1)%60 )
			t = datetime.time(hour,minute,second,microsecond)
		else:
			t = datetime.time(hour,minute,second)
		return t


class Marker:
	def __init__(self, verbose=False, filename=None):
		if filename:
			self.file=open(filename, 'w')
			self.file.write('set x2tics axis ( ')
		self.start_time = self.timenow()
		self.times = []
		self.verbose=verbose

	def __del__(self):
		if hasattr(self,'file'):
			self.file.seek(-1,1)  #Eliminate last commma
			self.file.write(' )')
			self.file.close()     #And close file
			
	def __str__(self):
		info_str = "\n".join( [ "%s: %s %s" % (name, self.pretty_time(t),self.pretty_elapsed(t,micro=True)) 
							for t,name in self.times] )
		return info_str

	def __call__(self, name):
		self.mark(name)
	
	def timenow(self):
		now = time.time()
		return now

	def mark(self, name=None):
		t = [self.timenow(),name]
		self.times += [t]   				#Save the current time and name
		self.save(t)							#Save in file, if it exists
		if self.verbose:
			print "%s: %s" % (t[1], self.pretty_elapsed(t[0]))

	def save(self,t):
		if hasattr(self,'file'):
			self.file.write('" %f",' % (t[0]-self.start_time))

	def pretty_time(self,t):
		return time.strftime("%D %H:%M:%S", time.localtime(t))
		
	def pretty_elapsed(self,endtime,micro=False):
		time = endtime-self.start_time
		
		second = int(round(time)%60)
		minute = int((time//60)%60)
		hour = int(time//3600)

		if micro:
			microsecond = int(round( 1e6*(time%1) ))
			second = int( (time//1)%60 )
			t = "%02d:%02d:%02d:%d" % (hour,minute,second,microsecond)
		else:
			t = "%02d:%02d:%02d" % (hour,minute,second)
		return t

#Optional timers
try:
	import resource
	class ResourceTimer(timer):
		def timenow(self):
			return resource.getrusage(resource.RUSAGE_SELF)[0]
except ImportError:
	pass



#Decorators:
def time_function(tick=None, prefix=""):
	if tick is None: #Use the root timer
		tick = root_timer
	def decorator(func):
		def timer_wrapper(*_args, **_kwargs):
			tick.start(prefix+func.__name__)
			retargs = func(*_args, **_kwargs)
			tick.stop(prefix+func.__name__)
			return retargs
		return timer_wrapper
	return decorator
	
def mark_function(tick=None):
	if tick is None: #Use the root timer
		tick = root_timer
 
	def decorator(func):
		def timer_wrapper(*_args, **_kwargs):
			tick.mark(func.__name__+" In")
			retargs = func(*_args, **_kwargs)
			tick.mark(func.__name__+" Out")
			return retargs
		return timer_wrapper
	return decorator

def report(tick=None):
	if tick is None: #Use the root timer
		tick = root_timer
	tick.report()

def reset(tick=None):
	if tick is None: #Use the root timer
		tick = root_timer
	tick.reset()

#Create root timer
root_timer = timer()

