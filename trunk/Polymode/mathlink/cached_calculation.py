# _*_ coding=utf-8 _*_
'''
Implementation of versitile cached calculation class, calculate once and access many times.

Todo:
* How do we know when to recalculate the cache?

---------------------------------------------------------------------------------
Copyright © 2008 Andrew Docherty

This program is part of ABCSolver.
ABCSolver is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''

class CachedCalculation(object):
	"A cached calculation using descriptors"
	def __init__(self, calc=None):
		self.calculate = calc
		self.name = calc.__name__
		self.cache = {}

	def __get__(self, obj, objtype=None):
		#Calculate data if no cache
		if obj in self.cache:
			return self.cache[obj]
		elif self.calculate is not None:
			self.cache[obj] = self.calculate(obj)
			return self.cache[obj]
		else:
			raise AttributeError, "No calculation specified"

	def __set__(self, obj, value):
		"Manually set cache"
		self.cache[obj] = value

	def __delete__(self, obj):
		if obj in self.cache:
			del self.cache[obj]

	def __getstate__(self):
		"Pickle all needed data, ignore cached data"
		state = self.__dict__.copy()
		if 'cache' in state: del state['cache']
		return state
	
	def __setstate__(self,state):
		"Restore data from pickler"
		self.__dict__.update(state)
		self.cache = {}
		

