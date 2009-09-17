# _*_ coding=utf-8 _*_
#
#---------------------------------------------------------------------------------
#Copyright © 2009 Andrew Docherty
#
#This program is part of Polymode.
#Polymode is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
#---------------------------------------------------------------------------------

"""
PMLSolver.py
===========

Solver class for PML implementation of Maxwell's equations
"""

from __future__ import division
import datetime

from .Solver import *
from .mathlink import blockarray,eigensolver,timer
from .difflounge import finitedifference, boundary

