from numpy import *
from .Solver import *


class Optimize(object):
	'''
	Optimize a given structure for a property that is evalutated by the evaluate function
	
	f = evaluate(solver)
	df/dp = evaluate.gradient(solver, modes, mode_gradients)

	returns a measure of the modes
	
	adjust(solver, p)
	object given a vector of paramters adjusts the structure to be evaluated
	
	solver
	solver to use
	'''
	def __init__(self, solver, adjust, evaluate):
		self.solver = solver
		self.evaluate = evaluate
		self.adjust = adjust

	def iteration(self, p):
		print "P:",p
		#Update structure
		self.adjust(self.solver, p)
		
		#Must setup again when wg changed - FIX!
		self.solver.setup()
		self.solver.clean_up()
		
		#Solve structure
		self.solver.calculate()
		
		#Calculate gradient
		
		#Evaluate results
		f = self.evaluate(self.solver)

		return f

	def __call__(self, p_start):
		from scipy import optimize
		
		output = optimize.fmin_cg(self.iteration, p_start, fprime=None, maxiter=None,\
			full_output=1, disp=1, callback=None)
		popt, fopt, func_calls, grad_calls, warnflag = output
		
		print "p:", popt, "f:", fopt
		
		return popt
		
