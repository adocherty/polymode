# -*- coding: utf-8 -*-
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
Classes for solving on MPI systems with mpi4py

"""

from __future__ import division
import os, sys, traceback, logging, pickle

import datetime as dt
import numpy as np

from numpy import *
from numpy.lib.scimath import sqrt

#Placeholder for MPI
MPI=None

class Status(object):
    pending = 0
    running = 1
    complete = 2
    senttonode = 3
    cancelled = 4
    nodepends = 5
    failed = -1

    @classmethod
    def asstring(cls, status):
        return cls.listall()[status]
    
    @classmethod
    def listall(cls):
        return { cls.failed:"Failed", cls.pending:"Pending", \
            cls.running:"Running", cls.senttonode:"Sent to node", \
            cls.complete:"Done", cls.complete:"Done", cls.cancelled:"Cancelled", \
            cls.nodepends:"AutoCancelled" }

class Queue(list):
    """
    Contains one or more queue items as a list, with extra queue information.
    """
    comment = ""
    def __init__(self, sequence=None):
        if sequence is None:
            list.__init__(self)
        else:
            list.__init__(self, sequence)
        
        #Copy comment attribute if there is one
        if hasattr(sequence, 'comment'):
            self.comment = sequence.comment
    

class QueueItem(object):
    """
    Contains all the information to run a job on one node, and to continue
    the job where it left off.
    """
    def __init__(self, solver=None, dependancies=[], parameters={}):
        self.solver = solver
        self.dependancies = dependancies
        self.id = id(self)
        self.parameters = parameters
        
        #Copy depenancy information from solver
        if hasattr(solver, 'dependancies'):
            self.dependancies = solver.dependancies
        if hasattr(solver, 'id'):
            self.id = solver.id

        #Vital task information
        self.status = Status.pending
        self.result = None
        
        #Set additional information
        self.started_at = None
        self.finished_at = None
        self.checkpoints = []

    def __eq__(self,other):
        return self.status == other

    def __repr__(self):
        return "<Job: %s>" % (Status.asstring(self.status))

    def __str__(self):
        return "Job %s: %s" % (Status.asstring(self.status), self.parameters)

    def calculate(self):
        if self.started_at is None:
            self.set_current_start_time()
            self.status = Status.running
            self.started_at = dt.datetime.now()
            self.initialize_external_data()

        #Wrap in try/except so we can continue with other jobs if one fails
        try:
            #Run calculation
            self.solver.calculate(**self.parameters)

            if self.solver.isfinished():
                #Update status if finished
                self.status = Status.complete
                self.finished_at = dt.datetime.now()
            else:
                #This is a checkpoint save otherwise
                self.checkpoints += [ dt.datetime.now() ]

            #Save item data externally, if supported
            self.save_external_data()
        except: #Catch all errors and save the information
            self.result = traceback.format_exc()
            self.status = Status.failed
            logging.error("Calculation of item failed")
            logging.error(self.result)
        
            #Finalize job after calculation if failed
            try:
                self.solver.finalize()
            except:
                logging.error("Finalization of item failed")
                logging.error(traceback.format_exc())

    def save_external_data(self):
        """
        Append new data to individual item datafile
        Requires solver to have two methods get_data() and clear_data() for this purpose
        """
        if hasattr(self, 'datafile'):
            #Open a file for appending
            savefile = open(self.datafile, 'ab')

            #Save data in separate file
            pickle.dump( self.solver.get_data(savefile), savefile)
            savefile.close()
            
            #Clear solver data once saved
            self.solver.clear_data()

    def initialize_external_data(self):
        "Initialize the external datafile"
        if hasattr(self, 'datafile'):
            open(self.datafile, 'wb').close()
    
    def info(self):
        print "Solver status: %s" % (Status.asstring(self.status))

        #Completion and timing information
        if self.started_at is not None:
            print "Solve started at: %s" % (self.started_at)
        if self.finished_at is not None:
            print "Solve ended at: %s, Solve took %s" % (self.finished_at,
                self.finished_at-self.started_at)
        
        #Checkpointing information
        if len(self.checkpoints):
            print "Checkpointed %d times. Last time at %s" % (len(self.checkpoints), self.checkpoints[-1])

        if self.isstatus(Status.nodepends):
            print "Solve dependancies not fulfilled."
        if self.isfailed():
            print "Solve failed."
            print self.result

    def current_job_runtime(self):
        return dt.datetime.now() - self.current_run_start_time
    
    def set_current_start_time(self):
        self.current_run_start_time = dt.datetime.now()
    
    def isfinished(self):
        return self.status in [Status.complete, Status.failed,
            Status.cancelled, Status.nodepends]
            
    def iscompleted(self):
        return self.status==Status.complete
    def ispending(self):
        return self.status==Status.pending
    def isfailed(self):
        return self.status==Status.failed
    def isrunning(self):
        return self.status==Status.running
    def isstatus(self, status):
        return self.status==status

class Worker:
    def __init__(self, comm, master):
        self.comm = comm
        self.master = master
    
    def run(self):
        """
        Main worker loop
        Wait for master to send a job, calculate it and return it to the master node
        """
        #Receive first job from master
        sts = MPI.Status()
        
        item = self.comm.recv(source=self.master, tag=MPI.ANY_TAG, status=sts)

        #Until we're sent a termination message
        tag_terminate = MPI.COMM_WORLD.Get_attr(MPI.TAG_UB)
        logging.info("%d: %s, %s, %s" % (self.comm.Get_rank(), item, sts.tag, tag_terminate))
        while sts.tag <> tag_terminate:
            #Calculate item and send back to master with the same tag
            item.calculate()
            logging.info("Finished calculation .. sending to master")
            self.comm.send(item, dest=self.master, tag=sts.tag)

            #If the item is finished then flag to wait for new job
            #Otherwise skip receiving a new job and continue with this one              
            if item.isfinished():
                item.finalize()
                item = self.comm.recv(source=self.master, tag=MPI.ANY_TAG, status=sts)
        
class Master:
    def __init__(self, comm, queue, savename):
        self.comm = comm
        self.queue = queue
        self.save_filename = savename
        self.last_saved = None
        self.min_save_interval = 120        #Only save once every 2 mins

        self.nworkers = comm.Get_size()-1
        if self.nworkers == 0:
            raise RuntimeError, "MPI solver should be run on more than one node."

        #Statistics
        self.mpi_calls = []
        
    def save_state(self, force_save=False):
        logging.info("Saving state")

        #Don't save to often
        if (not force_save) and (self.last_saved is not None):
            if (dt.datetime.now()-self.last_saved).seconds<self.min_save_interval:
                logging.info("Not saving, time since last save: %d" % (dt.datetime.now()-self.last_saved).seconds)
                return

        #Save the queue as a whole to a new file
        save_filename_temp = self.save_filename + '.tmp'
        pickle.dump(self.queue, open(save_filename_temp,'wb'), protocol=-1)
        
        logging.debug("Saved to %s" % save_filename_temp)
        
        #Remove the old file, if there is one
        try:
            os.remove(self.save_filename)
        except OSError:
            pass

        #Copy the new queue to the true name
        os.rename(save_filename_temp, self.save_filename)
        logging.debug("Now moved %s->%s" % (save_filename_temp, self.save_filename))
        
        self.last_saved = dt.datetime.now()

        #Save profiling info
        #pickle.dump(self.mpi_calls, open('mpi_call_info.prof','w'))

    def load_state(self):
        self.queue = pickle.load( open(self.save_filename,'rb') )
    
    def receive(self):
        sts = MPI.Status()
        #Blocking probe
        self.comm.Probe(MPI.ANY_SOURCE, MPI.ANY_TAG, sts)
        
        #Receive object data - note this uses 'comm.recv' to get a pickled
        #object, not the fast 'comm.Recv' which works on nump arrays
        data = self.comm.recv(source=sts.source, tag=sts.tag, status=sts)
        return sts.source, sts.tag, data

    def send_next(self):
        #Find next job that needs to be calculated & has all dependencies
        jobs_pending = False
        sent_job = False
        follow_neff=[]
        for tag, item in enumerate(self.queue):
            #Find next pending item that has dependencies fulfilled
            if item.ispending() and self.item_dependancies_fulfilled(item):
                jobs_pending = True
                break
            item = None

        #Send the item to the worker
        if item and len(self.available_workers)>0:
            #Choose a destination node from list of available
            dest = self.available_workers.pop(0)
            sent_job = True

            #Send job to node
            self.comm.send(item, dest=dest, tag=tag)
            #self.mpi_calls += mpi_log_info("send", dest, item)
            sent_job = True
            
            #Tag the master copy as running
            #Note:  the worker will seperately tag it as running so it knows
            #       that the job is new
            item.status = Status.running
            logging.info( "Sent job %d to node %d" % (tag, dest) )
        return jobs_pending, sent_job

    def item_dependancies_fulfilled(self, item):
        depends = []
        for ident in item.dependancies:
            #Find matching dependancy by id
            matching_item = filter(lambda i: i.id==ident, self.queue)
            
            if len(matching_item)==1:
                depends.append( matching_item[0].solver.depends() )
            else:
                print matching_item
                raise RuntimeError, "Couldn't find depenancy"
        
        return all(depends)

    def update_job_queue(self, autocancel=True):
        "Cancel jobs that are still pending but have no dependencies fulfilled"
        #for item in self.queue:
        for tag, item in enumerate(self.queue):
            if item.ispending() and not self.item_dependancies_fulfilled(item):
                item.status = Status.nodepends
    
    #Send signal to all workers to terminate
    def terminate_running_workers(self):
        print "Terminate running workers:", self.running_workers
        logging.info( "Terminating workers still running: %s" % (self.running_workers) )

        sts = MPI.Status()
        tag_terminate = MPI.COMM_WORLD.Get_attr(MPI.TAG_UB)
        while len(self.running_workers)>0:
            self.comm.Probe(MPI.ANY_SOURCE, MPI.ANY_TAG, sts)
            worker = sts.source

            logging.info( "Terminating node %d" % (worker) )
            
            self.comm.recv(source=worker, tag=sts.tag, status=sts) #Do we need this?
            self.comm.send(dest=worker, tag=tag_terminate)
            
            if worker in self.running_workers:
                self.running_workers.remove(worker)
        logging.info( "Finished terminating running nodes")
        
    def terminate_waiting_workers(self):
        logging.info( "Terminating waiting workers: %s" % (self.available_workers) )

        tag_terminate = MPI.COMM_WORLD.Get_attr(MPI.TAG_UB)
        while len(self.available_workers)>0:
            worker = self.available_workers.pop()
            logging.debug( "Terminating node %d" % (worker) )
            self.comm.send(dest=worker, tag=tag_terminate)
            
            if worker in self.running_workers:
                self.running_workers.remove(worker)
        logging.info( "Finished terminating available nodes")

    def run(self):
        """
        Main master loop
        Send jobs to workers, and wait for them to return them
        """
        #Startup & enumerate items correctly
        number_of_jobs = len(self.queue)
        number_calculated = sum([ item.isfinished() for item in self.queue ])

        item_update=None
        jobs_pending=True
        checkpoint=False
        logging.info( "MPISolver: Running %d simultaneous jobs" % (self.nworkers) )

        #Assume all workers are available and waiting
        self.available_workers = range(1,self.nworkers+1)
        self.running_workers = range(1,self.nworkers+1)

        #While we still have jobs to complete and something to run them on
        while number_calculated<number_of_jobs:
            #Send jobs to waiting workers
            sent_job = True
            while sent_job and not checkpoint:
                jobs_pending, sent_job = self.send_next()

            #if item_update is None or (not checkpoint and item_update.save_incremental):
            self.save_state()               #Checkpoint the queue

            #If there are no more jobs to be run then terminate waiting workers
            if not jobs_pending:
                self.terminate_waiting_workers()

            #If there are no jobs to wait for then quit
            if len(self.running_workers)==0:
                logging.error( "No jobs are running .. bailing master process" )
                break
            
            logging.info( "Waiting for response from %d workers" % len(self.running_workers) )
            worker_id, tag, item_update = self.receive()              #wait for any node to send us data
            self.queue[tag] = item_update                                   #Update the queue
            logging.info( "Received job %d from node %d" % (tag, worker_id) )
            
            if item_update.isfinished():        #Is the job complete? .. used for checkpointing
                number_calculated+=1
                logging.info( "Completed %d of %d jobs" % (number_calculated, number_of_jobs) )
                self.available_workers.append(worker_id)
                checkpoint = False
            else:                                                       #Just checkpointing
                #Check how long the job has been running
                runtime = item_update.current_job_runtime()
                if runtime>dt.timedelta(5,0,0):
                    logging.warning("Job has been running for a long time .. is there a problem?")
                logging.info( "Checkpointing job %d at X%%" % (tag) )
                checkpoint = True
        
        if jobs_pending and number_calculated>=number_of_jobs:
            logging.error("All jobs calculated but jobs still pending?")

        if not jobs_pending and number_calculated<number_of_jobs:
            self.update_job_queue(autocancel=True)
            logging.info("Jobs must have unfulfilled dependencies")

        #Final save of the queue, ensure we save it all
        self.save_state(True)
        
        #Terminate all waiting workers
        self.terminate_waiting_workers()

        #Wait for any active workers and kill them
        self.terminate_running_workers()    



def load_mpi():
    global MPI
    try:
        #If MPI not available this will fail
        #Save and restore working directory as mpi4py 0.5 stuffs it up
        oldcwd = os.getcwd()
        from mpi4py import MPI
        os.chdir(oldcwd)
        
    except ImportError:
        MPI = None

    return MPI

def is_mpi_run():
    load_mpi()

    if MPI is None or (MPI.COMM_WORLD.Get_size()<2):
        return False
    return True

def mpi_batch_continue(filename, **args):
    '''
    Continue solver with pre-existing queue file
    '''
    mpi_batch_solve(None, filename, restart=True, **args)

def mpi_direct_run_if_worker():
    """
    Run worker node directly, do not pass go if MPI is detected and not the master node
    """
    load_mpi()
    if MPI is None:
        logging.error("Error loading MPI: Cannot run MPI batch solve")
        #Print import error?
        return False

    master_id = 0
    comm = MPI.COMM_WORLD
    if comm.Get_rank() <> master_id:
        logging.debug( "Running worker on node %d." % comm.Get_rank() )
        Worker(comm, master_id).run()
        logging.debug( "Worker %d Finished." % comm.Get_rank() )
        return True
    else:
        return False

def mpi_is_master():
    """
    Simple test to check if we are running on the master node
    """
    load_mpi()
    if MPI is None:
        return True

    master_id = 0
    comm = MPI.COMM_WORLD
    if comm.Get_rank() == master_id:
        return True
    return False


def mpi_batch_solve(solvers, filename=None, savenumber=100, restart=False, try_again=False, comment=None):
    '''
    Top level solver for the MPI enabled solver
    '''
    load_mpi()
    if MPI is None:
        logging.error("Error loading MPI: Cannot run MPI batch solve")
        #Print import error?
        return False

    #Change jobs that are marked running when loaded to this status
    #usually pending if the job can be restarted, currently this is false
    #so change to done
    mpi_resume_running_status = Status.complete

    master_id = 0
    comm = MPI.COMM_WORLD

    #Note .. although solvers are supplied to all processes
    #the worker node doesn't and shouldn't use it
    if comm.Get_rank() <> master_id:
        logging.info( "Running worker on node %d." % comm.Get_rank() )
        Worker(comm, master_id).run()
        logging.info( "Worker %d Finished *" % comm.Get_rank() )
        return None
    
    #The master plan
    else:
        if filename is None:
            filename = "mpi_solve_%s.queue" % dt.date.today()
        logging.warning("Master node will save queue to '%s'" % filename)

        #if restart, check if an old queue exists
        queue = None
        if restart:
            try:
                queue = pickle.load(open(filename,'rb'))
            except IOError:
                logging.warning("Couldn't restart queue '%s' creating new queue" % filename)

        #Otherwise create a new queue
        if queue is None and solvers is None:
            logging.error("Cannout load queue and not given any solvers. Job will terminate.")
            queue = []
        
        elif queue is None:
            if isinstance(solvers,Queue):
                queue = solvers
            else:
                params = {'number':savenumber}
                queue = Queue([ QueueItem(s, parameters=params) for s in solvers ])
                queue.comment = comment
                
                logging.info("Created MPI job queue: %d items save every %d modes" \
                    % (len(queue),savenumber))
        else:
            #Change running jobs to pending for the calculation to restart
            for i in queue:
                i.set_current_start_time()          #Set new start time
                if i.isrunning() or i.isstatus(Status.senttonode):
                    i.status = mpi_resume_running_status
                if try_again and i.isfailed():
                    i.status = Status.pending

            nfinished = sum([i.isfinished() for i in queue])
            npending = sum([i.ispending() for i in queue])
            logging.info("Restarting with previous queue: " +
                "%d completed jobs, %d pending, %d total." % (nfinished, npending, len(queue)))

        if len(queue) < comm.Get_size()-1:
            logging.warning( "Queue has fewer items than MPI nodes, some workers will be idle" )

        amaster = Master(comm, queue, filename)
        try:
            amaster.run()
            print queue
            logging.info( "Master process finished." )
        except:
            logging.error( "Error in master process:\n%s" % traceback.format_exc() )
            modes = None

        return None

def queue_compress_modes(filename, new_size=[200,20], save_filename=None):
    "Load an existing queue and compress the modes to new_size"
    if save_filename is None: save_filename = filename

    try:
        queue = pickle.load(open(filename,'rb'))
    except:
        logging.error("Couldn't load queue named %s" % filename)
        return

    #For every solver in the queue compress the modes
    for item in queue:
        for m in item.solver.modes:
            if new_size is None:
                m.left = m.right = None
            else:
                newcoord = m.coord.new(Nr=new_size[0], Naz=new_size[1])
                m.compress(newcoord, astype=complex64)
    
    filename_temp = filename+".tmp"
    try:
        pickle.dump(queue, open(filename_temp,'wb'), protocol=-1)
    except:
        logging.error("Couldn't save compressed queue named %s" % filename)
        return

    #Remove the old file, if there is one
    try:
        os.remove(save_filename)
    except OSError:
        pass

    #Copy the new queue to the true name
    try:
        os.rename(filename_temp, save_filename)
    except OSError:
        logging.error("Filesystem error renaming temporary queue to %s" % save_filename)
    
def queue_load(filename):
    "Load preexisting queue and return modes"
    try:
        queue = pickle.load(open(filename,'rb'))
    except:
        logging.error("Couldn't open queue file")
        queue = None
    return queue

def queue_load_solvers(filename):
    "Load preexisting queue and return solvers"
    try:
        queue = pickle.load(open(filename,'rb'))
    except:
        logging.error("Couldn't open queue file")
        queue = None

    solvers = []
    for item in queue:
        solvers += [item.solver]
    return solvers

def queue_load_modes(filename, groups=False):
    """
    Load saved queue file and return modes
    optionally wgs are returned if groups is true
    """
    try:
        queue = pickle.load(open(filename,'rb'))
    except:
        logging.error("Couldn't open queue file")
        queue = None

    #Extract modes from queue items
    modes = []
    wgs = []
    for item in queue:
        if groups:
            modes.append([item.solver.get_data()])
            wgs.append(item.solver.wg)
        else:
            modes.append(item.solver.get_data())

    #Return waveguides if requested
    if groups:
        return modes, wgs
    else:
        return modes


