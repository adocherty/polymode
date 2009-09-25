#test block lu
from numpy import *

from Polymode.mathlink import blockarray, timer, ublocklu

t=timer.timer()
set_printoptions(precision=4, suppress=1, linewidth=100)

Nr=150; Naz=35; Nblock=2*Naz

data_type = complex_

print "Constructing Matricies"
A= blockarray.BlockArray( (Nr,3), (2*Naz,2*Naz), dtype=data_type )
Aublas= blockarray.BlockArray( (Nr,3), (2*Naz,2*Naz), dtype=data_type )
Ablock= blockarray.BlockArray( (Nr,3), (2*Naz,2*Naz), dtype=data_type )

A[:] = data_type( floor(100*random.random(A.shape)) )
Aublas[:]=A
Ablock[:]=A

x = data_type( 10*random.random((Nr,2*Naz)) )
shift = 0

print "Running Tests"

tlu = blockarray.TriBlockLU(overwrite=True)
t.start("Python TriBlock LU")
Ablock[:]=A
tlu(Ablock, shift=shift)
t.stop()

t.start("Python TriBlock Solve")
yblock = tlu.solve(x)
yblockt = tlu.solve_transpose(x)
t.stop()

t.start("Python TriBlock Update")
Ablock.blockview()[-1] = A.blockview()[-1]
tlu.update(Ablock, 1)
yblock = tlu.solve(x)
t.stop()

#C++ block LU
if issubclass(data_type, complexfloating):
	clu = ublocklu.cblocklu()
else:
	clu = ublocklu.dblocklu()

t.start("C++ TriBlock LU")
Aublas[:]=A
print A.shape
clu(Aublas, Nblock, shift)
t.stop()

t.start("C++ TriBlock Solve")
yublast = x.copy().reshape(Nr,Nblock);
yublas = x.copy().reshape(Nr,Nblock);
clu.solve(yublas)
clu.solve_transpose(yublast)
t.stop()

t.start("C++ TriBlock Update")
Aublas.blockview()[-2:] = 0
clu.update(A, 2)
yublas = x.copy().reshape(Nr,Nblock);
clu.solve(yublas)
t.stop()

print t.report()

#Check:
#print Aublas[:,:4]
#print Ablock[:,:4]

error_ublas = absolute(A.matvec(yublas.ravel())-x.ravel()).max()
error_ublas = max( error_ublas, absolute(A.rmatvec(yublast.ravel())-x.ravel()).max() )
error_block = absolute(A.matvec(yblock.ravel())-x.ravel()).max()
error_block = max( error_block, absolute(A.rmatvec(yblockt.ravel())-x.ravel()).max() )

print "Error in block:", error_block
print "Error in ublas:", error_ublas


