# _*_ coding=utf-8 _*_
"""
Block matrix class & routines v2

---------------------------------------------------------------------------------
Copyright © 2007 Andrew Docherty

This program is part of the ABCSolver suite.
It is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

__all__ = [ "BlockArray", "TriBlockLU" ]

from numpy import *
import logging

#Type coercion from numpy/linalg/iterative.py
_coerce_rules = {('l','l'):'l', ('l','d'):'d', ('l','F'):'F',
                 ('l','D'):'D', ('l','f'):'f', ('d','l'):'d',
                 ('l','f'):'f', ('f','l'):'f', ('l','d'):'d',
                 ('F','l'):'F', ('D','l'):'D', ('d','l'):'d',
                ('f','f'):'f', ('f','d'):'d', ('f','F'):'F',
                 ('f','D'):'D', ('d','f'):'d', ('d','d'):'d',
                 ('d','F'):'D', ('d','D'):'D', ('F','f'):'F',
                 ('F','d'):'D', ('F','F'):'F', ('F','D'):'D',
                 ('D','f'):'D', ('D','d'):'D', ('D','F'):'D',
                 ('D','D'):'D'}

class TriBlockLU:
    def __init__(self, overwrite=True):
        self.overwrite=overwrite
        
    #Tridiagonal block LU decomposition - overwrites the entries of A
    def __call__(self, A, shift=0):
        from scipy.linalg import lu_factor, lu_solve
        nrows = A.mshape[0]
        nblock = A.blockshape[0] #assume a square block!

        assert A.bandwidth==3, "Matrix bust be tridiagonal block matrix"

        #Overwrite current matrix?
        if self.overwrite:
            self.Alu = A
            
        #Create new internal matrix if none
        elif not hasattr(self, 'Alu'):
            logging.debug( "Creating new internal LU matrix" )
            self.Alu=BlockArray(A.mshape, A.blockshape, dtype=A.dtype)

        #Create new internal matrix if A is different shape
        elif (self.Alu.mshape!=A.mshape) or (self.Alu.blockshape!=A.blockshape):
            logging.debug( "Internal LU matrix incorrect shape; creating new one" )
            del self.Alu
            self.Alu = BlockArray(A.mshape, A.blockshape, dtype=A.dtype)

        #Vector for pivots
        self.shift = shift
        self.pivots = zeros((nrows, nblock), dtype=A.dtype)
        piv=0
                
        Ishift = shift*eye(nblock)
        self.Alu.set_raw_block(0,1,A.get_raw_block(0, 1) - Ishift)
        bnew = self.Alu.get_raw_block(0, 1)
        for row in range(0,nrows-1):
            a = A.get_raw_block(row+1, 0)
            b = bnew
            c = A.get_raw_block(row, 2)
            
            b, self.pivots[row] = lu_factor(b)
            # anew = a inv(b)
            anew = lu_solve((b,self.pivots[row]), a.T, trans=1).T
            bnew = A.get_raw_block(row+1, 1) - Ishift - dot(anew, c) 

            self.Alu.set_raw_block(row,1,b)                 #blu
            self.Alu.set_raw_block(row+1,0,anew)    #b anew = a

            if not self.overwrite:                      #Copy over block c, if not overwriting
                                        self.Alu.set_raw_block(row,2,c)

        #Final LU decomp of last block
        b, self.pivots[nrows-1] = lu_factor(bnew)
        self.Alu.set_raw_block(nrows-1,1,b)


    #Update last n rows of LU decomposition stored in A from Aupdate
    def update(self, Aupdate, uprows=None):
        from scipy.linalg import lu_factor, lu_solve
        nrowsup = Aupdate.mshape[0]
        nrows = self.Alu.mshape[0]
        nblock = self.Alu.blockshape[0]
        Ishift = self.shift*eye(nblock)

        #Default to update all rows
        if uprows is None: uprows=Aupdate.mshape[0]

        #Just check the block shapes
        assert Aupdate.blockshape==self.Alu.blockshape, "Block shapes not the same!"

        #b[i-1] is already LU decomposed from original decomposition
        bnew = self.Alu.get_raw_block(nrows-uprows-1, 1)
        pivot = self.pivots[nrows-uprows-1]
        for row in range(-uprows,0):
            a = Aupdate.get_raw_block(nrowsup+row,0)
            b = bnew
            if row==-uprows:    #Get c[i-1] from original matrix at start row
                                c = self.Alu.get_raw_block(nrows+row-1,2)
            else:                           #Then from update at sucessive rows
                                c = Aupdate.get_raw_block(nrowsup+row-1,2)
            
            #anew[i] bnew[i-1] = a[i]
            anew = lu_solve((b,pivot), a.T, trans=1).T

            #bnew[i] =  b[i] - a[i-1]c[i]
            bnew = Aupdate.get_raw_block(nrowsup+row, 1)  - Ishift - dot(anew, c) 
            bnew, pivot = lu_factor(bnew)
            
            #Store calculated blocks
            self.pivots[nrows+row] = pivot
            self.Alu.set_raw_block(nrows+row,1,bnew)
            self.Alu.set_raw_block(nrows+row,0,anew)
            self.Alu.set_raw_block(nrows+row-1,2, c)

    #Solve with block LU decomposed A
    def solve(self, din):
        from scipy.linalg import lu_solve
        nrows = self.Alu.mshape[0]
        nblock = self.Alu.blockshape[0] #assume a square block!

        d = din.reshape((nrows,nblock))
        x = zeros(d.shape, dtype=self.Alu.dtype)
        
        #Forward substitution pass
        # dnew[0] = d[0]
        # dnew[i] = d[i]-anew[i].dnew[i-1]
        x[0] = d[0]
        for row in range(0,nrows-1):
            x[row+1] = d[row+1] - dot(self.Alu.get_raw_block(row+1,0), x[row])

        #Backward substitution
        for row in range(nrows-1,-1,-1):
            if row<nrows-1:
                x[row] -= dot(self.Alu.get_raw_block(row,2), x[row+1])
            x[row] = lu_solve((self.Alu.get_raw_block(row,1),self.pivots[row]), x[row])

        return x.reshape(din.shape)

    #Solve transpose problem with block LU decomposed A
    def solve_transpose(self, din):
        from scipy.linalg import lu_solve
        nrows = self.Alu.mshape[0]
        nblock = self.Alu.blockshape[0] #assume a square block!

        d = din.reshape((nrows,nblock))
        x = zeros(d.shape, dtype=self.Alu.dtype)
        
        #Forward substitution pass
        # b[0].T dnew[0] = d[0]
        # b[i].T dnew[i] = d[i] - c[i-1].T dnew[i-1]
        x[0] = d[0]
        for row in range(0,nrows):
            if row>0:
                x[row] = d[row] - dot(self.Alu.get_raw_block(row-1,2).T, x[row-1])
            if any(isnan(x[row])) or any(isinf(x[row])):
                print(row, x[row])
            x[row] = lu_solve((self.Alu.get_raw_block(row,1),self.pivots[row]),\
                     x[row], trans=1)

        #Backward substitution
        # x[i] = d[i] - anew[i+1] x[i+1]
        for row in range(nrows-2,-1,-1):
            x[row] -= dot(self.Alu.get_raw_block(row+1,0).T, x[row+1])

        return x.reshape(din.shape)

class BlockArray(ndarray):
    '''
        Class for block matricies, a subclass of numpy.ndarray
        new_block_array  =  BlockArray(shape, blockshape = (1,1), dtype=float)
        create a block array with shape=(Nx,Ny) and blockshape is the shape of the block.
    '''
    memtotal = 0
    def __new__(subtype, shape, blockshape = (1,1), dtype=float):
        #Allow square shapes easily
        if isscalar(blockshape): blockshape = (blockshape,blockshape)
        newshape = shape[0]*blockshape[0], shape[1]*blockshape[1]
        
        self = ndarray.__new__(subtype, newshape, dtype=dtype)
        return self

    def __init__(self, shape, blockshape = (1,1), dtype=float):
        #Allow square shapes easily
        if isscalar(blockshape): blockshape = (blockshape,blockshape)
        if shape[1]%2==0:
            raise RuntimeError("Currently BlockArray only takes odd bandwidths")
                
        self.blockshape = blockshape
        self.mshape = shape
        
        #Set all elements to zero!
        self.__array__()[:]=0

#   def __array_finalize__(self, obj):
#       if isinstance(obj, BlockArray):
#           self.mshape=obj.mshape
#           self.blockshape=obj.blockshape
#           
#       elif obj is not None:
#           raise RuntimeError, "Allocating block matrix from object", obj
#   
#       #Calculate size of array, first get size of data, in bytes (hopefully)
#       self.memsize = (self.itemsize*prod(self.shape)/1024.**2)
#       BlockArray.memtotal += self.memsize
#       logging.debug( "Allocated %.2fMb for a total of %.2fMb for block matricies" \
#           % (self.memsize, BlockArray.memtotal) )

#   def __del__(self):
#       try:
#           BlockArray.memtotal-=self.memsize
#           logging.debug( "Deallocating %.2fMb; total now used %.2fMb" \
#               % (self.memsize, BlockArray.memtotal) )
#       except AttributeError:
#           pass
    
    def __repr__(self):
        return "<%dx%d block matrix with block size %dx%d of type '%s'>" % \
                (self.get_full_shape() + (self.dtype.type,))
    def __str__(self):
        return ("BlockArray{%dx%d, %dx%d}\n" % self.get_full_shape()) + str(self.__array__())

    #----------------------------------------------------------------------------
    #       Element access routines - for block, row, col, elementwise
    #----------------------------------------------------------------------------

    # Get blocks of banded matrix
    def set_raw_block(self,i,j,value):
        self.blockview().__setitem__((i,j),value)

    def get_raw_block(self,i,j):
        return self.blockview().__getitem__((i,j))

    def row_block_offsets(self,blockrow):
        '''
        Calculate the limits in block coords of the current row (also given in block coords)
        '''
        nrowblocks = self.mshape[0]
        block_loffset = max(0,blockrow - self.bandwidth//2)
        block_roffset = min(blockrow + (self.bandwidth+1)//2, nrowblocks)
        return block_loffset, block_roffset

    def row_offsets(self,row):
        '''
        Calculate the limits of the current row, namely the bandwidth range
        '''
        block_loffset, block_roffset = self.row_block_offsets(row//self.blockshape[0])
        return block_loffset*self.blockshape[1], block_roffset*self.blockshape[1]

    def map_raw_block(self,blockrow):
        nrowblocks = self.mshape[0]
        loffset = max(0,self.bandwidth//2 - blockrow)
        roffset = min(0,(nrowblocks-blockrow-1) - self.bandwidth//2)
        return loffset,roffset

    #Routines to help access rows & columns
    def __row_index__(self,row):
        #Kock off the first & last blocks which aren't in the true matrix
        nrowblocks = self.mshape[0]
        loffset = max(0,self.bandwidth//2 - row//self.blockshape[0])
        roffset = min(0,(nrowblocks-row//self.blockshape[0]-1) - self.bandwidth//2)
        
        if loffset==0 and roffset==0:
            return row
        else:
            return (row,slice(loffset*self.blockshape[1],\
                    (self.bandwidth+roffset)*self.blockshape[1]))

    def __column_index__(self,col):
        nrowblocks = self.mshape[0]
        loffset = max(0,self.bandwidth//2 - col//self.blockshape[0])
        roffset = min(0,(nrowblocks-col//self.blockshape[0]-1) - self.bandwidth//2)
        
        indicies = []
        bw = self.bandwidth; bs0,bs1 = self.blockshape
        startrow = col//bs1 - bw//2         #Starting blockrow
        minorcol = col%bs1                          #Column within block
        for i in range(loffset,bw+roffset):
            indicies += [ (slice((startrow+i)*bs0,(startrow+i+1)*bs0), (bw-i-1)*bs1+minorcol) ]
        return indicies
    
    #Return row as an array: not a blockarray
    def get_row(self,row):
        return self.array()[self.__row_index__(row)]

    #Set row as an array, not a blockarray. Otherwise blockarray objects get instanciated
    def set_row(self,row,value):
        self.array()[self.__row_index__(row)] = value

    def get_column(self,col):
        "Return column col of represented matrix with matrix bandwidth"
        col_indicies = self.__column_index__(col)
        col = zeros(0, dtype=self.dtype)
        for ind in col_indicies:
            col = append(col, self[ind])
        return col

    def set_diag(self, value):
        #Iterate over diagonal length or value length, whichever is smallest
        if not iterable(value):
            diag_length = self.shape[0]
        else:
            diag_length = min(size(value), self.shape[0])

        #Iterate over rows
        for row in range(diag_length):
            col = row%self.blockshape[1] + self.blockshape[1]*(self.bandwidth//2)
            if iterable(value):
                self[row,col] = value.flat[row]
            else:
                self[row,col] = value

    def get_diag(self):
        #Iterate over diagonal length or value length, whichever is smallest
        diag = empty(self.shape[0], dtype=self.dtype)

        #Iterate over rows
        for row in range(self.shape[0]):
            col = row%self.blockshape[1] + self.blockshape[1]*(self.bandwidth//2)
            diag[row] = self[row,col]
        return diag


    #----------------------------------------------------------------------------
    #       Utilities to access the underlying dense matrix
    #----------------------------------------------------------------------------

    #Caclulate banded block coord from dense block coord
    def     translate_block_coords(self,i,j):
            mapi=i
            mapj=j-mapi+self.mshape[1]//2
            return mapi, mapj

    # Get blocks of matrix represented here
    def set_block(self,i,j,value):
        mapi,mapj = self.translate_block_coords(i,j)
        self.blockview().__setitem__((mapi,mapj),value)

    def get_block(self,i,j):
        mapi,mapj = self.translate_block_coords(i,j)
        return self.blockview().__getitem__((mapi,mapj))

    def get_dense_row(self,row):
        "Returns a row of the represented dense matrix"
        row_length = self.mshape[0]*self.blockshape[0]
        
        loff,roff = self.row_offsets(row)
        matrix_row = zeros(row_length, dtype=self.dtype)
        matrix_row[loff:roff] = self.get_row(row)
        return matrix_row

    def get_dense_column(self,row):
        "Returns a column of the represented dense matrix"
        raise NotImplementedError("No right matvec yet")

    def toarray(self):
        "Converts the block matrix to a dense matrix"
        ncols = self.mshape[0]*self.blockshape[0]
        nrows = self.mshape[0]*self.blockshape[1]

        dense = zeros((ncols,nrows), dtype=self.dtype)
        for row in arange(nrows):
            loff,roff = self.row_offsets(row)
            dense[row,loff:roff] = self.get_row(row)
        return dense

    #----------------------------------------------------------------------------
    #       Shape routines
    #----------------------------------------------------------------------------
    def blockview(self):
        #Returns a view of the matrix that allows the blocks to be accessed
        #as the first two indices
        return self.array().reshape(self.viewshape).swapaxes(1,2)

    def get_full_shape(self):
        return self.mshape + self.blockshape
    
    def get_viewshape(self):
                                            return self.mshape[0],self.blockshape[0],self.mshape[1],self.blockshape[1]
    viewshape = property(get_viewshape)     
    
    def get_bandwidth(self):
        return self.mshape[1]
    bandwidth=property(get_bandwidth)

    def array(self):
        return self.__array__()

    #----------------------------------------------------------------------------
    #       Other routines
    #----------------------------------------------------------------------------

    #Slices should be arrays, not BlockArrays
    def __getslice__(self, i, j):
        return self.__array__().__getslice__(i, j)

    def __getitem__(self, index):
        return self.__array__().__getitem__(index)

    #matrix vector dot product
    def matvec(self, vec):
        "Multiply matrix by vector, A.x"
        nrows = self.shape[0]
        
        typ = typeDict[_coerce_rules[vec.dtype.char,self.dtype.char]]
        x = zeros(nrows, dtype=typ)
        d = vec.ravel()         
        for row in range(0,nrows):
            loffset, roffset = self.row_offsets(row)
            x[row] = dot(self.get_row(row), d[loffset:roffset])
        return x.reshape(vec.shape)

    def rmatvec(self, vec):
        "Multiply transposed matrix by vector, A^T.x"
        nrows = self.shape[0]

        typ = typeDict[_coerce_rules[vec.dtype.char,self.dtype.char]]
        x = zeros(nrows, dtype=typ)
        d = conj(vec).ravel()       
        for row in range(0,nrows):
            loffset, roffset = self.row_offsets(row)
            x[loffset:roffset] += self.get_row(row)*d[row]
        return conj(x).reshape(vec.shape)





