from igakit.io import PetIGA,VTK
from numpy import linspace
import glob
from multiprocessing import Pool
import time
	
	# read in discretization info and potentially geometry
nrb = PetIGA().read("igaphase.dat")
	
	# enter the refinement factor
refinement = 1
	
	# write a function to sample the nrbs object
uniform = lambda U: linspace(U[0], U[-1], int(len(U)*refinement))
	
	# function to print the fields into VTK files
def print_file(infile):
	sol = PetIGA().read_vec(infile,nrb)
	outfile = infile.split(".")[0] + ".vtk"
	VTK().write(outfile,       # output filename
		nrb,                    # igakit NURBS object
		fields=sol,             # sol is the numpy array to plot
		scalars={'phi':0})
	
if __name__ == '__main__':
	list_of_files = glob.glob("ch2d*.dat")
	t0 = time.time()
	p = Pool(24)
	p.map(print_file, list_of_files)
	t1 = time.time()
	print('Total post-processing time = %f secs'%(t1 - t0))
