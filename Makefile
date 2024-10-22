src=./src
targets=$(src)/grid.for $(src)/interpolation.for $(src)/M120_cla.for $(src)/M410_bc_new.for $(src)/M801_pvoutputs.for $(src)/main_mini.for    
flags=

main:	
	mpifort -o main.exe $(targets) $(flags)
