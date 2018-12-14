program gmrestest
	use sdc_gmres
	use problem
	implicit none
	integer :: i, kiter, pic, steps
	double precision :: time, dt
	character(len=*), parameter :: outfile1 = "energy", outfile2 = "trajectory"
	
	nfilter = 200			! defines the data writing step
	
	time = 0.01D0			! runtime in seconds
	steps = 10000000		! number of steps
	dt = time/steps		! time step
	kiter = 1				! GMRES-SDC iteration
	pic = 3				! Picard iteration
	
	call InitParticle
	call sdcgmresrun(steps, dt, kiter, pic, outfile1, outfile2)	
	
end program gmrestest




	





