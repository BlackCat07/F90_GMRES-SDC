program gmrestest
	use sdc_gmres
	use problem
	implicit none
	integer :: i, kiter, pic, steps
	double precision :: time, dt
	character(len=*), parameter :: mth = "gm_", fl1 = "nrg", fl2 = "ptc"
	character(len=8) :: fmt1 = '(I2.2)', fmt2 = '(I1)'
	character(len=2) :: spic, sit   
	character(len=1) :: snodes
	character(len=50):: name1, name2
	
	nfilter = 50			! defines the data writing step
	
	time = 0.001D0			! runtime in seconds
	steps = 500000			! number of steps
	dt = time/steps			! time step
	kiter = 1				! GMRES-SDC iteration
	pic = 3					! Picard iteration
	
	write (spic,    fmt1) pic
    write (sit,     fmt1) kiter
    write (snodes,  fmt2) colnodes
    	
	call InitParticle
	
    ! Full output names for energy and trajectory:
	name1 = mth//snodes//'_'//sit//spic//pname//fl1
	name2 = mth//snodes//'_'//sit//spic//pname//fl2
	
    call sdcgmresrun(steps, dt, kiter, pic, trim(name1), trim(name2))
	
end program gmrestest




	





