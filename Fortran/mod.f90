module sdc_gmres
	implicit none

#if defined (NODES)
	integer, parameter 	:: colnodes   = NODES				
#endif

	integer 			:: nfilter
	integer, parameter 	:: dimm = colnodes*6		
	double precision	:: Ekin0
	
	double precision, dimension(colnodes)			:: delta_tau, q
	double precision, dimension(colnodes,colnodes)	:: Qd_E, Qd_E2, Qd_T, Smat
	double precision, dimension(3,colnodes)			:: I_m_mp1, IV_m_mp1
	
	
	contains
	
	
	subroutine progress_bar(k, steps)
		implicit none
		integer :: k, steps
		write(*,FMT="(A1,A,t21,F6.2,A)",ADVANCE="NO") achar(13), &
		"Percent Complete: ", (real(k)/real(steps))*100.0, "%"
	end subroutine progress_bar
	
	
	subroutine savefiles(dt, steps, energy, pos, vel, filename1, filename2)
		use problem
		implicit none
		integer :: steps, i
		character(len=*) :: filename1, filename2
		double precision :: energy(steps), pos(3,steps+1), vel(3, steps+1)
		double precision :: R, Phi, Z, dt
		
		write(*,*) "Writing files..."
		
		open(5, file=filename1)	
		open(6, file=filename2)	
		do i = 1, steps
			write(5,*) dt*i*nfilter, energy(i)
			call CartesToCylinCoord(pos(:,i+1), R, Phi, Z)
			write(6,*)  R, Phi, Z
		end do		
		close(5)
		close(6)
						
	end subroutine savefiles
	
	
	subroutine getMatrices(dt)
		implicit none
		integer			:: i, m, n
		double precision		:: dt

#if(NODES==3)
		delta_tau(1) = 0.D0		
		delta_tau(2) = 0.5D0 * dt		
		delta_tau(3) = 0.5D0 * dt
			
		q(1) = 0.1666666666666663D0 * dt	
		q(2) = 0.66666666666666741D0 * dt	
		q(3) = 0.1666666666666663D0 * dt
				
		Smat(1,1) = 0.0D0				
		Smat(1,2) = 0.0D0				
		Smat(1,3) = 0.0D0
		Smat(2,1) = dt*0.2083333333333332D0 	
		Smat(2,2) = dt*0.33333333333333348D0	
		Smat(2,3) =-dt*0.41666666666666741D-1
		Smat(3,1) =-dt*0.41666666666666907D-1	
		Smat(3,2) = dt*0.33333333333333393D0	
		Smat(3,3) = dt*0.20833333333333304D0
		
#elif(NODES==5)

		delta_tau(1) = 0.D0		
		delta_tau(2) = 1.7267316464601151D-1 * dt		
		delta_tau(3) = 3.2732683535398849D-1 * dt
		delta_tau(4) = 3.2732683535398843D-1 * dt		
		delta_tau(5) = 1.7267316464601157D-1 * dt
				
		q(1) = 5.0000000000000051D-2 * dt	
		q(2) = 2.7222222222222225D-1 * dt	
		q(3) = 3.5555555555555529D-1 * dt
		q(4) = 2.7222222222222225D-1 * dt	
		q(5) = 5.0000000000000093D-2 * dt
				
		Smat(1,1) = 0.0D0				
		Smat(1,2) = 0.0D0				
		Smat(1,3) = 0.0D0
		Smat(1,4) = 0.0D0				
		Smat(1,5) = 0.0D0
		
		Smat(2,1) = dt*6.7728432186156914D-2	
		Smat(2,2) = dt*1.1974476934341179D-1	
		Smat(2,3) =-dt*2.1735721866558148D-2
		Smat(2,4) = dt*1.0635824225415503D-2	
		Smat(2,5) =-dt*3.7001392424145306D-3
		
		Smat(3,1) =-dt*2.7103432186156858D-2	
		Smat(3,2) = dt*1.8343941397963098D-1	
		Smat(3,3) = dt*1.9951349964433582D-1
		Smat(3,4) =-dt*4.1597785326236023D-2	
		Smat(3,5) = dt*1.3075139242414505D-2
	
		Smat(4,1) = dt*1.3075139242414519D-2	
		Smat(4,2) =-dt*4.1597785326235981D-2	
		Smat(4,3) = dt*1.9951349964433587D-1
		Smat(4,4) = dt*1.8343941397963098D-1	
		Smat(4,5) =-dt*2.7103432186156830D-2
		
		Smat(5,1) =-dt*3.7001392424145241D-3	
		Smat(5,2) = dt*1.0635824225415480D-2	
		Smat(5,3) =-dt*2.1735721866558255D-2
		Smat(5,4) = dt*1.1974476934341180D-1	
		Smat(5,5) = dt*6.7728432186156956D-2
#endif
	
		do i = 1, colnodes
			Qd_E(i,1:i) = delta_tau(2:i)
		end do
	
		do m = 1, colnodes
		  	do n = 1,colnodes
		    	Qd_E2(m,n) = 0.5D0*Qd_E(m,n)*Qd_E(m,n)
		    	end do
		end do
		
		do m = 1,colnodes
		  	Qd_T(m,m)  = delta_tau(m)
		  	do n = 1, m - 1
		    	Qd_T(m,n) = Qd_T(m,n) + delta_tau(n) + delta_tau(n+1)
		    	end do
		end do
		Qd_T = 0.5D0 * Qd_T
	end subroutine getMatrices
	
	
	subroutine sdcgmresrun(steps, dt, kiter, pic, filename1, filename2)
		use problem
		implicit none
		integer :: k, i, kiter, pic, steps, outsteps, jnum
		double precision :: dt
		double precision, allocatable, dimension(:)		:: energy
		double precision, allocatable, dimension(:,:) 	:: pos, vel
		character(len=*) :: filename1, filename2
		
		double precision, dimension(dimm)		:: u0, u_new, u_old
		double precision, dimension(3, colnodes)	:: x_old, v_old
		character(len=*), parameter  			:: fmx = "(A, 3f18.12)", fme = "(A, 3e17.10)", fmv = "(A, 3e18.8)"
		double precision :: tmp_pos(3,2), tmp_vel(3,2)

		if (mod(steps, nfilter) /= 0) then
		write(*,"(A)") "steps / nfilter should be integer,"
		write(*,"(A)") "program has stopped working."
		stop
		end if
		
		write(*,"(A, e13.6)") "dt =", dt		
		Ekin0 = getEnergy(x0,v0)
		call getMatrices(dt)

		outsteps = steps/nfilter
		allocate(energy(outsteps), pos(3, outsteps + 1), vel(3, outsteps + 1))
		
		pos(:,1) = x0
		vel(:,1) = v0
		tmp_pos(:,1) = x0  
		tmp_vel(:,1) = v0
		
		! Main cycle:		
		do k = 1, steps	
		
			call progress_bar(k, steps)			
			call InitData(tmp_pos(:,1), tmp_vel(:,1), u0)	
			u_old = sweep(u0)	
			x_old = unpackX(u_old)
			v_old = unpackV(u_old)    

			! SDC iterarions:
			u_new = mygmres(u0, u_old, kiter, x_old)

			! Picard iterations:
			do i = 1, pic
				u_new = u0 + quad(u_new)
			end do

			x_old = unpackX(u_new)
			v_old = unpackV(u_new)	

			call finalUpdateStep(x_old, v_old, tmp_pos(:,1), tmp_vel(:,1), tmp_pos(:,2), tmp_vel(:,2))

			if (mod(k, nfilter) == 0) then
				jnum = k/nfilter
				pos(:,jnum + 1) = tmp_pos(:,2)
				vel(:,jnum + 1) = tmp_vel(:,2)
				energy(jnum) = abs(Ekin0 - getEnergy(pos(:,jnum + 1),vel(:,jnum + 1)))/Ekin0
			end if
			
			tmp_pos(:,1) = tmp_pos(:,2)  
			tmp_vel(:,1) = tmp_vel(:,2)
				
		end do
				
		write(*,*) ''
		write(*,*) ''
		write(*,fmx) 'Final X:', pos(:,outsteps + 1)
		write(*,fmv) 'Final V:', vel(:,outsteps + 1)
		write(*,fme) 'Energy error:',  energy(outsteps)
		write(*,*) ''
		
		call savefiles(dt, outsteps, energy, pos, vel, filename1, filename2)
		deallocate(energy, pos, vel)
		
	end subroutine sdcgmresrun
	
	
	subroutine InitData(x0,v0,u0)
		implicit none
		integer						:: i, j
		double precision, dimension(3)			:: x0, v0
		double precision, dimension(dimm)		:: u0
		double precision, dimension(3, colnodes)	:: x0_, v0_
	
		do i = 1,colnodes
			do j = 1,3
				x0_(j, i) = x0(j)
				v0_(j, i) = v0(j)
			end do
		end do
	
		do i = 1, colnodes		
			u0(1 + 3*(i - 1):3 + 3*(i - 1)) = x0_(:,i)			
			u0(dimm/2 + 1 + 3*(i - 1) : dimm/2 + 3 + 3*(i - 1)) = v0_(:,i)
		end do	
	end subroutine InitData
	
	
	subroutine updateIntegrals(x, v)
		use problem
		implicit none
		integer			:: i, j, k
		double precision, dimension(3, colnodes)	:: F
		double precision, dimension(3, colnodes)	:: x, v

	    	F = 0.D0; 	I_m_mp1 = 0.D0;	IV_m_mp1 = 0.D0
	    	
    		do i = 1,colnodes
      		F(:,i) = getF(x(:,i),v(:,i)) 
		end do
		
    		do j = 1,colnodes
      		do k = 1,colnodes
        		I_m_mp1(:,j) = I_m_mp1(:,j) + Smat(j,k)*F(:,k)
        		IV_m_mp1(:,j) = IV_m_mp1(:,j) + Smat(j,k)*v(:,k)
        		end do
        	end do        	
      end subroutine updateIntegrals
	
	
	function packU(x,v) result (u)
		implicit none
		integer	:: i
		double precision, dimension(3, colnodes), intent(in)	:: x, v
		double precision, dimension(dimm)				:: u	
			
		do i = 1, colnodes		
			u(1 + 3*(i - 1):3 + 3*(i - 1)) = x(:,i)			
			u(dimm/2 + 1 + 3*(i - 1) : dimm/2 + 3 + 3*(i - 1)) = v(:,i)
		end do	
	end function packU
	
		
	function unpackX(u) result (x)
		implicit none
		integer	:: i, j
		double precision, dimension(dimm), intent(in)	:: u
		double precision, dimension(3, colnodes)		:: x			
		do j = 1,colnodes
			x(1:3,j) = u(1 + 3*(j - 1) : 3 + 3*(j - 1))			
		end do	
	end function unpackX
	
	
	function unpackV(u) result (v)
		implicit none
		integer	:: i
		double precision, dimension(dimm), intent(in)	:: u
		double precision, dimension(3, colnodes)		:: v	
		do i = 1,colnodes
			v(1:3,i) = u(dimm/2 + 1 + 3*(i - 1) : dimm/2 + 3 + 3*(i - 1))			
		end do
	end function unpackV
	
	
	function  boris_trick(B_new, v_old, c_i, coeff) result (v_new)
		use problem
		implicit none
		integer	:: i
		double precision, intent(in)				:: coeff
		double precision, dimension(3), intent(in)	:: B_new, v_old, c_i
		double precision, dimension(3)			:: v_new		
		double precision, dimension(3) :: E, t, s, v_min, v_plu, v_star, cross		
		
	    	E      = 0.D0
	    	t      = coeff*B_new*alpha	    	
	    	s      = 2.D0*t/(1.D0 + (t(1)*t(1) + t(2)*t(2) + t(3)*t(3)))
	    	v_min  = v_old + E + 0.5D0*c_i
	    	
	    	cross  = 	[v_min(2) * t(3) - v_min(3) * t(2), &
	    			 v_min(3) * t(1) - v_min(1) * t(3), &
	    			 v_min(1) * t(2) - v_min(2) * t(1)]	    	
	    	v_star = v_min + cross	    	
	    	cross  = 	[v_star(2) * s(3) - v_star(3) * s(2), &
	    			 v_star(3) * s(1) - v_star(1) * s(3), &
	    			 v_star(1) * s(2) - v_star(2) * s(1)]	    			 
	    	v_plu  = v_min + cross
	    	v_new  = v_plu + E + 0.5D0*c_i
	end function boris_trick
	
	
	function sweep(u) result (u_new)
		use problem
		implicit none
		integer	:: i, j
		double precision, dimension(dimm), intent(in)		:: u	
		double precision, dimension(dimm)	:: u_new
		double precision, dimension(3, colnodes)	:: x, v, xnew, vnew
		double precision, dimension(3)		:: b, bprime, cross, Bf, Func
		double precision					:: coeff

    		x = unpackX(u);	v = unpackV(u)    		
    		xnew = 0.D0;	vnew = 0.0D0;
    		
   		do j = 1, colnodes
   			b = 0.D0
     			xnew(:,j) = x(:,j)
      		do i = 1, j-1
       			xnew(:,j) = xnew(:,j) + Qd_E(j,i)*vnew(:,i)
       			Func = getF(xnew(:,i),vnew(:,i))
        			xnew(:,j) = xnew(:,j) + Qd_E2(j,i) * Func
        			b         = b + Qd_T(j,i) * Func
      		end do
      		
      		coeff     = Qd_T(j,j)
      		Bf	    = getB(xnew(:,j))            						
      		cross     =	[v(2,j) * Bf(3) - v(3,j) * Bf(2), &
	    			 	 v(3,j) * Bf(1) - v(1,j) * Bf(3), &
	    			 	 v(1,j) * Bf(2) - v(2,j) * Bf(1)] 
      		bprime    = b - coeff*cross*alpha
      		vnew(:,j) = boris_trick(Bf, v(:,j), bprime, coeff)
    		end do    
    		u_new = packU(xnew, vnew)      				
	end function sweep
	
	
	function sweep_lin(u, x_0) result (u_new)
		use problem
		implicit none
		integer	:: i, j
		double precision, dimension(3, colnodes), intent(in)	:: x_0
		double precision, dimension(dimm), intent(in)		:: u	
		double precision, dimension(dimm)	:: u_new
		double precision, dimension(3, colnodes)	:: x, v, xnew, vnew
		double precision, dimension(3)		:: b, bprime, cross, Bf, Func
		double precision					:: coeff

    		x = unpackX(u);	v = unpackV(u)
    		xnew = 0.D0;	vnew = 0.0D0;

    		do j = 1, colnodes
    		   	b = 0.D0
     			xnew(:,j) = x(:,j)
      		do i = 1, j-1
       			xnew(:,j) = xnew(:,j) + delta_tau(i+1)*vnew(:,i)
       			Func = getF(x_0(:,i),vnew(:,i))
        			xnew(:,j) = xnew(:,j) + 0.5D0*delta_tau(i+1)**2 * Func
        			b         = b + 0.5D0*(delta_tau(i) + delta_tau(i+1)) * Func
        		end do
        		
        		coeff     = 0.5D0*delta_tau(j)
      		Bf	    = getB(x_0(:,j))            						
      		cross     =	[v(2,j) * Bf(3) - v(3,j) * Bf(2), &
	    			 	 v(3,j) * Bf(1) - v(1,j) * Bf(3), &
	    			 	 v(1,j) * Bf(2) - v(2,j) * Bf(1)]	    			 		    			 		
      		bprime    = b - coeff*cross*alpha
      		vnew(:,j) = boris_trick(Bf, v(:,j), bprime, coeff)
    		end do    		
    		u_new = packU(xnew, vnew)      				
	end function sweep_lin
	
	
	function quad(u, x0) result (Q)
		implicit none
		integer	:: i, j, uu, nn
		double precision, dimension(dimm), intent(in)				:: u
		double precision, dimension(3, colnodes), optional, intent(in) 	:: x0
		
		double precision, dimension(dimm)		:: Q
		double precision, dimension(3, colnodes)	:: x, v, Qx, Qv

    		x = unpackX(u);	v = unpackV(u)
    		
    		if (present(x0)) then
    			call updateIntegrals(x0, v)
    		else
    			call updateIntegrals(x, v)
    		end if

    		Qx = 0.D0    		
    		do uu = 1, colnodes
      		do nn = 1, uu
        		Qx(:,uu) = Qx(:,uu) + IV_m_mp1(:,nn)
        		end do
        	end do
    		
    		Qv = 0.D0
    		do uu = 1, colnodes
      		do nn = 1, uu
        		Qv(:,uu) = Qv(:,uu) + I_m_mp1(:,nn)
        		end do
        	end do
    
    		Q = packU(Qx, Qv)    				
	end function quad
	
	
	function mygmres(b, x, m, x_old) result (u_new)
		implicit none
		integer :: n, k
		integer, intent(in) :: m
		double precision, dimension(dimm), intent(in)		:: b, x
		double precision, dimension(3, colnodes), intent(in)	:: x_old	
		double precision, dimension(dimm)				:: u_new, r, diff, s, c, A
		
		double precision, dimension(m)		:: y	
		double precision, dimension(dimm + 1)	:: e1, g
		double precision, dimension(dimm, m + 1)	:: Q
		double precision, dimension(m + 1, m)	:: H
		double precision 					:: b_norm, r_norm, error, nu, sm
    		
    		e1 = 0.D0;	g = 0.D0;	s = 0.D0;	c = 0.D0;    Q = 0.D0;	H = 0.D0;	

		b_norm = getnorm(b,dimm)
		if (b_norm == 0) b_norm = 1.D0
		
		A 	 = x - quad(x, x_old)		
		r      = sweep_lin(b - A, x_old)
	 
		r_norm = getnorm(r,dimm)
		error  = r_norm/b_norm
		
		Q(:,1) = r*r_norm**(-1.D0)
		e1(1)  = 1.D0
		g      = r_norm*e1

    		do n = 1, m
    		    		
    			call arnoldi(Q, n, Q(:,n+1), H(1:n+1,n), x_old, m)    			
    			if (n > 1) H(1:n+1,n) = apply_qn(H(1:n+1,n), c, s, n)
    			     			    			
			nu   = sqrt(H(n,n)**2 + H(n+1,n)**2)
			c(n) = H(n,n)/nu
			s(n) = -H(n+1,n)/nu
			H(n,n)   = c(n)*H(n,n) - s(n)*H(n+1,n)
			H(n+1,n) = 0.D0
			g(n+1) = s(n)*g(n)
			g(n)   = g(n)*c(n)
			
    		end do

    		! Solve linear system (back substitution algorithm):
    		y(m) = g(m)/H(m,m)
    		do n = m-1, 1, -1
    			sm = 0.D0
    			do k = n+1, m
    				sm = sm + y(k)*H(n,k)
    			end do
    			y(n) = (g(n) - sm)/H(n,n)
    		end do

    		u_new = x + matmul(Q(:,1:m),y)
    		
	end function mygmres
	
	
	subroutine arnoldi(Q, n, v, h, x_old, m)
		implicit none
		integer :: n, j, m
		double precision, dimension(n + 1)		:: h
		double precision, dimension(dimm)		:: v, A
		double precision, dimension(dimm, m + 1)	:: Q	
		double precision, dimension(3, colnodes)	:: x_old
			
		A = Q(:,n) - quad(Q(:,n), x_old)
		v = sweep_lin(A, x_old)
		
		do j = 1,n
			h(j) = dot_product(Q(:,j), v)
			v = v - h(j)*Q(:,j)
		end do
				
		h(n+1) = getnorm(v,dimm)

		v = v * 1.D0/h(n+1)
		
    	end subroutine arnoldi
    	
    	
    	function apply_qn(h, c, s, n) result (h_new)
		implicit none
    		integer :: i
		integer, intent(in) :: n
    		double precision, dimension(dimm), intent(in)	:: c, s
    		double precision, dimension(n + 1), intent(in)	:: h
        	double precision, dimension(n + 1)			:: h_new
        	double precision :: temp
    		
        	h_new = h    		
		do i = 1, n - 1
			temp   	 = c(i)*h_new(i) - s(i)*h_new(i + 1)
			h_new(i + 1) = s(i)*h_new(i) + c(i)*h_new(i + 1)
			h_new(i)   	 = temp
		end do

  	end function apply_qn
  	
  	
  	function getnorm(vector, ln) result (res)
		implicit none
    		integer :: i
    		integer, intent(in) 					:: ln
         	double precision, dimension(ln), intent(in)	:: vector   		
        	double precision 						:: res, sm
    		sm = 0.D0
		do i = 1, ln
			sm = sm + vector(i)*vector(i)
		end do
		res = sqrt(sm)
  	end function getnorm
  	
  	
  	subroutine finalUpdateStep(x, v, x_0, v_0, x_fin, v_fin)
  		use problem
		implicit none
		integer :: j
		double precision, dimension(3)		:: x_0, v_0, x_fin, v_fin
		double precision, dimension(3, colnodes)	:: F, x, v
		
		F = 0.D0
		
		do j = 1, colnodes
			F(:,j) = getF(x(:,j), v(:,j))
		end do

		x_fin = x_0
		v_fin = v_0
		do j = 1, colnodes
			x_fin = x_fin + q(j)*v(:,j)
			v_fin = v_fin + q(j)*F(:,j)
		end do      
	end subroutine finalUpdateStep
	
end module sdc_gmres
