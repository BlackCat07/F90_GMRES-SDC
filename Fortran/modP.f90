module problem
	implicit none 
	integer ::	icount
	double precision, parameter	:: alpha = 47918787.60368732D0
    double precision, parameter	:: R0 = 0.30004580D1, Z0 = 0.30397317D0
    double precision, parameter	:: Er0 = 0.5D5
    double precision, parameter	:: Rad0 = 1.5D0        		
	double precision, dimension(3):: x0, v0
	character(len = 4) :: pname

	contains
	
	
	subroutine InitParticle
	implicit none 
	double precision :: cs, sn
	
#if	    (PTCL==1)

	!Particle: passing
	double precision, dimension(3) :: posCyl = (/ 0.3085263955235203D1,  0.D0, -0.7329976002629984D-1 /)
	double precision, dimension(3) :: velCyl = (/ 0.8141583106593527D6,  0.9313541839057506D6,  0.1793580549387772D7 /)
    pname = "pas_"
	
#elif	(PTCL==2)
   
	!Particle: trapped
	double precision, dimension(3) :: posCyl = (/ 0.2188964117276106D1,  0.D0,  0.8635434778595501D0 /)
	double precision, dimension(3) :: velCyl = (/ 0.2269604314340626D7,  0.2922640610865196D6, -0.3385260666089379D6 /)
    pname = "trp_"
		
#endif
		cs = cos(posCyl(2))
		sn = sin(posCyl(2))    
		x0 = [posCyl(1)*cs, posCyl(1)*sn, posCyl(3) ]
		v0 = [velCyl(1)*cs - velCyl(2)*sn, velCyl(1)*sn + velCyl(2)*cs, velCyl(3)]

	end subroutine InitParticle	
	
	
	function CylinToCartesCoord(Cylind) result (Cartesian)
		implicit none 
		double precision, dimension(3), intent(in)  :: Cylind
		double precision, dimension(3)		        :: Cartesian
		
    		Cartesian(1) = Cylind(1)*cos( Cylind(2) )
    		Cartesian(2) = Cylind(1)*sin( Cylind(2) )
    		Cartesian(3) = Cylind(3)

	end function CylinToCartesCoord
	
	
	subroutine CartesToCylinCoord(x, R, Phi, Z)
		implicit none
		double precision :: R, Phi, Z, x(3)
		
		R = sqrt(x(1)**2 + x(2)**2)
		Phi = atan2(x(2),x(1))
		Z = x(3)
		
	end subroutine CartesToCylinCoord
	
	
!   Output: (r, phi, theta), Input: (R, Phi, Z)	
	function CylinToToroidCoord(Cylind) result (Toroid)
		implicit none 
		double precision, dimension(3), intent(in)  :: Cylind
		double precision, dimension(3)		        :: Toroid
		
        Toroid(3) = atan2(Cylind(3) - Z0, Cylind(1) - R0)
        Toroid(2) = Cylind(2)
        Toroid(1) = (Cylind(1) - R0)/cos(Toroid(3))
        
	end function CylinToToroidCoord
	
	
	function getB(x) result (B)
		implicit none
		double precision, dimension(3), intent(in) :: x
		double precision, dimension(3)		 :: B	
		double precision :: sia_slv, eps_slv, taa_slv, psi_slv, rma_slv
		double precision :: rmi_slv, zmg_slv, RB_phi, R, Phi, Z, xx, yy
		double precision :: B_R, B_Z, B_T	
		icount = icount + 1

		call CartesToCylinCoord(x, R, Phi, Z)

		sia_slv = 0.146387369075D1
		eps_slv = 0.226156682140D0
		taa_slv = 0.143320389205D1
		psi_slv = 0.113333149039D1
		rma_slv = 0.383120489000D1
		rmi_slv = 0.196085203000D1
		zmg_slv = 0.303973168000D0
		RB_phi  =-0.996056843000D1

		xx = 2.D0*(R - rmi_slv)/(rma_slv - rmi_slv) - 1.D0
		yy = Z - zmg_slv
                
		B_R  	= -(2.D0*yy/sia_slv**2)*(1.D0 - 0.25D0*eps_slv**2)*(1.D0 &
			+ taa_slv*eps_slv*xx*(2.D0 + eps_slv*xx))/(psi_slv*R)
              
		B_Z  	= 4.D0*(1.D0 + eps_slv*xx)*(xx - 0.5D0*eps_slv*(1.D0 - xx**2) &
			+ (1.D0 - 0.25D0*eps_slv**2)*yy**2*taa_slv*eps_slv/sia_slv**2)/(psi_slv*R*(rma_slv - rmi_slv))

		B_T  = RB_phi/R
      
		B(1) = B_R*cos(Phi) - B_T*sin(Phi)
		B(2) = B_R*sin(Phi) + B_T*cos(Phi)
		B(3) = B_Z

	end function getB
	
	
	function getE(x) result (E)
		implicit none
		double precision, dimension(3), intent(in)  :: x
		double precision, dimension(3)		        :: E		
		double precision :: Etor(3), Ecyl(3), tor(3), cyl(3)
	    double precision :: R, Phi, Z
	    
 		call CartesToCylinCoord(x, R, Phi, Z)
	    cyl = [R, Phi, Z]
	    tor = CylinToToroidCoord(cyl)
        
!       Toroidal components (r, phi, theta):
        Etor(1) = Er0*tor(1)**2/Rad0**2
        Etor(2) = 0.0
        Etor(3) = 0.0
        
!       Cylindrical components (R, Phi, Z):
        Ecyl(1) = Etor(1)*cos(tor(3))
        Ecyl(2) = 0.0
        Ecyl(3) = Etor(1)*sin(tor(3))

!       Cartesian components (X, Y, Z):
    	E(1) = Ecyl(1)*cos(Phi)
    	E(2) = Ecyl(1)*sin(Phi)
    	E(3) = Ecyl(3)

	end function getE
	
	
	function getF(x,v) result (F)
		implicit none
		double precision, dimension(3), intent(in) :: x, v
		double precision, dimension(3)		 :: F	
		double precision, dimension(3)		 :: cross, Bf

		Bf	    = getB(x)            						
      	cross     =	[v(2) * Bf(3) - v(3) * Bf(2), &
	    			 v(3) * Bf(1) - v(1) * Bf(3), &
	    			 v(1) * Bf(2) - v(2) * Bf(1)]	 
		F =  alpha * (getE(x) + cross)
	end function getF
	
	
	function getEnergy(x,v) result (energy)
		implicit none
		double precision, dimension(3), intent(in)  :: x,v
		double precision :: ekin, epot, energy
		double precision :: tor(3), cyl(3), R, Phi, Z
		 		
		call CartesToCylinCoord(x, R, Phi, Z)
	    cyl = [R, Phi, Z]
	    tor = CylinToToroidCoord(cyl)		
			
        ekin    = 0.5D0*(v(1)*v(1) + v(2)*v(2) + v(3)*v(3))
        epot    = - alpha*Er0*tor(1)**3/Rad0**2/3.D0
        energy  = ekin + epot

	end function getEnergy

end module problem
