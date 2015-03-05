module navier_stokes_solvers
use config
implicit none

contains
	subroutine mass_flow_rates(mu,mv,u,v,P,density,viscosity,dx,dy)
		!	-------------------------------------
		!	This subroutine calculates the mass
		!	flow rates of the system
		!	-------------------------------------
		
		!	INPUTS
		real(wp),intent(in)					::	density		!	Fluid density
		real(wp),intent(in)					::	viscosity	!	Fluid viscosity (not used)
		real(wp),intent(in)					::	dx,dy		!	Mesh parameters
		real(wp),dimension(:,:),intent(in)	::	u,v			!	Velocities
		real(wp),dimension(:,:),intent(in)	::	P			!	Pressures (not used)
		
		!	OUTPUTS
		real(wp),dimension(:,:),intent(out)	::	mu,mv		!	Mass flow rates
		
		!	INOUTS
		
		!	INTERNAL
		integer								::	i,j,i_count,j_count
		integer								::	x_steps,y_steps
		
		!	CALCULATE MASS FLOW RATES
		!	--	u mass flow rates
		x_steps = size(u,1)-1
		y_steps = size(u,2)-1
		do j_count=0,y_steps
			do i_count=0,x_steps
				i = i_count + 1
				j = j_count + 1
				mu(i,j) = density * dx * u(i,j)
			enddo
		enddo
		!	--	v mass flow rates
		x_steps = size(v,1)-1
		y_steps = size(v,2)-1
		y_steps = size(v,2)-1
		do j_count=0,y_steps
			do i_count=0,x_steps
				i = i_count + 1
				j = j_count + 1
				mv(i,j) = density * dx * v(i,j)
			enddo
		enddo
	end subroutine mass_flow_rates

	subroutine umomentum(mu,mv,u,v,P,density,viscosity,dx,dy,rel)
		!	-------------------------------------
		!	This subroutine solves the u-momentum
		!	component of the Navier-Stokes Equations.
		!	Uses an SOR scheme
		!	-------------------------------------
		
		!	INPUTS
		real(wp),intent(in)					::	density		!	Fluid density
		real(wp),intent(in)					::	viscosity	!	Fluid viscosity
		real(wp),intent(in)					::	dx,dy		!	Mesh parameters
		real(wp),dimension(:,:),intent(in)	::	v			!	v-velocities (not used)
		real(wp),dimension(:,:),intent(in)	::	P			!	Pressures
		real(wp),dimension(:,:),intent(in)	::	mu,mv		!	Mass flow rates
		real(wp),intent(in)					::	rel			!	Relaxation factor
		
		!	OUTPUTS
		
		!	INOUTS
		real(wp),dimension(:,:),intent(inout)	::	u			!	u-velocities
		
		!	INTERNAL
		integer								::	i,j,i_count,j_count
		integer								::	k
		integer,parameter					::	max_iter = 10
		integer								::	x_steps,y_steps
		real(wp)							::	AN,AE,AS,AW,AP
		real(wp)							::	mn_l,me_l,ms_l,mw_l		!	Local mass flow rates
		real(wp)							::	Pe,Pw					!	Local pressures
		real(wp)							::	UN,UE,US,UW				!	Local velocities
		
		!	CALCULATE NEW U_MOMENTUM
		x_steps = size(u,1)-1
		y_steps = size(u,2)-1
		do k = 1,max_iter
			do j_count = 1,y_steps-1
				do i_count = 2,x_steps-1
					i = i_count+1
					j = j_count+1
					
					if (j_count >= y_steps) then
						mn_l = 0.0_wp
					else
						mn_l = 0.5_wp * (mv(i,j+1) + mv(i-1,j+1))
					endif
					if (i_count >= x_steps) then
						me_l = 0.0_wp
					else
						me_l = 0.5_wp * (mu(i,j) + mu(i+1,j))
					endif
					
					ms_l = 0.5_wp * (mv(i,j) + mv(i-1,j))
					mw_l = 0.5_wp * (mu(i,j) + mu(i-1,j))
					
					AN = max(-mn_l,0.0_wp) + viscosity*dx/dy
					AE = max(-me_l,0.0_wp) + viscosity*dy/dx
					AS = max(ms_l,0.0_wp) + viscosity*dy/dx
					AW = max(mw_l,0.0_wp) + viscosity*dy/dx
					AP = (AN + AE + AS + AW)/rel
					
					UN = u(i,j+1)
					UE = u(i+1,j)
					US = u(i,j-1)
					UW = u(i-1,j)
					
					Pe = P(i_count,j_count)
					Pw = P(i_count-1,j_count)
					
					u(i,j) = (1.0_wp-rel)*u(i,j) + (1.0_wp/AP)*(AE*UE + AW*UW + AN*UN + AS*US + (Pw-Pe)*dy)
				enddo
			enddo
		enddo
	end subroutine umomentum

	subroutine vmomentum(mu,mv,u,v,P,density,viscosity,dx,dy,rel)
		!	-------------------------------------
		!	This subroutine solves the v-momentum
		!	component of the Navier-Stokes Equations.
		!	Uses an SOR scheme
		!	-------------------------------------
		
		!	INPUTS
		real(wp),intent(in)					::	density		!	Fluid density
		real(wp),intent(in)					::	viscosity	!	Fluid viscosity
		real(wp),intent(in)					::	dx,dy		!	Mesh parameters
		real(wp),dimension(:,:),intent(in)	::	u			!	u-velocities (not used)
		real(wp),dimension(:,:),intent(in)	::	P			!	Pressures
		real(wp),dimension(:,:),intent(in)	::	mu,mv		!	Mass flow rates
		real(wp),intent(in)					::	rel			!	Relaxation factor
		
		!	OUTPUTS
		
		!	INOUTS
		real(wp),dimension(:,:),intent(inout)	::	v			!	v-velocities
		
		!	INTERNAL
		integer								::	i,j,i_count,j_count
		integer								::	k
		integer,parameter					::	max_iter = 10
		integer								::	x_steps,y_steps
		real(wp)							::	AN,AE,AS,AW,AP
		real(wp)							::	mn_l,me_l,ms_l,mw_l		!	Local mass flow rates
		real(wp)							::	Pn,Ps					!	Local pressures
		real(wp)							::	VN,VE,VS,VW				!	Local velocities
		
		!	CALCULATE NEW U_MOMENTUM
		x_steps = size(v,1)-1
		y_steps = size(v,2)-1
		do k = 1,max_iter
			do j_count = 2,y_steps-1
				do i_count = 1,x_steps-1
					i = i_count+1
					j = j_count+1
					
					if (j_count >= y_steps) then
						mn_l = 0.0_wp
					else
						mn_l = 0.5_wp * (mv(i,j) + mv(i,j+1))
					endif
					if (i_count >= x_steps) then
						me_l = 0.0_wp
					else
						me_l = 0.5_wp * (mu(i+1,j-1) + mu(i+1,j))
					endif
					
					ms_l = 0.5_wp * (mv(i,j) + mv(i,j-1))
					mw_l = 0.5_wp * (mu(i,j-1) + mu(i,j))
					
					AN = max(-mn_l,0.0_wp) + viscosity*dx/dy
					AE = max(-me_l,0.0_wp) + viscosity*dy/dx
					AS = max(ms_l,0.0_wp) + viscosity*dy/dx
					AW = max(mw_l,0.0_wp) + viscosity*dy/dx
					AP = (AN + AE + AS + AW)/rel
					
					VN = v(i,j+1)
					VE = v(i+1,j)
					VS = v(i,j-1)
					VW = v(i-1,j)
					
					Pn = P(i_count,j_count)
					Ps = P(i_count,j_count-1)
					
					v(i,j) = (1.0_wp-rel)*v(i,j) + (1.0_wp/AP)*(AE*VE + AW*VW + AN*VN + AS*VS + (Ps-Pn)*dx)
				enddo
			enddo
		enddo
	end subroutine vmomentum






















end module navier_stokes_solvers