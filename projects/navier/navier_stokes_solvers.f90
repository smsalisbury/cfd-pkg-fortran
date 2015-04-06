module navier_stokes_solvers
use config
use array
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

	subroutine umomentum(mu,mv,u,v,P,density,viscosity,dx,dy,rel,bc_list,store_ap)
		!	-------------------------------------
		!	This subroutine solves the u-momentum
		!	component of the Navier-Stokes Equations.
		!	Uses an SOR scheme
		!	-------------------------------------
		
		!	INPUTS
		real(wp),intent(in)						::	density		!	Fluid density
		real(wp),intent(in)						::	viscosity	!	Fluid viscosity
		real(wp),intent(in)						::	dx,dy		!	Mesh parameters
		real(wp),dimension(:,:),intent(in)		::	v			!	v-velocities (not used)
		real(wp),dimension(:,:),intent(in)		::	P			!	Pressures
		real(wp),dimension(:,:),intent(in)		::	mu,mv		!	Mass flow rates
		real(wp),intent(in)						::	rel			!	Relaxation factor
		character(*),dimension(:),intent(in)	::	bc_list		!	A vector list of bc_types
		
		!	OUTPUTS
		real(wp),dimension(:,:),intent(out)	::	store_ap	!	an array to store the ap values
		
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
		
		integer								::	input_side,output_side
		real(wp)							::	input_m_dot,output_m_dot
		
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
					
					store_ap(i,j) = AP
					
					UN = u(i,j+1)
					UE = u(i+1,j)
					US = u(i,j-1)
					UW = u(i-1,j)
					
					Pe = P(i_count,j_count)
					Pw = P(i_count-1,j_count)
					
					!	-- Handle outlet boundary conditions (if applicable)
					if (j == y_steps .AND. trim(bc_list(1)) == 'outlet') then
						UN = u(i,j)
						u(i,j) = (1.0_wp-rel)*u(i,j) + (1.0_wp/AP)*(AE*UE + AW*UW + AN*UN + AS*US + (Pw-Pe)*dy)
					elseif (i == x_steps .AND. trim(bc_list(2)) == 'outlet') then
						UE = u(i,j)
						u(i,j) = (1.0_wp-rel)*u(i,j) + (1.0_wp/AP)*(AE*UE + AW*UW + AN*UN + AS*US + (Pw-Pe)*dy)
					elseif (j == 0 .AND. trim(bc_list(3)) == 'outlet') then
						US = u(i,j)
						u(i,j) = (1.0_wp-rel)*u(i,j) + (1.0_wp/AP)*(AE*UE + AW*UW + AN*UN + AS*US + (Pw-Pe)*dy)
					elseif (i == 0 .AND. trim(bc_list(4)) == 'outlet') then
						UW = u(i,j)
						u(i,j) = (1.0_wp-rel)*u(i,j) + (1.0_wp/AP)*(AE*UE + AW*UW + AN*UN + AS*US + (Pw-Pe)*dy)
					else
						u(i,j) = (1.0_wp-rel)*u(i,j) + (1.0_wp/AP)*(AE*UE + AW*UW + AN*UN + AS*US + (Pw-Pe)*dy)
					endif
				enddo
			enddo
		enddo
		
		!	SATISFY GLOBAL CONSERVATION
		input_side = array_search_char('velocity',bc_list)
		output_side = array_search_char('output',bc_list)
		
		select case (input_side)
			case (1)
				input_m_dot = sum(mv(:,y_steps))/size(mv(:,y_steps))
			case (2)
				input_m_dot = sum(mu(x_steps,:))/size(mu(x_steps,:))
			case (3)
				input_m_dot = sum(mv(:,1))/size(mv(:,1))
			case (4)
				input_m_dot = sum(mu(1,:))/size(mu(1,:))
			case default
		end select
		
		select case (output_side)
			case (1)
				output_m_dot = sum(mv(:,y_steps))/size(mv(:,y_steps))
				u(:,y_steps) = (input_m_dot/output_m_dot)*u(:,y_steps)
			case (2)
				output_m_dot = sum(mu(x_steps,:))/size(mu(x_steps,:))
				u(x_steps,:) = (input_m_dot/output_m_dot)*u(x_steps,:)
			case (3)
				output_m_dot = sum(mv(:,1))/size(mv(:,1))
				u(:,1) = (input_m_dot/output_m_dot)*u(:,1)
			case (4)
				output_m_dot = sum(mu(1,:))/size(mu(1,:))
				u(1,:) = (input_m_dot/output_m_dot)*u(1,:)
			case default
		end select
	end subroutine umomentum

	subroutine vmomentum(mu,mv,u,v,P,density,viscosity,dx,dy,rel,bc_list,store_ap)
		!	-------------------------------------
		!	This subroutine solves the v-momentum
		!	component of the Navier-Stokes Equations.
		!	Uses an SOR scheme
		!	-------------------------------------
		
		!	INPUTS
		real(wp),intent(in)						::	density		!	Fluid density
		real(wp),intent(in)						::	viscosity	!	Fluid viscosity
		real(wp),intent(in)						::	dx,dy		!	Mesh parameters
		real(wp),dimension(:,:),intent(in)		::	u			!	u-velocities (not used)
		real(wp),dimension(:,:),intent(in)		::	P			!	Pressures
		real(wp),dimension(:,:),intent(in)		::	mu,mv		!	Mass flow rates
		real(wp),intent(in)						::	rel			!	Relaxation factor
		character(*),dimension(:),intent(in)	::	bc_list		!	A vector list of bc_types
		
		!	OUTPUTS
		real(wp),dimension(:,:),intent(out)	::	store_ap	!	an array to store the ap values
		
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
					
					store_ap(i,j) = AP
					
					VN = v(i,j+1)
					VE = v(i+1,j)
					VS = v(i,j-1)
					VW = v(i-1,j)
					
					Pn = P(i_count,j_count)
					Ps = P(i_count,j_count-1)
					
					!	-- Handle outlet boundary conditions (if applicable)
					if (j == y_steps .AND. trim(bc_list(1)) == 'outlet') then
						VN = v(i,j)
						v(i,j) = (1.0_wp-rel)*v(i,j) + (1.0_wp/AP)*(AE*VE + AW*VW + AN*VN + AS*VS + (Ps-Pn)*dx)
					elseif (i == x_steps .AND. trim(bc_list(2)) == 'outlet') then
						VE = v(i,j)
						v(i,j) = (1.0_wp-rel)*v(i,j) + (1.0_wp/AP)*(AE*VE + AW*VW + AN*VN + AS*VS + (Ps-Pn)*dx)
					elseif (j == 0 .AND. trim(bc_list(3)) == 'outlet') then
						VS = v(i,j)
						v(i,j) = (1.0_wp-rel)*v(i,j) + (1.0_wp/AP)*(AE*VE + AW*VW + AN*VN + AS*VS + (Ps-Pn)*dx)
					elseif (i == 0 .AND. trim(bc_list(4)) == 'outlet') then
						VW = v(i,j)
						v(i,j) = (1.0_wp-rel)*v(i,j) + (1.0_wp/AP)*(AE*VE + AW*VW + AN*VN + AS*VS + (Ps-Pn)*dx)
					else
						v(i,j) = (1.0_wp-rel)*v(i,j) + (1.0_wp/AP)*(AE*VE + AW*VW + AN*VN + AS*VS + (Ps-Pn)*dx)
					endif
				enddo
			enddo
		enddo
	end subroutine vmomentum

	subroutine pressure(mu,mv,u,v,P,density,viscosity,dx,dy,rel_c,rel,AP_u,AP_v,continuity)
		!	-------------------------------------
		!	This subroutine solves the pressure
		!	component of the Navier-Stokes Equations.
		!	Uses an SOR scheme
		!	-------------------------------------
		
		!	INPUTS
		real(wp),intent(in)					::	density		!	Fluid density
		real(wp),intent(in)					::	viscosity	!	Fluid viscosity
		real(wp),intent(in)					::	dx,dy		!	Mesh parameters
		real(wp),dimension(:,:),intent(in)	::	mu,mv		!	Mass flow rates
		real(wp),intent(in)					::	rel,rel_c	!	Relaxation factor
		real(wp),dimension(:,:),intent(in)	::	AP_u,AP_v	!	APs from momentum solvers
		
		!	OUTPUTS
		real(wp),dimension(:,:)				::	continuity				!	Source terms
		
		!	INOUTS
		real(wp),dimension(:,:),intent(inout)	::	P			!	Pressures
		real(wp),dimension(:,:),intent(inout)	::	u			!	u-velocities
		real(wp),dimension(:,:),intent(inout)	::	v			!	v-velocities
		
		!	INTERNAL
		integer								::	i,j,i_offset,j_offset
		integer								::	k
		integer,parameter					::	max_iter = 10
		integer								::	x_steps,y_steps
		real(wp)							::	an_v,ae_u,as_v,aw_u		!	Relative velocity coefficients
		real(wp)							::	PN,PS,PE,PW,PP			!	Local pressures
		real(wp)							::	vn,ue,vs,uw				!	Local velocities
		real(wp)							::	AN,AE,AS,AW,AP
		real(wp)							::	tmp
		real(wp),dimension(:,:),allocatable	::	P_c						!	Pressure correction
		real(wp),dimension(:,:),allocatable	::	S						!	Source terms
		
		!	CALCULATE NEW U_MOMENTUM
		x_steps = size(P,1)
		y_steps = size(P,2)
		
		allocate(P_c(x_steps,y_steps))
		P_c = 0.0_wp
		allocate(S(x_steps,y_steps))
		
		do j = 1,y_steps
			do i = 1,x_steps
				!	Counter offsets for u and v-momentum
				i_offset = i+1
				j_offset = j+1
			
				!	Relative velocity coefficients
				an_v = AP_v(i_offset,j_offset+1)
				ae_u = AP_u(i_offset+1,j_offset)
				as_v = AP_v(i_offset,j_offset)
				aw_u = AP_u(i_offset,j_offset)
				
				!	Pressure correction
				AN = (density*dx**2)/an_v
				AE = (density*dy**2)/ae_u
				AS = (density*dx**2)/as_v
				AW = (density*dy**2)/aw_u
				if (abs(an_v) < 1.0E-5_wp) AN = 0.0_wp
				if (abs(ae_u) < 1.0E-5_wp) AE = 0.0_wp
				if (abs(as_v) < 1.0E-5_wp) AS = 0.0_wp
				if (abs(aw_u) < 1.0E-5_wp) AW = 0.0_wp
				vn = v(i_offset,j_offset+1)
				ue = u(i_offset+1,j_offset)
				vs = v(i_offset,j_offset)
				uw = u(i_offset,j_offset)
				
				S(i,j) = density*((ue-uw)*dy + (vn-vs)*dx)
				do k = 1,max_iter
					
					! 	Local pressures
					if (i .EQ. 1) then
						PW = P_c(i,j)
						AW = 0.0_wp
					else
						PW = P_c(i-1,j)
					endif
					if (j .EQ. 1) then
						PS = P_c(i,j)
						AS = 0.0_wp
					else
						PS = P_c(i,j-1)
					endif
					if (i .EQ. x_steps) then
						PE = P_c(i,j)
						AE = 0.0_wp
					else
						PE = P_c(i+1,j)
					endif
					if (j .EQ. y_steps) then
						PN = P_c(i,j)
						AN = 0.0_wp
					else
						PN = P_c(i,j+1)
					endif
					AP = (AN + AE + AS + AW)
					
					
					if (abs(AP) > 1.0E-5_wp) then
						P_c(i,j) = P_c(i,j) + (rel_c/AP)*(AN*PN + AE*PE + AS*PS + AW*PW - S(i,j) - AP*P_c(i,j))
					endif
				enddo
				PP = P_c(i,j)
				
				P(i,j) = P(i,j) + rel*P_c(i,j)
				
			enddo
		enddo
		
		!	UPDATE U-MOMENTUM
		do i=2,x_steps
			do j=1,y_steps
				!	Counter offsets for u and v-momentum
				i_offset = i+1
				j_offset = j+1
				
				!	Velocity Corrections
				if (abs(ae_u) > 1.0E-5_wp) u(i_offset+1,j_offset) = u(i_offset+1,j_offset) + (dy/ae_u)*(PP - PE)
				if (abs(aw_u) > 1.0E-5_wp) u(i_offset,j_offset) = u(i_offset,j_offset) + (dy/aw_u)*(PW - PP)
			enddo
		enddo
		
		!	UPDATE V-MOMENTUM
		do i=1,x_steps
			do j=2,y_steps
				!	Counter offsets for u and v-momentum
				i_offset = i+1
				j_offset = j+1
				
				!	Velocity Corrections
				if (abs(an_v) > 1.0E-5_wp) v(i_offset,j_offset+1) = v(i_offset,j_offset+1) + (dx/an_v)*(PP - PN)
				if (abs(as_v) > 1.0E-5_wp) v(i_offset,j_offset) = v(i_offset,j_offset) + (dx/as_v)*(PS - PP)
			enddo
		enddo
		
		continuity = S
		
		deallocate(P_c)
	end subroutine pressure




















end module navier_stokes_solvers