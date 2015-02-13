module config
implicit none
private

integer,parameter::sp=selected_real_kind(6)
integer,parameter::dp=selected_real_kind(15)
integer,parameter::ep=selected_real_kind(18)
integer,parameter::qp=selected_real_kind(30)
integer,parameter::wp=dp

public wp,disp_wp

contains
    subroutine disp_wp
        write(*,*)wp
    end subroutine
end module config
