module file_mgmt
use config
implicit none
private

!DATA DICTIONARY
real(wp)::A,B,k,q,error,bottom,top,left,right
integer::nx,ny,iter,method


public read_namelist,write_convergance,write_array,A,B,k,q,error,bottom,&
    top,left,right,nx,ny,iter,method

contains
    subroutine read_namelist
        
        !INTERNAL DATA
        character(40)::filename
        integer::i_err
        
        namelist/param/A,B,nx,ny,k,q,iter,error,method
        namelist/boundary/bottom,top,right,left
        
        !Get filename from user
        write(*,*)'Enter the input file name:'
        read(*,*)filename

        !Open the file
        open(unit=10,file=trim(adjustl(filename)),status='old',action='read',&
            & iostat=i_err)

        !Check the open
        if (i_err /= 0) then
            write(*,*)'Error opening file.'
            stop
        end if

        !Read the file
        read(10,param)
        read(10,boundary)

        !Close the file
        close(10)
    end subroutine read_namelist
    
    subroutine write_convergance(iter,error)
        !~~NOTE: Must open and close file outside of subroutine
        
        !INPUT VARIABLES
        integer,intent(in)::iter
        real(wp),intent(in)::error
        
        write(35,*)iter,error
    end subroutine write_convergance
    
    subroutine write_array(A,filename)
        !INPUT VARIABLES
        real(wp),intent(in),dimension(:,:)::A
        character(*)::filename
        
        !INTERNAL VARIABLES
        integer::xdim,ydim,j
        
        xdim=size(A,1)
        ydim=size(A,2)
        
        open(unit=33,file=trim(adjustl(filename)),status='replace',action='write')
        do j=ydim,1,-1
            write(33,*)A(:,j)
        end do
        close(33)        
    end subroutine write_array


end module file_mgmt
