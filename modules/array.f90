module array
use config

contains
	function array_search_char(needle,haystack)
		!	INPUTS
		character(*)				::	needle
		character(*),dimension(:)	::	haystack
		
		!	OUTPUTS
		integer						::	array_search_char
		
		!	INTERNAL
		integer						::	i
		
		!	FIND needle IN haystack
		array_search_char = -1
		
		do i=1,size(haystack)
			if (trim(needle) .EQ. trim(haystack(i))) then
				array_search_char = i
				exit
			endif
		enddo	
	end function array_search_char
	
	
end module array