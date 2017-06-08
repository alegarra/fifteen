module kinds
  integer, parameter :: single = SELECTED_REAL_KIND( 6, 37 )
  integer, parameter :: double = SELECTED_REAL_KIND( 15, 307 )
  integer, parameter :: extended = SELECTED_REAL_KIND( 18, 4931 )
  integer, parameter :: r4 = SELECTED_REAL_KIND( 6, 37 )
  integer, parameter :: r8 = SELECTED_REAL_KIND( 15, 307 )
  integer, parameter :: r16 = SELECTED_REAL_KIND( 18, 4931 )
  integer, parameter :: i8 = SELECTED_INT_KIND( 8 )
  integer, parameter :: i16 = SELECTED_INT_KIND( 16 )
  ! integer, parameter :: i8 = SELECTED_INT_KIND( 32 )
  ! current precison for hash storage
  integer, parameter :: rh=r8
end module kinds

