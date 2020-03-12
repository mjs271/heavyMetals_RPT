module FD_hMetal_mod
implicit none

! these can be used to get the component list and selected output headers
! from the commented code blocks near line 175 in dolo_ADRE.f90
type component_list
    character(:), allocatable :: comp
end type
type selectout_list
    character(:), allocatable :: head
end type

! global constants
double precision, parameter :: pi = 4.0d0 * atan(1.0d0)
integer, parameter          :: sp = kind(1.0), dp = kind(1.d0)
integer, parameter          :: nspec = 20

contains

! advection subroutine, using explicit upwinding
subroutine advect(mat, v, dx, dt, ncell)
    double precision, intent(inout) :: mat(:, :)
    double precision, intent(in   ) :: v(:), dx, dt
    integer,          intent(in   ) :: ncell
    double precision                :: vmat(ncell - 1, nspec)
    integer                         :: i

    do i = 1, nspec
        vmat(:, i) = v(1 : ncell - 1)
    enddo

    mat(2 : ncell, :) = mat(2 : ncell, :) - ((vmat * dt)/dx) *&
                        (mat(2 : ncell, :) - mat(1 : ncell - 1, :))
end subroutine advect

! diffusion subroutine using explicit Euler
subroutine diffuse(mat, D, dx, dt, phi, ncell)
    double precision, intent(inout) :: mat(:, :)
    double precision, intent(in   ) :: D(:), dx, dt, phi
    integer,          intent(in   ) :: ncell
    double precision                :: Dmat(ncell - 2, nspec), temp(nspec)
    integer                         :: i

    ! stop concentrations from diffusing in across upper boundary by eliminating
        ! concentration gradient between cells end and end - 1, then replace after
    temp = mat(ncell, :)
    mat(ncell, :) = mat(ncell - 1, :)
    do i = 1, nspec
        Dmat(:, i) = D(2 : ncell - 1)
    enddo

    mat(2 : ncell - 1, :) = mat(2 : ncell - 1, :) + ((Dmat * dt)/(phi * dx**2)) *&
                            (mat(3 : ncell, :) - 2 * mat(2 : ncell - 1, :)&
                                + mat(1 : ncell - 2, :))
end subroutine diffuse

end module FD_hMetal_mod
