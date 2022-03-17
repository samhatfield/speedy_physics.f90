module params
    implicit none

    public

    integer, parameter :: isc = 1

    integer, parameter :: ntrun = 30
    integer, parameter :: mtrun = 30
    integer, parameter :: nlon = 96
    integer, parameter :: nlat = 48
    integer, parameter :: ngp = nlon*nlat

    integer, parameter :: nx = ntrun + 2
    integer, parameter :: mx = mtrun + 1
    integer, parameter :: mxnx = mx*nx
    integer, parameter :: mx2 = 2*mx
    integer, parameter :: ntrun1 = ntrun + 1
    integer, parameter :: nxp = nx + 1
    integer, parameter :: mxp = isc*mtrun + 1
    integer, parameter :: lmax = mxp + nx - 2

    integer, parameter :: nlev = 8

    integer, parameter :: kx2 = 2*nlev
    integer, parameter :: kxm = nlev - 1
    integer, parameter :: kxp = nlev + 1
    integer, parameter :: ntr = 1
end module

