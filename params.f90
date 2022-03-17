module params
    implicit none

    public

    integer, parameter :: isc = 1

    integer, parameter :: ntrun = 30
    integer, parameter :: mtrun = 30
    integer, parameter :: ix = 96
    integer, parameter :: iy = 24

    integer, parameter :: nx = ntrun + 2
    integer, parameter :: mx = mtrun + 1
    integer, parameter :: mxnx = mx*nx
    integer, parameter :: mx2 = 2*mx
    integer, parameter :: il = 2*iy
    integer, parameter :: ntrun1 = ntrun + 1
    integer, parameter :: nxp = nx + 1
    integer, parameter :: mxp = isc*mtrun + 1
    integer, parameter :: lmax = mxp + nx - 2

    integer, parameter :: kx = 8

    integer, parameter :: kx2 = 2*kx
    integer, parameter :: kxm = kx - 1
    integer, parameter :: kxp = kx + 1
    integer, parameter :: ntr = 1
end module

