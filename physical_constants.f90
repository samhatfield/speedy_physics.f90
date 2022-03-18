!  For storing and initializing physical constants.
module physical_constants

    implicit none

    public

    real, parameter :: grav   = 9.81 !! Gravitational acceleration (m/s/s)

    ! Physical constants for thermodynamics
    real, parameter :: p0   = 1.e+5   !! Reference pressure (Pa)
    real, parameter :: cp   = 1004.0  !! Specific heat at constant pressure (J/K/kg)
    real, parameter :: alhc = 2501.0  !! Latent heat of condensation, in J/g for consistency with
end module
