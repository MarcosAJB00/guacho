!=================================

 module network

   use parameters, only : n_spec
   use exoplanet,  only : fHII, fHeII, alpha
   implicit none

   ! number of equilibrium equations
   integer, parameter :: n_equi = 3

   ! number of non-equilibrium equations
   integer, parameter :: n_nequ = n_spec - n_equi

   ! number of total elements
   integer, parameter :: n_elem = 2

   ! indexes of the different species
   integer, parameter :: iHI    = 1
   integer, parameter :: iHII   = 2
   integer, parameter :: iHeIS  = 3
   integer, parameter :: iHeIM  = 4
   integer, parameter :: iHeII  = 5
   integer, parameter :: ie     = 6

   ! indexes of the equilibrium species
   integer, parameter :: iH  = 1
   integer, parameter :: iHe = 2

   ! number of reaction rates
   integer, parameter :: n_reac = 18

   ! indexes of the different rates
   integer, parameter :: ichi      = 1
   integer, parameter :: iche13a   = 2
   integer, parameter :: iche31a   = 3
   integer, parameter :: iche31b   = 4
   integer, parameter :: iche31    = 5
   integer, parameter :: icheIS    = 6
   integer, parameter :: icheIM    = 7
   integer, parameter :: iahii_b   = 8
   integer, parameter :: iaheii_1  = 9
   integer, parameter :: iaheii_a  = 10
   integer, parameter :: iaheii_b  = 11
   integer, parameter :: iaheiim_b = 12
   integer, parameter :: iphiHI    = 13
   integer, parameter :: iphiHeIS  = 14
   integer, parameter :: iphiHeIM  = 15
   integer, parameter :: irheii    = 16
   integer, parameter :: iceheIS   = 17
   integer, parameter :: iceheII   = 18

   ! names of the equation to be solved
   character(len=10)  :: name_eqs(n_spec)

 contains

   !=======

   subroutine derv(y,rate,dydt,y0)

     implicit none
     real (kind=8), intent(in)  ::   y0(n_elem)
     real (kind=8), intent(in)  ::    y(n_spec)
     real (kind=8), intent(out) :: dydt(n_spec)
     real (kind=8), intent(in)  :: rate(n_reac)
     real (kind=8)              :: yy, p
     !cross section coeficitens value at 24.6 eV
     real (kind=8),parameter    :: ahI_246=1.24e-18, aheIS_246=2.42e-19, aheIM_246=4.26e-19

     !fraction of photons from helium recombination to ground state that can ionize hydrogen (ostebrook)
     yy  = ahI_246*y(iHI) / (y(iHI)*ahI_246 + y(iHeIS)*aheIS_246 + y(iHeIM)*aheIM_246)
     !probability  of photons that can ionize hydrogen from recombination to excited levels of helium
     p   = 0.67

     dydt(iHI) = - rate(ichi)    * y(iHI)  * y(ie)                             &
                 + rate(iahii_b) * y(iHII) * y(ie)                             &
                 - rate(iphiHI)  * y(iHI)                                      &
                 - rate(iaheii_b)* y(iHeII)* y(ie)*p                           &
                 - rate(iaheii_1)* y(iHeII)* y(ie)*yy

     dydt(iHeIS)= + rate(irheii)  * y(iHeIM)                                   &
                  + rate(iche31)  * y(iHI)* y(iHeIM)                           &
                  - rate(iche13a) * y(ie )* y(iHeIS)                           &
                  + rate(iche31a) * y(ie )* y(iHeIM)                           &
                  + rate(iche31b) * y(ie )* y(iHeIM)                           &
                  - rate(iaheii_1)* y(ie )* y(iHeII) * (1.-yy)                 &
                  + rate(iaheii_a)* y(ie )* y(iHeII)                           &
                  - rate(iphiHeIS)* y(iHeIS )                                  &
                  - rate(icheIS)  * y(iHeIS)*y(ie)                             &
                  - rate(iceheIS) * y(iHI)*y(iHeII)                            &
                  + rate(iceheII) * y(iHII)*y(iHeIS)

     dydt(iHeIM)= - rate(irheii)   * y(iHeIM)                                  &
                  - rate(iche31)   * y(iHI)* y(iHeIM)                          &
                  + rate(iche13a)  * y(ie )* y(iHeIS)                          &
                  - rate(iche31a)  * y(ie )* y(iHeIM)                          &
                  - rate(iche31b)  * y(ie )* y(iHeIM)                          &
                  + rate(iaheiim_b)* y(ie )* y(iHeII)                          &
                  - rate(iphiHeIM) * y(iHeIM )                                 &
                  - rate(icheIM)   * y(ie)*y(iHeIM)

     !  "conservation" equations

     dydt(iHII) =  -y0(iH) + y(iHI) + y(iHII)

     dydt(iHeII)=  -y0(iHe) + y(iHeIS) + y(iHeIM) + y(iHeII)

     dydt(ie)   =  -y(ie) + y(iHII) + y(iHeII)


   end subroutine derv

   !=======

   subroutine get_jacobian(y,jacobian,rate)

     implicit none
     real (kind=8), intent(in)  :: y(n_spec)
     real (kind=8), intent(out) :: jacobian(n_spec,n_spec)
     real (kind=8), intent(in)  :: rate(n_reac)
     real (kind=8)              :: dydhI, dydheIS, dydheIM
     real (kind=8)              :: yy, p
     !cross section coeficitens value at 24.6 eV
     real (kind=8),parameter    :: ahI_246=1.24e-18, aheIS_246=2.42e-19, aheIM_246=4.26e-19

     !fraction of photons from helium recombination to ground state that can ionize hydrogen (ostebrook)
     yy  = ahI_246*y(iHI) / (y(iHI)*ahI_246 + y(iHeIS)*aheIS_246 + y(iHeIM)*aheIM_246)
     !probability  of photons that can ionize hydrogen from recombination to excited levels of helium
     p   = 0.67

     dydhI   = (ahI_246*(y(iHI)*ahI_246 + y(iHeIS)*aheIS_246 + y(iHeIM)*aheIM_246) - ahI_246**2*y(iHI)) &
             / (y(iHI)*ahI_246 + y(iHeIS)*aheIS_246 + y(iHeIM)*aheIM_246)**2

     dydheIS = (-ahI_246*y(iHI)*aheIS_246) &
             / (y(iHI)*ahI_246 + y(iHeIS)*aheIS_246 + y(iHeIM)*aheIM_246)**2

     dydheIM = (-ahI_246*y(iHI)*aheIM_246) &
             / (y(iHI)*ahI_246 + y(iHeIS)*aheIS_246 + y(iHeIM)*aheIM_246)**2

    jacobian(iHI,iHI)      = - rate(ichi)*y(ie)                                &
                             - rate(iphiHI)                                    &
                             - rate(iaheii_1)*y(ie)*y(iHeII)*dydhI
    jacobian(iHI,iHII)     = + rate(iahii_b)*y(ie)
    jacobian(iHI,iHeIS)    = - rate(iaheii_1)*y(ie)*y(iHeII)*dydheIS
    jacobian(iHI,iHeIM)    = - rate(iaheii_1)*y(ie)*y(iHeII)*dydheIM
    jacobian(iHI,iHeII)    = - rate(iaheii_1)*y(ie)*yy                         &
                             - rate(iaheii_b)*y(ie)*p
    jacobian(iHI,ie)       = - rate(ichi)*y(iHI)                               &
                             + rate(iahii_b)*y(iHII)                           &
                             - rate(iaheii_1)*y(iHeII)*yy                      &
                             - rate(iaheii_b)*y(iHeII)*p

    jacobian(iHII,iHI)     =   1.
    jacobian(iHII,iHII)    =   1.
    jacobian(iHII,iHeIS)   =   0.
    jacobian(iHII,iHeIM)   =   0.
    jacobian(iHII,iHeII)   =   0.
    jacobian(iHII,ie)      =   0.

    jacobian(iHeIS,iHI)    = + rate(iaheii_1)*y(ie)*y(iHeII)*dydhI             &
                             + rate(iche31)*y(iHeIM)                           &
                             - rate(iceheIS)*y(iHeII)
    jacobian(iHeIS,iHII)   =   rate(iceheII)*y(iHeIS)
    jacobian(iHeIS,iHeIS)  = - rate(iche13a)*y(ie)                             &
                             - rate(iphiHeIS)                                  &
                             - rate(icheIS)*y(ie)                              &
                             + rate(iaheii_1)*y(ie)*y(iHeII)*dydheIS           &
                             + rate(iceheII)*y(iHII)

    jacobian(iHeIS,iHeIM)  = + rate(irheii)                                    &
                             + rate(iche31)*y(iHI)                             &
                             + rate(iche31a)*y(ie)                             &
                             + rate(iche31b)*y(ie)                             &
                             + rate(iaheii_1)*y(ie)*y(iHeII)*dydheIM

    jacobian(iHeIS,iHeII)  = + rate(iaheii_a)*y(ie)                            &
                             - rate(iaheii_1)*y(ie)*(1 - yy)                   &
                             - rate(iceheIS)*y(iHI)
    jacobian(iHeIS,ie)     = - rate(iche13a)*y(iHeIS)                          &
                             + rate(iche31a)*y(iHeIM)                          &
                             + rate(iche31b)*y(iHeIM)                          &
                             + rate(iaheii_a)*y(iHeII)                         &
                             - rate(iaheii_1)*y(iHeII)*(1 - yy)                &
                             - rate(icheIS) * y(iHeIS)

    jacobian(iHeIM,iHI)    = - rate(iche31)*y(iHeIM)
    jacobian(iHeIM,iHII)   =   0.
    jacobian(iHeIM,iHeIS)  = + rate(iche13a)*y(ie)
    jacobian(iHeIM,iHeIM)  = - rate(irheii)                                    &
                             - rate(iche31)*y(iHI)                             &
                             - rate(iche31a)*y(ie)                             &
                             - rate(iche31b)*y(ie)                             &
                             - rate(iphiHeIM)                                  &
                             - rate(icheIM) * y(ie)
    jacobian(iHeIM,iHeII)  = + rate(iaheiim_b)*y(ie)
    jacobian(iHeIM,ie)     = + rate(iche13a)*y(iHeIS)                          &
                             - rate(iche31a)*y(iHeIM)                          &
                             - rate(iche31b)*y(iHeIM)                          &
                             + rate(iaheiim_b)*y(iHeII)                        &
                             - rate(icheIM)*y(iHeIM)

    jacobian(iHeII,iHI)    =   0.
    jacobian(iHeII,iHII)   =   0.
    jacobian(iHeII,iHeIS)  =   1.
    jacobian(iHeII,iHeIM)  =   1.
    jacobian(iHeII,iHeII)  =   1.
    jacobian(iHeII,ie)     =   0.

    jacobian(ie,iHI)       =   0.
    jacobian(ie,iHII)      =   1.
    jacobian(ie,iHeIS)     =   0.
    jacobian(ie,iHeIM)     =   0.
    jacobian(ie,iHeII)     =   1.
    jacobian(ie,ie)        =  -1.

   end subroutine get_jacobian

   !=======================================================================

   subroutine get_reaction_rates(rate,T,phHI, phHeIS, phHeIM)
     implicit none
     real (kind=8), intent(in) :: T, phHI, phHeIS, phHeIM
     real (kind=8),intent(out) :: rate(n_reac)
     real (kind=8) :: gamma13, gamma31a, gamma31b,logT

     ! fitted values for gammaij from bray2000 only valid 3.75<logT<5.75 [K]
     logT = log10(T)
     if (logT >= 3.75 .and. logT<=5.75) then
       gamma13  = - 1.5699 + 1.303*log10(T) - 0.388*log10(T)**2 &
                    + 0.0516*log10(T)**3 - 0.0026*log10(T)**4
       gamma31a = - 209.42 + 170.03*log10(T)- 50.073*log10(T)**2  &
                    + 6.42*log10(T)**3 - 0.3045*log10(T)**4
       gamma31b = - 36.477 + 22.988*log10(T) - 4.5935*log10(T)**2 &
                    + 0.2968*log10(T)**3
     else
       gamma13  = 0.
       gamma31a = 0.
       gamma31b = 0.
     endif
     ! collisional ionisation rates for H (raga) and He (Lampon 2020)
     rate(ichi   ) = 5.830e-11*sqrt(T)*exp(-157828.0/T)    !colissional ioniz. of HI [cm3/s]
     rate(iche13a) = 2.10e-8*sqrt(13.6/(8.617333262145e-5*T))*exp(-19.81/(8.617333262145e-5*T))*gamma13    !colissional ioniz. of HeI [cm3/s]
     rate(iche31a) = 2.10e-8*sqrt(13.6/(8.617333262145e-5*T))*exp( -0.80/(8.617333262145e-5*T))*gamma31a/3 !colissional ioniz. of HeII [cm3/s]
     rate(iche31b) = 2.10e-8*sqrt(13.6/(8.617333262145e-5*T))*exp( -1.40/(8.617333262145e-5*T))*gamma31b/3 ![cm3/s]
     rate(iche31)  = 5.e-10                         !HeI colisions with HI
     rate(icheIS)  = 2.379e-11*sqrt(T)*exp(-285335.4/T)     !colisional ionisation of helium 11S (cen 1992)
     rate(icheIM)  = 8.335e-10*sqrt(T)*exp(-55338.0/T)      !colisional ionisation of helium 23S (Black 1981)

     ! recombination rates for He case B from Benjamin+1999
     rate(iahii_b)   = 2.59e-13*(1.0e4/T)**0.700     !recombination of HII [cm3/s]
     rate(iaheii_a)  = 4.27e-13*(1.0e4/T)**0.678     !recombination of HeII case A tot [cm3/s]
     rate(iaheii_b)  = 2.72e-13*(1.0e4/T)**0.789     !recombination of HeII case B tot [cm3/s]
     rate(iaheii_1)  = 1.54e-13*(1.0e4/T)**0.486     !recombination of HeII to ground only [cm3/s]
     rate(iaheiim_b) = 2.10e-13*(1.0e4/T)**0.778     !recombination of HeII to 2s3 case B [cm3/s]

     rate(irheii ) = 1.272e-4           !radiative transition He(23s) to He(11s)

     !charge exchange (okclopick, lampon)
     rate(iceheIS) = 1.25e-15*(300./T)**(-0.25)
     rate(iceheII) = 1.75e-11*(300./T)**0.75*exp(-128000.0/T)

     rate(iphiHI  ) = phHI             !photoionization of HI
     rate(iphiHeIS) = phHeIS           !photoionization of HeI from ground
     rate(iphiHeIM) = phHeIM           !photoionization of HeI from triplet

   end subroutine get_reaction_rates

   !=======================================================================

   subroutine nr_init(y,y0)
     implicit none
     real, intent(out) :: y(n_spec)
     real, intent(in ) :: y0(n_elem)

     y(iHI  ) = 0.5*y0(iH)!(1. - fHII)  * y0(iH )
     y(iHII ) = 0.5*y0(iH)!fHII         * y0(iH )
     y(iHeIS) = 0.0 !(1. - fHeII) * y0(iHe)
     y(iHeIM) = 0.0 !(1. - fHeII) * y0(iHe)/(1+alpha)
     y(iHeII) = 0.0 !fHeII        * y0(iHe)
     y(ie   ) = y(iHII) + y(iHeII)

     return
   end subroutine nr_init

   !=======================================================================

   logical function check_no_conservation(y,y0_in)
     implicit none
     real, intent(in)  :: y(n_spec)
     real, intent(in ) :: y0_in  (n_elem)
     real              :: y0_calc(n_elem)
     integer           :: i

     check_no_conservation = .false.

     y0_calc(iH )= y(iHI ) + y(iHII )
     y0_calc(iHe)= y(iHeIS) + y(iHeIM) + y(iHeII)

     do i = 1, n_elem
       if (y0_calc(i) < 1.01*y0_in(i)) check_no_conservation = .true.
       if (y0_calc(i) > 1.01*y0_in(i)) check_no_conservation = .true.
     end do

   end function check_no_conservation

   !=======================================================================

subroutine gsub(y, q, d, time, T, phHI, phHeIS, phHeIM)

    implicit none
    real (kind=8), intent(in)  :: y(n_spec)
    real (kind=8), intent(out) :: q(n_spec)     ! Producción
    real (kind=8), intent(out) :: d(n_spec)     ! Tasa de pérdida
    real (kind=8), intent(in)  :: time, T, phHI, phHeIS, phHeIM
    real (kind=8)              :: rate(n_reac)
    real (kind=8)              :: yy, prob_p
    !cross section coeficitens value at 24.6 eV
    real (kind=8),parameter    :: ahI_246=1.24e-18, aheIS_246=2.42e-19, aheIM_246=4.26e-19

    !fraction of photons from helium recombination to ground state that can ionize hydrogen (ostebrook)
    yy  = ahI_246*y(iHI) / (y(iHI)*ahI_246 + y(iHeIS)*aheIS_246 + y(iHeIM)*aheIM_246)
    !probability  of photons that can ionize hydrogen from recombination to excited levels of helium
    prob_p   = 0.67

    !  Obtain rates as function of current temperature
    call get_reaction_rates(rate,T, phHI, phHeIS, phHeIM)

    ! ================= HI =================
    ! Producción
    q(iHI) = + rate(iahii_b) * y(iHII) * y(ie)      &
             + rate(iceheII) * y(iHII) *y(iHeIS)

    ! Pérdida d(i) = p(i) * y(i)
    d(iHI) = + rate(ichi)    * y(iHI)  * y(ie)           &
             + rate(iphiHI)  * y(iHI)                    &
             + rate(iaheii_b)* y(iHeII)* y(ie)*prob_p    &
             + rate(iaheii_1)* y(iHeII)* y(ie)*yy        &
             + rate(iceheIS) * y(iHI)  *y(iHeII)


    ! ================= HII ==================
    ! Producción
    q(iHII) = + rate(ichi)    * y(iHI)  * y(ie)         &
              + rate(iphiHI)  * y(iHI)                  &
              + rate(iaheii_b)* y(iHeII)* y(ie)*prob_p  &
              + rate(iaheii_1)* y(iHeII)* y(ie)*yy      &
              + rate(iceheIS) * y(iHI)  *y(iHeII)

    ! Pérdida
    d(iHII) = + rate(iahii_b) * y(iHII) * y(ie)     &
              + rate(iceheII) * y(iHII) *y(iHeIS)


    ! ================= HeIS =================
    ! Producción
    q(iHeIS) = + rate(irheii)  * y(iHeIM)           &
               + rate(iche31)  * y(iHI)* y(iHeIM)   &
               + rate(iche31a) * y(ie )* y(iHeIM)   &
               + rate(iche31b) * y(ie )* y(iHeIM)   &
               + rate(iaheii_a)* y(ie )* y(iHeII)   &
               + rate(iceheIS) * y(iHI)* y(iHeII)
    ! Pérdida
    d(iHeIS) = + rate(iche13a) * y(ie )* y(iHeIS)    &
               + rate(iphiHeIS)* y(iHeIS )           &
               + rate(icheIS)  * y(iHeIS)*y(ie)      &
               + rate(iceheII) * y(iHII)*y(iHeIS)    &
               + rate(iaheii_1)* y(ie )* y(iHeII) * (1.-yy)

    ! ================ HeIM =================
    ! Producción
    q(iHeIM) = + rate(iche13a)  * y(ie )* y(iHeIS)   &
               + rate(iaheiim_b)* y(ie )* y(iHeII)


    ! Pérdida
    d(iHeIM) = + rate(irheii)   * y(iHeIM)           &
               + rate(iche31)   * y(iHI)* y(iHeIM)   &
               + rate(iche31a)  * y(ie )* y(iHeIM)   &
               + rate(iche31b)  * y(ie )* y(iHeIM)   &
               + rate(iphiHeIM) * y(iHeIM )          &
               + rate(icheIM)   * y(ie)*y(iHeIM)

    ! ================ HeII =================
    ! Producción
    q(iHeII) = + rate(iaheii_1)* y(ie )* y(iHeII) * (1.-yy)  &
               + rate(iphiHeIS)* y(iHeIS )                   &
               + rate(icheIS)  * y(iHeIS)*y(ie)              &
               + rate(iphiHeIM) * y(iHeIM )                  &
               + rate(icheIM)   * y(ie)*y(iHeIM)             &
               + rate(iceheII) * y(iHII)*y(iHeIS) !?? Estan bien estos?? Signo correcto en network ??

    ! Pérdida
    d(iHeIM) = + rate(iaheii_a)* y(ie )* y(iHeII)   &
               + rate(iceheIS) * y(iHI)*y(iHeII) !?? Estan bien estos?? Signo correcto en network ??


    ! ================= e ==================
    ! Producción
    q(ie) = + rate(ichi)    * y(iHI)  * y(ie)    &
            + rate(iphiHI)  * y(iHI)             &
            + rate(icheIS)  * y(iHeIS)*y(ie)     &
            + rate(iphiHeIS)* y(iHeIS )          &
            + rate(icheIM)   * y(ie)*y(iHeIM)    &
            + rate(iphiHeIM) * y(iHeIM )

    ! Pérdida
    d(ie) = + rate(iaheii_a)* y(ie )* y(iHeII)   &
            + rate(iahii_b) * y(iHII) * y(ie)


    ! =======================================

end subroutine gsub




 end module network

!=================================
