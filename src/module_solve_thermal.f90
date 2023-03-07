!======================================================================================================================
!> @file        module_solve_thermal.f90
!> @brief       This file contains a module of energy solver for PaScaL_TCS.
!>              MPI communication.
!> @author      
!>              - Kiha Kim (k-kiha@yonsei.ac.kr), Department of Computational Science & Engineering, Yonsei University
!>              - Ji-Hoon Kang (jhkang@kisti.re.kr), Korea Institute of Science and Technology Information
!>              - Jung-Il Choi (jic@yonsei.ac.kr), Department of Computational Science & Engineering, Yonsei University
!>
!> @date        October 2022
!> @version     1.0
!> @par         Copyright
!>              Copyright (c) 2022 Kiha Kim and Jung-Il choi, Yonsei University and 
!>              Ji-Hoon Kang, Korea Institute of Science and Technology Information, All rights reserved.
!> @par         License     
!>              This project is release under the terms of the MIT License (see LICENSE in )
!======================================================================================================================

!>
!> @brief       Module for energy solver.
!> @details     This module solves energy equation and returns temperature. It uses PaScaL_TDMA library to solve 
!>              the decoupled tri-diagonal systems of equations.
!>
module mpi_thermal

    implicit none

    !> @{ Main unknown array of temperature
    double precision, allocatable, dimension(:,:,:) :: T
    !> @}
    !> @{ Boundary condition of temperature at upper and lower walls in y-direction
    double precision, allocatable, dimension(:,:)   :: TBCbt_sub,TBCup_sub
    !> @}
    !> @{ Local arrays for thermal coefficients
    double precision, allocatable, dimension(:,:,:) :: KPP,dKPP,invRhoCp,dinvRhoCp
    !> @}

    contains

    !>
    !> @brief       Assign variables for energy solver
    !>
    subroutine mpi_thermal_allocation()

        use mpi_subdomain,  only : n1sub, n2sub, n3sub

        implicit none

        allocate(T(0:n1sub,0:n2sub,0:n3sub))
        allocate(TBCup_sub(0:n1sub,0:n3sub),TBCbt_sub(0:n1sub,0:n3sub))

        T(:,:,:) = 0.d0
        TBCup_sub(:,:) = 0.d0
        TBCbt_sub(:,:) = 0.d0

    end subroutine mpi_thermal_allocation

    !>
    !> @brief       Deallocate variables for energy solver
    !>
    subroutine mpi_thermal_clean()

        implicit none
        
        deallocate(T)
        deallocate(TBCup_sub,TBCbt_sub)
        
    end subroutine mpi_thermal_clean

    !>
    !> @brief       Initialize temperature filed
    !>
    subroutine mpi_thermal_initial()

        use mpi_subdomain,  only : n1sub, n2sub, n3sub

        implicit none

        integer :: i,j,k
        
        T(0:n1sub,0:n2sub,0:n3sub)=0.d0

    end subroutine mpi_thermal_initial

    !>
    !> @brief       Initialize temperature boundary condition
    !>
    subroutine mpi_thermal_boundary()

        use global,         only : Tup, Tbt
        use mpi_topology,   only : comm_1d_x1, comm_1d_x2, comm_1d_x3
        use mpi_subdomain,  only : n1sub, n2sub, n3sub

        implicit none

        integer :: k
        
       ! Example problem : Constant temperature for upper and lower walls in y-direction
        TBCup_sub(:,:)=Tup
        TBCbt_sub(:,:)=Tbt
        
        if(comm_1d_x2%myrank==0) then
            do k=0, n3sub
                T(0:n1sub,0,k)=TBCbt_sub(:,k)
            enddo
        endif
        if(comm_1d_x2%myrank==comm_1d_x2%nprocs-1) then
            do k=0, n3sub
                T(0:n1sub,n2sub,k)=TBCup_sub(:,k)
            enddo
        endif    
    end subroutine mpi_thermal_boundary

    !>
    !> @brief       Assign and calculate the thermal coefficients
    !>
    subroutine mpi_thermal_coeffi()

        use global
        use mpi_subdomain,  only : n1sub, n2sub, n3sub

        implicit none

        integer :: i,j,k
        double precision :: tmp1, tmp2, tmp3, tmp4, tmp5

        allocate(KPP(0:n1sub,0:n2sub,0:n3sub),dKPP(0:n1sub,0:n2sub,0:n3sub))
        allocate(invRhoCp(0:n1sub,0:n2sub,0:n3sub),dinvRhoCp(0:n1sub,0:n2sub,0:n3sub))
        
        do k = 0, n3sub
            do j = 0, n2sub
                do i = 0, n1sub
                tmp1 = DeltaT*T(i,j,k)
                tmp2 = tmp1 * tmp1
                tmp3 = tmp2 * tmp1
                tmp4 = tmp3 * tmp1
                tmp5 = tmp4 * tmp1
                KPP(i,j,k)  = c10*(1.d0 + c11*tmp1 +      c12*tmp2 +      c13*tmp3 +      c14*tmp4 +      c15*tmp5 )/KPP0            
                dKPP(i,j,k) = c10*(       c11      + 2.d0*c12*tmp1 + 3.d0*c13*tmp2 + 4.d0*c14*tmp3 + 5.d0*c15*tmp4 )/KPP0
                    
                invRhoCp(i,j,k) = Rho0*Cp0  &
                                 /(a10 + a10*(a11*tmp1 + a12*tmp2 + a13*tmp3 + a14*tmp4 + a15*tmp5)) &
                                 /(b10 + b10*(b11*tmp1 + b12*tmp2 + b13*tmp3 + b14*tmp4 + b15*tmp5))
                
                dinvRhoCp(i,j,k) = -( b11*(1.d0 + 3.d0*a12*tmp2 + 4.d0*a13*tmp3)                  &
                                     +a11*(1.d0 + 2.d0*b11*tmp1 + 3.d0*b12*tmp2 + 4.d0*b13*tmp3)  &
                                     +tmp1*( b12*(2.d0 + 5.d0*a13*tmp3)                           &
                                            +a12*(2.d0 + 4.d0*b12*tmp2 + 5.d0*b13*tmp3)           &
                                            +3.d0*tmp1*(a13 + b13 +  2.d0*a13*b13*tmp3)    )      &
                                    )/(1.d0 + a11*tmp1 + a12*tmp2 + a13*tmp3)**2.d0               &
                                     /(1.d0 + b11*tmp1 + b12*tmp2 + b13*tmp3)**2.d0
                
                enddo   
            enddo
        enddo

    end subroutine mpi_thermal_coeffi

    !>
    !> @brief       deallocate the thermal coefficients
    !>
    subroutine mpi_thermal_coeffi_clean()

        implicit none
        
        deallocate(KPP,dKPP,invRhoCp,dinvRhoCp)
        
    end subroutine mpi_thermal_coeffi_clean
    
    !>
    !> @brief       Main energy solver for temperature
    !> @param       U       Velocity field in x-direction
    !> @param       V       Velocity field in y-direction
    !> @param       W       Velocity field in z-direction
    !>
    subroutine mpi_thermal_solver(U,V,W)

        use global,         only : Ct, dt
        use mpi_topology,   only : comm_1d_x1, comm_1d_x2, comm_1d_x3
        use mpi_subdomain,  only : n1sub, n2sub, n3sub, n1msub, n2msub, n3msub
        use mpi_subdomain,  only : iC_BC, jC_BC, kC_BC
        use mpi_subdomain,  only : dx1_sub,dx2_sub,dx3_sub,dmx1_sub,dmx2_sub,dmx3_sub

        use PaScaL_TDMA

        implicit none
        
        double precision, dimension(0:n1sub,0:n2sub,0:n3sub) ::  U,V,W

        type(ptdma_plan_many)     :: ptdma_plan
        
        !> @{ Local pointer for subdomain variables
        double precision, dimension(:), pointer :: dx1, dx2, dx3
        double precision, dimension(:), pointer :: dmx1, dmx2, dmx3
        !> @}

        !> @{ Local indexing variables
        integer :: im, i,  ip
        integer :: jm, j,  jp
        integer :: km, k,  kp
        integer :: iem,iec,iep
        integer :: jem,jec,jep
        integer :: kem,kec,kep
        !> @}

        !> @{ Local arrays for matrix and RHS coefficients 
        double precision, allocatable, dimension(:,:,:) :: dTc,RHS,RHSI,RHSJ
        double precision, allocatable, dimension(:,:,:) :: AMK, ACK, APK, AMI, ACI, API, AMJ, ACJ, APJ
        !> @}

        !> @{ Local variables for matrix and RHS coefficients calculation
        double precision :: e1,e2,e3,e4,e5,e6
        double precision :: dedx1,dedx2,dedy3,dedy4,dedz5,dedz6
        double precision :: convect_e1,convect_e2,convect_e3
        double precision :: KPP1,KPP2,KPP3,KPP4,KPP5,KPP6
        double precision :: dKPP1,dKPP2,dKPP3,dKPP4,dKPP5,dKPP6
        double precision :: viscous_e1,viscous_e2,viscous_e3
        double precision :: ebc_up,ebc_down
        double precision :: eACK,eAPK,eAMK,eACI,eAPI,eAMI,eACJ,eAPJ,eAMJ    ! 1D-32
        double precision :: RHS_ijk
    
        double precision :: Tijk, Tip, Tim, Tjp, Tjm, Tkp, Tkm
        double precision :: Uijk, Uip, Uim, Ujp, Ujm, Ukp, Ukm
        double precision :: Vijk, Vip, Vim, Vjp, Vjm, Vkp, Vkm
        double precision :: Wijk, Wip, Wim, Wjp, Wjm, Wkp, Wkm
        double precision ::  KPPijk,  KPPip,  KPPim,  KPPjp,  KPPjm,  KPPkp,  KPPkm
        double precision :: dKPPijk, dKPPip, dKPPim, dKPPjp, dKPPjm, dKPPkp, dKPPkm
        double precision ::  invRhoCpijk
        double precision :: dinvRhoCpijk
        !> @}

        ! Pointer of grid information
        dx1 => dx1_sub
        dx2 => dx2_sub
        dx3 => dx3_sub
        dmx1 => dmx1_sub
        dmx2 => dmx2_sub
        dmx3 => dmx3_sub
    
        ! Allocate arrays for matrix and RHS coefficients
        allocate(RHS(1:n1msub,1:n2msub,1:n3msub))        
        allocate(AMK(1:n1msub,1:n2msub,1:n3msub), ACK(1:n1msub,1:n2msub,1:n3msub), APK(1:n1msub,1:n2msub,1:n3msub))
        allocate(AMI(1:n2msub,1:n3msub,1:n1msub), ACI(1:n2msub,1:n3msub,1:n1msub), API(1:n2msub,1:n3msub,1:n1msub))
        allocate(AMJ(1:n3msub,1:n1msub,1:n2msub), ACJ(1:n3msub,1:n1msub,1:n2msub), APJ(1:n3msub,1:n1msub,1:n2msub))

        ! 1st stage : build a matrix and RHS term to solve the tridagonal systems in z-direction
        do k = 1, n3msub
            km = k-1
            kp = k+1
            kem = kC_BC(k-1); kec = kC_BC(k); kep = kC_BC(k+1) 
            do j = 1, n2msub
                jm = j-1
                jp = j+1
                jem = jC_BC(j-1); jec = jC_BC(j); jep = jC_BC(j+1)
                do i = 1, n1msub
                    im = i-1
                    ip = i+1
                    iem = iC_BC(i-1); iec = iC_BC(i); iep = iC_BC(i+1)

                    Tijk = T(i,j,k); Tip  = T(ip,j,k); Tim  = T(im,j,k); Tjp  = T(i,jp,k); Tjm  = T(i,jm,k); Tkp  = T(i,j,kp); Tkm  = T(i,j,km);
                    Uijk = U(i,j,k); Uip  = U(ip,j,k); Uim  = U(im,j,k); Ujp  = U(i,jp,k); Ujm  = U(i,jm,k); Ukp  = U(i,j,kp); Ukm  = U(i,j,km);
                    Vijk = V(i,j,k); Vip  = V(ip,j,k); Vim  = V(im,j,k); Vjp  = V(i,jp,k); Vjm  = V(i,jm,k); Vkp  = V(i,j,kp); Vkm  = V(i,j,km);
                    Wijk = W(i,j,k); Wip  = W(ip,j,k); Wim  = W(im,j,k); Wjp  = W(i,jp,k); Wjm  = W(i,jm,k); Wkp  = W(i,j,kp); Wkm  = W(i,j,km);
                     KPPijk =  KPP(i,j,k);  KPPip =  KPP(ip,j,k);  KPPim =  KPP(im,j,k);  KPPjp =  KPP(i,jp,k);  KPPjm =  KPP(i,jm,k);  KPPkp  =  KPP(i,j,kp);  KPPkm  =  KPP(i,j,km);
                    dKppijk = dKpp(i,j,k); dKppip = dKpp(ip,j,k); dKppim = dKpp(im,j,k); dKppjp = dKpp(i,jp,k); dKppjm = dKpp(i,jm,k); dKppkp  = dKpp(i,j,kp); dKppkm  = dKpp(i,j,km);
                     invRhoCpijk =  invRhoCp(i,j,k);
                    dinvRhoCpijk = dinvRhoCp(i,j,k); 
                
                    ! Convection term
                    e1 = (dx1(i )*Tim  + dx1(im)*Tijk)*(0.5d0/dmx1(i ))
                    e2 = (dx1(ip)*Tijk + dx1(i )*Tip )*(0.5d0/dmx1(ip))
                    e3 = (dx2(j )*Tjm  + dx2(jm)*Tijk)*(0.5d0/dmx2(j ))
                    e4 = (dx2(jp)*Tijk + dx2(j )*Tjp )*(0.5d0/dmx2(jp))  
                    e5 = (dx3(k )*Tkm  + dx3(km)*Tijk)*(0.5d0/dmx3(k ))
                    e6 = (dx3(k )*Tkp  + dx3(kp)*Tijk)*(0.5d0/dmx3(kp))
                
                    dedx1 = (Tijk - Tim )/dmx1(i )
                    dedx2 = (Tip  - Tijk)/dmx1(ip)                      
                    dedy3 = (Tijk - Tjm )/dmx2(j )
                    dedy4 = (Tjp  - Tijk)/dmx2(jp)
                    dedz5 = (Tijk - Tkm )/dmx3(k )
                    dedz6 = (Tkp  - Tijk)/dmx3(kp) 
                
                    ! Diffusion term
                    KPP1 = 0.5d0/dmx1(i )*(dx1(i )*KPPim  + dx1(im)*KPPijk)
                    KPP2 = 0.5d0/dmx1(ip)*(dx1(ip)*KPPijk + dx1(i )*KPPip ) 
                    KPP3 = 0.5d0/dmx2(j )*(dx2(j )*KPPjm  + dx2(jm)*KPPijk)
                    KPP4 = 0.5d0/dmx2(jp)*(dx2(jp)*KPPijk + dx2(j )*KPPjp ) 
                    KPP5 = 0.5d0/dmx3(k )*(dx3(k )*KPPkm  + dx3(km)*KPPijk)
                    KPP6 = 0.5d0/dmx3(kp)*(dx3(kp)*KPPijk + dx3(k )*KPPkp )
                    
                    dKPP1 = 0.5d0/dmx1(i )*(dx1(i )*dKPPim  + dx1(im)*dKPPijk)
                    dKPP2 = 0.5d0/dmx1(ip)*(dx1(ip)*dKPPijk + dx1(i )*dKPPip )
                    dKPP3 = 0.5d0/dmx2(j )*(dx2(j )*dKPPjm  + dx2(jm)*dKPPijk)
                    dKPP4 = 0.5d0/dmx2(jp)*(dx2(jp)*dKPPijk + dx2(j )*dKPPjp )
                    dKPP5 = 0.5d0/dmx3(k )*(dx3(k )*dKPPkm  + dx3(km)*dKPPijk)
                    dKPP6 = 0.5d0/dmx3(kp)*(dx3(kp)*dKPPijk + dx3(k )*dKPPkp )
                
                    !MTn--------------------------------------------------------------------------------------------
                    !1   : X-direction
                    !1-1 : Viscous
                    eAMI = - 0.5d0* Ct*invRhoCpijk /dx1(i)*(                   KPP1/dmx1(i)                                          - dx1(i )*(0.5d0/dmx1(i))*dKPP1*dedx1 )
                    eACI = - 0.5d0* Ct*invRhoCpijk /dx1(i)*( - KPP2/dmx1(ip) - KPP1/dmx1(i) + dx1(ip)*(0.5d0/dmx1(ip)) * dKPP2*dedx2 - dx1(im)*(0.5d0/dmx1(i))*dKPP1*dedx1 )
                    eAPI = - 0.5d0* Ct*invRhoCpijk /dx1(i)*(   KPP2/dmx1(ip)                + dx1(i )*(0.5d0/dmx1(ip)) * dKPP2*dedx2                                       )
                    !1-2 : Convection  
                    eAMI = eAMI + (                         - 0.25d0*Uijk/dmx1(i) )
                    eACI = eACI + ( - 0.25d0*Uip/dmx1(ip) + 0.25d0*Uijk/dmx1(i) )
                    eAPI = eAPI + ( + 0.25d0*Uip/dmx1(ip)                         )
                    
                    !2   : Y-direction
                    !2-1 : Viscous 
                    eAMJ = - 0.5d0* Ct*invRhoCpijk /dx2(j)*(                   KPP3/dmx2(j)                                                     - dx2(j )*(0.5d0/dmx2(j))*dKPP3*dedy3            )
                    eACJ = - 0.5d0* Ct*invRhoCpijk /dx2(j)*( - KPP4/dmx2(jp) - KPP3/dmx2(j) + dx2(jp)*(0.5d0/dmx2(jp)) * dKPP4*dedy4 *dble(jep) - dx2(jm)*(0.5d0/dmx2(j))*dKPP3*dedy3 *dble(jem) )
                    eAPJ = - 0.5d0* Ct*invRhoCpijk /dx2(j)*(   KPP4/dmx2(jp)                + dx2(j )*(0.5d0/dmx2(jp)) * dKPP4*dedy4                                                             )
                    !2-2 : Convection 
                    eAMJ = eAMJ + (                       - 0.25d0*Vijk/dmx2(j)  )
                    eACJ = eACJ + ( - 0.25d0*Vjp/dmx2(jp) + 0.25d0*Vijk/dmx2(j)  )
                    eAPJ = eAPJ + ( + 0.25d0*Vjp/dmx2(jp)                        )
                    !2-3 : BC check
                    eAMJ = eAMJ*dble(jem)
                    eAPJ = eAPJ*dble(jep)    
                
                
                    !3   : Z-direction
                    !3-1 : Viscous
                    eAMK = - 0.5d0* Ct*invRhoCpijk /dx3(k)*(                   KPP5/dmx3(k)                                          - dx3(k )*(0.5d0/dmx3(k))*dKPP5*dedz5 )
                    eACK = - 0.5d0* Ct*invRhoCpijk /dx3(k)*( - KPP6/dmx3(kp) - KPP5/dmx3(k) + dx3(kp)*(0.5d0/dmx3(kp)) * dKPP6*dedz6 - dx3(km)*(0.5d0/dmx3(k))*dKPP5*dedz5 )
                    eAPK = - 0.5d0* Ct*invRhoCpijk /dx3(k)*(   KPP6/dmx3(kp)                + dx3(k )*(0.5d0/dmx3(kp)) * dKPP6*dedz6                                       )
                    !3-2 : Convection
                    eAMK = eAMK + (                         - 0.25d0*Wijk/dmx3(k) )
                    eACK = eACK + ( - 0.25d0*Wkp/dmx3(kp) + 0.25d0*Wijk/dmx3(k) )
                    eAPK = eAPK + ( + 0.25d0*Wkp/dmx3(kp)                         )
                    
                    !4   : RHS (RHSijk = -convect + viscous + ebc -RHS_e)
                    convect_e1 = (Uijk*dedx1 + Uip*dedx2)*0.5d0
                    convect_e2 = (Vijk*dedy3 + Vjp*dedy4)*0.5d0 
                    convect_e3 = (Wijk*dedz5 + Wkp*dedz6)*0.5d0
                    !4-1 : rn (* 1/dt*un are offseted in Au^n term.)
                    RHS_ijk  = - 0.5d0 * (convect_e1 + convect_e2 + convect_e3)
                    
                    viscous_e1 = (KPP2*dedx2 - KPP1*dedx1 - e2*dKPP2*dedx2 + e1*dKPP1*dedx1)/dx1(i)
                    viscous_e2 = (KPP4*dedy4 - KPP3*dedy3 - e4*dKPP4*dedy4 + e3*dKPP3*dedy3)/dx2(j)
                    viscous_e3 = (KPP6*dedz6 - KPP5*dedz5 - e6*dKPP6*dedz6 + e5*dKPP5*dedz5)/dx3(k)
                    !4-1 : rn (* 1/dt*un are offseted in Au^n term.)
                    RHS_ijk  = RHS_ijk + 0.5d0* Ct * invRhoCpijk      * ( viscous_e1 + viscous_e2 + viscous_e3 ) 
                    RHS_ijk  = RHS_ijk - 0.5d0* Ct *dinvRhoCpijk*Tijk * ( (KPP2*dedx2 - KPP1*dedx1)/dx1(i) + (KPP4*dedy4 - KPP3*dedy3)/dx2(j) + (KPP6*dedz6 - KPP5*dedz5)/dx3(k) )  
                    
                    !4-2A:    ! MTn
                    RHS_ijk  = RHS_ijk - ( +eAMI*Tim + eACI*Tijk + eAPI*Tip &
                                           +eAMJ*Tjm + eACJ*Tijk + eAPJ*Tjp &
                                           +eAMK*Tkm + eACJ*Tijk + eAPK*Tkp )
                
                    !4-3 : BC
                    ! From Convection Terms
                    ebc_up   = -0.25d0*Vjp /dmx2(jp)*TBCup_sub(i,k)
                    ebc_down = +0.25d0*Vijk/dmx2(j )*TBCbt_sub(i,k)
                    
                    ! From Diffusion Terms
                    ebc_up   = ebc_up                                                          &
                             + 0.5d0*Ct*invRhoCpijk/dx2(j)* KPP4/dmx2(jp)*TBCup_sub(i,k) & 
                             + 0.5d0*Ct*invRhoCpijk/dx2(j)*dKPP4*dedy4   *TBCup_sub(i,k)   
                    ebc_down = ebc_down                                                        &
                             + 0.5d0*Ct*invRhoCpijk/dx2(j)* KPP3/dmx2(j )*TBCbt_sub(i,k) &
                             - 0.5d0*Ct*invRhoCpijk/dx2(j)*dKPP3*dedy3   *TBCbt_sub(i,k)
                    
                    ! ebc For Y-direction
                    ! RHS_dijk = -convect + viscous + ebc -RHS_e
                    RHS_ijk  = RHS_ijk + ( dble(1-jem) ) *ebc_down   &
                                       + ( dble(1-jep) ) *ebc_up
                    
                
                    RHS_ijk  = RHS_ijk + 0.5d0* Ct*dinvRhoCpijk*Tijk *( +(KPP2*dedx2 - KPP1*dedx1)/dx1(i) + (KPP4*dedy4 - KPP3*dedy3)/dx2(j) + (KPP6*dedz6 - KPP5*dedz5)/dx3(k) )
                    
                    ! Note that eAC in RHS and AC are different.
                    eACI = eACI - 0.25d0*Ct*dinvRhoCpijk*( +(KPP2*dedx2 - KPP1*dedx1)/dx1(i) +(KPP4*dedy4 - KPP3*dedy3)/dx2(j) +(KPP6*dedz6 - KPP5*dedz5)/dx3(k))
                    eACJ = eACJ - 0.25d0*Ct*dinvRhoCpijk*( +(KPP2*dedx2 - KPP1*dedx1)/dx1(i) +(KPP4*dedy4 - KPP3*dedy3)/dx2(j) +(KPP6*dedz6 - KPP5*dedz5)/dx3(k))
                    eACK = eACK - 0.25d0*Ct*dinvRhoCpijk*( +(KPP2*dedx2 - KPP1*dedx1)/dx1(i) +(KPP4*dedy4 - KPP3*dedy3)/dx2(j) +(KPP6*dedz6 - KPP5*dedz5)/dx3(k))
                
                    ! MTn: Z-direction
                    AMK(i,j,k) = eAMK     * dt
                    ACK(i,j,k) = eACK     * dt + 1.0d0
                    APK(i,j,k) = eAPK     * dt
                    RHS(i,j,k) = RHS_ijk  * dt
                
                    ! MTn: X-direction
                    AMI(j,k,i) = eAMI * dt
                    ACI(j,k,i) = eACI * dt + 1.0d0
                    API(j,k,i) = eAPI * dt
                
                    ! MTn: Y-direction
                    AMJ(k,i,j) = eAMJ * dt
                    ACJ(k,i,j) = eACJ * dt + 1.0d0
                    APJ(k,i,j) = eAPJ * dt
                enddo
            enddo
        enddo

        ! 1st stage : solve TDM in z-direction
        call PaScaL_TDMA_plan_many_create(ptdma_plan, n1msub*n2msub, comm_1d_x3%myrank, comm_1d_x3%nprocs, comm_1d_x3%mpi_comm)
        call PaScaL_TDMA_many_solve_cycle(ptdma_plan, AMK, ACK, APK, RHS, n1msub*n2msub,n3msub)
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x3%nprocs)

        ! Deallocate A-matrix for the 2nd stage
        deallocate(AMK, ACK, APK)

        ! Transpose the array (i,j,k) to array (j,k,i) to solve the tridagonal systems in x-direction
        allocate(RHSI(1:n2msub,1:n3msub,1:n1msub))
        do i = 1, n1msub
        do k = 1, n3msub
        do j = 1, n2msub
            RHSI(j,k,i)=RHS(i,j,k)
        enddo
        enddo
        enddo

        ! Deallocate the original RHS array which is not required any more
        deallocate(RHS) 

        ! 2nd stage : solve TDM in x-direction
        call PaScaL_TDMA_plan_many_create(ptdma_plan, n2msub*n3msub, comm_1d_x1%myrank, comm_1d_x1%nprocs, comm_1d_x1%mpi_comm)
        call PaScaL_TDMA_many_solve_cycle(ptdma_plan, AMI, ACI, API, RHSI, n2msub*n3msub,n1msub)
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x1%nprocs)

        ! Deallocate A-matrix for the 3rd stage
        deallocate(AMI, ACI, API)

        ! Transpose the array (j,k,i) to array (i,k,j) to solve the tridagonal systems in y-direction
        allocate(RHSJ(1:n3msub,1:n1msub,1:n2msub))
        do j = 1, n2msub
        do i = 1, n1msub
        do k = 1, n3msub
            RHSJ(k,i,j)=RHSI(j,k,i)
        enddo
        enddo
        enddo

        ! Deallocate the RHSI array which is not required any more
        deallocate(RHSI)

        ! 3rd stage : solve TDM in y-direction
        call PaScaL_TDMA_plan_many_create(ptdma_plan, n3msub*n1msub, comm_1d_x2%myrank, comm_1d_x2%nprocs, comm_1d_x2%mpi_comm)
        call PaScaL_TDMA_many_solve(ptdma_plan, AMJ, ACJ, APJ, RHSJ, n3msub*n1msub, n2msub)
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x2%nprocs)

        ! Deallocate A-matrix for the 3rd stage
        deallocate(AMJ, ACJ, APJ)

        ! Update the temperature with the transpose of the solution array (i,k,j)
        do k = 1, n3msub
        do j = 1, n2msub
        do i = 1, n1msub
            T(i,j,k) = T(i,j,k) + RHSJ(k,i,j)
        enddo
        enddo
        enddo

        ! Deallocate the RHSJ array which is not required any more
        deallocate(RHSJ)

        ! Nullify grid information pointer
        nullify(dx1, dx2, dx3, dmx1, dmx2, dmx3)

    end subroutine mpi_thermal_solver

end module mpi_thermal
