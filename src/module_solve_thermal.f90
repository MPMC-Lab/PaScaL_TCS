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
        else if(comm_1d_x2%myrank==comm_1d_x2%nprocs-1) then
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
        use mpi_subdomain,  only : n1sub, n2sub, n3sub, n1msub, n2msub, n3msub, jC_BC
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
        integer :: i,j,k
        integer :: im,jm,km
        integer :: ip,jp,kp
        integer :: jep,jem
        !> @}

        !> @{ Local arrays for matrix and RHS coefficients 
        double precision, allocatable, dimension(:,:,:) :: dTc,RHS,RHSI,RHSJ
        double precision, allocatable, dimension(:,:,:) :: am,ac,ap
        !> @}

        !> @{ Local variables for matrix and RHS coefficients calculation
        double precision :: e1,e2,e3,e4,e5,e6
        double precision :: dedx1,dedx2,dedy3,dedy4,dedz5,dedz6
        double precision :: convect_e1,convect_e2,convect_e3
        double precision :: viscous_e1,viscous_e2,viscous_e3
        double precision :: KPP1,KPP2,KPP3,KPP4,KPP5,KPP6
        double precision :: dKPP1,dKPP2,dKPP3,dKPP4,dKPP5,dKPP6
        double precision :: ebc_up,ebc_down
        double precision :: eACK,eAPK,eAMK
        double precision :: eACI,eAPI,eAMI
        double precision :: eACJ,eAPJ,eAMJ
        double precision :: invRhoCpCt_half, dinvRhoCpCt_half
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
        allocate(am(1:n1msub,1:n2msub,1:n3msub))
        allocate(ac(1:n1msub,1:n2msub,1:n3msub))
        allocate(ap(1:n1msub,1:n2msub,1:n3msub))

        ! 1st stage : build a matrix and RHS term to solve the tridagonal systems in z-direction
        do k = 1, n3msub

            km=k-1
            kp=k+1

            do j = 1, n2msub

                jm=j-1
                jp=j+1
                jem = jC_BC(jm)
                jep = jC_BC(jp)

                do i = 1, n1msub

                    ip = i + 1
                    im = i - 1

                    invRhoCpCt_half  = 0.5d0 * Ct * invRhoCp(i,j,k)
                    dinvRhoCpCt_half = 0.5d0 * Ct * dinvRhoCp(i,j,k)

                    ! Convection term
                    e1 = (dx1(i )*T(im,j, k ) + dx1(im)*T(i,j,k )) *(0.5d0/dmx1(i))
                    e2 = (dx1(ip)*T(i, j, k ) + dx1(i )*T(ip,j,k)) *(0.5d0/dmx1(ip))
                    e3 = (dx2(j )*T(i, jm,k ) + dx2(jm)*T(i,j,k )) *(0.5d0/dmx2(j))
                    e4 = (dx2(jp)*T(i, j, k ) + dx2(j )*T(i,jp,k)) *(0.5d0/dmx2(jp))  
                    e5 = (dx3(k )*T(i, j, km) + dx3(km)*T(i,j,k )) *(0.5d0/dmx3(k))
                    e6 = (dx3(kp)*T(i, j, k ) + dx3(k )*T(i,j,kp)) *(0.5d0/dmx3(kp))
                    dedx1 = (T(i, j, k) - T(im,j, k)) /dmx1(i)
                    dedx2 = (T(ip,j, k) - T(i, j, k)) /dmx1(ip)                      
                    dedy3 = (T(i, j, k) - T(i, jm,k)) /dmx2(j)
                    dedy4 = (T(i, jp,k) - T(i, j, k)) /dmx2(jp)
                    dedz5 = (T(i, j, k )- T(i, j, km))/dmx3(k)
                    dedz6 = (T(i, j, kp)- T(i, j, k)) /dmx3(kp)

                    ! 1/2 ( U(Tp)x + V(Tp)y + W(Tp)z )
                    convect_e1 = (U(i,j,k)*dedx1 + U(ip,j,k)*dedx2)*0.5d0
                    convect_e2 = (V(i,j,k)*dedy3 + V(i,jp,k)*dedy4)*0.5d0 
                    convect_e3 = (W(i,j,k)*dedz5 + W(i,j,kp)*dedz6)*0.5d0
                    RHS(i,j,k) = -0.5d0*(convect_e1 + convect_e2 + convect_e3)
                    
                    ! Diffusion term
                    kPP1 = 0.5d0/dmx1(i) *(dx1(i) *KPP(im,j,k)+ dx1(im)*KPP(i,j ,k ))
                    kPP2 = 0.5d0/dmx1(ip)*(dx1(ip)*KPP(i,j,k) + dx1(i)*KPP(ip,j ,k )) 
                    kPP3 = 0.5d0/dmx2(j) *(dx2(j) *KPP(i,jm,k)+ dx2(jm)*KPP(i,j ,k ))
                    kPP4 = 0.5d0/dmx2(jp)*(dx2(jp)*KPP(i,j,k) + dx2(j )*KPP(i,jp,k )) 
                    kPP5 = 0.5d0/dmx3(k) *(dx3(k) *KPP(i,j,km)+ dx3(km)*KPP(i,j ,k ))
                    kPP6 = 0.5d0/dmx3(kp)*(dx3(kp)*KPP(i,j,k) + dx3(k )*KPP(i,j ,kp))

                    dkPP1 = 0.5d0/dmx1(i)*(dx1(i)*dKPP(im,j ,k ) + dx1(im)*dKPP(i,j ,k ))
                    dkPP2 = 0.5d0/dmx1(ip)*(dx1(ip)*dKPP(i,j ,k ) + dx1(i)*dKPP(ip,j ,k )) 
                    dkPP3 = 0.5d0/dmx2(j )*(dx2(j )*dKPP(i,jm,k ) + dx2(jm)*dKPP(i,j ,k ))
                    dkPP4 = 0.5d0/dmx2(jp)*(dx2(jp)*dKPP(i,j ,k ) + dx2(j )*dKPP(i,jp,k )) 
                    dkPP5 = 0.5d0/dmx3(k )*(dx3(k )*dKPP(i,j ,km) + dx3(km)*dKPP(i,j ,k ))
                    dkPP6 = 0.5d0/dmx3(kp)*(dx3(kp)*dKPP(i,j ,k ) + dx3(k )*dKPP(i,j ,kp))

                    viscous_e1 = (kPP2*dedx2 - kPP1*dedx1 - e2*dkPP2*dedx2 + e1*dkPP1*dedx1)/dx1(i)
                    viscous_e2 = (kPP4*dedy4 - kPP3*dedy3 - e4*dkPP4*dedy4 + e3*dkPP3*dedy3)/dx2(j)
                    viscous_e3 = (kPP6*dedz6 - kPP5*dedz5 - e6*dkPP6*dedz6 + e5*dkPP5*dedz5)/dx3(k)
                    
                    RHS(i,j,k) = RHS(i,j,k) &
                               + invRhoCpCt_half*(viscous_e1 + viscous_e2 + viscous_e3) &
                               - dinvRhoCpCt_half*T(i,j,k)*( (kPP2*dedx2 - kPP1*dedx1)/dx1(i) &
                                                            +(kPP4*dedy4 - kPP3*dedy3)/dx2(j) &
                                                            +(kPP6*dedz6 - kPP5*dedz5)/dx3(k) )  

                    ! ebc for y-direction- From convection terms
                    ebc_up   =-0.25d0*V(i,jp,k)/dmx2(jp)*TBCup_sub(i,k)
                    ebc_down = 0.25d0*V(i,j ,k)/dmx2(j )*TBCbt_sub(i,k)
                    ! ebc for y-direction- From diffusion terms
                    ebc_up = ebc_up &
                           + invRhoCpCt_half/dx2(j)* kPP4/dmx2(jp)*TBCup_sub(i,k) & 
                           + invRhoCpCt_half/dx2(j)*dkPP4*dedy4*TBCup_sub(i,k)   
                    ebc_down = ebc_down &
                             + invRhoCpCt_half/dx2(j)*kPP3/dmx2(j)*TBCbt_sub(i,k) &
                            - invRhoCpCt_half/dx2(j)*dkPP3*dedy3*TBCbt_sub(i,k)

                    ! ebc for y-direction
                    RHS(i,j,k) = RHS(i,j,k) + dble(1-jem)*ebc_down + dble(1-jep)*ebc_up
                    
                    ! Z-direction term
                    eACK = invRhoCpCt_half/dx3(k)*( kPP6/dmx3(kp) + kPP5/dmx3(k) &
                                                    -dx3(kp)*(0.5d0/dmx3(kp))*dkPP6*dedz6 &
                                                    +dx3(km)*(0.5d0/dmx3(k ))*dkPP5*dedz5 ) &
                                + (-0.25d0*W(i,j,kp)/dmx3(kp) + 0.25d0*W(i,j,k)/dmx3(k ))
                    eAPK =-invRhoCpCt_half/dx3(k)*(kPP6/dmx3(kp) + dx3(k)*(0.5d0/dmx3(kp))*dkPP6*dedz6) + 0.25d0*W(i,j,kp)/dmx3(kp)
                    eAMK =-invRhoCpCt_half/dx3(k)*(kPP5/dmx3(k ) - dx3(k)*(0.5d0/dmx3(k ))*dkPP5*dedz5) - 0.25d0*W(i,j,k )/dmx3(k )

                    ! X-direction term
                    eACI = invRhoCpCt_half/dx1(i)*( kPP2/dmx1(ip) + kPP1/dmx1(i) &
                                                                        -dx1(ip)*(0.5d0/dmx1(ip))*dkPP2*dedx2 &
                                                                        +dx1(im)*(0.5d0/dmx1(i))*dkPP1*dedx1 ) &
                                + (0.25d0*U(i,j,k)/dmx1(i) - 0.25d0*U(ip,j,k)/dmx1(ip))   
                    eAPI =-invRhoCpCt_half/dx1(i)*(kPP2/dmx1(ip) + dx1(i)*(0.5d0/dmx1(ip))*dkPP2*dedx2) + 0.25d0*U(ip,j,k)/dmx1(ip)  
                    eAMI =-invRhoCpCt_half/dx1(i)*(kPP1/dmx1(i)  - dx1(i)*(0.5d0/dmx1(i))*dkPP1*dedx1)  - 0.25d0*U(i,j,k)/dmx1(i)
                            
                    ! Y-direction term
                    eACJ = invRhoCpCt_half/dx2(j)*( kPP4/dmx2(jp) + kPP3/dmx2(j) &
                                                                    -dx2(jp)*(0.5d0/dmx2(jp))*dkPP4*dedy4*dble(jep) &
                                                                    +dx2(jm)*(0.5d0/dmx2(j ))*dkPP3*dedy3*dble(jem) ) &
                                + (0.25d0*V(i,j,k)/dmx2(j ) - 0.25d0*V(i,jp,k)/dmx2(jp))   
                    eAPJ =-invRhoCpCt_half/dx2(j)*(kPP4/dmx2(jp) + dx2(j)*(0.5d0/dmx2(jp))*dkPP4*dedy4) + 0.25d0*V(i,jp,k)/dmx2(jp)
                    eAMJ =-invRhoCpCt_half/dx2(j)*(kPP3/dmx2(j ) - dx2(j)*(0.5d0/dmx2(j ))*dkPP3*dedy3) - 0.25d0*V(i,j ,k)/dmx2(j )                           
                    eAPJ = eAPJ*dble(jep)    
                    eAMJ = eAMJ*dble(jem)

                    ! RHS
                    RHS(i,j,k) = RHS(i,j,k) &
                                -( eAPK*T(i,j ,kp) + eACK*T(i,j ,k ) + eAMK*T(i,j ,km) &
                                +eAPJ*T(i,jp,k ) + eACJ*T(i,j ,k ) + eAMJ*T(i,jm,k ) &
                                +eAPI*T(ip,j ,k ) + eACI*T(i,j ,k ) + eAMI*T(im,j ,k ) &
                                -dinvRhoCpCt_half*T(i,j ,k )*( (kPP2*dedx2 - kPP1*dedx1)/dx1(i) &
                                                                +(kPP4*dedy4 - kPP3*dedy3)/dx2(j) &
                                                                +(kPP6*dedz6 - kPP5*dedz5)/dx3(k) ) )

                    ! Coefficients of tridagonal A-matrix
                    ac(i,j,k)= eACK*dt &
                                - 0.5d0 * dinvRhoCpCt_half*( (kPP2*dedx2 - kPP1*dedx1)/dx1(i) & 
                                                            +(kPP4*dedy4 - kPP3*dedy3)/dx2(j) &
                                                            +(kPP6*dedz6 - kPP5*dedz5)/dx3(k) )*dt &
                                + 1.d0
                    ap(i,j,k)=eAPK*dt
                    am(i,j,k)=eAMK*dt
                    RHS(i,j,k)=RHS(i,j,k)*dt
                enddo
            enddo
        enddo

        ! 1st stage : solve TDM in z-direction
        call PaScaL_TDMA_plan_many_create(ptdma_plan, n1msub*n2msub, comm_1d_x3%myrank, comm_1d_x3%nprocs, comm_1d_x3%mpi_comm)
        call PaScaL_TDMA_many_solve_cycle(ptdma_plan, am, ac, ap, RHS,n1msub*n2msub,n3msub)
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x3%nprocs)

        ! Deallocate A-matrix for the 2nd stage
        deallocate(am,ac,ap)

        ! Transpose the array (i,j,k) to array (j,k,i) to solve the tridagonal systems in x-direction
        allocate(RHSI(1:n2msub,1:n3msub,1:n1msub))
        do k = 1, n3msub
            do i = 1, n1msub
                do j = 1, n2msub
                    RHSI(j,k,i)=RHS(i,j,k)
                enddo
            enddo
        enddo

        ! Deallocate the original RHS array which is not required any more
        deallocate(RHS) 

        ! Re-allocate arrays for A matrix coefficient in the second stage
        allocate(am(1:n2msub,1:n3msub,1:n1msub),ac(1:n2msub,1:n3msub,1:n1msub),ap(1:n2msub,1:n3msub,1:n1msub))

        ! 2nd stage : build A matrix and RHS to solve the tridagonal systems in x-direction
        do k = 1, n3msub
            do j = 1, n2msub
                do i = 1, n1msub

                    ip = i + 1
                    im = i - 1

                    invRhoCpCt_half  = 0.5d0 * Ct * invRhoCp(i,j,k)

                    kPP1  = 0.5d0/dmx1(i )*(dx1(i )* KPP(im,j ,k ) + dx1(im)* KPP(i, j ,k ))
                    kPP2  = 0.5d0/dmx1(ip)*(dx1(ip)* KPP(i, j ,k ) + dx1(i )* KPP(ip,j ,k )) 
                    KPP3  = 0.5d0/dmx2(j )*(dx2(j )* KPP(i, jm,k ) + dx2(jm)* KPP(i, j ,k ))
                    KPP4  = 0.5d0/dmx2(jp)*(dx2(jp)* KPP(i, j ,k ) + dx2(j )* KPP(i, jp,k )) 
                    KPP5  = 0.5d0/dmx3(k )*(dx3(k )* KPP(i, j ,km) + dx3(km)* KPP(i, j ,k ))
                    KPP6  = 0.5d0/dmx3(kp)*(dx3(kp)* KPP(i, j ,k ) + dx3(k )* KPP(i, j ,kp))
    
                    dkPP1 = 0.5d0/dmx1(i )*(dx1(i )*dKPP(im,j ,k ) + dx1(im)*dKPP(i, j ,k ))
                    dkPP2 = 0.5d0/dmx1(ip)*(dx1(ip)*dKPP(i, j ,k ) + dx1(i )*dKPP(ip,j ,k )) 

                    dedx1 = (T(i, j, k ) - T(im,j, k ))/dmx1(i )
                    dedx2 = (T(ip,j, k ) - T(i, j, k ))/dmx1(ip)                      
                    dedy3 = (T(i, j ,k ) - T(i, jm,k ))/dmx2(j )
                    dedy4 = (T(i, jp,k ) - T(i, j, k ))/dmx2(jp)
                    dedz5 = (T(i, j, k ) - T(i, j, km))/dmx3(k )
                    dedz6 = (T(i, j, kp) - T(i, j, k ))/dmx3(kp)
                    ! X-dirction 
                    eACI = invRhoCpCt_half/dx1(i)*( kPP2/dmx1(ip) + kPP1/dmx1(i) &
                                                   -dx1(ip)*(0.5d0/dmx1(ip))*dkPP2*dedx2 &
                                                   +dx1(im)*(0.5d0/dmx1(i ))*dkPP1*dedx1 ) &
                            + (0.25d0*U(i,j,k)/dmx1(i) - 0.25d0*U(ip,j,k)/dmx1(ip)) &  
                            - 0.5d0*invRhoCpCt_half*((kPP2*dedx2 - kPP1*dedx1)/dx1(i) & 
                                                    +(kPP4*dedy4 - kPP3*dedy3)/dx2(j) &
                                                    +(kPP6*dedz6 - kPP5*dedz5)/dx3(k) )
                    eAPI =-invRhoCpCt_half/dx1(i)*(kPP2/dmx1(ip) + dx1(i)*(0.5d0/dmx1(ip))*dkPP2*dedx2) + 0.25d0*U(ip,j,k)/dmx1(ip)
                    eAMI =-invRhoCpCt_half/dx1(i)*(kPP1/dmx1(i) - dx1(i)*(0.5d0/dmx1(i))*dkPP1*dedx1)   - 0.25d0*U(i,j,k)/dmx1(i)

                    ac(j,k,i)=eACI*dt + 1.d0
                    ap(j,k,i)=eAPI*dt
                    am(j,k,i)=eAMI*dt
                enddo
            enddo
        enddo

        ! 2nd stage : solve TDM in x-direction
        call PaScaL_TDMA_plan_many_create(ptdma_plan, n2msub*n3msub, comm_1d_x1%myrank, comm_1d_x1%nprocs, comm_1d_x1%mpi_comm)
        call PaScaL_TDMA_many_solve_cycle(ptdma_plan, am, ac, ap, RHSI,n2msub*n3msub,n1msub)
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x1%nprocs)

        ! Deallocate A-matrix for the 3rd stage
        deallocate(am,ac,ap)

        ! Transpose the array (j,k,i) to array (i,k,j) to solve the tridagonal systems in y-direction
        allocate(RHSJ(1:n1msub,1:n3msub,1:n2msub))
        do j = 1, n2msub
            do k = 1, n3msub
                do i = 1, n1msub
                    RHSJ(i,k,j)=RHSI(j,k,i)
                enddo
            enddo
        enddo

        ! Deallocate the RHSI array which is not required any more
        deallocate(RHSI)

        ! Re-allocate arrays for A matrix coefficient in the third stage
        allocate(am(1:n1msub,1:n3msub,1:n2msub),ac(1:n1msub,1:n3msub,1:n2msub),ap(1:n1msub,1:n3msub,1:n2msub))

        ! 3rd stage : build A matrix and RHS to solve the tridagonal systems in y-direction  
        do k = 1, n3msub
            do j = 1, n2msub

                jm=j-1
                jp=j+1
                jem = jC_BC(jm)
                jep = jC_BC(jp)

                do i = 1, n1msub

                    ip = i + 1
                    im = i - 1

                    invRhoCpCt_half  = 0.5d0 * Ct * invRhoCp(i,j,k)
                
                    KPP1  = 0.5d0/dmx1(i )*(dx1(i )* KPP(im,j ,k ) + dx1(im)* KPP(i, j ,k ))
                    KPP2  = 0.5d0/dmx1(ip)*(dx1(ip)* KPP(i, j ,k ) + dx1(i )* KPP(ip,j ,k )) 
                    kPP3  = 0.5d0/dmx2(j )*(dx2(j )* KPP(i, jm,k ) + dx2(jm)* KPP(i, j ,k ))
                    kPP4  = 0.5d0/dmx2(jp)*(dx2(jp)* KPP(i, j ,k ) + dx2(j )* KPP(i, jp,k ))
                    KPP5  = 0.5d0/dmx3(k )*(dx3(k )* KPP(i, j ,km) + dx3(km)* KPP(i, j ,k ))
                    KPP6  = 0.5d0/dmx3(kp)*(dx3(kp)* KPP(i, j ,k ) + dx3(k )* KPP(i, j ,kp))
    
                    dkPP3 = 0.5d0/dmx2(j )*(dx2(j )*dKPP(i, jm,k ) + dx2(jm)*dKPP(i, j ,k ))
                    dkPP4 = 0.5d0/dmx2(jp)*(dx2(jp)*dKPP(i, j ,k ) + dx2(j )*dKPP(i, jp,k )) 

                    dedx1 = (T(i,j,k)  - T(im,j,k))/dmx1(i)
                    dedx2 = (T(ip,j,k) - T(i,j,k))/dmx1(ip)          
                    dedy3 = (T(i,j ,k) - T(i,jm,k))/dmx2(j )
                    dedy4 = (T(i,jp,k) - T(i,j ,k))/dmx2(jp)
                    dedz5 = (T(i,j,k ) - T(i,j,km))/dmx3(k )
                    dedz6 = (T(i,j,kp) - T(i,j,k ))/dmx3(kp)
                    
                    ! Y-DIRECTION   
                    eACJ = invRhoCpCt_half/dx2(j)*( kPP4/dmx2(jp) + kPP3/dmx2(j) &
                                                   -dx2(jp)*(0.5d0/dmx2(jp))*dkPP4*dedy4*dble(jep) &
                                                   +dx2(jm)*(0.5d0/dmx2(j ))*dkPP3*dedy3*dble(jem) ) &
                            + (0.25d0*V(i,j,k)/dmx2(j ) - 0.25d0*V(i,jp,k)/dmx2(jp))  &
                            - 0.5d0*invRhoCpCt_half*((KPP2*dedx2 - KPP1*dedx1)/dx1(i) & 
                                                    +(KPP4*dedy4 - KPP3*dedy3)/dx2(j) &
                                                    +(KPP6*dedz6 - KPP5*dedz5)/dx3(k)   ) 

                    eAPJ =-invRhoCpCt_half/dx2(j)*(kPP4/dmx2(jp) + dx2(j)*(0.5d0/dmx2(jp))*dkPP4*dedy4) + 0.25d0*V(i,jp,k)/dmx2(jp)
                    eAMJ =-invRhoCpCt_half/dx2(j)*(kPP3/dmx2(j ) - dx2(j)*(0.5d0/dmx2(j ))*dkPP3*dedy3) - 0.25d0*V(i,j ,k)/dmx2(j )                           
                    eAPJ = eAPJ*dble(jep)    
                    eAMJ = eAMJ*dble(jem)

                    ac(i,k,j)=eACJ*dt + 1.d0
                    ap(i,k,j)=eAPJ*dt
                    am(i,k,j)=eAMJ*dt
                enddo
            enddo
        enddo
        
        ! 3rd stage : solve TDM in y-direction
        call PaScaL_TDMA_plan_many_create(ptdma_plan, n3msub*n1msub, comm_1d_x2%myrank, comm_1d_x2%nprocs, comm_1d_x2%mpi_comm)
        call PaScaL_TDMA_many_solve(ptdma_plan, am, ac, ap, RHSJ, n3msub*n1msub, n2msub)
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x2%nprocs)

        ! Deallocate A-matrix for the 3rd stage
        deallocate(am,ac,ap)

        ! Update the temperature with the transpose of the solution array (i,k,j)
        do k = 1, n3msub
            do j = 1, n2msub
                do i = 1, n1msub
                    T(i,j,k) = T(i,j,k) + RHSJ(i,k,j)
                enddo
            enddo
        enddo

        ! Deallocate the RHSJ array which is not required any more
        deallocate(RHSJ)

        ! Nullify grid information pointer
        nullify(dx1, dx2, dx3, dmx1, dmx2, dmx3)

    end subroutine mpi_thermal_solver

end module mpi_thermal
