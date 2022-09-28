!======================================================================================================================
!> @file        module_solve_momentum.f90
!> @brief       This file contains a module of velocity solver for PaScaL_TCS.
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
!> @brief       Module for momentum solver
!> @details     This module solves momentum equation and returns velocity fileds in each direction.
!>              It uses PaScaL_TDMA library to solve the decoupled tri-diagonal systems of equations.
!>
module mpi_momentum

    use global

    implicit none

    !> @{ Main unknown arrays of velocity and its indrement
    double precision, allocatable, dimension(:,:,:) :: U,V,W
    double precision, allocatable, dimension(:,:,:) :: dU,dV,dW
    !> }@
    !> @{ Boundary condition of velocity at upper and lower walls in y-direction
    double precision, allocatable, dimension(:,:) :: UBCup_sub,VBCup_sub,WBCup_sub
    double precision, allocatable, dimension(:,:) :: UBCbt_sub,VBCbt_sub,WBCbt_sub
    !> }@

    !> @{ Array of pressure field : it is updated at mpi_pressure module
    double precision, allocatable, dimension(:,:,:) :: P
    !>@
    
    !> @{ Local arrays of momentum coefficients
    double precision, allocatable, dimension(:,:,:) :: invRho,Mu
    !>@

    contains    

    !>
    !> @brief       Assign variables for momentum solver
    !>
    subroutine mpi_momentum_allocation()
        
        use mpi_subdomain,  only : n1sub, n2sub, n3sub

        implicit none

        allocate(U(0:n1sub,0:n2sub,0:n3sub),V(0:n1sub,0:n2sub,0:n3sub),W(0:n1sub,0:n2sub,0:n3sub))
        allocate(P(0:n1sub,0:n2sub,0:n3sub))

        allocate( UBCup_sub(0:n1sub,0:n3sub),VBCup_sub(0:n1sub,0:n3sub),WBCup_sub(0:n1sub,0:n3sub))
        allocate( UBCbt_sub(0:n1sub,0:n3sub),VBCbt_sub(0:n1sub,0:n3sub),WBCbt_sub(0:n1sub,0:n3sub))

    end subroutine mpi_momentum_allocation

    !>
    !> @brief       Deallocate variables for momentum solver
    !>
    subroutine mpi_momentum_clean()

        implicit none
        
        deallocate(U,V,W)
        deallocate(P)

        deallocate( UBCup_sub,VBCup_sub,WBCup_sub)
        deallocate( UBCbt_sub,VBCbt_sub,WBCbt_sub)
        
    end subroutine mpi_momentum_clean

    !>
    !> @brief       Initialize velocity field
    !>
    subroutine mpi_momentum_initial()
        
        use mpi_subdomain,  only : n1sub, n2sub, n3sub

        implicit none

        double precision:: x1(0:n1sub),x2(0:n2sub),x3(0:n3sub)
        integer :: i,j,k

        U(0:n1sub,0:n2sub,0:n3sub)=0.d0
        V(0:n1sub,0:n2sub,0:n3sub)=0.d0
        W(0:n1sub,0:n2sub,0:n3sub)=0.d0
        P(0:n1sub,0:n2sub,0:n3sub)=0.d0

    end subroutine mpi_momentum_initial

    !>
    !> @brief       Initialize velocity boundary condition
    !>
    subroutine mpi_momentum_boundary()

        use mpi_topology,   only : comm_1d_x1, comm_1d_x2, comm_1d_x3
        use mpi_subdomain,  only : n1sub, n2sub, n3sub

        implicit none

        integer :: i, k
        
       ! Example problem : Constant temperature for upper and lower walls in y-direction
        UBCup_sub(:,:)=Uup
        VBCup_sub(:,:)=Vup
        WBCup_sub(:,:)=Wup
        
        UBCbt_sub(:,:)=Ubt
        VBCbt_sub(:,:)=Vbt
        WBCbt_sub(:,:)=Wbt

        if(comm_1d_x2%myrank==0) then
            do k=0, n3sub
                do i=0, n1sub
                    U(i,0,k)=UBCbt_sub(i,k)
                    W(i,0,k)=WBCbt_sub(i,k)
                    V(i,0,k)=VBCbt_sub(i,k)
                    V(i,1,k)=VBCbt_sub(i,k)
                enddo
            enddo
        else if(comm_1d_x2%myrank==comm_1d_x2%nprocs-1) then
            do k=0, n3sub
                do i=0, n1sub
                    U(i,n2sub,k)=UBCup_sub(i,k)
                    W(i,n2sub,k)=WBCup_sub(i,k)
                    V(i,n2sub,k)=VBCup_sub(i,k)
                enddo
            enddo
        endif    

    end subroutine mpi_momentum_boundary

    !>
    !> @brief       Assign and calculate the momentum coefficients
    !>
    subroutine mpi_momentum_coeffi(T)

        use mpi_subdomain,  only : n1sub, n2sub, n3sub

        implicit none

        double precision :: T(0:n1sub,0:n2sub,0:n3sub)
        double precision :: tmp1, tmp2, tmp3, tmp4, tmp5
        integer :: i,j,k

        allocate(invRho(0:n1sub,0:n2sub,0:n3sub),Mu(0:n1sub,0:n2sub,0:n3sub))
        do k=0,n3sub
            do j=0,n2sub
                do i=0,n1sub
                    tmp1 = DeltaT*T(i,j,k)
                    tmp2 = tmp1 * tmp1
                    tmp3 = tmp2 * tmp1
                    tmp4 = tmp3 * tmp1
                    tmp5 = tmp4 * tmp1

                    Mu(:,j,k) = (a10 + a10*(a11*tmp1 + a12*tmp2 + a13*tmp3 + a14*tmp4 + a15*tmp5)) &
                               *(d10 + d10*(d11*tmp1 + d12*tmp2 + d13*tmp3 + d14*tmp4 + d15*tmp5))/Mu0

                    invRho(:,j,k) = Rho0/(a10 + a10*(a11*tmp1 + a12*tmp2 + a13*tmp3 + a14*tmp4 + a15*tmp5))
                enddo
            enddo
        enddo

    end subroutine mpi_momentum_coeffi

    !>
    !> @brief       Update velocity field with the solved incremental velocity fields
    !>
    subroutine mpi_momentum_pseudoupdateUVW()

        use mpi_subdomain,  only : n1sub, n2sub, n3sub, n1msub, n2msub, n3msub
        use mpi_subdomain,  only : j_indexS

        implicit none
        integer :: i,j,k

        do k = 1, n3msub
            do j = 2, n2msub
                do i = 1, n1msub
                    U(i,j,k)=U(i,j,k)+dU(i,j,k)
                    V(i,j,k)=V(i,j,k)+dV(i,j,k)
                    W(i,j,k)=W(i,j,k)+dW(i,j,k)
                enddo
            enddo

            ! Special treatment for velocity field in case of lower wall. Valid for the current example case
            j=1
            do i = 1, n1msub
                U(i,j,k)=U(i,j,k)+dU(i,j,k)
                V(i,j,k)=V(i,j,k)+dV(i,j,k)*dble(2-j_indexS)
                W(i,j,k)=W(i,j,k)+dW(i,j,k)
            enddo
        enddo

        deallocate(dU)
        deallocate(dV)
        deallocate(dW)
        
    end subroutine mpi_momentum_pseudoupdateUVW

    !>
    !> @brief       deallocate the momentum coefficients
    !>
    subroutine mpi_momentum_coeffi_clean()

        implicit none

        deallocate(invRho,Mu)
        
    end subroutine mpi_momentum_coeffi_clean

    !>
    !> @brief       Main momentum solver for incremental velocity in x-direction (du)
    !> @param       T       Temperature field
    !>
    subroutine mpi_momentum_solvedU(T)

        use mpi_topology,   only : comm_1d_x1, comm_1d_x2, comm_1d_x3
        use mpi_subdomain,  only : n1sub, n2sub, n3sub, n1msub, n2msub, n3msub
        use mpi_subdomain,  only : dx1_sub,dx2_sub,dx3_sub,dmx1_sub,dmx2_sub,dmx3_sub
        use mpi_subdomain,  only : i_indexS, jC_BC

        use PaScaL_TDMA

        implicit none

        double precision, dimension(0:n1sub,0:n2sub,0:n3sub) ::  T

        type(ptdma_plan_many)     :: ptdma_plan

        !> @{ Local pointer for subdomain variables
        double precision, dimension(:), pointer :: dx1, dx2, dx3
        double precision, dimension(:), pointer :: dmx1, dmx2, dmx3
        !> @}

        !> @{ Local indexing variables
        integer :: i,j,k
        integer :: im,jm,km
        integer :: ip,jp,kp
        integer :: jup,jum
        !> @}

        !> @{ Local arrays for matrix RHS coefficients
        double precision, allocatable, dimension(:,:,:) :: ddU,RHS,RHSI,RHSJ
        double precision, allocatable, dimension(:,:,:) :: am,ac,ap
        !> @}

        !> @{ Local variables for matrix and RHS coefficients calculation
        double precision :: u1,u2,u5,u6,v3,v4,w5,w6,Tc
        double precision :: dudx1,dudx2,dudy3,dudy4,dudz5,dudz6,dvdx3,dvdx4,dwdx5,dwdx6
        double precision :: mua,mub,muc,mu3,mu4,mu5,mu6
        double precision :: invRhoc,invRhocCmu_half,viscous_u1,viscous_u2,viscous_u3,viscous_u12,viscous_u13
        double precision :: ubc_up,ubc_down
        double precision :: mACI,mAPI,mAMI,mACJ,mAPJ,mAMJ,mACK,mAPK,mAMK
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
                jm = j-1
                jp = j+1
                jum = jC_BC(jm)
                jup = jC_BC(jp)

                do i = i_indexS, n1msub

                    im = i-1
                    ip = i+1
                            
                    u1 = 0.5d0*(U(im,j,k) + U(i,j,k))            
                    u2 = 0.5d0*(U(ip,j,k) + U(i,j,k))             

                    v3 = 0.5d0*(dx1(im)*V(i,j ,k) + dx1(i)*V(im,j ,k))/dmx1(i)            
                    v4 = 0.5d0*(dx1(im)*V(i,jp,k) + dx1(i)*V(im,jp,k))/dmx1(i)
                    
                    w5 = 0.5d0*(dx1(im)*W(i,j,k ) + dx1(i)*W(im,j,k ))/dmx1(i)
                    w6 = 0.5d0*(dx1(im)*W(i,j,kp) + dx1(i)*W(im,j,kp))/dmx1(i)

                    dudx1 = (U(i,j ,k ) - U(im,j ,k ))/dx1(im)
                    dudx2 = (U(ip,j ,k ) - U(i,j ,k ))/dx1(i)
                    dudy3 = (U(i,j ,k ) - U(i,jm,k ))/dmx2(j )
                    dudy4 = (U(i,jp,k ) - U(i,j ,k ))/dmx2(jp)
                    dudz5 = (U(i,j ,k ) - U(i,j ,km))/dmx3(k )
                    dudz6 = (U(i,j ,kp) - U(i,j ,k ))/dmx3(kp)

                    dvdx3 = (V(i,j ,k) - V(im,j ,k))/dmx1(i)
                    dvdx4 = (V(i,jp,k) - V(im,jp,k))/dmx1(i)

                    dwdx5 = (W(i,j,k ) - W(im,j,k ))/dmx1(i)
                    dwdx6 = (W(i,j,kp) - W(im,j,kp))/dmx1(i)
                                        
                    mua = 0.5d0*(dx1(im)*Mu(i,jm,k ) + dx1(i)*Mu(im,jm,k ))/dmx1(i)
                    muc = 0.5d0*(dx1(im)*Mu(i,j ,k ) + dx1(i)*Mu(im,j ,k ))/dmx1(i)
                    mub = 0.5d0*(dx1(im)*Mu(i,jp,k ) + dx1(i)*Mu(im,jp,k ))/dmx1(i)
                    
                    mu3 = 0.5d0*(dx2(jm)*muc + dx2(j)*mua)/dmx2(j )
                    mu4 = 0.5d0*(dx2(jp)*muc + dx2(j)*mub)/dmx2(jp)
                    
                    mua = 0.5d0*(dx1(im)*Mu(i,j ,km) + dx1(i)*Mu(im,j ,km))/dmx1(i)
                    muc = 0.5d0*(dx1(im)*Mu(i,j ,k ) + dx1(i)*Mu(im,j ,k ))/dmx1(i)
                    mub = 0.5d0*(dx1(im)*Mu(i,j ,kp) + dx1(i)*Mu(im,j ,kp))/dmx1(i)
                    
                    mu5 = 0.5d0*(dx3(km)*muc + dx3(k)*mua)/dmx3(k )
                    mu6 = 0.5d0*(dx3(kp)*muc + dx3(k)*mub)/dmx3(kp)
                    
                    invRhoc = 0.5d0*(dx1(im)*invRho(i,j ,k )+dx1(i)*invRho(im,j ,k ) )/dmx1(i)
                    invRhocCmu_half = 0.5d0*Cmu*invRhoc
    
                    ! Viscous term
                    viscous_u1 = 1.d0*(Mu(i,j ,k )*dudx2 - Mu(im,j ,k )*dudx1)/dmx1(i)
                    viscous_u2 = 1.d0*(mu4*dudy4 - mu3*dudy3)/dx2(j )
                    viscous_u3 = 1.d0*(mu6*dudz6 - mu5*dudz5)/dx3(k )            
                    viscous_u12 = 1.d0*(mu4*dvdx4 - mu3*dvdx3)/dx2(j)
                    viscous_u13 = 1.d0*(mu6*dwdx6 - mu5*dwdx5)/dx3(k)

                    RHS(i,j,k) = invRhocCmu_half*(2.d0*viscous_u1 + viscous_u2 + viscous_u3 + viscous_u12 + viscous_u13)

                    ! Pressure 
                    RHS(i,j,k) = RHS(i,j,k) &
                                - Cmp*invRhoc*(P(i,j,k) - P(im,j,k))/dmx1(i)

                    ! Buoyancy term (we use Cmt*theta for the buoyancy term)
                    Tc = 0.5d0*(dx1(im)*T(i,j,k) + dx1(i)*T(im,j,k))/dmx1(i)
                    !THETAy = Cmt*((Tc(i) - Thetam))**b
                    !THETAy = Cmt*(Tc(i))
                    RHS(i,j,k) = RHS(i,j,k) &
                                + Cmt*(Tc + a12pera11*Tc**2.*DeltaT)*invRhoc

                    ! Upper B.C.
                    ! From convection term
                    ubc_down = 0.25d0*v3/dmx2(j)*UBCbt_sub(i,k) - 0.25d0*dudy3*(dx1(i)*VBCbt_sub(im,k) + dx1(im)*VBCbt_sub(i,k))/dmx1(i)*0.5d0
                    ubc_up  = -0.25d0*v4/dmx2(jp)*UBCup_sub(i,k)- 0.25d0*dudy4*(dx1(i)*VBCup_sub(im,k) + dx1(im)*VBCup_sub(i,k))/dmx1(i)*0.5d0

                    ! From diffusion term
                    ubc_down = ubc_down + invRhocCmu_half/dx2(j)*mu3/dmx2(j )*UBCbt_sub(i,k)
                    ubc_down = ubc_down - invRhocCmu_half/dx2(j)*mu3/dmx1(i)*(VBCbt_sub(i,k) - VBCbt_sub(im,k))

                    ubc_up   = ubc_up   + invRhocCmu_half/dx2(j)*mu4/dmx2(jp)*UBCup_sub(i,k)
                    ubc_up   = ubc_up   + invRhocCmu_half/dx2(j)*mu4/dmx1(i)*(VBCup_sub(i,k) - VBCup_sub(im,k))

                    RHS(i,j,k) = RHS(i,j,k) &
                                + dble(1.d0 - jum)*ubc_down + dble(1.d0 - jup)*ubc_up

                    !M11Un
                    !X-direction
                    mACI =  invRhocCmu_half/dmx1(i)*(Mu(i,j ,k )/dx1(i) + Mu(i,j ,k )/dx1(im))* 2.d0 &
                            + 0.25d0*(dudx2*0.5d0 + dudx1*0.5d0 - u2/dx1(i) + u1/dx1(im) )
                    mAPI = -invRhocCmu_half/dmx1(i)*Mu(i,j ,k )/dx1(i)* 2.d0 &
                            + 0.25d0*( u2/dx1(i) + dudx2*0.5d0) 
                    mAMI = -invRhocCmu_half/dmx1(i)*Mu(i,j ,k )/dx1(im)* 2.d0 &
                            + 0.25d0*(-u1/dx1(im) + dudx1*0.5d0)

                    !Y-direction
                    mACJ =  invRhocCmu_half/dx2(j)*(mu4/dmx2(jp) + mu3/dmx2(j)) &
                            + 0.25d0*(-v4/dmx2(jp) + v3/dmx2(j ))
                    mAPJ = -invRhocCmu_half/dx2(j)*mu4/dmx2(jp) &
                            + 0.25d0*( v4/dmx2(jp))
                    mAMJ = -invRhocCmu_half/dx2(j)*mu3/dmx2(j)  &
                            + 0.25d0*(-v3/dmx2(j))
                    mAPJ = mAPJ*dble(jup)
                    mAMJ = mAMJ*dble(jum)
                    
                    !Z-direction
                    mACK =  invRhocCmu_half/dx3(k)*(mu6/dmx3(kp) + mu5/dmx3(k)) &
                            + 0.25d0*(-w6/dmx3(kp) + w5/dmx3(k))
                    mAPK = -invRhocCmu_half/dx3(k)*mu6/dmx3(kp) &
                            + 0.25d0*( w6/dmx3(kp))
                    mAMK = -invRhocCmu_half/dx3(k)*mu5/dmx3(k ) &
                            + 0.25d0*(-w5/dmx3(k))

                    RHS(i,j,k) = RHS(i,j,k) &
                                - ( mAPI*U(ip,j ,k ) + mACI*U(i,j,k) + mAMI*U(im,j ,k ) &
                                  + mAPJ*U(i, jp,k ) + mACJ*U(i,j,k) + mAMJ*U(i, jm,k ) &
                                  + mAPK*U(i, j ,kp) + mACK*U(i,j,k) + mAMK*U(i, j, km) )

                    RHS(i,j,k) = RHS(i,j,k) &
                                - ( 0.25d0*(v4*dble(jup)*dudy4 + v3*dble(jum)*dudy3) &
                                - invRhocCmu_half/dx2(j)*(mu4*dvdx4*dble(jup) - mu3*dvdx3*dble(jum)) )

                    RHS(i,j,k) = RHS(i,j,k) &
                                - ( 0.25d0*(w6*dudz6 + w5*dudz5) &
                                - invRhocCmu_half/dx3(k)*(mu6*dwdx6 - mu5*dwdx5) )
                        
                    ! Coefficients of tridagonal A-matrix
                    ac(i,j,k) = mACK*dt + 1.d0
                    ap(i,j,k) = mAPK*dt
                    am(i,j,k) = mAMK*dt
                    RHS(i,j,k)= RHS(i,j,k)*dt
                enddo
            enddo
        enddo

        ! 1st stage : solve TDM in z-direction
        call PaScaL_TDMA_plan_many_create(ptdma_plan, (n1msub)*(n2msub), comm_1d_x3%myrank, comm_1d_x3%nprocs, comm_1d_x3%mpi_comm)
        call PaScaL_TDMA_many_solve_cycle(ptdma_plan, am, ac, ap, RHS,(n1msub)*(n2msub),(n3msub))
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x3%nprocs)

        ! Deallocate A-matrix for the 2nd stage
        deallocate(am,ac,ap)

        ! Transpose the array (i,j,k) to array (j,k,i) to solve the tridagonal systems in x-direction
        allocate(RHSI(1:n2msub,1:n3msub,1:n1msub))
        do k = 1, n3msub
            do j = 1, n2msub
                do i = i_indexS, n1msub
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
                do i = i_indexS, n1msub

                    im = i- 1
                    ip = i+1

                    u1 = 0.5d0*(U(im,j,k) + U(i,j,k))            
                    u2 = 0.5d0*(U(ip,j,k) + U(i,j,k))             

                    dudx1 = (U(i,j ,k ) - U(im,j ,k ))/dx1(im)
                    dudx2 = (U(ip,j ,k ) - U(i,j ,k ))/dx1(i)
                    
                    invRhoc = 0.5d0*(dx1(im)*invRho(i,j ,k )+dx1(i)*invRho(im,j ,k ) )/dmx1(i)
                    invRhocCmu_half = 0.5d0*Cmu*invRhoc

                    !X-direction
                    mACI =  invRhocCmu_half/dmx1(i)*(Mu(i,j ,k )/dx1(i) + Mu(i,j ,k )/dx1(im))* 2.d0 &
                            + 0.25d0*(dudx2*0.5d0 + dudx1*0.5d0 - u2/dx1(i) + u1/dx1(im) )
                    mAPI = -invRhocCmu_half/dmx1(i)*Mu(i,j ,k )/dx1(i)* 2.d0 &
                            + 0.25d0*( u2/dx1(i) + dudx2*0.5d0) 
                    mAMI = -invRhocCmu_half/dmx1(i)*Mu(i,j ,k )/dx1(im)* 2.d0 &
                            + 0.25d0*(-u1/dx1(im) + dudx1*0.5d0)

                    ac(j,k,i)=mACI*dt + 1.d0
                    ap(j,k,i)=mAPI*dt
                    am(j,k,i)=mAMI*dt
                enddo
            enddo
        enddo

        ! 2nd stage : solve TDM in x-direction
        call PaScaL_TDMA_plan_many_create(ptdma_plan, (n2msub)*(n3msub), comm_1d_x1%myrank, comm_1d_x1%nprocs, comm_1d_x1%mpi_comm)
        call PaScaL_TDMA_many_solve_cycle(ptdma_plan, am, ac, ap, RHSI,(n2msub)*(n3msub),(n1msub))
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x1%nprocs)

        ! Deallocate A-matrix for the 3rd stage
        deallocate(am,ac,ap)

        ! Transpose the array (j,k,i) to array (i,k,j) to solve the tridagonal systems in y-direction
        allocate(RHSJ(1:n1msub,1:n3msub,1:n2msub))
        do j = 1, n2msub
            do i = i_indexS, n1msub
                do k = 1, n3msub
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
                jum = jC_BC(jm)
                jup = jC_BC(jp)
                
                do i = i_indexS, n1msub
                    im = i- 1
                    ip = i+1

                    v3 = 0.5d0*(dx1(im)*V(i,j ,k) + dx1(i)*V(im,j ,k))/dmx1(i)
                    v4 = 0.5d0*(dx1(im)*V(i,jp,k) + dx1(i)*V(im,jp,k))/dmx1(i)
                
                    mua = 0.5d0*(dx1(im)*Mu(i,jm,k ) + dx1(i)*Mu(im,jm,k ))/dmx1(i)
                    muc = 0.5d0*(dx1(im)*Mu(i,j ,k ) + dx1(i)*Mu(im,j ,k ))/dmx1(i)
                    mub = 0.5d0*(dx1(im)*Mu(i,jp,k ) + dx1(i)*Mu(im,jp,k ))/dmx1(i)
                    mu3 = 0.5d0*(dx2(jm)*muc + dx2(j)*mua)/dmx2(j )
                    mu4 = 0.5d0*(dx2(jp)*muc + dx2(j)*mub)/dmx2(jp)
                    
                    invRhoc = 0.5d0*(dx1(im)*invRho(i,j ,k )+dx1(i)*invRho(im,j ,k ) )/dmx1(i)
                    invRhocCmu_half = 0.5d0*Cmu*invRhoc

                    !Y-direction
                    mACJ =  invRhocCmu_half/dx2(j)*(mu4/dmx2(jp) + mu3/dmx2(j)) &
                            + 0.25d0*(-v4/dmx2(jp) + v3/dmx2(j ))
                    mAPJ = -invRhocCmu_half/dx2(j)*mu4/dmx2(jp) &
                            + 0.25d0*( v4/dmx2(jp))
                    mAMJ = -invRhocCmu_half/dx2(j)*mu3/dmx2(j) &
                            + 0.25d0*(-v3/dmx2(j))
                    mAPJ = mAPJ*dble(jup)
                    mAMJ = mAMJ*dble(jum)

                    ac(i,k,j)=mACJ*dt + 1.d0
                    ap(i,k,j)=mAPJ*dt
                    am(i,k,j)=mAMJ*dt
                enddo
            enddo
        enddo

        ! 3rd stage : solve TDM in y-direction
        call PaScaL_TDMA_plan_many_create(ptdma_plan, (n3msub)*(n1msub), comm_1d_x2%myrank, comm_1d_x2%nprocs, comm_1d_x2%mpi_comm)
        call PaScaL_TDMA_many_solve(ptdma_plan, am, ac, ap, RHSJ,(n3msub)*(n1msub),(n2msub))
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x2%nprocs)        

        ! Deallocate A-matrix for the 3rd stage
        deallocate(am,ac,ap)

        ! Update the velocity increments with the transpose of the solution array (i,k,j)
        allocate(dU(0:n1sub,0:n2sub,0:n3sub))
        do k = 1, n3msub
            do j = 1, n2msub
                do i = 1, n1msub
                    dU(i,j,k)=RHSJ(i,k,j)
                enddo
            enddo
        enddo

        ! Deallocate the RHSJ array which is not required any more
        deallocate(RHSJ)

        ! Nullify grid information pointer
        nullify(dx1, dx2, dx3, dmx1, dmx2, dmx3)

    end subroutine mpi_momentum_solvedU

    !>
    !> @brief       Main momentum solver for incremental velocity in x-direction (dv)
    !> @param       T       Temperature field
    !>
    subroutine mpi_momentum_solvedV(T)

        use mpi_topology,   only : comm_1d_x1, comm_1d_x2, comm_1d_x3
        use mpi_subdomain,  only : j_indexS, jS_BC
        use mpi_subdomain,  only : n1sub, n2sub, n3sub, n1msub, n2msub, n3msub
        use mpi_subdomain,  only : dx1_sub,dx2_sub,dx3_sub,dmx1_sub,dmx2_sub,dmx3_sub

        use PaScaL_TDMA

        implicit none

        double precision, dimension(0:n1sub,0:n2sub,0:n3sub) ::  T

        type(ptdma_plan_many)     :: ptdma_plan

        !> @{ Local pointer for subdomain variables
        double precision, dimension(:), pointer :: dx1, dx2, dx3
        double precision, dimension(:), pointer :: dmx1, dmx2, dmx3
        !> @}

        !> @{ Local indexing variables
        integer :: i,j,k
        integer :: im,jm,km
        integer :: ip,jp,kp
        integer :: jvp,jvm
        !> @}

        !> @{ Local arrays for matrix and RHS coefficients
        double precision, allocatable, dimension(:,:,:) :: ddV,RHS,RHSI,RHSJ
        double precision, allocatable, dimension(:,:,:) :: am,ac,ap
        !> @}

        !> @{ Local variables for matrix and RHS coefficients calculation
        double precision :: u1,u2,v3,v4,w5,w6,Tc
        double precision :: dudy1,dudy2,dvdx1,dvdx2,dvdy3,dvdy4,dvdz5,dvdz6,dwdy5,dwdy6
        double precision :: ddu1,ddu2,ddudy1,ddudy2
        double precision :: mua,muc,mub,mu1,mu2,mu5,mu6
        double precision :: invRhoc,invRhocCmu_half,viscous_v1,viscous_v2,viscous_v3,viscous_v21,viscous_v23
        double precision :: vbc_up,vbc_down
        double precision :: mACI,mAPI,mAMI,mACJ,mAPJ,mAMJ,mACK,mAPK,mAMK
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

            kp = k+1
            km = k-1

            do j =  j_indexS, n2msub

                jp = j + 1
                jm = j - 1
                jvp = jS_BC(jm)
                jvm = jS_BC(jp)

                do i = 1, n1msub

                    ip = i + 1
                    im = i - 1

                    v3 = 0.5d0*(V(i,jm,k) + V(i,j ,k))
                    v4 = 0.5d0*(V(i,j ,k) + V(i,jp,k))

                    u1 = (dx2(jm)*U(i, j,k) + dx2(j)*U(i, jm,k))/dmx2(j)*0.5d0
                    u2 = (dx2(jm)*U(ip,j,k) + dx2(j)*U(ip,jm,k))/dmx2(j)*0.5d0

                    w5 = (dx2(jm)*W(i,j,k ) + dx2(j)*W(i,jm,k ))/dmx2(j)*0.5d0
                    w6 = (dx2(jm)*W(i,j,kp) + dx2(j)*W(i,jm,kp))/dmx2(j)*0.5d0
                    
                    ddu1 = (dx2(jm)*dU(i,j,k) + dx2(j)*dU(i,jm,k))/dmx2(j)*0.5d0
                    ddu2 = (dx2(jm)*dU(ip,j,k) + dx2(j)*dU(ip,jm,k))/dmx2(j)*0.5d0

                    
                    dvdx1 = (V(i, j,k) - V(im,j,k))/dmx1(i)
                    dvdx2 = (V(ip,j,k) - V(i, j,k))/dmx1(ip)
                    dvdy3 = (V(i,j ,k) - V(i,jm,k))/dx2(jm)
                    dvdy4 = (V(i,jp,k) - V(i,j ,k))/dx2(j )
                    dvdz5 = (V(i,j,k ) - V(i,j,km))/dmx3(k )
                    dvdz6 = (V(i,j,kp) - V(i,j,k ))/dmx3(kp)
                    

                    dudy1 = (U(i, j,k) - U(i, jm,k))/dmx2(j)
                    dudy2 = (U(ip,j,k) - U(ip,jm,k))/dmx2(j)

                    dwdy5 = (W(i,j,k ) - W(i,jm,k ))/dmx2(j)
                    dwdy6 = (W(i,j,kp) - W(i,jm,kp))/dmx2(j)

                    dvdx1 = (V(i, j,k) - V(im,j,k))/dmx1(i)
                    dvdx2 = (V(ip,j,k) - V(i, j,k))/dmx1(ip)
                    ddudy1 = (dU(i, j,k) - dU(i, jm,k))/dmx2(j)
                    ddudy2 = (dU(ip,j,k) - dU(ip,jm,k))/dmx2(j) 

                    mua = 0.5d0*(dx2(jm)*Mu(im,j ,k ) + dx2(j)*Mu(im,jm,k ))/dmx2(j)
                    muc = 0.5d0*(dx2(jm)*Mu(i,j ,k ) + dx2(j)*Mu(i,jm,k ))/dmx2(j)
                    mub = 0.5d0*(dx2(jm)*Mu(ip,j ,k ) + dx2(j)*Mu(ip,jm,k ))/dmx2(j)
                    mu1 = 0.5d0*(dx1(im)*muc + dx1(i)*mua)/dmx1(i)
                    mu2 = 0.5d0*(dx1(ip)*muc + dx1(i)*mub)/dmx1(ip)

                    mua = 0.5d0*(dx2(jm)*Mu(i,j ,km) + dx2(j)*Mu(i,jm,km))/dmx2(j)
                    muc = 0.5d0*(dx2(jm)*Mu(i,j ,k ) + dx2(j)*Mu(i,jm,k ))/dmx2(j)
                    mub = 0.5d0*(dx2(jm)*Mu(i,j ,kp) + dx2(j)*Mu(i,jm,kp))/dmx2(j)
                    mu5 = 0.5d0*(dx3(km)*muc + dx3(k)*mua)/dmx3(k )
                    mu6 = 0.5d0*(dx3(kp)*muc + dx3(k)*mub)/dmx3(kp)

                    invRhoc = 0.5d0*(dx2(jm)*invRho(i,j ,k ) + dx2(j )*invRho(i,jm,k ))/dmx2(j)
                    invRhocCmu_half = 0.5d0*Cmu*invRhoc

                    ! Viscous term
                    viscous_v1 = 1.d0*(mu2*dvdx2 - mu1*dvdx1)/dx1(i)
                    viscous_v2 = 1.d0*(Mu(i,j ,k )*dvdy4 - Mu(i,jm,k )*dvdy3)/dmx2(j)
                    viscous_v3 = 1.d0*(mu6*dvdz6 - mu5*dvdz5)/dx3(k )
                    viscous_v21 = 1.d0*(mu2*dudy2 - mu1*dudy1)/dx1(i)
                    viscous_v23 = 1.d0*(mu6*dwdy6 - mu5*dwdy5)/dx3(k)

                    RHS(i,j,k) = invRhocCmu_half*(viscous_v1 + 2.*viscous_v2 + viscous_v3 + viscous_v21 + viscous_v23)          ! GE

                    ! Pressure term (invRhoc(i) is included for the NOB case)
                    RHS(i,j,k) = RHS(i,j,k)     &
                                - Cmp*invRhoc*(P(i,j,k) - P(i,jm,k))/dmx2(j)


                    ! Buoyancy term (we use Cmt*theta for the buoyancy term)
                    ! Tc = 0.5d0*(dx2(jm)*T(i,j,k) + dx2(j)*T(i,jm,k))/dmx2(j)
                    ! !THETAy = Cmt*((Tc - Thetam))**b
                    ! !THETAy = Cmt*(Tc)
                    ! RHS(i,j,k) = RHS(i,j,k)     &
                    !             + Cmt*(Tc + a12pera11*Tc**2.*DeltaT)*invRhoc
            
                    ! Upper B.C.
                    vbc_down =  ( 0.25d0/dx2(jm)*v3*VBCbt_sub(i,k) &
                                - 0.25d0*0.5d0*dvdy3*VBCbt_sub(i,k) )

                    vbc_up =    (-0.25d0/dx2(j )*v4*VBCup_sub(i,k) &
                                - 0.25d0*0.5d0*dvdy4*VBCup_sub(i,k) )

                    vbc_down = vbc_down + invRhocCmu_half/dmx2(j)*Mu(i,jm,k )/dx2(jm)*VBCbt_sub(i,k)*2.d0           !'*2.' is needed since the modified GE
                    vbc_up   = vbc_up   + invRhocCmu_half/dmx2(j)*Mu(i,j ,k )/dx2(j )*VBCup_sub(i,k)*2.d0

                    RHS(i,j,k) = RHS(i,j,k) &
                                + dble(1.d0 - jvm)*vbc_down + dble(1.d0 - jvp)*vbc_up

                    !----M22Vn
                    !X-direction
                    mACI =  invRhocCmu_half/dx1(i)*(mu2/dmx1(ip) + mu1/dmx1(i)) &
                            + 0.25d0*(-u2/dmx1(ip)+ u1/dmx1(i))
                    mAPI = -invRhocCmu_half/dx1(i)*mu2/dmx1(ip) &
                            + 0.25d0*( u2/dmx1(ip))
                    mAMI = -invRhocCmu_half/dx1(i)*mu1/dmx1(i)  &
                            + 0.25d0*(-u1/dmx1(i))

                    !Y-direction
                    mACJ = invRhocCmu_half/dmx2(j)*(Mu(i,j ,k )/dx2(j) + Mu(i,jm,k )/dx2(jm))*2.d0 &
                            + 0.25d0*(dvdy4*0.5d0 + dvdy3*0.5d0-v4/dx2(j) + v3/dx2(jm) )
                    mAPJ = -invRhocCmu_half/dmx2(j)*Mu(i,j ,k )/dx2(j )*2.d0 &
                            + 0.25d0*( v4/dx2(j ) + dvdy4*0.5d0)
                    mAPJ = mAPJ*dble(jvp)
                    mAMJ = -invRhocCmu_half/dmx2(j)*Mu(i,jm,k )/dx2(jm)*2.d0 &
                            + 0.25d0*(-v3/dx2(jm) + dvdy3*0.5d0)
                    mAMJ = mAMJ*dble(jvm)

                    !Z-direction
                    mACK = invRhocCmu_half/dx3(k)*(mu6/dmx3(kp) + mu5/dmx3(k)) &
                            + 0.25d0*(-w6/dmx3(kp) + w5/dmx3(k))
                    mAPK = -invRhocCmu_half/dx3(k)*mu6/dmx3(kp) &
                            + 0.25d0*( w6/dmx3(kp))
                    mAMK = -invRhocCmu_half/dx3(k)*mu5/dmx3(k ) &
                            + 0.25d0*(-w5/dmx3(k ))

                    RHS(i,j,k) = RHS(i,j,k) &
                                -(mAPI*V(ip,j, k ) + mACI*V(i,j,k) + mAMI*V(im,j, k ) &
                                + mAPJ*V(i, jp,k ) + mACJ*V(i,j,k) + mAMJ*V(i, jm,k ) &
                                + mAPK*V(i, j, kp) + mACK*V(i,j,k) + mAMK*V(i, j, km) )

                    RHS(i,j,k) = RHS(i,j,k) &
                                -(0.25d0*(u2*dvdx2 + u1*dvdx1) &
                                - invRhocCmu_half/dx1(i)*(mu2*dudy2 - mu1*dudy1) )
                        
                    RHS(i,j,k) = RHS(i,j,k) &
                                -(0.25d0*(w6*dvdz6 + w5*dvdz5) &
                                - invRhocCmu_half/dx3(k)*(mu6*dwdy6 - mu5*dwdy5))
                
                    RHS(i,j,k) = RHS(i,j,k) &
                                -(0.25d0*(ddu2*dvdx2+ddu1*dvdx1) &
                                - invRhocCmu_half/dx1(i)*(mu2*ddudy2 - mu1*ddudy1) )

                    ! Coefficients of tridagonal A-matrix
                    ac(i,j,k) =mACK*dt + 1.d0
                    ap(i,j,k) =mAPK*dt
                    am(i,j,k) =mAMK*dt
                    RHS(i,j,k)=RHS(i,j,k)*dt
                enddo
            end do
        end do

        ! 1st stage : solve TDM in z-direction
        call PaScaL_TDMA_plan_many_create(ptdma_plan, (n1msub)*(n2msub), comm_1d_x3%myrank, comm_1d_x3%nprocs, comm_1d_x3%mpi_comm)
        call PaScaL_TDMA_many_solve_cycle(ptdma_plan, am, ac, ap, RHS,(n1msub)*(n2msub),(n3msub))
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x3%nprocs)

        ! Deallocate A-matrix for the 2nd stage
        deallocate(am,ac,ap)

        ! Transpose the array (i,j,k) to array (j,k,i) for the tridagonal systems in x-direction
        allocate(RHSI(1:n2msub,1:n3msub,1:n1msub))
        do k = 1, n3msub
            do j = 1, n2msub
                do i = 1, n1msub
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
            do j = j_indexS, n2msub
                jp = j + 1
                jm = j - 1
                jvp = jS_BC(jm)
                jvm = jS_BC(jp)

                do i = 1, n1msub

                    ip = i + 1
                    im = i - 1

                    u1 = (dx2(jm)*U(i, j,k) + dx2(j)*U(i, jm,k))/dmx2(j)*0.5d0
                    u2 = (dx2(jm)*U(ip,j,k) + dx2(j)*U(ip,jm,k))/dmx2(j)*0.5d0
                
                    invRhoc = 0.5d0*(dx2(jm)*invRho(i,j ,k ) + dx2(j )*invRho(i,jm,k ))/dmx2(j)
                    invRhocCmu_half = 0.5d0*Cmu*invRhoc

                    ! M22Vn:X-Direction
                    mACI =  invRhocCmu_half/dx1(i)*(mu2/dmx1(ip) + mu1/dmx1(i)) &
                            + 0.25d0*(-u2/dmx1(ip)+ u1/dmx1(i))
                    mAPI = -invRhocCmu_half/dx1(i)*mu2/dmx1(ip) &
                            + 0.25d0*( u2/dmx1(ip))
                    mAMI = -invRhocCmu_half/dx1(i)*mu1/dmx1(i) &
                            + 0.25d0*(-u1/dmx1(i))
                            
                    ac(j,k,i)=mACI*dt + 1.d0
                    ap(j,k,i)=mAPI*dt
                    am(j,k,i)=mAMI*dt
                enddo
            enddo
        enddo
        
        ! 2nd stage : solve TDM in x-direction
        call PaScaL_TDMA_plan_many_create(ptdma_plan, (n2msub)*(n3msub), comm_1d_x1%myrank, comm_1d_x1%nprocs, comm_1d_x1%mpi_comm)
        call PaScaL_TDMA_many_solve_cycle(ptdma_plan, am, ac, ap, RHSI,(n2msub)*(n3msub),(n1msub))
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x1%nprocs)

        ! Deallocate A-matrix for the 3rd stage
        deallocate(am,ac,ap)

        ! Transpose the array (j,k,i) to array (i,k,j) to solve the tridagonal systems in y-direction
        allocate(RHSJ(1:n1msub,1:n3msub,1:n2sub-j_indexS))
        do k = 1, n3msub
            do j = j_indexS, n2msub
                do i = 1, n1msub
                    RHSJ(i,k,j-j_indexS+1)=RHSI(j,k,i)
                enddo
            enddo
        enddo

        ! Deallocate the RHSI array which is not required any more
        deallocate(RHSI)

        ! Re-allocate arrays for A matrix coefficient in the third stage
        allocate(am(1:n1msub,1:n3msub,1:n2sub-j_indexS),ac(1:n1msub,1:n3msub,1:n2sub-j_indexS),ap(1:n1msub,1:n3msub,1:n2sub-j_indexS))

        ! 3rd stage : build A matrix and RHS to solve the tridagonal systems in y-direction  
        do k = 1, n3msub
            do j = j_indexS, n2msub
                jm=j-1
                jp=j+1
                jvp = jS_BC(jm)
                jvm = jS_BC(jp)

                do i = 1, n1msub

                    ip = i + 1
                    im = i - 1

                    v3 = 0.5d0*(V(i,jm,k) + V(i,j ,k))
                    v4 = 0.5d0*(V(i,j ,k) + V(i,jp,k))

                    dvdy3 = (V(i,j ,k) - V(i,jm,k))/dx2(jm)
                    dvdy4 = (V(i,jp,k) - V(i,j ,k))/dx2(j )

                    invRhoc = 0.5d0*(dx2(jm)*invRho(i,j ,k ) + dx2(j )*invRho(i,jm,k ))/dmx2(j)
                    invRhocCmu_half = 0.5d0*Cmu*invRhoc

                    ! M22Vn:Y-DIRECTION
                    mACJ =  invRhocCmu_half/dmx2(j)*(Mu(i,j ,k )/dx2(j) + Mu(i,jm,k )/dx2(jm))*2.d0 &
                            + 0.25d0*(dvdy4*0.5d0 + dvdy3*0.5d0-v4/dx2(j) + v3/dx2(jm) )
                    mAPJ = -invRhocCmu_half/dmx2(j)*Mu(i,j ,k )/dx2(j )*2.d0 &
                            + 0.25d0*( v4/dx2(j ) + dvdy4*0.5d0)
                    mAPJ = mAPJ*dble(jvp)
                    mAMJ = -invRhocCmu_half/dmx2(j)*Mu(i,jm,k )/dx2(jm)*2.d0 &
                            + 0.25d0*(-v3/dx2(jm) + dvdy3*0.5d0)
                    mAMJ = mAMJ*dble(jvm)

                    ac(i,k,j-j_indexS+1)=mACJ*dt + 1.d0
                    ap(i,k,j-j_indexS+1)=mAPJ*dt
                    am(i,k,j-j_indexS+1)=mAMJ*dt
                enddo
            enddo
        enddo

        ! 3rd stage : solve TDM in y-direction
        call PaScaL_TDMA_plan_many_create(ptdma_plan, (n3msub)*(n1msub), comm_1d_x2%myrank, comm_1d_x2%nprocs, comm_1d_x2%mpi_comm)
        call PaScaL_TDMA_many_solve(ptdma_plan, am, ac, ap, RHSJ,(n3msub)*(n1msub),(n2sub-j_indexS))
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x2%nprocs)

        ! Deallocate A-matrix for the 3rd stage
        deallocate(am,ac,ap)

        ! Update the velocity increments with the transpose of the solution array (i,k,j)
        allocate(dV(0:n1sub,0:n2sub,0:n3sub))
        dV(0:n1sub,0:1,0:n3sub)=0.d0
        dV(0:n1sub,n2sub,0:n3sub)=0.d0
        do k = 1, n3msub
            do j = j_indexS, n2msub
                do i = 1, n1msub
                    dV(i,j,k)=RHSJ(i,k,j-j_indexS+1)
                enddo
            enddo
        enddo

        ! Deallocate the RHSJ array which is not required any more
        deallocate(RHSJ)

        ! Nullify grid information pointer
        nullify(dx1, dx2, dx3, dmx1, dmx2, dmx3)

    end subroutine mpi_momentum_solvedV
    
    !>
    !> @brief       Main momentum solver for incremental velocity in z-direction (dw)
    !> @param       T       Temperature field
    !>
    subroutine mpi_momentum_solvedW(T)

        use mpi_topology,   only : comm_1d_x1, comm_1d_x2, comm_1d_x3
        use mpi_subdomain,  only : k_indexS, jC_BC
        use mpi_subdomain,  only : n1sub, n2sub, n3sub, n1msub, n2msub, n3msub
        use mpi_subdomain,  only : dx1_sub,dx2_sub,dx3_sub,dmx1_sub,dmx2_sub,dmx3_sub

        use PaScaL_TDMA

        implicit none

        double precision, dimension(0:n1sub,0:n2sub,0:n3sub) ::  T

        type(ptdma_plan_many)     :: ptdma_plan

        !> @{ Local pointer for subdomain variables
        double precision, dimension(:), pointer :: dx1, dx2, dx3
        double precision, dimension(:), pointer :: dmx1, dmx2, dmx3
        !> @}

        !> @{ Local indexing variables
        integer :: i,j,k
        integer :: im,jm,km
        integer :: ip,jp,kp
        integer :: jwp,jwm
        !> @}

        !> @{ Local arrays for matrix RHS coefficients
        double precision, allocatable, dimension(:,:,:) :: ddW,RHS,RHSI,RHSJ
        double precision, allocatable, dimension(:,:,:) :: am,ac,ap
        !> @}

        !> @{ Local variables for matrix and RHS coefficients calculation
        double precision :: u1,u2,v3,v4,w6,w5
        double precision :: dudz1,dudz2,dvdz3,dvdz4,dwdx1,dwdx2,dwdy3,dwdy4,dwdz5,dwdz6
        double precision :: ddu1,ddu2,ddudz1,ddudz2
        double precision :: ddv3,ddv4,ddvdz3,ddvdz4
        double precision :: mua,muc,mub,mu1,mu2,mu3,mu4,invRhoc,invRhocCmu_half
        double precision :: viscous_w1,viscous_w2,viscous_w3,viscous_w31,viscous_w32
        double precision :: zbc_up,zbc_down
        double precision :: mACI,mAPI,mAMI,mACJ,mAPJ,mAMJ,mACK,mAPK,mAMK
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
        do k = k_indexS, n3msub
            kp = k+1
            km = k-1
            do j = 1, n2msub
                jp = j + 1
                jm = j - 1
                jwp = jC_BC(jp)
                jwm = jC_BC(jm)

                do i = 1, n1msub

                    ip = i + 1
                    im = i - 1

                    w5 = 0.5d0*(W(i,j,k) + W(i,j,km))
                    w6 = 0.5d0*(W(i,j,k) + W(i,j,kp))
                    
                    u1 = (dx3(km)*U(i, j,k) + dx3(k)*U(i, j,km))/dmx3(k)*0.5d0
                    u2 = (dx3(km)*U(ip,j,k) + dx3(k)*U(ip,j,km))/dmx3(k)*0.5d0
                    v3 = (dx3(km)*V(i,j, k) + dx3(k)*V(i,j, km))/dmx3(k)*0.5d0
                    v4 = (dx3(km)*V(i,jp,k) + dx3(k)*V(i,jp,km))/dmx3(k)*0.5d0
                    
                    dwdx1 = (W(i, j,k) - W(im,j,k))/dmx1(i)
                    dwdx2 = (W(ip,j,k) - W(i, j,k))/dmx1(ip)
                    dwdy3 = (W(i,j, k) - W(i,jm,k))/dmx2(j )
                    dwdy4 = (W(i,jp,k) - W(i,j, k))/dmx2(jp)
                    dwdz5 = (W(i,j,k ) - W(i,j,km))/dx3(km)
                    dwdz6 = (W(i,j,kp) - W(i,j,k ))/dx3(k )
                    
                    dudz1 = (U(i, j,k) - U(i, j,km))/dmx3(k)
                    dudz2 = (U(ip,j,k) - U(ip,j,km))/dmx3(k)
                                        
                    dvdz3 = (V(i,j ,k) - V(i,j ,km))/dmx3(k)
                    dvdz4 = (V(i,jp,k) - V(i,jp,km))/dmx3(k)
                    
                    ddudz1 = (dU(i, j,k) - dU(i, j,km))/dmx3(k)
                    ddudz2 = (dU(ip,j,k) - dU(ip,j,km))/dmx3(k) 
                    ddvdz3 = (dV(i,j ,k) - dV(i,j ,km))/dmx3(k)
                    ddvdz4 = (dV(i,jp,k) - dV(i,jp,km))/dmx3(k)
                    
                    ddu1 = (dx3(km)*dU(i, j,k) + dx3(k)*dU(i, j,km))/dmx3(k)*0.5d0
                    ddu2 = (dx3(km)*dU(ip,j,k) + dx3(k)*dU(ip,j,km))/dmx3(k)*0.5d0
                    ddv3 = (dx3(km)*dV(i,j, k) + dx3(k)*dV(i,j, km))/dmx3(k)*0.5d0*dble(jwm)
                    ddv4 = (dx3(km)*dV(i,jp,k) + dx3(k)*dV(i,jp,km))/dmx3(k)*0.5d0*dble(jwp)
                    
                    mua = 0.5d0*(dx3(km)*Mu(im,j ,k ) + dx3(k)*Mu(im,j ,km))/dmx3(k)
                    muc = 0.5d0*(dx3(km)*Mu(i, j ,k ) + dx3(k)*Mu(i, j ,km))/dmx3(k)
                    mub = 0.5d0*(dx3(km)*Mu(ip,j ,k ) + dx3(k)*Mu(ip,j ,km))/dmx3(k)
                    mu1 = 0.5d0*(dx1(im)*muc + dx1(i)*mua)/dmx1(i)
                    mu2 = 0.5d0*(dx1(ip)*muc + dx1(i)*mub)/dmx1(ip)
                    
                    mua = 0.5d0*(dx3(km)*Mu(i,jm,k ) + dx3(k)*Mu(i,jm,km))/dmx3(k)
                    muc = 0.5d0*(dx3(km)*Mu(i,j ,k ) + dx3(k)*Mu(i,j ,km))/dmx3(k)
                    mub = 0.5d0*(dx3(km)*Mu(i,jp,k ) + dx3(k)*Mu(i,jp,km))/dmx3(k)
                    mu3 = 0.5d0*(dx2(jm)*muc + dx2(j)*mua)/dmx2(j )
                    mu4 = 0.5d0*(dx2(jp)*muc + dx2(j)*mub)/dmx2(jp)
                    
                    invRhoc = 0.5d0*(dx3(km)*invRho(i,j ,k ) + dx3(k )*invRho(i,j ,km) )/ dmx3(k)
                    invRhocCmu_half = 0.5d0*Cmu*invRhoc
    
                    !---Viscous term
                    viscous_w1 = 1.d0*(mu2*dwdx2 - mu1*dwdx1)/dx1(i)
                    viscous_w2 = 1.d0*(mu4*dwdy4 - mu3*dwdy3)/dx2(j )
                    viscous_w3 = 1.d0*(Mu(i,j ,k )*dwdz6 - Mu(i,j ,km)*dwdz5) / dmx3(k)
                    viscous_w31 = 1.d0*(mu2*dudz2 - mu1*dudz1)/dx1(i)
                    viscous_w32 = 1.d0*(mu4*dvdz4 - mu3*dvdz3)/dx2(j)

                    RHS(i,j,k) =  invRhocCmu_half * (viscous_w1 + viscous_w2 + 2.*viscous_w3 + viscous_w31 + viscous_w32)

                    !---Pressure term
                    RHS(i,j,k) = RHS(i,j,k) &
                                - Cmp*invRhoc*(P(i,j,k) - P(i,j,km))/dmx3(k)

                    !---wbc
                    zbc_down = 0.25d0*v3/dmx2(j)*WBCbt_sub(i,k) &
                             - 0.25d0*dwdy3*(dx3(km)*VBCbt_sub(i,k) + dx3(k)*VBCbt_sub(i,km))/dmx3(k)*0.5d0     

                    zbc_up = - 0.25d0*v4/dmx2(jp)*WBCup_sub(i,k) &
                             - 0.25d0*dwdy4*(dx3(km)*VBCup_sub(i,k) + dx3(k)*VBCup_sub(i,km))/dmx3(k)*0.5d0

                    zbc_down = zbc_down + invRhocCmu_half/dx2(j)*mu3/dmx2(jp)*WBCbt_sub(i,k)
                    zbc_up   = zbc_up   + invRhocCmu_half/dx2(j)*mu4/dmx2(j )*WBCup_sub(i,k)

                    zbc_down = zbc_down - invRhocCmu_half/dx2(j)*mu3/dmx3(k)*(VBCbt_sub(i,k) - VBCbt_sub(i,km))
                    zbc_up   = zbc_up   + invRhocCmu_half/dx2(j)*mu4/dmx3(k)*(VBCup_sub(i,k) - VBCup_sub(i,km))

                    RHS(i,j,k) = RHS(i,j,k) &
                                + dble(1.d0 - jwm)*zbc_down + dble(1.d0 - jwp)*zbc_up

                    !---M33Wn
                    mACI =  invRhocCmu_half/dx1(i)*(mu2/dmx1(ip) + mu1/dmx1(i)) &
                            + 0.25d0*(-u2/dmx1(ip) + u1/dmx1(i))
                    mAPI = -invRhocCmu_half/dx1(i)*mu2/dmx1(ip) &
                            + 0.25d0*( u2/dmx1(ip))
                    mAMI = -invRhocCmu_half/dx1(i)*mu1/dmx1(i) &
                            + 0.25d0*(-u1/dmx1(i))

                    mACJ =  invRhocCmu_half/dx2(j)*(mu4/dmx2(jp) + mu3/dmx2(j)) &
                            + 0.25d0*(-v4/dmx2(jp) + v3/dmx2(j ))
                    mAPJ = -invRhocCmu_half/dx2(j)*mu4/dmx2(jp) &
                            + 0.25d0*( v4/dmx2(jp))
                    mAPJ = mAPJ*dble(jwp)
                    mAMJ = -invRhocCmu_half/dx2(j)*mu3/dmx2(j) &
                            + 0.25d0*(-v3/dmx2(j))
                    mAMJ = mAMJ*dble(jwm)
                    
                    mACK =  invRhocCmu_half/dmx3(k)*(Mu(i,j ,k )/dx3(k) + Mu(i,j ,km)/dx3(km))*2.d0 &
                                + 0.25d0*(-w6/dx3(k) + w5/dx3(km) + 0.5d0*dwdz6 + 0.5d0*dwdz5)
                    mAPK = -invRhocCmu_half/dmx3(k)*Mu(i,j ,k )/dx3(k)*2.d0 &
                                + 0.25d0*( w6/dx3(k) + 0.5d0*dwdz6)
                    mAMK = -invRhocCmu_half/dmx3(k)*Mu(i,j ,km)/dx3(km)*2.d0 &
                                + 0.25d0*(-w5/dx3(km) + 0.5d0*dwdz5)

                    RHS(i,j,k) = RHS(i,j,k) &
                            -  ( mAPI*W(ip,j, k ) + mACI*W(i,j,k) + mAMI*W(im,j, k ) &
                                +mAPJ*W(i, jp,k ) + mACJ*W(i,j,k) + mAMJ*W(i, jm,k ) &
                                +mAPK*W(i, j, kp) + mACK*W(i,j,k) + mAMK*W(i, j, km) )
                            
                    RHS(i,j,k) = RHS(i,j,k) &
                                - (0.25d0*(u2*dwdx2 + u1*dwdx1) & 
                                -invRhocCmu_half/dx1(i)*(mu2*dudz2 - mu1*dudz1) )

                    RHS(i,j,k) = RHS(i,j,k) &
                                -(0.25d0*(v3*dble(jwm)*dwdy3 + v4*dble(jwp)*dwdy4) &
                                -invRhocCmu_half/dx2(j)*(mu4*dvdz4*dble(jwp) - mu3*dvdz3*dble(jwm)) )

                    RHS(i,j,k) = RHS(i,j,k) &
                                -(0.25d0*(ddu2*dwdx2 + ddu1*dwdx1) &
                                -invRhocCmu_half/dx1(i)*(mu2*ddudz2 - mu1*ddudz1) )

                    RHS(i,j,k) = RHS(i,j,k) &
                                -(0.25d0*(ddv3*dble(jwm)*dwdy3 + ddv4*dble(jwp)*dwdy4) &
                                -invRhocCmu_half/dx2(j)*(mu4*ddvdz4*dble(jwp) - mu3*ddvdz3*dble(jwm)) )

                    ! M33WVn: Z-direction
                    ac(i,j,k) = mACK*dt + 1.d0
                    ap(i,j,k) = mAPK*dt
                    am(i,j,k) = mAMK*dt
                    RHS(i,j,k)= RHS(i,j,k)*dt
                enddo
            enddo  
        enddo

        ! 1st stage : solve TDM in z-direction
        call PaScaL_TDMA_plan_many_create(ptdma_plan, (n1msub)*(n2msub), comm_1d_x3%myrank, comm_1d_x3%nprocs, comm_1d_x3%mpi_comm)
        call PaScaL_TDMA_many_solve_cycle(ptdma_plan, am, ac, ap, RHS,(n1msub)*(n2msub),(n3msub))
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x3%nprocs)

        ! Deallocate A-matrix for the 2nd stage
        deallocate(am,ac,ap)

        ! Transpose the array (i,j,k) to array (j,k,i) for the tridagonal systems in x-direction
        allocate(RHSI(1:n2msub,1:n3msub,1:n1msub))
        do k = k_indexS, n3msub
            do j = 1, n2msub
                do i = 1, n1msub
                    RHSI(j,k,i)=RHS(i,j,k)
                enddo
            enddo
        enddo

        ! Deallocate the original RHS array which is not required any more
        deallocate(RHS) 

        ! Re-allocate arrays for A matrix coefficient in the second stage
        allocate(am(1:n2msub,1:n3msub,1:n1msub),ac(1:n2msub,1:n3msub,1:n1msub),ap(1:n2msub,1:n3msub,1:n1msub))

        ! 2nd stage : build A matrix and RHS to solve the tridagonal systems in x-direction
        do k = k_indexS, n3msub
            kp = k+1
            km = k-1
            do j = 1, n2msub
                jp = j + 1
                jm = j - 1
                jwp = jC_BC(jp)
                jwm = jC_BC(jm)

                do i = 1, n1msub

                    ip = i + 1
                    im = i - 1

                    u1 = (dx3(km)*U(i,j,k) + dx3(k)*U(i, j,km))/dmx3(k)*0.5d0
                    u2 = (dx3(km)*U(ip,j,k) + dx3(k)*U(ip,j,km))/dmx3(k)*0.5d0                
                
                    invRhoc = 0.5d0*(dx3(km)*invRho(i,j ,k ) + dx3(k )*invRho(i,j ,km) )/dmx3(k)
                    invRhocCmu_half = 0.5d0*Cmu*invRhoc
   
                    ! M33Wn: X-direction
                    mACI =  invRhocCmu_half/dx1(i)*(mu2/dmx1(ip) + mu1/dmx1(i)) &
                            + 0.25d0*(-u2/dmx1(ip) + u1/dmx1(i))
                    mAPI = -invRhocCmu_half/dx1(i)*mu2/dmx1(ip) &
                            + 0.25d0*( u2/dmx1(ip))
                    mAMI = -invRhocCmu_half/dx1(i)*mu1/dmx1(i) &
                            + 0.25d0*(-u1/dmx1(i))
                            
                    ac(j,k,i)=mACI*dt + 1.d0
                    ap(j,k,i)=mAPI*dt
                    am(j,k,i)=mAMI*dt
                enddo
            enddo
        enddo

        ! 2nd stage : solve TDM in x-direction
        call PaScaL_TDMA_plan_many_create(ptdma_plan, (n2msub)*(n3msub), comm_1d_x1%myrank, comm_1d_x1%nprocs, comm_1d_x1%mpi_comm)
        call PaScaL_TDMA_many_solve_cycle(ptdma_plan, am, ac, ap, RHSI,(n2msub)*(n3msub),(n1msub))
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x1%nprocs)

        ! Deallocate A-matrix for the 3rd stage
        deallocate(am,ac,ap)

        ! Transpose the array (j,k,i) to array (i,k,j) to solve the tridagonal systems in y-direction
        allocate(RHSJ(1:n1msub,1:n3msub,1:n2msub))
        do k = 1, n3msub
            do j = 1, n2msub
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
        do k = k_indexS, n3msub
            do j = 1, n2msub
                jm=j-1
                jp=j+1
                jwp = jC_BC(jp)
                jwm = jC_BC(jm)

                do i = 1, n1msub

                    v3 = (dx3(km)*V(i,j, k) + dx3(k)*V(i,j, km))/dmx3(k)*0.5d0
                    v4 = (dx3(km)*V(i,jp,k) + dx3(k)*V(i,jp,km))/dmx3(k)*0.5d0                
                
                    invRhoc = 0.5d0*(dx3(km)*invRho(i,j ,k ) + dx3(k )*invRho(i,j ,km) )/dmx3(k)
                    invRhocCmu_half = 0.5d0*Cmu*invRhoc

                ! M33Wn: Y-direction
                    mACJ =  invRhocCmu_half/dx2(j)*(mu4/dmx2(jp) + mu3/dmx2(j)) &
                            + 0.25d0*(-v4/dmx2(jp) + v3/dmx2(j ))
                    mAPJ = -invRhocCmu_half/dx2(j)*mu4/dmx2(jp) &
                            + 0.25d0*( v4/dmx2(jp))
                    mAPJ = mAPJ*dble(jwp)
                    mAMJ = -invRhocCmu_half/dx2(j)*mu3/dmx2(j)  &
                            + 0.25d0*(-v3/dmx2(j))
                    mAMJ = mAMJ*dble(jwm)

                    ac(i,k,j)=mACJ*dt + 1.d0
                    ap(i,k,j)=mAPJ*dt
                    am(i,k,j)=mAMJ*dt
                enddo
            enddo
        enddo

        ! 3rd stage : solve TDM in y-direction
        call PaScaL_TDMA_plan_many_create(ptdma_plan, (n3msub)*(n1msub), comm_1d_x2%myrank, comm_1d_x2%nprocs, comm_1d_x2%mpi_comm)
        call PaScaL_TDMA_many_solve(ptdma_plan, am, ac, ap, RHSJ,(n3msub)*(n1msub),(n2msub))
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x2%nprocs)

        ! Deallocate A-matrix for the 3rd stage
        deallocate(am,ac,ap)

        ! Update the velocity increments with the transpose of the solution array (i,k,j)
        allocate(dW(0:n1sub,0:n2sub,0:n3sub))
        do k = k_indexS, n3msub
            do j = 1, n2msub
                do i = 1, n1msub
                    dW(i,j,k)=RHSJ(i,k,j)
                enddo
            enddo
        enddo

        ! Deallocate the RHSJ array which is not required any more
        deallocate(RHSJ)

        ! Nullify grid information pointer
        nullify(dx1, dx2, dx3, dmx1, dmx2, dmx3)

    end subroutine mpi_momentum_solvedW

    !>
    !> @brief       dV update from the intermediate velocity increments
    !> @param       T       Temperature field
    !>
    subroutine mpi_momentum_blockLdV(T)

        use mpi_subdomain,  only : j_indexS, jS_BC
        use mpi_subdomain,  only : n1sub, n2sub, n3sub, n1msub, n2msub, n3msub
        use mpi_subdomain,  only : dx2_sub,dx3_sub,dmx2_sub,dmx3_sub

        implicit none

        double precision, dimension(0:n1sub,0:n2sub,0:n3sub) ::  T

        !> @{ Local pointer for subdomain variables
        double precision, dimension(:), pointer :: dx2, dx3
        double precision, dimension(:), pointer :: dmx2, dmx3
        !> @}

        !> @{ Local indexing variables
        integer :: i,j,k
        integer :: im,jm,km
        integer :: ip,jp,kp
        integer :: jvp,jvm
        !> @}

        !> @{ Local variables for dV update
        double precision :: dwm5,dwm6
        double precision :: dvdz5,dvdz6,ddwdy5,ddwdy6
        double precision :: mua,muc,mub,mu5,mu6
        double precision :: invRhoc
        !> @}

        ! Pointer of grid information
        dx2 => dx2_sub
        dx3 => dx3_sub
        dmx2 => dmx2_sub
        dmx3 => dmx3_sub

        ! Update
        do k = 1, n3msub
            kp = k+1
            km = k-1
            do j =  j_indexS, n2msub
                jp = j + 1
                jm = j - 1
                jvp = jS_BC(jm)
                jvm = jS_BC(jp)

                do i = 1, n1msub
                
                    dwm5 = (dx2(jm)*dW(i,j,k ) + dx2(j)*dW(i,jm,k ))/dmx2(j)*0.5d0
                    dwm6 = (dx2(jm)*dW(i,j,kp) + dx2(j)*dW(i,jm,kp))/dmx2(j)*0.5d0

                    ddwdy5 = (dW(i,j,k ) - dW(i,jm,k ))/dmx2(j)
                    ddwdy6 = (dW(i,j,kp) - dW(i,jm,kp))/dmx2(j)  
                    
                    dvdz5 = (V(i,j,k ) - V(i,j,km))/dmx3(k )
                    dvdz6 = (V(i,j,kp) - V(i,j,k ))/dmx3(kp)

                    invRhoc = 0.5d0*(dx2(jm)*invRho(i,j ,k ) + dx2(j )*invRho(i,jm,k ))/dmx2(j)
                    
                    mua = 0.5d0*(dx2(jm)*Mu(i,j ,km) + dx2(j)*Mu(i,jm,km))/dmx2(j)
                    muc = 0.5d0*(dx2(jm)*Mu(i,j ,k ) + dx2(j)*Mu(i,jm,k ))/dmx2(j)
                    mub = 0.5d0*(dx2(jm)*Mu(i,j ,kp) + dx2(j)*Mu(i,jm,kp))/dmx2(j)

                    mu5 = 0.5d0*(dx3(km)*muc + dx3(k)*mua)/dmx3(k )
                    mu6 = 0.5d0*(dx3(kp)*muc + dx3(k)*mub)/dmx3(kp)

                    invRhoc = 0.5d0*(dx2(jm)*invRho(i,j ,k ) + dx2(j )*invRho(i,jm,k ))/dmx2(j)
                    !>  dV(i,j,k) = dV(i,j,k) - dt*M23dWm        
                    dV(i,j,k) =  dV(i,j,k) &
                                - dt*(0.25d0*(dwm5*dvdz5 + dwm6*dvdz6) &
                                - 0.5d0*Cmu*invRhoc/dx3(k)*(mu6*ddwdy6 - mu5*ddwdy5) )

                enddo
            enddo
        enddo

        ! Nullify grid information pointer
        nullify(dx2, dx3, dmx2, dmx3)

    end subroutine mpi_momentum_blockLdV

    !>
    !> @brief       dU update from the intermediate velocity increments
    !> @param       T       Temperature field
    !>
    subroutine mpi_momentum_blockLdU(T)

        use mpi_subdomain,  only : i_indexS, jC_BC
        use mpi_subdomain,  only : n1sub, n2sub, n3sub, n1msub, n2msub, n3msub
        use mpi_subdomain,  only : dx1_sub,dx2_sub,dx3_sub,dmx1_sub,dmx2_sub,dmx3_sub

        implicit none

        double precision, dimension(0:n1sub,0:n2sub,0:n3sub) ::  T

        !> @{ Local pointer for subdomain variables
        double precision, dimension(:), pointer :: dx1, dx2, dx3
        double precision, dimension(:), pointer :: dmx1, dmx2, dmx3
        !> @}
        
        !> @{ Local indexing variables
        integer :: i,j,k
        integer :: im,jm,km
        integer :: ip,jp,kp
        integer :: jup,jum
        !> @}

        !> @{ Local variables for dV update
        double precision :: dvm3,dvm4
        double precision :: dudy3,dudy4,dvmdx3,dvmdx4
        double precision :: mua,mub,muc,mu3,mu4
        double precision :: invRhoc,viscous_u1

        double precision :: dwm5,dwm6
        double precision :: dudz5,dudz6,dwmdx5,dwmdx6
        double precision :: mu5,mu6
        !> @}

        ! Pointer of grid information
        dx1 => dx1_sub
        dx2 => dx2_sub
        dx3 => dx3_sub
        dmx1 => dmx1_sub
        dmx2 => dmx2_sub
        dmx3 => dmx3_sub
    
        ! Update
        do k = 1, n3msub
            km=k-1
            kp=k+1
            do j = 1, n2msub
                jm = j-1
                jp = j+1
                jum = jC_BC(jm)
                jup = jC_BC(jp)

                do i = 1, n1msub

                    im = i - 1
               
                    invRhoc = 0.5d0*(dx1(im)*invRho(i,j ,k ) + dx1(i)*invRho(im,j ,k ) )/dmx1(i)

                    dvm3 = (dx1(i)*dV(im,j, k) + dx1(im)*dV(i,j, k))/dmx1(i)*0.5d0*dble(jum)
                    dvm4 = (dx1(i)*dV(im,jp,k) + dx1(im)*dV(i,jp,k))/dmx1(i)*0.5d0*dble(jup)
                
                    dudy3 = (U(i,j, k) - U(i,jm,k))/dmx2(j )
                    dudy4 = (U(i,jp,k) - U(i,j, k))/dmx2(jp)
                
                    dvmdx3 = (dV(i,j ,k) - dV(im,j ,k))/dmx1(i)*dble(jum)
                    dvmdx4 = (dV(i,jp,k) - dV(im,jp,k))/dmx1(i)*dble(jup)
                
                    mua = 0.5d0*(dx1(im)*Mu(i,jm,k ) + dx1(i)*Mu(im,jm,k ))/dmx1(i)
                    muc = 0.5d0*(dx1(im)*Mu(i,j ,k ) + dx1(i)*Mu(im,j ,k ))/dmx1(i)
                    mub = 0.5d0*(dx1(im)*Mu(i,jp,k ) + dx1(i)*Mu(im,jp,k ))/dmx1(i)
                    mu3 = 0.5d0*(dx2(jm)*muc + dx2(j)*mua)/dmx2(j )
                    mu4 = 0.5d0*(dx2(jp)*muc + dx2(j)*mub)/dmx2(jp)
                
                    dwm5 = (dx1(i)*dW(im,j,k ) + dx1(im)*dW(i,j,k ))/dmx1(i)*0.5d0 
                    dwm6 = (dx1(i)*dW(im,j,kp) + dx1(im)*dW(i,j,kp))/dmx1(i)*0.5d0 
                
                    dudz5 = (U(i,j,k ) - U(i,j,km))/dmx3(k )
                    dudz6 = (U(i,j,kp) - U(i,j,k ))/dmx3(kp)
                
                    dwmdx5 = (dW(i,j,k ) - dW(im,j,k ))/dmx1(i)
                    dwmdx6 = (dW(i,j,kp) - dW(im,j,kp))/dmx1(i)
                
                    mua = 0.5d0*(dx1(im)*Mu(i,j ,km) + dx1(i)*Mu(im,j ,km))/dmx1(i)
                    muc = 0.5d0*(dx1(im)*Mu(i,j ,k ) + dx1(i)*Mu(im,j ,k ))/dmx1(i)
                    mub = 0.5d0*(dx1(im)*Mu(i,j ,kp) + dx1(i)*Mu(im,j ,kp))/dmx1(i)            
                    mu5 = 0.5d0*(dx3(km)*muc + dx3(k)*mua)/dmx3(k )
                    mu6 = 0.5d0*(dx3(kp)*muc + dx3(k)*mub)/dmx3(kp)
                    
                    !>  dU(i,j,k) = dU(i,j,k) - dt*M12dVm - dt*M13dWm
                    dU(i,j,k) = dU(i,j,k) &
                                -dt*( 0.25d0*(dvm4*dudy4 + dvm3*dudy3) - 0.5d0*Cmu*invRhoc/dx2(j)*(mu4*dvmdx4 - mu3*dvmdx3) ) &
                                -dt*( 0.25d0*(dwm6*dudz6 + dwm5*dudz5) - 0.5d0*Cmu*invRhoc/dx3(k)*(mu6*dwmdx6 - mu5*dwmdx5) )
                enddo
            enddo
        enddo

        ! Nullify grid information pointer
        nullify(dx1, dx2, dx3, dmx1, dmx2, dmx3)

    end subroutine mpi_momentum_blockLdU

end module mpi_momentum