!======================================================================================================================
!> @file        module_solve_pressure.f90
!> @brief       This file contains a module of pressure solver for PaScaL_TCS.
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
!> @brief       Module for pressure solver
!> @details     This module solves (incremental) pressure Poisson equation and update pressure fileds.
!>              It uses 2D FFT in x-z directions to decouple the discretized equation and and 
!>              solve the remaining equation using 1D PaScaL_TDMA.
!>
module mpi_pressure

    use global

    implicit none

    !> @{ Main unknown arrays of incremental pressure and intermediate terms
    double precision, allocatable, dimension(:,:,:) :: dPhat,dP0,dP
    !> @}

    !> @{ RHS coefficient
    double precision, allocatable, dimension(:,:,:) :: PRHS
    !> @}

    !> @{ Local array for wave number
    double precision, allocatable, dimension(:) :: dzk
    double precision, allocatable, dimension(:) :: dxk1, dxk2
    !> @}

    contains    

    !>
    !> @brief       Assign variables for pressure solver
    !>
    subroutine mpi_pressure_allocation()

        use mpi_subdomain,  only : n1sub, n2sub, n3sub 
        
        implicit none

        allocate(dPhat(0:n1sub,0:n2sub,0:n3sub),dP0(0:n1sub,0:n2sub,0:n3sub))

        dPhat(:,:,:) = 0.d0
        dP0(:,:,:) = 0.d0

    end subroutine mpi_pressure_allocation

    !>
    !> @brief       Deallocate variables for pressure solver
    !>
    subroutine mpi_pressure_clean()

        implicit none

        deallocate(dPhat,dP0)
        deallocate(dzk, dxk1, dxk2)

    end subroutine mpi_pressure_clean

    !>
    !> @brief       Initialize intermediate pressure field
    !>
    subroutine mpi_pressure_initial()

        use mpi_subdomain,  only : n1sub, n2sub, n3sub 

        implicit none

        integer :: i,j,k
        
        dPhat(0:n1sub,0:n2sub,0:n3sub)=0.d0
        dP0(0:n1sub,0:n2sub,0:n3sub) = 0.d0

    end subroutine mpi_pressure_initial

    !>
    !> @brief       Calculate wave number for decoupled pressure Poisson equation
    !>
    subroutine mpi_pressure_wave_number()

        use mpi_subdomain,  only : h3p, h1psub_Ksub, h1psub_Ksub_ista, h1pKsub, h1pKsub_ista

        integer ::i, k, im, km
        double precision :: ddx,ddz

        ddx=L1/dble(n1m)
        ddz=L3/dble(n3m)

        allocate(dzk(1:n3m))
        allocate(dxk1(1:h1pKsub))   
        allocate(dxk2(1:h1psub_Ksub))

        ! Wave number in the z-direction
        do k = 1, n3m
            if(k <= h3p)then
                km = k - 1
            else
                km = N3m - k + 1
            end if
            dzk(k) = 2.0d0*(1.0d0-cos(2.0d0*PI*dble(km)*ddz/L3))/(ddz*ddz) 
        enddo 

        ! Wave number in the x-direction for transper scheme 1
        do i = 1, h1pKsub
            im = (h1pKsub_ista+i-1) - 1
            dxk1(i) = 2.0d0*(1.0d0-cos(2.0d0*PI*dble(im)*ddx/L1))/(ddx*ddx)
        enddo

        ! Wave number in the x-direction for transper scheme 2
        do i = 1, h1psub_Ksub
            im = (h1psub_Ksub_ista+i-1) - 1
            dxk2(i) = 2.0d0*(1.0d0-cos(2.0d0*PI*dble(im)*ddx/L1))/(ddx*ddx)
        enddo

    end subroutine mpi_pressure_wave_number

    !>
    !> @brief       Calculate wave number for decoupled pressure Poisson equation
    !> @param       invRho  I   nverse of density field
    !> @param       U           Velocity field in x-direction
    !> @param       V           Velocity field in y-direction
    !> @param       W           Velocity field in z-direction
    !> @param       VBCup_sub   y-velocity boundary condition at upper wall
    !> @param       VBCbt_sub   y-velocity boundary condition at lower wall
    !>
    subroutine mpi_pressure_RHS(invRho,U,V,W,VBCup_sub,VBCbt_sub)

        use mpi_subdomain,  only : n1sub, n2sub, n3sub 
        use mpi_subdomain,  only : i_indexS, jC_BC
        use mpi_subdomain,  only : n1sub, n2sub, n3sub, n1msub, n2msub, n3msub
        use mpi_subdomain,  only : dx1_sub,dx2_sub,dx3_sub,dmx1_sub,dmx2_sub,dmx3_sub

        implicit none

        double precision, dimension(0:n1sub,0:n2sub,0:n3sub) :: invRho,U,V,W
        double precision, dimension(0:n1sub,0:n3sub) :: VBCup_sub,VBCbt_sub

        !> @{ Local pointer for subdomain variables
        double precision, dimension(:), pointer :: dx1, dx2, dx3
        double precision, dimension(:), pointer :: dmx1, dmx2, dmx3
        !> @}

        !> @{ Local indexing variables
        integer :: i,j,k
        integer :: kp,km  
        integer :: jp,jm,jvm,jvp
        integer :: ip,im
        !> @}

        !> @{ Local variables for matrix and RHS coefficients calculation
        double precision :: DivUm,cbc
        double precision :: ddpdx1,ddpdx2,ddpdy3,ddpdy4,ddpdz5,ddpdz6
        double precision :: invRho1,invRho2,invRho3,invRho4,invRho5,invRho6
        double precision :: ExtraTerm
        !> @}

        ! Pointer of grid information
        dx1 => dx1_sub
        dx2 => dx2_sub
        dx3 => dx3_sub
        dmx1 => dmx1_sub
        dmx2 => dmx2_sub
        dmx3 => dmx3_sub
    
        ! Allocate arrays for RHS coefficients
        allocate(PRHS(1:n1msub,1:n2msub,1:n3msub))
        
        do k = 1, n3msub

            kp = k+1
            km = k-1

            do j = 1, n2msub

                jp = j + 1
                jm = j - 1
                jvm = jC_BC(jm)
                jvp = jC_BC(jp)

                do i = i_indexS, n1msub

                    ip = i + 1
                    im = i - 1
                    DivUm = (            U(ip,j, k) -           U(i,j,k)  )/dx1(i) &
                          + ( dble(jvp) *V(i, jp,k) - dble(jvm)*V(i,j,k)  )/dx2(j) &
                          + (            W(i, j, kp) -          W(i,j,k)  )/dx3(k)

                    cbc = dble(1.d0 - jvm)*VBCbt_sub(i,k)/dx2(j) &
                        - dble(1.d0 - jvp)*VBCup_sub(i,k)/dx2(j)

                    ddpdx1 = (dPhat(i, j, k) - dPhat(im,j, k))/dmx1(i)
                    ddpdx2 = (dPhat(ip,j, k) - dPhat(i, j, k))/dmx1(ip)

                    ddpdy3 = (dPhat(i,j , k) - dPhat(i,jm, k))/dmx2(j )
                    ddpdy4 = (dPhat(i,jp, k) - dPhat(i,j , k))/dmx2(jp)
                    
                    ddpdz5 = (dPhat(i,j , k) - dPhat(i,j, km))/dmx3(k )
                    ddpdz6 = (dPhat(i,j, kp) - dPhat(i,j , k))/dmx3(kp)

                    invRho1 = 0.5d0/dmx1(i) *(dx1(im)*invRho(i,j,k ) + dx1(i)*invRho(im,j ,k ))
                    invRho2 = 0.5d0/dmx1(ip)*(dx1(ip)*invRho(i,j,k ) + dx1(i)*invRho(ip,j ,k ))
                    invRho3 = 0.5d0/dmx2(j) *(dx2(jm)*invRho(i,j,k ) + dx2(j)*invRho(i, jm,k ))
                    invRho4 = 0.5d0/dmx2(jp)*(dx2(jp)*invRho(i,j,k ) + dx2(j)*invRho(i, jp,k ))
                    invRho5 = 0.5d0/dmx3(k) *(dx3(km)*invRho(i,j,k ) + dx3(K)*invRho(i, j ,km))
                    invRho6 = 0.5d0/dmx3(kp)*(dx3(kp)*invRho(i,j,k ) + dx3(K)*invRho(i, j ,kp))

                    ExtraTerm = ((1.d0-invRho2)*ddpdx2 - (1.d0-invRho1)*ddpdx1)/dx1(i) &
                              + ((1.d0-invRho4)*ddpdy4 - (1.d0-invRho3)*ddpdy3)/dx2(j) &
                              + ((1.d0-invRho6)*ddpdz6 - (1.d0-invRho5)*ddpdz5)/dx3(k) 

                    PRHS(i,j,k) = (DivUm - cbc)/dt/Cmp + ExtraTerm
                enddo
                
            end do
        end do
        
    end subroutine mpi_pressure_RHS

    !>
    !> @brief       Main pressure solver with transpose scheme 1
    !> @param       dx2         Grid spacing in y-direction for the example problem
    !> @param       dmx2        Mesh length in y-direction for the example problem
    !>
    subroutine mpi_pressure_Poisson_FFT1(dx2,dmx2)

        use MPI
        use PaScaL_TDMA
        use mpi_topology,   only : comm_1d_x1, comm_1d_x2, comm_1d_x3, comm_1d_x1n2, myrank
        use mpi_subdomain,  only : n1sub, n2sub, n3sub, n1msub, n2msub, n3msub
        use mpi_subdomain,  only : h1p, h3p
        use mpi_subdomain,  only : h1pKsub, n2mIsub, n2mKsub, h1pKsub_ista, h1pKsub_iend, n2mKsub_jsta, n2mKsub_jend
        use mpi_subdomain,  only : ddtype_dble_C_in_C2I, ddtype_dble_I_in_C2I
        use mpi_subdomain,  only : ddtype_cplx_I_in_I2K, ddtype_cplx_K_in_I2K
        use mpi_subdomain,  only : countsendI, countdistI, countsendK, countdistK

        implicit none

        include 'fftw3.f'

        double precision :: dx2(0:n2sub),dmx2(0:n2sub)

        type(ptdma_plan_many)     :: ptdma_plan

        integer :: i,j,k,ierr,jp        !< Local index variables 

        integer(8) :: plan1             !< Plan for FFTW

        !> @{ Local variables for matrix coefficient calculation
        double precision :: fft_amj,fft_apj,fft_acj 
        !> @}

        !> @{ Local arrays for matrix coefficients
        double precision, allocatable, dimension(:,:,:) :: Am_r,Ac_r,Ap_r,Be_r
        double precision, allocatable, dimension(:,:,:) :: Am_c,Ac_c,Ap_c,Be_c
        !> @}

        !> @{ Local variables calculation
        double precision :: Pbc_a,Pbc_b
        !> @}

        !> @{ Local variables to calculate average of solution
        double precision :: AVERsub, AVERmpi
        !> @}

        integer     :: buffer_cd_size, buffer_dp_size       !< Local variables for ommunication buffer size

        !> @{ Local variables for communication buffers
        double precision, allocatable, target, dimension(:) :: buffer_dp1, buffer_dp2
        complex(8), allocatable, target, dimension(:)       :: buffer_cd1, buffer_cd2
        !> @}

        !> @{ Local variables for RHS transpose
        double precision, pointer, dimension(:,:,:)         :: RHS_Iline
        double precision, pointer, dimension(:,:,:)         :: tmp
        complex(8), pointer, dimension(:,:,:)               :: RHSIhat_Iline
        complex(8), pointer, dimension(:,:,:)               :: RHSIhat_Kline
        complex(8), pointer, dimension(:,:,:)               :: RHSIKhat_Kline
        !> @}

        buffer_dp_size = max( n1m * n2mIsub * n3msub, n1msub*n2msub*n3msub )
        buffer_cd_size = max( h1pKsub * n2mKsub * n3m, h1p * n2mIsub * n3msub )

        allocate( buffer_dp1( buffer_dp_size ) )
        allocate( buffer_dp2( buffer_dp_size ) )
        allocate( buffer_cd1( buffer_cd_size ) )
        allocate( buffer_cd2( buffer_cd_size ) )

        ! Alltoall C to I
        RHS_Iline(1:n1m,1:n2mIsub,1:n3msub) => buffer_dp1
        call MPI_ALLTOALLW(PRHS,      countsendI, countdistI, ddtype_dble_C_in_C2I &
                          ,RHS_Iline, countsendI, countdistI, ddtype_dble_I_in_C2I &
                          ,comm_1d_x1%mpi_comm, ierr )
        deallocate(PRHS)

        ! Forward FFT in x-direction
        RHSIhat_Iline(1:h1p,1:n2mIsub,1:n3msub) => buffer_cd1
        call dfftw_plan_many_dft_r2c(plan1, 1, n1m, n2mIsub*n3msub, RHS_Iline, n1m, 1, n1m, RHSIhat_Iline, n1m, 1, h1p, FFTW_ESTIMATE)
        call dfftw_execute_dft_r2c(plan1,RHS_Iline, RHSIhat_Iline)
        call dfftw_destroy_plan(plan1)
        nullify(RHS_Iline)
                        
        ! Alltoall I to K
        RHSIhat_Kline( 1:h1pKsub,1:n2mKsub,1:n3m ) => buffer_cd2
        call MPI_ALLTOALLW(RHSIhat_Iline,   countsendK, countdistK, ddtype_cplx_I_in_I2K &
                            ,RHSIhat_Kline, countsendK, countdistK, ddtype_cplx_K_in_I2K &
                            ,comm_1d_x3%mpi_comm, ierr )
        nullify(RHSIhat_Iline)

        ! Forward FFT in z-direction
        RHSIKhat_Kline(1:n3m,1:h1pKsub,1:n2mKsub) => buffer_cd1
        call dfftw_plan_many_dft(plan1, 1, n3m, h1pKsub*n2mKsub, RHSIhat_Kline, n3m, h1pKsub*n2mKsub, 1, RHSIKhat_Kline, n3m, 1, n3m, fftw_forward, FFTW_ESTIMATE)
        call dfftw_execute_dft(plan1, RHSIhat_Kline, RHSIKhat_Kline)
        call dfftw_destroy_plan(plan1)
        nullify(RHSIhat_Kline)

        ! Allocate arrays for matrix coefficients
        allocate(Am_r(1:n3m,1:h1pKsub,n2mKsub),Ac_r(1:n3m,1:h1pKsub,n2mKsub))
        allocate(Ap_r(1:n3m,1:h1pKsub,n2mKsub),Be_r(1:n3m,1:h1pKsub,n2mKsub))
        allocate(Am_c(1:n3m,1:h1pKsub,n2mKsub),Ac_c(1:n3m,1:h1pKsub,n2mKsub))
        allocate(Ap_c(1:n3m,1:h1pKsub,n2mKsub),Be_c(1:n3m,1:h1pKsub,n2mKsub))

        ! Build matrix coefficients for direct TDM solver in y-direction
        ! Calculate real and imaginary coefficients seperately.
        do j = 1, n2mKsub     
            jp = j + 1
            fft_amj = 1.0d0/dx2(j+n2mKsub_jsta-1)/dmx2(j +n2mKsub_jsta-1)
            fft_apj = 1.0d0/dx2(j+n2mKsub_jsta-1)/dmx2(jp+n2mKsub_jsta-1)

            ! Lower boundary in y-direction
            if(comm_1d_x1n2%myrank==0.and.j==1  )  fft_amj = 0.d0
            ! Upper boundary in y-direction
            if(comm_1d_x1n2%myrank==comm_1d_x1n2%nprocs-1.and.j==n2mKsub )  fft_apj = 0.d0

            fft_acj =  - fft_amj - fft_apj
    
            do i = 1, h1pKsub
                do k = 1, n3m

                    ! Define the RHS for both 'real' and 'complex' 
                    Be_r(k,i,j) = dble(RHSIKhat_Kline(k,i,j))
                    Be_c(k,i,j) = aimag(RHSIKhat_Kline(k,i,j))
        
                    Am_r(k,i,j) = fft_amj
                    Ac_r(k,i,j) = fft_acj - dxk1(i) - dzk(k)
                    Ap_r(k,i,j) = fft_apj
        
                    Am_c(k,i,j) = fft_amj
                    Ac_c(k,i,j) = fft_acj - dxk1(i) - dzk(k)
                    Ap_c(k,i,j) = fft_apj
                enddo
            enddo
        enddo
    
        ! Special treatment for global 1st element at (1,1,1)
        if(comm_1d_x1n2%myrank==0.and.h1pKsub_ista==1.and.comm_1d_x3%myrank==0) then
            Am_r(1,1,1) = 0.d0
            Ac_r(1,1,1) = 1.d0
            Ap_r(1,1,1) = 0.d0
            Be_r(1,1,1) = 0.d0
        endif

        ! Solve TDM in y-direction : Obtain solutions for real and imaginary part seperately.
        call PaScaL_TDMA_plan_many_create(ptdma_plan, n3m*h1pKsub, comm_1d_x1n2%myrank, comm_1d_x1n2%nprocs, comm_1d_x1n2%mpi_comm)
        call PaScaL_TDMA_many_solve(ptdma_plan,Am_r,Ac_r,Ap_r,Be_r,n3m*h1pKsub,n2mKsub)
        call PaScaL_TDMA_many_solve(ptdma_plan,Am_c,Ac_c,Ap_c,Be_c,n3m*h1pKsub,n2mKsub)
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x1n2%nprocs)

        ! Build complex array from the TDM solution
        do j = 1, n2mKsub
            if(h1pKsub_ista==1) then
                Be_c(1,1,j)= 0.d0;
                Be_c(h3p,1,j)= 0.d0;
            else if(h1pKsub_iend==h1p) then
                Be_c(1,h1pKsub,j)= 0.d0;
                Be_c(h3p,h1pKsub,j)= 0.d0;
            endif

            do i = 1, h1pKsub
                do k = 1, n3m
                    RHSIKhat_Kline(k,i,j) = DCMPLX(Be_r(k,i,j), Be_c(k,i,j))
                enddo
            enddo
        end do

        ! Deallocate A-matrix
        deallocate(Am_r,Ac_r,Ap_r,Be_r)
        deallocate(Am_c,Ac_c,Ap_c,Be_c)
                
        ! Backward FFT in z-direction
        RHSIhat_Kline( 1:h1pKsub,1:n2mKsub,1:n3m ) => buffer_cd2
        call dfftw_plan_many_dft(plan1, 1, n3m, h1pKsub*n2mKsub, RHSIKhat_Kline, n3m, 1, n3m, RHSIhat_Kline, n3m, h1pKsub*n2mKsub, 1, fftw_backward, FFTW_ESTIMATE)
        call dfftw_execute_dft(plan1, RHSIKhat_Kline, RHSIhat_Kline)
        call dfftw_destroy_plan(plan1)
        nullify(RHSIKhat_Kline)

        !== Alltoall K to I
        RHSIhat_Iline(1:h1p,1:n2mIsub,1:n3msub) => buffer_cd1
        call MPI_ALLTOALLW(RHSIhat_Kline, countsendK, countdistK, ddtype_cplx_K_in_I2K &
                          ,RHSIhat_Iline, countsendK, countdistK, ddtype_cplx_I_in_I2K &
                          ,comm_1d_x3%mpi_comm, ierr )
        nullify(RHSIhat_Kline)

        ! Backward FFT in x-direction
        RHS_Iline(1:n1m,1:n2mIsub,1:n3msub) => buffer_dp1
        call dfftw_plan_many_dft_c2r(plan1, 1, n1m, n2mIsub*n3msub, RHSIhat_Iline, n1m, 1, h1p, RHS_Iline, n1m, 1, n1m, FFTW_ESTIMATE)
        call dfftw_execute_dft_c2r(plan1, RHSIhat_Iline, RHS_Iline)
        call dfftw_destroy_plan(plan1)
        nullify(RHSIhat_Iline)

        !== Alltoall I to C
        tmp(1:n1msub,1:n2msub,1:n3msub) => buffer_dp2
        call MPI_ALLTOALLW(RHS_Iline, countsendI, countdistI, ddtype_dble_I_in_C2I &
                          ,tmp      , countsendI, countdistI, ddtype_dble_C_in_C2I &
                          ,comm_1d_x1%mpi_comm, ierr )
        nullify(RHS_Iline)

        ! Update the temporary solution
        do k = 1, n3msub
            do j = 1, n2msub
                do i = 1, n1msub
                    tmp(i,j,k) = tmp(i,j,k) / dble(n1m) / dble(n3m)
                enddo
            enddo
        enddo

        ! Calculate thee average of obtained incremental pressure (dP)
        AVERsub=0.d0
        do k=1,n3msub
            AVERsub = AVERsub+sum(tmp(1:n1msub,1:n2msub,k))
        enddo

        AVERmpi=0.d0
        call MPI_ALLREDUCE(AVERsub  , AVERmpi, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
        AVERmpi=AVERmpi/dble(n1m)/dble(n2m)/dble(n3m)

        ! Remove the average
        allocate(dP(0:n1sub,0:n2sub,0:n3sub))
        dP(1:n1msub,1:n2msub,1:n3msub)=tmp(1:n1msub,1:n2msub,1:n3msub)- AVERmpi
        nullify(tmp)
        
        ! Boundary condition treatment in y-direction for the example problem
        ! Lower wall
        if(comm_1d_x2%myrank == 0) then
            Pbc_a = (dmx2(1)+dmx2(2))**2.0d0/((dmx2(1)+dmx2(2))**2.0d0-dmx2(1)**2.0d0)
            Pbc_b = dmx2(1)**2.0d0          /((dmx2(1)+dmx2(2))**2.0d0-dmx2(1)**2.0d0)
            
            do k=1,n3sub
                dP(1:n1msub,0, k) = Pbc_a*dP(1:n1msub,1,   k) - Pbc_b*dP(1:n1msub,2,   k)
            enddo
        endif
        ! Upper wall
        if(comm_1d_x2%myrank==comm_1d_x2%nprocs-1) then
            Pbc_a = (dmx2(n2sub)+dmx2(n2msub))**2.0d0/((dmx2(n2sub)+dmx2(n2msub))**2.0d0-dmx2(n2sub)**2.0d0)
            Pbc_b = dmx2(n2sub)**2.0d0               /((dmx2(n2sub)+dmx2(n2msub))**2.0d0-dmx2(n2sub)**2.0d0)
            
            do k=1,n3sub
                dP(1:n1msub,n2sub,k) = Pbc_a*dP(1:n1msub,n2msub,k) - Pbc_b*dP(1:n1msub,n2sub-2,k)
            enddo
            
        endif

        ! Deallocate communication buffer
        deallocate( buffer_dp1 )
        deallocate( buffer_dp2 )
        deallocate( buffer_cd1 )
        deallocate( buffer_cd2 )

    end subroutine mpi_pressure_Poisson_FFT1

    !>
    !> @brief       Main pressure solver with transpose scheme 2
    !> @param       dx2         Grid spacing in y-direction for the example problem
    !> @param       dmx2        Mesh length in y-direction for the example problem
    !>
    subroutine mpi_pressure_Poisson_FFT2(dx2,dmx2)
        
        use MPI
        use PaScaL_TDMA
        use mpi_topology,   only : comm_1d_x1, comm_1d_x2, comm_1d_x3
        use mpi_subdomain,  only : n1sub, n2sub, n3sub, n1msub, n2msub, n3msub
        use mpi_subdomain,  only : h1p, h3p
        use mpi_subdomain,  only : n3msub_Isub, h1psub, h1psub_Ksub, h1psub_Ksub_ista, h1psub_Ksub_iend
        use mpi_subdomain,  only : ddtype_dble_C_in_C2I, ddtype_dble_I_in_C2I
        use mpi_subdomain,  only : ddtype_cplx_C_in_C2I, ddtype_cplx_I_in_C2I
        use mpi_subdomain,  only : ddtype_cplx_C_in_C2K, ddtype_cplx_K_in_C2K
        use mpi_subdomain,  only : countsendI, countdistI, countsendK, countdistK

        implicit none
        
        include 'fftw3.f'

        double precision :: dx2(0:n2sub),dmx2(0:n2sub)
        type(ptdma_plan_many)     :: ptdma_plan

        integer :: i,j,k,jp,ierr        !< Local index variables 

        integer(8) :: plan1             !< Plan for FFTW
        
        !> @{ Local variables for matrix coefficient calculation
        double precision :: fft_amj,fft_apj,fft_acj
        !> @}

        !> @{ Local arrays for matrix coefficients
        double precision, allocatable, dimension(:,:,:) :: Am_r,Ac_r,Ap_r,Be_r
        double precision, allocatable, dimension(:,:,:) :: Am_c,Ac_c,Ap_c,Be_c
        !> @}

        !> @{ Local variables calculation
        double precision :: Pbc_a,Pbc_b
        !> @}

        !> @{ Local variables to calculate average of solution
        double precision :: AVERsub, AVERmpi
        !> @}

        integer     :: buffer_cd_size, buffer_dp_size       !< Local variables for ommunication buffer size

        !> @{ Local variables for communication buffers
        double precision, allocatable, target, dimension(:) :: buffer_dp1, buffer_dp2
        complex(8), allocatable, target, dimension(:)       :: buffer_cd1, buffer_cd2
        !> @}

        !> @{ Local variables for RHS transpose
        double precision, pointer, dimension(:,:,:)     :: RHS_Iline
        double precision, pointer, dimension(:,:,:)     :: tmp
        complex(8), pointer, dimension(:,:,:)           :: RHSIhat_Iline
        complex(8), pointer, dimension(:,:,:)           :: RHSIhat
        complex(8), pointer, dimension(:,:,:)           :: RHSIhat_Kline
        complex(8), pointer, dimension(:,:,:)           :: RHSIKhat_Kline
        !> @}

        buffer_dp_size = max( n1m * n2msub * n3msub_Isub, n1msub*n2msub*n3msub )
        buffer_cd_size = max( h1p * n2msub * n3msub_Isub, h1psub * n2msub * n3msub, h1psub_Ksub * n2msub * n3m)

        allocate( buffer_dp1( buffer_dp_size ) )
        allocate( buffer_dp2( buffer_dp_size ) )

        allocate( buffer_cd1( buffer_cd_size ) )
        allocate( buffer_cd2( buffer_cd_size ) )

        RHS_Iline(1:n1m,1:n2msub,1:n3msub_Isub) => buffer_dp1

        !== Alltoall C to I
        call MPI_ALLTOALLW(PRHS     , countsendI, countdistI, ddtype_dble_C_in_C2I &
                          ,RHS_Iline, countsendI, countdistI, ddtype_dble_I_in_C2I &
                          ,comm_1d_x1%mpi_comm, ierr )
        deallocate(PRHS)

        ! Forward FFT in x-direction
        RHSIhat_Iline(1:h1p,1:n2msub,1:n3msub_Isub) => buffer_cd1
        call dfftw_plan_many_dft_r2c(plan1, 1, n1m, n2msub*n3msub_Isub, RHS_Iline, n1m, 1, n1m, RHSIhat_Iline, n1m, 1, h1p, FFTW_ESTIMATE)
        call dfftw_execute_dft_r2c(plan1,RHS_Iline, RHSIhat_Iline)
        call dfftw_destroy_plan(plan1)
        nullify(RHS_Iline)

        !== Alltoall I to C
        RHSIhat(1:h1psub,1:n2msub,1:n3msub) => buffer_cd2
        call MPI_ALLTOALLW(RHSIhat_Iline, countsendI, countdistI, ddtype_cplx_I_in_C2I &
                          ,RHSIhat      , countsendI, countdistI, ddtype_cplx_C_in_C2I &
                          ,comm_1d_x1%mpi_comm, ierr )
        nullify(RHSIhat_Iline)

        !== Alltoall C to K
        RHSIhat_Kline(1:h1psub_Ksub,1:n2msub,1:n3m) => buffer_cd1
        call MPI_ALLTOALLW(RHSIhat      , countsendK, countdistK, ddtype_cplx_C_in_C2K &
                          ,RHSIhat_Kline, countsendK, countdistK, ddtype_cplx_K_in_C2K &
                          ,comm_1d_x3%mpi_comm, ierr )
        nullify(RHSIhat)

        ! Forward FFT in z-direction
        RHSIKhat_Kline(1:n3m,1:h1psub_Ksub,1:n2msub) => buffer_cd2
        call dfftw_plan_many_dft(plan1, 1, n3m, h1psub_Ksub*n2msub, RHSIhat_Kline, n3m, h1psub_Ksub*n2msub, 1, RHSIKhat_Kline, n3m, 1, n3m, fftw_forward, FFTW_ESTIMATE)
        call dfftw_execute_dft(plan1, RHSIhat_Kline, RHSIKhat_Kline)
        call dfftw_destroy_plan(plan1)
        nullify(RHSIhat_Kline)

        ! Allocate arrays for matrix coefficients
        allocate(Am_r(1:n3m,1:h1psub_Ksub,n2msub),Ac_r(1:n3m,1:h1psub_Ksub,n2msub))
        allocate(Ap_r(1:n3m,1:h1psub_Ksub,n2msub),Be_r(1:n3m,1:h1psub_Ksub,n2msub))
        allocate(Am_c(1:n3m,1:h1psub_Ksub,n2msub),Ac_c(1:n3m,1:h1psub_Ksub,n2msub))
        allocate(Ap_c(1:n3m,1:h1psub_Ksub,n2msub),Be_c(1:n3m,1:h1psub_Ksub,n2msub))

        ! Build matrix coefficients for direct TDM solver in y-direction
        ! Calculate real and imaginary coefficients seperately.
        do j = 1, n2msub     
            jp = j + 1   
            fft_amj = 1.d0/dx2(j)/dmx2(j )
            fft_apj = 1.d0/dx2(j)/dmx2(jp)

            ! Lower boundary in y-direction
            if(comm_1d_x2%myrank==0.and.j==1  )  fft_amj = 0.d0
            ! Upper boundary in y-direction
            if(comm_1d_x2%myrank==comm_1d_x2%nprocs-1.and.j==n2msub )  fft_apj = 0.d0

            fft_acj =  - fft_amj - fft_apj
    
            do i = 1, h1psub_Ksub
                    
                ! Define the RHS for both 'real' and 'complex' 
                Be_r(1:n3m,i,j) = dble(RHSIKhat_Kline(1:n3m,i,j))
                Be_c(1:n3m,i,j) = aimag(RHSIKhat_Kline(1:n3m,i,j))
    
                Am_r(1:n3m,i,j) = fft_amj
                Ac_r(1:n3m,i,j) = fft_acj - dxk2(i) - dzk(1:n3m)
                Ap_r(1:n3m,i,j) = fft_apj
    
                Am_c(1:n3m,i,j) = fft_amj
                Ac_c(1:n3m,i,j) = fft_acj - dxk2(i) - dzk(1:n3m)
                Ap_c(1:n3m,i,j) = fft_apj
                
            end do
        end do

        ! Special treatment for global 1st element at (1,1,1)
        if(comm_1d_x1%myrank==0.and.comm_1d_x2%myrank==0.and.comm_1d_x3%myrank==0) then
            Am_r(1,1,1) = 0.d0
            Ac_r(1,1,1) = 1.d0
            Ap_r(1,1,1) = 0.d0
            Be_r(1,1,1) = 0.d0
        endif

        ! Solve TDM in y-direction : Obtain solutions for real and imaginary part seperately.
        call PaScaL_TDMA_plan_many_create(ptdma_plan, n3m*h1psub_Ksub, comm_1d_x2%myrank, comm_1d_x2%nprocs, comm_1d_x2%mpi_comm)
        call PaScaL_TDMA_many_solve(ptdma_plan,Am_r,Ac_r,Ap_r,Be_r,n3m*h1psub_Ksub,n2msub)
        call PaScaL_TDMA_many_solve(ptdma_plan,Am_c,Ac_c,Ap_c,Be_c,n3m*h1psub_Ksub,n2msub)
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x2%nprocs)

        ! Build complex array from the TDM solution
        do j = 1, n2msub   
            if(h1psub_Ksub_ista==1) then
                Be_c(1,1,j)= 0.d0;
                Be_c(h3p,1,j)= 0.d0;
            else if(h1psub_Ksub_iend==h1p) then
                Be_c(1,h1psub_Ksub,j)= 0.d0;
                Be_c(h3p,h1psub_Ksub,j)= 0.d0;
            endif
            RHSIKhat_Kline(1:n3m,1:h1psub_Ksub,j) = DCMPLX(Be_r(1:n3m,1:h1psub_Ksub,j), Be_c(1:n3m,1:h1psub_Ksub,j))
        end do

        ! Deallocate A-matrix
        deallocate(Am_r,Ac_r,Ap_r,Be_r)
        deallocate(Am_c,Ac_c,Ap_c,Be_c)

        ! Backward FFT in z-direction
        RHSIhat_Kline(1:h1psub_Ksub,1:n2msub,1:n3m) => buffer_cd1
        call dfftw_plan_many_dft(plan1, 1, n3m, h1psub_Ksub*n2msub, RHSIKhat_Kline, n3m, 1, n3m, RHSIhat_Kline, n3m, h1psub_Ksub*n2msub, 1, fftw_backward, FFTW_ESTIMATE)
        call dfftw_execute_dft(plan1, RHSIKhat_Kline, RHSIhat_Kline)
        call dfftw_destroy_plan(plan1)
        nullify(RHSIKhat_Kline)

        !== Alltoall K to C
        RHSIhat(1:h1psub,1:n2msub,1:n3msub) => buffer_cd2
        call MPI_ALLTOALLW(RHSIhat_Kline, countsendK, countdistK, ddtype_cplx_K_in_C2K &
                          ,RHSIhat      , countsendK, countdistK, ddtype_cplx_C_in_C2K &
                          ,comm_1d_x3%mpi_comm, ierr )
        nullify(RHSIhat_Kline)

        !== Alltoall C to I
        RHSIhat_Iline(1:h1p,1:n2msub,1:n3msub_Isub) => buffer_cd1
        call MPI_ALLTOALLW(RHSIhat      , countsendI, countdistI, ddtype_cplx_C_in_C2I &
                          ,RHSIhat_Iline, countsendI, countdistI, ddtype_cplx_I_in_C2I &
                          ,comm_1d_x1%mpi_comm, ierr )
        nullify(RHSIhat)

        ! Backward FFT in x-direction
        RHS_Iline(1:n1m,1:n2msub,1:n3msub_Isub) => buffer_dp1
        call dfftw_plan_many_dft_c2r(plan1, 1, n1m, n2msub*n3msub_Isub, RHSIhat_Iline, n1m, 1, h1p, RHS_Iline, n1m, 1, n1m, FFTW_ESTIMATE)
        call dfftw_execute_dft_c2r(plan1, RHSIhat_Iline, RHS_Iline)
        call dfftw_destroy_plan(plan1)
        nullify(RHSIhat_Iline)

        !== Alltoall I to C
        tmp(1:n1msub,1:n2msub,1:n3msub) => buffer_dp2
        call MPI_ALLTOALLW(RHS_Iline, countsendI, countdistI, ddtype_dble_I_in_C2I &
                          ,tmp      , countsendI, countdistI, ddtype_dble_C_in_C2I &
                          ,comm_1d_x1%mpi_comm, ierr )
        nullify(RHS_Iline)

        ! Update the temporary solution
        do k = 1, n3msub
            do j = 1, n2msub
                do i = 1, n1msub
                    tmp(i,j,k) = tmp(i,j,k) / dble(n1m) / dble(n3m)
                enddo
            enddo
        enddo
        
        ! Calculate thee average of obtained incremental pressure (dP)
        AVERsub=0.d0
        do k=1,n3msub
            AVERsub = AVERsub+sum(tmp(1:n1msub,1:n2msub,k))
        enddo

        AVERmpi=0.d0
        call MPI_ALLREDUCE(AVERsub  , AVERmpi, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierr)
        AVERmpi=AVERmpi/dble(n1m)/dble(n2m)/dble(n3m)

        ! Remove the average
        allocate(dP(0:n1sub,0:n2sub,0:n3sub))
        dP(1:n1msub,1:n2msub,1:n3msub)=tmp(1:n1msub,1:n2msub,1:n3msub)- AVERmpi
        nullify(tmp)
        
        ! Boundary condition treatment in y-direction for the example problem
        ! Lower wall
        if(comm_1d_x2%myrank == 0) then
            Pbc_a = (dmx2(1)+dmx2(2))**2.d0/((dmx2(1)+dmx2(2))**2.d0-dmx2(1)**2.d0)
            Pbc_b = dmx2(1)**2.d0          /((dmx2(1)+dmx2(2))**2.d0-dmx2(1)**2.d0)
            
            do k=1,n3sub
                dP(1:n1msub,0, k) = Pbc_a*dP(1:n1msub,1, k) - Pbc_b*dP(1:n1msub,2, k)
            enddo
        endif
        ! Upper wall
        if(comm_1d_x2%myrank==comm_1d_x2%nprocs-1) then
            Pbc_a = (dmx2(n2sub)+dmx2(n2msub))**2.d0/((dmx2(n2sub)+dmx2(n2msub))**2.d0-dmx2(n2sub)**2.d0)
            Pbc_b = dmx2(n2sub)**2.d0               /((dmx2(n2sub)+dmx2(n2msub))**2.d0-dmx2(n2sub)**2.d0)
            
            do k=1,n3sub
                dP(1:n1msub,n2sub,k) = Pbc_a*dP(1:n1msub,n2msub,k) - Pbc_b*dP(1:n1msub,n2sub-2,k)
            enddo
            
        endif

        ! Deallocate communication buffer
        deallocate( buffer_dp1 )
        deallocate( buffer_dp2 )
        deallocate( buffer_cd1 )
        deallocate( buffer_cd2 )

    end subroutine mpi_pressure_Poisson_FFT2
    
    !>
    !> @brief       Update the projected incremental pressure
    !> @param       invRho  I   nverse of density field
    !> @param       U           Velocity field in x-direction
    !> @param       V           Velocity field in y-direction
    !> @param       W           Velocity field in z-direction
    !> @param       P           Pressure field
    !> @param       Jrank       MPI rank in y-direction
    !> @param       Jnprocs     MPI size in y-direction
    !>
    subroutine mpi_pressure_Projection(invRho,U,V,W,P,Jrank,Jnprocs)

        use mpi_subdomain,  only : j_indexS
        use mpi_subdomain,  only : n1sub, n2sub, n3sub, n1msub, n2msub, n3msub
        use mpi_subdomain,  only : dx1_sub,dx2_sub,dx3_sub,dmx1_sub,dmx2_sub,dmx3_sub

        implicit none

        double precision, dimension(0:n1sub,0:n2sub,0:n3sub) :: invRho,U,V,W,P
        integer :: Jrank,Jnprocs

        !> @{ Local pointer for subdomain variables
        double precision, dimension(:), pointer :: dx1, dx2, dx3
        double precision, dimension(:), pointer :: dmx1, dmx2, dmx3
        !> @}

        !> @{ Local indexing variables
        integer :: i,j,k
        integer :: kp,km  
        integer :: jp,jm,jvm,jvp
        integer :: im
        !> @}

        !> @{ Local variables for calculation
        double precision :: invRhoc        
        double precision :: Pbc_a,Pbc_b
        !> @}

        ! Pointer of grid information
        dx1 => dx1_sub
        dx2 => dx2_sub
        dx3 => dx3_sub
        dmx1 => dmx1_sub
        dmx2 => dmx2_sub
        dmx3 => dmx3_sub
    
        ! Update pressure and Keep dP0 and dPhat for NOB treatment
        do k = 1, n3msub
            kp = k+1
            km = k-1

            do j = 2, n2msub
                jp = j + 1
                jm = j - 1

                do i = 1, n1msub
                    im = i - 1

                    invRhoc = 0.5d0/dmx1(i)*(dx1(i)*invRho(im,j ,k ) + dx1(im)*invRho(i,j ,k ))
                    U(i,j ,k ) = U(i,j ,k ) &
                                - dt*Cmp*( (dP(i,j ,k ) - dP(im,j , k))/dmx1(i) &
                                +(invRhoc - 1.d0)*(dPhat(i,j ,k )- dPhat(im,j ,k ))/dmx1(i) )
                    
                    invRhoc = 0.5d0/dmx2(j)*(dx2(j)*invRho(i,jm,k ) + dx2(jm)*invRho(i,j ,k ))
                    V(i,j ,k ) = V(i,j ,k ) &
                                - dt*Cmp*( (dP(i,j ,k ) - dP(i,jm,k ))/dmx2(j) &
                                +(invRhoc - 1.d0)*(dPhat(i,j ,k ) - dPhat(i,jm,k ))/dmx2(j) )
                    
                    invRhoc = 0.5d0/dmx3(k)*(dx3(k)*invRho(i,j ,km) + dx3(km)*invRho(i,j ,k ))
                    W(i,j ,k ) = W(i,j ,k ) &
                                - dt*Cmp*( (dP(i,j ,k ) - dP(i, j, km))/dmx3(k) &
                                +(invRhoc - 1.d0)*(dPhat(i,j ,k )- dPhat(i,j ,km))/dmx3(k) )
                    P(i,j,k) = P(i,j,k) + dP(i,j,k)

                    dP0(i,j,k)=(dPhat(i,j,k)+dP0(i,j,k))/2.d0
                    dPhat(i,j,k)=2.d0*dP(i,j,k)-dP0(i,j,k)
                enddo
            enddo

            ! Lower plane treatment
            j=1
            jp = j + 1
            jm = j - 1

            do i = 1, n1msub
                im = i - 1

                invRhoc = 0.5d0/dmx1(i)*(dx1(i)*invRho(im,j ,k ) + dx1(im)*invRho(i,j ,k ))
                U(i,j ,k )  = U(i,j ,k ) &
                            - dt*Cmp*( (dP(i,j ,k ) - dP(im,j , k))/dmx1(i) &
                            +(invRhoc - 1.d0)*(dPhat(i,j ,k )- dPhat(im,j ,k ))/dmx1(i) )

                invRhoc = 0.5d0/dmx2(j)*(dx2(j)*invRho(i,jm,k ) + dx2(jm)*invRho(i,j ,k ))
                V(i,j ,k )  = V(i,j ,k ) &
                            - dt*Cmp*( (dP(i,j ,k ) - dP(i,jm,k ))/dmx2(j) &
                            +(invRhoc - 1.d0)*(dPhat(i,j ,k ) - dPhat(i,jm,k ))/dmx2(j) )*dble(2-j_indexS)
                
                invRhoc = 0.5d0/dmx3(k)*(dx3(k)*invRho(i,j ,km) + dx3(km)*invRho(i,j ,k ))
                W(i,j ,k )  = W(i,j ,k ) &
                            - dt*Cmp*( (dP(i,j ,k ) - dP(i, j, km))/dmx3(k) &
                            +(invRhoc - 1.d0)*(dPhat(i,j ,k )- dPhat(i,j ,km))/dmx3(k) )
                P(i,j,k) = P(i,j,k) + dP(i,j,k)

                dP0(i,j,k)=(dPhat(i,j,k)+dP0(i,j,k))/2.d0
                dPhat(i,j,k)=2.d0*dP(i,j,k)-dP0(i,j,k)
            enddo

            ! Lower boundary treatment
            if(Jrank == 0) then
                Pbc_a = (dmx2(1)+dmx2(2))**2.0d0/((dmx2(1)+dmx2(2))**2.0d0-dmx2(1)**2.0d0)
                Pbc_b = dmx2(1)**2.0d0          /((dmx2(1)+dmx2(2))**2.0d0-dmx2(1)**2.0d0)
                do i = 1, n1sub - 1
                    P(i,0, k) = Pbc_a*P(i,1, k) - Pbc_b*P(i,2, k)
                enddo
            endif

            ! Upper boundary treatment
            if(Jrank==Jnprocs-1) then
                Pbc_a = (dmx2(n2sub)+dmx2(n2msub))**2.0d0/((dmx2(n2sub)+dmx2(n2msub))**2.0d0-dmx2(n2sub)**2.0d0)
                Pbc_b = dmx2(n2sub)**2.0d0               /((dmx2(n2sub)+dmx2(n2msub))**2.0d0-dmx2(n2sub)**2.0d0)
                do i = 1, n1sub - 1
                    P(i,n2sub,k) = Pbc_a*P(i,n2msub,k) - Pbc_b*P(i,n2sub-2,k)
                enddo
            endif
        enddo

        deallocate(dP)

    end subroutine mpi_pressure_Projection

end module mpi_pressure
