!======================================================================================================================
!> @file        module_mpi_subdomain.f90
!> @brief       This file contains a module of subdomains for PaScaL_TCS.
!> @details     The module includes the informations on the partitioned domains and corresponding derived datatype for
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
!> @brief       Module for building subdomains from the physical domain.
!> @details     This module has simulation parameters for subdomains and communication between the subdomains.
!>

module mpi_subdomain

    use mpi
    use global

    implicit none
    
    !> @{ Grid numbers in the subdomain
    integer, public :: n1sub,n2sub,n3sub
    integer, public :: n1msub,n2msub,n3msub
    !> @}
    !> @{ Grid indices of the assigned range
    integer, public :: ista, iend, jsta, jend, ksta, kend
    !> @}

    !> @{ Coordinates of grid points in the subdomain
    double precision, allocatable, dimension(:), public, target :: x1_sub, x2_sub, x3_sub
    !> @}
    !> @{ Grid lengths in the subdomain
    double precision, allocatable, dimension(:), public, target ::  dx1_sub,  dx2_sub,  dx3_sub
    double precision, allocatable, dimension(:), public, target :: dmx1_sub, dmx2_sub, dmx3_sub
    !> @}

    integer, allocatable, dimension(:), public              :: iC_BC, iS_BC, jC_BC, jS_BC, kC_BC, kS_BC
    integer :: i_indexC, i_indexS
    integer :: j_indexC, j_indexS
    integer :: k_indexC, k_indexS
    ! Staggered grid

    !> @{ Derived datatype for communication between x-neighbor subdomains
    integer :: ddtype_sendto_E, ddtype_recvfrom_W, ddtype_sendto_W, ddtype_recvfrom_E
    !> @}
    !> @{ Derived datatype for communication between y-neighbor subdomains
    integer :: ddtype_sendto_N, ddtype_recvfrom_S, ddtype_sendto_S, ddtype_recvfrom_N
    !> @}
    !> @{ Derived datatype for communication between z-neighbor subdomains
    integer :: ddtype_sendto_F, ddtype_recvfrom_B, ddtype_sendto_B, ddtype_recvfrom_F
    !> @}

    !> @{ Half grid numbers plus one in the global domain
    integer :: h1p, h3p
    !> @}
    !> @{ Partitioned grid numbers in the subdomain for transpose scheme 1
    integer :: h1pKsub, n2mIsub, n2mKsub, h1pKsub_ista, h1pKsub_iend, n2mKsub_jsta, n2mKsub_jend
    !> @}

    !> @{ Partitioned grid numbers in the subdomain for transpose scheme 2
    integer :: n3msub_Isub, h1psub, h1psub_Ksub, h1psub_Ksub_ista, h1psub_Ksub_iend
    !> @}

    !> @{ Derived datatype for transpose scheme 1 and transpose scheme 2. C: cubic (org.), I: x-aligned, K: z-aligned.
    integer, allocatable, dimension(:) :: ddtype_dble_C_in_C2I, ddtype_dble_I_in_C2I    ! Scheme 1 & 2
    integer, allocatable, dimension(:) :: ddtype_cplx_C_in_C2I, ddtype_cplx_I_in_C2I    ! Scheme 1 & 2
    integer, allocatable, dimension(:) :: ddtype_cplx_I_in_I2K, ddtype_cplx_K_in_I2K    ! Scheme 1
    integer, allocatable, dimension(:) :: ddtype_cplx_C_in_C2K, ddtype_cplx_K_in_C2K    ! Scheme 2
    !> @}
    !> @{ Array for MPI_Alltoallw
    integer, allocatable, dimension(:) :: countsendI, countdistI
    integer, allocatable, dimension(:) :: countsendK, countdistK
    !> @}

    public  :: subdomain_para_range

    public  :: mpi_subdomain_make
    public  :: mpi_subdomain_clean
    public  :: mpi_subdomain_mesh
    public  :: mpi_subdomain_DDT_ghostcell
    public  :: mpi_subdomain_DDT_transpose1
    public  :: mpi_subdomain_DDT_transpose2
    public  :: mpi_subdomain_ghostcell_update

    contains

    !>
    !> @brief       Prepare the subdomain and determine the size of the subdomain.
    !>
    subroutine mpi_subdomain_make()

        use mpi_topology,   only : comm_1d_x1, comm_1d_x2, comm_1d_x3

        implicit none

        ! Assigning grid numbers and grid indices of my subdomain.
        call subdomain_para_range(1, n1-1, comm_1d_x1%nprocs, comm_1d_x1%myrank, ista, iend)
        n1sub = iend - ista + 2
        call subdomain_para_range(1, n2-1, comm_1d_x2%nprocs, comm_1d_x2%myrank, jsta, jend)
        n2sub = jend - jsta + 2
        call subdomain_para_range(1, n3-1, comm_1d_x3%nprocs, comm_1d_x3%myrank, ksta, kend)
        n3sub = kend - ksta + 2

        n1msub=n1sub-1
        n2msub=n2sub-1
        n3msub=n3sub-1

        ! Allocate subdomain variables.
        allocate( x1_sub(0:n1sub), dmx1_sub(0:n1sub), dx1_sub(0:n1sub))
        allocate( x2_sub(0:n2sub), dmx2_sub(0:n2sub), dx2_sub(0:n2sub))
        allocate( x3_sub(0:n3sub), dmx3_sub(0:n3sub), dx3_sub(0:n3sub))
        
    end subroutine mpi_subdomain_make

    !>
    !> @brief       Deallocate subdomain variables
    !>
    subroutine mpi_subdomain_clean
    
        implicit none
    
        deallocate(x1_sub, dmx1_sub, dx1_sub)
        deallocate(x2_sub, dmx2_sub, dx2_sub)
        deallocate(x3_sub, dmx3_sub, dx3_sub)

        if(allocated(ddtype_dble_C_in_C2I)) deallocate(ddtype_dble_C_in_C2I)
        if(allocated(ddtype_dble_I_in_C2I)) deallocate(ddtype_dble_I_in_C2I)
        if(allocated(ddtype_cplx_C_in_C2I)) deallocate(ddtype_cplx_C_in_C2I)
        if(allocated(ddtype_cplx_I_in_C2I)) deallocate(ddtype_cplx_I_in_C2I)
        if(allocated(ddtype_cplx_C_in_C2K)) deallocate(ddtype_cplx_C_in_C2K)
        if(allocated(ddtype_cplx_K_in_C2K)) deallocate(ddtype_cplx_K_in_C2K)
        if(allocated(ddtype_cplx_I_in_I2K)) deallocate(ddtype_cplx_I_in_I2K)
        if(allocated(ddtype_cplx_K_in_I2K)) deallocate(ddtype_cplx_K_in_I2K)

        if(allocated(countsendI)) deallocate(countsendI)
        if(allocated(countdistI)) deallocate(countdistI)
        if(allocated(countsendK)) deallocate(countsendK)
        if(allocated(countdistK)) deallocate(countdistK)

    end subroutine mpi_subdomain_clean

    !>
    !> @brief       Assign grid coordinates and lengths of subdomains.
    !>
    subroutine mpi_subdomain_mesh()

        use mpi_topology,   only : comm_1d_x1, comm_1d_x2, comm_1d_x3

        implicit none

        integer :: i,j,k
        integer :: request_S2E, request_S2W, ierr
        integer :: status(MPI_STATUS_SIZE)
        
        x1_sub(0:n1sub)=0.d0; dmx1_sub(0:n1sub)=0.d0; dx1_sub(0:n1sub)=0.d0;
        x2_sub(0:n2sub)=0.d0; dmx2_sub(0:n2sub)=0.d0; dx2_sub(0:n2sub)=0.d0;
        x3_sub(0:n3sub)=0.d0; dmx3_sub(0:n3sub)=0.d0; dx3_sub(0:n3sub)=0.d0;
        
        ! X-direction
        ! x_sub is for coordinates.
        if(UNIFORM1 == 1) then 
            do i = ista-1, iend+1
                x1_sub(i-ista+1) = dble(i-1)*L1/dble(N1m) + x1_start
            end do
        else
            do i = ista-1, iend+1
                x1_sub(i-ista+1) = L1*0.5d0*(1.d0 + tanh(0.5d0*GAMMA1*(2.d0*dble(i-1)/dble(N1m)-1.0d0))/tanh(GAMMA1*0.5d0) ) + x1_start
            end do

        end if
        if(pbc1==.false. .and. comm_1d_x1%myrank==0) x1_sub(0)=x1_sub(1)

        ! dx_sub is for the length between two grid points.
        do i = 1, n1msub
            dx1_sub(i) = x1_sub(i+1) - x1_sub(i)
        end do
        ! Communication for grid spacing involved with a ghost grid point
        call MPI_Isend(dx1_sub(n1msub),1, MPI_DOUBLE_PRECISION, comm_1d_x1%east_rank, 111, comm_1d_x1%mpi_comm, request_S2E, ierr)
        call MPI_Irecv(dx1_sub(0),     1, MPI_DOUBLE_PRECISION, comm_1d_x1%west_rank, 111, comm_1d_x1%mpi_comm, request_S2E, ierr)
        call MPI_Wait(request_S2E,STATUS,ierr)
        call MPI_Isend(dx1_sub(1),     1, MPI_DOUBLE_PRECISION, comm_1d_x1%west_rank, 111, comm_1d_x1%mpi_comm, request_S2W, ierr)
        call MPI_Irecv(dx1_sub(n1sub), 1, MPI_DOUBLE_PRECISION, comm_1d_x1%east_rank, 111, comm_1d_x1%mpi_comm, request_S2W, ierr)
        call MPI_Wait(request_S2W,STATUS,ierr)

        if(pbc1==.false. .and. comm_1d_x1%myrank==0) dx1_sub(0)=0.d0
        if(pbc1==.false. .and. comm_1d_x1%myrank==comm_1d_x1%nprocs-1) dx1_sub(n1sub)=0.d0

        ! dmx_sub is for the length of mesh, half of two grid-point spacing.
        do i = 1, n1msub
            dmx1_sub(i) = 0.5d0*(dx1_sub(i-1)+dx1_sub(i))
        end do

        ! Communication for mesh length involved with a ghost grid.
        call MPI_Isend(dmx1_sub(n1msub),1, MPI_DOUBLE_PRECISION, comm_1d_x1%east_rank, 111, comm_1d_x1%mpi_comm, request_S2E, ierr)
        call MPI_Irecv(dmx1_sub(0),     1, MPI_DOUBLE_PRECISION, comm_1d_x1%west_rank, 111, comm_1d_x1%mpi_comm, request_S2E, ierr)
        call MPI_Wait(request_S2E,STATUS,ierr)
        call MPI_Isend(dmx1_sub(1),     1, MPI_DOUBLE_PRECISION, comm_1d_x1%west_rank, 111, comm_1d_x1%mpi_comm, request_S2W, ierr)
        call MPI_Irecv(dmx1_sub(n1sub), 1, MPI_DOUBLE_PRECISION, comm_1d_x1%east_rank, 111, comm_1d_x1%mpi_comm, request_S2W, ierr)
        call MPI_Wait(request_S2W,STATUS,ierr)

        if(pbc1==.false. .and. comm_1d_x1%myrank==comm_1d_x1%nprocs-1) dmx1_sub(n1sub)=0.5d0*(dx1_sub(n1msub)+dx1_sub(n1sub))

        ! Y-direction
        ! y_sub is for coordinates.
        if(UNIFORM2 == 1) then 
            do j = jsta-1, jend+1
                x2_sub(j-jsta+1) = dble(j-1)*L2/dble(N2m) + x2_start
            end do
        else
            do j = jsta-1, jend+1
                x2_sub(j-jsta+1) = L2*0.5*(1. + tanh(0.5*GAMMA2*(2.*dble(j-1)/dble(N2m)-1.0))/tanh(GAMMA2*0.5) ) + x2_start
            end do
        end if
        if(pbc2==.false. .and. comm_1d_x2%myrank==0) x2_sub(0)=x2_sub(1)

        ! dy_sub is for the length between two grid points.
        do j = 1, n2msub
            dx2_sub(j) = x2_sub(j+1) - x2_sub(j)
        end do

        ! Communication for grid spacing involved with a ghost grid point
        call MPI_Isend(dx2_sub(n2msub),1, MPI_DOUBLE_PRECISION, comm_1d_x2%east_rank, 111, comm_1d_x2%mpi_comm, request_S2E, ierr)
        call MPI_Irecv(dx2_sub(0),     1, MPI_DOUBLE_PRECISION, comm_1d_x2%west_rank, 111, comm_1d_x2%mpi_comm, request_S2E, ierr)
        call MPI_Wait(request_S2E,STATUS,ierr)
        call MPI_Isend(dx2_sub(1),     1, MPI_DOUBLE_PRECISION, comm_1d_x2%west_rank, 111, comm_1d_x2%mpi_comm, request_S2W, ierr)
        call MPI_Irecv(dx2_sub(n2sub), 1, MPI_DOUBLE_PRECISION, comm_1d_x2%east_rank, 111, comm_1d_x2%mpi_comm, request_S2W, ierr)
        call MPI_Wait(request_S2W,STATUS,ierr)

        if(pbc2==.false. .and. comm_1d_x2%myrank==0) dx2_sub(0)=0.d0
        if(pbc2==.false. .and. comm_1d_x2%myrank==comm_1d_x2%nprocs-1) dx2_sub(n2sub)=0.d0

        ! dmy_sub is for the length of mesh, half of two grid-point spacing.
        do j = 1, n2msub
            dmx2_sub(j) = 0.5*(dx2_sub(j-1)+dx2_sub(j))
        end do

        ! Communication for mesh length involved with a ghost grid.
        call MPI_Isend(dmx2_sub(n2msub),1, MPI_DOUBLE_PRECISION, comm_1d_x2%east_rank, 111, comm_1d_x2%mpi_comm, request_S2E, ierr)
        call MPI_Irecv(dmx2_sub(0),     1, MPI_DOUBLE_PRECISION, comm_1d_x2%west_rank, 111, comm_1d_x2%mpi_comm, request_S2E, ierr)
        call MPI_Wait(request_S2E,STATUS,ierr)
        call MPI_Isend(dmx2_sub(1),     1, MPI_DOUBLE_PRECISION, comm_1d_x2%west_rank, 111, comm_1d_x2%mpi_comm, request_S2W, ierr)
        call MPI_Irecv(dmx2_sub(n2sub), 1, MPI_DOUBLE_PRECISION, comm_1d_x2%east_rank, 111, comm_1d_x2%mpi_comm, request_S2W, ierr)
        call MPI_Wait(request_S2W,STATUS,ierr)

        if(pbc2==.false. .and. comm_1d_x2%myrank==comm_1d_x2%nprocs-1) dmx2_sub(n2sub)=0.5*(dx2_sub(n2msub)+dx2_sub(n2sub))

        ! Z-direction
        ! z_sub is for coordinates.
        if(UNIFORM3 == 1) then 
            do k = ksta-1, kend+1
                x3_sub(k-ksta+1) = dble(k-1)*L3/dble(N3m) + x3_start
            end do
        else
            do k = ksta-1, kend+1
                x3_sub(k-ksta+1) = L3*0.5*(1. + tanh(0.5*GAMMA3*(2.*dble(k-1)/dble(N3m)-1.0))/tanh(GAMMA3*0.5) ) + x3_start
            end do
        end if
        if(pbc3==.false. .and. comm_1d_x3%myrank==0) x3_sub(0)=x3_sub(1)

        ! dz_sub is for the length between two grid points.
        do k = 1, n3msub
            dx3_sub(k) = x3_sub(k+1) - x3_sub(k)
        end do

        ! Communication for grid spacing involved with a ghost grid point
        call MPI_Isend(dx3_sub(n3msub),1, MPI_DOUBLE_PRECISION, comm_1d_x3%east_rank, 111, comm_1d_x3%mpi_comm, request_S2E, ierr)
        call MPI_Irecv(dx3_sub(0),     1, MPI_DOUBLE_PRECISION, comm_1d_x3%west_rank, 111, comm_1d_x3%mpi_comm, request_S2E, ierr)
        call MPI_Wait(request_S2E,STATUS,ierr)
        call MPI_Isend(dx3_sub(1),     1, MPI_DOUBLE_PRECISION, comm_1d_x3%west_rank, 111, comm_1d_x3%mpi_comm, request_S2W, ierr)
        call MPI_Irecv(dx3_sub(n3sub), 1, MPI_DOUBLE_PRECISION, comm_1d_x3%east_rank, 111, comm_1d_x3%mpi_comm, request_S2W, ierr)
        call MPI_Wait(request_S2W,STATUS,ierr)

        if(pbc3==.false. .and. comm_1d_x3%myrank==0) dx3_sub(0)=0.d0
        if(pbc3==.false. .and. comm_1d_x3%myrank==comm_1d_x3%nprocs-1) dx3_sub(n3sub)=0.d0

        ! dmz_sub is for the length of mesh, half of two grid-point spacing.
        do k = 1, n3msub
            dmx3_sub(k) = 0.5*(dx3_sub(k-1)+dx3_sub(k))
        end do

        ! Communication for mesh length involved with a ghost grid.
        call MPI_Isend(dmx3_sub(n3msub),1, MPI_DOUBLE_PRECISION, comm_1d_x3%east_rank, 111, comm_1d_x3%mpi_comm, request_S2E, ierr)
        call MPI_Irecv(dmx3_sub(0),     1, MPI_DOUBLE_PRECISION, comm_1d_x3%west_rank, 111, comm_1d_x3%mpi_comm, request_S2E, ierr)
        call MPI_Wait(request_S2E,STATUS,ierr)
        call MPI_Isend(dmx3_sub(1),     1, MPI_DOUBLE_PRECISION, comm_1d_x3%west_rank, 111, comm_1d_x3%mpi_comm, request_S2W, ierr)
        call MPI_Irecv(dmx3_sub(n3sub), 1, MPI_DOUBLE_PRECISION, comm_1d_x3%east_rank, 111, comm_1d_x3%mpi_comm, request_S2W, ierr)
        call MPI_Wait(request_S2W,STATUS,ierr)

        if(pbc3==.false. .and. comm_1d_x3%myrank==comm_1d_x3%nprocs-1) dmx3_sub(n3sub)=0.5*(dx3_sub(n3msub)+dx3_sub(n3sub))

    end subroutine mpi_subdomain_mesh

    !>
    !> @brief       Determine whether the previous and next grids are empty(0) or filled(1) in each direction.
    !>
    subroutine mpi_subdomain_indices()

        use mpi_topology,   only : comm_1d_x1, comm_1d_x2, comm_1d_x3

        implicit none

        allocate(iC_BC(0:n1sub),iS_BC(0:n1sub))
        allocate(jC_BC(0:n2sub),jS_BC(0:n2sub))
        allocate(kC_BC(0:n3sub),kS_BC(0:n3sub))

        ! All values are initialized with 1, meaning all grids are not empty and effective.
        i_indexC=1; i_indexS=1
        j_indexC=1; j_indexS=1
        k_indexC=1; k_indexS=1
        
        iC_BC(0:n1sub)=1;iS_BC(0:n1sub)=1;
        jC_BC(0:n2sub)=1;jS_BC(0:n2sub)=1;
        kC_BC(0:n3sub)=1;kS_BC(0:n3sub)=1;

        if(pbc1==.false.) then

            ! Set the first iC_BC and iS_BC to 0, having it empty in case of boundary domain.
            if(comm_1d_x1%myrank==0) then
                i_indexS=2
                iC_BC(0)=0
                iS_BC(0:1)=0
            endif

            ! Set the last iC_BC and iS_BC to 0, having it empty in case of boundary domain.
            if(comm_1d_x1%myrank==comm_1d_x1%nprocs-1) then
                iC_BC(n1sub)=0
                iS_BC(n1sub)=0
            endif

        endif

        if(pbc2==.false.) then

            ! Set the first jC_BC and jS_BC to 0, having it empty in case of boundary domain.
            if(comm_1d_x2%myrank==0) then
                j_indexS=2
                jC_BC(0)=0
                jS_BC(0:1)=0
            endif

            ! Set the last jC_BC and jS_BC to 0, having it empty in case of boundary domain.
            if(comm_1d_x2%myrank==comm_1d_x2%nprocs-1) then
                jC_BC(n2sub)=0
                jS_BC(n2sub)=0
            endif

        endif

        if(pbc3==.false.) then

            ! Set the first kC_BC and kS_BC to 0, having it empty in case of boundary domain.
            if(comm_1d_x3%myrank==0) then
                k_indexS=2
                kC_BC(0)=0
                kS_BC(0:1)=0
            endif

            ! Set the last kC_BC and kS_BC to 0, having it empty in case of boundary domain.
            if(comm_1d_x3%myrank==comm_1d_x3%nprocs-1) then
                kC_BC(n3sub)=0
                kS_BC(n3sub)=0
            endif

        endif
        
    end subroutine mpi_subdomain_indices

    !>
    !> @brief       Deallocate indexing variables
    !>
    subroutine mpi_subdomain_indices_clean()
        implicit none

        deallocate(iC_BC,iS_BC)
        deallocate(jC_BC,jS_BC)
        deallocate(kC_BC,kS_BC)
        
    end subroutine mpi_subdomain_indices_clean

    !>
    !> @brief       Build derived datatypes for subdomain communication using ghostcells.
    !>
    subroutine mpi_subdomain_DDT_ghostcell

        implicit none
        integer :: sizes(0:2), subsizes(0:2), starts(0:2), ierr     ! Local variables for MPI_Type_create_subarray

        ! ddtype sending data to east MPI process (x+ neighbor)
        sizes    = (/n1sub+1,n2sub+1,n3sub+1/)
        subsizes = (/      1,n2sub+1,n3sub+1/)
        starts   = (/n1sub-1,      0,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_sendto_E, ierr)
        call MPI_Type_commit(ddtype_sendto_E,ierr)

        ! ddtype receiving data from west MPI process (x- neighbor)
        starts   = (/      0,      0,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_recvfrom_W, ierr)
        call MPI_Type_commit(ddtype_recvfrom_W,ierr)

        ! ddtype sending data to west MPI process (x- neighbor)
        starts   = (/      1,      0,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_sendto_W, ierr)
        call MPI_Type_commit(ddtype_sendto_W,ierr)

        ! ddtype receiving data from east MPI process (x+ neighbor)
        starts   = (/  n1sub,      0,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_recvfrom_E, ierr)
        call MPI_Type_commit(ddtype_recvfrom_E,ierr)

        ! ddtype sending data to north MPI process (y+ neighbor)
        sizes    = (/n1sub+1,n2sub+1,n3sub+1/)
        subsizes = (/n1sub+1,      1,n3sub+1/)
        starts   = (/      0,n2sub-1,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_sendto_N, ierr)
        call MPI_Type_commit(ddtype_sendto_N,ierr)

        ! ddtype receiving data from south MPI process (y- neighbor)
        starts   = (/      0,      0,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_recvfrom_S, ierr)
        call MPI_Type_commit(ddtype_recvfrom_S,ierr)

        ! ddtype sending data to south MPI process (y- neighbor)
        starts   = (/      0,      1,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_sendto_S, ierr)
        call MPI_Type_commit(ddtype_sendto_S,ierr)

        ! ddtype receiving data from north MPI process (y+ neighbor)
        starts   = (/      0,  n2sub,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_recvfrom_N, ierr)
        call MPI_Type_commit(ddtype_recvfrom_N,ierr)

        ! ddtype sending data to forth MPI process (z+ neighbor)
        sizes    = (/n1sub+1,n2sub+1,n3sub+1/)
        subsizes = (/n1sub+1,n2sub+1,      1/)
        starts   = (/      0,      0,n3sub-1/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_sendto_F, ierr)
        call MPI_Type_commit(ddtype_sendto_F,ierr)

        ! ddtype receiving data from back MPI process (z- neighbor)
        starts   = (/      0,      0,      0/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_recvfrom_B, ierr)
        call MPI_Type_commit(ddtype_recvfrom_B,ierr)

        ! ddtype sending data to back MPI process (z- neighbor)
        starts   = (/      0,      0,      1/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_sendto_B, ierr)
        call MPI_Type_commit(  ddtype_sendto_B,ierr)

        ! ddtype receiving data from forth MPI process (z+ neighbor)
        starts   = (/      0,      0,  n3sub/)
        call MPI_Type_create_subarray(3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                    MPI_DOUBLE_PRECISION, ddtype_recvfrom_F, ierr)
        call MPI_Type_commit(ddtype_recvfrom_F,ierr)

    end subroutine mpi_subdomain_DDT_ghostcell

    !>
    !> @brief       Build derived datatypes for FFT with transpose scheme 1.
    !>
    subroutine mpi_subdomain_DDT_transpose1()

        use mpi_topology,   only : comm_1d_x1, comm_1d_x2, comm_1d_x3
        
        implicit none

        integer :: i,j,k
        integer :: bigsize(3), subsize(3), start(3), ierr 
        integer :: indexA,indexB

        integer, allocatable, dimension(:) :: n1msubAll,n3msubAll,h1pKsubAll,n2mIsubAll
        
        ! C means the partitioned domain in cubic shape as in original decomposition.
        ! I means the partitioned domain in x-aligned shaped for FFT in x-direction
        ! K means the partitioned domain in z-aligned shaped for FFT in z-direction
        allocate(ddtype_dble_C_in_C2I(0:comm_1d_x1%nprocs-1),ddtype_dble_I_in_C2I(0:comm_1d_x1%nprocs-1))
        allocate(ddtype_cplx_I_in_I2K(0:comm_1d_x3%nprocs-1),ddtype_cplx_K_in_I2K(0:comm_1d_x3%nprocs-1))

        allocate(countsendI(0:comm_1d_x1%nprocs-1), countdistI(0:comm_1d_x1%nprocs-1))
        allocate(countsendK(0:comm_1d_x3%nprocs-1), countdistK(0:comm_1d_x3%nprocs-1))

        countsendI(:)=1
        countdistI(:)=0
        countsendK(:)=1
        countdistK(:)=0

        h1p = int(n1m/2)+1
        h3p = int(n3m/2)+1

        call subdomain_para_range(1, h1p, comm_1d_x3%nprocs, comm_1d_x3%myrank, h1pKsub_ista, h1pKsub_iend)
        h1pKsub = h1pKsub_iend - h1pKsub_ista + 1   

        call subdomain_para_range(1, n2msub, comm_1d_x1%nprocs, comm_1d_x1%myrank, indexA, indexB)
        n2mIsub = indexB - indexA + 1   
        n2mKsub = n2mIsub

        allocate(n1msubAll(0:comm_1d_x1%nprocs-1),n3msubAll(0:comm_1d_x3%nprocs-1))
        allocate(h1pKsubAll(0:comm_1d_x3%nprocs-1))
        allocate(n2mIsubAll(0:comm_1d_x1%nprocs-1))

        do i=0,comm_1d_x1%nprocs-1
            call subdomain_para_range(1, n1m, comm_1d_x1%nprocs, i, indexA, indexB)
            n1msubAll(i)= indexB - indexA + 1
            call subdomain_para_range(1, n2msub, comm_1d_x1%nprocs, i, indexA, indexB)
            n2mIsubAll(i)= indexB - indexA + 1
        enddo

        do i=0,comm_1d_x3%nprocs-1      
            call subdomain_para_range(1, n3m, comm_1d_x3%nprocs, i, indexA, indexB)
            n3msubAll(i)= indexB - indexA + 1   
            call subdomain_para_range(1, h1p, comm_1d_x3%nprocs, i, indexA, indexB)
            h1pKsubAll(i)= indexB - indexA + 1                  
        enddo

        call subdomain_para_range(1, n2msub, comm_1d_x1%nprocs, comm_1d_x1%myrank, n2mKsub_jsta, n2mKsub_jend)

        ! DDT for I-C transpose (real type)
        do i=0,comm_1d_x1%nprocs-1
            bigsize(1) = n1msub
            bigsize(2) = n2msub
            bigsize(3) = n3msub
            subsize(1) = n1msub
            subsize(2) = n2mIsubAll(i)
            subsize(3) = n3msub
            start(1) = 0
            start(2) = sum(n2mIsubAll(0:i)) - n2mIsubAll(i)
            start(3) = 0
                        
            call MPI_TYPE_CREATE_SUBARRAY( 3, bigsize, subsize, start, MPI_ORDER_FORTRAN &
                                         , MPI_DOUBLE_PRECISION, ddtype_dble_C_in_C2I(i), ierr )
            call MPI_TYPE_COMMIT(ddtype_dble_C_in_C2I(i),ierr)
                        
            bigsize(1) = n1m
            bigsize(2) = n2mIsub
            bigsize(3) = n3msub
            subsize(1) = n1msubAll(i)
            subsize(2) = n2mIsub
            subsize(3) = n3msub
            start(1) = sum(n1msubAll(0:i)) - n1msubAll(i)
            start(2) = 0
            start(3) = 0
                        
            call MPI_TYPE_CREATE_SUBARRAY( 3, bigsize, subsize, start, MPI_ORDER_FORTRAN &
                                         , MPI_DOUBLE_PRECISION, ddtype_dble_I_in_C2I(i), ierr )
            call MPI_TYPE_COMMIT(ddtype_dble_I_in_C2I(i),ierr)
        enddo

        ! DDT for I-K transpose (complex type)
        do k=0,comm_1d_x3%nprocs-1
            bigsize(1) = h1p
            bigsize(2) = n2mIsub
            bigsize(3) = n3msub
            subsize(1) = h1pKsubAll(k)
            subsize(2) = n2mIsub
            subsize(3) = n3msub
            start(1) = sum(h1pKsubAll(0:k)) - h1pKsubAll(k)
            start(2) = 0
            start(3) = 0
                        
            call MPI_TYPE_CREATE_SUBARRAY( 3, bigsize, subsize, start, MPI_ORDER_FORTRAN &
                                         , MPI_DOUBLE_COMPLEX, ddtype_cplx_I_in_I2K(k), ierr )
            call MPI_TYPE_COMMIT(ddtype_cplx_I_in_I2K(k),ierr)
                        
            bigsize(1) = h1pKsub
            bigsize(2) = n2mIsub
            bigsize(3) = n3m
            subsize(1) = h1pKsub
            subsize(2) = n2mIsub
            subsize(3) = n3msubAll(k)
            start(1) = 0
            start(2) = 0
            start(3) = sum(n3msubAll(0:k)) - n3msubAll(k)
                        
            call MPI_TYPE_CREATE_SUBARRAY( 3, bigsize, subsize, start, MPI_ORDER_FORTRAN &
                                         , MPI_DOUBLE_COMPLEX, ddtype_cplx_K_in_I2K(k), ierr )
            call MPI_TYPE_COMMIT(ddtype_cplx_K_in_I2K(k),ierr)
        enddo

        deallocate(n1msubAll,n3msubAll)
        deallocate(h1pKsubAll)
        deallocate(n2mIsubAll)

    end subroutine mpi_subdomain_DDT_transpose1

    !>
    !> @brief       Build derived datatypes for FFT with transpose scheme 1.
    !>
    subroutine mpi_subdomain_DDT_transpose2()
                                    
        use mpi_topology,   only : comm_1d_x1, comm_1d_x2, comm_1d_x3, myrank

        implicit none

        integer :: i
        integer :: bigsize(3), subsize(3), start(3), ierr 
        integer :: indexA, indexB
        integer, allocatable, dimension(:) :: n3msub_IsubAll,n1msubAll,h1psubAll,h1psub_KsubAll,n3msubAll

        ! C means the partitioned domain in cubic shape as in original decomposition.
        ! I means the partitioned domain in x-aligned shaped for FFT in x-direction
        ! K means the partitioned domain in z-aligned shaped for FFT in z-direction
        allocate(ddtype_dble_C_in_C2I(0:comm_1d_x1%nprocs-1),ddtype_dble_I_in_C2I(0:comm_1d_x1%nprocs-1))
        allocate(ddtype_cplx_C_in_C2I(0:comm_1d_x1%nprocs-1),ddtype_cplx_I_in_C2I(0:comm_1d_x1%nprocs-1))
        allocate(ddtype_cplx_C_in_C2K(0:comm_1d_x3%nprocs-1),ddtype_cplx_K_in_C2K(0:comm_1d_x3%nprocs-1))

        allocate(countsendI(0:comm_1d_x1%nprocs-1), countdistI(0:comm_1d_x1%nprocs-1))
        allocate(countsendK(0:comm_1d_x3%nprocs-1), countdistK(0:comm_1d_x3%nprocs-1))

        countsendI(:)=1
        countdistI(:)=0
        countsendK(:)=1
        countdistK(:)=0

        allocate(n3msub_IsubAll(0:comm_1d_x1%nprocs-1),n1msubAll(0:comm_1d_x1%nprocs-1),h1psubAll(0:comm_1d_x1%nprocs-1))
        allocate(h1psub_KsubAll(0:comm_1d_x3%nprocs-1),n3msubAll(0:comm_1d_x3%nprocs-1))
        
        call subdomain_para_range(1, n3msub, comm_1d_x1%nprocs, comm_1d_x1%myrank, indexA, indexB)
        n3msub_Isub= indexB - indexA + 1     
        
        h1p = int(n1m/2)+1
        h3p = int(n3m/2)+1

        call subdomain_para_range(1, h1p, comm_1d_x1%nprocs, comm_1d_x1%myrank, indexA, indexB)
        h1psub= indexB - indexA + 1

        call subdomain_para_range(indexA, indexB, comm_1d_x3%nprocs, comm_1d_x3%myrank, h1psub_Ksub_ista, h1psub_Ksub_iend)
        h1psub_Ksub= h1psub_Ksub_iend - h1psub_Ksub_ista + 1

        do i=0,comm_1d_x1%nprocs-1
            call subdomain_para_range(1, n3msub, comm_1d_x1%nprocs, i, indexA, indexB)
            n3msub_IsubAll(i)= indexB - indexA + 1        
            call subdomain_para_range(1, n1m, comm_1d_x1%nprocs, i, indexA, indexB)
            n1msubAll(i)= indexB - indexA + 1        
            call subdomain_para_range(1, h1p, comm_1d_x1%nprocs, i, indexA, indexB)
            h1psubAll(i)= indexB - indexA + 1
        enddo
        
        do i=0,comm_1d_x3%nprocs-1             
            call subdomain_para_range(1, h1psub, comm_1d_x3%nprocs, i, indexA, indexB)
            h1psub_KsubAll(i) =indexB - indexA + 1 
            call subdomain_para_range(1, n3-1, comm_1d_x3%nprocs, i, indexA, indexB)
            n3msubAll(i) =indexB - indexA + 1 
        enddo

        ! DDT for I-C transpose (real type)
        do i=0,comm_1d_x1%nprocs-1
            bigsize(1) = n1msub
            bigsize(2) = n2msub
            bigsize(3) = n3msub
            subsize(1) = n1msub
            subsize(2) = n2msub
            subsize(3) = n3msub_IsubAll(i)
            start(1) = 0
            start(2) = 0
            start(3) = sum(n3msub_IsubAll(0:i)) - n3msub_IsubAll(i)
                        
            call MPI_TYPE_CREATE_SUBARRAY( 3, bigsize, subsize, start, MPI_ORDER_FORTRAN &
                                         , MPI_DOUBLE_PRECISION, ddtype_dble_C_in_C2I(i), ierr )
            call MPI_TYPE_COMMIT(ddtype_dble_C_in_C2I(i),ierr)
                        
            bigsize(1) = n1-1
            bigsize(2) = n2msub
            bigsize(3) = n3msub_Isub
            subsize(1) = n1msubAll(i)
            subsize(2) = n2msub
            subsize(3) = n3msub_Isub
            start(1) = sum(n1msubAll(0:i)) - n1msubAll(i)
            start(2) = 0
            start(3) = 0
                        
            call MPI_TYPE_CREATE_SUBARRAY( 3, bigsize, subsize, start, MPI_ORDER_FORTRAN &
                                         , MPI_DOUBLE_PRECISION, ddtype_dble_I_in_C2I(i), ierr )
            call MPI_TYPE_COMMIT(ddtype_dble_I_in_C2I(i),ierr)
        enddo

        ! DDT for I-C transpose (complex type)
        do i=0,comm_1d_x1%nprocs-1
            bigsize(1) = h1psub
            bigsize(2) = n2msub
            bigsize(3) = n3msub
            subsize(1) = h1psub
            subsize(2) = n2msub
            subsize(3) = n3msub_IsubAll(i)
            start(1) = 0
            start(2) = 0
            start(3) = sum(n3msub_IsubAll(0:i)) - n3msub_IsubAll(i)
            call MPI_TYPE_CREATE_SUBARRAY( 3, bigsize, subsize, start, MPI_ORDER_FORTRAN &
                                         , MPI_DOUBLE_COMPLEX, ddtype_cplx_C_in_C2I(i), ierr )
            call MPI_TYPE_COMMIT(ddtype_cplx_C_in_C2I(i),ierr)
                                                
            bigsize(1) = h1p
            bigsize(2) = n2msub
            bigsize(3) = n3msub_Isub
            subsize(1) = h1psubAll(i)
            subsize(2) = n2msub
            subsize(3) = n3msub_Isub
            start(1) = sum(h1psubAll(0:i)) - h1psubAll(i)
            start(2) = 0
            start(3) = 0
            call MPI_TYPE_CREATE_SUBARRAY( 3, bigsize, subsize, start, MPI_ORDER_FORTRAN &
                                         , MPI_DOUBLE_COMPLEX, ddtype_cplx_I_in_C2I(i), ierr )
            call MPI_TYPE_COMMIT(ddtype_cplx_I_in_C2I(i),ierr)
        enddo

        ! DDT for I-K transpose (complex type)
        do i=0,comm_1d_x3%nprocs-1
            bigsize(1) = h1psub
            bigsize(2) = n2msub
            bigsize(3) = n3msub
            subsize(1) = h1psub_KsubAll(i)
            subsize(2) = n2msub
            subsize(3) = n3msub
            start(1) = sum(h1psub_KsubAll(0:i)) - h1psub_KsubAll(i)
            start(2) = 0
            start(3) = 0
            call MPI_TYPE_CREATE_SUBARRAY( 3, bigsize, subsize, start, MPI_ORDER_FORTRAN &
                                         , MPI_DOUBLE_COMPLEX, ddtype_cplx_C_in_C2K(i), ierr )
            call MPI_TYPE_COMMIT(ddtype_cplx_C_in_C2K(i),ierr)
                                                
            bigsize(1) = h1psub_Ksub
            bigsize(2) = n2msub
            bigsize(3) = n3-1
            subsize(1) = h1psub_Ksub
            subsize(2) = n2msub
            subsize(3) = n3msubAll(i)
            start(1) = 0
            start(2) = 0
            start(3) = sum(n3msubAll(0:i)) - n3msubAll(i)
            call MPI_TYPE_CREATE_SUBARRAY( 3, bigsize, subsize, start, MPI_ORDER_FORTRAN &
                                         , MPI_DOUBLE_COMPLEX, ddtype_cplx_K_in_C2K(i), ierr )
            call MPI_TYPE_COMMIT(ddtype_cplx_K_in_C2K(i),ierr)
        enddo
        deallocate(n3msub_IsubAll,n1msubAll,h1psubAll,h1psub_KsubAll,n3msubAll)

    end subroutine mpi_subdomain_DDT_transpose2

    !>
    !> @brief       Update the values of boundary ghostcells through communication in all directions.
    !>
    subroutine mpi_subdomain_ghostcell_update(Value_sub)

        use mpi_topology,   only : comm_1d_x1, comm_1d_x2, comm_1d_x3

        implicit none

        double precision, dimension(0:n1sub, 0:n2sub, 0:n3sub), intent(inout) :: Value_sub

        integer :: ierr
        integer :: request(4)

        ! Update the ghostcells in the x-direction using derived datatypes and subcommunicator.
        call MPI_Isend(Value_sub,1, ddtype_sendto_E  , comm_1d_x1%east_rank, 111, comm_1d_x1%mpi_comm, request(1), ierr)
        call MPI_Irecv(Value_sub,1, ddtype_recvfrom_W, comm_1d_x1%west_rank, 111, comm_1d_x1%mpi_comm, request(2), ierr)
        call MPI_Isend(Value_sub,1, ddtype_sendto_W  , comm_1d_x1%west_rank, 222, comm_1d_x1%mpi_comm, request(3), ierr)
        call MPI_Irecv(Value_sub,1, ddtype_recvfrom_E, comm_1d_x1%east_rank, 222, comm_1d_x1%mpi_comm, request(4), ierr)
        call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierr)
        
        ! Update the ghostcells in the y-direction using derived datatypes and subcommunicator.
        call MPI_Isend(Value_sub,1, ddtype_sendto_N  , comm_1d_x2%east_rank, 111, comm_1d_x2%mpi_comm, request(1), ierr)
        call MPI_Irecv(Value_sub,1, ddtype_recvfrom_S, comm_1d_x2%west_rank, 111, comm_1d_x2%mpi_comm, request(2), ierr)
        call MPI_Isend(Value_sub,1, ddtype_sendto_S  , comm_1d_x2%west_rank, 222, comm_1d_x2%mpi_comm, request(3), ierr)
        call MPI_Irecv(Value_sub,1, ddtype_recvfrom_N, comm_1d_x2%east_rank, 222, comm_1d_x2%mpi_comm, request(4), ierr)
        call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierr)

        ! Update the ghostcells in the z-direction using derived datatypes and subcommunicator.
        call MPI_Isend(Value_sub,1, ddtype_sendto_F  , comm_1d_x3%east_rank, 111, comm_1d_x3%mpi_comm, request(1), ierr)
        call MPI_Irecv(Value_sub,1, ddtype_recvfrom_B, comm_1d_x3%west_rank, 111, comm_1d_x3%mpi_comm, request(2), ierr)
        call MPI_Isend(Value_sub,1, ddtype_sendto_B  , comm_1d_x3%west_rank, 222, comm_1d_x3%mpi_comm, request(3), ierr)
        call MPI_Irecv(Value_sub,1, ddtype_recvfrom_F, comm_1d_x3%east_rank, 222, comm_1d_x3%mpi_comm, request(4), ierr)
        call MPI_Waitall(4, request, MPI_STATUSES_IGNORE, ierr)
        
        
    end subroutine mpi_subdomain_ghostcell_update

    !> @brief       Compute the indices of the assigned range for each MPI process .
    !> @param       nsta        First index of total range
    !> @param       nend        Last index of total range
    !> @param       nprocs      Number of MPI process
    !> @param       myrank      Rank ID of my MPI process
    !> @param       index_sta   First index of assigned range for myrank
    !> @param       index_end   Last index of assigned range for myrank
    subroutine subdomain_para_range(nsta, nend, nprocs, myrank, index_sta, index_end)

        implicit none

        integer, intent(in)     :: nsta, nend, nprocs, myrank
        integer, intent(out)    :: index_sta, index_end
        integer :: iwork1, iwork2

        iwork1 = int((nend - nsta + 1) / nprocs)
        iwork2 = mod(nend - nsta + 1, nprocs)
        index_sta = myrank * iwork1 + nsta + min(myrank, iwork2)
        index_end = index_sta + iwork1 - 1
        if (iwork2 > myrank) index_end = index_end + 1

    end subroutine subdomain_para_range

end module mpi_subdomain