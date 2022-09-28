!======================================================================================================================
!> @file        module_mpi_topology.f90
!> @brief       This file contains a module of communication topology for PaScaL_TCS.
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
!> @brief       Module for creating the cartesian topology of the MPI processes and subcommunicators.
!> @details     This module has three subcommunicators in each-direction and related subroutines.
!>
module mpi_topology

    use mpi

    implicit none

    integer :: mpi_world_cart       !< Communicator for cartesian topology
    integer :: nprocs               !< World size
    integer :: myrank               !< World rank

    !> @brief   Type variable for the information of 1D communicator
    type, public :: cart_comm_1d
        integer :: myrank           !< Rank ID in current communicator
        integer :: nprocs                   !< Number of processes in current communicator
        integer :: west_rank                !< Previous rank ID in current communicator
        integer :: east_rank                !< Next rank ID in current communicator
        integer :: mpi_comm                 !< Current communicator
    end type cart_comm_1d

    type(cart_comm_1d)  :: comm_1d_x1       !< Subcommunicator information in x-direction
    type(cart_comm_1d)  :: comm_1d_x2       !< Subcommunicator information in x-direction
    type(cart_comm_1d)  :: comm_1d_x3       !< Subcommunicator information in x-direction
    type(cart_comm_1d)  :: comm_1d_x1n2     !< Subcommunicator information in x-y direction for data transpose

    public  :: mpi_topology_make
    public  :: mpi_topology_clean

    contains

    !>
    !> @brief       Destroy the communicator for cartesian topology.
    !>
    subroutine mpi_topology_clean()

        implicit none
        integer :: ierr

        call MPI_Comm_free(mpi_world_cart, ierr)

    end subroutine mpi_topology_clean

    !>
    !> @brief       Create the cartesian topology for the MPI processes and subcommunicators.
    !>
    subroutine mpi_topology_make()
    
        use global, only : np1, np2, np3, pbc1, pbc2, pbc3

        implicit none

        integer :: np_dim(0:2)
        logical :: period(0:2)
        logical :: remain(0:2)
        integer :: ierr

        np_dim(0) = np1;    np_dim(1) = np2;    np_dim(2) = np3
        period(0) = pbc1;   period(1) = pbc2;   period(2) = pbc3 !: period(:) <- 'mpi_topology'

        ! Create the cartesian topology.
        call MPI_Cart_create( MPI_COMM_WORLD    &!  input  | integer      | Input communicator (handle).
                            , 3                 &!  input  | integer      | Number of dimensions of Cartesian grid (integer).
                            , np_dim            &!  input  | integer(1:3) | Integer array of size ndims specifying the number of processes in each dimension.
                            , period            &!  input  | logical(1:3) | Logical array of size ndims specifying whether the grid is periodic (true=1) or not (false=0) in each dimension.
                            , .false.           &!  input  | logical      | Ranking may be reordered (true=1) or not (false=0) (logical).
                            , mpi_world_cart    &! *output | integer      | Communicator with new Cartesian topology (handle).
                            , ierr              &!  output | integer      | Fortran only: Error status
                            )

        ! Create subcommunicators and assign two neighboring processes in the x-direction.
        remain(0) = .true.
        remain(1) = .false.
        remain(2) = .false.
        call MPI_Cart_sub( mpi_world_cart, remain, comm_1d_x1%mpi_comm, ierr)
        call MPI_Comm_rank(comm_1d_x1%mpi_comm, comm_1d_x1%myrank, ierr)
        call MPI_Comm_size(comm_1d_x1%mpi_comm, comm_1d_x1%nprocs, ierr)
        call MPI_Cart_shift(comm_1d_x1%mpi_comm, 0, 1, comm_1d_x1%west_rank, comm_1d_x1%east_rank, ierr)

        ! Create subcommunicators and assign two neighboring processes in the y-direction
        remain(0) = .false.
        remain(1) = .true.
        remain(2) = .false.
        call MPI_Cart_sub( mpi_world_cart, remain, comm_1d_x2%mpi_comm, ierr)
        call MPI_Comm_rank(comm_1d_x2%mpi_comm, comm_1d_x2%myrank, ierr)
        call MPI_Comm_size(comm_1d_x2%mpi_comm, comm_1d_x2%nprocs, ierr)
        call MPI_Cart_shift(comm_1d_x2%mpi_comm, 0, 1, comm_1d_x2%west_rank, comm_1d_x2%east_rank, ierr)

        ! Create subcommunicators and assign two neighboring processes in the z-direction
        remain(0) = .false.
        remain(1) = .false.
        remain(2) = .true.
        call MPI_Cart_sub( mpi_world_cart, remain, comm_1d_x3%mpi_comm, ierr)
        call MPI_Comm_rank(comm_1d_x3%mpi_comm, comm_1d_x3%myrank, ierr)
        call MPI_Comm_size(comm_1d_x3%mpi_comm, comm_1d_x3%nprocs, ierr)
        call MPI_Cart_shift(comm_1d_x3%mpi_comm, 0, 1, comm_1d_x3%west_rank, comm_1d_x3%east_rank, ierr)

        ! For Possion
        call MPI_Comm_split(mpi_world_cart, comm_1d_x3%myrank, comm_1d_x1%myrank+comm_1d_x2%myrank*comm_1d_x1%nprocs, comm_1d_x1n2%mpi_comm, ierr)
        call MPI_Comm_rank(comm_1d_x1n2%mpi_comm, comm_1d_x1n2%myrank, ierr)
        call MPI_Comm_size(comm_1d_x1n2%mpi_comm, comm_1d_x1n2%nprocs, ierr)

    end subroutine mpi_topology_make

end module mpi_topology
