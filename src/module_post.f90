!======================================================================================================================
!> @file        module_post.f90
!> @brief       This file contains a module of post-treatment for PaScaL_TCS.
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
!> @brief       Module for post-treatment for PaScaL_TCS.
!> @details     This module contains the followings : 
!>              1. Divergency and CFL number
!>              2. Monitoring the output parameters of interest
!>              3. File IO with explicit aggregation
!>
module mpi_Post

    use mpi_topology
    use global
    use mpi_subdomain

    implicit none

    private

    !> @{ Control flag for output print
    logical :: FirstPrint1=.true.,FirstPrint2=.true.,FirstPrint3=.true.
    !> @}
    
    !> @{ Function for divergence, CFL, and output print
    public :: mpi_Post_Div
    public :: mpi_Post_CFL
    public :: mpi_Post_MonitorOut
    !> @}

    !> @{ MPI derived datatype for file IO
    integer :: ddtype_inner_domain
    integer :: ddtype_global_domain_pr_IO
    integer :: ddtype_global_domain_MPI_IO
    integer :: ddtype_global_domain_MPI_IO_aggregation
    integer :: ddtype_aggregation_chunk
    !> @}

    !> @{ Temporary variables for building DDT for file IO
    integer, allocatable, DIMENSION(:)  :: cnts_pr, disp_pr, cnts_aggr, disp_aggr
    !> @}

    !> @{ Communicator type for aggregated file IO
    type, private :: comm_IO
        integer :: myrank
        integer :: nprocs
        integer :: mpi_comm
    end type comm_IO
    !> @}

    type(comm_IO)  :: comm_IO_aggregation   !> Communicator for aggregation
    type(comm_IO)  :: comm_IO_master        !> Communicator for master process having aggregated data

    !> @{ Function for file IO
    public :: mpi_Post_allocation

    public :: mpi_Post_FileOut_InstanField              ! Instanteneous filed print in ascii format
    
    public :: mpi_Post_FileIn_Continue_Single_Process_IO
    public :: mpi_Post_FileOut_Continue_Single_Process_IO

    public :: mpi_Post_FileIn_Continue_Single_Process_IO_with_Aggregation
    public :: mpi_Post_FileOut_Continue_Single_Process_IO_with_Aggregation

    public :: mpi_Post_FileIn_Continue_MPI_IO
    public :: mpi_Post_FileOut_Continue_MPI_IO

    public :: mpi_Post_FileIn_Continue_MPI_IO_with_Aggregation
    public :: mpi_Post_FileOut_Continue_MPI_IO_with_Aggregation

    public :: mpi_Post_FileIn_Continue_Post_Reassembly_IO
    public :: mpi_Post_FileOut_Continue_Post_Reassembly_IO

    public :: mpi_Post_WallShearStress
    !> @}

    character(len=17) :: xyzrank            !> Rank information in filename

    contains

    !>
    !> @brief       Calculate maximum divergence of velocity
    !> @param       U           Velocity field in x-direction
    !> @param       V           Velocity field in y-direction
    !> @param       W           Velocity field in z-direction
    !> @param       maxDivU     Maximum divergence of velocity
    !>
    subroutine mpi_Post_Div(U,V,W,maxDivU)

        use MPI

        implicit none

        double precision, dimension(0:n1sub,0:n2sub,0:n3sub) :: U,V,W
        double precision :: maxDivU

        !> @{ Local indexing variable
        integer :: i,j,k
        integer :: im,jm,km
        integer :: ip,jp,kp
        integer :: Sc,Ec
        integer :: Sm,Em
        integer :: Sp,Ep
        !> @}
        
        !> @{ Temporary reduce variable for (maybe) better allreduce communication
        double precision :: maxI,maxIJ,maxIJK
        double precision, allocatable, dimension(:) :: findmaxk
        double precision, allocatable, dimension(:,:) :: findmaxIJ
        !> @}

        integer :: ierr

        allocate(findmaxIJ(1:n1sub-1,1:n2sub-1),findmaxk(1:n3sub-1))
        
        ! Define range
        Sc=1  ;Ec=n1sub-1
        Sm=1-1;Em=n1sub-1-1
        Sp=1+1;Ep=n1sub-1+1
        
        ! Find maximum divergence in the subdomain with plane-by-plane comparison
        do k = 1, n3sub-1
            kp = k+1
            km = k-1
            do j = 1, n2sub-1
                jp = j + 1
                jm = j - 1
            
                findmaxIJ(Sc:Ec,j) = ( U(Sp:Ep,j ,k ) - U(Sc:Ec,j,k)  )/dx1_sub(Sc:Ec)    &
                                + ( V(Sc:Ec,jp,k ) - V(Sc:Ec,j,k)  )/dx2_sub(j)        &
                                + ( W(Sc:Ec,j ,kp) - W(Sc:Ec,j,k)  )/dx3_sub(k)

            end do
            findmaxk(k)=maxval(findmaxIJ)
        end do

        ! Find maximum divergence in the global domain with three step communication
        maxDivU=maxval(findmaxk)
        call MPI_ALLREDUCE(maxDivU, maxI  , 1, MPI_DOUBLE, MPI_MAX, comm_1d_x1%mpi_comm, ierr)
        call MPI_ALLREDUCE(maxI   , maxIJ , 1, MPI_DOUBLE, MPI_MAX, comm_1d_x2%mpi_comm, ierr)
        call MPI_ALLREDUCE(maxIJ  , maxIJK, 1, MPI_DOUBLE, MPI_MAX, comm_1d_x3%mpi_comm, ierr)
        maxDivU=maxIJK

        deallocate(findmaxIJ,findmaxk)

    end subroutine mpi_Post_Div

    !>
    !> @brief       Calculate maximum CFL number and modify time step
    !> @param       U           Velocity field in x-direction
    !> @param       V           Velocity field in y-direction
    !> @param       W           Velocity field in z-direction
    !> @param       newCFL      Maximum CFL number
    !> @param       newdt       Modified dt
    !>
    subroutine mpi_Post_CFL(U,V,W,newCFL,newdt)

        use MPI

        implicit none
        double precision, dimension(0:n1sub,0:n2sub,0:n3sub) :: U,V,W
        double precision :: newCFL,newdt

        !> @{ Local indexing variable
        integer :: i,j,k
        integer :: im,jm,km
        integer :: ip,jp,kp
        integer :: Sc,Ec
        integer :: Sm,Em
        integer :: Sp,Ep
        !> @}

        !> @{ Temporary reduce variable for (maybe) better allreduce communication
        double precision :: maxI,maxIJ,maxIJK,maxUovX
        double precision, allocatable, dimension(:) :: findmaxk
        double precision, allocatable, dimension(:,:) :: findmaxIJ
        !> @}
        integer :: ierr

        allocate(findmaxIJ(1:n1sub-1,1:n2sub-1),findmaxk(1:n3sub-1))
        
        ! Define range
        Sc=1  ;Ec=n1sub-1
        Sm=1-1;Em=n1sub-1-1
        Sp=1+1;Ep=n1sub-1+1
        
        ! Find maximum CFL number in the subdomain with plane-by-plane comparison
        do k = 1, n3sub-1
            kp = k+1
            km = k-1
            do j = 1, n2sub-1
                jp = j + 1
                jm = j - 1
            
                findmaxIJ(Sc:Ec,j) = (ABS( U(Sp:Ep,j ,k ) + U(Sc:Ec,j,k)  )/2.d0)/dx1_sub(Sc:Ec)    &
                                + (ABS( V(Sc:Ec,jp,k ) + V(Sc:Ec,j,k)  )/2.d0)/dx2_sub(j)        &
                                + (ABS( W(Sc:Ec,j ,kp) + W(Sc:Ec,j,k)  )/2.d0)/dx3_sub(k)

            end do
            findmaxk(k)=maxval(findmaxIJ)
        end do

        ! Find maximum CFL number in the global domain with three step communication
        maxUovX=maxval(findmaxk)
        call MPI_ALLREDUCE(maxUovX, maxI  , 1, MPI_DOUBLE, MPI_MAX, comm_1d_x1%mpi_comm, ierr)
        call MPI_ALLREDUCE(maxI   , maxIJ , 1, MPI_DOUBLE, MPI_MAX, comm_1d_x2%mpi_comm, ierr)
        call MPI_ALLREDUCE(maxIJ  , maxIJK, 1, MPI_DOUBLE, MPI_MAX, comm_1d_x3%mpi_comm, ierr)
        maxUovX=maxIJK

        ! Calculate CFL number and modify dt
        newCFL=newdt*maxUovX
        newdt=MaxCFL/maxUovX
        if(newdt>dtMax) newdt=dtMax
        newCFL=newdt*maxUovX

        deallocate(findmaxIJ,findmaxk)

    end subroutine mpi_Post_CFL

    !>
    !> @brief       Print the output parameters for monitoring
    !> @param       myrank          Rank for print
    !> @param       outTimeStep     Current time step number
    !> @param       outtime         Dimensionless time
    !> @param       outdt           dt
    !> @param       outCFL          CFL number
    !> @param       outmaxDivU      Maximum divergence
    !> @param       outtimer        Wall-clock time of current time step
    !>
    subroutine mpi_Post_MonitorOut(myrank,outTimeStep,outtime,outdt,outCFL,outmaxDivU,wss,outtimer)

        implicit none

        integer :: outTimeStep,myrank
        double precision :: outtime,outdt,outCFL,outmaxDivU,wss,outtimer

        ! Print output paramters every 20 steps
        if(mod(outTimeStep,20)==1) then
            if(myrank==0) write(*,*)
            if(myrank==0) write(*,'(7A15)') 'Timestep', 'Time', 'dt', 'CFL', 'max_DivU', 'WSS', 'WTime/step'
        endif
        if(myrank==0) write(*,'(1I15,6E15.5)') outTimeStep,outtime,outdt,outCFL,outmaxDivU,wss,outtimer

    end subroutine mpi_Post_MonitorOut

    !>
    !> @brief       Write instantaneous field variables for post-processing in ascii format to file
    !> @detail      An example subroutine
    !> @param       myrank          Rank for print
    !> @param       tstepin         Current time step number for file write
    !> @param       timein          Current dimensionless time for file write
    !> @param       Uin             Velocity field in x-direction for file write
    !> @param       Vin             Velocity field in y-direction for file write
    !> @param       Win             Velocity field in z-direction for file write
    !> @param       Pin             Pressure field for file write
    !> @param       Tin             Temperature field for file write
    !>
    subroutine mpi_Post_FileOut_InstanField(myrank,tstepin,timein,Uin,Vin,Win,Pin,Tin)

        use mpi
        
        implicit none
        
        integer :: myrank,tstepin
        double precision :: timein
        double precision, dimension(0:n1sub, 0:n2sub, 0:n3sub) :: Uin,Vin,Win,Pin,Tin

        !> @{ File name
        character(len=22) :: filename_instantfieldXY
        character(len=27) :: filename_instantfieldXZ_wall
        !> @}

        integer :: i,j,k

        ! Define range
        filename_instantfieldXY='Output_instantfield_XY'
        filename_instantfieldXZ_wall='Output_instantfield_XZ_wall'

        ! Time step for file write
        if(tstepin==0.or.(tstepin>=print_start_step.and.mod(tstepin-print_start_step,print_interval_step)==0)) then

        ! Write field variables in X-Y plane
            if(comm_1d_x3%myrank==comm_1d_x3%nprocs-1) then
                open(unit=myrank,file=dir_instantfield//filename_instantfieldXY//xyzrank//'.plt', position='append')
                    if(FirstPrint1) then
                        write(myrank,*) 'VARIABLES="X","Y","Z","U","V","W","P","T"'  !--
                        write(myrank,*) 'zone t="',tstepin,'"','i=',n1sub+1,'j=',n2sub+1,'k=',1
                        write(myrank,*) 'STRANDID=', tstepin
                        write(myrank,*) 'SOLUTIONTIME=', timein
                        k= (n3sub-1)/2
                        do j=0,n2sub
                        do i=0,n1sub
                            write(myrank,'(3D20.10,5D30.20)') x1_sub(i),x2_sub(j),x3_sub(k)    &
                                                        &, Uin(i,j,k),Vin(i,j,k),Win(i,j,k),Pin(i,j,k),Tin(i,j,k)
                        enddo
                        enddo

                        FirstPrint1=.false.
                    else
                        write(myrank,*) 'zone t="',tstepin,'"','i=',n1sub+1,'j=',n2sub+1,'k=',1
                        write(myrank,*) 'VARSHARELIST= ([1-3]=1)'
                        write(myrank,*) 'STRANDID=', tstepin
                        write(myrank,*) 'SOLUTIONTIME=', timein
                        k= (n3sub-1)/2
                        do j=0,n2sub
                        do i=0,n1sub
                            write(myrank,'(5D30.20)') Uin(i,j,k),Vin(i,j,k),Win(i,j,k),Pin(i,j,k),Tin(i,j,k)
                        enddo
                        enddo

                    endif
                close(myrank)
            endif

        ! Write field variables in X-Z plane
            if(comm_1d_x2%myrank==0) then
                open(unit=myrank,file=dir_instantfield//filename_instantfieldXZ_wall//xyzrank//'.plt', position='append')
                    if(FirstPrint2) then
                        write(myrank,*) 'VARIABLES="X","Y","Z","U","V","W","P","T"'  !--
                        write(myrank,*) 'zone t="',tstepin,'"','i=',n1sub+1,'j=',1,'k=',n3sub+1
                        write(myrank,*) 'STRANDID=', tstepin
                        write(myrank,*) 'SOLUTIONTIME=', timein
                        j=3
                        do k=0,n3sub
                        do i=0,n1sub
                            write(myrank,'(3D20.10,5D30.20)') x1_sub(i),x2_sub(j),x3_sub(k)    &
                                                        &, Uin(i,j,k),Vin(i,j,k),Win(i,j,k),Pin(i,j,k),Tin(i,j,k)
                        enddo
                        enddo

                        FirstPrint2=.false.
                    else
                        write(myrank,*) 'zone t="',tstepin,'"','i=',n1sub+1,'j=',1,'k=',n3sub+1
                        write(myrank,*) 'VARSHARELIST= ([1-3]=1)'
                        write(myrank,*) 'STRANDID=', tstepin
                        write(myrank,*) 'SOLUTIONTIME=', timein
                        j=3
                        do k=0,n3sub
                        do i=0,n1sub
                            write(myrank,'(5D30.20)') Uin(i,j,k),Vin(i,j,k),Win(i,j,k),Pin(i,j,k),Tin(i,j,k)
                        enddo
                        enddo

                    endif
                close(myrank)
            endif
        endif

    end subroutine mpi_Post_FileOut_InstanField

    !>
    !> @brief       Initialize the required variables for file IO of field variable
    !> @detail      File IO of field variables is for restart, so it writes whole field variables
    !>              to binary file.
    !> @param       chunk_size  Aggregation size of MPI processes
    !>
    subroutine mpi_Post_allocation(chunk_size)

        implicit none

        integer :: chunk_size

        !> @{ Temporary local variables for communicator setup
        integer :: color
        integer :: i, j, k, icur, ierr
        integer :: sizes(3), subsizes(3), starts(3)

        integer :: ddtype_temp
        integer :: r8size
        integer :: mpi_world_group, mpi_master_group
        !> @}

        !> @{ Temporary local variables for communicator setup
        integer, allocatable, dimension(:,:)    :: cart_coord
        integer, allocatable, dimension(:)      :: n1msub_cnt,  n2msub_cnt,  n3msub_cnt
        integer, allocatable, dimension(:)      :: n1msub_disp, n2msub_disp, n3msub_disp
        integer, allocatable, dimension(:)      :: maste_group_rank
        integer(kind=MPI_ADDRESS_KIND)          :: extent, lb
        !> @}

        ! String of rank ID in x, y, and z-direction
        write(xyzrank, '(I5.5,1A1,I5.5,1A1,I5.5)' ) comm_1d_x1%myrank,'_',   &
                                                    comm_1d_x2%myrank,'_',   &
                                                    comm_1d_x3%myrank

        
        !>>>>>>>>>> MPI communicator for data aggregation (comm_IO_aggregation)
        ! comm_IO_aggregation includes MPI processes of chunk_size (mpi_size_aggregation) in z-direction.
        comm_IO_aggregation%nprocs = chunk_size

        if( mod(comm_1d_x3%nprocs, comm_IO_aggregation%nprocs).ne.0 ) then
            if( myrank.eq.0) print *, '[Error] Chunk_size for IO aggregation should be a measure of mpisize in z-direction'
            call MPI_Abort(MPI_COMM_WORLD, 11, ierr)
        endif

        color = myrank / comm_IO_aggregation%nprocs

        call MPI_Comm_split(MPI_COMM_WORLD, color, myrank, comm_IO_aggregation%mpi_comm, ierr )
        call MPI_Comm_rank(comm_IO_aggregation%mpi_comm, comm_IO_aggregation%myrank, ierr)

        !>>>>>>>>>> MPI communicator for aggregation_master (comm_IO_master)
        ! comm_IO_master includes rank0 processes in comm_IO_aggregation
        comm_IO_master%mpi_comm = MPI_COMM_NULL
        comm_IO_master%nprocs = nprocs / comm_IO_aggregation%nprocs
        allocate( maste_group_rank(comm_IO_master%nprocs) )
        do i = 1, comm_IO_master%nprocs
            maste_group_rank(i) = (i-1) * comm_IO_aggregation%nprocs
        enddo
        
        call MPI_Comm_group(MPI_COMM_WORLD, mpi_world_group, ierr)
        call MPI_Group_incl(mpi_world_group, comm_IO_master%nprocs, maste_group_rank, mpi_master_group, ierr)
        call MPI_Comm_create_group(MPI_COMM_WORLD, mpi_master_group, 0, comm_IO_master%mpi_comm, ierr)

        if(comm_IO_master%mpi_comm.ne.MPI_COMM_NULL) then
            call MPI_Comm_rank(comm_IO_master%mpi_comm, comm_IO_master%myrank, ierr)
        endif
        deallocate( maste_group_rank )
                                            
        !>>>>>>>>>> Derived datatype for inner domain without ghost cells
        sizes    = (/ n1msub+2, n2msub+2, n3msub+2 /)
        subsizes = (/ n1msub,   n2msub,   n3msub   /)
        starts   = (/ 1, 1, 1 /)

        call MPI_Type_create_subarray(  3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                        MPI_DOUBLE_PRECISION, ddtype_inner_domain, ierr)
        call MPI_Type_commit(ddtype_inner_domain, ierr)

        !>>>>>>>>>> Derived datatype for post reassembly IO
        ! Post reassembly can be used when a single node is capable of memory allocation for the global domain.
        ! All data is gathered into rank 0 and scattered from rank 0 using MPI_Scatterv and MPI_Gatherv
        sizes    = (/ n1m,    n2m,    n3m    /)
        subsizes = (/ n1msub, n2msub, n3msub /)
        starts   = (/ 0, 0, 0 /)

        call MPI_Type_create_subarray(  3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                        MPI_DOUBLE_PRECISION, ddtype_temp, ierr)

        CALL MPI_Type_size(MPI_DOUBLE_PRECISION, r8size, ierr)
        lb = 0
        extent = r8size

        call MPI_Type_create_resized(ddtype_temp, lb, extent, ddtype_global_domain_pr_IO, ierr)
        call MPI_Type_commit(ddtype_global_domain_pr_IO, ierr)

        ! Counts and displacements for MPI_Gatherv and MPI_Scatterv 
        allocate( cart_coord(3,0:nprocs-1) )
        do i = 0, nprocs-1
            call MPI_Cart_coords(mpi_world_cart, i, 3, cart_coord(:,i), ierr )
        enddo

        allocate(n1msub_cnt(0:comm_1d_x1%nprocs-1), n1msub_disp(0:comm_1d_x1%nprocs-1))
        allocate(n2msub_cnt(0:comm_1d_x2%nprocs-1), n2msub_disp(0:comm_1d_x2%nprocs-1))
        allocate(n3msub_cnt(0:comm_1d_x3%nprocs-1), n3msub_disp(0:comm_1d_x3%nprocs-1))

        call MPI_Allgather(n1msub, 1, MPI_INTEGER, n1msub_cnt, 1, MPI_INTEGER, comm_1d_x1%mpi_comm, ierr)
        call MPI_Allgather(n2msub, 1, MPI_INTEGER, n2msub_cnt, 1, MPI_INTEGER, comm_1d_x2%mpi_comm, ierr)
        call MPI_Allgather(n3msub, 1, MPI_INTEGER, n3msub_cnt, 1, MPI_INTEGER, comm_1d_x3%mpi_comm, ierr)
    
        n1msub_disp(0) = 0
        do i = 1, comm_1d_x1%nprocs-1
            n1msub_disp(i)=sum(n1msub_cnt(0:i-1))
        enddo
    
        n2msub_disp(0) = 0
        do i = 1, comm_1d_x2%nprocs-1
            n2msub_disp(i)=sum(n2msub_cnt(0:i-1))
        enddo
    
        n3msub_disp(0) = 0
        do i = 1, comm_1d_x3%nprocs-1
            n3msub_disp(i)=sum(n3msub_cnt(0:i-1))
        enddo

        allocate( cnts_pr(0:nprocs-1) )
        allocate( disp_pr(0:nprocs-1) )
        
        do i = 0, nprocs-1
            cnts_pr(i) = 1
            disp_pr(i) =  n1msub_disp(cart_coord(1,i))  &
                        + n2msub_disp(cart_coord(2,i)) * n1m  &
                        + n3msub_disp(cart_coord(3,i)) * n1m * n2m
        enddo

        deallocate(cart_coord)
        deallocate(n1msub_cnt, n2msub_cnt, n3msub_cnt)
        deallocate(n1msub_disp, n2msub_disp, n3msub_disp)

        !>>>>>>>>>> Derived datatype for MPI IO
        sizes    = (/ n1m,    n2m,    n3m    /)
        subsizes = (/ n1msub, n2msub, n3msub /)
        starts   = (/ ista-1, jsta-1, ksta-1 /)

        call MPI_Type_create_subarray(  3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                        MPI_DOUBLE_PRECISION, ddtype_global_domain_MPI_IO, ierr)
        call MPI_Type_commit(ddtype_global_domain_MPI_IO, ierr)


        !>>>>>>>>>> Derived datatype for data aggregation.
        sizes    = (/ n1msub, n2msub, n3msub * comm_IO_aggregation%nprocs /)
        subsizes = (/ n1msub, n2msub, n3msub /)
        starts   = (/ 0, 0, n3msub * comm_IO_aggregation%myrank /)
        
        call MPI_Type_create_subarray(  3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                        MPI_DOUBLE_PRECISION, ddtype_temp, ierr)

        CALL MPI_Type_size(MPI_DOUBLE_PRECISION, r8size, ierr)
        lb = 0
        extent = r8size

        call MPI_Type_create_resized(ddtype_temp, lb, extent, ddtype_aggregation_chunk, ierr)
        call MPI_Type_commit(ddtype_aggregation_chunk, ierr)

        allocate( cnts_aggr(0:comm_IO_aggregation%nprocs-1) )
        allocate( disp_aggr(0:comm_IO_aggregation%nprocs-1) )

        do i = 0, comm_IO_aggregation%nprocs-1
            cnts_aggr(i) = 1
            disp_aggr(i) = n1msub * n2msub * n3msub * i
        enddo

        !>>>>>>>>>> Derived datatype for aggregated MPI IO.
        sizes    = (/ n1m,    n2m,    n3m    /)
        subsizes = (/ n1msub, n2msub, n3msub * comm_IO_aggregation%nprocs /)
        starts   = (/ ista-1, jsta-1, ksta-1 - n3msub * comm_IO_aggregation%myrank /)

        call MPI_Type_create_subarray(  3, sizes, subsizes, starts, MPI_ORDER_FORTRAN, &
                                        MPI_DOUBLE_PRECISION, ddtype_global_domain_MPI_IO_aggregation, ierr)
        call MPI_Type_commit(ddtype_global_domain_MPI_IO_aggregation, ierr)

    end subroutine mpi_Post_allocation

    !>
    !> @brief       Read field variables using single process IO
    !> @param       myrank          Rank for print
    !> @param       timein          Current dimensionless time for file write
    !> @param       dt              dt
    !> @param       Uin             Velocity field in x-direction for file write
    !> @param       Vin             Velocity field in y-direction for file write
    !> @param       Win             Velocity field in z-direction for file write
    !> @param       Pin             Pressure field for file write
    !> @param       Tin             Temperature field for file write
    !>
    subroutine mpi_Post_FileIn_Continue_Single_Process_IO(myrank,timein,dt,Uin,Vin,Win,Pin,Tin)

        implicit none

        integer :: myrank
        double precision :: timein,dt
        double precision, dimension(0:n1sub, 0:n2sub, 0:n3sub) :: Uin,Vin,Win,Pin,Tin

        integer :: i,j,k
        
        ! Read time step information
        open(myrank, FILE=trim(dir_cont_filein)//'cont_time_'//xyzrank//'.bin', FORM='UNFORMATTED', &
            STATUS='OLD', ACTION='READ', ACCESS='DIRECT', RECL=sizeof(timein)*2 )
        read(myrank, REC=1) timein,dt
        close(myrank)

        ! Read U
        open(myrank, FILE=trim(dir_cont_filein)//'cont_U_'//xyzrank//'.bin', FORM='UNFORMATTED', &
            STATUS='OLD', ACTION='READ', ACCESS='DIRECT', RECL=sizeof(Uin) )
        read(myrank, REC=1) Uin
        close(myrank)

        ! Read V
        open(myrank, FILE=trim(dir_cont_filein)//'cont_V_'//xyzrank//'.bin', FORM='UNFORMATTED', &
            STATUS='OLD', ACTION='READ', ACCESS='DIRECT', RECL=sizeof(Vin) )
        read(myrank, REC=1) Vin
        close(myrank)

        ! Read W
        open(myrank, FILE=trim(dir_cont_filein)//'cont_W_'//xyzrank//'.bin', FORM='UNFORMATTED', &
            STATUS='OLD', ACTION='READ', ACCESS='DIRECT', RECL=sizeof(Win) )
        read(myrank, REC=1) Win
        close(myrank)

        ! Read P
        open(myrank, FILE=trim(dir_cont_filein)//'cont_P_'//xyzrank//'.bin', FORM='UNFORMATTED', &
            STATUS='OLD', ACTION='READ', ACCESS='DIRECT', RECL=sizeof(Pin) )
        read(myrank, REC=1) Pin
        close(myrank)

        ! Read T
        open(myrank, FILE=trim(dir_cont_filein)//'cont_THETA_'//xyzrank//'.bin', FORM='UNFORMATTED', &
            STATUS='OLD', ACTION='READ', ACCESS='DIRECT', RECL=sizeof(Tin) )
        read(myrank, REC=1) Tin
        close(myrank)

        if(myrank.eq.0) print '(a)', 'Read continue file using single process IO'

    end subroutine mpi_Post_FileIn_Continue_Single_Process_IO

    !>
    !> @brief       write field variables using single process IO
    !> @param       myrank          Rank for print
    !> @param       timein          Current dimensionless time for file write
    !> @param       dt              dt
    !> @param       Uin             Velocity field in x-direction for file write
    !> @param       Vin             Velocity field in y-direction for file write
    !> @param       Win             Velocity field in z-direction for file write
    !> @param       Pin             Pressure field for file write
    !> @param       Tin             Temperature field for file write
    !>
    subroutine mpi_Post_FileOut_Continue_Single_Process_IO(myrank,timein,dt,Uin,Vin,Win,Pin,Tin)

        implicit none

        integer :: myrank
        double precision :: timein,dt
        double precision, dimension(0:n1sub, 0:n2sub, 0:n3sub) :: Uin,Vin,Win,Pin,Tin

        integer :: i,j,k
        
        ! Write time step information
        open(myrank, FILE=trim(dir_cont_fileout)//'cont_time_'//xyzrank//'.bin', FORM='UNFORMATTED', &
            STATUS='REPLACE', ACTION='WRITE', ACCESS='DIRECT', RECL=sizeof(timein)*2 )
        write(myrank, REC=1) timein,dt
        close(myrank)

        ! Write U
        open(myrank, FILE=trim(dir_cont_fileout)//'cont_U_'//xyzrank//'.bin', FORM='UNFORMATTED', &
            STATUS='REPLACE', ACTION='WRITE', ACCESS='DIRECT', RECL=sizeof(Uin) )
        write(myrank, REC=1) Uin
        close(myrank)

        ! Write V
        open(myrank, FILE=trim(dir_cont_fileout)//'cont_V_'//xyzrank//'.bin', FORM='UNFORMATTED', &
            STATUS='REPLACE', ACTION='WRITE', ACCESS='DIRECT', RECL=sizeof(Vin) )
        write(myrank, REC=1) Vin
        close(myrank)

        ! Write W
        open(myrank, FILE=trim(dir_cont_fileout)//'cont_W_'//xyzrank//'.bin', FORM='UNFORMATTED', &
            STATUS='REPLACE', ACTION='WRITE', ACCESS='DIRECT', RECL=sizeof(Win) )
        write(myrank, REC=1) Win
        close(myrank)

        ! Write P
        open(myrank, FILE=trim(dir_cont_fileout)//'cont_P_'//xyzrank//'.bin', FORM='UNFORMATTED', &
            STATUS='REPLACE', ACTION='WRITE', ACCESS='DIRECT', RECL=sizeof(Pin) )
        write(myrank, REC=1) Pin
        close(myrank)

        ! Write T
        open(myrank, FILE=trim(dir_cont_fileout)//'cont_THETA_'//xyzrank//'.bin', FORM='UNFORMATTED', &
            STATUS='REPLACE', ACTION='WRITE', ACCESS='DIRECT', RECL=sizeof(Tin) )
        write(myrank, REC=1) Tin
        close(myrank)

        if(myrank.eq.0) print '(a)', 'Write continue file using single process IO'

    end subroutine mpi_Post_FileOut_Continue_Single_Process_IO

    !>
    !> @brief       Read field variables using single process IO with explicit aggregation
    !> @param       myrank          Rank for print
    !> @param       timein          Current dimensionless time for file write
    !> @param       dt              dt
    !> @param       Uin             Velocity field in x-direction for file write
    !> @param       Vin             Velocity field in y-direction for file write
    !> @param       Win             Velocity field in z-direction for file write
    !> @param       Pin             Pressure field for file write
    !> @param       Tin             Temperature field for file write
    !>
    subroutine mpi_Post_FileIn_Continue_Single_Process_IO_with_Aggregation(myrank,timein,dt,Uin,Vin,Win,Pin,Tin)
        use mpi
        implicit none
        integer :: myrank
        double precision :: timein,dt
        double precision, dimension(0:n1sub, 0:n2sub, 0:n3sub) :: Uin,Vin,Win,Pin,Tin

        integer :: ierr
        double precision, allocatable, dimension(:,:,:) :: var_chunk
        
        ! The master rank for data aggregation allocates chunk array which have aggregated data
        if(comm_IO_aggregation%myrank.eq.0) then
            allocate( var_chunk(n1msub, n2msub, n3msub * comm_IO_aggregation%nprocs) )
        endif

        ! The master rank for data aggregation reads time step information
        if(comm_IO_aggregation%myrank.eq.0) then
            open(myrank, FILE=trim(dir_cont_filein)//'cont_time_'//xyzrank//'.bin', FORM='UNFORMATTED', &
                STATUS='OLD', ACTION='READ', ACCESS='DIRECT', RECL=sizeof(timein)*2 )
            read(myrank, REC=1) timein,dt
            close(myrank)
        endif
        ! The master rank for data aggregation broadcast its local communicator for aggregation
        call MPI_Bcast(timein, 1, MPI_DOUBLE_PRECISION, 0, comm_IO_aggregation%mpi_comm, ierr)
        call MPI_Bcast(dt,     1, MPI_DOUBLE_PRECISION, 0, comm_IO_aggregation%mpi_comm, ierr)

        ! The master rank for data aggregation reads U
        if(comm_IO_aggregation%myrank.eq.0) then
            open(myrank, FILE=trim(dir_cont_filein)//'cont_U_'//xyzrank//'.bin', FORM='UNFORMATTED', &
                STATUS='OLD', ACTION='READ', ACCESS='DIRECT', RECL=sizeof(var_chunk) )
            read(myrank, REC=1) var_chunk
            close(myrank)
        endif

        ! The master rank for data aggregation scatter the readed U to its local communicator using the defined DDT
        call MPI_Scatterv( var_chunk, cnts_aggr, disp_aggr, ddtype_aggregation_chunk, Uin, 1, ddtype_inner_domain, 0, &
                        comm_IO_aggregation%mpi_comm, ierr)

        ! The master rank for data aggregation reads V
        if(comm_IO_aggregation%myrank.eq.0) then
            open(myrank, FILE=trim(dir_cont_filein)//'cont_V_'//xyzrank//'.bin', FORM='UNFORMATTED', &
                STATUS='OLD', ACTION='READ', ACCESS='DIRECT', RECL=sizeof(var_chunk) )
            read(myrank, REC=1) var_chunk
            close(myrank)
        endif

        ! The master rank for data aggregation scatter the readed V to its local communicator using the defined DDT
        call MPI_Scatterv( var_chunk, cnts_aggr, disp_aggr, ddtype_aggregation_chunk, Vin, 1, ddtype_inner_domain, 0, &
                        comm_IO_aggregation%mpi_comm, ierr)

        ! The master rank for data aggregation reads W
        if(comm_IO_aggregation%myrank.eq.0) then
            open(myrank, FILE=trim(dir_cont_filein)//'cont_W_'//xyzrank//'.bin', FORM='UNFORMATTED', &
                STATUS='OLD', ACTION='READ', ACCESS='DIRECT', RECL=sizeof(var_chunk) )
            read(myrank, REC=1) var_chunk
            close(myrank)
        endif

        ! The master rank for data aggregation scatter the readed W to its local communicator using the defined DDT
        call MPI_Scatterv( var_chunk, cnts_aggr, disp_aggr, ddtype_aggregation_chunk, Win, 1, ddtype_inner_domain, 0, &
                        comm_IO_aggregation%mpi_comm, ierr)

        ! The master rank for data aggregation reads P
        if(comm_IO_aggregation%myrank.eq.0) then
            open(myrank, FILE=trim(dir_cont_filein)//'cont_P_'//xyzrank//'.bin', FORM='UNFORMATTED', &
                STATUS='OLD', ACTION='READ', ACCESS='DIRECT', RECL=sizeof(var_chunk) )
            read(myrank, REC=1) var_chunk
            close(myrank)
        endif

        ! The master rank for data aggregation scatter the readed P to its local communicator using the defined DDT
        call MPI_Scatterv( var_chunk, cnts_aggr, disp_aggr, ddtype_aggregation_chunk, Pin, 1, ddtype_inner_domain, 0, &
                        comm_IO_aggregation%mpi_comm, ierr)

        ! The master rank for data aggregation reads T
        if(comm_IO_aggregation%myrank.eq.0) then
            open(myrank, FILE=trim(dir_cont_filein)//'cont_THETA_'//xyzrank//'.bin', FORM='UNFORMATTED', &
                STATUS='OLD', ACTION='READ', ACCESS='DIRECT', RECL=sizeof(var_chunk) )
            read(myrank, REC=1) var_chunk
            close(myrank)
        endif

        ! The master rank for data aggregation scatter the readed T to its local communicator using the defined DDT
        call MPI_Scatterv( var_chunk, cnts_aggr, disp_aggr, ddtype_aggregation_chunk, Tin, 1, ddtype_inner_domain, 0, &
                        comm_IO_aggregation%mpi_comm, ierr)

        ! The master rank for data aggregation deallocates the chunk array for aggregation
        if(comm_IO_aggregation%myrank.eq.0) then
            deallocate( var_chunk )
        endif

        ! Ghostcell update
        call mpi_subdomain_ghostcell_update(Uin)
        call mpi_subdomain_ghostcell_update(Vin)
        call mpi_subdomain_ghostcell_update(Win)
        call mpi_subdomain_ghostcell_update(Pin)
        call mpi_subdomain_ghostcell_update(Tin)

        if(myrank.eq.0) print '(a)', 'Read continue file using single process IO with aggregation'

    end subroutine mpi_Post_FileIn_Continue_Single_Process_IO_with_Aggregation

    !>
    !> @brief       Write field variables using single process IO with explicit aggregation
    !> @param       myrank          Rank for print
    !> @param       timein          Current dimensionless time for file write
    !> @param       dt              dt
    !> @param       Uin             Velocity field in x-direction for file write
    !> @param       Vin             Velocity field in y-direction for file write
    !> @param       Win             Velocity field in z-direction for file write
    !> @param       Pin             Pressure field for file write
    !> @param       Tin             Temperature field for file write
    !>
    subroutine mpi_Post_FileOut_Continue_Single_Process_IO_with_Aggregation(myrank,timein,dt,Uin,Vin,Win,Pin,Tin)

        use mpi

        implicit none

        integer :: myrank
        double precision :: timein,dt
        double precision, dimension(0:n1sub, 0:n2sub, 0:n3sub) :: Uin,Vin,Win,Pin,Tin

        double precision, allocatable, dimension(:,:,:) :: var_chunk
        integer :: ierr
        
        ! The master rank for data aggregation allocates chunk array which have aggregated data
        if(comm_IO_aggregation%myrank.eq.0) then
            allocate( var_chunk(n1msub, n2msub, n3msub * comm_IO_aggregation%nprocs) )
        endif

        ! The master rank for data aggregation writes time step information
        if(comm_IO_aggregation%myrank.eq.0) then
            open(myrank, FILE=trim(dir_cont_fileout)//'cont_time_'//xyzrank//'.bin', FORM='UNFORMATTED', &
                STATUS='REPLACE', ACTION='WRITE', ACCESS='DIRECT', RECL=sizeof(timein)*2 )
            write(myrank, REC=1) timein,dt
            close(myrank)
        endif

        ! The master rank for data aggregation gathers U to its chunk array
        call MPI_Gatherv( Uin, 1, ddtype_inner_domain, var_chunk, cnts_aggr, disp_aggr, ddtype_aggregation_chunk, 0, &
                        comm_IO_aggregation%mpi_comm, ierr)

        ! The master rank for data aggregation writes U
        if(comm_IO_aggregation%myrank.eq.0) then
            open(myrank, FILE=trim(dir_cont_fileout)//'cont_U_'//xyzrank//'.bin', FORM='UNFORMATTED', &
                STATUS='REPLACE', ACTION='WRITE', ACCESS='DIRECT', RECL=sizeof(var_chunk) )
            write(myrank, REC=1) var_chunk
            close(myrank)
        endif

        ! The master rank for data aggregation gathers V to its chunk array
        call MPI_Gatherv( Vin, 1, ddtype_inner_domain, var_chunk, cnts_aggr, disp_aggr, ddtype_aggregation_chunk, 0, &
                        comm_IO_aggregation%mpi_comm, ierr)

        ! The master rank for data aggregation writes V
        if(comm_IO_aggregation%myrank.eq.0) then
            open(myrank, FILE=trim(dir_cont_fileout)//'cont_V_'//xyzrank//'.bin', FORM='UNFORMATTED', &
                STATUS='REPLACE', ACTION='WRITE', ACCESS='DIRECT', RECL=sizeof(var_chunk) )
            write(myrank, REC=1) var_chunk
            close(myrank)
        endif

        ! The master rank for data aggregation gathers W to its chunk array
        call MPI_Gatherv( Win, 1, ddtype_inner_domain, var_chunk, cnts_aggr, disp_aggr, ddtype_aggregation_chunk, 0, &
                        comm_IO_aggregation%mpi_comm, ierr)

        ! The master rank for data aggregation writes W
        if(comm_IO_aggregation%myrank.eq.0) then
            open(myrank, FILE=trim(dir_cont_fileout)//'cont_W_'//xyzrank//'.bin', FORM='UNFORMATTED', &
                STATUS='REPLACE', ACTION='WRITE', ACCESS='DIRECT', RECL=sizeof(var_chunk) )
            write(myrank, REC=1) var_chunk
            close(myrank)
        endif

        ! The master rank for data aggregation gathers P to its chunk array
        call MPI_Gatherv( Pin, 1, ddtype_inner_domain, var_chunk, cnts_aggr, disp_aggr, ddtype_aggregation_chunk, 0, &
                        comm_IO_aggregation%mpi_comm, ierr)

        ! The master rank for data aggregation writes P
        if(comm_IO_aggregation%myrank.eq.0) then
            open(myrank, FILE=trim(dir_cont_fileout)//'cont_P_'//xyzrank//'.bin', FORM='UNFORMATTED', &
                STATUS='REPLACE', ACTION='WRITE', ACCESS='DIRECT', RECL=sizeof(var_chunk) )
            write(myrank, REC=1) var_chunk
            close(myrank)
        endif

        ! The master rank for data aggregation gathers T to its chunk array
        call MPI_Gatherv( Tin, 1, ddtype_inner_domain, var_chunk, cnts_aggr, disp_aggr, ddtype_aggregation_chunk, 0, &
                        comm_IO_aggregation%mpi_comm, ierr)

        ! The master rank for data aggregation writes T
        if(comm_IO_aggregation%myrank.eq.0) then
            open(myrank, FILE=trim(dir_cont_fileout)//'cont_THETA_'//xyzrank//'.bin', FORM='UNFORMATTED', &
                STATUS='REPLACE', ACTION='WRITE', ACCESS='DIRECT', RECL=sizeof(var_chunk) )
            write(myrank, REC=1) var_chunk
            close(myrank)
        endif

        ! The master rank for data aggregation deallocates the chunk array for aggregation
        if(comm_IO_aggregation%myrank.eq.0) then
            deallocate( var_chunk )
        endif

        if(myrank.eq.0) print '(a)', 'Write continue file using single process IO with aggregation'

    end subroutine mpi_Post_FileOut_Continue_Single_Process_IO_with_Aggregation

    !>
    !> @brief       Read field variables from a single file using MPI IO
    !> @param       myrank          Rank for print
    !> @param       timein          Current dimensionless time for file write
    !> @param       dt              dt
    !> @param       Uin             Velocity field in x-direction for file write
    !> @param       Vin             Velocity field in y-direction for file write
    !> @param       Win             Velocity field in z-direction for file write
    !> @param       Pin             Pressure field for file write
    !> @param       Tin             Temperature field for file write
    !>
    subroutine mpi_Post_FileIn_Continue_MPI_IO(myrank,timein,dt,Uin,Vin,Win,Pin,Tin)

        use mpi

        implicit none

        integer :: myrank
        double precision :: timein,dt
        double precision, dimension(0:n1sub, 0:n2sub, 0:n3sub) :: Uin,Vin,Win,Pin,Tin

        integer :: filep, ierr
        integer(kind=MPI_OFFSET_KIND) :: disp
        
        ! Rank 0 reads time step information and broadcast to every rank
        if(myrank.eq.0) then
            open(myrank, FILE=trim(dir_cont_filein)//'cont_time.bin', FORM='UNFORMATTED', &
                STATUS='OLD', ACTION='READ', ACCESS='DIRECT', RECL=sizeof(timein)*2 )
            read(myrank, REC=1) timein,dt
            close(myrank)
        endif
        call MPI_Bcast(timein, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(dt,     1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

        disp = 0

        ! Read U from a single file using MPI IO
        call MPI_File_open(MPI_COMM_WORLD, trim(dir_cont_filein)//'cont_U.bin', MPI_MODE_RDONLY, MPI_INFO_NULL, filep, ierr)
        call MPI_File_set_view(filep, disp, MPI_DOUBLE_PRECISION, ddtype_global_domain_MPI_IO, 'native', MPI_INFO_NULL, ierr)
        call MPI_File_read(filep, Uin, 1, ddtype_inner_domain, MPI_STATUS_IGNORE, ierr)
        call MPI_File_close(filep, ierr)

        ! Read V from a single file using MPI IO
        call MPI_File_open(MPI_COMM_WORLD, trim(dir_cont_filein)//'cont_V.bin', MPI_MODE_RDONLY, MPI_INFO_NULL, filep, ierr)
        call MPI_File_set_view(filep, disp, MPI_DOUBLE_PRECISION, ddtype_global_domain_MPI_IO, 'native', MPI_INFO_NULL, ierr)
        call MPI_File_read(filep, Vin, 1, ddtype_inner_domain, MPI_STATUS_IGNORE, ierr)
        call MPI_File_close(filep, ierr)

        ! Read W from a single file using MPI IO
        call MPI_File_open(MPI_COMM_WORLD, trim(dir_cont_filein)//'cont_W.bin', MPI_MODE_RDONLY, MPI_INFO_NULL, filep, ierr)
        call MPI_File_set_view(filep, disp, MPI_DOUBLE_PRECISION, ddtype_global_domain_MPI_IO, 'native', MPI_INFO_NULL, ierr)
        call MPI_File_read(filep, Win, 1, ddtype_inner_domain, MPI_STATUS_IGNORE, ierr)
        call MPI_File_close(filep, ierr)

        ! Read P from a single file using MPI IO
        call MPI_File_open(MPI_COMM_WORLD, trim(dir_cont_filein)//'cont_P.bin', MPI_MODE_RDONLY, MPI_INFO_NULL, filep, ierr)
        call MPI_File_set_view(filep, disp, MPI_DOUBLE_PRECISION, ddtype_global_domain_MPI_IO, 'native', MPI_INFO_NULL, ierr)
        call MPI_File_read(filep, Pin, 1, ddtype_inner_domain, MPI_STATUS_IGNORE, ierr)
        call MPI_File_close(filep, ierr)

        ! Read T from a single file using MPI IO
        call MPI_File_open(MPI_COMM_WORLD, trim(dir_cont_filein)//'cont_THETA.bin', MPI_MODE_RDONLY, MPI_INFO_NULL, filep, ierr)
        call MPI_File_set_view(filep, disp, MPI_DOUBLE_PRECISION, ddtype_global_domain_MPI_IO, 'native', MPI_INFO_NULL, ierr)
        call MPI_File_read(filep, Tin, 1, ddtype_inner_domain, MPI_STATUS_IGNORE, ierr)
        call MPI_File_close(filep, ierr)

        ! Ghostcell update
        call mpi_subdomain_ghostcell_update(Uin)
        call mpi_subdomain_ghostcell_update(Vin)
        call mpi_subdomain_ghostcell_update(Win)
        call mpi_subdomain_ghostcell_update(Pin)
        call mpi_subdomain_ghostcell_update(Tin)

        if(myrank.eq.0) print '(a)', 'Read continue file using MPI IO'

    end subroutine mpi_Post_FileIn_Continue_MPI_IO

    !>
    !> @brief       Write field variables to a single file using MPI IO
    !> @param       myrank          Rank for print
    !> @param       timein          Current dimensionless time for file write
    !> @param       dt              dt
    !> @param       Uin             Velocity field in x-direction for file write
    !> @param       Vin             Velocity field in y-direction for file write
    !> @param       Win             Velocity field in z-direction for file write
    !> @param       Pin             Pressure field for file write
    !> @param       Tin             Temperature field for file write
    !>
    subroutine mpi_Post_FileOut_Continue_MPI_IO(myrank,timein,dt,Uin,Vin,Win,Pin,Tin)

        use mpi

        implicit none

        integer :: myrank
        double precision :: timein,dt
        double precision, dimension(0:n1sub, 0:n2sub, 0:n3sub) :: Uin,Vin,Win,Pin,Tin

        integer :: filep, ierr
        integer(kind=MPI_OFFSET_KIND) :: disp

        ! Rank 0 writes time step information
        if(myrank.eq.0) then
            open(myrank, FILE=trim(dir_cont_fileout)//'cont_time.bin', FORM='UNFORMATTED', &
            STATUS='REPLACE', ACTION='WRITE', ACCESS='DIRECT', RECL=sizeof(timein)*2 )
            write(myrank, REC=1) timein,dt
            close(myrank)
        endif

        disp = 0

        ! Write U to a single file using MPI IO
        call MPI_File_open(MPI_COMM_WORLD, trim(dir_cont_fileout)//'cont_U.bin', MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL, filep, ierr)
        call MPI_File_set_view(filep, disp, MPI_DOUBLE_PRECISION, ddtype_global_domain_MPI_IO, 'native', MPI_INFO_NULL, ierr)
        call MPI_File_write(filep, Uin, 1, ddtype_inner_domain, MPI_STATUS_IGNORE, ierr)
        call MPI_File_close(filep, ierr)

        ! Write V to a single file using MPI IO
        call MPI_File_open(MPI_COMM_WORLD, trim(dir_cont_fileout)//'cont_V.bin', MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL, filep, ierr)
        call MPI_File_set_view(filep, disp, MPI_DOUBLE_PRECISION, ddtype_global_domain_MPI_IO, 'native', MPI_INFO_NULL, ierr)
        call MPI_File_write(filep, Vin, 1, ddtype_inner_domain, MPI_STATUS_IGNORE, ierr)
        call MPI_File_close(filep, ierr)

        ! Write W to a single file using MPI IO
        call MPI_File_open(MPI_COMM_WORLD, trim(dir_cont_fileout)//'cont_W.bin', MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL, filep, ierr)
        call MPI_File_set_view(filep, disp, MPI_DOUBLE_PRECISION, ddtype_global_domain_MPI_IO, 'native', MPI_INFO_NULL, ierr)
        call MPI_File_write(filep, Win, 1, ddtype_inner_domain, MPI_STATUS_IGNORE, ierr)
        call MPI_File_close(filep, ierr)

        ! Write P to a single file using MPI IO
        call MPI_File_open(MPI_COMM_WORLD, trim(dir_cont_fileout)//'cont_P.bin', MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL, filep, ierr)
        call MPI_File_set_view(filep, disp, MPI_DOUBLE_PRECISION, ddtype_global_domain_MPI_IO, 'native', MPI_INFO_NULL, ierr)
        call MPI_File_write(filep, Pin, 1, ddtype_inner_domain, MPI_STATUS_IGNORE, ierr)
        call MPI_File_close(filep, ierr)

        ! Write T to a single file using MPI IO
        call MPI_File_open(MPI_COMM_WORLD, trim(dir_cont_fileout)//'cont_THETA.bin', MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL, filep, ierr)
        call MPI_File_set_view(filep, disp, MPI_DOUBLE_PRECISION, ddtype_global_domain_MPI_IO, 'native', MPI_INFO_NULL, ierr)
        call MPI_File_write(filep, Tin, 1, ddtype_inner_domain, MPI_STATUS_IGNORE, ierr)
        call MPI_File_close(filep, ierr)

        if(myrank.eq.0) print '(a)', 'Write continue file using MPI IO'

    end subroutine mpi_Post_FileOut_Continue_MPI_IO

    !>
    !> @brief       Read field variables from a single file using MPI IO with explicit aggregation
    !> @param       myrank          Rank for print
    !> @param       timein          Current dimensionless time for file write
    !> @param       dt              dt
    !> @param       Uin             Velocity field in x-direction for file write
    !> @param       Vin             Velocity field in y-direction for file write
    !> @param       Win             Velocity field in z-direction for file write
    !> @param       Pin             Pressure field for file write
    !> @param       Tin             Temperature field for file write
    !>
    subroutine mpi_Post_FileIn_Continue_MPI_IO_with_Aggregation(myrank,timein,dt,Uin,Vin,Win,Pin,Tin)

        use mpi

        implicit none

        integer :: myrank
        double precision :: timein,dt
        double precision, dimension(0:n1sub, 0:n2sub, 0:n3sub) :: Uin,Vin,Win,Pin,Tin

        integer :: filep, ierr
        integer(kind=MPI_OFFSET_KIND) :: disp
        
        double precision, allocatable, dimension(:,:,:) :: var_chunk

        ! Rank 0 reads time step information and broadcast to every rank
        if(myrank.eq.0) then
            open(myrank, FILE=trim(dir_cont_filein)//'cont_time.bin', FORM='UNFORMATTED', &
                STATUS='OLD', ACTION='READ', ACCESS='DIRECT', RECL=sizeof(timein)*2 )
            read(myrank, REC=1) timein,dt
            close(myrank)
        endif
        call MPI_Bcast(timein, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(dt,     1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

        ! The master rank for data aggregation allocates chunk array which have aggregated data
        if(comm_IO_aggregation%myrank.eq.0) then
            allocate( var_chunk(n1msub, n2msub, n3msub * comm_IO_aggregation%nprocs) )
        endif

        disp = 0

        ! The master rank for data aggregation reads U from a single file using MPI IO
        if(comm_IO_aggregation%myrank.eq.0) then
            call MPI_File_open(comm_IO_master%mpi_comm, trim(dir_cont_filein)//'cont_U.bin', MPI_MODE_RDONLY, MPI_INFO_NULL, filep, ierr)
            call MPI_File_set_view(filep, disp, MPI_DOUBLE_PRECISION, ddtype_global_domain_MPI_IO_aggregation, 'native', MPI_INFO_NULL, ierr)
            call MPI_File_read(filep, var_chunk, size(var_chunk), MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
            call MPI_File_close(filep, ierr)
        endif

        ! The master rank for data aggregation scatter the readed U to its local communicator using the defined DDT
        call MPI_Scatterv( var_chunk, cnts_aggr, disp_aggr, ddtype_aggregation_chunk, Uin, 1, ddtype_inner_domain, 0, comm_IO_aggregation%mpi_comm, ierr)

        ! The master rank for data aggregation reads V from a single file using MPI IO
        if(comm_IO_aggregation%myrank.eq.0) then
            call MPI_File_open(comm_IO_master%mpi_comm, trim(dir_cont_filein)//'cont_V.bin', MPI_MODE_RDONLY, MPI_INFO_NULL, filep, ierr)
            call MPI_File_set_view(filep, disp, MPI_DOUBLE_PRECISION, ddtype_global_domain_MPI_IO_aggregation, 'native', MPI_INFO_NULL, ierr)
            call MPI_File_read(filep, var_chunk, size(var_chunk), MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
            call MPI_File_close(filep, ierr)
        endif

        ! The master rank for data aggregation scatter the readed V to its local communicator using the defined DDT
        call MPI_Scatterv( var_chunk, cnts_aggr, disp_aggr, ddtype_aggregation_chunk, Vin, 1, ddtype_inner_domain, 0, comm_IO_aggregation%mpi_comm, ierr)

        ! The master rank for data aggregation reads W from a single file using MPI IO
        if(comm_IO_aggregation%myrank.eq.0) then
            call MPI_File_open(comm_IO_master%mpi_comm, trim(dir_cont_filein)//'cont_W.bin', MPI_MODE_RDONLY, MPI_INFO_NULL, filep, ierr)
            call MPI_File_set_view(filep, disp, MPI_DOUBLE_PRECISION, ddtype_global_domain_MPI_IO_aggregation, 'native', MPI_INFO_NULL, ierr)
            call MPI_File_read(filep, var_chunk, size(var_chunk), MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
            call MPI_File_close(filep, ierr)
        endif

        ! The master rank for data aggregation scatter the readed W to its local communicator using the defined DDT
        call MPI_Scatterv( var_chunk, cnts_aggr, disp_aggr, ddtype_aggregation_chunk, Win, 1, ddtype_inner_domain, 0, comm_IO_aggregation%mpi_comm, ierr)

        ! The master rank for data aggregation reads P from a single file using MPI IO
        if(comm_IO_aggregation%myrank.eq.0) then
            call MPI_File_open(comm_IO_master%mpi_comm, trim(dir_cont_filein)//'cont_P.bin', MPI_MODE_RDONLY, MPI_INFO_NULL, filep, ierr)
            call MPI_File_set_view(filep, disp, MPI_DOUBLE_PRECISION, ddtype_global_domain_MPI_IO_aggregation, 'native', MPI_INFO_NULL, ierr)
            call MPI_File_read(filep, var_chunk, size(var_chunk), MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
            call MPI_File_close(filep, ierr)
        endif

        ! The master rank for data aggregation scatter the readed P to its local communicator using the defined DDT
        call MPI_Scatterv( var_chunk, cnts_aggr, disp_aggr, ddtype_aggregation_chunk, Pin, 1, ddtype_inner_domain, 0, comm_IO_aggregation%mpi_comm, ierr)

        ! The master rank for data aggregation reads T from a single file using MPI IO
        if(comm_IO_aggregation%myrank.eq.0) then
            call MPI_File_open(comm_IO_master%mpi_comm, trim(dir_cont_filein)//'cont_THETA.bin', MPI_MODE_RDONLY, MPI_INFO_NULL, filep, ierr)
            call MPI_File_set_view(filep, disp, MPI_DOUBLE_PRECISION, ddtype_global_domain_MPI_IO_aggregation, 'native', MPI_INFO_NULL, ierr)
            call MPI_File_read(filep, var_chunk, size(var_chunk), MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
            call MPI_File_close(filep, ierr)
        endif

        ! The master rank for data aggregation scatter the readed T to its local communicator using the defined DDT
        call MPI_Scatterv( var_chunk, cnts_aggr, disp_aggr, ddtype_aggregation_chunk, Tin, 1, ddtype_inner_domain, 0, comm_IO_aggregation%mpi_comm, ierr)

        if(comm_IO_aggregation%myrank.eq.0) then
            deallocate( var_chunk )
        endif

        ! Ghostcell update
        call mpi_subdomain_ghostcell_update(Uin)
        call mpi_subdomain_ghostcell_update(Vin)
        call mpi_subdomain_ghostcell_update(Win)
        call mpi_subdomain_ghostcell_update(Pin)
        call mpi_subdomain_ghostcell_update(Tin)
    
        if(myrank.eq.0) print '(a)', 'Read continue file using MPI IO with aggregation'

    end subroutine mpi_Post_FileIn_Continue_MPI_IO_with_Aggregation

    !>
    !> @brief       Write field variables to a single file using MPI IO with explicit aggregation
    !> @param       myrank          Rank for print
    !> @param       timein          Current dimensionless time for file write
    !> @param       dt              dt
    !> @param       Uin             Velocity field in x-direction for file write
    !> @param       Vin             Velocity field in y-direction for file write
    !> @param       Win             Velocity field in z-direction for file write
    !> @param       Pin             Pressure field for file write
    !> @param       Tin             Temperature field for file write
    !>
    subroutine mpi_Post_FileOut_Continue_MPI_IO_with_Aggregation(myrank,timein,dt,Uin,Vin,Win,Pin,Tin)

        use mpi

        implicit none

        integer :: myrank
        double precision :: timein,dt
        double precision, dimension(0:n1sub, 0:n2sub, 0:n3sub) :: Uin,Vin,Win,Pin,Tin

        integer :: filep, ierr
        integer(kind=MPI_OFFSET_KIND) :: disp
        
        double precision, allocatable, dimension(:,:,:) :: var_chunk

        ! Rank 0 writes time step information
        if(myrank.eq.0) then
            open(myrank, FILE=trim(dir_cont_fileout)//'cont_time.bin', FORM='UNFORMATTED', &
                STATUS='REPLACE', ACTION='WRITE', ACCESS='DIRECT', RECL=sizeof(timein)*2 )
            write(myrank, REC=1) timein,dt
            close(myrank)
        endif

        ! The master rank for data aggregation allocates chunk array which have aggregated data
        if(comm_IO_aggregation%myrank.eq.0) then
            allocate( var_chunk(n1msub, n2msub, n3msub * comm_IO_aggregation%nprocs) )
        endif

        disp = 0

        ! The master rank for data aggregation gathers U to its chunk array
        call MPI_Gatherv( Uin, 1, ddtype_inner_domain, var_chunk, cnts_aggr, disp_aggr, ddtype_aggregation_chunk, 0, &
                        comm_IO_aggregation%mpi_comm, ierr)

        ! The master rank for data aggregation writes U to a single file using MPI IO
        if(comm_IO_aggregation%myrank.eq.0) then
            call MPI_File_open(comm_IO_master%mpi_comm, trim(dir_cont_fileout)//'cont_U.bin', MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL, filep, ierr)
            call MPI_File_set_view(filep, disp, MPI_DOUBLE_PRECISION, ddtype_global_domain_MPI_IO_aggregation, 'native', MPI_INFO_NULL, ierr)
            call MPI_File_write(filep, var_chunk, size(var_chunk), MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
            call MPI_File_close(filep, ierr)
        endif

        ! The master rank for data aggregation gathers V to its chunk array
        call MPI_Gatherv( Vin, 1, ddtype_inner_domain, var_chunk, cnts_aggr, disp_aggr, ddtype_aggregation_chunk, 0, &
                        comm_IO_aggregation%mpi_comm, ierr)

        ! The master rank for data aggregation writes V to a single file using MPI IO
        if(comm_IO_aggregation%myrank.eq.0) then
            call MPI_File_open(comm_IO_master%mpi_comm, trim(dir_cont_fileout)//'cont_V.bin', MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL, filep, ierr)
            call MPI_File_set_view(filep, disp, MPI_DOUBLE_PRECISION, ddtype_global_domain_MPI_IO_aggregation, 'native', MPI_INFO_NULL, ierr)
            call MPI_File_write(filep, var_chunk, size(var_chunk), MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
            call MPI_File_close(filep, ierr)
        endif

        ! The master rank for data aggregation gathers W to its chunk array
        call MPI_Gatherv( Win, 1, ddtype_inner_domain, var_chunk, cnts_aggr, disp_aggr, ddtype_aggregation_chunk, 0, &
                        comm_IO_aggregation%mpi_comm, ierr)

        ! The master rank for data aggregation writes W to a single file using MPI IO
        if(comm_IO_aggregation%myrank.eq.0) then
            call MPI_File_open(comm_IO_master%mpi_comm, trim(dir_cont_fileout)//'cont_W.bin', MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL, filep, ierr)
            call MPI_File_set_view(filep, disp, MPI_DOUBLE_PRECISION, ddtype_global_domain_MPI_IO_aggregation, 'native', MPI_INFO_NULL, ierr)
            call MPI_File_write(filep, var_chunk, size(var_chunk), MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
            call MPI_File_close(filep, ierr)
        endif

        ! The master rank for data aggregation gathers P to its chunk array
        call MPI_Gatherv( Pin, 1, ddtype_inner_domain, var_chunk, cnts_aggr, disp_aggr, ddtype_aggregation_chunk, 0, &
                        comm_IO_aggregation%mpi_comm, ierr)

        ! The master rank for data aggregation writes P to a single file using MPI IO
        if(comm_IO_aggregation%myrank.eq.0) then
            call MPI_File_open(comm_IO_master%mpi_comm, trim(dir_cont_fileout)//'cont_P.bin', MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL, filep, ierr)
            call MPI_File_set_view(filep, disp, MPI_DOUBLE_PRECISION, ddtype_global_domain_MPI_IO_aggregation, 'native', MPI_INFO_NULL, ierr)
            call MPI_File_write(filep, var_chunk, size(var_chunk), MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
            call MPI_File_close(filep, ierr)
        endif

        ! The master rank for data aggregation gathers T to its chunk array
        call MPI_Gatherv( Tin, 1, ddtype_inner_domain, var_chunk, cnts_aggr, disp_aggr, ddtype_aggregation_chunk, 0, &
                        comm_IO_aggregation%mpi_comm, ierr)

        ! The master rank for data aggregation writes T to a single file using MPI IO
        if(comm_IO_aggregation%myrank.eq.0) then
            call MPI_File_open(comm_IO_master%mpi_comm, trim(dir_cont_fileout)//'cont_THETA.bin', MPI_MODE_WRONLY+MPI_MODE_CREATE, MPI_INFO_NULL, filep, ierr)
            call MPI_File_set_view(filep, disp, MPI_DOUBLE_PRECISION, ddtype_global_domain_MPI_IO_aggregation, 'native', MPI_INFO_NULL, ierr)
            call MPI_File_write(filep, var_chunk, size(var_chunk), MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr)
            call MPI_File_close(filep, ierr)
        endif

        ! The master rank for data aggregation deallocates the chunk array for aggregation
        if(comm_IO_aggregation%myrank.eq.0) then
            deallocate( var_chunk )
        endif

        if(myrank.eq.0) print '(a)', 'Write continue file using MPI IO with aggregation'

    end subroutine mpi_Post_FileOut_Continue_MPI_IO_with_Aggregation

    !>
    !> @brief       Read field variables using post-reassembly IO from a single file
    !> @param       myrank          Rank for print
    !> @param       timein          Current dimensionless time for file write
    !> @param       dt              dt
    !> @param       Uin             Velocity field in x-direction for file write
    !> @param       Vin             Velocity field in y-direction for file write
    !> @param       Win             Velocity field in z-direction for file write
    !> @param       Pin             Pressure field for file write
    !> @param       Tin             Temperature field for file write
    !>
    subroutine mpi_Post_FileIn_Continue_Post_Reassembly_IO(myrank,timein,dt,Uin,Vin,Win,Pin,Tin)

        use mpi

        implicit none

        integer :: myrank
        double precision :: timein,dt
        double precision, dimension(0:n1sub, 0:n2sub, 0:n3sub) :: Uin,Vin,Win,Pin,Tin

        double precision, allocatable, dimension(:,:,:) :: var_all
        integer :: ierr

        ! Rank 0 reads time step information and broadcast to every rank
        if(myrank.eq.0) then
            open(myrank, FILE=trim(dir_cont_filein)//'cont_time.bin', FORM='UNFORMATTED', &
                STATUS='OLD', ACTION='READ', ACCESS='DIRECT', RECL=sizeof(timein)*2 )
            read(myrank, REC=1) timein,dt
            close(myrank)
        endif
        call MPI_Bcast(timein, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
        call MPI_Bcast(dt,     1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

        ! Rank 0 allocates chunk array which have all field array data
        if(myrank.eq.0) allocate( var_all( n1m, n2m, n3m) )

        ! Rank 0 reads U from a single file
        if(myrank.eq.0) then
            open(myrank, FILE=trim(dir_cont_filein)//'cont_U.bin', FORM='UNFORMATTED', &
                STATUS='OLD', ACTION='READ', ACCESS='DIRECT', RECL=sizeof(var_all) )
            read(myrank, REC=1) var_all
            close(myrank)
        endif

        ! Rank 0 scatter U using the defined DDT
        call MPI_Scatterv(var_all, cnts_pr, disp_pr, ddtype_global_domain_pr_IO, Uin, 1, ddtype_inner_domain, 0, MPI_COMM_WORLD, ierr)

        ! Rank 0 reads V from a single file
        if(myrank.eq.0) then
            open(myrank, FILE=trim(dir_cont_filein)//'cont_V.bin', FORM='UNFORMATTED', &
                STATUS='OLD', ACTION='READ', ACCESS='DIRECT', RECL=sizeof(var_all) )
            read(myrank, REC=1) var_all
            close(myrank)
        endif

        ! Rank 0 scatter V using the defined DDT
        call MPI_Scatterv(var_all, cnts_pr, disp_pr, ddtype_global_domain_pr_IO, Vin, 1, ddtype_inner_domain, 0, MPI_COMM_WORLD, ierr)

        ! Rank 0 reads W from a single file
        if(myrank.eq.0) then
            open(myrank, FILE=trim(dir_cont_filein)//'cont_W.bin', FORM='UNFORMATTED', &
                STATUS='OLD', ACTION='READ', ACCESS='DIRECT', RECL=sizeof(var_all) )
            read(myrank, REC=1) var_all
            close(myrank)
        endif

        ! Rank 0 scatter W using the defined DDT
        call MPI_Scatterv(var_all, cnts_pr, disp_pr, ddtype_global_domain_pr_IO, Win, 1, ddtype_inner_domain, 0, MPI_COMM_WORLD, ierr)

        ! Rank 0 reads P from a single file
        if(myrank.eq.0) then
            open(myrank, FILE=trim(dir_cont_filein)//'cont_P.bin', FORM='UNFORMATTED', &
                STATUS='OLD', ACTION='READ', ACCESS='DIRECT', RECL=sizeof(var_all) )
            read(myrank, REC=1) var_all
            close(myrank)
        endif

        ! Rank 0 scatter P using the defined DDT
        call MPI_Scatterv(var_all, cnts_pr, disp_pr, ddtype_global_domain_pr_IO, Pin, 1, ddtype_inner_domain, 0, MPI_COMM_WORLD, ierr)

        ! Rank 0 reads T from a single file
        if(myrank.eq.0) then
            open(myrank, FILE=trim(dir_cont_filein)//'cont_THETA.bin', FORM='UNFORMATTED', &
                STATUS='OLD', ACTION='READ', ACCESS='DIRECT', RECL=sizeof(var_all) )
            read(myrank, REC=1) var_all
            close(myrank)
        endif

        ! Rank 0 scatter T using the defined DDT
        call MPI_Scatterv(var_all, cnts_pr, disp_pr, ddtype_global_domain_pr_IO, Tin, 1, ddtype_inner_domain, 0, MPI_COMM_WORLD, ierr)

        ! Rank 0 deallocates chunk array
        if(myrank.eq.0) deallocate( var_all )

        ! Ghostcell update
        call mpi_subdomain_ghostcell_update(Uin)
        call mpi_subdomain_ghostcell_update(Vin)
        call mpi_subdomain_ghostcell_update(Win)
        call mpi_subdomain_ghostcell_update(Pin)
        call mpi_subdomain_ghostcell_update(Tin)
    
        if(myrank.eq.0) print '(a)', 'Read continue file using post reassembly IO'

    end subroutine mpi_Post_FileIn_Continue_Post_Reassembly_IO

    !>
    !> @brief       Read field variables using post-reassembly IO from a single file
    !> @param       myrank          Rank for print
    !> @param       timein          Current dimensionless time for file write
    !> @param       dt              dt
    !> @param       Uin             Velocity field in x-direction for file write
    !> @param       Vin             Velocity field in y-direction for file write
    !> @param       Win             Velocity field in z-direction for file write
    !> @param       Pin             Pressure field for file write
    !> @param       Tin             Temperature field for file write
    !>
    subroutine mpi_Post_FileOut_Continue_Post_Reassembly_IO(myrank,timein,dt,Uin,Vin,Win,Pin,Tin)

        use mpi

        implicit none

        integer :: myrank
        double precision :: timein,dt
        double precision, dimension(0:n1sub, 0:n2sub, 0:n3sub) :: Uin,Vin,Win,Pin,Tin

        double precision, allocatable, dimension(:,:,:) :: var_all
        integer :: ierr

        ! Rank 0 writes time step information
        if(myrank.eq.0) then
            open(myrank, FILE=trim(dir_cont_fileout)//'cont_time.bin', FORM='UNFORMATTED', &
            STATUS='REPLACE', ACTION='WRITE', ACCESS='DIRECT', RECL=sizeof(timein)*2 )
            write(myrank, REC=1) timein,dt
            close(myrank)
        endif

        ! Rank 0 allocates chunk array which have all field array data
        if(myrank.eq.0) allocate( var_all( n1m, n2m, n3m) )

        ! Rank 0 gathers U using the defined DDT
        call MPI_Gatherv(Uin, 1, ddtype_inner_domain, var_all, cnts_pr, disp_pr, ddtype_global_domain_pr_IO, 0, MPI_COMM_WORLD, ierr)

        ! Rank 0 writes U to a single file
        if(myrank.eq.0) then
            open(myrank, FILE=trim(dir_cont_fileout)//'cont_U.bin', FORM='UNFORMATTED', &
                STATUS='REPLACE', ACTION='WRITE', ACCESS='DIRECT', RECL=sizeof(var_all) )
            write(myrank, REC=1) var_all
            close(myrank)
        endif

        ! Rank 0 gathers V using the defined DDT
        call MPI_Gatherv(Vin, 1, ddtype_inner_domain, var_all, cnts_pr, disp_pr, ddtype_global_domain_pr_IO, 0, MPI_COMM_WORLD, ierr)

        ! Rank 0 writes V to a single file
        if(myrank.eq.0) then
            open(myrank, FILE=trim(dir_cont_fileout)//'cont_V.bin', FORM='UNFORMATTED', &
                STATUS='REPLACE', ACTION='WRITE', ACCESS='DIRECT', RECL=sizeof(var_all) )
            write(myrank, REC=1) var_all
            close(myrank)
        endif

        ! Rank 0 gathers W using the defined DDT
        call MPI_Gatherv(Win, 1, ddtype_inner_domain, var_all, cnts_pr, disp_pr, ddtype_global_domain_pr_IO, 0, MPI_COMM_WORLD, ierr)

        ! Rank 0 writes W to a single file
        if(myrank.eq.0) then
            open(myrank, FILE=trim(dir_cont_fileout)//'cont_W.bin', FORM='UNFORMATTED', &
                STATUS='REPLACE', ACTION='WRITE', ACCESS='DIRECT', RECL=sizeof(var_all) )
            write(myrank, REC=1) var_all
            close(myrank)
        endif

        ! Rank 0 gathers P using the defined DDT
        call MPI_Gatherv(Pin, 1, ddtype_inner_domain, var_all, cnts_pr, disp_pr, ddtype_global_domain_pr_IO, 0, MPI_COMM_WORLD, ierr)

        ! Rank 0 writes P to a single file
        if(myrank.eq.0) then
            open(myrank, FILE=trim(dir_cont_fileout)//'cont_P.bin', FORM='UNFORMATTED', &
                STATUS='REPLACE', ACTION='WRITE', ACCESS='DIRECT', RECL=sizeof(var_all) )
            write(myrank, REC=1) var_all
            close(myrank)
        endif

        ! Rank 0 gathers T using the defined DDT
        call MPI_Gatherv(Tin, 1, ddtype_inner_domain, var_all, cnts_pr, disp_pr, ddtype_global_domain_pr_IO, 0, MPI_COMM_WORLD, ierr)

        ! Rank 0 writes T to a single file
        if(myrank.eq.0) then
            open(myrank, FILE=trim(dir_cont_fileout)//'cont_THETA.bin', FORM='UNFORMATTED', &
                STATUS='REPLACE', ACTION='WRITE', ACCESS='DIRECT', RECL=sizeof(var_all) )
            write(myrank, REC=1) var_all
            close(myrank)
        endif

        ! Rank 0 deallocates chunk array
        if(myrank.eq.0) deallocate( var_all )

        if(myrank.eq.0) print '(a)', 'Write continue file using post reassembly IO'

    end subroutine mpi_Post_FileOut_Continue_Post_Reassembly_IO


    subroutine mpi_Post_WallShearStress(myrank,U,V,W,TH,outTimeStep,wss)
        
        use MPI
        use mpi_topology,   only : comm_1d_x1, comm_1d_x2, comm_1d_x3
        use mpi_subdomain,  only : dx2_sub
        use mpi_subdomain,  only : n1sub, n2sub, n3sub

        implicit none
        integer :: outTimeStep,myrank
        double precision, dimension(0:n1sub,0:n2sub,0:n3sub) :: U,V,W,TH
        double precision, dimension(:), pointer :: dx2
        double precision :: dUdy_loc,dUdy_I,dUdy_K,wss,Tss
        double precision :: dTdy_loc,dTdy_I,dTdy_K
        character(len=23) :: filename_wallshearstress

        integer :: i,j,k
        integer :: im,jm,km
        integer :: ip,jp,kp
        integer :: Sc,Ec
        integer :: Sm,Em
        integer :: Sp,Ep
        integer :: ierr

        dx2  => dx2_sub

        dUdy_loc=0.d0; dUdy_I=0.d0; dUdy_K=0.d0
        dTdy_loc=0.d0; dTdy_I=0.d0; dTdy_K=0.d0
        wss=0.d0     ; Tss=0.d0

        if(comm_1d_x2%myrank==0) then
            do k=1,n3sub
            do i=1,n1sub
                dudy_loc=dudy_loc+( U(i,2,k)- U(i,1,k))/(dx2(1)*0.5d0)
                dTdy_loc=dTdy_loc+(TH(i,2,k)-TH(i,1,k))/(dx2(1)*0.5d0)
            enddo
            enddo
        endif
        if(comm_1d_x2%myrank==comm_1d_x2%nprocs-1) then
            do k=1,n3sub
            do i=1,n1sub
                dudy_loc=dudy_loc+( U(i,n2sub-1,k)- U(i,n2sub,k))/(dx2_sub(n2sub-1)*0.5d0)
                dTdy_loc=dTdy_loc+(TH(i,      2,k)-TH(i,    1,k))/(dx2_sub(1      )*0.5d0)
            enddo
            enddo
        endif

        dUdy_loc=dUdy_loc/real(n1sub*n3sub,8)
        dTdy_loc=dTdy_loc/real(n1sub*n3sub,8)

        call MPI_ALLREDUCE(dUdy_loc, dUdy_I, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_1d_x1%mpi_comm, ierr)
        call MPI_ALLREDUCE(dUdy_I  , dUdy_K, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_1d_x3%mpi_comm, ierr)
        call MPI_ALLREDUCE(dUdy_K  , wss   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_1d_x2%mpi_comm, ierr)

        call MPI_ALLREDUCE(dTdy_loc, dTdy_I, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_1d_x1%mpi_comm, ierr)
        call MPI_ALLREDUCE(dTdy_I  , dTdy_K, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_1d_x3%mpi_comm, ierr)
        call MPI_ALLREDUCE(dTdy_K  , Tss   , 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_1d_x2%mpi_comm, ierr)

        !wss=wss/real(np_dim(0)*np_dim(2),8)/2.d0
        !Tss=Tss/real(np_dim(0)*np_dim(2),8)/2.d0

        wss=wss/real(np1*np3,8)/2.d0
        Tss=Tss/real(np1*np3,8)/2.d0

        filename_wallshearstress='Monitor_WallShearStress'

        if(myrank==0) then
            open(unit=myrank,file=dir_instantfield//filename_wallshearstress//'.plt',position='append')
                if(outTimeStep==1) write(myrank,*) 'VARIABLES="Timestep","wss","Tss"'
                write(myrank,*) outTimeStep, wss, Tss
            close(myrank)
        endif

    end subroutine mpi_Post_WallShearStress

end module mpi_Post