!!======================================================================================================================
!> @file        main.f90
!> @brief       This file contains the main subroutines for an example problem of PaScaL_TCS.
!> @details     The target example problem is a three-dimensional(3D) Rayleigh-Bernard natural convection problem
!>              in a channel domain with the boundary conditions of vertically different wall temperature and 
!>              horizontally periodic boundaries. The non-Oberbeck–Boussinesq effect is considered. 
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
!> @brief       Main execution program for the example problem of PaScaL_TCS.
!>

program main
  
    use mpi
    use mpi_topology
    use global
    use mpi_subdomain    
    use mpi_thermal
    use mpi_momentum
    use mpi_pressure

    use mpi_Post
    use timer

    implicit none

    integer :: i,j,k,TimeStep
    integer :: ierr
    double precision :: maxDivU

    character(len=64)   :: timer_str(64)
    double precision    :: timer_a,timer_b
    integer, parameter  :: stamp_main = 1
    
    call MPI_Init(ierr)
    call MPI_Comm_size( MPI_COMM_WORLD, nprocs, ierr)
    call MPI_Comm_rank( MPI_COMM_WORLD, myrank, ierr)

    if(myrank==0) write(*,*) '[Main] The simulation starts!'

    ! Timer string setup. 
    timer_str(1)  = '[Main] Thermal coeff             '
    timer_str(2)  = '[Main] Thermal solve             '
    timer_str(3)  = '[Main] Momentum coeff            '
    timer_str(4)  = '[Main] Momentum solvedU          '
    timer_str(5)  = '[Main] Momentum solvedV          '
    timer_str(6)  = '[Main] Momentum solvedW          '
    timer_str(7)  = '[Main] Momentum blockLdV         '
    timer_str(8)  = '[Main] Momentum blockLdU         '
    timer_str(9)  = '[Main] Momentum updatePseudoUVW  '
    timer_str(10) = '[Main] Pressure RHS              '
    timer_str(11) = '[Main] Pressure PoissonSolve     '
    timer_str(12) = '[Main] Pressure projection       '
    timer_str(13) = '[Main] Update GhostCell          '
    timer_str(14) = '[Main] Clean coeff               '
    timer_str(15) = '[Main] Post div                  '
    timer_str(16) = '[Main] Post CFL                  '
    timer_str(17) = '[Main] Misc operation            '
    timer_str(18) = '[Main] File write                '

    call timer_init(18,timer_str)
    if(myrank==0) write(*,*) '[Main] Timer initialized!'

    if(myrank==0) call system('mkdir -p ./data')
    if(myrank==0) call system('mkdir -p ./data/1_continue')
    if(myrank==0) call system('mkdir -p ./data/2_instanfield')
    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    
    call global_inputpara()
    if(myrank==0) write(*,*) '[Main] Read input parameters!'

    call mpi_topology_make()
    call mpi_subdomain_make()
    call mpi_subdomain_mesh()
    call mpi_subdomain_indices()

    call mpi_subdomain_DDT_ghostcell()

    call mpi_thermal_allocation()
    call mpi_thermal_initial()
    call mpi_thermal_boundary()
    
    call mpi_momentum_allocation()
    call mpi_momentum_initial()
    call mpi_momentum_boundary()

    call mpi_pressure_allocation()
    call mpi_pressure_initial()
    call mpi_subdomain_DDT_transpose2()
    call mpi_pressure_wave_number()

    call mpi_subdomain_ghostcell_update(T)
    call mpi_subdomain_ghostcell_update(U)
    call mpi_subdomain_ghostcell_update(V)
    call mpi_subdomain_ghostcell_update(W)
    call mpi_subdomain_ghostcell_update(P)

    call MPI_Barrier(MPI_COMM_WORLD,ierr)

    call mpi_Post_allocation(1)

    if(myrank==0) write(*,*) '[Main] Simulation setup completed!'

    CFL = MaxCFL
    time = tStart 
    dt = dtStart

    if ( ContinueFilein==1 ) then
        call mpi_Post_FileIn_Continue_Post_Reassembly_IO(myrank,tStart,dt,U,V,W,P,T)
        call mpi_Post_CFL(U,V,W,CFL,dt)
        time = tStart
    end if

    if(myrank==0) write(*,*) '[Main] Iteration starts!'

    do TimeStep = 1, Timestepmax
        timer_a=MPI_WTIME()
        call timer_stamp0(stamp_main)

        call mpi_thermal_coeffi()
        call timer_stamp(1, stamp_main)

        call mpi_thermal_solver(U,V,W)
        call timer_stamp(2, stamp_main)

        call mpi_thermal_coeffi_clean()
        call timer_stamp(14, stamp_main)

        call mpi_subdomain_ghostcell_update(T)
        call timer_stamp(13, stamp_main)
        
        call mpi_momentum_coeffi(T)
        call timer_stamp(3, stamp_main)

        call mpi_momentum_solvedU(T)
        call timer_stamp(4, stamp_main)

        call mpi_subdomain_ghostcell_update(dU)
        call timer_stamp(13, stamp_main)
        
        call mpi_momentum_solvedV(T)
        call timer_stamp(5, stamp_main)

        call mpi_subdomain_ghostcell_update(dV)
        call timer_stamp(13, stamp_main)

        call mpi_momentum_solvedW(T)
        call timer_stamp(6, stamp_main)

        call mpi_subdomain_ghostcell_update(dW)
        call timer_stamp(13, stamp_main)

        call mpi_momentum_blockLdV(T)
        call timer_stamp(7, stamp_main)

        call mpi_subdomain_ghostcell_update(dV)        
        call timer_stamp(13, stamp_main)

        call mpi_momentum_blockLdU(T)
        call timer_stamp(8, stamp_main)

        call mpi_momentum_pseudoupdateUVW()
        call timer_stamp(9, stamp_main)

        call mpi_subdomain_ghostcell_update(U)
        call mpi_subdomain_ghostcell_update(V)
        call mpi_subdomain_ghostcell_update(W)
        call timer_stamp(13, stamp_main)

        call mpi_pressure_RHS(invRho,U,V,W,VBCup_sub,VBCbt_sub)
        call timer_stamp(10, stamp_main)

        call mpi_pressure_Poisson_FFT2(dx2_sub,dmx2_sub)
        call timer_stamp(11, stamp_main)

        call mpi_subdomain_ghostcell_update(dP)
        call timer_stamp(13, stamp_main)

        call mpi_pressure_Projection(invRho,U,V,W,P,comm_1d_x2%myrank,comm_1d_x2%nprocs) ! dealloc :: dP
        call timer_stamp(12, stamp_main)

        call mpi_subdomain_ghostcell_update(U)
        call mpi_subdomain_ghostcell_update(V)
        call mpi_subdomain_ghostcell_update(W)
        call mpi_subdomain_ghostcell_update(P)
        call mpi_subdomain_ghostcell_update(dPhat)
        call timer_stamp(13, stamp_main)
        
        call mpi_momentum_coeffi_clean()
        call timer_stamp(14, stamp_main)

        call mpi_Post_Div(U,V,W,maxDivU)
        call timer_stamp(15, stamp_main)
        
        if(maxDivU>=1.d-3) then
            call mpi_Post_FileOut_InstanField(myrank,TimeStep,time,U,V,W,P,T)
            if(myrank==0) write(*,*) 'BLOW UP',TimeStep
            exit
        endif

        timer_b=MPI_WTIME()
        call mpi_Post_MonitorOut(myrank,TimeStep,time,dt,CFL,maxDivU,timer_b-timer_a)
        call timer_stamp(17, stamp_main)
        
        call mpi_Post_CFL(U,V,W,CFL,dt)
        call timer_stamp(16, stamp_main)

        time=time+dt

    enddo
    if(myrank==0) write(*,*) '[Main] Iteration ends!'

    if ( ContinueFileout==1 ) then
        call mpi_Post_FileOut_Continue_Post_Reassembly_IO(myrank,time,dt,U,V,W,P,T)
    endif
    call timer_stamp(18, stamp_main)

    call timer_reduction()
    call timer_output(myrank, nprocs)

    call mpi_pressure_clean()
    call mpi_momentum_clean()    
    call mpi_thermal_clean()

    call mpi_subdomain_indices_clean()
    call mpi_subdomain_clean()

    call mpi_topology_clean()
    call MPI_FINALIZE(ierr)

    if(myrank==0) write(*,*) '[Main] The main simulation complete! '
end
