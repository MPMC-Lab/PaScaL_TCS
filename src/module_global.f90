!======================================================================================================================
!> @file        module_global.f90
!> @brief       This file contains a module of global input parameters for the example problem of PaScaL_TCS.
!> @details     The input parameters include global domain information, boundary conditions, fluid properties, 
!>              flow conditions, and simulation control parameters.
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
!> @brief       Module for global parameters.
!> @details     This global module has simulation parameters and a subroutine to initialize the parameters. 
!>
module global
    implicit none
    double precision, parameter :: PI = acos(-1.d0)

    ! Physical parameters
    double precision :: Pr, Ra, Re_c, Re_u, Gr
    double precision :: Cmu, Cmt, Ct, Cmp

    ! Computational size for the physical domain and time discretization
    integer :: n1,n2,n3
    integer :: n1m,n2m,n3m
    integer :: n1p,n2p,n3p
    integer :: np1,np2,np3

    logical :: pbc1, pbc2, pbc3
    integer :: UNIFORM1,UNIFORM2,UNIFORM3
    double precision :: GAMMA1, GAMMA2, GAMMA3

    integer :: Timestepmax, problem
    double precision :: time,dt
    double precision :: CFL,MaxCFL, dtStart, tStart, dtMax

    ! Physical size of the computational domain
    double precision :: L1,L2,L3, Volume
    double precision :: H,Aspect1,Aspect3
    double precision :: x1_start,x2_start,x3_start
    double precision :: x1_end,x2_end,x3_end

    double precision :: Uup,Vup,Wup,Tup
    double precision :: Ubt,Vbt,Wbt,Tbt

    double precision :: DeltaT
    double precision :: a10, a11, a12, a13, a14, a15, a12pera11
    double precision :: b10, b11, b12, b13, b14, b15
    double precision :: c10, c11, c12, c13, c14, c15
    double precision :: d10, d11, d12, d13, d14, d15
    double precision :: KPP0, Rho0, Cp0, Mu0

    integer  :: les_model

    character(len=128) dir_cont_fileout
    character(len=128) dir_instantfield
    character(len=128) dir_cont_filein
    character(len=128) dir_statistics ! Statistics from YMG

    integer :: print_start_step,print_interval_step,print_j_index_wall,print_j_index_bulk,print_k_index
    integer :: ContinueFilein
    integer :: ContinueFileout

    ! CJY : Parameters for channel flow
    double precision :: presgrad1, wss, vper !YMG : adding vper for initial channel random number
    
    contains
    !>
    !> @brief       Assign global parameters.
    !> @param       np_dim      Number of MPI processes in 3D topology
    !>

    subroutine global_inputpara()

        implicit none
        double precision :: tttmp(1:3)
        character(len=512) :: temp_char
        integer :: i

        ! Namelist variables for file input
        namelist /meshes/               n1m, n2m, n3m
        namelist /MPI_procs/            np1, np2, np3
        namelist /periodic_boundary/    pbc1, pbc2, pbc3
        namelist /uniform_mesh/         uniform1, uniform2, uniform3
        namelist /mesh_stretch/         gamma1, gamma2, gamma3
        namelist /aspect_ratio/         Aspect1, H, Aspect3
        namelist /flow_style/           problem
        namelist /setting_RBC/          Ra, Pr
        namelist /setting_channel/      Re_c, Pr
        namelist /setting_urban/        Re_u, Pr
        namelist /les/                  les_model
        namelist /sim_parameter/        DeltaT, MaxCFL, vper
        namelist /sim_control/          dtStart, tStart, Timestepmax, print_start_step, print_interval_step
        namelist /sim_continue/         ContinueFilein, ContinueFileout, dir_cont_filein, dir_cont_fileout, dir_instantfield

        ! Using file input
        open(unit = 1, file = "PARA_INPUT.dat")
        read(1, sim_continue)
        read(1, meshes)
        read(1, MPI_procs)
        read(1, periodic_boundary)
        read(1, uniform_mesh)
        read(1, mesh_stretch)
        read(1, aspect_ratio)
        read(1, flow_style)
        read(1, setting_RBC)
        read(1, setting_channel)
        read(1, setting_urban)
        read(1, les)
        read(1, sim_parameter)
        read(1, sim_control)
        close(1)


        Uup = 0.d0;Ubt = 0.d0
        Vup = 0.d0;Vbt = 0.d0
        Wup = 0.d0;Wbt = 0.d0
        Tup = -0.5d0;Tbt =0.5d0
        
        ! Coefficients for Water
        ! a10 = 0.9922d+3; a11 = -3.736d-4;  a12 = -3.98d-6;  a13 = 0.;         a14 = 0.;  a15 = 0.
        ! b10 = 4.1690d+3; b11 =  0.084d-4;  b12 =  4.60d-6;  b13 = 0.;         b14 = 0.;  b15 = 0.
        ! c10 = 0.6297;    c11 = 21.99d-4;   c12 = -17.8d-6;  c13 = 0.;         c14 = 0.;  c15 = 0. 
        ! d10 = 0.6690d-6; d11 = -175.9d-4;  d12 = 295.8d-6;  d13 = -460.d-8;   d14 = 0.;  d15 = 0.
        a10 = 1.d0; a11 =0.d0; a12 = 0.d0;  a13 = 0.d0; a14 = 0.;  a15 = 0.
        b10 = 1.d0; b11 =0.d0; b12 = 0.d0;  b13 = 0.d0; b14 = 0.;  b15 = 0.
        c10 = 1.d0; c11 =0.d0; c12 = 0.d0;  c13 = 0.d0; c14 = 0.;  c15 = 0. 
        d10 = 1.d0; d11 =0.d0; d12 = 0.d0;  d13 = 0.d0; d14 = 0.;  d15 = 0.
        
        if(a12<1.0D-14) then
            a12pera11=0.d0
        else
            a12pera11=a12/a11
        endif
        
        ! Coefficients for Glycerol
        ! a10 = 1.2477d+3; a11 = -4.789d-4;  a12 = -0.3795d-6;  a13 = 0.;         a14 = 0.;  a15 = 0.
        ! b10 = 2.5108d+3; b11 = 22.511d-4;  b12 =  0.       ;  b13 = 0.;         b14 = 0.;  b15 = 0.
        ! c10 = 2.9351d-3; c11 = 3.863d-4;   c12 =  0.       ;  c13 = 0.;         c14 = 0.;  c15 = 0. 
        ! d10 = 238.71d-6; d11 =-702.83d-4; d12 = 2393.1d-6;    d13 = -6923.0d-8; d14 = 33131.3d-10; d15 = -71517.5d-12

        Rho0  = a10
        Cp0   = b10
        KPP0  = c10
        Mu0   = a10*d10            

        if(problem.eq.0) then ! RBC
            Cmu = sqrt(Pr/Ra);    Cmt = 1.0d0;   Ct  = 1.0d0/sqrt(Ra*Pr)
            Gr  = Ra/Pr
            Cmp = 1.0d0
        else if(problem.eq.1) then ! Channel
            Cmu = 1.0d0/Re_c; Cmt = 0.0d0; Ct  = 1.0d0/(Re_c*Pr)
            Gr  = Ra/Pr
            Cmp = 1.0d0
        else if(problem.eq.2) then  ! Urban
            Cmu = 1.0d0/Re_u; Cmt = 0.0d0; Ct  = 1.0d0/(Re_u*Pr)
            Gr  = Ra/Pr
            Cmp = 1.0d0
        endif

        ! Computational size for the physical domain and time discretization
        n1=n1m+1;n1p=n1+1;
        n2=n2m+1;n2p=n2+1;
        n3=n3m+1;n3p=n3+1;

        ! Physical size of the computational domain
        L1=H*Aspect1
        L2=2.0d0*H
        L3=H*Aspect3

        Volume = L1*L2*L3

        x1_start = 0.0d0
        x2_start = 0.0d0
        x3_start = 0.0d0

        x1_end=x1_start+L1
        x2_end=x2_start+L2
        x3_end=x3_start+L3
        
        tttmp(1)=L1/dble(n1-1)
        tttmp(2)=L2/dble(n2-1)
        tttmp(3)=L3/dble(n3-1)

        dtMax=minval(tttmp)*100
        
    end subroutine global_inputpara

end module global