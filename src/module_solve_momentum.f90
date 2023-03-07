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
!    double precision, allocatable, dimension(:,:) :: UBCup_sub,VBCup_sub,WBCup_sub
 !   double precision, allocatable, dimension(:,:) :: UBCbt_sub,VBCbt_sub,WBCbt_sub
    double precision, allocatable, dimension(:,:,:,:) :: XMBC, YMBC, ZMBC ! 3: 1,2,3(U,V,W) 4: 1,2(bt,up)
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
        integer :: i, j, k

        allocate(U(0:n1sub,0:n2sub,0:n3sub),V(0:n1sub,0:n2sub,0:n3sub),W(0:n1sub,0:n2sub,0:n3sub))
        allocate(P(0:n1sub,0:n2sub,0:n3sub))

        ! allocate( UBCup_sub(0:n1sub,0:n3sub),VBCup_sub(0:n1sub,0:n3sub),WBCup_sub(0:n1sub,0:n3sub))
        ! allocate( UBCbt_sub(0:n1sub,0:n3sub),VBCbt_sub(0:n1sub,0:n3sub),WBCbt_sub(0:n1sub,0:n3sub))

        ! do k = 0, n3sub
        ! do i = 0, n1sub
        !     UBCup_sub(i,k) = 0.d0
        !     UBCbt_sub(i,k) = 0.d0 
        !     VBCup_sub(i,k) = 0.d0
        !     VBCbt_sub(i,k) = 0.d0 
        !     WBCup_sub(i,k) = 0.d0
        !     WBCbt_sub(i,k) = 0.d0
        ! enddo
        ! enddo

        allocate(XMBC(0:n2sub,0:n3sub, 1:3, 1:2), YMBC(0:n1sub,0:n3sub, 1:3, 1:2), ZMBC(0:n1sub,0:n2sub, 1:3, 1:2))
        XMBC=0.0d0
        YMBC=0.0d0
        ZMBC=0.0d0

    end subroutine mpi_momentum_allocation

    !>
    !> @brief       Deallocate variables for momentum solver
    !>
    subroutine mpi_momentum_clean()

        implicit none
        
        deallocate(U,V,W)
        deallocate(P)

       ! deallocate( UBCup_sub,VBCup_sub,WBCup_sub)
       ! deallocate( UBCbt_sub,VBCbt_sub,WBCbt_sub)
        deallocate(XMBC,YMBC,ZMBC)
        
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
    !> @brief       Initialize velocity field for channel flow
    !>
    subroutine mpi_momentum_initial_channel() !cell
        use MPI
        use mpi_topology,   only : myrank, comm_1d_x1, comm_1d_x2, comm_1d_x3
        use mpi_subdomain,  only : x1_sub, x2_sub, x3_sub, dx1_sub, dx2_sub, dx3_sub
        use mpi_subdomain,  only : n1sub, n2sub, n3sub, n1msub, n2msub, n3msub

        implicit none
    
        double precision, dimension(:), pointer ::  x1, x2, x3
        double precision, dimension(:), pointer :: dx1, dx2, dx3
        !integer, dimension(-1:n1sub+1,-1:n2sub+1,-1:n3sub+1), intent(in) :: cell
        double precision, allocatable, dimension(:,:,:) :: RanNum1, RanNum2, RanNum3
        
        integer :: i,j,k
        integer :: kp,km  
        integer :: jp,jm
        integer :: ip,im
        integer :: ierr
        !integer :: cellijk, cellip, celljp, cellkp, cellijp, cellikp, celljkp
    
        double precision:: yh, Umx
        double precision:: Um, Um_I, Um_K, Um_total
        double precision:: Vm, Vm_I, Vm_K, Vm_total
        double precision:: Wm, Wm_I, Wm_K, Wm_total
    
        allocate(RanNum1(0:n1sub,0:n2sub,0:n3sub))
        allocate(RanNum2(0:n1sub,0:n2sub,0:n3sub))
        allocate(RanNum3(0:n1sub,0:n2sub,0:n3sub))

        ! Pointer of grid information
        x1 => x1_sub
        x2 => x2_sub
        x3 => x3_sub
        dx1 => dx1_sub
        dx2 => dx2_sub
        dx3 => dx3_sub

        U(0:n1sub, 0:n2sub, 0:n3sub)=0.d0
        V(0:n1sub, 0:n2sub, 0:n3sub)=0.d0
        W(0:n1sub, 0:n2sub, 0:n3sub)=0.d0
        P(0:n1sub, 0:n2sub, 0:n3sub)=0.d0
   
        RanNum1(0:n1sub,0:n2sub,0:n3sub)=0.d0
        RanNum2(0:n1sub,0:n2sub,0:n3sub)=0.d0
        RanNum3(0:n1sub,0:n2sub,0:n3sub)=0.d0

        Umx  = 1.5d0
        call random_seed() 
        call random_number(RanNum1)
        call random_seed() 
        call random_number(RanNum2)
        call random_seed() 
        call random_number(RanNum3)

        ! Eliminate mean quantities
        Um=0.d0; Vm=0.d0; Wm=0.d0

        RanNum1(1:n1msub, 1:n2msub, 1:n3msub) = RanNum1(1:n1msub, 1:n2msub, 1:n3msub) - 0.5d0
        RanNum2(1:n1msub, 1:n2msub, 1:n3msub) = RanNum2(1:n1msub, 1:n2msub, 1:n3msub) - 0.5d0
        RanNum3(1:n1msub, 1:n2msub, 1:n3msub) = RanNum3(1:n1msub, 1:n2msub, 1:n3msub) - 0.5d0

        do k=1,n3msub
        do j=1,n2msub
        do i=1,n1msub       
            Um = Um + RanNum1(i,j,k) * dx1(i) * dx2(j) * dx3(k)
            Vm = Vm + RanNum2(i,j,k) * dx1(i) * dx2(j) * dx3(k)
            Wm = Wm + RanNum3(i,j,k) * dx1(i) * dx2(j) * dx3(k)
        enddo
        enddo
        enddo
        
        Um = Um / (x1(n1sub)-x1(1)); Um = Um / (x2(n2sub)-x2(1)); Um = Um / (x3(n3sub)-x3(1)) 
        Vm = Vm / (x1(n1sub)-x1(1)); Vm = Vm / (x2(n2sub)-x2(1)); Vm = Vm / (x3(n3sub)-x3(1)) 
        Wm = Wm / (x1(n1sub)-x1(1)); Wm = Wm / (x2(n2sub)-x2(1)); Wm = Wm / (x3(n3sub)-x3(1))

        RanNum1(1:n1msub, 1:n2msub, 1:n3msub) = RanNum1(1:n1msub, 1:n2msub, 1:n3msub) - Um
        RanNum2(1:n1msub, 1:n2msub, 1:n3msub) = RanNum2(1:n1msub, 1:n2msub, 1:n3msub) - Vm
        RanNum3(1:n1msub, 1:n2msub, 1:n3msub) = RanNum3(1:n1msub, 1:n2msub, 1:n3msub) - Wm 
    
        ! For Test
        Um=0.d0; Vm=0.d0; Wm=0.d0

        do k=1,n3msub
        do j=1,n2msub
        do i=1,n1msub       
            Um = Um + RanNum1(i,j,k) * dx1(i) * dx2(j) * dx3(k)
            Vm = Vm + RanNum2(i,j,k) * dx1(i) * dx2(j) * dx3(k)
            Wm = Wm + RanNum3(i,j,k) * dx1(i) * dx2(j) * dx3(k)
        enddo
        enddo
        enddo
    
        Um = Um / (x1(n1sub)-x1(1)); Um = Um / (x2(n2sub)-x2(1)); Um = Um / (x3(n3sub)-x3(1)) 
        Vm = Vm / (x1(n1sub)-x1(1)); Vm = Vm / (x2(n2sub)-x2(1)); Vm = Vm / (x3(n3sub)-x3(1)) 
        Wm = Wm / (x1(n1sub)-x1(1)); Wm = Wm / (x2(n2sub)-x2(1)); Wm = Wm / (x3(n3sub)-x3(1))
    
        Um=0.d0; Um_I=0.d0;  Um_K=0.d0;  Um_total=0.d0
        Vm=0.d0; Vm_I=0.d0;  Vm_K=0.d0;  Vm_total=0.d0
        Wm=0.d0; Wm_I=0.d0;  Wm_K=0.d0;  Wm_total=0.d0

        do k=1,n3msub
        do j=1,n2msub
        do i=1,n1msub
            yh = 0.5d0*( x2(j)+x2(j+1) )
            ! For channel flow
            U(i,j,k) = Umx*(1.d0-(yh/H-1.d0)**2.d0) + vper*Umx*RanNum1(i,j,k)
            ! For uniform
            ! U(i,j,k) = 1.d0 + vper*Umx*RanNum1(i,j,k)
            V(i,j,k) =        vper*Umx*RanNum2(i,j,k)
            W(i,j,k) =        vper*Umx*RanNum3(i,j,k)
    
            Um = Um + U(i,j,k) * dx1(i) * dx2(j) * dx3(k)
            Vm = Vm + V(i,j,k) * dx1(i) * dx2(j) * dx3(k)
            Wm = Wm + W(i,j,k) * dx1(i) * dx2(j) * dx3(k)
        enddo
        enddo
        enddo

        !Cell classification
            ! do k=1,n3msub
            !     kp = k + 1
            ! do j=1,n2msub
            !     jp = j + 1
            ! do i=1,n1msub
            !     ip = i + 1
            !     cellijk=cell(i ,j ,k )
            !     cellip =cell(ip,j ,k )
            !     celljp =cell(i ,jp,k )
            !     cellkp =cell(i ,j ,kp)
            !     cellijp=cell(ip,jp,k )
            !     cellikp=cell(ip,j ,kp)
            !     celljkp=cell(i ,jp,kp)
        
            !     if( MOD((cellijk+celljp+cellkp+celljkp),10) .eq. 0 ) then
            !         Um = Um - u(i,j,k) * dx1(i) * dx2(j) * dx3(k)
            !         U(i,j,k)=0.d0
            !     endif
            !     if( MOD((cellijk+cellip+cellkp+cellikp),10) .eq. 0 ) then
            !         Vm = Vm - v(i,j,k) * dx1(i) * dx2(j) * dx3(k)
            !         V(i,j,k)=0.d0
            !     endif
            !     if( MOD((cellijk+cellip+celljp+cellijp),10) .eq. 0 ) then
            !         Wm = Wm - w(i,j,k) * dx1(i) * dx2(j) * dx3(k)
            !         W(i,j,k)=0.d0
            !     endif
            ! enddo
            ! enddo
            ! enddo

            ! call MPI_ALLREDUCE(Um    , Um_I,     1, MPI_real_type, MPI_SUM, comm_1d_x1%mpi_comm, ierr)
            ! call MPI_ALLREDUCE(Um_I  , Um_K,     1, MPI_real_type, MPI_SUM, comm_1d_x3%mpi_comm, ierr)
            ! call MPI_ALLREDUCE(Um_K  , Um_total, 1, MPI_real_type, MPI_SUM, comm_1d_x2%mpi_comm, ierr)
            ! call MPI_ALLREDUCE(Vm    , Vm_I,     1, MPI_real_type, MPI_SUM, comm_1d_x1%mpi_comm, ierr)
            ! call MPI_ALLREDUCE(Vm_I  , Vm_K,     1, MPI_real_type, MPI_SUM, comm_1d_x3%mpi_comm, ierr)
            ! call MPI_ALLREDUCE(Vm_K  , Vm_total, 1, MPI_real_type, MPI_SUM, comm_1d_x2%mpi_comm, ierr)
            ! call MPI_ALLREDUCE(Wm    , Wm_I,     1, MPI_real_type, MPI_SUM, comm_1d_x1%mpi_comm, ierr)
            ! call MPI_ALLREDUCE(Wm_I  , Wm_K,     1, MPI_real_type, MPI_SUM, comm_1d_x3%mpi_comm, ierr)
            ! call MPI_ALLREDUCE(Wm_K  , Wm_total, 1, MPI_real_type, MPI_SUM, comm_1d_x2%mpi_comm, ierr)
        
            ! Um_total = Um_total / Volume
            ! Vm_total = Vm_total / Volume 
            ! Wm_total = Wm_total / Volume
        
            ! Um=0.d0; Vm=0.d0; Wm=0.d0
        
            ! do k=0,n3sub
            !     kp = k + 1
            ! do j=0,n2sub
            !     jp = j + 1
            ! do i=0,n1sub
            !     ip = i + 1
            !     cellijk=cell(i ,j ,k )
            !     cellip =cell(ip,j ,k )
            !     celljp =cell(i ,jp,k )
            !     cellkp =cell(i ,j ,kp)
            !     cellijp=cell(ip,jp,k )
            !     cellikp=cell(ip,j ,kp)
            !     celljkp=cell(i ,jp,kp)
        
            !     if( MOD((cellijk+celljp+cellkp+celljkp),10) .ne. 0 ) then
            !         u(i,j,k)=u(i,j,k)*(1.d0/Um_total)
            !     endif
            !     if( MOD((cellijk+cellip+cellkp+cellikp),10) .ne. 0 ) then
            !         !  v(i,j,k)=v(i,j,k)
            !     endif
            !     if( MOD((cellijk+cellip+celljp+cellijp),10) .ne. 0 ) then
            !         ! w(i,j,k)=w(i,j,k)*(1.0d0/Wm_total)
            !     endif
        
            !     Um = Um + U(i,j,k) * dx1(i) * dx2(j) * dx3(k)
            !     Vm = Vm + V(i,j,k) * dx1(i) * dx2(j) * dx3(k)
            !     Wm = Wm + W(i,j,k) * dx1(i) * dx2(j) * dx3(k)
            ! enddo
            ! enddo
            ! enddo
        !Cell classification end

        call MPI_ALLREDUCE(Um    , Um_I,     1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_1d_x1%mpi_comm, ierr)
        call MPI_ALLREDUCE(Um_I  , Um_K,     1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_1d_x3%mpi_comm, ierr)
        call MPI_ALLREDUCE(Um_K  , Um_total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_1d_x2%mpi_comm, ierr)
        call MPI_ALLREDUCE(Vm    , Vm_I,     1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_1d_x1%mpi_comm, ierr)
        call MPI_ALLREDUCE(Vm_I  , Vm_K,     1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_1d_x3%mpi_comm, ierr)
        call MPI_ALLREDUCE(Vm_K  , Vm_total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_1d_x2%mpi_comm, ierr)
        call MPI_ALLREDUCE(Wm    , Wm_I,     1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_1d_x1%mpi_comm, ierr)
        call MPI_ALLREDUCE(Wm_I  , Wm_K,     1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_1d_x3%mpi_comm, ierr)
        call MPI_ALLREDUCE(Wm_K  , Wm_total, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm_1d_x2%mpi_comm, ierr)
    
        Um_total = Um_total / Volume
        Vm_total = Vm_total / Volume 
        Wm_total = Wm_total / Volume
    
        if(myrank==0) Write(*,'(1A20,1A5,1F8.3,1A5,1F8.3,1A5,1F8.3)') "Initial Condition","Um:", Um_total,"Vm:", Vm_total, "Wm:", Wm_total
    
        deallocate(RanNum1,RanNum2,RanNum3)
    
    end subroutine mpi_momentum_initial_channel
    
    !>
    !> @brief       Initialize velocity boundary condition
    !>
    subroutine mpi_momentum_boundary()

        use mpi_topology,   only : comm_1d_x1, comm_1d_x2, comm_1d_x3
        use mpi_subdomain,  only : n1sub, n2sub, n3sub

        implicit none

        integer :: i, j, k
        
       ! Example problem : Constant temperature for upper and lower walls in y-direction
        ! UBCup_sub(:,:)=Uup
        ! VBCup_sub(:,:)=Vup
        ! WBCup_sub(:,:)=Wup
        
        ! UBCbt_sub(:,:)=Ubt
        ! VBCbt_sub(:,:)=Vbt
        ! WBCbt_sub(:,:)=Wbt

        do k=0,n3sub
        do j=0,n2sub
            XMBC(j,k,1,2)=Uup
            XMBC(j,k,2,2)=Vup
            XMBC(j,k,3,2)=Wup
            
            XMBC(j,k,1,1)=Ubt
            XMBC(j,k,2,1)=Vbt
            XMBC(j,k,3,1)=Wbt
        enddo
        enddo

        do k=0,n3sub
        do i=0,n1sub
            YMBC(i,k,1,2)=0.0d0
            YMBC(i,k,2,2)=0.0d0
            YMBC(i,k,3,2)=0.0d0
            
            YMBC(i,k,1,1)=0.0d0
            YMBC(i,k,2,1)=0.0d0
            YMBC(i,k,3,1)=0.0d0
        enddo
        enddo

        do j=0,n2sub
        do i=0,n1sub
            ZMBC(i,j,1,2)=0.0d0
            ZMBC(i,j,2,2)=0.0d0
            ZMBC(i,j,3,2)=0.0d0
            
            ZMBC(i,j,1,1)=0.0d0
            ZMBC(i,j,2,1)=0.0d0
            ZMBC(i,j,3,1)=0.0d0
        enddo
        enddo

        ! if(comm_1d_x2%myrank==0) then
        !     do k=0, n3sub
        !         do i=0, n1sub
        !             U(i,0,k)=UBCbt_sub(i,k)
        !             W(i,0,k)=WBCbt_sub(i,k)
        !             V(i,0,k)=VBCbt_sub(i,k)
        !             V(i,1,k)=VBCbt_sub(i,k)
        !         enddo
        !     enddo
        ! else if(comm_1d_x2%myrank==comm_1d_x2%nprocs-1) then
        !     do k=0, n3sub
        !         do i=0, n1sub
        !             U(i,n2sub,k)=UBCup_sub(i,k)
        !             W(i,n2sub,k)=WBCup_sub(i,k)
        !             V(i,n2sub,k)=VBCup_sub(i,k)
        !         enddo
        !     enddo
        ! endif    


        if(pbc1==.False.) then
            if(comm_1d_x1%myrank==0) then
                do k=0, n3sub
                do j=0, n2sub
                    u(0,j,k)=XMBC(j,k,1,1)
                    u(1,j,k)=XMBC(j,k,1,1)
                    v(0,j,k)=XMBC(j,k,2,1)
                    w(0,j,k)=XMBC(j,k,3,1)
                enddo
                enddo
            endif 
            if(comm_1d_x1%myrank==comm_1d_x1%nprocs-1) then
                do k=0, n3sub
                do j=0, n2sub
                    u(n1sub,j,k)=XMBC(j,k,1,2)
                    v(n1sub,j,k)=XMBC(j,k,2,2)
                    w(n1sub,j,k)=XMBC(j,k,3,2)   
                enddo
                enddo
            endif
        endif
        if(pbc2==.False.) then
            if(comm_1d_x2%myrank==0) then
                do k=0, n3sub
                do i=0, n1sub
                    u(i,0,k)=YMBC(i,k,1,1)
                    v(i,0,k)=YMBC(i,k,2,1)
                    v(i,1,k)=YMBC(i,k,2,1)
                    w(i,0,k)=YMBC(i,k,3,1)   
                enddo
                enddo
            endif 
            if(comm_1d_x2%myrank==comm_1d_x2%nprocs-1) then
                do k=0, n3sub
                do i=0, n1sub
                    u(i,n2sub,k)=YMBC(i,k,1,2)
                    v(i,n2sub,k)=YMBC(i,k,2,2)
                    w(i,n2sub,k)=YMBC(i,k,3,2)   
                enddo
                enddo
            endif    
        endif
        if(pbc3==.False.) then
            if(comm_1d_x3%myrank==0) then   !No i,j Dirichlet 
                do j=0, n2sub
                do i=0, n1sub
                    u(i,j,0)=ZMBC(i,j,1,1)
                    v(i,j,0)=ZMBC(i,j,2,1)
                    w(i,j,0)=ZMBC(i,j,3,1)
                    w(i,j,1)=ZMBC(i,j,3,1)   
                enddo
                enddo
            endif 
            if(comm_1d_x3%myrank==comm_1d_x3%nprocs-1) then
                do j=0, n2sub
                do i=0, n1sub
                    u(i,j,n3sub)=ZMBC(i,j,1,2)
                    v(i,j,n3sub)=ZMBC(i,j,2,2)
                    w(i,j,n3sub)=ZMBC(i,j,3,2)   
                enddo
                enddo
            endif    
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

        if(problem.eq.0) then
            do k=0,n3sub
            do j=0,n2sub
            do i=0,n1sub
                tmp1 = DeltaT*T(i,j,k)
                tmp2 = tmp1 * tmp1
                tmp3 = tmp2 * tmp1
                tmp4 = tmp3 * tmp1
                tmp5 = tmp4 * tmp1

                Mu(i,j,k) = (a10 + a10*(a11*tmp1 + a12*tmp2 + a13*tmp3 + a14*tmp4 + a15*tmp5)) &
                            *(d10 + d10*(d11*tmp1 + d12*tmp2 + d13*tmp3 + d14*tmp4 + d15*tmp5))/Mu0

                invRho(i,j,k) = Rho0/(a10 + a10*(a11*tmp1 + a12*tmp2 + a13*tmp3 + a14*tmp4 + a15*tmp5))
            enddo
            enddo
            enddo
        else if(problem.eq.1) then
            do k=0,n3sub
            do j=0,n2sub
            do i=0,n1sub
                Mu(i,j,k) = 1.0d0
                invRho(i,j,k) = 1.0d0
            enddo
            enddo
                enddo
        else if(problem.eq.2) then
            do k=0,n3sub
            do j=0,n2sub
            do i=0,n1sub
                Mu(i,j,k) = 1.0d0
                invRho(i,j,k) = 1.0d0
            enddo
            enddo
            enddo
        endif

    end subroutine mpi_momentum_coeffi

    !>
    !> @brief       Large eddy simulation - Eddy viscosity model:  Smagorisky model 
    !>
    subroutine mpi_momentum_LES_constant_sm(Mu)

        use mpi_subdomain,  only : n1sub, n2sub, n3sub
        use mpi_subdomain,  only : x2_sub
        use mpi_subdomain,  only : dx1_sub,dx2_sub,dx3_sub,dmx1_sub,dmx2_sub,dmx3_sub

        implicit none

        !> @{ Local indexing variables
        integer :: i,j,k
        integer :: im,jm,km
        integer :: ip,jp,kp
        integer :: jup,jum
        !> @}

        !> @{ Local variables for LES calc 
        double precision :: delta, damping, abs_strain
        double precision :: SR11, SR22, SR33, SR12, SR13, SR23
        double precision :: U1,U2,V1,V2,W1,W2 
        double precision :: dis1,dis2
        !> @}

        !> @{ Local pointer for subdomain variables
        double precision, dimension(:), pointer :: Y      
        double precision, dimension(:), pointer :: dx1, dx2, dx3
        double precision, dimension(:), pointer :: dmx1, dmx2, dmx3
        !> @}

        double precision :: Es_Retau, Ap, Cs, Cmu, Re
        double precision, dimension(0:n1sub,0:n2sub,0:n3sub) :: Mu 

        ! Pointer of grid information
        Y => x2_sub
        dx1 => dx1_sub
        dx2 => dx2_sub
        dx3 => dx3_sub
        dmx1 => dmx1_sub
        dmx2 => dmx2_sub
        dmx3 => dmx3_sub

        Ap= 0.1d0
        Cs= 0.1d0
        Es_Retau=(2.0d0*Re_c)**0.88d0*0.09d0

        ! Get strain rate at cell center
        do k=1,n3sub-1
            kp=k+1
            km=k-1
        do j=1,n2sub-1
            jp=j+1
            jm=j-1
        do i=1,n1sub-1
            ip=i+1
            im=i-1

            delta = (dx1(i)*dx2(j)*dx3(k))**(1./3.)
            damping = 1.0d0

            dis1=abs(Y(j)-x2_start)
            dis2=abs(Y(j)-x2_end)

            if(dis1.gt.dis2) damping=1.0d0 - exp(-(x2_end-Y(j))*(Es_Retau/Ap))
            if(dis1.le.dis2) damping=1.0d0 - exp(-Y(j)*(Es_Retau/Ap))

            SR11= ( U(ip,j ,k )-U(i ,j ,k ) )/dx1(i)
            SR22= ( V(i ,jp,k )-V(i ,j ,k ) )/dx2(j)
            SR33= ( W(i ,j ,kp)-W(i ,j ,k ) )/dx3(k)

            U1= ( 0.5d0/dx1(i)*( dmx1(i)*U(ip,j ,k) + dmx1(ip)*U(i,j ,k) )*0.5d0*dx2(jm) + &
                & 0.5d0/dx1(i)*( dmx1(i)*U(ip,jm,k) + dmx1(ip)*U(i,jm,k) )*0.5d0*dx2(j ) )/dmx2(j )
            U2= ( 0.5d0/dx1(i)*( dmx1(i)*U(ip,jp,k) + dmx1(ip)*U(i,jp,k) )*0.5d0*dx2(j ) + &
                & 0.5d0/dx1(i)*( dmx1(i)*U(ip,j ,k) + dmx1(ip)*U(i,j ,k) )*0.5d0*dx2(jp) )/dmx2(jp)
            V1= ( 0.5d0/dx2(j)*( dmx2(j)*V(i ,jp,k) + dmx2(jp)*V(i ,j,k) )*0.5d0*dx1(im) + &
                & 0.5d0/dx2(j)*( dmx2(j)*V(im,jp,k) + dmx2(jp)*V(im,j,k) )*0.5d0*dx1(i ) )/dmx1(i )
            V2= ( 0.5d0/dx2(j)*( dmx2(j)*V(ip,jp,k) + dmx2(jp)*V(ip,j,k) )*0.5d0*dx1(i ) + &
                & 0.5d0/dx2(j)*( dmx2(j)*V(i ,jp,k) + dmx2(jp)*V(i ,j,k) )*0.5d0*dx1(ip) )/dmx1(ip)
            SR12= 0.5d0*( (U2-U1)/dx2(j) + (V2-V1)/dx1(i) )

            U1= ( 0.5d0/dx1(i)*( dmx1(i)*U(ip,j,k ) + dmx1(ip)*U(i,j,k ) )*0.5d0*dx3(km) + &
                & 0.5d0/dx1(i)*( dmx1(i)*U(ip,j,km) + dmx1(ip)*U(i,j,km) )*0.5d0*dx3(k ) )/dmx3(k )
            U2= ( 0.5d0/dx1(i)*( dmx1(i)*U(ip,j,kp) + dmx1(ip)*U(i,j,kp) )*0.5d0*dx3(k ) + &
                & 0.5d0/dx1(i)*( dmx1(i)*U(ip,j,k ) + dmx1(ip)*U(i,j,k ) )*0.5d0*dx3(kp) )/dmx3(kp)
            W1= ( 0.5d0/dx3(k)*( dmx3(k)*W(i ,j,kp) + dmx3(kp)*W(i ,j,k) )*0.5d0*dx1(im) + &
                & 0.5d0/dx3(k)*( dmx3(k)*W(im,j,kp) + dmx3(kp)*W(im,j,k) )*0.5d0*dx1(i ) )/dmx1(i )
            W2= ( 0.5d0/dx3(k)*( dmx3(k)*W(ip,j,kp) + dmx3(kp)*W(ip,j,k) )*0.5d0*dx1(i ) + &
                & 0.5d0/dx3(k)*( dmx3(k)*W(i ,j,kp) + dmx3(kp)*W(i ,j,k) )*0.5d0*dx1(ip) )/dmx1(ip)
            SR13= 0.5d0*( (U2-U1)/dx3(k) + (W2-W1)/dx1(i) )

            V1= ( 0.5d0/dx2(j)*( dmx2(j)*V(i,jp,k ) + dmx2(jp)*V(i,j,k ) )*0.5d0*dx3(km) + &
                & 0.5d0/dx2(j)*( dmx2(j)*V(i,jp,km) + dmx2(jp)*V(i,j,km) )*0.5d0*dx3(k ) )/dmx3(k )
            V2= ( 0.5d0/dx2(j)*( dmx2(j)*V(i,jp,kp) + dmx2(jp)*V(i,j,kp) )*0.5d0*dx3(k ) + &
                & 0.5d0/dx2(j)*( dmx2(j)*V(i,jp,k ) + dmx2(jp)*V(i,j,k ) )*0.5d0*dx3(kp) )/dmx3(kp)
            W1= ( 0.5d0/dx3(k)*( dmx3(k)*W(i,j ,kp) + dmx3(kp)*W(i,j ,k) )*0.5d0*dx2(jm) + &
                & 0.5d0/dx3(k)*( dmx3(k)*W(i,jm,kp) + dmx3(kp)*W(i,jm,k) )*0.5d0*dx2(j ) )/dmx2(j )
            W2= ( 0.5d0/dx3(k)*( dmx3(k)*W(i,jp,kp) + dmx3(kp)*W(i,jp,k) )*0.5d0*dx2(j ) + &
                & 0.5d0/dx3(k)*( dmx3(k)*W(i,j ,kp) + dmx3(kp)*W(i,j ,k) )*0.5d0*dx2(jp) )/dmx2(jp)
            SR23= 0.5d0*( (V2-V1)/dx3(k) + (W2-W1)/dx2(j) )

            abs_strain= sqrt( 2.0d0*SR11**2.0d0 + 2.0d0*SR22**2.0d0 + 2.0d0*SR33**2.0d0 + &
                          & + 4.0d0*SR12**2.0d0 + 4.0d0*SR13**2.0d0 + 4.0d0*SR23**2.0d0 )

            ! SR(j)=SR(j) + ( Cs*((delta*damping)**2.0d0*abs_strain )/dble(m1msb*m3msb)*dt
            ! SR(j)=SR(j) + ( (Cs**2.0d0*((delta*damping)**2.0d0)*abs_strain)/CmU + 1.0d0) / dble(m1msb*m3msb)*dt

            Mu(i,j,k)     = (Cs**2.0d0*((delta*damping)**2.0d0)*abs_strain)/Cmu + 1.0d0
        enddo
        enddo
        enddo

    end subroutine mpi_momentum_LES_constant_sm

    !>
    !> @brief       Large eddy simulation - Eddy viscosity model:  Vreman model 
    !>
    subroutine mpi_momentum_LES_constant_vr(Mu)

        use mpi_subdomain,  only : n1sub, n2sub, n3sub
        use mpi_subdomain,  only : x2_sub
        use mpi_subdomain,  only : dx1_sub,dx2_sub,dx3_sub,dmx1_sub,dmx2_sub,dmx3_sub

        implicit none

        !> @{ Local indexing variables
        integer :: i,j,k
        integer :: im,jm,km
        integer :: ip,jp,kp
        !> @}

        !> @{ Local variables for LES calc 
        double precision :: BB11, BB22, BB33, BB12, BB13, BB23, BB, AIJAIJ     
        double precision :: DV1DX1, DV1DX2, DV1DX3, DV2DX1, DV2DX2, DV2DX3, DV3DX1, DV3DX2, DV3DX3
        double precision :: U1,U2,V1,V2,W1,W2 
        !> @}

        !> @{ Local pointer for subdomain variables
        double precision, dimension(:), pointer :: Y      
        double precision, dimension(:), pointer :: dx1, dx2, dx3
        double precision, dimension(:), pointer :: dmx1, dmx2, dmx3
        !> @}

        !> @{ Coefficients for LES calc 
        double precision :: Es_Retau, Ap, Cs, Re
        double precision :: Cv, temp
        double precision :: damping,delta,dis1,dis2 
        !> @}

        double precision, dimension(0:n1sub,0:n2sub,0:n3sub) :: Mu 


        ! Pointer of grid information
        Y => x2_sub
        dx1 => dx1_sub
        dx2 => dx2_sub
        dx3 => dx3_sub
        dmx1 => dmx1_sub
        dmx2 => dmx2_sub
        dmx3 => dmx3_sub

        Ap= 0.1d0
        Cs= 0.1d0
        Es_Retau=(2.0d0*Re_c)**0.88d0*0.09d0
        Cv=2.5d0*Cs*Cs

        ! Get strain rate at cell center
        do k=1,n3sub-1
            kp=k+1
            km=k-1
        do j=1,n2sub-1
            jp=j+1
            jm=j-1
        do i=1,n1sub-1
            ip=i+1
            im=i-1

            delta = (dx1(i)*dx2(j)*dx3(k))**(1.d0/3.d0)
            damping = 1.d0
            dis1=dabs(Y(j)-x2_start)
            dis2=dabs(Y(j)-x2_end)

            if(dis1.gt.dis2) damping=1.0d0 - exp(-(x2_end-Y(j))*(Es_Retau/Ap))
            if(dis1.le.dis2) damping=1.0d0 - exp(-Y(j)*(Es_Retau/Ap))

            DV1DX1= ( U(ip,j ,k )-U(i ,j ,k ) )/dx1(i)
            DV2DX2= ( V(i ,jp,k )-V(i ,j ,k ) )/dx2(j)
            DV3DX3= ( W(i ,j ,kp)-W(i ,j ,k ) )/dx3(k)

            U1= ( 0.5d0/dx1(i)*( dmx1(i)*U(ip,j ,k) + dmx1(ip)*U(i,j ,k) )*0.5d0*dx2(jm) + &
                & 0.5d0/dx1(i)*( dmx1(i)*U(ip,jm,k) + dmx1(ip)*U(i,jm,k) )*0.5d0*dx2(j ) )/dmx2(j )
            U2= ( 0.5d0/dx1(i)*( dmx1(i)*U(ip,jp,k) + dmx1(ip)*U(i,jp,k) )*0.5d0*dx2(j ) + &
                & 0.5d0/dx1(i)*( dmx1(i)*U(ip,j ,k) + dmx1(ip)*U(i,j ,k) )*0.5d0*dx2(jp) )/dmx2(jp)
            DV1DX2= (U2-U1)/dx2(j)

            U1= ( 0.5d0/dx1(i)*( dmx1(i)*U(ip,j,k ) + dmx1(ip)*U(i,j,k ) )*0.5d0*dx3(km) + &
                & 0.5d0/dx1(i)*( dmx1(i)*U(ip,j,km) + dmx1(ip)*U(i,j,km) )*0.5d0*dx3(k ) )/dmx3(k )
            U2= ( 0.5d0/dx1(i)*( dmx1(i)*U(ip,j,kp) + dmx1(ip)*U(i,j,kp) )*0.5d0*dx3(k ) + &
                & 0.5d0/dx1(i)*( dmx1(i)*U(ip,j,k ) + dmx1(ip)*U(i,j,k ) )*0.5d0*dx3(kp) )/dmx3(kp)
            DV1DX3= (U2-U1)/dx3(k)

            V1= ( 0.5d0/dx2(j)*( dmx2(j)*V(i ,jp,k) + dmx2(jp)*V(i ,j,k) )*0.5d0*dx1(im) + &
                & 0.5d0/dx2(j)*( dmx2(j)*V(im,jp,k) + dmx2(jp)*V(im,j,k) )*0.5d0*dx1(i ) )/dmx1(i )
            V2= ( 0.5d0/dx2(j)*( dmx2(j)*V(ip,jp,k) + dmx2(jp)*V(ip,j,k) )*0.5d0*dx1(i ) + &
                & 0.5d0/dx2(j)*( dmx2(j)*V(i ,jp,k) + dmx2(jp)*V(i ,j,k) )*0.5d0*dx1(ip) )/dmx1(ip)
            DV2DX1= (V2-V1)/dx1(i)

            V1= ( 0.5d0/dx2(j)*( dmx2(j)*V(i,jp,k ) + dmx2(jp)*V(i,j,k ) )*0.5d0*dx3(km) + &
                & 0.5d0/dx2(j)*( dmx2(j)*V(i,jp,km) + dmx2(jp)*V(i,j,km) )*0.5d0*dx3(k ) )/dmx3(k )
            V2= ( 0.5d0/dx2(j)*( dmx2(j)*V(i,jp,kp) + dmx2(jp)*V(i,j,kp) )*0.5d0*dx3(k ) + &
                & 0.5d0/dx2(j)*( dmx2(j)*V(i,jp,k ) + dmx2(jp)*V(i,j,k ) )*0.5d0*dx3(kp) )/dmx3(kp)
            DV2DX3= (V2-V1)/dx3(k)

            W1= ( 0.5d0/dx3(k)*( dmx3(k)*W(i ,j,kp) + dmx3(kp)*W(i ,j,k) )*0.5d0*dx1(im) + &
                & 0.5d0/dx3(k)*( dmx3(k)*W(im,j,kp) + dmx3(kp)*W(im,j,k) )*0.5d0*dx1(i ) )/dmx1(i )
            W2= ( 0.5d0/dx3(k)*( dmx3(k)*W(ip,j,kp) + dmx3(kp)*W(ip,j,k) )*0.5d0*dx1(i ) + &
                & 0.5d0/dx3(k)*( dmx3(k)*W(i ,j,kp) + dmx3(kp)*W(i ,j,k) )*0.5d0*dx1(ip) )/dmx1(ip)             
            DV3DX1= (W2-W1)/dx1(i)

            W1= ( 0.5d0/dx3(k)*( dmx3(k)*W(i,j ,kp) + dmx3(kp)*W(i,j ,k) )*0.5d0*dx2(jm) + &
                & 0.5d0/dx3(k)*( dmx3(k)*W(i,jm,kp) + dmx3(kp)*W(i,jm,k) )*0.5d0*dx2(j ) )/dmx2(j )
            W2= ( 0.5d0/dx3(k)*( dmx3(k)*W(i,jp,kp) + dmx3(kp)*W(i,jp,k) )*0.5d0*dx2(j ) + &
                & 0.5d0/dx3(k)*( dmx3(k)*W(i,j ,kp) + dmx3(kp)*W(i,j ,k) )*0.5d0*dx2(jp) )/dmx2(jp)
            DV3DX2= (W2-W1)/dx2(j)

            ! vreman paper 
            !BB11 = DV1DX1*DV1DX1*(dx1(i)*dx1(i)) + DV1DX2*DV1DX2*(dx2(j)*dx2(j)) + DV1DX3*DV1DX3*(dx3(k)*dx3(k))
            !BB22 = DV2DX1*DV2DX1*(dx1(i)*dx1(i)) + DV2DX2*DV2DX2*(dx2(j)*dx2(j)) + DV2DX3*DV2DX3*(dx3(k)*dx3(k))
            !BB33 = DV3DX1*DV3DX1*(dx1(i)*dx1(i)) + DV3DX2*DV3DX2*(dx2(j)*dx2(j)) + DV3DX3*DV3DX3*(dx3(k)*dx3(k))
            !BB12 = DV1DX1*DV2DX1*(dx1(i)*dx1(i)) + DV1DX2*DV2DX2*(dx2(j)*dx2(j)) + DV1DX3*DV2DX3*(dx3(k)*dx3(k))
            !BB13 = DV1DX1*DV3DX1*(dx1(i)*dx1(i)) + DV1DX2*DV3DX2*(dx2(j)*dx2(j)) + DV1DX3*DV3DX3*(dx3(k)*dx3(k))
            !BB23 = DV2DX1*DV3DX1*(dx1(i)*dx1(i)) + DV2DX2*DV3DX2*(dx2(j)*dx2(j)) + DV2DX3*DV3DX3*(dx3(k)*dx3(k))
            
            ! +) considering damping 
            BB11 = DV1DX1*DV1DX1*(damping*dx1(i)) + DV1DX2*DV1DX2*(damping*dx2(j)) + DV1DX3*DV1DX3*(damping*dx3(k))
            BB22 = DV2DX1*DV2DX1*(damping*dx1(i)) + DV2DX2*DV2DX2*(damping*dx2(j)) + DV2DX3*DV2DX3*(damping*dx3(k))
            BB33 = DV3DX1*DV3DX1*(damping*dx1(i)) + DV3DX2*DV3DX2*(damping*dx2(j)) + DV3DX3*DV3DX3*(damping*dx3(k))
            BB12 = DV1DX1*DV2DX1*(damping*dx1(i)) + DV1DX2*DV2DX2*(damping*dx2(j)) + DV1DX3*DV2DX3*(damping*dx3(k))
            BB13 = DV1DX1*DV3DX1*(damping*dx1(i)) + DV1DX2*DV3DX2*(damping*dx2(j)) + DV1DX3*DV3DX3*(damping*dx3(k))
            BB23 = DV2DX1*DV3DX1*(damping*dx1(i)) + DV2DX2*DV3DX2*(damping*dx2(j)) + DV2DX3*DV3DX3*(damping*dx3(k))

            BB = BB11*BB22 - BB12*BB12 + BB11*BB33 - BB13*BB13 + BB22*BB33 - BB23*BB23
        
            AIJAIJ = DV1DX1*DV1DX1 + DV2DX1*DV2DX1 + DV3DX1*DV3DX1    &
                   + DV1DX2*DV1DX2 + DV2DX2*DV2DX2 + DV3DX2*DV3DX2    & 
                   + DV1DX3*DV1DX3 + DV2DX3*DV2DX3 + DV3DX3*DV3DX3    
        
            ! For Single Precision
            temp = BB / AIJAIJ

            if(temp.lt.1e-12) then
                temp = 0.d0
            endif
            
            Mu(i,j,k)     = ( Cv * SQRT ( temp ) ) / Cmu + 1.d0

        enddo
        enddo
        enddo

    end subroutine mpi_momentum_LES_constant_vr

    !>
    !> @brief       Update velocity field with the solved incremental velocity fields
    !>
    subroutine mpi_momentum_pseudoupdateUVW()

        use mpi_subdomain,  only : n1sub, n2sub, n3sub, n1msub, n2msub, n3msub
        use mpi_subdomain,  only : iS_BC, jS_BC, kS_BC

        implicit none
        integer :: i,j,k
        integer :: iuc,jvc,kwc

        do k = 1, n3msub
        do j = 1, n2msub
        do i = 1, n1msub
            iuc = iS_BC(i)
            jvc = jS_BC(j)
            kwc = kS_BC(k)
            
            u(i,j,k)=u(i,j,k) + du(i,j,k)*dble(iuc)
            v(i,j,k)=v(i,j,k) + dv(i,j,k)*dble(jvc)
            w(i,j,k)=w(i,j,k) + dw(i,j,k)*dble(kwc)
        enddo
        enddo
        enddo

    end subroutine mpi_momentum_pseudoupdateUVW

    subroutine mpi_momentum_deallocate_dUVW()
        implicit none
        
        deallocate(dU)
        deallocate(dV)
        deallocate(dW)

    end subroutine mpi_momentum_deallocate_dUVW

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
        use mpi_subdomain,  only : i_indexS, iS_BC, jC_BC, kC_BC

        use PaScaL_TDMA

        implicit none

        double precision, dimension(0:n1sub,0:n2sub,0:n3sub) ::  T

        type(ptdma_plan_many)     :: ptdma_plan

        !> @{ Local pointer for subdomain variables
        double precision, dimension(:), pointer :: dx1, dx2, dx3
        double precision, dimension(:), pointer :: dmx1, dmx2, dmx3
        !> @}

        !> @{ Local indexing variables
        integer :: im, i, ip, ium, iuc, iup
        integer :: jm, j, jp, jum, juc, jup
        integer :: km, k, kp, kum, kuc, kup
        !> @}

        !> @{ Local arrays for matrix RHS coefficients
        double precision, allocatable, dimension(:,:,:) :: ddU, RHS, RHSI, RHSJ
        double precision, allocatable, dimension(:,:,:) :: AMK, ACK, APK, AMI, ACI, API, AMJ, ACJ, APJ
        !> @}

        !> @{ Local variables for matrix and RHS coefficients calculation
        double precision :: u1,u2,v3,v4,w5,w6
        double precision :: dudx1,dudx2,dudy3,dudy4,dudz5,dudz6
        double precision :: dvdx3,dvdx4,dwdx5,dwdx6
        double precision :: mua, mub, muc, mu3, mu4, mu5, mu6
        double precision :: viscous_u1,viscous_u2,viscous_u3,viscous_u12,viscous_u13,Tc
        double precision :: ubc_down,ubc_up
        double precision :: M11MI, M11CI, M11PI, M11MJ, M11CJ, M11PJ, M11MK, M11CK, M11PK
        double precision :: M12Vn, M13Wn
        double precision :: M12Vall_down, M12Vall_up, M13Wall_down, M13Wall_up
        double precision :: invRhoc, invRhocCmu_half

        double precision :: RHS_ijk
        double precision :: Tijk, Tim
        double precision :: Uijk, Uip, Uim, Ujp, Ujm, Ukp, Ukm
        double precision :: Vijk, Vip, Vim, Vjp, Vjm, Vkp, Vkm, Vimjp
        double precision :: Wijk, Wip, Wim, Wjp, Wjm, Wkp, Wkm, Wimkp
        double precision :: Muijk, Mujm, Muimjm, Mujp, Muimjp, Mukm, Muimkm, Mukp, Muimkp, Muim
        double precision :: UBCbt, UBCup, VBCbtim, VBCbt, VBCupim, VBCup, WBCbtim, WBCbt, WBCupim, WBCup
    
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
            km=k-1
            kp=k+1
            kum = kC_BC(km); kuc = kC_BC(k); kup = kC_BC(kp) 
            do j = 1, n2msub
                jm = j-1
                jp = j+1
                jum = jC_BC(jm); juc = jC_BC(j); jup = jC_BC(jp)
                do i = 1, n1msub
                    im = max(i,i_indexS)-1 ! To prevent NaN values
                    ip = i+1
                    ium = iS_BC(im); iuc = iS_BC(i); iup = iS_BC(ip)
                            
                    Uijk = u(i ,j ,k );
                    Uip  = u(ip,j ,k ); Uim  = u(im,j ,k )
                    Ujp  = u(i ,jp,k ); Ujm  = u(i ,jm,k )
                    Ukp  = u(i ,j ,kp); Ukm  = u(i ,j ,km)
                
                    u1 = 0.5d0*(Uim + Uijk)            
                    u2 = 0.5d0*(Uip + Uijk)  
                    dudx1 = (Uijk - Uim )/ dx1(im)
                    dudx2 = (Uip  - Uijk)/ dx1(i )
                    dudy3 = (Uijk - Ujm )/dmx2(j )
                    dudy4 = (Ujp  - Uijk)/dmx2(jp)
                    dudz5 = (Uijk - Ukm )/dmx3(k )
                    dudz6 = (Ukp  - Uijk)/dmx3(kp)
                
                    Vijk = v(i,j ,k) ; Vim   = v(im,j ,k)
                    Vjp  = v(i,jp,k) ; Vimjp = v(im,jp,k)
                    v3 = 0.5d0*(dx1(im)*Vijk + dx1(i)*Vim  )/dmx1(i)
                    v4 = 0.5d0*(dx1(im)*Vjp  + dx1(i)*Vimjp)/dmx1(i)
                    dvdx3 = (Vijk - Vim  )/dmx1(i)
                    dvdx4 = (Vjp  - Vimjp)/dmx1(i)
                
                    Wijk = w(i,j,k ) ; Wim  = w(im,j,k )
                    Wkp  = w(i,j,kp) ; Wimkp= w(im,j,kp)
                    w5 = 0.5d0*(dx1(im)*Wijk + dx1(i)*Wim  )/dmx1(i)
                    w6 = 0.5d0*(dx1(im)*Wkp  + dx1(i)*Wimkp)/dmx1(i)
                    dwdx5 = (Wijk - Wim  )/dmx1(i)
                    dwdx6 = (Wkp  - Wimkp)/dmx1(i)
                
                    Muijk = Mu(i ,j,k); Mujm   = Mu(i ,jm,k); Mujp   = Mu(i ,jp,k); Mukm  =Mu(i ,j,km); Mukp  =Mu(i ,j,kp)
                    Muim  = Mu(im,j,k); Muimjm = Mu(im,jm,k); Muimjp = Mu(im,jp,k); Muimkm=Mu(im,j,km); Muimkp=Mu(im,j,kp)    
                    muc = 0.5d0*(dx1(im)*Muijk + dx1(i)*Muim  )/dmx1(i )
                    mua = 0.5d0*(dx1(im)*Mujm  + dx1(i)*Muimjm)/dmx1(i )
                    mub = 0.5d0*(dx1(im)*Mujp  + dx1(i)*Muimjp)/dmx1(i )
                    mu3 = 0.5d0*(dx2(jm)*muc   + dx2(j)*mua   )/dmx2(j )
                    mu4 = 0.5d0*(dx2(jp)*muc   + dx2(j)*mub   )/dmx2(jp)
                
                    !muc = 0.5d0*(dx1(im)*Muijk + dx1(i)*Muim  )/dmx1(i )
                    mua = 0.5d0*(dx1(im)*Mukm  + dx1(i)*Muimkm)/dmx1(i )
                    mub = 0.5d0*(dx1(im)*Mukp  + dx1(i)*Muimkp)/dmx1(i )
                    mu5 = 0.5d0*(dx3(km)*muc   + dx3(k)*mua   )/dmx3(k )
                    mu6 = 0.5d0*(dx3(kp)*muc   + dx3(k)*mub   )/dmx3(kp)
                
                    invRhoc        = 0.5d0*(dx1(im)*invRho(i,j ,k )+dx1(i)*invRho(im,j ,k ) )/dmx1(i)
                    invRhocCmu_half= 0.5d0* Cmu*invRhoc

                    !M11Un--------------------------------------------------------------------------------------------
                    !1   : X-direction (NoBC)
                    !1-1 : Viscous              ! mu * d^2( )/dx^2
                    M11MI = - invRhocCmu_half *  (                  + Muim / dx1(im) ) / dmx1(i) * 2.0d0 
                    M11CI = - invRhocCmu_half *  ( - Muijk / dx1(i) - Muim / dx1(im) ) / dmx1(i) * 2.0d0 
                    M11PI = - invRhocCmu_half *  ( + Muijk / dx1(i)                  ) / dmx1(i) * 2.0d0 
                    !1-2 : Convection              ! Un*d( )/dx                   ! (dUn/dx)*( )       
                    M11MI = M11MI + 0.25d0* (             - u1/dx1(im) + 0.5d0 * (         + dudx1 ) )
                    M11CI = M11CI + 0.25d0* ( - u2/dx1(i) + u1/dx1(im) + 0.5d0 * ( + dudx2 + dudx1 ) )
                    M11PI = M11PI + 0.25d0* ( + u2/dx1(i)              + 0.5d0 * ( + dudx2         ) )

                    !2   : Y-direction (NoBC)
                    !2-1 : Viscous               ! mu * d^2( )/dy^2 
                    M11MJ = - invRhocCmu_half * (                + mu3/dmx2(j) ) / dx2(j)
                    M11CJ = - invRhocCmu_half * ( - mu4/dmx2(jp) - mu3/dmx2(j) ) / dx2(j)  
                    M11PJ = - invRhocCmu_half * ( + mu4/dmx2(jp)               ) / dx2(j)               
                    !2-2 : Convection         ! Vn*d( )/dy 
                    M11MJ = M11MJ + 0.25d0* (               - v3/dmx2(j) )
                    M11CJ = M11CJ + 0.25d0* ( - v4/dmx2(jp) + v3/dmx2(j) )
                    M11PJ = M11PJ + 0.25d0* ( + v4/dmx2(jp)              )

                    !3   : Z-direction (NoBC)
                    !3-1 : Viscous             ! mu * d^2( )/dz^2 
                    M11MK = - invRhocCmu_half * (                + mu5/dmx3(k) ) / dx3(k) 
                    M11CK = - invRhocCmu_half * ( - mu6/dmx3(kp) - mu5/dmx3(k) ) / dx3(k) 
                    M11PK = - invRhocCmu_half * ( + mu6/dmx3(kp)               ) / dx3(k) 
                    !3-2 : Convection            ! Wn*d( )/dz 
                    M11MK = M11MK + 0.25d0* (               - w5/dmx3(k) )
                    M11CK = M11CK + 0.25d0* ( - w6/dmx3(kp) + w5/dmx3(k) )
                    M11PK = M11PK + 0.25d0* ( + w6/dmx3(kp)              )
    
                    ! M11Un: Z-direction
                    AMK(i,j,k) = ( M11MK * dble(kum) * dt ) * dble(iuc)
                    ACK(i,j,k) = ( M11CK             * dt ) * dble(iuc) + 1.0d0
                    APK(i,j,k) = ( M11PK * dble(kup) * dt ) * dble(iuc)
                    ! M11Un: X-direction
                    AMI(j,k,i) = ( M11MI * dble(ium) * dt ) * dble(iuc)
                    ACI(j,k,i) = ( M11CI             * dt ) * dble(iuc) + 1.0d0
                    API(j,k,i) = ( M11PI * dble(iup) * dt ) * dble(iuc)
                    ! M11Un: Y-direction
                    AMJ(k,i,j) = ( M11MJ * dble(jum) * dt ) * dble(iuc)
                    ACJ(k,i,j) = ( M11CJ             * dt ) * dble(iuc) + 1.0d0
                    APJ(k,i,j) = ( M11PJ * dble(jup) * dt ) * dble(iuc)
                
                    !4   : RHS
                    viscous_u1  = (Muijk*dudx2 - Muim*dudx1) / dmx1(i) 
                    viscous_u2  = (mu4  *dudy4 - mu3 *dudy3) /  dx2(j)
                    viscous_u3  = (mu6  *dudz6 - mu5 *dudz5) /  dx3(k)            
                    viscous_u12 = (mu4  *dvdx4 - mu3 *dvdx3) /  dx2(j) 
                    viscous_u13 = (mu6  *dwdx6 - mu5 *dwdx5) /  dx3(k)
                    !4-1 : rn (* 1/dt*un are offseted in Au^n term.) (NoBC+BC)
                    RHS_ijk  = invRhocCmu_half * ( (viscous_u1 + viscous_u2 + viscous_u3) + (viscous_u1 + viscous_u12 + viscous_u13) )
                
                    !4-2A:       ! M11Un (NoBC)
                    RHS_ijk  = RHS_ijk - ( + M11MI * dble(ium) * Uim + M11MJ * dble(jum) * Ujm + M11MK * dble(kum) * Ukm  &
                                           +(M11CI                   + M11CJ                   + M11CK           ) * Uijk &
                                           + M11PI * dble(iup) * Uip + M11PJ * dble(jup) * Ujp + M11PK * dble(kup) * Ukp  ) 
                
                    !4-2B:       ! M12Vn (NoBC)
                    M12Vn = 0.25d0*( v3*dble(jum) *dudy3 + v4*dble(jup) *dudy4 ) - invRhocCmu_half * ( mu4* dvdx4*dble(jup) - mu3* dvdx3*dble(jum) ) / dx2(j)
                    RHS_ijk   = RHS_ijk - M12Vn
                
                    !4-2C:       ! M13Wn (NoBC)
                    M13Wn = 0.25d0*( w5*dble(kum) *dudz5 + w6*dble(kup) *dudz6 ) - invRhocCmu_half * ( mu6* dwdx6*dble(kup) - mu5* dwdx5*dble(kum) ) / dx3(k)
                    RHS_ijk   = RHS_ijk - M13Wn
                    
                    !4-3 : Pressure
                    RHS_ijk  = RHS_ijk - Cmp*invRhoc * ( p(i,j,k) - p(im,j,k) )/dmx1(i)
                 
                    !4-4 : Buoyancy(we use Cmt*theta for the buoyancy term)
                    !Tijk= T(i ,j,k)
                    !Tim = T(im,j,k)
                    !Tc  = 0.5d0*(dx1(im)*Tijk + dx1(i)*Tim)/dmx1(i)
                    !           !Buoyancy
                    !RHS_ijk = RHS_ijk + Cmt*(Tc + a12pera11*Tc**2.0d0*DeltaT)*invRhoc
                
                    !4-5 : BC
                    ! X direction
                    UBCbt = XMBC(j,k,1,1)
                    UBCup = XMBC(j,k,1,2)
                
                    ubc_down = - ( M11MI ) * UBCbt
                    ubc_up   = - ( M11PI ) * UBCup
                
                    RHS_ijk  = RHS_ijk + dble(1-ium)*ubc_down + dble(1-iup)*ubc_up
                               
                    ! Y direction
                    UBCbt  = YMBC(i ,k,1,1)   
                    UBCup  = YMBC(i ,k,1,2)   
                    VBCbtim= YMBC(im,k,2,1)   
                    VBCbt  = YMBC(i ,k,2,1) 
                    VBCupim= YMBC(im,k,2,2) 
                    VBCup  = YMBC(i ,k,2,2)
                
                    v3    = 0.5d0*(dx1(im)*VBCbt + dx1(i)*VBCbtim) /dmx1(i)
                    v4    = 0.5d0*(dx1(im)*VBCup + dx1(i)*VBCupim) /dmx1(i)
                    dvdx3 = (VBCbt - VBCbtim) /dmx1(i)
                    dvdx4 = (VBCup - VBCupim) /dmx1(i)
                
                    M12Vall_down = 0.25d0* ( v3 * dudy3 ) - invRhocCmu_half * ( - mu3 * dvdx3 ) /dx2(j)
                    M12Vall_up   = 0.25d0* ( v4 * dudy4 ) - invRhocCmu_half * ( + mu4 * dvdx4 ) /dx2(j)
                
                    ubc_down = - ( M11MJ ) * UBCbt - M12Vall_down
                    ubc_up   = - ( M11PJ ) * UBCup - M12Vall_up
                
                    RHS_ijk = RHS_ijk  + dble(1-jum)*ubc_down + dble(1-jup)*ubc_up
                
                    ! Z direction
                    UBCbt  = ZMBC(i ,j,1,1)   
                    UBCup  = ZMBC(i ,j,1,2)   
                    WBCbtim= ZMBC(im,j,3,1)   
                    WBCbt  = ZMBC(i ,j,3,1) 
                    WBCupim= ZMBC(im,j,3,2) 
                    WBCup  = ZMBC(i ,j,3,2) 
                
                    w5    = 0.5d0*(dx1(im)*WBCbt + dx1(i)*WBCbtim) /dmx1(i)
                    w6    = 0.5d0*(dx1(im)*WBCup + dx1(i)*WBCupim) /dmx1(i)
                    dwdx5 = (WBCbt - WBCbtim) /dmx1(i)
                    dwdx6 = (WBCup - WBCupim) /dmx1(i)
                
                    M13Wall_down = 0.25d0* ( w5 * dudz5 ) - invRhocCmu_half * ( - mu5 * dwdx5 ) /dx3(k)
                    M13Wall_up   = 0.25d0* ( w6 * dudz6 ) - invRhocCmu_half * ( + mu6 * dwdx6 ) /dx3(k)
                
                    ubc_down = - ( M11MK ) * UBCbt - M13Wall_down 
                    ubc_up   = - ( M11PK ) * UBCup - M13Wall_up 
                
                    RHS_ijk = RHS_ijk  + dble(1-kum)*ubc_down + dble(1-kup)*ubc_up
                
                    RHS(i,j,k) = RHS_ijk * dt * dble(iuc)
                enddo
            enddo
        enddo

        ! 1st stage : solve TDM in z-direction
        call PaScaL_TDMA_plan_many_create(ptdma_plan, (n1msub)*(n2msub), comm_1d_x3%myrank, comm_1d_x3%nprocs, comm_1d_x3%mpi_comm)
        call PaScaL_TDMA_many_solve_cycle(ptdma_plan, AMK, ACK, APK, RHS,(n1msub)*(n2msub),(n3msub))
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x3%nprocs)

        ! Deallocate A-matrix for the 2nd stage
        deallocate(AMK,ACK,APK)

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
        call PaScaL_TDMA_plan_many_create(ptdma_plan, (n2msub)*(n3msub), comm_1d_x1%myrank, comm_1d_x1%nprocs, comm_1d_x1%mpi_comm)
        call PaScaL_TDMA_many_solve_cycle(ptdma_plan, AMI, ACI, API, RHSI,(n2msub)*(n3msub),(n1msub))
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
        call PaScaL_TDMA_plan_many_create(ptdma_plan, (n3msub)*(n1msub), comm_1d_x2%myrank, comm_1d_x2%nprocs, comm_1d_x2%mpi_comm)
        call PaScaL_TDMA_many_solve(ptdma_plan, AMJ, ACJ, APJ, RHSJ,(n3msub)*(n1msub),(n2msub))
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x2%nprocs)        

        ! Deallocate A-matrix for the 3rd stage
        deallocate(AMJ,ACJ,APJ)

        ! Update the velocity increments with the transpose of the solution array (i,k,j)
        allocate(dU(0:n1sub,0:n2sub,0:n3sub))
        dU=0.0d0
        do k = 1, n3msub
        do j = 1, n2msub
        do i = 1, n1msub
            iuc = iS_BC(i)
            dU(i,j,k)=RHSJ(k,i,j)  * dble(iuc)
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
        use mpi_subdomain,  only : n1sub, n2sub, n3sub, n1msub, n2msub, n3msub
        use mpi_subdomain,  only : dx1_sub,dx2_sub,dx3_sub,dmx1_sub,dmx2_sub,dmx3_sub
        use mpi_subdomain,  only : j_indexS, iC_BC, jS_BC, kC_BC

        use PaScaL_TDMA

        implicit none

        double precision, dimension(0:n1sub,0:n2sub,0:n3sub) ::  T

        type(ptdma_plan_many)     :: ptdma_plan

        !> @{ Local pointer for subdomain variables
        double precision, dimension(:), pointer :: dx1, dx2, dx3
        double precision, dimension(:), pointer :: dmx1, dmx2, dmx3
        !> @}

        !> @{ Local indexing variables
        integer :: im, i, ip, ivm, ivc, ivp
        integer :: jm, j, jp, jvm, jvc, jvp
        integer :: km, k, kp, kvm, kvc, kvp
        !> @}

        !> @{ Local arrays for matrix and RHS coefficients
        double precision, allocatable, dimension(:,:,:) :: ddV, RHS, RHSI, RHSJ
        double precision, allocatable, dimension(:,:,:) :: AMK, ACK, APK, AMI, ACI, API, AMJ, ACJ, APJ
        !> @}

        !> @{ Local variables for matrix and RHS coefficients calculation
        double precision :: u1,u2,v3,v4,w5,w6
        double precision :: ddu1,ddu2,dvdx1,dvdx2,dvdy3,dvdy4,dvdz5,dvdz6
        double precision :: dudy1,dudy2,dwdy5,dwdy6,ddudy1,ddudy2
        double precision :: mua,muc,mub,mu1,mu2,mu5,mu6
        double precision :: viscous_v1,viscous_v2,viscous_v3,viscous_v21,viscous_v23
        double precision :: vbc_down,vbc_up
        double precision :: M22MI, M22CI, M22PI, M22MJ, M22CJ, M22PJ, M22MK, M22CK, M22PK
        double precision :: M21Un, M23Wn
        double precision :: M21Uall_down, M21Uall_up, M23Wall_down, M23Wall_up
        double precision :: invRhoc, invRhocCmu_half
    
        double precision :: RHS_ijk
        double precision :: Uijk, Uip, Uim, Ujp, Ujm, Ukp, Ukm, Uipjm
        double precision :: Vijk, Vip, Vim, Vjp, Vjm, Vkp, Vkm
        double precision :: Wijk, Wip, Wim, Wjp, Wjm, Wkp, Wkm, Wjmkp
        double precision :: dUijk, dUjm, dUip, dUipjm    
        double precision :: Muijk, Muim,Muimjm,Mujm,Muip,Muipjm,Mukm,Mujmkm,Mukp,Mujmkp
        double precision :: VBCbt, VBCup, UBCbtjm, UBCbt, UBCupjm, UBCup, WBCbtjm, WBCbt, WBCupjm, WBCup
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
            km=k-1
            kp=k+1
            kvm = kC_BC(km); kvc = kC_BC(k); kvp = kC_BC(kp) 
            do j =  1, n2msub
                jm = max(j,j_indexS)-1 ! To prevent NaN values
                jp = j+1
                jvm = jS_BC(jm); jvc = jS_BC(j); jvp = jS_BC(jp)
                do i = 1, n1msub
                    im = i-1
                    ip = i+1
                    ivm = iC_BC(im); ivc = iC_BC(i); ivp = iC_BC(ip)

                    Vijk = v(i ,j ,k )
                    Vip  = v(ip,j ,k ); Vim  = v(im,j ,k )
                    Vjp  = v(i ,jp,k ); Vjm  = v(i ,jm,k )
                    Vkp  = v(i ,j ,kp); Vkm  = v(i ,j ,km)
                
                    v3 = 0.5d0*(Vjm + Vijk)
                    v4 = 0.5d0*(Vjp + Vijk)
                    dvdx1 = (Vijk - Vim )/dmx1(i )
                    dvdx2 = (Vip  - Vijk)/dmx1(ip)
                    dvdy3 = (Vijk - Vjm )/ dx2(jm)
                    dvdy4 = (Vjp  - Vijk)/ dx2(j ) 
                    dvdz5 = (Vijk - Vkm )/dmx3(k )
                    dvdz6 = (Vkp  - Vijk)/dmx3(kp)
                    
                    Wijk = w(i,j,k ); Wjm    = w(i,jm,k )
                    Wkp  = w(i,j,kp); Wjmkp  = w(i,jm,kp)
                    w5 = 0.5d0*(dx2(jm)*Wijk + dx2(j)*Wjm  )/dmx2(j)
                    w6 = 0.5d0*(dx2(jm)*Wkp  + dx2(j)*Wjmkp)/dmx2(j)
                    dwdy5 = (Wijk - Wjm  )/dmx2(j)
                    dwdy6 = (Wkp  - Wjmkp)/dmx2(j)
                
                    Uijk = u(i ,j,k); Ujm    = u(i ,jm,k)
                    Uip  = u(ip,j,k); Uipjm  = u(ip,jm,k)
                    u1 = 0.5d0*(dx2(jm)*Uijk + dx2(j)*Ujm  )/dmx2(j)
                    u2 = 0.5d0*(dx2(jm)*Uip  + dx2(j)*Uipjm)/dmx2(j)
                    dudy1 = (Uijk - Ujm  )/dmx2(j)
                    dudy2 = (Uip  - Uipjm)/dmx2(j)
                
                    Muijk = Mu(i,j ,k); Muim   = Mu(im,j ,k); Muip  =Mu(ip,j ,k); Mukm  =Mu(i ,j ,km); Mukp  =Mu(i ,j ,kp)
                    Mujm  = Mu(i,jm,k); Muimjm = Mu(im,jm,k); Muipjm=Mu(ip,jm,k); Mujmkm=Mu(i ,jm,km); Mujmkp=Mu(i ,jm,kp)
                    
                    muc = 0.5d0*(dx2(jm)*Muijk + dx2(j)*Mujm  )/dmx2(j)
                    mua = 0.5d0*(dx2(jm)*Muim  + dx2(j)*Muimjm)/dmx2(j)
                    mub = 0.5d0*(dx2(jm)*Muip  + dx2(j)*Muipjm)/dmx2(j)
                    mu1 = 0.5d0*(dx1(im)*muc   + dx1(i)*mua   )/dmx1(i)
                    mu2 = 0.5d0*(dx1(ip)*muc   + dx1(i)*mub   )/dmx1(ip)

                    !muc = 0.5d0*(dx2(jm)*Muijk + dx2(j)*Mujm  )/dmx2(j)
                    mua = 0.5d0*(dx2(jm)*Mukm  + dx2(j)*Mujmkm)/dmx2(j)
                    mub = 0.5d0*(dx2(jm)*Mukp  + dx2(j)*Mujmkp)/dmx2(j)
                    mu5 = 0.5d0*(dx3(km)*muc   + dx3(k)*mua   )/dmx3(k )
                    mu6 = 0.5d0*(dx3(kp)*muc   + dx3(k)*mub   )/dmx3(kp)
                
                    invRhoc = 0.5d0*(dx2(jm)*invRho(i,j ,k ) + dx2(j )*invRho(i,jm,k ))/dmx2(j)
                    invRhocCmu_half= 0.5d0* Cmu*invRhoc
                
                    !M22Vn--------------------------------------------------------------------------------------------
                    !2   : Y-direction (NoBC)
                    !2-1 : Viscous             ! mu * d^2( )/dy^2
                    M22MJ = - invRhocCmu_half * (                  + Mujm / dx2(jm) ) / dmx2(j) * 2.0d0      
                    M22CJ = - invRhocCmu_half * ( - Muijk / dx2(j) - Mujm / dx2(jm) ) / dmx2(j) * 2.0d0    
                    M22PJ = - invRhocCmu_half * ( + Muijk / dx2(j)                  ) / dmx2(j) * 2.0d0  
                    !2-2 : Convection             ! Vn*d( )/dy                  ! (dVn/dy)*( )
                    M22MJ = M22MJ + 0.25d0* (             - v3/dx2(jm) + 0.5d0 * (         + dvdy3 ) )
                    M22CJ = M22CJ + 0.25d0* ( - v4/dx2(j) + v3/dx2(jm) + 0.5d0 * ( + dvdy4 + dvdy3 ) )
                    M22PJ = M22PJ + 0.25d0* ( + v4/dx2(j)              + 0.5d0 * ( + dvdy4         ) )
                
                    !3   : Z-direction (NoBC)
                    !3-1 : Viscous                       ! mu * d^2( )/dz^2 
                    M22MK = - invRhocCmu_half * (                + mu5/dmx3(k) ) /dx3(k) 
                    M22CK = - invRhocCmu_half * ( - mu6/dmx3(kp) - mu5/dmx3(k) ) /dx3(k) 
                    M22PK = - invRhocCmu_half * ( + mu6/dmx3(kp)               ) /dx3(k) 
                    !3-2 : Convection             ! Wn*d( )/dz
                    M22MK = M22MK + 0.25d0* (               - w5/dmx3(k) )
                    M22CK = M22CK + 0.25d0* ( - w6/dmx3(kp) + w5/dmx3(k) )
                    M22PK = M22PK + 0.25d0* ( + w6/dmx3(kp)              )
                
                    !1   : X-direction (NoBC)
                    !1-1 : Viscous                        ! mu * d^2( )/dx^2               
                    M22MI = - invRhocCmu_half * (                + mu1/dmx1(i) ) / dx1(i)   
                    M22CI = - invRhocCmu_half * ( - mu2/dmx1(ip) - mu1/dmx1(i) ) / dx1(i)             
                    M22PI = - invRhocCmu_half * ( + mu2/dmx1(ip)               ) / dx1(i)        
                    !1-2 : Convection            ! Un*d( )/dx 
                    M22MI = M22MI + 0.25d0*(               - u1/dmx1(i) )
                    M22CI = M22CI + 0.25d0*( - u2/dmx1(ip) + u1/dmx1(i) )
                    M22PI = M22PI + 0.25d0*( + u2/dmx1(ip)              )
                
                    ! M22Vn: Z-direction
                    AMK(i,j,k) = ( M22MK * dble(kvm) * dt ) * dble(jvc)
                    ACK(i,j,k) = ( M22CK             * dt ) * dble(jvc) + 1.0d0
                    APK(i,j,k) = ( M22PK * dble(kvp) * dt ) * dble(jvc)
                    ! M22Vn: X-direction
                    AMI(j,k,i) = ( M22MI * dble(ivm) * dt ) * dble(jvc)
                    ACI(j,k,i) = ( M22CI             * dt ) * dble(jvc) + 1.0d0
                    API(j,k,i) = ( M22PI * dble(ivp) * dt ) * dble(jvc)
                    ! M22Vn: Y-direction
                    AMJ(k,i,j) = ( M22MJ * dble(jvm) * dt ) * dble(jvc)
                    ACJ(k,i,j) = ( M22CJ             * dt ) * dble(jvc) + 1.0d0
                    APJ(k,i,j) = ( M22PJ * dble(jvp) * dt ) * dble(jvc)
                
                    !4   : RHS 
                    viscous_v1  = (mu2  *dvdx2 - mu1 *dvdx1)/ dx1(i)
                    viscous_v2  = (Muijk*dvdy4 - Mujm*dvdy3)/dmx2(j)
                    viscous_v3  = (mu6  *dvdz6 - mu5 *dvdz5)/ dx3(k)
                    viscous_v21 = (mu2  *dudy2 - mu1 *dudy1)/ dx1(i)
                    viscous_v23 = (mu6  *dwdy6 - mu5 *dwdy5)/ dx3(k)
                    !4-1 : rn (* 1/dt*un are offseted in Au^n term.) (NoBC+BC)
                    RHS_ijk  = invRhocCmu_half * ( (viscous_v1 + viscous_v2 + viscous_v3) + (viscous_v21 + viscous_v2 + viscous_v23) )
                    
                    !4-2A:       ! M11Un (NoBC)
                    RHS_ijk  = RHS_ijk - ( + M22MI * dble(ivm) * Vim + M22MJ * dble(jvm) * Vjm + M22MK * dble(kvm) * Vkm  &
                                           +(M22CI                   + M22CJ                   + M22CK           ) * Vijk &
                                           + M22PI * dble(ivp) * Vip + M22PJ * dble(jvp) * Vjp + M22PK * dble(kvp) * Vkp  ) 
                                   
                    !4-2B:       ! M21Un (NoBC)
                    M21Un = + 0.25d0 * ( u1*dble(ivm) *dvdx1 + u2*dble(ivp) *dvdx2 ) - invRhocCmu_half * ( mu2* dudy2*dble(ivp)  - mu1* dudy1*dble(ivm) ) / dx1(i)
                    RHS_ijk   = RHS_ijk - M21Un
                    !4-2C:       ! M23Wn (NoBC)
                    M23Wn = + 0.25d0 * ( w5*dble(kvm) *dvdz5 + w6*dble(kvp) *dvdz6 ) - invRhocCmu_half * ( mu6* dwdy6*dble(kvp)  - mu5* dwdy5*dble(kvm) ) / dx3(k)
                    RHS_ijk   = RHS_ijk - M23Wn 
                
                    !4-3 : Pressure
                    RHS_ijk  = RHS_ijk - Cmp*invRhoc * ( p(i,j ,k) - p(i,jm,k) )/dmx2(j)
                
                    !4-4 : Buoyancy (we use Cmt*theta for the buoyancy term)
                    ! Tc = 0.5d0*(dx2(jm)*t(i,j,k) + dx2(j)*t(i,jm,k))/dmx2(j)
                    ! !THETAy = Cmt*((Tc - Thetam))**b
                    ! !THETAy = Cmt*(Tc)
                    ! RHS_ijk(i,j,k) = RHS_ijk(i,j,k)     &
                    !             + Cmt*(Tc + a12pera11*Tc**2.*DeltaT)*invRhoc
                
                    !4-5 : BC
                    ! Y direction
                    VBCbt = YMBC(i,k,2,1)
                    VBCup = YMBC(i,k,2,2)
                
                    vbc_down = - ( M22MJ ) * VBCbt
                    vbc_up   = - ( M22PJ ) * VBCup
                
                    RHS_ijk  = RHS_ijk + dble(1-jvm)*vbc_down + dble(1-jvp)*vbc_up
                
                    ! Z direction
                    VBCbt  = ZMBC(i,j ,2,1)   
                    VBCup  = ZMBC(i,j ,2,2)   
                    WBCbtjm= ZMBC(i,jm,3,1)   
                    WBCbt  = ZMBC(i,j ,3,1) 
                    WBCupjm= ZMBC(i,jm,3,2) 
                    WBCup  = ZMBC(i,j ,3,2) 
                
                    w5    = 0.5d0* ( dx2(jm)*WBCbt + dx2(j)*WBCbtjm ) /dmx2(j)
                    w6    = 0.5d0* ( dx2(jm)*WBCup + dx2(j)*WBCupjm ) /dmx2(j)
                    dwdy5 = (WBCbt - WBCbtjm) /dmx2(j)
                    dwdy6 = (WBCup - WBCupjm) /dmx2(j)
                
                    M23Wall_down = 0.25d0* ( w5 * dvdz5 ) - invRhocCmu_half * ( - mu5 * dwdy5 ) /dx3(k)
                    M23Wall_up   = 0.25d0* ( w6 * dvdz6 ) - invRhocCmu_half * ( + mu6 * dwdy6 ) /dx3(k)
                
                    vbc_down = - ( M22MK ) * VBCbt - M23Wall_down
                    vbc_up   = - ( M22PK ) * VBCup - M23Wall_up
                
                    RHS_ijk = RHS_ijk  + dble(1-kvm)*vbc_down + dble(1-kvp)*vbc_up
                
                    ! X direction
                    VBCbt  = XMBC(j ,k,2,1)   
                    VBCup  = XMBC(j ,k,2,2)   
                    UBCbtjm= XMBC(jm,k,1,1)   
                    UBCbt  = XMBC(j ,k,1,1) 
                    UBCupjm= XMBC(jm,k,1,2) 
                    UBCup  = XMBC(j ,k,1,2) 
                
                    u1    = 0.5d0* ( dx2(jm)*UBCbt + dx2(j)*UBCbtjm ) /dmx2(j)
                    u2    = 0.5d0* ( dx2(jm)*UBCup + dx2(j)*UBCupjm ) /dmx2(j)
                    dudy1 = (UBCbt - UBCbtjm)/dmx2(j)
                    dudy2 = (UBCup - UBCupjm)/dmx2(j)
                
                    M21Uall_down = 0.25d0* ( u1 * dvdx1 ) - invRhocCmu_half * ( - mu1 * dudy1 ) /dx1(i)
                    M21Uall_up   = 0.25d0* ( u2 * dvdx2 ) - invRhocCmu_half * ( + mu2 * dudy2 ) /dx1(i)
                
                    vbc_down = - ( M22MI ) * VBCbt - M21Uall_down
                    vbc_up   = - ( M22PI ) * VBCup - M21Uall_up
                
                    RHS_ijk = RHS_ijk  + dble(1-ivm)*vbc_down + dble(1-ivp)*vbc_up
                
                    !5: M21ddU
                    dUijk  = du(i ,j,k); dUjm   = du(i ,jm,k)
                    dUip   = du(ip,j,k); dUipjm = du(ip,jm,k)
                    ddudy1 = ( dUijk - dUjm  )/dmx2(j)
                    ddudy2 = ( dUip  - dUipjm)/dmx2(j)
                
                    ddu1   = ( dx2(jm)*dUijk  + dx2(j)*dUjm   ) / dmx2(j) *0.5d0
                    ddu2   = ( dx2(jm)*dUip   + dx2(j)*dUipjm ) / dmx2(j) *0.5d0
                    
                    RHS_ijk  = RHS_ijk - 0.25d0 * ( dble(ivm)*dvdx1*ddu1 + dble(ivp)*dvdx2*ddu2 ) + invRhocCmu_half * ( dble(ivp)*ddudy2*mu2 - dble(ivm)*ddudy1*mu1 ) / dx1(i)
                
                    RHS(i,j,k) = RHS_ijk * dt * dble(jvc)
                enddo
            end do
        end do

        ! 1st stage : solve TDM in z-direction
        call PaScaL_TDMA_plan_many_create(ptdma_plan, (n1msub)*(n2msub), comm_1d_x3%myrank, comm_1d_x3%nprocs, comm_1d_x3%mpi_comm)
        call PaScaL_TDMA_many_solve_cycle(ptdma_plan, AMK, ACK, APK, RHS,(n1msub)*(n2msub),(n3msub))
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x3%nprocs)

        ! Deallocate A-matrix for the 2nd stage
        deallocate(AMK, ACK, APK)

        ! Transpose the array (i,j,k) to array (j,k,i) for the tridagonal systems in x-direction
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
        call PaScaL_TDMA_plan_many_create(ptdma_plan, (n2msub)*(n3msub), comm_1d_x1%myrank, comm_1d_x1%nprocs, comm_1d_x1%mpi_comm)
        call PaScaL_TDMA_many_solve_cycle(ptdma_plan, AMI, ACI, API, RHSI,(n2msub)*(n3msub),(n1msub))
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
        call PaScaL_TDMA_plan_many_create(ptdma_plan, (n3msub)*(n1msub), comm_1d_x2%myrank, comm_1d_x2%nprocs, comm_1d_x2%mpi_comm)
        call PaScaL_TDMA_many_solve(ptdma_plan, AMJ, ACJ, APJ, RHSJ,(n3msub)*(n1msub),(n2msub))
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x2%nprocs)

        ! Deallocate A-matrix for the 3rd stage
        deallocate(AMJ,ACJ,APJ)

        ! Update the velocity increments with the transpose of the solution array (i,k,j)
        allocate(dV(0:n1sub,0:n2sub,0:n3sub))
        dV=0.0d0
        do k = 1, n3msub
        do j = 1, n2msub
            jvc = jS_BC(j)
            do i = 1, n1msub
                dV(i,j,k)=RHSJ(k,i,j) * dble(jvc)
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
        use mpi_subdomain,  only : n1sub, n2sub, n3sub, n1msub, n2msub, n3msub
        use mpi_subdomain,  only : dx1_sub,dx2_sub,dx3_sub,dmx1_sub,dmx2_sub,dmx3_sub
        use mpi_subdomain,  only : k_indexS, iC_BC, jC_BC, kS_BC

        use PaScaL_TDMA

        implicit none

        double precision, dimension(0:n1sub,0:n2sub,0:n3sub) ::  T

        type(ptdma_plan_many)     :: ptdma_plan

        !> @{ Local pointer for subdomain variables
        double precision, dimension(:), pointer :: dx1, dx2, dx3
        double precision, dimension(:), pointer :: dmx1, dmx2, dmx3
        !> @}

        !> @{ Local indexing variables
        integer :: im, i, ip, iwm, iwc, iwp
        integer :: jm, j, jp, jwm, jwc, jwp
        integer :: km, k, kp, kwm, kwc, kwp
        !> @}

        !> @{ Local arrays for matrix RHS coefficients
        double precision, allocatable, dimension(:,:,:) :: ddW, RHS, RHSI, RHSJ
        double precision, allocatable, dimension(:,:,:) :: AMK, ACK, APK, AMI, ACI, API, AMJ, ACJ, APJ
        !> @}

        !> @{ Local variables for matrix and RHS coefficients calculation
        double precision :: u1,u2,v3,v4,w5,w6,dwdx1,dwdx2,dwdy3,dwdy4,dwdz5,dwdz6
        double precision :: dudz1,dudz2,dvdz3,dvdz4,ddudz1,ddudz2,ddvdz3,ddvdz4,ddu1,ddu2,ddv3,ddv4
        double precision :: mua,muc,mub,mu1,mu2,mu3,mu4
        double precision :: viscous_w1,viscous_w2,viscous_w3,viscous_w31,viscous_w32
        double precision :: wbc_down, wbc_up
        double precision :: M33MI, M33CI, M33PI, M33MJ, M33CJ, M33PJ, M33MK, M33CK, M33PK
        double precision :: M31Un, M32Vn
        double precision :: M31Uall_down, M31Uall_up, M32Vall_down, M32Vall_up
        double precision :: invRhoc, invRhocCmu_half

        double precision :: RHS_ijk
        double precision :: Uijk, Uip, Uim, Ujp, Ujm, Ukp, Ukm, Uipkm
        double precision :: Vijk, Vip, Vim, Vjp, Vjm, Vkp, Vkm, Vjpkm
        double precision :: Wijk, Wip, Wim, Wjp, Wjm, Wkp, Wkm
        double precision :: dUijk, dUkm, dUip, dUipkm
        double precision :: dVijk, dVkm, dVjp, dVjpkm
        double precision :: Muijk, Muim, Muimkm, Mukm, Muip, Muipkm, Mujm, Mujmkm, Mujp, Mujpkm
        double precision :: Pijk, Pkm
        double precision :: UBCbt, UBCbtkm, UBCup, UBCupkm, VBCbt, VBCbtkm, VBCup, VBCupkm, WBCbt, WBCup
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
            km = max(k,k_indexS)-1 ! To prevent NaN values
            kp=k+1
            kwm = kS_BC(km); kwc = kS_BC(k); kwp = kS_BC(kp) 
            do j = 1, n2msub
                jm = j-1
                jp = j+1
                jwm = jC_BC(jm); jwc = jC_BC(j); jwp = jC_BC(jp)
                do i = 1, n1msub
                    im = i-1
                    ip = i+1
                    iwm = iC_BC(im); iwc = iC_BC(i); iwp = iC_BC(ip)

                    Wijk = w(i ,j ,k )
                    Wip  = w(ip,j ,k ); Wim  = w(im,j ,k ); 
                    Wjp  = w(i ,jp,k ); Wjm  = w(i ,jm,k );
                    Wkp  = w(i ,j ,kp); Wkm  = w(i ,j ,km);
                
                    w5 = 0.5d0*(Wijk + Wkm)
                    w6 = 0.5d0*(Wijk + Wkp)
                
                    dwdx1 = (Wijk - Wim )/dmx1(i )
                    dwdx2 = (Wip  - Wijk)/dmx1(ip)
                    dwdy3 = (Wijk - Wjm )/dmx2(j )
                    dwdy4 = (Wjp  - Wijk)/dmx2(jp)
                    dwdz5 = (Wijk - Wkm )/ dx3(km)
                    dwdz6 = (Wkp  - Wijk)/ dx3(k )
                
                    Uijk = u(i ,j,k);  Ukm    = u(i ,j,km)
                    Uip  = u(ip,j,k);  Uipkm  = u(ip,j,km)
                    u1 = 0.5d0*(dx3(km)*Uijk + dx3(k)*Ukm  )/dmx3(k)
                    u2 = 0.5d0*(dx3(km)*Uip  + dx3(k)*Uipkm)/dmx3(k)
                    dudz1 = (Uijk - Ukm  )/dmx3(k)
                    dudz2 = (Uip  - Uipkm)/dmx3(k)
                
                    Vijk = v(i,j ,k); Vkm    = v(i,j ,km)
                    Vjp  = v(i,jp,k); Vjpkm  = v(i,jp,km)
                    v3 = 0.5d0*(dx3(km)*Vijk + dx3(k)*Vkm  )/dmx3(k)
                    v4 = 0.5d0*(dx3(km)*Vjp  + dx3(k)*Vjpkm)/dmx3(k)
                    dvdz3 = (Vijk - Vkm  )/dmx3(k)
                    dvdz4 = (Vjp  - Vjpkm)/dmx3(k)
                
                    Muijk = Mu(i,j,k ); Muim   = Mu(im,j,k ); Muip   = Mu(ip,j,k ); Mujm   =Mu(i,jm,k ); Mujp   = Mu(i,jp,k )
                    Mukm  = Mu(i,j,km); Muimkm = Mu(im,j,km); Muipkm = Mu(ip,j,km); Mujmkm =Mu(i,jm,km); Mujpkm = Mu(i,jp,km)
                
                    muc = 0.5d0*(dx3(km)*Muijk + dx3(k)*Mukm  )/dmx3(k)
                    mua = 0.5d0*(dx3(km)*Muim  + dx3(k)*Muimkm)/dmx3(k)
                    mub = 0.5d0*(dx3(km)*Muip  + dx3(k)*Muipkm)/dmx3(k)
                    mu1 = 0.5d0*(dx1(im)*muc   + dx1(i)*mua   )/dmx1(i)
                    mu2 = 0.5d0*(dx1(ip)*muc   + dx1(i)*mub   )/dmx1(ip)

                    !muc = 0.5d0*(dx3(km)*Muijk + dx3(k)*Mukm  )/dmx3(k)
                    mua = 0.5d0*(dx3(km)*Mujm  + dx3(k)*Mujmkm)/dmx3(k)
                    mub = 0.5d0*(dx3(km)*Mujp  + dx3(k)*Mujpkm)/dmx3(k)
                    mu3 = 0.5d0*(dx2(jm)*muc   + dx2(j)*mua   )/dmx2(j )
                    mu4 = 0.5d0*(dx2(jp)*muc   + dx2(j)*mub   )/dmx2(jp)
                
                    invRhoc = 0.5d0*(dx3(km)*invRho(i,j,k ) + dx3(k )*invRho(i,j,km) )/dmx3(k)
                    invRhocCmu_half = 0.5d0* Cmu*invRhoc
                
                    !M33Wn--------------------------------------------------------------------------------------------
                    !3   : Z-direction (NoBC)
                    !3-1 : Viscous             ! mu * d^2( )/dz^2
                    M33MK = - invRhocCmu_half * (                  + Mukm / dx3(km) ) /dmx3(k) *2.0d0 
                    M33CK = - invRhocCmu_half * ( - Muijk / dx3(k) - Mukm / dx3(km) ) /dmx3(k) *2.0d0 
                    M33PK = - invRhocCmu_half * ( + Muijk / dx3(k)                  ) /dmx3(k) *2.0d0 
                    !3-2 : Convection              ! Wn*d( )/dz                 ! (dWn/dz)*( )         
                    M33MK = M33MK + 0.25d0* (             - w5/dx3(km) + 0.5d0* (         + dwdz5 ))
                    M33CK = M33CK + 0.25d0* ( - w6/dx3(k) + w5/dx3(km) + 0.5d0* ( + dwdz6 + dwdz5 ))
                    M33PK = M33PK + 0.25d0* ( + w6/dx3(k)              + 0.5d0* ( + dwdz6         ))
                
                    !1   : X-direction (NoBC)
                    !1-1 : Viscous              ! mu * d^2( )/dx^2           
                    M33MI = - invRhocCmu_half * (                + mu1/dmx1(i) ) /dx1(i) 
                    M33CI = - invRhocCmu_half * ( - mu2/dmx1(ip) - mu1/dmx1(i) ) /dx1(i) 
                    M33PI = - invRhocCmu_half * ( + mu2/dmx1(ip)               ) /dx1(i) 
                    !1-2 : Convection              ! Un*d( )/dx 
                    M33MI = M33MI + 0.25d0 * (               - u1/dmx1(i) )
                    M33CI = M33CI + 0.25d0 * ( - u2/dmx1(ip) + u1/dmx1(i) )
                    M33PI = M33PI + 0.25d0 * ( + u2/dmx1(ip)              )
                
                    !2   : Y-direction (NoBC)
                    !2-1 : Viscous                               ! mu * d^2( )/dy^2 
                    M33MJ = - invRhocCmu_half *(                + mu3/dmx2(j) ) /dx2(j)
                    M33CJ = - invRhocCmu_half *( - mu4/dmx2(jp) - mu3/dmx2(j) ) /dx2(j)
                    M33PJ = - invRhocCmu_half *( + mu4/dmx2(jp)               ) /dx2(j)
                    !2-2 : Convection            ! Vn*d( )/dy  
                    M33MJ = M33MJ + 0.25d0* (               - v3/dmx2(j) )
                    M33CJ = M33CJ + 0.25d0* ( - v4/dmx2(jp) + v3/dmx2(j) ) 
                    M33PJ = M33PJ + 0.25d0* ( + v4/dmx2(jp)              )
                
                    !M33Wn: Z-direction
                    AMK(i,j,k) = ( M33MK * dble(kwm) * dt ) * dble(kwc)
                    ACK(i,j,k) = ( M33CK             * dt ) * dble(kwc) + 1.0d0
                    APK(i,j,k) = ( M33PK * dble(kwp) * dt ) * dble(kwc)
                    !M33Wn: X-direction
                    AMI(j,k,i) = ( M33MI * dble(iwm) * dt ) * dble(kwc)
                    ACI(j,k,i) = ( M33CI             * dt ) * dble(kwc) + 1.0d0
                    API(j,k,i) = ( M33PI * dble(iwp) * dt ) * dble(kwc)
                    !M33Wn: Y-direction
                    AMJ(k,i,j) = ( M33MJ * dble(jwm) * dt ) * dble(kwc)
                    ACJ(k,i,j) = ( M33CJ             * dt ) * dble(kwc) + 1.0d0
                    APJ(k,i,j) = ( M33PJ * dble(jwp) * dt ) * dble(kwc)
                    
                    !4   : RHS
                    viscous_w1  = (mu2  *dwdx2 - mu1 *dwdx1)/ dx1(i)
                    viscous_w2  = (mu4  *dwdy4 - mu3 *dwdy3)/ dx2(j)
                    viscous_w3  = (Muijk*dwdz6 - Mukm*dwdz5)/dmx3(k)
                    viscous_w31 = (mu2  *dudz2 - mu1 *dudz1)/ dx1(i)
                    viscous_w32 = (mu4  *dvdz4 - mu3 *dvdz3)/ dx2(j)
                    !4-1 : rn (* 1/dt*un are offseted in Au^n term.) (NoBC+BC)    
                    RHS_ijk  =  invRhocCmu_half *( (viscous_w1 + viscous_w2 + viscous_w3) + (viscous_w31 + viscous_w32 + viscous_w3) )
                
                    !4-2A:       ! M33Wn (NoBC)
                    RHS_ijk  = RHS_ijk - ( + M33MI * dble(iwm) * Wim + M33MJ * dble(jwm) * Wjm + M33MK * dble(kwm) * Wkm  &
                                           +(M33CI                   + M33CJ                      + M33CK              ) * Wijk &
                                           + M33PI * dble(iwp) * Wip + M33PJ * dble(jwp) * Wjp + M33PK * dble(kwp) * Wkp  ) 
                    !4-2B:       ! M31Un (NoBC)
                    M31Un = + 0.25d0* ( + u1*dble(iwm) *dwdx1 + u2*dble(iwp) *dwdx2 ) - invRhocCmu_half *( mu2*dudz2*dble(iwp) - mu1*dudz1*dble(iwm) ) / dx1(i)
                    RHS_ijk   = RHS_ijk - M31Un
                    !4-2C:       ! M32Vn (NoBC)
                    M32Vn = + 0.25d0* ( + v3*dble(jwm) *dwdy3 + v4*dble(jwp) *dwdy4 ) - invRhocCmu_half *( mu4*dvdz4*dble(jwp) - mu3*dvdz3*dble(jwm) ) / dx2(j)
                    RHS_ijk   = RHS_ijk - M32Vn
                
                    !4-3 : Pressure
                    Pijk = p(i,j,k )
                    Pkm  = p(i,j,km)
                    RHS_ijk  = RHS_ijk - Cmp*invRhoc * ( Pijk - Pkm ) / dmx3(k)
                
                    !4-4 : Buoyancy term
                
                    !4-5 : BC
                    ! Z direction
                    WBCbt = ZMBC(i,j,3,1)
                    WBCup = ZMBC(i,j,3,2)
                
                    wbc_down = - ( M33MK ) * WBCbt
                    wbc_up   = - ( M33PK ) * WBCup
                
                    RHS_ijk  = RHS_ijk + dble(1-kwm)*wbc_down + dble(1-kwp)*wbc_up
                
                    ! X direction
                    WBCbt  =XMBC(j,k ,3,1)
                    WBCup  =XMBC(j,k ,3,2)
                    UBCbt  =XMBC(j,k ,1,1)
                    UBCbtkm=XMBC(j,km,1,1) 
                    UBCup  =XMBC(j,k ,1,2)
                    UBCupkm=XMBC(j,km,1,2)
                
                    u1 = 0.5d0* ( dx3(km)*UBCbt + dx3(k)*UBCbtkm ) /dmx3(k)
                    u2 = 0.5d0* ( dx3(km)*UBCup + dx3(k)*UBCupkm ) /dmx3(k)
                    dudz1 = (UBCbt - UBCbtkm)/dmx3(k)
                    dudz2 = (UBCup - UBCupkm)/dmx3(k)
                
                    M31Uall_down = + 0.25d0* ( u1 * dwdx1 ) - invRhocCmu_half * ( - mu1 * dudz1 ) /dx1(i)
                    M31Uall_up   = + 0.25d0* ( u2 * dwdx2 ) - invRhocCmu_half * ( + mu2 * dudz2 ) /dx1(i)
                
                    wbc_down = - ( M33MI ) * WBCbt - M31Uall_down 
                    wbc_up   = - ( M33PI ) * WBCup - M31Uall_up 
                
                    RHS_ijk  = RHS_ijk + dble(1-iwm)*wbc_down + dble(1-iwp)*wbc_up
                
                    ! Y direction
                    WBCbt  =YMBC(i,k ,3,1)
                    WBCup  =YMBC(i,k ,3,2)
                    VBCbt  =YMBC(i,k ,2,1)
                    VBCbtkm=YMBC(i,km,2,1) 
                    VBCup  =YMBC(i,k ,2,2)
                    VBCupkm=YMBC(i,km,2,2)
                    
                    v3    = (dx3(k)*VBCbtkm + dx3(km)*VBCbt)*0.5d0 /dmx3(k)
                    v4    = (dx3(k)*VBCupkm + dx3(km)*VBCup)*0.5d0 /dmx3(k)
                    dvdz3 = (VBCbt - VBCbtkm)/dmx3(k)
                    dvdz4 = (VBCup - VBCupkm)/dmx3(k)
                
                    M32Vall_down = + 0.25d0 * ( v3 * dwdy3 ) - invRhocCmu_half * ( - mu3 * dvdz3 ) /dx2(j)
                    M32Vall_up   = + 0.25d0 * ( v4 * dwdy4 ) - invRhocCmu_half * ( + mu4 * dvdz4 ) /dx2(j)
                
                    wbc_down = - ( M33MJ ) * WBCbt - M32Vall_down 
                    wbc_up   = - ( M33PJ ) * WBCup - M32Vall_up 
                
                    RHS_ijk  = RHS_ijk + dble(1-jwm)*wbc_down + dble(1-jwp)*wbc_up
                
                    !5-1 : M31ddU
                    dUijk  = du(i ,j,k); dUkm   = du(i ,j,km)
                    dUip   = du(ip,j,k); dUipkm = du(ip,j,km)
                
                    ddudz1 = (dUijk - dUkm  )/dmx3(k)
                    ddudz2 = (dUip  - dUipkm)/dmx3(k) 
                
                    ddu1 = (dx3(km)*dUijk + dx3(k)*dUkm  )/dmx3(k)*0.5d0
                    ddu2 = (dx3(km)*dUip  + dx3(k)*dUipkm)/dmx3(k)*0.5d0
                
                    RHS_ijk  = RHS_ijk - 0.25d0* ( ddu1*dble(iwm) *dwdx1 + ddu2*dble(iwp) *dwdx2 ) + invRhocCmu_half *( mu2* ddudz2*dble(iwp) - mu1* ddudz1*dble(iwm) ) /dx1(i)
                    
                    !5-2 : M32ddV
                    dVijk  = dv(i,j ,k); dVkm   = dv(i,j ,km)
                    dVjp   = dv(i,jp,k); dVjpkm = dv(i,jp,km)
                
                    ddvdz3 = (dVijk - dVkm  )/dmx3(k)
                    ddvdz4 = (dVjp  - dVjpkm)/dmx3(k)
                
                    ddv3 = (dx3(km)*dVijk + dx3(k)*dVkm  )/dmx3(k)*0.5d0
                    ddv4 = (dx3(km)*dVjp  + dx3(k)*dVjpkm)/dmx3(k)*0.5d0
                    
                    RHS_ijk  = RHS_ijk - 0.25d0* ( ddv3*dble(jwm) *dwdy3 + ddv4*dble(jwp) *dwdy4 ) + invRhocCmu_half *( mu4* ddvdz4*dble(jwp) - mu3* ddvdz3*dble(jwm) ) /dx2(j)
                    
                    RHS(i,j,k) = ( RHS_ijk  * dt ) * dble(kwc)
                enddo
            enddo  
        enddo

        ! 1st stage : solve TDM in z-direction
        call PaScaL_TDMA_plan_many_create(ptdma_plan, (n1msub)*(n2msub), comm_1d_x3%myrank, comm_1d_x3%nprocs, comm_1d_x3%mpi_comm)
        call PaScaL_TDMA_many_solve_cycle(ptdma_plan, AMK, ACK, APK, RHS,(n1msub)*(n2msub),(n3msub))
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x3%nprocs)

        ! Deallocate A-matrix for the 2nd stage
        deallocate(AMK,ACK,APK)

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
        call PaScaL_TDMA_plan_many_create(ptdma_plan, (n2msub)*(n3msub), comm_1d_x1%myrank, comm_1d_x1%nprocs, comm_1d_x1%mpi_comm)
        call PaScaL_TDMA_many_solve_cycle(ptdma_plan, AMI, ACI, API, RHSI,(n2msub)*(n3msub),(n1msub))
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
        call PaScaL_TDMA_plan_many_create(ptdma_plan, (n3msub)*(n1msub), comm_1d_x2%myrank, comm_1d_x2%nprocs, comm_1d_x2%mpi_comm)
        call PaScaL_TDMA_many_solve(ptdma_plan, AMJ, ACJ, APJ, RHSJ,(n3msub)*(n1msub),(n2msub))
        call PaScaL_TDMA_plan_many_destroy(ptdma_plan,comm_1d_x2%nprocs)        

        ! Deallocate A-matrix for the 3rd stage
        deallocate(AMJ,ACJ,APJ)

        ! Update the velocity increments with the transpose of the solution array (i,k,j)
        allocate(dW(0:n1sub,0:n2sub,0:n3sub))
        dW=0.0d0
        do k = 1, n3msub
            kwc = kS_BC(k)
        do j = 1, n2msub
        do i = 1, n1msub
            dW(i,j,k)=RHSJ(k,i,j) * dble(kwc)
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

        use mpi_subdomain,  only : j_indexS, iC_BC, jS_BC, kC_BC
        use mpi_subdomain,  only : n1sub, n2sub, n3sub, n1msub, n2msub, n3msub
        use mpi_subdomain,  only : dx2_sub,dx3_sub,dmx2_sub,dmx3_sub

        implicit none

        double precision, dimension(0:n1sub,0:n2sub,0:n3sub) ::  T

        !> @{ Local pointer for subdomain variables
        double precision, dimension(:), pointer :: dx2, dx3
        double precision, dimension(:), pointer :: dmx2, dmx3
        !> @}

        !> @{ Local indexing variables
        integer :: im, i,  ip
        integer :: jm, j,  jp
        integer :: km, k,  kp
        integer :: ivm,ivc,ivp
        integer :: jvm,jvc,jvp
        integer :: kvm,kvc,kvp
        !> @}

        !> @{ Local variables for dV update
        double precision :: dwm5,dwm6
        double precision :: dvdz5,dvdz6,ddwdy5,ddwdy6
        double precision :: mua,muc,mub,mu5,mu6
        double precision :: invRhoc, invRhocCmu_half
        double precision :: M23dWm
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
            kvm = kC_BC(k-1); kvc = kC_BC(k); kvp = kC_BC(k+1)
            do j =  1, n2msub
                jp = j + 1
                jm = max(j,j_indexS)-1 ! To prevent NaN values
                jvm = jS_BC(j-1); jvc = jS_BC(j); jvp = jS_BC(j+1)
                do i = 1, n1msub

                    dwm5 = (dx2(jm)*dw(i,j,k ) + dx2(j)*dw(i,jm,k ))/dmx2(j)*0.5d0
                    dwm6 = (dx2(jm)*dw(i,j,kp) + dx2(j)*dw(i,jm,kp))/dmx2(j)*0.5d0
                
                    ddwdy5 = (dw(i,j,k ) - dw(i,jm,k ))/dmx2(j)
                    ddwdy6 = (dw(i,j,kp) - dw(i,jm,kp))/dmx2(j)  
                
                    dvdz5 = (v(i,j,k ) - v(i,j,km))/dmx3(k )
                    dvdz6 = (v(i,j,kp) - v(i,j,k ))/dmx3(kp)
                
                    muc = 0.5d0*(dx2(jm)*Mu(i,j ,k ) + dx2(j)*Mu(i,jm,k ))/dmx2(j)
                    mua = 0.5d0*(dx2(jm)*Mu(i,j ,km) + dx2(j)*Mu(i,jm,km))/dmx2(j)
                    mub = 0.5d0*(dx2(jm)*Mu(i,j ,kp) + dx2(j)*Mu(i,jm,kp))/dmx2(j)
                    mu5 = 0.5d0*(dx3(km)*muc           + dx3(k)*mua          )/dmx3(k )
                    mu6 = 0.5d0*(dx3(kp)*muc           + dx3(k)*mub          )/dmx3(kp)
                
                    invRhoc = 0.5d0*(dx2(jm)*invRho(i,j,k ) + dx2(j )*invRho(i,jm,k ))/dmx2(j)
                    invRhocCmu_half = 0.5d0* Cmu*invRhoc
 
                    M23dWm = + dble(0.25) * ( dwm5*dble(kvm) *dvdz5 + dwm6*dble(kvp) *dvdz6 ) - invRhocCmu_half * ( mu6* ddwdy6*dble(kvp)  - mu5* ddwdy5*dble(kvm) ) / dx3(k)
                    dv(i,j,k) = dv(i,j,k) - dt * ( M23dWm ) * dble(jvc)
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

        use mpi_subdomain,  only : i_indexS, iS_BC, jC_BC, kC_BC
        use mpi_subdomain,  only : n1sub, n2sub, n3sub, n1msub, n2msub, n3msub
        use mpi_subdomain,  only : dx1_sub,dx2_sub,dx3_sub,dmx1_sub,dmx2_sub,dmx3_sub

        implicit none

        double precision, dimension(0:n1sub,0:n2sub,0:n3sub) ::  T

        !> @{ Local pointer for subdomain variables
        double precision, dimension(:), pointer :: dx1, dx2, dx3
        double precision, dimension(:), pointer :: dmx1, dmx2, dmx3
        !> @}
        
        !> @{ Local indexing variables
        integer :: im, i,  ip
        integer :: jm, j,  jp
        integer :: km, k,  kp
        integer :: ium,iuc,iup
        integer :: jum,juc,jup
        integer :: kum,kuc,kup
        !> @}

        !> @{ Local variables for dV update
        double precision :: dvm3,dvm4
        double precision :: dudy3,dudy4,dvmdx3,dvmdx4
        double precision :: mua,mub,muc,mu3,mu4
        double precision :: invRhoc, invRhocCmu_half
    
        double precision :: dwm5,dwm6
        double precision :: dudz5,dudz6,dwmdx5,dwmdx6
        double precision :: mu5,mu6
        double precision :: M12dVm, M13dWm
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
            km = k-1;  kp = k+1;
            kum = kC_BC(k-1); kuc = kC_BC(k); kup = kC_BC(k+1)
            do j = 1, n2msub
                jm = j-1;  jp = j+1;
                jum = jC_BC(j-1); juc = jC_BC(j); jup = jC_BC(j+1)
                do i = 1, n1msub
                    im = max(i,i_indexS)-1 ! To prevent NaN values
                    ip = i+1;
                    ium = iS_BC(i-1); iuc = iS_BC(i); iup = iS_BC(i+1)
               
                    invRhoc = 0.5d0*(dx1(im)*invRho(i,j ,k ) + dx1(i)*invRho(im,j ,k ) )/dmx1(i)
                    invRhocCmu_half = 0.5d0* Cmu*invRhoc
                
                    dvm3 = ( dx1(im)*dv(i,j, k) + dx1(i)*dv(im,j, k) )/dmx1(i)*0.5d0
                    dvm4 = ( dx1(im)*dv(i,jp,k) + dx1(i)*dv(im,jp,k) )/dmx1(i)*0.5d0
                
                    dudy3 = (u(i,j, k) - u(i,jm,k))/dmx2(j )
                    dudy4 = (u(i,jp,k) - u(i,j, k))/dmx2(jp)
                
                    dvmdx3 = (dv(i,j ,k) - dv(im,j ,k))/dmx1(i)
                    dvmdx4 = (dv(i,jp,k) - dv(im,jp,k))/dmx1(i)
                
                    mua = 0.5d0*(dx1(im)*Mu(i,jm,k ) + dx1(i)*Mu(im,jm,k ))/dmx1(i)
                    muc = 0.5d0*(dx1(im)*Mu(i,j ,k ) + dx1(i)*Mu(im,j ,k ))/dmx1(i)
                    mub = 0.5d0*(dx1(im)*Mu(i,jp,k ) + dx1(i)*Mu(im,jp,k ))/dmx1(i)
                    mu3 = 0.5d0*(dx2(jm)*muc         + dx2(j)*mua         )/dmx2(j )
                    mu4 = 0.5d0*(dx2(jp)*muc         + dx2(j)*mub         )/dmx2(jp)
                
                    dwm5 = (dx1(i)*dw(im,j,k ) + dx1(im)*dw(i,j,k ))/dmx1(i)*0.5d0 
                    dwm6 = (dx1(i)*dw(im,j,kp) + dx1(im)*dw(i,j,kp))/dmx1(i)*0.5d0 
                
                    dudz5 = (u(i,j,k ) - u(i,j,km))/dmx3(k )
                    dudz6 = (u(i,j,kp) - u(i,j,k ))/dmx3(kp)
                
                    dwmdx5 = (dw(i,j,k ) - dw(im,j,k ))/dmx1(i)
                    dwmdx6 = (dw(i,j,kp) - dw(im,j,kp))/dmx1(i)
                
                    muc = 0.5d0*(dx1(im)*Mu(i,j ,k ) + dx1(i)*Mu(im,j ,k ))/dmx1(i)
                    mua = 0.5d0*(dx1(im)*Mu(i,j ,km) + dx1(i)*Mu(im,j ,km))/dmx1(i)
                    mub = 0.5d0*(dx1(im)*Mu(i,j ,kp) + dx1(i)*Mu(im,j ,kp))/dmx1(i)            
                    mu5 = 0.5d0*(dx3(km)*muc         + dx3(k)*mua         )/dmx3(k )
                    mu6 = 0.5d0*(dx3(kp)*muc         + dx3(k)*mub         )/dmx3(kp)
                
                    !>  dU(i,j,k) = dU(i,j,k) - dt*M12dVm - dt*M13dWm
                    M12dVm = 0.25d0*( dvm3*dble(jum) *dudy3 + dvm4*dble(jup) *dudy4 ) - invRhocCmu_half *( mu4* dvmdx4*dble(jup) - mu3* dvmdx3*dble(jum) )/dx2(j)
                    M13dWm = 0.25d0*( dwm5*dble(kum) *dudz5 + dwm6*dble(kup) *dudz6 ) - invRhocCmu_half *( mu6* dwmdx6*dble(kup) - mu5* dwmdx5*dble(kum) )/dx3(k)
                    du(i,j,k) = du(i,j,k) - dt * ( M12dVm ) * dble(iuc) &
                                          - dt * ( M13dWm ) * dble(iuc)
                enddo
            enddo
        enddo

        ! Nullify grid information pointer
        nullify(dx1, dx2, dx3, dmx1, dmx2, dmx3)

    end subroutine mpi_momentum_blockLdU

    subroutine mpi_momentum_masscorrection(dU, dP) !cell
        
        use MPI
        use mpi_topology,   only : comm_1d_x1, comm_1d_x2, comm_1d_x3
        use mpi_subdomain,  only : dx1_sub, dx2_sub, dx3_sub, dmx1_sub, dmx2_sub, dmx3_sub
        use mpi_subdomain,  only : n1sub, n2sub, n3sub, n1msub, n2msub, n3msub

        implicit none

        double precision, dimension(:), pointer :: dx1, dx2, dx3
        double precision, dimension(:), pointer :: dmx1, dmx2, dmx3
        double precision, dimension(0:n1sub,0:n2sub,0:n3sub), intent(in) :: dU, dP

        integer :: i,j,k
        integer :: myrank
        ! !integer :: cellijk, celljp, cellkp, celljkp
        
        double precision :: flow1, flow1_total
        double precision :: DMpresg, DMpresg_total
        double precision, dimension(1:2) :: package, package_I, package_K, package_total
        integer :: ierr

        ! ! Pointer of grid information

        dx1  => dx1_sub
        dx2  => dx2_sub
        dx3  => dx3_sub
        dmx1 => dmx1_sub
        dmx2 => dmx2_sub
        dmx3 => dmx3_sub

        flow1=0.d0; DMpresg=0.d0
        flow1_total=0.d0; DMpresg_total=0.d0
        package=0.d0; package_I=0.d0; package_K=0.d0; package_total=0.d0

        do k=1,n3msub
        do j=1,n2msub
        do i=1,n1msub
            flow1   =   flow1 + (U(i,j,k)-dU(i,j,k))*dx2(j)*dmx1(i)*dx3(k)
            DMpresg = DMpresg + U(i,j,k)*dx2(j)*dmx1(i)*dx3(k) - dt*(dP(i,j,k)-dP(i-1,j,k))*dx2(j)*dx3(k)
        enddo
        enddo
        enddo

        package(1)=flow1
        package(2)=DMpresg
        call MPI_Allreduce(package  , package_I    ,2,MPI_DOUBLE_PRECISION,MPI_SUM,comm_1d_x1%mpi_comm, ierr)
        call MPI_Allreduce(package_I, package_K    ,2,MPI_DOUBLE_PRECISION,MPI_SUM,comm_1d_x3%mpi_comm, ierr)
        call MPI_Allreduce(package_K, package_total,2,MPI_DOUBLE_PRECISION,MPI_SUM,comm_1d_x2%mpi_comm, ierr)

        flow1_total   = package_total(1)
        DMpresg_total = package_total(2)
        DMpresg_total = (DMpresg_total - flow1_total)/Volume/dt
        presgrad1     = presgrad1 + DMpresg_total

        do k = 1, n3msub
        do j = 1, n2msub
        do i = 1, n1msub
            U(i,j,k) = U(i,j,k) - dt*DMpresg_total
        enddo
        enddo
        enddo

    end subroutine mpi_momentum_masscorrection

end module mpi_momentum