module mpi_statistics
    ! use debug
    use global

    implicit none

    double precision :: avg_time_length
	double precision :: u_tau, delta_nu, Temp_tau ! YMG: Non-dimensionlized Parameter Calculation

    ! Output File =======
    double precision, allocatable, dimension(:) :: Umean, Vmean, Wmean, Tmean, Pmean ! mean profile avged in x- and z- dir.
    double precision, allocatable, dimension(:) :: UUavg_t, VVavg_t, WWavg_t, TTavg_t   ! ( fluctuation component )**2
    double precision, allocatable, dimension(:) :: UVavg_t, UTavg_t, UWavg_t, VWavg_t

    ! For 1D-Turbulence Statistics =======
    double precision, allocatable, dimension(:) :: UUUavg_t, UUVavg_t, UUWavg_t ! YMG: For Turbulent Diffusion
	double precision, allocatable, dimension(:) :: VVUavg_t, VVVavg_t, VVWavg_t
	double precision, allocatable, dimension(:) :: WWUavg_t, WWVavg_t, WWWavg_t
	double precision, allocatable, dimension(:) :: UVUavg_t, UVVavg_t, UVWavg_t
	
	double precision, allocatable, dimension(:) :: dudxdudx, dudydudy, dudzdudz ! YMG: For Dissipation
	double precision, allocatable, dimension(:) :: dvdxdvdx, dvdydvdy, dvdzdvdz 
	double precision, allocatable, dimension(:) :: dwdxdwdx, dwdydwdy, dwdzdwdz 
	double precision, allocatable, dimension(:) :: dudxdvdx, dudydvdy, dudzdvdz 
	
	double precision, allocatable, dimension(:) :: udpdx, vdpdy, wdpdz, udpdy, vdpdx ! YMG: For V-P Gradient

    !! For 3D-Turbulence Statistics
    ! double precision, allocatable, dimension(:,:,:,:) :: UU_mean3D, dUdxdUdx_mean3D, UUU_mean3D, Udpdx_mean3D ! YMG Reynolds - Budget
    ! double precision, allocatable, dimension(:,:,:) :: Ugmean, Vgmean, Wgmean, Pgmean

    double precision, allocatable, dimension(:) :: Usq, Vsq, Wsq, Tsq, UV, UT ! (instant velocity )**2
    double precision, allocatable, dimension(:) :: Mumean

    character(len=17) :: xyzrank

contains

    subroutine mpi_statistics_xyzrank_allocation()

        use mpi_topology,   only : comm_1d_x1, comm_1d_x2, comm_1d_x3
        
        implicit none
        
        integer :: r1_0,r2_0,r3_0
        integer :: r1_1,r2_1,r3_1
        integer :: r1_2,r2_2,r3_2
        integer :: r1_3,r2_3,r3_3
        integer :: r1_4,r2_4,r3_4
        
        r1_4=int((comm_1d_x1%myrank)/10000)
        r1_3=int((comm_1d_x1%myrank-r1_4*10000)/1000 )
        r1_2=int((comm_1d_x1%myrank-r1_4*10000-r1_3*1000)/100 )
        r1_1=int((comm_1d_x1%myrank-r1_4*10000-r1_3*1000-r1_2*100)/10 )
        r1_0=int((comm_1d_x1%myrank-r1_4*10000-r1_3*1000-r1_2*100-r1_1*10)/1)

        r2_4=int((comm_1d_x2%myrank)/10000)
        r2_3=int((comm_1d_x2%myrank-r2_4*10000)/1000 )
        r2_2=int((comm_1d_x2%myrank-r2_4*10000-r2_3*1000)/100 )
        r2_1=int((comm_1d_x2%myrank-r2_4*10000-r2_3*1000-r2_2*100)/10 )
        r2_0=int((comm_1d_x2%myrank-r2_4*10000-r2_3*1000-r2_2*100-r2_1*10)/1)
        
        r3_4=int((comm_1d_x3%myrank)/10000)
        r3_3=int((comm_1d_x3%myrank-r3_4*10000)/1000 )
        r3_2=int((comm_1d_x3%myrank-r3_4*10000-r3_3*1000)/100 )
        r3_1=int((comm_1d_x3%myrank-r3_4*10000-r3_3*1000-r3_2*100)/10 )
        r3_0=int((comm_1d_x3%myrank-r3_4*10000-r3_3*1000-r3_2*100-r3_1*10)/1)
        
        write(xyzrank, '(5I1,1A1,5I1,1A1,5I1)' ) r1_4,r1_3,r1_2,r1_1,r1_0,'_'   &
                                                 ,r2_4,r2_3,r2_2,r2_1,r2_0,'_'   &
                                                 ,r3_4,r3_3,r3_2,r3_1,r3_0

    end subroutine mpi_statistics_xyzrank_allocation

    subroutine mpi_statistics_allocation()

        use mpi_subdomain,  only : n1sub, n2sub, n3sub

        implicit none

        ! Output File =======
        allocate( Umean(0:n2sub), Tmean(0:n2sub) )
        allocate( Vmean(0:n2sub), Wmean(0:n2sub), Pmean(0:n2sub) )
        allocate( Mumean(0:n2sub) )
        allocate( UUavg_t(0:n2sub), VVavg_t(0:n2sub) )
        allocate( WWavg_t(0:n2sub), TTavg_t(0:n2sub) )
		allocate( UVavg_t(0:n2sub), UTavg_t(0:n2sub) )
		allocate( UWavg_t(0:n2sub), VWavg_t(0:n2sub) )


        ! For 1D-Turbulence Statistics =======
        allocate( UUUavg_t(0:n2sub), UUVavg_t(0:n2sub), UUWavg_t(0:n2sub) ) ! YMG: For Turbulent Diffusion
		allocate( VVUavg_t(0:n2sub), VVVavg_t(0:n2sub), VVWavg_t(0:n2sub) )
		allocate( WWUavg_t(0:n2sub), WWVavg_t(0:n2sub), WWWavg_t(0:n2sub) ) 
		allocate( UVUavg_t(0:n2sub), UVVavg_t(0:n2sub), UVWavg_t(0:n2sub) )
		
		allocate( dudxdudx(0:n2sub), dudydudy(0:n2sub), dudzdudz(0:n2sub) ) ! YMG: For Dissipation
		allocate( dvdxdvdx(0:n2sub), dvdydvdy(0:n2sub), dvdzdvdz(0:n2sub) )
		allocate( dwdxdwdx(0:n2sub), dwdydwdy(0:n2sub), dwdzdwdz(0:n2sub) ) 
		allocate( dudxdvdx(0:n2sub), dudydvdy(0:n2sub), dudzdvdz(0:n2sub) )
	
		allocate( udpdx(0:n2sub), vdpdy(0:n2sub), wdpdz(0:n2sub), udpdy(0:n2sub), vdpdx(0:n2sub) )  ! YMG: For V-P Gradient

        ! ! For 3D-Turbulence Statistics
        ! allocate( UU_mean3D(6,0:s1sb,0:s2sb,0:s3sb), dUdxdUdx_mean3D(12,0:s1sb,0:s2sb,0:s3sb), UUU_mean3D(10,0:s1sb,0:s2sb,0:s3sb), Udpdx_mean3D(4,0:s1sb,0:s2sb,0:s3sb) )  ! YMG Reynolds - Budget
        ! allocate( Ugmean(0:s1sb,0:s2sb,0:s3sb), Vgmean(0:s1sb,0:s2sb,0:s3sb), Wgmean(0:s1sb,0:s2sb,0:s3sb), Pgmean(0:s1sb,0:s2sb,0:s3sb))

        ! Instant velocity **2
        allocate(Usq(0:n2sub), Vsq(0:n2sub), Wsq(0:n2sub), Tsq(0:n2sub) )
        allocate(UV(0:n2sub), UT(0:n2sub) )

        avg_time_length = 0.0d0

        ! Output File =======
        Umean=0.0d0;   Tmean=0.0d0
        Vmean=0.0d0;   Wmean=0.0d0;   Pmean=0.0d0
		
        UUavg_t=0.0d0; VVavg_t=0.0d0
        WWavg_t=0.0d0; TTavg_t=0.0d0
        UVavg_t=0.0d0; UTavg_t=0.0d0
		UWavg_t=0.0d0; VWavg_t=0.0d0

        ! For 1D-Turbulence Statistics =======
        UUUavg_t=0.0d0; UUVavg_t=0.0d0; UUWavg_t=0.0d0 ! YMG: For Turbulent Diffusion
		VVUavg_t=0.0d0; VVVavg_t=0.0d0; VVWavg_t=0.0d0
		WWUavg_t=0.0d0; WWVavg_t=0.0d0; WWWavg_t=0.0d0
		UVUavg_t=0.0d0; UVVavg_t=0.0d0; UVWavg_t=0.0d0
		
		dudxdudx=0.0d0; dudydudy=0.0d0; dudzdudz=0.0d0  ! YMG: For Dissipation
		dvdxdvdx=0.0d0; dvdydvdy=0.0d0; dvdzdvdz=0.0d0  
		dwdxdwdx=0.0d0; dwdydwdy=0.0d0; dwdzdwdz=0.0d0  
		dudxdvdx=0.0d0; dudydvdy=0.0d0; dudzdvdz=0.0d0 
		
		udpdx=0.0d0; vdpdy=0.0d0; wdpdz=0.0d0; udpdy=0.0d0; vdpdx=0.0d0  ! YMG: For V-P Gradient


        !! For 3D-Turbulence Statistics
        ! UU_mean3D=0.0d0; dUdxdUdx_mean3D=0.0d0; UUU_mean3D=0.0d0; Udpdx_mean3D=0.0d0
        ! Ugmean=0.0d0; Vgmean=0.0d0; Wgmean=0.0d0; Pgmean=0.0d0
        ! For device memory
        ! UU_mean3D=UU_mean3D; dUdxdUdx_mean3D=dUdxdUdx_mean3D; UUU_mean3D=UUU_mean3D; Udpdx_mean3D=Udpdx_mean3D
        ! Ugmean=Ugmean; Vgmean=Vgmean; Wgmean=Wgmean; Pgmean=Pgmean
		
        Usq=0.0d0; Vsq=0.0d0; Wsq=0.0d0; Tsq=0.0d0
        UV=0.0d0; UT=0.0d0
        Mumean=0.0d0

		u_tau=0.0d0; delta_nu=0.0d0; Temp_tau=0.0d0

    end subroutine mpi_statistics_allocation

    subroutine mpi_statistics_clean()
        implicit none
         ! Output File =======
        deallocate( Umean, Tmean )
        deallocate( Vmean, Wmean, Pmean )
		
        deallocate( UUavg_t, VVavg_t )
        deallocate( WWavg_t, TTavg_t )
        deallocate( UVavg_t, UTavg_t )
		deallocate( UWavg_t, VWavg_t )
      

        ! For 1D-Turbulence Statistics =======
		deallocate( UUUavg_t, UUVavg_t, UUWavg_t ) ! YMG: For Turbulent Diffusion
		deallocate( VVUavg_t, VVVavg_t, VVWavg_t )
		deallocate( WWUavg_t, WWVavg_t, WWWavg_t )
		deallocate( UVUavg_t, UVVavg_t, UVWavg_t )
		
		deallocate( dudxdudx, dudydudy, dudzdudz ) ! YMG: For Dissipation
		deallocate( dvdxdvdx, dvdydvdy, dvdzdvdz )
		deallocate( dwdxdwdx, dwdydwdy, dwdzdwdz ) 
		deallocate( dudxdvdx, dudydvdy, dudzdvdz )
		
		deallocate( udpdx, vdpdy, wdpdz, udpdy, vdpdx ) ! YMG: For V-P Gradient
      

        !! For 3D-Turbulence Statistics
        ! deallocate( UU_mean3D, dUdxdUdx_mean3D, UUU_mean3D, Udpdx_mean3D)  ! YMG Reynolds - Budget
        ! deallocate( Ugmean, Vgmean, Wgmean, Pgmean)

        deallocate( Usq, Vsq, Wsq, Tsq, UV, UT )
        deallocate( Mumean   )

    end subroutine mpi_statistics_clean

    subroutine mpi_statistics_avg_xzt(myrank,U,V,W,P,T,Mu,dtime) ! YMG 210407 !Cell
        
        use mpi_topology,   only : comm_1d_x1, comm_1d_x2, comm_1d_x3
        use mpi_subdomain,  only : dx1_sub, dx2_sub, dx3_sub, dmx1_sub, dmx2_sub, dmx3_sub
        use mpi_subdomain,  only : n1sub, n2sub, n3sub, n1msub, n2msub, n3msub
        
        implicit none
    
        integer :: i,j,k,count,myrank
        double precision, dimension(0:n1sub,0:n2sub,0:n3sub) :: U,V,W,P,T,Mu
        !integer,  dimension(-1:s1sb+1,-1:s2sb+1,-1:s3sb+1), intent(in) :: Cell
        double precision, dimension(:), pointer :: dx1,  dx2,  dx3
        double precision, dimension(:), pointer :: dmx1, dmx2, dmx3
        double precision :: dtime
        double precision :: u1,v1,w1
        double precision :: u2,v2,w2
        double precision :: t1,t2,t3,t4,t5,t6
        double precision :: p1,p2,p3,p4,p5,p6
        double precision :: ug,vg,wg,tg,pg
		double precision :: ugp, ugm, vgp, vgm, wgp, wgm, pgyp, pgym, pgxp, pgxm, pgzp, pgzm
        double precision :: temp1, temp2, temp3, temp4, temp5

        double precision :: Umean_tmp, Vmean_tmp, Wmean_tmp, Tmean_tmp, Pmean_tmp, Mumean_tmp

        double precision :: UUavg_t_tmp, VVavg_t_tmp, WWavg_t_tmp, TTavg_t_tmp   ! ( fluctuation component )**2
        double precision :: UVavg_t_tmp, UTavg_t_tmp, UWavg_t_tmp, VWavg_t_tmp
        double precision :: Usq_tmp, Vsq_tmp, Wsq_tmp, Tsq_tmp, UV_tmp, UT_tmp ! (instant velocity )**2

        double precision :: dudxdudx_tmp, dudydudy_tmp, dudzdudz_tmp ! YMG: For Dissipation
        double precision :: dvdxdvdx_tmp, dvdydvdy_tmp, dvdzdvdz_tmp
        double precision :: dwdxdwdx_tmp, dwdydwdy_tmp, dwdzdwdz_tmp
        double precision :: dudxdvdx_tmp, dudydvdy_tmp, dudzdvdz_tmp
        double precision :: udpdx_tmp, vdpdy_tmp, wdpdz_tmp, udpdy_tmp, vdpdx_tmp ! YMG: For V-P Gradient

        double precision :: UUUavg_t_tmp, UUVavg_t_tmp, UUWavg_t_tmp ! YMG: For Turbulent Diffusion
        double precision :: VVUavg_t_tmp, VVVavg_t_tmp, VVWavg_t_tmp
        double precision :: WWUavg_t_tmp, WWVavg_t_tmp, WWWavg_t_tmp
        double precision :: UVUavg_t_tmp, UVVavg_t_tmp, UVWavg_t_tmp

        ! s1msb = s1sb-1
        ! s2msb = s2sb-1
        ! s3msb = s3sb-1

        dx1  => dx1_sub
        dx2  => dx2_sub
        dx3  => dx3_sub
        dmx1 => dmx1_sub
        dmx2 => dmx2_sub
        dmx3 => dmx3_sub

        avg_time_length=avg_time_length+dt 

        do j=1, n2sub
            Umean_tmp =0.0d0
            Vmean_tmp =0.0d0
            Wmean_tmp =0.0d0
            Tmean_tmp =0.0d0
            Pmean_tmp =0.0d0
            Mumean_tmp=0.0d0
            do k= 1, n3msub
            do i =1, n1msub
                Umean_tmp= Umean_tmp + ( U(i,j,k)*dx2(j-1)+U(i,j-1,k)*dx2(j) )/(dx2(j)+dx2(j-1)) /dble(n3msub*n1msub)*dtime
                Vmean_tmp= Vmean_tmp +   V(i,j,k)/dble(n3msub*n1msub)*dtime
                Wmean_tmp= Wmean_tmp + ( W(i,j,k)*dx2(j-1)+W(i,j-1,k)*dx2(j) )/(dx2(j)+dx2(j-1)) /dble(n3msub*n1msub)*dtime
                ! For Building Array
                ! if(MOD(i,40).eq.1 .and. MOD(k,40).eq.21 ) then
                !     u1 = ( U(i,j,k-1)*dx2(j-1)+U(i,j-1,k-1)*dx2(j) )/( dx2(j)+dx2(j-1) )
                !     u2 = ( U(i,j,k  )*dx2(j-1)+U(i,j-1,k  )*dx2(j) )/( dx2(j)+dx2(j-1) )
                !     ug = ( u2*dx3(k-1)+u1*dx3(k) )/( dx3(k)+dx3(k-1) )

                !     v1 = ( V(i,j,k-1)*dx1(i-1)+V(i-1,j,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
                !     v2 = ( V(i,j,k  )*dx1(i-1)+V(i-1,j,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
                !     vg = ( v2*dx3(k-1)+v1*dx3(k) )/( dx3(k)+dx3(k-1) )
        
                !     w1 = ( W(i-1,j,k)*dx2(j-1)+W(i-1,j-1,k)*dx2(j) )/( dx2(j)+dx2(j-1) )
                !     w2 = ( W(i  ,j,k)*dx2(j-1)+W(i  ,j-1,k)*dx2(j) )/( dx2(j)+dx2(j-1) )
                !     wg = ( w2*dx1(i-1)+w1*dx1(i) )/( dx1(i)+dx1(i-1) )

                !     Umean_tmp= Umean_tmp + ug /dble(144)*dtime
                !     Vmean_tmp= Vmean_tmp + vg /dble(144)*dtime
                !     Wmean_tmp= Wmean_tmp + wg /dble(144)*dtime
                ! endif
            ! Vmean_tmp= Vmean_tmp +   V(i,j,k)/dble(s3msb*s1msb)*dtime
            ! Wmean_tmp= Wmean_tmp + ( W(i,j,k)*dx2(j-1)+W(i,j-1,k)*dx2(j) )/(dx2(j)+dx2(j-1)) /dble(s3msb*s1msb)*dtime
                Tmean_tmp  = Tmean_tmp + ( T(i,j,k)*dx2(j-1)+ T(i,j-1,k)*dx2(j) )/(dx2(j)+dx2(j-1)) /dble(n3msub*n1msub)*dtime
                Pmean_tmp  = Pmean_tmp + ( P(i,j,k)*dx2(j-1)+ P(i,j-1,k)*dx2(j) )/(dx2(j)+dx2(j-1)) /dble(n3msub*n1msub)*dtime
                Mumean_tmp = Mumean_tmp + (Mu(i,j,k)*dx2(j-1)+Mu(i,j-1,k)*dx2(j) )/(dx2(j)+dx2(j-1)) /dble(n3msub*n1msub)*dtime
            enddo
            enddo
            Umean(j)= Umean(j) + Umean_tmp
            Vmean(j)= Vmean(j) + Vmean_tmp
            Wmean(j)= Wmean(j) + Wmean_tmp
            Pmean(j)= Pmean(j) + Pmean_tmp
            Tmean(j)= Tmean(j) + Tmean_tmp
            Mumean(j)=Mumean(j) + Mumean_tmp
        enddo
    
        temp1=Umean(2)
        temp2=Umean(1)
        temp3=dx2(1) 
        temp4=Tmean(2)
        temp5=Tmean(1)

		! YMG: Non-dimensionlized Parameter Calculation
        u_tau =  SQRT( Cmu * ( temp1/avg_time_length - temp2/avg_time_length ) / ( temp3 ) )
		delta_nu =  Cmu / u_tau
		Temp_tau =  - (temp4 - temp5 ) * Ct /  u_tau

        ! !$acc parallel loop private(u1,u2,ug,v1,v2,vg,w1,w2,wg,t1,t2,t3,t4,t5,t6,tg,p1,p2,p3,p4,p5,p6,pg)
        ! do k=1,s3sb
        ! !$acc loop vector
        ! do j=1,s2sb
        ! !$acc loop vector
        ! do i=1,s1sb

        !     u1 = ( U(i,j,k-1)*dx2(j-1)+U(i,j-1,k-1)*dx2(j) )/( dx2(j)+dx2(j-1) )
        !     u2 = ( U(i,j,k  )*dx2(j-1)+U(i,j-1,k  )*dx2(j) )/( dx2(j)+dx2(j-1) )
        !     ug = ( u2*dx3(k-1)+u1*dx3(k) )/( dx3(k)+dx3(k-1) )
    
        !     v1 = ( V(i,j,k-1)*dx1(i-1)+V(i-1,j,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
        !     v2 = ( V(i,j,k  )*dx1(i-1)+V(i-1,j,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
        !     vg = ( v2*dx3(k-1)+v1*dx3(k) )/( dx3(k)+dx3(k-1) )

        !     w1 = ( W(i-1,j,k)*dx2(j-1)+W(i-1,j-1,k)*dx2(j) )/( dx2(j)+dx2(j-1) )
        !     w2 = ( W(i  ,j,k)*dx2(j-1)+W(i  ,j-1,k)*dx2(j) )/( dx2(j)+dx2(j-1) )
        !     wg = ( w2*dx1(i-1)+w1*dx1(i) )/( dx1(i)+dx1(i-1) )

        !     t1 = ( T(i  ,j  ,k  )*dx1(i-1)+T(i-1,j  ,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
        !     t2 = ( T(i  ,j-1,k  )*dx1(i-1)+T(i-1,j-1,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
        !     t3 = ( T(i  ,j  ,k-1)*dx1(i-1)+T(i-1,j  ,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
        !     t4 = ( T(i  ,j-1,k-1)*dx1(i-1)+T(i-1,j-1,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
    
        !     t5 = ( t1*dx2(j-1)+t2*dx2(j) )/( dx2(j)+dx2(j-1) )
        !     t6 = ( t3*dx2(j-1)+t4*dx2(j) )/( dx2(j)+dx2(j-1) )
        !     tg = ( t5*dx3(k-1)+t6*dx3(k) )/( dx3(k)+dx3(k-1) )
            
        !     p1 = ( P(i  ,j  ,k  )*dx1(i-1)+P(i-1,j  ,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
        !     p2 = ( P(i  ,j-1,k  )*dx1(i-1)+P(i-1,j-1,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
        !     p3 = ( P(i  ,j  ,k-1)*dx1(i-1)+P(i-1,j  ,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
        !     p4 = ( P(i  ,j-1,k-1)*dx1(i-1)+P(i-1,j-1,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
    
        !     p5 = ( p1*dx2(j-1)+p2*dx2(j) )/( dx2(j)+dx2(j-1) )
        !     p6 = ( p3*dx2(j-1)+p4*dx2(j) )/( dx2(j)+dx2(j-1) )
        !     pg = ( p5*dx3(k-1)+p6*dx3(k) )/( dx3(k)+dx3(k-1) )

        !     !$acc atomic update
        !     UUavg_t(j) =  UUavg_t(j)+ (ug-Umean(j)/avg_time_length)*(ug-Umean(j)/avg_time_length)/dble(s3sb*s1sb)*dtime
        !     !$acc atomic update
        !     VVavg_t(j) =  VVavg_t(j)+ (vg-Vmean(j)/avg_time_length)*(vg-Vmean(j)/avg_time_length)/dble(s3sb*s1sb)*dtime
        !     !$acc atomic update
        !     WWavg_t(j) =  WWavg_t(j)+ (wg-Wmean(j)/avg_time_length)*(wg-Wmean(j)/avg_time_length)/dble(s3sb*s1sb)*dtime
        !     !$acc atomic update
        !     UVavg_t(j) =  UVavg_t(j)+ (ug-Umean(j)/avg_time_length)*(vg-Vmean(j)/avg_time_length)/dble(s3sb*s1sb)*dtime
        !     !$acc atomic update
        !     TTavg_t(j) =  TTavg_t(j)+ (tg-Tmean(j)/avg_time_length)*(tg-Tmean(j)/avg_time_length)/dble(s3sb*s1sb)*dtime
        !     !$acc atomic update
        !     UTavg_t(j) =  UTavg_t(j)+ (ug-Umean(j)/avg_time_length)*(tg-Tmean(j)/avg_time_length)/dble(s3sb*s1sb)*dtime
        !     !$acc atomic update
        !     UWavg_t(j) =  UWavg_t(j)+ (ug-Umean(j)/avg_time_length)*(wg-Wmean(j)/avg_time_length)/dble(s3sb*s1sb)*dtime
        !     !$acc atomic update
        !     VWavg_t(j) =  VWavg_t(j)+ (vg-Vmean(j)/avg_time_length)*(wg-Wmean(j)/avg_time_length)/dble(s3sb*s1sb)*dtime
        !     !$acc atomic update
        !     Usq(j)=Usq(j) + ug*ug/dble(s3sb*s1sb)*dtime
        !     !$acc atomic update
        !     Vsq(j)=Vsq(j) + vg*vg/dble(s3sb*s1sb)*dtime
        !     !$acc atomic update
        !     Wsq(j)=Wsq(j) + wg*wg/dble(s3sb*s1sb)*dtime
        !     !$acc atomic update
        !     Tsq(j)=Tsq(j) + tg*tg/dble(s3sb*s1sb)*dtime
        !     !$acc atomic update
        !     UV(j) =UV(j)  + ug*vg/dble(s3sb*s1sb)*dtime
        !     !$acc atomic update
        !     UT(j) =UT(j)  + ug*tg/dble(s3sb*s1sb)*dtime
        ! enddo
        ! enddo
        ! enddo

        do j=1,n2sub
            UUavg_t_tmp=0.0d0
            VVavg_t_tmp=0.0d0
            WWavg_t_tmp=0.0d0
            TTavg_t_tmp=0.0d0
            UVavg_t_tmp=0.0d0
            UTavg_t_tmp=0.0d0
            UWavg_t_tmp=0.0d0
            VWavg_t_tmp=0.0d0
            Usq_tmp=0.0d0
            Vsq_tmp=0.0d0
            Wsq_tmp=0.0d0
            Tsq_tmp=0.0d0
            UV_tmp=0.0d0
            UT_tmp=0.0d0
            do k=1,n3msub
            do i=1,n1msub
                u1 = ( U(i,j,k-1)*dx2(j-1)+U(i,j-1,k-1)*dx2(j) )/( dx2(j)+dx2(j-1) )
                u2 = ( U(i,j,k  )*dx2(j-1)+U(i,j-1,k  )*dx2(j) )/( dx2(j)+dx2(j-1) )
                ug = ( u2*dx3(k-1)+u1*dx3(k) )/( dx3(k)+dx3(k-1) )
        
                v1 = ( V(i,j,k-1)*dx1(i-1)+V(i-1,j,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
                v2 = ( V(i,j,k  )*dx1(i-1)+V(i-1,j,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
                vg = ( v2*dx3(k-1)+v1*dx3(k) )/( dx3(k)+dx3(k-1) )
    
                w1 = ( W(i-1,j,k)*dx2(j-1)+W(i-1,j-1,k)*dx2(j) )/( dx2(j)+dx2(j-1) )
                w2 = ( W(i  ,j,k)*dx2(j-1)+W(i  ,j-1,k)*dx2(j) )/( dx2(j)+dx2(j-1) )
                wg = ( w2*dx1(i-1)+w1*dx1(i) )/( dx1(i)+dx1(i-1) )
    
                t1 = ( T(i  ,j  ,k  )*dx1(i-1)+T(i-1,j  ,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
                t2 = ( T(i  ,j-1,k  )*dx1(i-1)+T(i-1,j-1,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
                t3 = ( T(i  ,j  ,k-1)*dx1(i-1)+T(i-1,j  ,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
                t4 = ( T(i  ,j-1,k-1)*dx1(i-1)+T(i-1,j-1,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
        
                t5 = ( t1*dx2(j-1)+t2*dx2(j) )/( dx2(j)+dx2(j-1) )
                t6 = ( t3*dx2(j-1)+t4*dx2(j) )/( dx2(j)+dx2(j-1) )
                tg = ( t5*dx3(k-1)+t6*dx3(k) )/( dx3(k)+dx3(k-1) )
                
                p1 = ( P(i  ,j  ,k  )*dx1(i-1)+P(i-1,j  ,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
                p2 = ( P(i  ,j-1,k  )*dx1(i-1)+P(i-1,j-1,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
                p3 = ( P(i  ,j  ,k-1)*dx1(i-1)+P(i-1,j  ,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
                p4 = ( P(i  ,j-1,k-1)*dx1(i-1)+P(i-1,j-1,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
        
                p5 = ( p1*dx2(j-1)+p2*dx2(j) )/( dx2(j)+dx2(j-1) )
                p6 = ( p3*dx2(j-1)+p4*dx2(j) )/( dx2(j)+dx2(j-1) )
                pg = ( p5*dx3(k-1)+p6*dx3(k) )/( dx3(k)+dx3(k-1) )
    
            !    UUavg_t_tmp =  UUavg_t_tmp+ (ug-Umean(j)/avg_time_length)*(ug-Umean(j)/avg_time_length)/dble(s3msb*s1msb)*dtime
            !    VVavg_t_tmp =  VVavg_t_tmp+ (vg-Vmean(j)/avg_time_length)*(vg-Vmean(j)/avg_time_length)/dble(s3msb*s1msb)*dtime
            !    WWavg_t_tmp =  WWavg_t_tmp+ (wg-Wmean(j)/avg_time_length)*(wg-Wmean(j)/avg_time_length)/dble(s3msb*s1msb)*dtime
            !    UVavg_t_tmp =  UVavg_t_tmp+ (ug-Umean(j)/avg_time_length)*(vg-Vmean(j)/avg_time_length)/dble(s3msb*s1msb)*dtime
                if(MOD(i,40).eq.1 .and. MOD(k,40).eq.21 ) then ! In boundary layer, UW
                    UUavg_t_tmp =  UUavg_t_tmp+ (ug-Umean(j)/avg_time_length)*(ug-Umean(j)/avg_time_length)/dble(144)*dtime
                    VVavg_t_tmp =  VVavg_t_tmp+ (vg-Vmean(j)/avg_time_length)*(vg-Vmean(j)/avg_time_length)/dble(144)*dtime
                    WWavg_t_tmp =  WWavg_t_tmp+ (wg-Wmean(j)/avg_time_length)*(wg-Wmean(j)/avg_time_length)/dble(144)*dtime
                    UVavg_t_tmp =  UVavg_t_tmp+ (ug-Umean(j)/avg_time_length)*(vg-Vmean(j)/avg_time_length)/dble(144)*dtime
                endif
                TTavg_t_tmp =  TTavg_t_tmp+ (tg-Tmean(j)/avg_time_length)*(tg-Tmean(j)/avg_time_length)/dble(n3msub*n1msub)*dtime
                UTavg_t_tmp =  UTavg_t_tmp+ (ug-Umean(j)/avg_time_length)*(tg-Tmean(j)/avg_time_length)/dble(n3msub*n1msub)*dtime
                UWavg_t_tmp =  UWavg_t_tmp+ (ug-Umean(j)/avg_time_length)*(wg-Wmean(j)/avg_time_length)/dble(n3msub*n1msub)*dtime
                VWavg_t_tmp =  VWavg_t_tmp+ (vg-Vmean(j)/avg_time_length)*(wg-Wmean(j)/avg_time_length)/dble(n3msub*n1msub)*dtime

               ! Usq_tmp= Usq_tmp + ug*ug/dble(s3sb*s1sb)*dtime
               ! Vsq_tmp= Vsq_tmp + vg*vg/dble(s3sb*s1sb)*dtime
               ! Wsq_tmp= Wsq_tmp + wg*wg/dble(s3sb*s1sb)*dtime
                if(MOD(i,40).eq.1 .and. MOD(k,40).eq.21 ) then !In boundary layer, UV
                    Usq_tmp= Usq_tmp + ug*ug/dble(144)*dtime
                    Vsq_tmp= Vsq_tmp + vg*vg/dble(144)*dtime
                    Wsq_tmp= Wsq_tmp + wg*wg/dble(144)*dtime
               endif
                Tsq_tmp= Tsq_tmp + tg*tg/dble(n3msub*n1msub)*dtime
                UV_tmp = UV_tmp  + ug*vg/dble(n3msub*n1msub)*dtime
                UT_tmp = UT_tmp  + ug*tg/dble(n3msub*n1msub)*dtime
            enddo
            enddo
            UUavg_t(j) =  UUavg_t(j)+ UUavg_t_tmp 
            VVavg_t(j) =  VVavg_t(j)+ VVavg_t_tmp 
            WWavg_t(j) =  WWavg_t(j)+ WWavg_t_tmp 
            UVavg_t(j) =  UVavg_t(j)+ UVavg_t_tmp 
            TTavg_t(j) =  TTavg_t(j)+ TTavg_t_tmp 
            UTavg_t(j) =  UTavg_t(j)+ UTavg_t_tmp 
            UWavg_t(j) =  UWavg_t(j)+ UWavg_t_tmp 
            VWavg_t(j) =  VWavg_t(j)+ VWavg_t_tmp 
            Usq(j)=Usq(j) + Usq_tmp
            Vsq(j)=Vsq(j) + Vsq_tmp
            Wsq(j)=Wsq(j) + Wsq_tmp
            Tsq(j)=Tsq(j) + Tsq_tmp
            UV(j) =UV(j) + UV_tmp
            UT(j) =UT(j) + UT_tmp 
        enddo


        ! !$acc parallel loop private(u1,u2,ug,ugp,ugm,v1,v2,vg,vgp,vgm,w1,w2,wg,wgp,wgm,t1,t2,t3,t4,t5,t6,tg,p1,p2,p3,p4,p5,p6,pg,pgxp,pgxm,pgyp,pgym,pgzp,pgzm)
        ! do k=2,s3msb
        ! !$acc loop vector
        ! do j=2,s2msb
        ! !$acc loop vector
        ! do i=2,s1msb
        !     u1 = ( U(i,j,k-1)*dx2(j-1)+U(i,j-1,k-1)*dx2(j) )/( dx2(j)+dx2(j-1) )
        !     u2 = ( U(i,j,k  )*dx2(j-1)+U(i,j-1,k  )*dx2(j) )/( dx2(j)+dx2(j-1) )
        !     ug = ( u2*dx3(k-1)+u1*dx3(k) )/( dx3(k)+dx3(k-1) )

        !     u1 = ( U(i,j+1,k-1)*dx2(j)+U(i,j,k-1)*dx2(j+1) )/( dx2(j+1)+dx2(j) )
        !     u2 = ( U(i,j+1,k  )*dx2(j)+U(i,j,k  )*dx2(j+1) )/( dx2(j+1)+dx2(j) )
        !    ugp = ( u2*dx3(k-1)+u1*dx3(k) )/( dx3(k)+dx3(k-1) )
			
        !     u1 = ( U(i,j-1,k-1)*dx2(j-2)+U(i,j-2,k-1)*dx2(j-1) )/( dx2(j-1)+dx2(j-2) )
        !     u2 = ( U(i,j-1,k  )*dx2(j-2)+U(i,j-2,k  )*dx2(j-1) )/( dx2(j-1)+dx2(j-2) )
        !    ugm = ( u2*dx3(k-1)+u1*dx3(k) )/( dx3(k)+dx3(k-1) )			

        !     v1 = ( V(i,j,k-1)*dx1(i-1)+V(i-1,j,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
        !     v2 = ( V(i,j,k  )*dx1(i-1)+V(i-1,j,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
        !     vg = ( v2*dx3(k-1)+v1*dx3(k) )/( dx3(k)+dx3(k-1) )
			
        !     v1 = ( V(i,j+1,k-1)*dx1(i-1)+V(i-1,j+1,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
        !     v2 = ( V(i,j+1,k  )*dx1(i-1)+V(i-1,j+1,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
        !    vgp = ( v2*dx3(k-1)+v1*dx3(k) )/( dx3(k)+dx3(k-1) )
			
		! 	v1 = ( V(i,j-1,k-1)*dx1(i-1)+V(i-1,j-1,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
        !     v2 = ( V(i,j-1,k  )*dx1(i-1)+V(i-1,j-1,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
        !    vgm = ( v2*dx3(k-1)+v1*dx3(k) )/( dx3(k)+dx3(k-1) )
			
        !     w1 = ( W(i-1,j,k)*dx2(j-1)+W(i-1,j-1,k)*dx2(j) )/( dx2(j)+dx2(j-1) )
        !     w2 = ( W(i  ,j,k)*dx2(j-1)+W(i  ,j-1,k)*dx2(j) )/( dx2(j)+dx2(j-1) )
        !     wg = ( w2*dx1(i-1)+w1*dx1(i) )/( dx1(i)+dx1(i-1) )
			
		! 	w1 = ( W(i-1,j+1,k)*dx2(j)+W(i-1,j,k)*dx2(j+1) )/( dx2(j+1)+dx2(j) )
        !     w2 = ( W(i  ,j+1,k)*dx2(j)+W(i  ,j,k)*dx2(j+1) )/( dx2(j+1)+dx2(j) )
        !    wgp = ( w2*dx1(i-1)+w1*dx1(i) )/( dx1(i)+dx1(i-1) )
			
		! 	w1 = ( W(i-1,j-1,k)*dx2(j-2)+W(i-1,j-2,k)*dx2(j-1) )/( dx2(j-1)+dx2(j-2) )
        !     w2 = ( W(i  ,j-1,k)*dx2(j-2)+W(i  ,j-2,k)*dx2(j-1) )/( dx2(j-1)+dx2(j-2) )
        !    wgm = ( w2*dx1(i-1)+w1*dx1(i) )/( dx1(i)+dx1(i-1) )

        !     t1 = ( T(i  ,j  ,k  )*dx1(i-1)+T(i-1,j  ,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
        !     t2 = ( T(i  ,j-1,k  )*dx1(i-1)+T(i-1,j-1,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
        !     t3 = ( T(i  ,j  ,k-1)*dx1(i-1)+T(i-1,j  ,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
        !     t4 = ( T(i  ,j-1,k-1)*dx1(i-1)+T(i-1,j-1,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )

        !     t5 = ( t1*dx2(j-1)+t2*dx2(j) )/( dx2(j)+dx2(j-1) )
        !     t6 = ( t3*dx2(j-1)+t4*dx2(j) )/( dx2(j)+dx2(j-1) )
        !     tg = ( t5*dx3(k-1)+t6*dx3(k) )/( dx3(k)+dx3(k-1) )
			
        !     p1 = ( P(i  ,j  ,k  )*dx1(i-1)+P(i-1,j  ,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
        !     p2 = ( P(i  ,j-1,k  )*dx1(i-1)+P(i-1,j-1,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
        !     p3 = ( P(i  ,j  ,k-1)*dx1(i-1)+P(i-1,j  ,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
        !     p4 = ( P(i  ,j-1,k-1)*dx1(i-1)+P(i-1,j-1,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )

        !     p5 = ( p1*dx2(j-1)+p2*dx2(j) )/( dx2(j)+dx2(j-1) )
        !     p6 = ( p3*dx2(j-1)+p4*dx2(j) )/( dx2(j)+dx2(j-1) )
        !     pg = ( p5*dx3(k-1)+p6*dx3(k) )/( dx3(k)+dx3(k-1) )
			
		! 	! Y-direction 
		! 	p1 = ( P(i  ,j+1,k  )*dx1(i-1)+P(i-1,j+1,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
        !     p2 = ( P(i  ,j  ,k  )*dx1(i-1)+P(i-1,j  ,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
        !     p3 = ( P(i  ,j+1,k-1)*dx1(i-1)+P(i-1,j+1,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
        !     p4 = ( P(i  ,j  ,k-1)*dx1(i-1)+P(i-1,j  ,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )

        !     p5 = ( p1*dx2(j-1)+p2*dx2(j+1) )/( dx2(j)+dx2(j) )
        !     p6 = ( p3*dx2(j-1)+p4*dx2(j+1) )/( dx2(j)+dx2(j) )
        !   pgyp = ( p5*dx3(k-1)+p6*dx3(k) )/( dx3(k)+dx3(k-1) )
			
		! 	p1 = ( P(i  ,j-1,k  )*dx1(i-1)+P(i-1,j-1,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
        !     p2 = ( P(i  ,j-2,k  )*dx1(i-1)+P(i-1,j-2,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
        !     p3 = ( P(i  ,j-1,k-1)*dx1(i-1)+P(i-1,j-1,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
        !     p4 = ( P(i  ,j-2,k-1)*dx1(i-1)+P(i-1,j-2,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )

        !     p5 = ( p1*dx2(j-2)+p2*dx2(j-1) )/( dx2(j-1)+dx2(j-2) )
        !     p6 = ( p3*dx2(j-2)+p4*dx2(j-1) )/( dx2(j-1)+dx2(j-2) )
        !   pgym = ( p5*dx3(k-1)+p6*dx3(k) )/( dx3(k)+dx3(k-1) )
		   
		!     ! Z-direction
        !     p1 = ( P(i  ,j  ,k+1)*dx1(i-1)+P(i-1,j  ,k+1)*dx1(i) )/( dx1(i)+dx1(i-1) )
        !     p2 = ( P(i  ,j-1,k+1)*dx1(i-1)+P(i-1,j-1,k+1)*dx1(i) )/( dx1(i)+dx1(i-1) )
        !     p3 = ( P(i  ,j  ,k  )*dx1(i-1)+P(i-1,j  ,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
        !     p4 = ( P(i  ,j-1,k  )*dx1(i-1)+P(i-1,j-1,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )

        !     p5 = ( p1*dx2(j-1)+p2*dx2(j  ) )/( dx2(j  )+dx2(j-1) )
        !     p6 = ( p3*dx2(j-1)+p4*dx2(j  ) )/( dx2(j  )+dx2(j-1) )
        !   pgzp = ( p5*dx3(k  )+p6*dx3(k+1) )/( dx3(k+1)+dx3(k) )
			
        !     p1 = ( P(i  ,j  ,k-1)*dx1(i-1)+P(i-1,j  ,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
        !     p2 = ( P(i  ,j-1,k-1)*dx1(i-1)+P(i-1,j-1,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
        !     p3 = ( P(i  ,j  ,k-2)*dx1(i-1)+P(i-1,j  ,k-2)*dx1(i) )/( dx1(i)+dx1(i-1) )
        !     p4 = ( P(i  ,j-1,k-2)*dx1(i-1)+P(i-1,j-1,k-2)*dx1(i) )/( dx1(i)+dx1(i-1) )

        !     p5 = ( p1*dx2(j-1)+p2*dx2(j  ) )/( dx2(j  )+dx2(j-1) )
        !     p6 = ( p3*dx2(j-1)+p4*dx2(j  ) )/( dx2(j  )+dx2(j-1) )
        !   pgzm = ( p5*dx3(k-2)+p6*dx3(k-1) )/( dx3(k-1)+dx3(k-2) )
		   
		!    	! X-direction
        !     p1 = ( P(i+1,j  ,k  )*dx1(i)+P(i,j  ,k  )*dx1(i+1) )/( dx1(i+1)+dx1(i) )
        !     p2 = ( P(i+1,j-1,k  )*dx1(i)+P(i,j-1,k  )*dx1(i+1) )/( dx1(i+1)+dx1(i) )
        !     p3 = ( P(i+1,j  ,k-1)*dx1(i)+P(i,j  ,k-1)*dx1(i+1) )/( dx1(i+1)+dx1(i) )
        !     p4 = ( P(i+1,j-1,k-1)*dx1(i)+P(i,j-1,k-1)*dx1(i+1) )/( dx1(i+1)+dx1(i) )

        !     p5 = ( p1*dx2(j-1)+p2*dx2(j) )/( dx2(j)+dx2(j-1) )
        !     p6 = ( p3*dx2(j-1)+p4*dx2(j) )/( dx2(j)+dx2(j-1) )
        !   pgxp = ( p5*dx3(k-1)+p6*dx3(k) )/( dx3(k)+dx3(k-1) )
			
        !     p1 = ( P(i-1,j  ,k  )*dx1(i-1)+P(i-2,j  ,k  )*dx1(i) )/( dx1(i-2)+dx1(i-1) )
        !     p2 = ( P(i-1,j-1,k  )*dx1(i-1)+P(i-2,j-1,k  )*dx1(i) )/( dx1(i-2)+dx1(i-1) )
        !     p3 = ( P(i-1,j  ,k-1)*dx1(i-1)+P(i-2,j  ,k-1)*dx1(i) )/( dx1(i-2)+dx1(i-1) )
        !     p4 = ( P(i-1,j-1,k-1)*dx1(i-1)+P(i-2,j-1,k-1)*dx1(i) )/( dx1(i-2)+dx1(i-1) )

        !     p5 = ( p1*dx2(j-1)+p2*dx2(j) )/( dx2(j)+dx2(j-1) )
        !     p6 = ( p3*dx2(j-1)+p4*dx2(j) )/( dx2(j)+dx2(j-1) )
        !   pgxm = ( p5*dx3(k-1)+p6*dx3(k) )/( dx3(k)+dx3(k-1) )

		! 	! YMG: For Turbulent Diffusion
		! 	! UUUavg_t(j) = 0.0d0
		! 	! UUVavg_t(j) = UUVavg_t(j)+ (ug-Umean(j)/avg_time_length)*(ug-Umean(j)/avg_time_length)*(vg-Vmean(j)/avg_time_length)/dble((s3msb-1)*(s1msb-1))*dtime
		! 	! UUWavg_t(j) = 0.0d0
		! 	! VVUavg_t(j) = 0.0d0
		! 	! VVVavg_t(j) = VVVavg_t(j)+ (vg-Vmean(j)/avg_time_length)*(vg-Vmean(j)/avg_time_length)*(vg-Vmean(j)/avg_time_length)/dble((s3msb-1)*(s1msb-1))*dtime
		! 	! VVWavg_t(j) = 0.0d0
		! 	! WWUavg_t(j) = 0.0d0
		! 	! WWVavg_t(j) = WWVavg_t(j)+ (wg-Wmean(j)/avg_time_length)*(wg-Wmean(j)/avg_time_length)*(vg-Vmean(j)/avg_time_length)/dble((s3msb-1)*(s1msb-1))*dtime
		! 	! WWWavg_t(j) = 0.0d0
		! 	! UVUavg_t(j) = 0.0d0
		! 	! UVVavg_t(j) = UVVavg_t(j)+ (ug-Umean(j)/avg_time_length)*(vg-Vmean(j)/avg_time_length)*(vg-Vmean(j)/avg_time_length)/dble((s3msb-1)*(s1msb-1))*dtime
		! 	! UVWavg_t(j) = 0.0d0

		! 	! YMG: For Dissipation
		! 	dudxdudx(j) = 0.0d0
        !     !$acc atomic update
		! 	dudydudy(j) = dudydudy(j)+ ( ((ugp-Umean(j+1)/avg_time_length) - (ugm-Umean(j-1)/avg_time_length) ) / ( dx2(j) + dx2(j-1) ))**dble(2.0) /dble((s3msb-1)*(s1msb-1))*dtime
		! 	dudzdudz(j) = 0.0d0
		! 	dvdxdvdx(j) = 0.0d0
        !     !$acc atomic update
		! 	dvdydvdy(j) = dvdydvdy(j)+ ( ((vgp-Vmean(j+1)/avg_time_length) - (vgm-Vmean(j-1)/avg_time_length) ) / ( dx2(j) + dx2(j-1) ))**dble(2.0) /dble((s3msb-1)*(s1msb-1))*dtime
		! 	dvdzdvdz(j) = 0.0d0
		! 	dwdxdwdx(j) = 0.0d0
        !     !$acc atomic update
		! 	dwdydwdy(j) = dwdydwdy(j)+ ( ((wgp-Wmean(j+1)/avg_time_length) - (wgm-Wmean(j-1)/avg_time_length) ) / ( dx2(j) + dx2(j-1) ))**dble(2.0) /dble((s3msb-1)*(s1msb-1))*dtime
		! 	dwdzdwdz(j) = 0.0d0
		! 	dudxdvdx(j) = 0.0d0
        !     !$acc atomic update
		! 	dudydvdy(j) = dudydvdy(j)+   (((ugp-Umean(j+1)/avg_time_length) - (ugm-Umean(j-1)/avg_time_length) ) / ( dx2(j) + dx2(j-1) ))  &
		! 	                              & *(((vgp-Vmean(j+1)/avg_time_length) - (vgm-Vmean(j-1)/avg_time_length) ) / ( dx2(j) + dx2(j-1) ))/dble((s3msb-1)*(s1msb-1))*dtime
		! 	dudzdvdz(j) = 0.0d0
			
		! 	! YMG: For V-P Gradient
        !     !$acc atomic update
		! 	udpdx(j)=udpdx(j) + (ug-Umean(j)/avg_time_length)* ( (pgxp-Pmean(j  )/avg_time_length) - (pgxm-Pmean(j  )/avg_time_length) ) / ( dx1(i) + dx1(i-1) ) /dble((s3msb-1)*(s1msb-1))*dtime
        !     !$acc atomic update
        !     vdpdy(j)=vdpdy(j) + (vg-Vmean(j)/avg_time_length)* ( (pgyp-Pmean(j+1)/avg_time_length) - (pgym-Pmean(j-1)/avg_time_length) ) / ( dx2(j) + dx2(j-1) ) /dble((s3msb-1)*(s1msb-1))*dtime
        !     !$acc atomic update
        !     wdpdz(j)=wdpdz(j) + (wg-Wmean(j)/avg_time_length)* ( (pgzp-Pmean(j  )/avg_time_length) - (pgzm-Pmean(j  )/avg_time_length) ) / ( dx3(k) + dx3(k-1) ) /dble((s3msb-1)*(s1msb-1))*dtime
        !     !$acc atomic update
        !     udpdy(j)=udpdy(j) + (ug-Umean(j)/avg_time_length)* ( (pgyp-Pmean(j+1)/avg_time_length) - (pgym-Pmean(j-1)/avg_time_length) ) / ( dx2(j) + dx2(j-1) ) /dble((s3msb-1)*(s1msb-1))*dtime
        !     !$acc atomic update
        !     vdpdx(j)=vdpdx(j) + (vg-Vmean(j)/avg_time_length)* ( (pgxp-Pmean(j  )/avg_time_length) - (pgxm-Pmean(j  )/avg_time_length) ) / ( dx1(i) + dx1(i-1) ) /dble((s3msb-1)*(s1msb-1))*dtime
        ! enddo
        ! enddo
        ! enddo

        do j=2,n2msub
            dudxdudx_tmp=0.0d0; dudydudy_tmp=0.0d0;  dudzdudz_tmp=0.0d0
            dvdxdvdx_tmp=0.0d0; dvdydvdy_tmp=0.0d0;  dvdzdvdz_tmp=0.0d0
            dwdxdwdx_tmp=0.0d0; dwdydwdy_tmp=0.0d0;  dwdzdwdz_tmp=0.0d0
            dudxdvdx_tmp=0.0d0; dudydvdy_tmp=0.0d0;  dudzdvdz_tmp=0.0d0
            udpdx_tmp=0.0d0; vdpdy_tmp=0.0d0;  wdpdz_tmp=0.0d0;  udpdy_tmp=0.0d0;  vdpdx_tmp=0.0d0
            do k=2,n3msub
            do i=2,n1msub
                u1 = ( U(i,j,k-1)*dx2(j-1)+U(i,j-1,k-1)*dx2(j) )/( dx2(j)+dx2(j-1) )
                u2 = ( U(i,j,k  )*dx2(j-1)+U(i,j-1,k  )*dx2(j) )/( dx2(j)+dx2(j-1) )
                ug = ( u2*dx3(k-1)+u1*dx3(k) )/( dx3(k)+dx3(k-1) )
    
                u1 = ( U(i,j+1,k-1)*dx2(j)+U(i,j,k-1)*dx2(j+1) )/( dx2(j+1)+dx2(j) )
                u2 = ( U(i,j+1,k  )*dx2(j)+U(i,j,k  )*dx2(j+1) )/( dx2(j+1)+dx2(j) )
               ugp = ( u2*dx3(k-1)+u1*dx3(k) )/( dx3(k)+dx3(k-1) )
                
                u1 = ( U(i,j-1,k-1)*dx2(j-2)+U(i,j-2,k-1)*dx2(j-1) )/( dx2(j-1)+dx2(j-2) )
                u2 = ( U(i,j-1,k  )*dx2(j-2)+U(i,j-2,k  )*dx2(j-1) )/( dx2(j-1)+dx2(j-2) )
               ugm = ( u2*dx3(k-1)+u1*dx3(k) )/( dx3(k)+dx3(k-1) )			
    
                v1 = ( V(i,j,k-1)*dx1(i-1)+V(i-1,j,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
                v2 = ( V(i,j,k  )*dx1(i-1)+V(i-1,j,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
                vg = ( v2*dx3(k-1)+v1*dx3(k) )/( dx3(k)+dx3(k-1) )
                
                v1 = ( V(i,j+1,k-1)*dx1(i-1)+V(i-1,j+1,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
                v2 = ( V(i,j+1,k  )*dx1(i-1)+V(i-1,j+1,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
               vgp = ( v2*dx3(k-1)+v1*dx3(k) )/( dx3(k)+dx3(k-1) )
                
                v1 = ( V(i,j-1,k-1)*dx1(i-1)+V(i-1,j-1,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
                v2 = ( V(i,j-1,k  )*dx1(i-1)+V(i-1,j-1,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
               vgm = ( v2*dx3(k-1)+v1*dx3(k) )/( dx3(k)+dx3(k-1) )
                
                w1 = ( W(i-1,j,k)*dx2(j-1)+W(i-1,j-1,k)*dx2(j) )/( dx2(j)+dx2(j-1) )
                w2 = ( W(i  ,j,k)*dx2(j-1)+W(i  ,j-1,k)*dx2(j) )/( dx2(j)+dx2(j-1) )
                wg = ( w2*dx1(i-1)+w1*dx1(i) )/( dx1(i)+dx1(i-1) )
                
                w1 = ( W(i-1,j+1,k)*dx2(j)+W(i-1,j,k)*dx2(j+1) )/( dx2(j+1)+dx2(j) )
                w2 = ( W(i  ,j+1,k)*dx2(j)+W(i  ,j,k)*dx2(j+1) )/( dx2(j+1)+dx2(j) )
               wgp = ( w2*dx1(i-1)+w1*dx1(i) )/( dx1(i)+dx1(i-1) )
                
                w1 = ( W(i-1,j-1,k)*dx2(j-2)+W(i-1,j-2,k)*dx2(j-1) )/( dx2(j-1)+dx2(j-2) )
                w2 = ( W(i  ,j-1,k)*dx2(j-2)+W(i  ,j-2,k)*dx2(j-1) )/( dx2(j-1)+dx2(j-2) )
               wgm = ( w2*dx1(i-1)+w1*dx1(i) )/( dx1(i)+dx1(i-1) )
    
                t1 = ( T(i  ,j  ,k  )*dx1(i-1)+T(i-1,j  ,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
                t2 = ( T(i  ,j-1,k  )*dx1(i-1)+T(i-1,j-1,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
                t3 = ( T(i  ,j  ,k-1)*dx1(i-1)+T(i-1,j  ,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
                t4 = ( T(i  ,j-1,k-1)*dx1(i-1)+T(i-1,j-1,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
    
                t5 = ( t1*dx2(j-1)+t2*dx2(j) )/( dx2(j)+dx2(j-1) )
                t6 = ( t3*dx2(j-1)+t4*dx2(j) )/( dx2(j)+dx2(j-1) )
                tg = ( t5*dx3(k-1)+t6*dx3(k) )/( dx3(k)+dx3(k-1) )
                
                p1 = ( P(i  ,j  ,k  )*dx1(i-1)+P(i-1,j  ,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
                p2 = ( P(i  ,j-1,k  )*dx1(i-1)+P(i-1,j-1,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
                p3 = ( P(i  ,j  ,k-1)*dx1(i-1)+P(i-1,j  ,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
                p4 = ( P(i  ,j-1,k-1)*dx1(i-1)+P(i-1,j-1,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
    
                p5 = ( p1*dx2(j-1)+p2*dx2(j) )/( dx2(j)+dx2(j-1) )
                p6 = ( p3*dx2(j-1)+p4*dx2(j) )/( dx2(j)+dx2(j-1) )
                pg = ( p5*dx3(k-1)+p6*dx3(k) )/( dx3(k)+dx3(k-1) )
                
                ! Y-direction 
                p1 = ( P(i  ,j+1,k  )*dx1(i-1)+P(i-1,j+1,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
                p2 = ( P(i  ,j  ,k  )*dx1(i-1)+P(i-1,j  ,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
                p3 = ( P(i  ,j+1,k-1)*dx1(i-1)+P(i-1,j+1,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
                p4 = ( P(i  ,j  ,k-1)*dx1(i-1)+P(i-1,j  ,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
    
                p5 = ( p1*dx2(j-1)+p2*dx2(j+1) )/( dx2(j)+dx2(j) )
                p6 = ( p3*dx2(j-1)+p4*dx2(j+1) )/( dx2(j)+dx2(j) )
              pgyp = ( p5*dx3(k-1)+p6*dx3(k) )/( dx3(k)+dx3(k-1) )
                
                p1 = ( P(i  ,j-1,k  )*dx1(i-1)+P(i-1,j-1,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
                p2 = ( P(i  ,j-2,k  )*dx1(i-1)+P(i-1,j-2,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
                p3 = ( P(i  ,j-1,k-1)*dx1(i-1)+P(i-1,j-1,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
                p4 = ( P(i  ,j-2,k-1)*dx1(i-1)+P(i-1,j-2,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
    
                p5 = ( p1*dx2(j-2)+p2*dx2(j-1) )/( dx2(j-1)+dx2(j-2) )
                p6 = ( p3*dx2(j-2)+p4*dx2(j-1) )/( dx2(j-1)+dx2(j-2) )
              pgym = ( p5*dx3(k-1)+p6*dx3(k) )/( dx3(k)+dx3(k-1) )
               
                ! Z-direction
                p1 = ( P(i  ,j  ,k+1)*dx1(i-1)+P(i-1,j  ,k+1)*dx1(i) )/( dx1(i)+dx1(i-1) )
                p2 = ( P(i  ,j-1,k+1)*dx1(i-1)+P(i-1,j-1,k+1)*dx1(i) )/( dx1(i)+dx1(i-1) )
                p3 = ( P(i  ,j  ,k  )*dx1(i-1)+P(i-1,j  ,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
                p4 = ( P(i  ,j-1,k  )*dx1(i-1)+P(i-1,j-1,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
    
                p5 = ( p1*dx2(j-1)+p2*dx2(j  ) )/( dx2(j  )+dx2(j-1) )
                p6 = ( p3*dx2(j-1)+p4*dx2(j  ) )/( dx2(j  )+dx2(j-1) )
              pgzp = ( p5*dx3(k  )+p6*dx3(k+1) )/( dx3(k+1)+dx3(k) )
                
                p1 = ( P(i  ,j  ,k-1)*dx1(i-1)+P(i-1,j  ,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
                p2 = ( P(i  ,j-1,k-1)*dx1(i-1)+P(i-1,j-1,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
                p3 = ( P(i  ,j  ,k-2)*dx1(i-1)+P(i-1,j  ,k-2)*dx1(i) )/( dx1(i)+dx1(i-1) )
                p4 = ( P(i  ,j-1,k-2)*dx1(i-1)+P(i-1,j-1,k-2)*dx1(i) )/( dx1(i)+dx1(i-1) )
    
                p5 = ( p1*dx2(j-1)+p2*dx2(j  ) )/( dx2(j  )+dx2(j-1) )
                p6 = ( p3*dx2(j-1)+p4*dx2(j  ) )/( dx2(j  )+dx2(j-1) )
              pgzm = ( p5*dx3(k-2)+p6*dx3(k-1) )/( dx3(k-1)+dx3(k-2) )
               
                   ! X-direction
                p1 = ( P(i+1,j  ,k  )*dx1(i)+P(i,j  ,k  )*dx1(i+1) )/( dx1(i+1)+dx1(i) )
                p2 = ( P(i+1,j-1,k  )*dx1(i)+P(i,j-1,k  )*dx1(i+1) )/( dx1(i+1)+dx1(i) )
                p3 = ( P(i+1,j  ,k-1)*dx1(i)+P(i,j  ,k-1)*dx1(i+1) )/( dx1(i+1)+dx1(i) )
                p4 = ( P(i+1,j-1,k-1)*dx1(i)+P(i,j-1,k-1)*dx1(i+1) )/( dx1(i+1)+dx1(i) )
    
                p5 = ( p1*dx2(j-1)+p2*dx2(j) )/( dx2(j)+dx2(j-1) )
                p6 = ( p3*dx2(j-1)+p4*dx2(j) )/( dx2(j)+dx2(j-1) )
              pgxp = ( p5*dx3(k-1)+p6*dx3(k) )/( dx3(k)+dx3(k-1) )
                
                p1 = ( P(i-1,j  ,k  )*dx1(i-1)+P(i-2,j  ,k  )*dx1(i) )/( dx1(i-2)+dx1(i-1) )
                p2 = ( P(i-1,j-1,k  )*dx1(i-1)+P(i-2,j-1,k  )*dx1(i) )/( dx1(i-2)+dx1(i-1) )
                p3 = ( P(i-1,j  ,k-1)*dx1(i-1)+P(i-2,j  ,k-1)*dx1(i) )/( dx1(i-2)+dx1(i-1) )
                p4 = ( P(i-1,j-1,k-1)*dx1(i-1)+P(i-2,j-1,k-1)*dx1(i) )/( dx1(i-2)+dx1(i-1) )
    
                p5 = ( p1*dx2(j-1)+p2*dx2(j) )/( dx2(j)+dx2(j-1) )
                p6 = ( p3*dx2(j-1)+p4*dx2(j) )/( dx2(j)+dx2(j-1) )
              pgxm = ( p5*dx3(k-1)+p6*dx3(k) )/( dx3(k)+dx3(k-1) )
    
                ! YMG: For Dissipation
                dudxdudx_tmp = 0.0d0
                dudydudy_tmp = dudydudy_tmp+ ( ((ugp-Umean(j+1)/avg_time_length) - (ugm-Umean(j-1)/avg_time_length) ) / ( dx2(j) + dx2(j-1) ))**dble(2.0) /dble((n3msub-1)*(n1msub-1))*dtime
                dudzdudz_tmp = 0.0d0
                dvdxdvdx_tmp = 0.0d0
                dvdydvdy_tmp = dvdydvdy_tmp+ ( ((vgp-Vmean(j+1)/avg_time_length) - (vgm-Vmean(j-1)/avg_time_length) ) / ( dx2(j) + dx2(j-1) ))**dble(2.0) /dble((n3msub-1)*(n1msub-1))*dtime
                dvdzdvdz_tmp = 0.0d0
                dwdxdwdx_tmp = 0.0d0
                dwdydwdy_tmp = dwdydwdy_tmp+ ( ((wgp-Wmean(j+1)/avg_time_length) - (wgm-Wmean(j-1)/avg_time_length) ) / ( dx2(j) + dx2(j-1) ))**dble(2.0) /dble((n3msub-1)*(n1msub-1))*dtime
                dwdzdwdz_tmp = 0.0d0
                dudxdvdx_tmp = 0.0d0
                dudydvdy_tmp = dudydvdy_tmp+   (((ugp-Umean(j+1)/avg_time_length) - (ugm-Umean(j-1)/avg_time_length) ) / ( dx2(j) + dx2(j-1) ))  &
                                                & *(((vgp-Vmean(j+1)/avg_time_length) - (vgm-Vmean(j-1)/avg_time_length) ) / ( dx2(j) + dx2(j-1) ))/dble((n3msub-1)*(n1msub-1))*dtime
                dudzdvdz_tmp = 0.0d0
                
                ! YMG: For V-P Gradient
                udpdx_tmp=udpdx_tmp + (ug-Umean(j)/avg_time_length)* ( (pgxp-Pmean(j  )/avg_time_length) - (pgxm-Pmean(j  )/avg_time_length) ) / ( dx1(i) + dx1(i-1) ) /dble((n3msub-1)*(n1msub-1))*dtime
                vdpdy_tmp=vdpdy_tmp + (vg-Vmean(j)/avg_time_length)* ( (pgyp-Pmean(j+1)/avg_time_length) - (pgym-Pmean(j-1)/avg_time_length) ) / ( dx2(j) + dx2(j-1) ) /dble((n3msub-1)*(n1msub-1))*dtime
                wdpdz_tmp=wdpdz_tmp + (wg-Wmean(j)/avg_time_length)* ( (pgzp-Pmean(j  )/avg_time_length) - (pgzm-Pmean(j  )/avg_time_length) ) / ( dx3(k) + dx3(k-1) ) /dble((n3msub-1)*(n1msub-1))*dtime
                udpdy_tmp=udpdy_tmp + (ug-Umean(j)/avg_time_length)* ( (pgyp-Pmean(j+1)/avg_time_length) - (pgym-Pmean(j-1)/avg_time_length) ) / ( dx2(j) + dx2(j-1) ) /dble((n3msub-1)*(n1msub-1))*dtime
                vdpdx_tmp=vdpdx_tmp + (vg-Vmean(j)/avg_time_length)* ( (pgxp-Pmean(j  )/avg_time_length) - (pgxm-Pmean(j  )/avg_time_length) ) / ( dx1(i) + dx1(i-1) ) /dble((n3msub-1)*(n1msub-1))*dtime
            enddo
            enddo
			dudxdudx(j) = 0.0d0
			dudydudy(j) = dudydudy(j)+ dudydudy_tmp
			dudzdudz(j) = 0.0d0
			dvdxdvdx(j) = 0.0d0
			dvdydvdy(j) = dvdydvdy(j)+ dvdydvdy_tmp
			dvdzdvdz(j) = 0.0d0
			dwdxdwdx(j) = 0.0d0
			dwdydwdy(j) = dwdydwdy(j)+ dwdydwdy_tmp
			dwdzdwdz(j) = 0.0d0
			dudxdvdx(j) = 0.0d0
			dudydvdy(j) = dudydvdy(j)+ dudydvdy_tmp
			dudzdvdz(j) = 0.0d0
			
			! YMG: For V-P Gradient
			udpdx(j)=udpdx(j) + udpdx_tmp
            vdpdy(j)=vdpdy(j) + vdpdy_tmp
            wdpdz(j)=wdpdz(j) + wdpdz_tmp
            udpdy(j)=udpdy(j) + udpdy_tmp
            vdpdx(j)=vdpdx(j) + vdpdx_tmp
        enddo

        ! !$acc parallel loop private(u1,u2,ug,v1,v2,vg,w1,w2,wg)
        ! do k=2,s3msb
        ! !$acc loop vector
        ! do j=1,s2msb
        ! !$acc loop vector
        ! do i=2,s1msb
        !     u1 = ( U(i,j,k-1)*dx2(j-1)+U(i,j-1,k-1)*dx2(j) )/( dx2(j)+dx2(j-1) )
        !     u2 = ( U(i,j,k  )*dx2(j-1)+U(i,j-1,k  )*dx2(j) )/( dx2(j)+dx2(j-1) )
        !     ug = ( u2*dx3(k-1)+u1*dx3(k) )/( dx3(k)+dx3(k-1) )

        !     v1 = ( V(i,j,k-1)*dx1(i-1)+V(i-1,j,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
        !     v2 = ( V(i,j,k  )*dx1(i-1)+V(i-1,j,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
        !     vg = ( v2*dx3(k-1)+v1*dx3(k) )/( dx3(k)+dx3(k-1) )
            
        !     w1 = ( W(i-1,j,k)*dx2(j-1)+W(i-1,j-1,k)*dx2(j) )/( dx2(j)+dx2(j-1) )
        !     w2 = ( W(i  ,j,k)*dx2(j-1)+W(i  ,j-1,k)*dx2(j) )/( dx2(j)+dx2(j-1) )
        !     wg = ( w2*dx1(i-1)+w1*dx1(i) )/( dx1(i)+dx1(i-1) )
            

        !     ! YMG: For Turbulent Diffusion
        !     UUUavg_t(j) = 0.0d0
        !     !$acc atomic update
        !     UUVavg_t(j) = UUVavg_t(j)+ (ug-Umean(j)/avg_time_length)*(ug-Umean(j)/avg_time_length)*(vg-Vmean(j)/avg_time_length)/dble((s3msb-1)*(s1msb-1))*dtime
        !     UUWavg_t(j) = 0.0d0
        !     VVUavg_t(j) = 0.0d0
        !     !$acc atomic update
        !     VVVavg_t(j) = VVVavg_t(j)+ (vg-Vmean(j)/avg_time_length)*(vg-Vmean(j)/avg_time_length)*(vg-Vmean(j)/avg_time_length)/dble((s3msb-1)*(s1msb-1))*dtime
        !     VVWavg_t(j) = 0.0d0
        !     WWUavg_t(j) = 0.0d0
        !     !$acc atomic update
        !     WWVavg_t(j) = WWVavg_t(j)+ (wg-Wmean(j)/avg_time_length)*(wg-Wmean(j)/avg_time_length)*(vg-Vmean(j)/avg_time_length)/dble((s3msb-1)*(s1msb-1))*dtime
        !     WWWavg_t(j) = 0.0d0
        !     UVUavg_t(j) = 0.0d0
        !     !$acc atomic update
        !     UVVavg_t(j) = UVVavg_t(j)+ (ug-Umean(j)/avg_time_length)*(vg-Vmean(j)/avg_time_length)*(vg-Vmean(j)/avg_time_length)/dble((s3msb-1)*(s1msb-1))*dtime
        !     UVWavg_t(j) = 0.0d0
        ! enddo
        ! enddo
        ! enddo

        do j=1,n2msub
            UUUavg_t_tmp=0.0d0; UUVavg_t_tmp=0.0d0; UUWavg_t_tmp=0.0d0; 
            VVUavg_t_tmp=0.0d0; VVVavg_t_tmp=0.0d0; VVWavg_t_tmp=0.0d0; 
            WWUavg_t_tmp=0.0d0; WWVavg_t_tmp=0.0d0; WWWavg_t_tmp=0.0d0; 
            UVUavg_t_tmp=0.0d0; UVVavg_t_tmp=0.0d0; UVWavg_t_tmp=0.0d0; 
            do k=2,n3msub
            do i=2,n1msub
                u1 = ( U(i,j,k-1)*dx2(j-1)+U(i,j-1,k-1)*dx2(j) )/( dx2(j)+dx2(j-1) )
                u2 = ( U(i,j,k  )*dx2(j-1)+U(i,j-1,k  )*dx2(j) )/( dx2(j)+dx2(j-1) )
                ug = ( u2*dx3(k-1)+u1*dx3(k) )/( dx3(k)+dx3(k-1) )
    
                v1 = ( V(i,j,k-1)*dx1(i-1)+V(i-1,j,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
                v2 = ( V(i,j,k  )*dx1(i-1)+V(i-1,j,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
                vg = ( v2*dx3(k-1)+v1*dx3(k) )/( dx3(k)+dx3(k-1) )
                
                w1 = ( W(i-1,j,k)*dx2(j-1)+W(i-1,j-1,k)*dx2(j) )/( dx2(j)+dx2(j-1) )
                w2 = ( W(i  ,j,k)*dx2(j-1)+W(i  ,j-1,k)*dx2(j) )/( dx2(j)+dx2(j-1) )
                wg = ( w2*dx1(i-1)+w1*dx1(i) )/( dx1(i)+dx1(i-1) )
                
                UUUavg_t_tmp = 0.0d0
                UUVavg_t_tmp = UUVavg_t_tmp+ (ug-Umean(j)/avg_time_length)*(ug-Umean(j)/avg_time_length)*(vg-Vmean(j)/avg_time_length)/dble((n3msub-1)*(n1msub-1))*dtime
                UUWavg_t_tmp = 0.0d0
                VVUavg_t_tmp = 0.0d0
                VVVavg_t_tmp = VVVavg_t_tmp+ (vg-Vmean(j)/avg_time_length)*(vg-Vmean(j)/avg_time_length)*(vg-Vmean(j)/avg_time_length)/dble((n3msub-1)*(n1msub-1))*dtime
                VVWavg_t_tmp = 0.0d0
                WWUavg_t_tmp = 0.0d0
                WWVavg_t_tmp = WWVavg_t_tmp+ (wg-Wmean(j)/avg_time_length)*(wg-Wmean(j)/avg_time_length)*(vg-Vmean(j)/avg_time_length)/dble((n3msub-1)*(n1msub-1))*dtime
                WWWavg_t_tmp = 0.0d0
                UVUavg_t_tmp = 0.0d0
                UVVavg_t_tmp = UVVavg_t_tmp+ (ug-Umean(j)/avg_time_length)*(vg-Vmean(j)/avg_time_length)*(vg-Vmean(j)/avg_time_length)/dble((n3msub-1)*(n1msub-1))*dtime
                UVWavg_t_tmp = 0.0d0
            enddo
            enddo
            ! YMG: For Turbulent Diffusion
            UUUavg_t(j) = 0.0d0
            UUVavg_t(j) = UUVavg_t(j)+ UUVavg_t_tmp
            UUWavg_t(j) = 0.0d0
            VVUavg_t(j) = 0.0d0
            VVVavg_t(j) = VVVavg_t(j)+ VVVavg_t_tmp
            VVWavg_t(j) = 0.0d0
            WWUavg_t(j) = 0.0d0
            WWVavg_t(j) = WWVavg_t(j)+ WWVavg_t_tmp
            WWWavg_t(j) = 0.0d0
            UVUavg_t(j) = 0.0d0
            UVVavg_t(j) = UVVavg_t(j)+ UVVavg_t_tmp
            UVWavg_t(j) = 0.0d0

            enddo

    end subroutine mpi_statistics_avg_xzt

    subroutine mpi_statistics_Reynolds_Budget_out(myrank) ! YMG 210407
        
        use mpi
        use mpi_topology,   only : comm_1d_x1, comm_1d_x2, comm_1d_x3
        use mpi_subdomain,  only : x2_sub
        use mpi_subdomain,  only : dx1_sub, dx2_sub, dx3_sub, dmx1_sub, dmx2_sub, dmx3_sub
        use mpi_subdomain,  only : n1sub, n2sub, n3sub, n1msub, n2msub, n3msub
 
        implicit none
        
        integer :: i,j,k
        integer :: myrank
        integer :: onebyone,ierr
		
        double precision, dimension(:), pointer:: dx1, dx2, dx3
        double precision, dimension(:), pointer :: dmx1, dmx2, dmx3
        double precision :: RS_All(0:n2sub,4), RS_P(0:n2sub,4), RS(0:n2sub,4), RS_T(0:n2sub,4), RS_Pi(0:n2sub,4), RS_V(0:n2sub,4) ! 1:UU 2:VV 3:WW 4:UV

        ! s1msb = s1sb-1
        ! s2msb = s2sb-1
        ! s3msb = s3sb-1

        dx1  => dx1_sub
        dx2  => dx2_sub
        dx3  => dx3_sub
        dmx1 => dmx1_sub
        dmx2 => dmx2_sub
        dmx3 => dmx3_sub
        
		! Reynolds_Stress_Budgets UU
        do j=1,n2sub
			RS_P (j,1)   =  -dble(2.0)  *  (UVavg_t(j)/avg_time_length)  *  ( (Umean(j+1)/avg_time_length) - (Umean(j-1)/avg_time_length) ) / (  dx2(j) + dx2(j-1)  )
			RS (j,1)   =  -dble(2.0)  * Cmu * dudydudy(j)/avg_time_length
			RS_T (j,1)   =  -( (UUVavg_t(j+1)/avg_time_length) - (UUVavg_t(j-1)/avg_time_length) ) / (  dx2(j) + dx2(j-1)  ) 
			RS_Pi(j,1)   =  -dble(2.0)  * udpdx(j)/avg_time_length
			RS_V (j,1)   =   Cmu    * (  ( (UUavg_t(j+1)/avg_time_length) - (UUavg_t(j  )/avg_time_length) ) / dx2(j)  -             &
		                            &    ( (UUavg_t(j  )/avg_time_length) - (UUavg_t(j-1)/avg_time_length) ) / dx2(j-1)    ) / dmx2(j)
										   
			RS_ALL(j,1)  =  RS_P(j,1)  +  RS(j,1)  +  RS_T(j,1)  +  RS_Pi(j,1)  +  RS_V(j,1)	
        enddo
		
		! Reynolds_Stress_Budgets VV
        do j=1,n2sub
			RS_P (j,2)   =  -dble(2.0)  *  (VVavg_t(j)/avg_time_length)  *  ( (Vmean(j+1)/avg_time_length) - (Vmean(j-1)/avg_time_length) ) / (  dx2(j) + dx2(j-1)  )
			RS (j,2)   =  -dble(2.0)  * Cmu * dvdydvdy(j)/avg_time_length
			RS_T (j,2)   =  -( (VVVavg_t(j+1)/avg_time_length) - (VVVavg_t(j-1)/avg_time_length) ) / (  dx2(j) + dx2(j-1)  ) 
			RS_Pi(j,2)   =  -dble(2.0)  * vdpdy(j)/avg_time_length
			RS_V (j,2)   =   Cmu    * (  ( (VVavg_t(j+1)/avg_time_length) - (VVavg_t(j  )/avg_time_length) ) / dx2(j)  -             &
		                            &    ( (VVavg_t(j  )/avg_time_length) - (VVavg_t(j-1)/avg_time_length) ) / dx2(j-1)    ) / dmx2(j)
										   
			RS_ALL(j,2)  =  RS_P(j,2)  +  RS(j,2)  +  RS_T(j,2)  +  RS_Pi(j,2)  +  RS_V(j,2)	
        enddo
		
		! Reynolds_Stress_Budgets WW
        do j=1,n2sub
			RS_P (j,3)   =  -dble(2.0)  *  (VWavg_t(j)/avg_time_length)  *  ( (Wmean(j+1)/avg_time_length) - (Wmean(j-1)/avg_time_length) ) / (  dx2(j) + dx2(j-1)  )
			RS (j,3)   =  -dble(2.0)  * Cmu * dwdydwdy(j)/avg_time_length
			RS_T (j,3)   =  -( (WWVavg_t(j+1)/avg_time_length) - (WWVavg_t(j-1)/avg_time_length) ) / (  dx2(j) + dx2(j-1)  ) 
			RS_Pi(j,3)   =  -dble(2.0)  * wdpdz(j)/avg_time_length
			RS_V (j,3)   =   Cmu    * (  ( (WWavg_t(j+1)/avg_time_length) - (WWavg_t(j  )/avg_time_length) ) / dx2(j)  -             &
		                            &    ( (WWavg_t(j  )/avg_time_length) - (WWavg_t(j-1)/avg_time_length) ) / dx2(j-1)    ) / dmx2(j)
										   
			RS_ALL(j,3)  =  RS_P(j,3)  +  RS(j,3)  +  RS_T(j,3)  +  RS_Pi(j,3)  +  RS_V(j,3)	
        enddo

		! Reynolds_Stress_Budgets UV
        do j=1,n2sub
			RS_P (j,4)   =  - (UVavg_t(j)/avg_time_length)  *  ( (Vmean(j+1)/avg_time_length) - (Vmean(j-1)/avg_time_length) ) / (  dx2(j) + dx2(j-1)  )  &
			              & - (VVavg_t(j)/avg_time_length)  *  ( (Umean(j+1)/avg_time_length) - (Umean(j-1)/avg_time_length) ) / (  dx2(j) + dx2(j-1)  )
						  
			RS (j,4)   =  -dble(2.0)  * Cmu * dudydvdy(j)/avg_time_length
			RS_T (j,4)   =  -( (UVVavg_t(j+1)/avg_time_length) - (UVVavg_t(j-1)/avg_time_length) ) / (  dx2(j) + dx2(j-1)  ) 
			RS_Pi(j,4)   =  -1.0d0  * udpdy(j)/avg_time_length - 1.0d0  * vdpdx(j)/avg_time_length
			RS_V (j,4)   =   Cmu    * (  ( (UVavg_t(j+1)/avg_time_length) - (UVavg_t(j  )/avg_time_length) ) / dx2(j)  -             &
		                            &    ( (UVavg_t(j  )/avg_time_length) - (UVavg_t(j-1)/avg_time_length) ) / dx2(j-1)    ) / dmx2(j)
										   
			RS_ALL(j,4)  =  RS_P(j,4)  +  RS(j,4)  +  RS_T(j,4)  +  RS_Pi(j,4)  +  RS_V(j,4)	
        enddo


		
		! Reynolds_Stress_Budgets Fileout : It is normalized by "u_tau^4/nu"
        do onebyone=0,comm_1d_x2%nprocs-1
            if(onebyone==comm_1d_x2%myrank) then
                open(myrank, file=dir_statistics//'051_Reynolds_Budget_UU'//xyzrank//'.plt')  
                    write(myrank,*) 'variables="y","RS_All","RS_P","RS","RS_T","RS_Pi","RS_V"'
                    do j=2,n2msub
                        write(myrank,'(7E16.8)') x2_sub(j),            RS_All(j,1)  &
																	&, RS_P  (j,1)  &
																	&, RS  (j,1)  &
																	&, RS_T  (j,1)  &
																	&, RS_Pi (j,1)  &
																	&, RS_V  (j,1) 
                    enddo
                close(myrank)
                open(myrank, file=dir_statistics//'052_Reynolds_Budget_VV'//xyzrank//'.plt')
                    write(myrank,*) 'variables="y","RS_All","RS_P","RS","RS_T","RS_Pi","RS_V"'
                    do j=2,n2msub
                        write(myrank,'(7E16.8)') x2_sub(j),            RS_All(j,2)  &
																	&, RS_P  (j,2)  &
																	&, RS  (j,2)  &
																	&, RS_T  (j,2)  &
																	&, RS_Pi (j,2)  &
																	&, RS_V  (j,2) 
                    enddo	
                close(myrank)
                open(myrank, file=dir_statistics//'053_Reynolds_Budget_WW'//xyzrank//'.plt')
                    write(myrank,*) 'variables="y","RS_All","RS_P","RS","RS_T","RS_Pi","RS_V"'
                    do j=2,n2msub
                        write(myrank,'(7E16.8)') x2_sub(j),            RS_All(j,3)  &
																	&, RS_P  (j,3)  &
																	&, RS  (j,3)  &
																	&, RS_T  (j,3)  &
																	&, RS_Pi (j,3)  &
																	&, RS_V  (j,3) 
                    enddo
                close(myrank)
                open(myrank, file=dir_statistics//'054_Reynolds_Budget_UV'//xyzrank//'.plt')
                    write(myrank,*) 'variables="y","-RS_All","-RS_P","-RS","-RS_T","-RS_Pi","-RS_V"'
                    do j=2,n2msub
                        write(myrank,'(7E16.8)') x2_sub(j),           -RS_All(j,4)  &
																	&,-RS_P  (j,4)  &
																	&,-RS  (j,4)  &
																	&,-RS_T  (j,4)  &
																	&,-RS_Pi (j,4)  &
																	&,-RS_V  (j,4) 
                    enddo
                close(myrank)				
			endif
            call MPI_Barrier(MPI_COMM_WORLD,ierr)
        enddo		
		

    end subroutine mpi_statistics_Reynolds_Budget_out
	
    subroutine mpi_statistics_fileout(myrank, SR)
        use mpi
        use mpi_topology,   only : comm_1d_x2
        use mpi_subdomain,  only : x2_sub
        use mpi_subdomain,  only : n2sub
        use mpi_subdomain,  only : n2msub

        implicit none
        integer :: i,j,k
        integer :: myrank
        integer :: onebyone,ierr
        double precision :: SR(0:n2sub)

        do onebyone=0,comm_1d_x2%nprocs-1
            if(onebyone==comm_1d_x2%myrank) then
                open(myrank, file=dir_statistics//'01_mean_profiles_'//xyzrank//'.plt')
                    write(myrank,*) 'variables="y","Um","Vm","Wm","Tm","Pm"'
                    do j=1,n2sub
                        write(myrank,'(6E16.8)') x2_sub(j),            Umean(j)/avg_time_length    &
																	&, Vmean(j)/avg_time_length    &
																	&, Wmean(j)/avg_time_length    &
																	&, Tmean(j)/avg_time_length    &
																	&, Pmean(j)/avg_time_length
                        !write(myrank,'(5D16.8)') x2_sub(j),Umean(j) &
                        !                                 &,Vmean(j) &
                        !                                 &,Wmean(j) &
                        !                                 &,Tmean(j)
                    enddo
                close(myrank)
                open(myrank, file=dir_statistics//'02_Stresses_profiles_'//xyzrank//'.plt')
                    write(myrank,*) 'variables="y","u_std","v_std","w_std","uv","TT","UT"'
                    do j=1,n2sub
                        write(myrank,'(7E16.8)') x2_sub(j),            sqrt(abs(UUavg_t(j))/avg_time_length) &
																	&, sqrt(abs(VVavg_t(j))/avg_time_length) &
																	&, sqrt(abs(WWavg_t(j))/avg_time_length) &
																	&,      abs(UVavg_t(j))/avg_time_length  &
																	&, sqrt(abs(TTavg_t(j))/avg_time_length) &
																	&,         (UTavg_t(j))/avg_time_length
                        !write(myrank,'(7D16.8)') x2_sub(j),UUavg_t(j) &
                        !                                 &,VVavg_t(j) &
                        !                                 &,WWavg_t(j) &
                        !                                 &,UVavg_t(j) &
                        !                                 &,TTavg_t(j) &
                        !                                 &,UTavg_t(j)
                    enddo
                close(myrank)
                open(myrank, file=dir_statistics//'03_Stresses_profiles_02_'//xyzrank//'.plt')
                    write(myrank,*) 'variables="y","Usq","Vsq","Wsq","UV","TT","UT"'
                    do j=1,n2sub
                        write(myrank,'(7E16.8)') x2_sub(j),            Usq(j)/avg_time_length    &
                                                                    &, Vsq(j)/avg_time_length &
                                                                    &, Wsq(j)/avg_time_length   &
																	&, UV(j) /avg_time_length   &
																	&, Tsq(j)/avg_time_length   &
																	&, UT(j) /avg_time_length
                    enddo
                close(myrank)
                open(myrank, file=dir_statistics//'04_Strainrate_profile_'//xyzrank//'.plt')
                    write(myrank,*) 'variables="y","absSR"'
                    do j=1,n2msub
                        write(myrank,'(2E16.8)') 0.5*(x2_sub(j)+x2_sub(j+1)), SR(j)/avg_time_length 
                    enddo
                close(myrank)
            endif
            call MPI_Barrier(MPI_COMM_WORLD,ierr)
        enddo

    end subroutine mpi_statistics_fileout

    ! subroutine mpi_statistics_avg_xzt_3D(U,V,W,P,T,dx1,dx2,dx3,dmx1,dmx2,dmx3,dtime)
    !     implicit none
    !     integer :: i,j,k
    !     double precision, device, intent(in), dimension(0:s1sb,0:s2sb,0:s3sb) :: U, V, W, P, T
    !     double precision, device ::  dx1(0:s1sb), dx2(0:s2sb), dx3(0:s3sb)
    !     double precision, device :: dmx1(0:s1sb),dmx2(0:s2sb),dmx3(0:s3sb)
    !     double precision :: dtime
    !     double precision :: u1,v1,w1
    !     double precision :: u2,v2,w2
    !     double precision :: t1,t2,t3,t4,t5,t6
    !     double precision :: p1,p2,p3,p4,p5,p6
    !     double precision :: ug,vg,wg,tg,pg
	! 	double precision :: ugp, ugm, vgp, vgm, wgp, wgm, pgyp, pgym, pgxp, pgxm, pgzp, pgzm

    !     s1msb = s1sb-1
    !     s2msb = s2sb-1
    !     s3msb = s3sb-1

    !     avg_time_length=avg_time_length+dtime

    !     $acc kernels
    !     $acc loop
    !     do k=1,s3msb
    !     do j=1,s2msb
    !     $acc loop
    !     do j=1,s2sb
    !     $acc loop
    !     do i=1,s1msb
        
    !         Umean(1:s1msb,j,k)=Umean(1:s1msb,j,k)+U(1:s1msb,j,k)*dtime/dble(s2msb*s3msb)
    !             $acc atomic update
    !             Umean(j) = Umean(j)    + ( U(i,j,k)*dx2(j-1)+U(i,j-1,k)*dx2(j) )/(dx2(j)+dx2(j-1)) /dble(s3msb*s1msb)*dtime
    !             $acc atomic update
    !             Vmean(j) = Vmean(j)    +   V(i,j,k)/dble(s3msb*s1msb)*dtime
    !             $acc atomic update
    !             Wmean(j) = Wmean(j)    + ( W(i,j,k)*dx2(j-1)+W(i,j-1,k)*dx2(j) )/(dx2(j)+dx2(j-1)) /dble(s3msb*s1msb)*dtime
    !             $acc atomic update
    !             Tmean(j) = Tmean(j)    + ( T(i,j,k)*dx2(j-1)+T(i,j-1,k)*dx2(j) )/(dx2(j)+dx2(j-1)) /dble(s3msb*s1msb)*dtime
    !             $acc atomic update
	! 			Pmean(j) = Pmean(j)    + ( P(i,j,k)*dx2(j-1)+P(i,j-1,k)*dx2(j) )/(dx2(j)+dx2(j-1)) /dble(s3msb*s1msb)*dtime
    !     enddo
    !     enddo
    !     enddo
    !     $acc end kernels
		

	! 	YMG: Non-dimensionlized Parameter Calculation
	!        u_tau =  SQRT( Cmu * ( Umean(2)/avg_time_length - Umean(1)/avg_time_length ) / ( dx2(1) ) )
	! 	delta_nu =  Cmu / u_tau
	! 	Temp_tau =  - ( Tmean(2) - Tmean(1) ) * Ct /  u_tau


    !     $acc kernels
    !     $acc loop private(u1,u2,ug,v1,v2,vg,w1,w2,wg,t1,t2,t3,t4,t5,t6,tg,p1,p2,p3,p4,p5,p6,pg)
    !     do k=1,s3sb
    !     $acc loop
    !     do j=1,s2sb
    !     $acc loop
    !     do i=1,s1sb
    !         u1 = ( U(i,j,k-1)*dx2(j-1)+U(i,j-1,k-1)*dx2(j) )/( dx2(j)+dx2(j-1) )
    !         u2 = ( U(i,j,k  )*dx2(j-1)+U(i,j-1,k  )*dx2(j) )/( dx2(j)+dx2(j-1) )
    !         ug = ( u2*dx3(k-1)+u1*dx3(k) )/( dx3(k)+dx3(k-1) )
    
    !         v1 = ( V(i,j,k-1)*dx1(i-1)+V(i-1,j,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
    !         v2 = ( V(i,j,k  )*dx1(i-1)+V(i-1,j,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
    !         vg = ( v2*dx3(k-1)+v1*dx3(k) )/( dx3(k)+dx3(k-1) )

    !         w1 = ( W(i-1,j,k)*dx2(j-1)+W(i-1,j-1,k)*dx2(j) )/( dx2(j)+dx2(j-1) )
    !         w2 = ( W(i  ,j,k)*dx2(j-1)+W(i  ,j-1,k)*dx2(j) )/( dx2(j)+dx2(j-1) )
    !         wg = ( w2*dx1(i-1)+w1*dx1(i) )/( dx1(i)+dx1(i-1) )

    !         t1 = ( T(i  ,j  ,k  )*dx1(i-1)+T(i-1,j  ,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
    !         t2 = ( T(i  ,j-1,k  )*dx1(i-1)+T(i-1,j-1,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
    !         t3 = ( T(i  ,j  ,k-1)*dx1(i-1)+T(i-1,j  ,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
    !         t4 = ( T(i  ,j-1,k-1)*dx1(i-1)+T(i-1,j-1,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
    
    !         t5 = ( t1*dx2(j-1)+t2*dx2(j) )/( dx2(j)+dx2(j-1) )
    !         t6 = ( t3*dx2(j-1)+t4*dx2(j) )/( dx2(j)+dx2(j-1) )
    !         tg = ( t5*dx3(k-1)+t6*dx3(k) )/( dx3(k)+dx3(k-1) )
            
    !         p1 = ( P(i  ,j  ,k  )*dx1(i-1)+P(i-1,j  ,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
    !         p2 = ( P(i  ,j-1,k  )*dx1(i-1)+P(i-1,j-1,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
    !         p3 = ( P(i  ,j  ,k-1)*dx1(i-1)+P(i-1,j  ,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
    !         p4 = ( P(i  ,j-1,k-1)*dx1(i-1)+P(i-1,j-1,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
    
    !         p5 = ( p1*dx2(j-1)+p2*dx2(j) )/( dx2(j)+dx2(j-1) )
    !         p6 = ( p3*dx2(j-1)+p4*dx2(j) )/( dx2(j)+dx2(j-1) )
    !         pg = ( p5*dx3(k-1)+p6*dx3(k) )/( dx3(k)+dx3(k-1) )

    !         $acc atomic update
    !         UUavg_t(j) =  UUavg_t(j)+ (ug-Umean(j)/avg_time_length)*(ug-Umean(j)/avg_time_length)/dble(s3sb*s1sb)*dtime
    !         $acc atomic update
    !         VVavg_t(j) =  VVavg_t(j)+ (vg-Vmean(j)/avg_time_length)*(vg-Vmean(j)/avg_time_length)/dble(s3sb*s1sb)*dtime
    !         $acc atomic update
    !         WWavg_t(j) =  WWavg_t(j)+ (wg-Wmean(j)/avg_time_length)*(wg-Wmean(j)/avg_time_length)/dble(s3sb*s1sb)*dtime
    !         $acc atomic update
    !         UVavg_t(j) =  UVavg_t(j)+ (ug-Umean(j)/avg_time_length)*(vg-Vmean(j)/avg_time_length)/dble(s3sb*s1sb)*dtime
    !         $acc atomic update
    !         TTavg_t(j) =  TTavg_t(j)+ (tg-Tmean(j)/avg_time_length)*(tg-Tmean(j)/avg_time_length)/dble(s3sb*s1sb)*dtime
    !         $acc atomic update
    !         UTavg_t(j) =  UTavg_t(j)+ (ug-Umean(j)/avg_time_length)*(tg-Tmean(j)/avg_time_length)/dble(s3sb*s1sb)*dtime
    !         $acc atomic update
    !         UWavg_t(j) =  UWavg_t(j)+ (ug-Umean(j)/avg_time_length)*(wg-Wmean(j)/avg_time_length)/dble(s3sb*s1sb)*dtime
    !         $acc atomic update
    !         VWavg_t(j) =  VWavg_t(j)+ (vg-Vmean(j)/avg_time_length)*(wg-Wmean(j)/avg_time_length)/dble(s3sb*s1sb)*dtime
			
    !         Usq(j)=Usq(j) + ug*ug/dble(s3sb*s1sb)*dtime
    !         Vsq(j)=Vsq(j) + vg*vg/dble(s3sb*s1sb)*dtime
    !         Wsq(j)=Wsq(j) + wg*wg/dble(s3sb*s1sb)*dtime
    !         Tsq(j)=Tsq(j) + tg*tg/dble(s3sb*s1sb)*dtime
    !         UV(j) =UV(j)  + ug*vg/dble(s3sb*s1sb)*dtime
    !         UT(j) =UT(j)  + ug*tg/dble(s3sb*s1sb)*dtime

    !     enddo
    !     enddo
    !     enddo
    !     $acc end kernels

    ! end subroutine mpi_statistics_avg_xzt_3D
    
    ! subroutine mpi_statistics_Reynolds_budget_3D(U,V,W,P,T,dx1,dx2,dx3,dmx1,dmx2,dmx3,dtime)
    !     use MPI
    !     use mpi_subdomain

    !     implicit none
    !     integer :: i,j,k, iter
    !     double precision, device, dimension(0:s1sb,0:s2sb,0:s3sb) :: U,V,W,P,T
    !     double precision, device :: dx1(0:s1sb),dx2(0:s2sb),dx3(0:s3sb)
    !     double precision, device :: dmx1(0:s1sb),dmx2(0:s2sb),dmx3(0:s3sb)
    !     double precision :: dtime
    !     double precision :: u1, u2, v1, v2, w1, w2
    !     double precision :: p1, p2, p3, p4, p5, p6
    !     double precision :: Ug, Vg, Wg, Pg

    !     double precision, device, allocatable, dimension (:,:,:) :: Uf, Vf, Wf, Pf, Temp

    !     allocate(Uf(0:s1sb,0:s2sb,0:s3sb), Vf(0:s1sb,0:s2sb,0:s3sb), Wf(0:s1sb,0:s2sb,0:s3sb), Pf(0:s1sb,0:s2sb,0:s3sb))
    !     allocate(Temp(0:s1sb,0:s2sb,0:s3sb))
    !     Uf=0.0d0; Vf=0.0d0; Wf=0.0d0; Pf=0.0d0; Temp=0.0d0

    !     s1msb = s1sb-1
    !     s2msb = s2sb-1
    !     s3msb = s3sb-1

    !     avg_time_length=avg_time_length+dtime

    !     $acc parallel loop collapse(3) reduction(+:Ugmean, Vgmean, Wgmean, Pgmean, Uf, Vf, Wf, Pf) private(u1,u2,ug,v1,v2,vg,w1,w2,wg,p1,p2,p3,p4,p5,p6,pg)
    !     do k=1,s3sb
    !     do j=1,s2sb
    !     do i=1,s1sb
    !         u1 = ( U(i,j,k-1)*dx2(j-1)+U(i,j-1,k-1)*dx2(j) )/( dx2(j)+dx2(j-1) )
    !         u2 = ( U(i,j,k  )*dx2(j-1)+U(i,j-1,k  )*dx2(j) )/( dx2(j)+dx2(j-1) )
    !         ug = ( u2*dx3(k-1)+u1*dx3(k) )/( dx3(k)+dx3(k-1) )
    
    !         v1 = ( V(i,j,k-1)*dx1(i-1)+V(i-1,j,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
    !         v2 = ( V(i,j,k  )*dx1(i-1)+V(i-1,j,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
    !         vg = ( v2*dx3(k-1)+v1*dx3(k) )/( dx3(k)+dx3(k-1) )

    !         w1 = ( W(i-1,j,k)*dx2(j-1)+W(i-1,j-1,k)*dx2(j) )/( dx2(j)+dx2(j-1) )
    !         w2 = ( W(i  ,j,k)*dx2(j-1)+W(i  ,j-1,k)*dx2(j) )/( dx2(j)+dx2(j-1) )
    !         wg = ( w2*dx1(i-1)+w1*dx1(i) )/( dx1(i)+dx1(i-1) )

    !         p1 = ( P(i  ,j  ,k  )*dx1(i-1)+P(i-1,j  ,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
    !         p2 = ( P(i  ,j-1,k  )*dx1(i-1)+P(i-1,j-1,k  )*dx1(i) )/( dx1(i)+dx1(i-1) )
    !         p3 = ( P(i  ,j  ,k-1)*dx1(i-1)+P(i-1,j  ,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
    !         p4 = ( P(i  ,j-1,k-1)*dx1(i-1)+P(i-1,j-1,k-1)*dx1(i) )/( dx1(i)+dx1(i-1) )
    
    !         p5 = ( p1*dx2(j-1)+p2*dx2(j) )/( dx2(j)+dx2(j-1) )
    !         p6 = ( p3*dx2(j-1)+p4*dx2(j) )/( dx2(j)+dx2(j-1) )
    !         pg = ( p5*dx3(k-1)+p6*dx3(k) )/( dx3(k)+dx3(k-1) )
            
    !         Ugmean(i,j,k) = Ugmean(i,j,k) + ( Ug )*dtime
    !         Vgmean(i,j,k) = Vgmean(i,j,k) + ( Vg )*dtime
    !         Wgmean(i,j,k) = Wgmean(i,j,k) + ( Wg )*dtime
    !         Pgmean(i,j,k) = Pgmean(i,j,k) + ( Pg )*dtime

    !         Uf(i,j,k) = Ug - ( Ugmean(i,j,k) / avg_time_length )
    !         Vf(i,j,k) = Vg - ( Vgmean(i,j,k) / avg_time_length )
    !         Wf(i,j,k) = Wg - ( Wgmean(i,j,k) / avg_time_length )
    !         Pf(i,j,k) = Pg - ( Pgmean(i,j,k) / avg_time_length )

    !     enddo
    !     enddo
    !     enddo
    !     $acc end parallel

    !     call cuda_mpi_subdomain_ghostcell_update(Uf)
    !     call cuda_mpi_subdomain_ghostcell_update(Vf)
    !     call cuda_mpi_subdomain_ghostcell_update(Wf)
    !     call cuda_mpi_subdomain_ghostcell_update(Pf)

    !     Mean Calculation

    !     $acc kernels
    !     $acc loop
    !     do k=1,s3msb
    !     $acc loop
    !     do j=1,s2msb
    !     $acc loop
    !     do i=1,s1msb

    !         $acc atomic update
    !         UU_mean3D(1,i,j,k) = UU_mean3D(1,i,j,k) + ( Uf(i,j,k) ) * ( Uf(i,j,k) ) * dtime  ! UU
    !         $acc atomic update
    !         UU_mean3D(2,i,j,k) = UU_mean3D(2,i,j,k) + ( Vf(i,j,k) ) * ( Vf(i,j,k) ) * dtime  ! VV
    !         $acc atomic update
    !         UU_mean3D(3,i,j,k) = UU_mean3D(3,i,j,k) + ( Wf(i,j,k) ) * ( Wf(i,j,k) ) * dtime  ! WW
    !         $acc atomic update
    !         UU_mean3D(4,i,j,k) = UU_mean3D(4,i,j,k) + ( Uf(i,j,k) ) * ( Vf(i,j,k) ) * dtime  ! UV
    !         $acc atomic update
    !         UU_mean3D(5,i,j,k) = UU_mean3D(5,i,j,k) + ( Uf(i,j,k) ) * ( Wf(i,j,k) ) * dtime  ! UW
    !         $acc atomic update
    !         UU_mean3D(6,i,j,k) = UU_mean3D(6,i,j,k) + ( Vf(i,j,k) ) * ( Wf(i,j,k) ) * dtime  ! VW

    !         $acc atomic update
    !         dUdxdUdx_mean3D(1 ,i,j,k) =  dUdxdUdx_mean3D(1 ,i,j,k) + ( ((Uf(i+1,j,k)-Uf(i-1,j,k))/(dx1(i)+dx1(i-1))) **dble(2.0) ) * dtime  ! dUdxdUdx
    !         $acc atomic update
    !         dUdxdUdx_mean3D(2 ,i,j,k) =  dUdxdUdx_mean3D(2 ,i,j,k) + ( ((Uf(i,j+1,k)-Uf(i,j-1,k))/(dx2(j)+dx2(j-1))) **dble(2.0) ) * dtime  ! dUdydUdy
    !         $acc atomic update
    !         dUdxdUdx_mean3D(3 ,i,j,k) =  dUdxdUdx_mean3D(3 ,i,j,k) + ( ((Uf(i,j,k+1)-Uf(i,j,k-1))/(dx3(k)+dx3(k-1))) **dble(2.0) ) * dtime  ! dUdzdUdz
    !         $acc atomic update
    !         dUdxdUdx_mean3D(4 ,i,j,k) =  dUdxdUdx_mean3D(4 ,i,j,k) + ( ((Vf(i+1,j,k)-Vf(i-1,j,k))/(dx1(i)+dx1(i-1))) **dble(2.0) ) * dtime  ! dVdxdVdx
    !         $acc atomic update
    !         dUdxdUdx_mean3D(5 ,i,j,k) =  dUdxdUdx_mean3D(5 ,i,j,k) + ( ((Vf(i,j+1,k)-Vf(i,j-1,k))/(dx2(j)+dx2(j-1))) **dble(2.0) ) * dtime  ! dVdydVdy
    !         $acc atomic update
    !         dUdxdUdx_mean3D(6 ,i,j,k) =  dUdxdUdx_mean3D(6 ,i,j,k) + ( ((Vf(i,j,k+1)-Vf(i,j,k-1))/(dx3(k)+dx3(k-1))) **dble(2.0) ) * dtime  ! dVdzdVdz
    !         $acc atomic update
    !         dUdxdUdx_mean3D(7 ,i,j,k) =  dUdxdUdx_mean3D(7 ,i,j,k) + ( ((Wf(i+1,j,k)-Wf(i-1,j,k))/(dx1(i)+dx1(i-1))) **dble(2.0) ) * dtime  ! dWdxdWdx
    !         $acc atomic update
    !         dUdxdUdx_mean3D(8 ,i,j,k) =  dUdxdUdx_mean3D(8 ,i,j,k) + ( ((Wf(i,j+1,k)-Wf(i,j-1,k))/(dx2(j)+dx2(j-1))) **dble(2.0) ) * dtime  ! dWdydWdy
    !         $acc atomic update
    !         dUdxdUdx_mean3D(9 ,i,j,k) =  dUdxdUdx_mean3D(9 ,i,j,k) + ( ((Wf(i,j,k+1)-Wf(i,j,k-1))/(dx3(k)+dx3(k-1))) **dble(2.0) ) * dtime  ! dWdzdWdz
    !         $acc atomic update
    !         dUdxdUdx_mean3D(10,i,j,k) =  dUdxdUdx_mean3D(10,i,j,k) + ( ((Uf(i+1,j,k)-Uf(i-1,j,k))/(dx1(i)+dx1(i-1))) * ((Vf(i+1,j,k)-Vf(i-1,j,k))/(dx1(i)+dx1(i-1))) ) * dtime  ! dUdxdVdx
    !         $acc atomic update
    !         dUdxdUdx_mean3D(11,i,j,k) =  dUdxdUdx_mean3D(11,i,j,k) + ( ((Uf(i,j+1,k)-Uf(i,j-1,k))/(dx2(j)+dx2(j-1))) * ((Vf(i,j+1,k)-Vf(i,j-1,k))/(dx2(j)+dx2(j-1))) ) * dtime  ! dUdydVdy
    !         $acc atomic update
    !         dUdxdUdx_mean3D(12,i,j,k) =  dUdxdUdx_mean3D(12,i,j,k) + ( ((Uf(i,j,k+1)-Uf(i,j,k-1))/(dx3(k)+dx3(k-1))) * ((Vf(i,j,k+1)-Vf(i,j,k-1))/(dx3(k)+dx3(k-1))) ) * dtime  ! dUdzdVdz

    !         $acc atomic update
    !         UUU_mean3D(1 ,i,j,k) = UUU_mean3D(1 ,i,j,k) + ( Uf(i,j,k) ) * ( Uf(i,j,k) ) * ( Uf(i,j,k) ) * dtime  ! UUU
    !         $acc atomic update
    !         UUU_mean3D(2 ,i,j,k) = UUU_mean3D(2 ,i,j,k) + ( Uf(i,j,k) ) * ( Uf(i,j,k) ) * ( Vf(i,j,k) ) * dtime  ! UUV
    !         $acc atomic update
    !         UUU_mean3D(3 ,i,j,k) = UUU_mean3D(3 ,i,j,k) + ( Uf(i,j,k) ) * ( Uf(i,j,k) ) * ( Wf(i,j,k) ) * dtime  ! UUW
    !         $acc atomic update
    !         UUU_mean3D(4 ,i,j,k) = UUU_mean3D(4 ,i,j,k) + ( Vf(i,j,k) ) * ( Vf(i,j,k) ) * ( Uf(i,j,k) ) * dtime  ! VVU
    !         $acc atomic update
    !         UUU_mean3D(5 ,i,j,k) = UUU_mean3D(5 ,i,j,k) + ( Vf(i,j,k) ) * ( Vf(i,j,k) ) * ( Vf(i,j,k) ) * dtime  ! VVV
    !         $acc atomic update
    !         UUU_mean3D(6 ,i,j,k) = UUU_mean3D(6 ,i,j,k) + ( Vf(i,j,k) ) * ( Vf(i,j,k) ) * ( Wf(i,j,k) ) * dtime  ! VVW
    !         $acc atomic update
    !         UUU_mean3D(7 ,i,j,k) = UUU_mean3D(7 ,i,j,k) + ( Wf(i,j,k) ) * ( Wf(i,j,k) ) * ( Uf(i,j,k) ) * dtime  ! WWU
    !         $acc atomic update
    !         UUU_mean3D(8 ,i,j,k) = UUU_mean3D(8 ,i,j,k) + ( Wf(i,j,k) ) * ( Wf(i,j,k) ) * ( Vf(i,j,k) ) * dtime  ! WWV
    !         $acc atomic update
    !         UUU_mean3D(9 ,i,j,k) = UUU_mean3D(9 ,i,j,k) + ( Wf(i,j,k) ) * ( Wf(i,j,k) ) * ( Wf(i,j,k) ) * dtime  ! WWW
    !         $acc atomic update
    !         UUU_mean3D(10,i,j,k) = UUU_mean3D(10,i,j,k) + ( Uf(i,j,k) ) * ( Vf(i,j,k) ) * ( Wf(i,j,k) ) * dtime  ! UVW

    !         $acc atomic update
    !         Udpdx_mean3D(1,i,j,k) = Udpdx_mean3D(1,i,j,k) +   Uf(i,j,k) * ((Pf(i+1,j,k)-Pf(i-1,j,k))/(dx1(i)+dx1(i-1))) * dble(2.0) *dtime  ! UdPdx
    !         $acc atomic update
    !         Udpdx_mean3D(2,i,j,k) = Udpdx_mean3D(2,i,j,k) +   Vf(i,j,k) * ((Pf(i,j+1,k)-Pf(i,j-1,k))/(dx2(j)+dx2(j-1))) * dble(2.0) *dtime  ! VdPdy
    !         $acc atomic update
    !         Udpdx_mean3D(3,i,j,k) = Udpdx_mean3D(3,i,j,k) +   Wf(i,j,k) * ((Pf(i,j,k+1)-Pf(i,j,k-1))/(dx3(k)+dx3(k-1))) * dble(2.0) *dtime  ! WdPdz
    !         $acc atomic update
    !         Udpdx_mean3D(4,i,j,k) = Udpdx_mean3D(4,i,j,k) + ( Uf(i,j,k) * ((Pf(i,j+1,k)-Pf(i,j-1,k))/(dx2(j)+dx2(j-1))) + Vf(i,j,k) * ((Pf(i+1,j,k)-Pf(i-1,j,k))/(dx1(i)+dx1(i-1))) ) *dtime  ! UdPdy + VdPdy
            
    !     enddo
    !     enddo
    !     enddo
    !     $acc end kernels

    !     ! Ghost cell update
    !     do iter=1,6
    !         Temp(:,:,:)=UU_mean3D(iter,:,:,:)

    !         call cuda_mpi_subdomain_ghostcell_update(Temp)

    !         UU_mean3D(iter,:,:,:)=Temp(:,:,:)
    !     enddo

    !     do iter=1,12
    !         Temp(:,:,:)=dUdxdUdx_mean3D(iter,:,:,:)

    !         call cuda_mpi_subdomain_ghostcell_update(Temp)

    !         dUdxdUdx_mean3D(iter,:,:,:)=Temp(:,:,:)
    !     enddo

    !     do iter=1,10
    !         Temp(:,:,:)=UUU_mean3D(iter,:,:,:)

    !         call cuda_mpi_subdomain_ghostcell_update(Temp)

    !         UUU_mean3D(iter,:,:,:)=Temp(:,:,:)
    !     enddo
        
    !     do iter=1,4
    !         Temp(:,:,:)=Udpdx_mean3D(iter,:,:,:)

    !         call cuda_mpi_subdomain_ghostcell_update(Temp)

    !         Udpdx_mean3D(iter,:,:,:)=Temp(:,:,:)
    !     enddo

    !     deallocate(Uf, Vf, Wf, Pf, Temp)

    ! end subroutine mpi_statistics_Reynolds_budget_3D

    ! subroutine mpi_statistics_Reynolds_budget_out_3D(myrank,dx1,dx2,dx3,dmx1,dmx2,dmx3) ! YMG
    !     use mpi
    !     implicit none
    !     integer :: i,j,k
    !     integer :: myrank
    !     integer :: onebyone,ierr
		
    !     double precision :: dx1(0:s1sb),dx2(0:s2sb),dx3(0:s3sb)
    !     double precision :: dmx1(0:s1sb),dmx2(0:s2sb),dmx3(0:s3sb)
    !     double precision :: RS_All(0:s2sb,4), RS_P(0:s2sb,4), RS(0:s2sb,4), RS_T(0:s2sb,4), RS_Pi(0:s2sb,4), RS_V(0:s2sb,4) ! 1:UU 2:VV 3:WW 4:UV

    !     s1msb = s1sb-1
    !     s2msb = s2sb-1
    !     s3msb = s3sb-1

    !     RS_All=0.0d0; RS_P=0.0d0; RS=0.0d0; RS_T=0.0d0; RS_Pi=0.0d0; RS_V=0.0d0 

    !     must be divided by dble(s3sb*s1sb)
    !     do k=1,s3msb
    !     do j=2,s2msb
    !     do i=1,s1msb

    !         Reynolds_Stress_Budgets UU
	! 		RS_P (j,1)   =  RS_P (j,1) - (dble(2.0)  * (  ( UU_mean3D(4,i,j,k)/avg_time_length )  *  ( (Umean(j+1)/avg_time_length) - (Umean(j-1)/avg_time_length) ) / (  dx2(j) + dx2(j-1)  ) ) ) /dble(s3msb*s1msb)

	! 		RS (j,1)   =  RS (j,1) - (dble(2.0)  * Cmu * ( dUdxdUdx_mean3D(1 ,i,j,k)/avg_time_length + dUdxdUdx_mean3D(2 ,i,j,k)/avg_time_length + dUdxdUdx_mean3D(3 ,i,j,k)/avg_time_length ) ) /dble(s3msb*s1msb)

	! 		RS_T (j,1)   =  RS_T (j,1) - ( ( ( UUU_mean3D(1 ,i+1,j,k)/avg_time_length - UUU_mean3D(1 ,i-1,j,k)/avg_time_length ) / (dx1(i)+dx1(i-1)) ) + &
    !                                      & ( ( UUU_mean3D(2 ,i,j+1,k)/avg_time_length - UUU_mean3D(2 ,i,j-1,k)/avg_time_length ) / (dx2(j)+dx2(j-1)) ) + &
    !                                      & ( ( UUU_mean3D(3 ,i,j,k+1)/avg_time_length - UUU_mean3D(3 ,i,j,k-1)/avg_time_length ) / (dx3(k)+dx3(k-1)) ) ) /dble(s3msb*s1msb)

	! 		RS_Pi(j,1)   =  RS_Pi(j,1) - Udpdx_mean3D(1,i,j,k) / avg_time_length /dble(s3msb*s1msb)

	! 		RS_V (j,1)   =  RS_V (j,1) +   Cmu  *  (  (   ( (UU_mean3D(1,i+1,j,k)/avg_time_length) - (UU_mean3D(1,i  ,j,k)/avg_time_length) ) / dx1(i)  -                &
	! 	                                             &    ( (UU_mean3D(1,i  ,j,k)/avg_time_length) - (UU_mean3D(1,i-1,j,k)/avg_time_length) ) / dx1(i-1)    ) / dmx1(i)+ &
    !                                                   (   ( (UU_mean3D(1,i,j+1,k)/avg_time_length) - (UU_mean3D(1,i,j  ,k)/avg_time_length) ) / dx2(j)  -                &
	! 	                                             &    ( (UU_mean3D(1,i,j  ,k)/avg_time_length) - (UU_mean3D(1,i,j-1,k)/avg_time_length) ) / dx2(j-1)    ) / dmx2(j)+ &
    !                                                   (   ( (UU_mean3D(1,i,j,k+1)/avg_time_length) - (UU_mean3D(1,i,j,k  )/avg_time_length) ) / dx3(k)  -                &
	! 	                                             &    ( (UU_mean3D(1,i,j,k  )/avg_time_length) - (UU_mean3D(1,i,j,k-1)/avg_time_length) ) / dx3(k-1)    ) / dmx3(k) )  /dble(s3msb*s1msb)
										   
				
    !         Reynolds_Stress_Budgets VV
    !         RS_P (j,2)   =  RS_P (j,2) - ( 0.0d0 ) /dble(s3msb*s1msb)

    !         RS (j,2)   =  RS (j,2) - (dble(2.0)  * Cmu * ( dUdxdUdx_mean3D(4 ,i,j,k)/avg_time_length + dUdxdUdx_mean3D(5 ,i,j,k)/avg_time_length + dUdxdUdx_mean3D(6 ,i,j,k)/avg_time_length ) ) /dble(s3msb*s1msb)

    !         RS_T (j,2)   =  RS_T (j,2) - ( ( ( UUU_mean3D(4 ,i+1,j,k)/avg_time_length - UUU_mean3D(4 ,i-1,j,k)/avg_time_length ) / (dx1(i)+dx1(i-1)) ) + &
    !                                      & ( ( UUU_mean3D(5 ,i,j+1,k)/avg_time_length - UUU_mean3D(5 ,i,j-1,k)/avg_time_length ) / (dx2(j)+dx2(j-1)) ) + &
    !                                      & ( ( UUU_mean3D(6 ,i,j,k+1)/avg_time_length - UUU_mean3D(6 ,i,j,k-1)/avg_time_length ) / (dx3(k)+dx3(k-1)) ) ) /dble(s3msb*s1msb)

    !         RS_Pi(j,2)   =  RS_Pi(j,2) - Udpdx_mean3D(2,i,j,k) / avg_time_length /dble(s3msb*s1msb)

    !         RS_V (j,2)   =  RS_V (j,2) +   Cmu  *  (  (   ( (UU_mean3D(2,i+1,j,k)/avg_time_length) - (UU_mean3D(2,i  ,j,k)/avg_time_length) ) / dx1(i)  -                &
    !                                                  &    ( (UU_mean3D(2,i  ,j,k)/avg_time_length) - (UU_mean3D(2,i-1,j,k)/avg_time_length) ) / dx1(i-1)    ) / dmx1(i)+ &
    !                                                   (   ( (UU_mean3D(2,i,j+1,k)/avg_time_length) - (UU_mean3D(2,i,j  ,k)/avg_time_length) ) / dx2(j)  -                &
    !                                                  &    ( (UU_mean3D(2,i,j  ,k)/avg_time_length) - (UU_mean3D(2,i,j-1,k)/avg_time_length) ) / dx2(j-1)    ) / dmx2(j)+ &
    !                                                   (   ( (UU_mean3D(2,i,j,k+1)/avg_time_length) - (UU_mean3D(2,i,j,k  )/avg_time_length) ) / dx3(k)  -                &
    !                                                  &    ( (UU_mean3D(2,i,j,k  )/avg_time_length) - (UU_mean3D(2,i,j,k-1)/avg_time_length) ) / dx3(k-1)    ) / dmx3(k) )  /dble(s3msb*s1msb)


    !         Reynolds_Stress_Budgets WW
    !         RS_P (j,3)   =  RS_P (j,3) - ( 0.0d0 ) /dble(s3msb*s1msb)

    !         RS (j,3)   =  RS (j,3) - (dble(2.0)  * Cmu * ( dUdxdUdx_mean3D(7 ,i,j,k)/avg_time_length + dUdxdUdx_mean3D(8 ,i,j,k)/avg_time_length + dUdxdUdx_mean3D(9 ,i,j,k)/avg_time_length ) ) /dble(s3msb*s1msb)

    !         RS_T (j,3)   =  RS_T (j,3) - ( ( ( UUU_mean3D(7 ,i+1,j,k)/avg_time_length - UUU_mean3D(7 ,i-1,j,k)/avg_time_length ) / (dx1(i)+dx1(i-1)) ) + &
    !                                      & ( ( UUU_mean3D(8 ,i,j+1,k)/avg_time_length - UUU_mean3D(8 ,i,j-1,k)/avg_time_length ) / (dx2(j)+dx2(j-1)) ) + &
    !                                      & ( ( UUU_mean3D(9 ,i,j,k+1)/avg_time_length - UUU_mean3D(9 ,i,j,k-1)/avg_time_length) / (dx3(k)+dx3(k-1)) ) ) /dble(s3msb*s1msb)

    !         RS_Pi(j,3)   =  RS_Pi(j,3) - Udpdx_mean3D(3,i,j,k) / avg_time_length /dble(s3msb*s1msb)

    !         RS_V (j,3)   =  RS_V (j,3) +   Cmu  *  (  (   ( (UU_mean3D(3,i+1,j,k)/avg_time_length) - (UU_mean3D(3,i  ,j,k)/avg_time_length) ) / dx1(i)  -                &
    !                                                  &    ( (UU_mean3D(3,i  ,j,k)/avg_time_length) - (UU_mean3D(3,i-1,j,k)/avg_time_length) ) / dx1(i-1)    ) / dmx1(i)+ &
    !                                                   (   ( (UU_mean3D(3,i,j+1,k)/avg_time_length) - (UU_mean3D(3,i,j  ,k)/avg_time_length) ) / dx2(j)  -                &
    !                                                  &    ( (UU_mean3D(3,i,j  ,k)/avg_time_length) - (UU_mean3D(3,i,j-1,k)/avg_time_length) ) / dx2(j-1)    ) / dmx2(j)+ &
    !                                                   (   ( (UU_mean3D(3,i,j,k+1)/avg_time_length) - (UU_mean3D(3,i,j,k  )/avg_time_length) ) / dx3(k)  -                &
    !                                                  &    ( (UU_mean3D(3,i,j,k  )/avg_time_length) - (UU_mean3D(3,i,j,k-1)/avg_time_length) ) / dx3(k-1)    ) / dmx3(k) )  /dble(s3msb*s1msb)


    !         Reynolds_Stress_Budgets UV
    !         RS_P (j,4)   =  RS_P (j,4) - ( (  ( UU_mean3D(2,i,j,k)/avg_time_length )  *  ( (Umean(j+1)/avg_time_length) - (Umean(j-1)/avg_time_length) ) / (  dx2(j) + dx2(j-1)  ) ) ) /dble(s3msb*s1msb)

    !         RS (j,4)   =  RS (j,4) - (dble(2.0)  * Cmu * ( dUdxdUdx_mean3D(10,i,j,k)/avg_time_length + dUdxdUdx_mean3D(11,i,j,k)/avg_time_length + dUdxdUdx_mean3D(12,i,j,k)/avg_time_length ) ) /dble(s3msb*s1msb)

    !         RS_T (j,4)   =  RS_T (j,4) - ( ( ( UUU_mean3D(2 ,i+1,j,k)/avg_time_length - UUU_mean3D(2 ,i-1,j,k)/avg_time_length ) / (dx1(i)+dx1(i-1)) ) + &
    !                                      & ( ( UUU_mean3D(4 ,i,j+1,k)/avg_time_length - UUU_mean3D(4 ,i,j-1,k)/avg_time_length ) / (dx2(j)+dx2(j-1)) ) + &
    !                                      & ( ( UUU_mean3D(10,i,j,k+1)/avg_time_length - UUU_mean3D(10,i,j,k-1)/avg_time_length ) / (dx3(k)+dx3(k-1)) ) ) /dble(s3msb*s1msb)

    !         RS_Pi(j,4)   =  RS_Pi(j,4) - Udpdx_mean3D(4,i,j,k) / avg_time_length /dble(s3msb*s1msb)

    !         RS_V (j,4)   =  RS_V (j,4) +   Cmu  *  (  (   ( (UU_mean3D(4,i+1,j,k)/avg_time_length) - (UU_mean3D(4,i  ,j,k)/avg_time_length) ) / dx1(i)  -                &
    !                                                  &    ( (UU_mean3D(4,i  ,j,k)/avg_time_length) - (UU_mean3D(4,i-1,j,k)/avg_time_length) ) / dx1(i-1)    ) / dmx1(i)+ &
    !                                                   (   ( (UU_mean3D(4,i,j+1,k)/avg_time_length) - (UU_mean3D(4,i,j  ,k)/avg_time_length) ) / dx2(j)  -                &
    !                                                  &    ( (UU_mean3D(4,i,j  ,k)/avg_time_length) - (UU_mean3D(4,i,j-1,k)/avg_time_length) ) / dx2(j-1)    ) / dmx2(j)+ &
    !                                                   (   ( (UU_mean3D(4,i,j,k+1)/avg_time_length) - (UU_mean3D(4,i,j,k  )/avg_time_length) ) / dx3(k)  -                &
    !                                                  &    ( (UU_mean3D(4,i,j,k  )/avg_time_length) - (UU_mean3D(4,i,j,k-1)/avg_time_length) ) / dx3(k-1)    ) / dmx3(k) )  /dble(s3msb*s1msb)

    !     enddo
    !     enddo
    !     enddo

    !     do j=1,s2msb
    !         RS_ALL(j,1)  =  RS_P(j,1)  +  RS(j,1)  +  RS_T(j,1)  +  RS_Pi(j,1)  +  RS_V(j,1)
    !         RS_ALL(j,2)  =  RS_P(j,2)  +  RS(j,2)  +  RS_T(j,2)  +  RS_Pi(j,2)  +  RS_V(j,2)	
    !         RS_ALL(j,3)  =  RS_P(j,3)  +  RS(j,3)  +  RS_T(j,3)  +  RS_Pi(j,3)  +  RS_V(j,3)
    !         RS_ALL(j,4)  =  RS_P(j,4)  +  RS(j,4)  +  RS_T(j,4)  +  RS_Pi(j,4)  +  RS_V(j,4)
	! 	enddo

	! 	Reynolds_Stress_Budgets Fileout
    !     do onebyone=0,comm_1d_x2%nprocs-1
    !         if(onebyone==comm_1d_x2%myrank) then
    !             open(myrank, file=dir_statistics//'051_Reynolds_Budget_UU'//xyzrank//'.plt')  
    !                 write(myrank,*) 'variables="y","RS_All","RS_P","RS","RS_T","RS_Pi","RS_V"'
    !                 do j=2,s2msb
    !                     write(myrank,'(7E16.8)') x2_sub(j),    RS_All(j,1)  &
	! 													    &, RS_P  (j,1)  &
	! 													    &, RS  (j,1)  &
	! 													    &, RS_T  (j,1)  &
	! 													    &, RS_Pi (j,1)  &
	! 													    &, RS_V  (j,1) 
    !                 enddo
    !             close(myrank)
    !             open(myrank, file=dir_statistics//'052_Reynolds_Budget_VV'//xyzrank//'.plt')
    !                 write(myrank,*) 'variables="y","RS_All","RS_P","RS","RS_T","RS_Pi","RS_V"'
    !                 do j=2,s2msb
    !                     write(myrank,'(7E16.8)') x2_sub(j),    RS_All(j,2)  &
	! 														&, RS_P  (j,2)  &
	! 														&, RS  (j,2)  &
	! 														&, RS_T  (j,2)  &
	! 														&, RS_Pi (j,2)  &
	! 														&, RS_V  (j,2) 
    !                 enddo	
    !             close(myrank)
    !             open(myrank, file=dir_statistics//'053_Reynolds_Budget_WW'//xyzrank//'.plt')
    !                 write(myrank,*) 'variables="y","RS_All","RS_P","RS","RS_T","RS_Pi","RS_V"'
    !                 do j=2,s2msb
    !                     write(myrank,'(7E16.8)') x2_sub(j),    RS_All(j,3)  &
	! 														&, RS_P  (j,3)  &
	! 														&, RS  (j,3)  &
	! 														&, RS_T  (j,3)  &
	! 														&, RS_Pi (j,3)  &
	! 														&, RS_V  (j,3) 
    !                 enddo
    !             close(myrank)
    !             open(myrank, file=dir_statistics//'054_Reynolds_Budget_UV'//xyzrank//'.plt')
    !                 write(myrank,*) 'variables="y","-RS_All","-RS_P","-RS","-RS_T","-RS_Pi","-RS_V"'
    !                 do j=2,s2msb
    !                     write(myrank,'(7E16.8)') x2_sub(j),       -RS_All(j,4)  &
	! 															&,-RS_P  (j,4)  &
	! 															&,-RS  (j,4)  &
	! 													    	&,-RS_T  (j,4)  &
	! 															&,-RS_Pi (j,4)  &
	! 															&,-RS_V  (j,4) 
    !                 enddo
    !             close(myrank)				
	! 		endif
    !         call MPI_Barrier(MPI_COMM_WORLD,ierr)
    !     enddo		
		

    ! end subroutine mpi_statistics_Reynolds_budget_out_3D

end module mpi_Statistics

