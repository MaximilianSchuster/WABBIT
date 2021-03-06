
!-----------------------------------------------------------------
!! Module of 2D/3D gas turbine plenum
!> \version 18/10/18
!> \author M. Schuster
!-----------------------------------------------------------------

module module_plenum

  !-------------------------------------------------------
  ! modules
  use module_navier_stokes_params
  use module_precision
  use module_ns_penalization
  use module_ini_files_parser_mpi, only : read_param_mpi
  use mpi
  !--------------------------------------------------------
  ! variables

  implicit none

  ! I usually find it helpful to use the private keyword by itself initially, which specifies
  ! that everything within the module is private unless explicitly marked public.
  PRIVATE

  !**********************************************************************************************
  ! These are the important routines that are visible to WABBIT:
  !**********************************************************************************************
  PUBLIC :: read_params_plenum,mean_quantity_plenum,draw_plenum, &
            set_inicond_plenum,plenum_penalization2D!,plenum_penalization3D,integrate_over_pump_area,
  !**********************************************************************************************

!  real(kind=rk),    allocatable,     save        :: mask(:,:,:)
  character(len=80),save    :: mask_geometry!273.15_rk
  logical      ,save        :: smooth_mask, use_sponge
  !real(kind=rk),save        :: C_eta_inv, C_sp_inv, L_sponge
  real(kind=rk),save        :: domain_size(3)=0.0_rk
  ! radius of domain (Ly/2)
  real(kind=rk),save        :: R_domain
  real(kind=rk),save        :: Rs,gamma_


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
! identifyers of the different parts of the plenum
! they are used in the array mask_color
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++
  integer(kind=2),parameter :: color_inlet            =6
  integer(kind=2),parameter :: color_blades           =5
  integer(kind=2),parameter :: color_rad_sponge       =4
  integer(kind=2),parameter :: color_walls            =3
  integer(kind=2),parameter :: color_west_sponge      =2
  integer(kind=2),parameter :: color_east_wall        =1
!++++++++++++++++++++++++++++++++++++++++++++++++


  type :: type_plenum
      character(len=80)   :: name,blade_position  
      real(kind=rk)       :: diameter_ple          ! Plenum diameter
      real(kind=rk)       :: length_ple 
      real(kind=rk), dimension(9)   :: lower_coefficients,upper_coefficients
      real(kind=rk)       :: sponge_density
      real(kind=rk)       :: sponge_pressure
      real(kind=rk)       :: temperature 
      real(kind=rk)       :: wall_thickness       
      real(kind=rk)       :: diameter_pip                                    
      real(kind=rk)       :: length_pip, boxwidth,boxheight
      real(kind=rk)       ::inlet_velocity,inlet_pressure, inlet_density
      
      !BLADE VARIABLES
      real(kind=rk)       :: L_blade
      real(kind=rk)       :: AoA
      integer(kind=ik)    :: n_blades
      integer(kind=ik)    :: set_inibox
      real(kind=rk)       :: gap_blade
      real(kind=rk)       :: t_blade,inibox_length

  end type type_plenum


  !------------------------------------------------
  type(type_plenum)   , save :: plenum
  !------------------------------------------------

contains

  include 'plenum2D.f90'
  !include 'plenum3D.f90'



    !> \brief reads parameters for mask function from file
    subroutine read_params_plenum(params,FILE)

      implicit none
      real(kind=rk)::L_blade, AoA, gap_blade,t_blade
      character(len=90) :: name
      integer(kind=ik)::n_blades
      ! character(len=*), intent(in) :: filename
      type(inifile) , intent(inout) :: FILE
      !> params structure of navier stokes
      type(type_params_ns),intent(inout)  :: params
      ! READ IN geometry
      ! ----------------
      call read_param_mpi(FILE, 'plenum', 'name', plenum%name,'')
      call read_param_mpi(FILE, 'plenum', 'blade_position', plenum%blade_position,'')
      call read_param_mpi(FILE, 'plenum', 'diameter_ple'  , plenum%diameter_ple, 0.1_rk)
      call read_param_mpi(FILE, 'plenum', 'length_ple'  , plenum%length_ple, 0.156_rk)
      call read_param_mpi(FILE, 'plenum', 'sponge_density', plenum%sponge_density,1.22_rk)
      call read_param_mpi(FILE, 'plenum', 'sponge_pressure', plenum%sponge_pressure, 101000.00_rk) 
      call read_param_mpi(FILE, 'plenum', 'temperature'  , plenum%temperature, 300.0_rk)
      call read_param_mpi(FILE, 'plenum', 'wall_thickness',plenum%wall_thickness, 0.00625_rk)
      call read_param_mpi(FILE, 'plenum', 'diameter_pip',plenum%diameter_pip, 0.1_rk)
      call read_param_mpi(FILE, 'plenum', 'length_pip',plenum%length_pip, 0.8_rk)

      
      
      
      call read_param_mpi(FILE, 'plenum', 'inlet_density', plenum%inlet_density,1.22_rk)
      call read_param_mpi(FILE, 'plenum', 'inlet_pressure', plenum%inlet_pressure, 101000.00_rk) 
      call read_param_mpi(FILE, 'plenum', 'inlet_velocity', plenum%inlet_velocity, 0.00_rk) 
      ! BLADE VARIABLES AND PARAMETERS
      call read_param_mpi(FILE, 'plenum', 'L_blade'  , plenum%L_blade, 0.1_rk)
      call read_param_mpi(FILE, 'plenum', 'AoA'  , plenum%AoA, 0.0_rk)
      call read_param_mpi(FILE, 'plenum', 'n_blades'  , plenum%n_blades, 1)
      call read_param_mpi(FILE, 'plenum', 'set_inibox'  , plenum%set_inibox, 0)
      call read_param_mpi(FILE, 'plenum', 'inibox_length'  , plenum%inibox_length, 0.1_rk)
      call read_param_mpi(FILE, 'plenum', 'gap_blade', plenum%gap_blade, 0.05_rk)
      call read_param_mpi(FILE, 'plenum', 't_blade', plenum%t_blade, 0.02_rk)
      call read_param_mpi(FILE, 'plenum', 'boxheight'  , plenum%boxheight, 0.05_rk)
      call read_param_mpi(FILE, 'plenum', 'boxwidth'  , plenum%boxwidth, 0.1_rk)
      ! POLYNOMIAL COEFFICEINTS)  
      call read_param_mpi(FILE, 'plenum', 'lower_coefficients'  , plenum%lower_coefficients(1:9))               
      call read_param_mpi(FILE, 'plenum', 'upper_coefficients'  , plenum%upper_coefficients(1:9))  
      ! these parameters are global in plenum module!
      Rs         =params%Rs
      gamma_     =params%gamma_
      domain_size=params%domain_size
      R_domain   =params%domain_size(2)*0.5_rk
      C_sp_inv   =1.0_rk/params%C_sp
      C_eta_inv   =1.0_rk/params%C_eta

end subroutine read_params_plenum


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!> Allocate and compute mask for 2D/3D plenum. The different parts of the mask can be
!> colored if the boolean mask_is_colored is true. If the boolean is false then mask returns
!> only the mask of the solid objects (like walls and plates)
subroutine  draw_plenum(x0, dx, Bs, g, mask, mask_is_colored)
  implicit none
  ! -----------------------------------------------------------------
  integer(kind=ik), intent(in)  :: g        !< grid parameter
  integer(kind=ik), dimension(3), intent(in) :: Bs
  real(kind=rk), intent(in)     :: x0(3), dx(3) !< coordinates of block and block spacinf
  real(kind=rk), intent(inout)  :: mask(:,:,:)    !< mask function
  logical, optional, intent(in) :: mask_is_colored
  integer(kind=2),allocatable   :: mask_color(:,:,:)!< identifyers of mask parts (plates etc)
  logical, save :: is_colored =.false.
  real(kind=rk), allocatable  :: mask_tmp(:,:,:,:)    !< mask function for the statevector
  ! -----------------------------------------------------------------
  if (size(mask,1) /= Bs(1)+2*g) call abort(127109,"wrong array size!")
  ! if variable is present the default (false) is overwritten by the input
  if( present(mask_is_colored)) is_colored=mask_is_colored
  ! allocate and compute mask and colored mask
!  if (params_ns%dim==3) then
!    if (.not. allocated(mask_color))  allocate(mask_color(1:Bs+2*g, 1:Bs+2*g, 1:Bs+2*g))
!    if (.not. allocated(mask_tmp))        allocate(mask_tmp(1:Bs+2*g, 1:Bs+2*g, 1:Bs+2*g,5))
!    mask_tmp    = 0.0_rk
!    mask_color  = 0
!    call  draw_plenum3D(x0, dx, Bs, g, mask_tmp, mask_color)
!    call  draw_sponge3D(x0,dx,Bs,g,mask_tmp,mask_color)
!  else
    if (.not. allocated(mask_color))  allocate(mask_color(1:Bs(1)+2*g, 1:Bs(2)+2*g, 1))
    if (.not. allocated(mask_tmp))        allocate(mask_tmp(1:Bs(1)+2*g, 1:Bs(2)+2*g, 1,4))
    mask_tmp    = 0.0_rk
    mask_color  = 0
    call  draw_plenum2D(x0, dx, Bs, g, mask_tmp(:,:,1,:), mask_color(:,:,1))
    call  draw_plenum_sponge2D(x0,dx,Bs,g,mask_tmp(:,:,1,:), mask_color(:,:,1))

!  endif

  ! mask coloring is optional, which is mainly used for plotting the different parts
  ! of the plenum in paraview
  if (.false.) then
    mask(:,:,:)= real(mask_color(:,:,:),kind=rk)
  else
    ! if the mask is not colored we use the mask of the solid obstacles
    mask(:,:,:) = mask_tmp(:,:,:,UxF)
  endif

end subroutine draw_plenum
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


  !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !> \brief Set the initial condition of a specific case
  !> \details
   subroutine set_inicond_plenum(x0, dx, Bs, g, u)
       implicit none
       ! -----------------------------------------------------------------
       integer(kind=ik), intent(in)  :: g        !< grid parameter
       integer(kind=ik), dimension(3), intent(in) :: Bs
       real(kind=rk), intent(in)     :: x0(3), dx(3) !< coordinates of block and block spacinf
       real(kind=rk), intent(inout)  :: u(:,:,:,:)    !< Statevector for t=0
       ! -----------------------------------------------------------------
       real(kind=rk),allocatable:: mask(:,:,:)
       real(kind=rk)            :: y_rel,p_init, rho_init,u_init(3),T_init,b
       real(kind=rk)            ::  r_pip
       integer(kind=ik)         :: iy,ix
       real(kind=rk)     :: x, y, r, h

       p_init    =params_ns%initial_pressure
       rho_init  =params_ns%initial_density
       u_init    =params_ns%initial_velocity
       T_init    =params_ns%initial_temp
       ! FU HARDCODE
       if (x< 0.01 .and. abs(r)<0.015) then
       p_init = 330000.0_rk
       rho_init=2.71_rk
       u_init = 323.0_rk
       endif
       ! FU HARDCODE
       allocate(mask(size(u,1), size(u,2), size(u,3)))
      ! set velocity field u(x)=1 for x in mask
      ! u(x)=(1-mask(x))*u0 to make sure that flow is zero at mask values
      call draw_plenum(x0, dx, Bs, g, mask)
      u( :, :, :, pF) = p_init
      u( :, :, :, rhoF) = sqrt(rho_init)
      u( :, :, :, UxF) = ( 1 - mask ) * u_init(1)*sqrt(rho_init) !flow in x
      u( :, :, :, UyF) = (1-mask)*u_init(2)*sqrt(rho_init) !flow in y
      if (params_ns%dim==3) then
        u( :, :, :, UzF) = (1-mask)*u_init(2)*sqrt(rho_init) !flow in z
      endif


    end subroutine set_inicond_plenum
    !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



 subroutine integrate_over_pump_area(u,g,Bs,x0,dx,integral,area)
      implicit none
      !---------------------------------------------------------------
      integer(kind=ik), intent(in):: g            !< grid parameter (g ghostnotes,Bs Bulk)
      integer(kind=ik), dimension(3), intent(in) :: Bs      
      real(kind=rk), intent(in)   :: u(:,:,:,:)       !< statevector in PURE VARIABLES \f$ (rho,u,v,w,p) \f$
      real(kind=rk),  intent(in)  :: x0(3), dx(3)     !< spacing and origin of block
      real(kind=rk),intent(out)   :: integral(5), area!< mean values
      !---------------------------------------------------------------!


     if ( params_ns%dim==2 ) then
        call integrate_over_pump_area2D(u(:,:,1,:),g,Bs,x0(1:2),dx(1:2),integral(1:4),area)
!      else
!        call integrate_over_pump_area3D(u,g,Bs,x0,dx,integral,area)
      end if

  end subroutine integrate_over_pump_area





 subroutine mean_quantity_plenum(integral,area)
      !> area of taking the mean
     real(kind=rk),intent(in)    :: area
      !> integral over the area
      real(kind=rk),intent(inout) :: integral(1:)

      ! temporary values
      real(kind=rk),allocatable,save :: tmp(:)
      real(kind=rk)                  :: A
      integer(kind=ik)               :: mpierr,Nq


      Nq = size(integral,1)
      if ( .not. allocated(tmp) ) allocate(tmp(Nq))

      tmp=integral

    !   integrate over all procs
      call MPI_ALLREDUCE(tmp  ,integral, Nq , MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)
      call MPI_ALLREDUCE(area ,A       , 1  , MPI_DOUBLE_PRECISION, MPI_SUM, WABBIT_COMM, mpierr)

      if ( .not. abs(A) > 0) then
        call abort(24636,"Error [plenum.f90]: only chuck norris can devide by zero!!")
      endif

      !devide by the area of the region
      integral = integral / A
   !   plenum%pump_density = integral(rhoF)
  !    plenum%pump_pressure = integral(pF)
  end subroutine mean_quantity_plenum

end module module_plenum
