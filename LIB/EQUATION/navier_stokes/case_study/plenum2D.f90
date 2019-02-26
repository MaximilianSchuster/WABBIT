!=========================================================================
!# SECTION SUBROUTINES
!=========================================================================
subroutine plenum_penalization2D(Bs, g, x0, dx, phi, mask, phi_ref)
   use module_helpers
   implicit none
   !----------------------------------------------------------------------
   integer(kind=ik), intent(in) :: g              !< grid parameter, blocksize and number of ghost nodes
   integer(kind=ik), dimension(3), intent(in) :: Bs
   real(kind=rk), intent(in)    :: x0(2), dx(2)      !<coordinates of block, and block spacing
   real(kind=rk), intent(in)     :: phi(:,:,:)     !< state vector
   real(kind=rk), intent(inout)  :: phi_ref(:,:,:) !< reference values of penalized volume
   real(kind=rk), intent(inout)  :: mask(:,:,:)    !< mask function
   integer(kind=2), allocatable,save:: mask_color(:,:)!< identifyers of mask parts (plates etc)
   logical                       :: mesh_was_adapted=.true.
!-------------------------------------------------------------------------
    if (.not. allocated(mask_color))  allocate(mask_color(1:Bs(1)+2*g, 1:Bs(2)+2*g))
    !!> todo implement function check_if_mesh_adapted (true/false) in adapt mesh
    if ( mesh_was_adapted .eqv. .true. ) then
      ! dont switch the order of draw_plenum3D and draw_sponge3D,
      ! because mask and color are reset in the draw_plenum
      call draw_plenum2D(x0, dx, Bs, g, mask, mask_color)
      call draw_plenum_sponge2D(x0, dx, Bs, g, mask, mask_color)
    end if

   call compute_plenum_penal2D(mask_color,mask,phi, x0, dx, Bs, g ,phi_ref)

end subroutine plenum_penalization2D

!========================================================================

subroutine draw_plenum2D(x0, dx, Bs, g, mask, mask_color)
    implicit none

!########
!#DECLARE VARIABLES
    ! -----------------------------------------------------------------
    integer(kind=ik), intent(in)             :: g          !< grid parameter
    integer(kind=ik), dimension(3), intent(in) :: Bs
    real(kind=rk), intent(in)                :: x0(2), dx(2)   !< coordinates of block and block spacinf
    real(kind=rk), intent(inout)             :: mask(:,:,:)    !< mask function
    integer(kind=2), intent(inout), optional :: mask_color(:,:)!< identifyers of mask parts (plates etc)

    ! -----------------------------------------------------------------
    real(kind=rk)     :: x, y, r, h
    real(kind=rk)     :: chi,Ly_eff,Ly_div
    integer(kind=ik)  :: ix, iy,n,blade_no ! loop variables
    integer(kind=ik) ::NoB
    real(kind=rk)     ::y_p_stag
  ! -----------------------------------------------------------------
NoB = plenum%n_blades
y_p_stag = 0.2_rk
h  = 1.5_rk*max(dx(1), dx(2))
Ly_eff = domain_size(2)-2.0_rk*plenum%wall_thickness
    ! smooth width in x and y direction
    do iy=g+1, Bs(2)+g
       y = dble(iy-(g+1)) * dx(2) + x0(2)
       r = abs(y-domain_size(2)*0.5_rk)
       do ix=g+1, Bs(1)+g
            x = dble(ix-(g+1)) * dx(1) + x0(1)

            !=============================
!     /\     reset the mask function!
!    /  \    caution! Reseting should only be done once
!   /    \   for optimal performance
!  / stop \
! +--------+
            mask_color(ix,iy)=0
            mask(ix,iy,:)=0.0_rk
            !=============================

            ! Walls
            ! -----
            chi = draw_plenum_walls(x,r,plenum,h)            
            if (chi>0.0_rk) then                       ! default values on the walls
              mask_color(ix,iy)  = color_walls
              mask(ix,iy,1:4)    = mask(ix,iy,1:4) + chi
            end if


        
 
 
 if (NoB>0) then
       do blade_no=1,NoB
                  Ly_div = (domain_size(2)-2.0_rk*plenum%wall_thickness)/(NoB+1.0_rk)
                  y_p_stag = blade_no*Ly_div+plenum%wall_thickness
                  chi = draw_stat_blades(x,y,r,plenum,y_p_stag)
                     if (chi>0.0_rk) then
                        mask_color(ix,iy)  = color_walls
                        mask(ix,iy,1:4)    = mask(ix,iy,1:4) + chi
                     endif
          
           
           end do
           end if
       end do
    end do

end subroutine draw_plenum2D

!=======================================================================

subroutine draw_plenum_sponge2D(x0, dx, Bs, g, mask, mask_color)
    implicit none
!############
!#DECLARE VARIABLES
!###########
    ! -----------------------------------------------------------------
    integer(kind=ik), intent(in)             :: g          !< grid parameter
    integer(kind=ik), dimension(3), intent(in) :: Bs
    real(kind=rk)                           :: sp_thickness
    real(kind=rk), intent(in)                :: x0(2), dx(2)   !< coordinates of block and block spacinf
    real(kind=rk), intent(inout)             :: mask(:,:,:)    !< mask function
    integer(kind=2), intent(inout), optional :: mask_color(:,:)!< identifyers of mask parts (plates etc)
    ! -----------------------------------------------------------------
    real(kind=rk)     :: x, y, r, h
    real(kind=rk)     :: chi
    integer(kind=ik)  :: ix, iy,n ! loop variables

  ! -----plenum%length_pip+plenum%wall_thickness------------------------------------------------------------
!INITIALISE PARAMETERS
!sp_thickness = plenum%sponge_thickness

sp_thickness=0.05_rk*domain_size(1)
    ! parameter for smoothing function (width)
    h  = 1.5_rk*max(dx(1), dx(2))

    ! smooth width in x and y direction
    do iy=g+1, Bs(2)+g
       y = dble(iy-(g+1)) * dx(2) + x0(2)
       r = abs(y-domain_size(2)*0.5_rk)
        do ix=g+1, Bs(1)+g
            x = dble(ix-(g+1)) * dx(1) + x0(1)

            chi = draw_plenum_rad_sponge(x,r,plenum,h)
            if(chi>0.0_rk) then
                mask_color(ix,iy) = color_walls
                mask(ix,iy,1) = mask(ix,iy,1)+chi
                mask(ix,iy,2) = mask(ix,iy,2)+chi
               mask(ix,iy,3:4) = mask(ix,iy,3:4)+chi
            endif
        end do
    end do
end subroutine draw_plenum_sponge2D

!=======================================================================

subroutine compute_plenum_penal2D(mask_color,mask,phi,x0,dx,Bs,g,phi_ref)
!##################
!#DECLARE VARIABLES
    implicit none
    ! -----------------------------------------------------------------
    integer(kind=ik), intent(in)  :: g          !< grid parameter
    integer(kind=ik), dimension(3), intent(in) :: Bs
    real(kind=rk), intent(in)     :: x0(2), dx(2)   !< coordinates of block and block spacinf
    integer(kind=2), intent(inout):: mask_color(:,:)!< identifyers of mask parts (plates etc)
    real(kind=rk), intent(in)     :: phi(:,:,:)     !< state vector
    real(kind=rk), intent(inout)  :: mask(:,:,:)     !< mask
    real(kind=rk), intent(inout)  :: phi_ref(:,:,:)  !< plenum penalization term
    ! -----------------------------------------------------------------
!##################
!# ABOUT PHI_REF:1-rho 2- u, 3-v, 4-p
!#initialise parameters
real(kind=rk)    :: x,y,r,h,velocity
real(kind=rk)   :: rho,chi, v_ref,dq,u,v,p,C_inv
real(kind=rk) :: u_inlet, v_inlet, rho_inlet, p_inlet
real(kind=rk) :: rho_outlet,p_outlet
real(kind=rk) :: rho_sponge, p_sponge,v_sponge,temperature
integer(kind=ik) :: ix,iy,n
real(kind=rk)    :: length_pip,diameter_pip,radius_pip
!real(kind=rk)    ::sponge_thickness
! INITIALIZE PARAMETERS
! -------------------------------------
phi_ref=0.0_rk
u_inlet=plenum%inlet_velocity(1)
v_inlet=plenum%inlet_velocity(2)
rho_inlet=plenum%inlet_density
p_inlet = plenum%inlet_pressure
!--------------------------------------
rho_outlet=plenum%outlet_density
p_outlet  =plenum%outlet_pressure
temperature = plenum%temperature
!--------------------------------------
rho_sponge = plenum%sponge_density
p_sponge = plenum%sponge_pressure
v_sponge = plenum%sponge_velocity(2)
length_pip = plenum%length_pip
diameter_pip = plenum%diameter_pip
radius_pip = 0.5_rk*diameter_pip

!--------------------------------------
 !set smoothing parameter
h  = 1.5_rk*max(dx(1), dx(2))
    do iy=1, Bs(2)+2*g
       y = dble(iy-(g+1)) * dx(2) + x0(2)
       r = abs(y-domain_size(2)*0.5_rk)
       do ix=1, Bs(1)+2*g
            x = dble(ix-(g+1)) * dx(1) + x0(1)

            rho = phi(ix,iy,rhoF)
            u   = phi(ix,iy,UxF)
            v   = phi(ix,iy,UyF)
            p   = phi(ix,iy,pF)
            C_inv=C_eta_inv

     !what happens in solid obstacles (walls) no slip wall is penalized by u=v=0, rho == arbitrary value and p == computed via ideal gas law
            if (mask_color(ix,iy)==color_walls) then
               Phi_ref(ix,iy,2) = 0.0_rk ! no velocity in x
               Phi_ref(ix,iy,3) = 0.0_rk ! no velocity in y
               Phi_ref(ix,iy, 4) = rho * Rs * temperature! set pressure with ideal gas law
               Phi_ref(ix,iy,1) = rho_sponge !
              C_inv = C_eta_inv
            end if

    mask(ix,iy,:) = C_inv*mask(ix,iy,:)
   end do
    end do

end subroutine compute_plenum_penal2D

!=======================================================================

!###SECTION FUNCTIONS

!=======================================================================

function draw_plenum_walls(x, r, plenum, h)
implicit none
!-----------------------------------------------------------------------!
  real(kind=rk),    intent(in)          :: x, r, h
  type(type_plenum),intent(in)          ::plenum

  real(kind=rk)                         ::mask, draw_plenum_walls
  real(kind=rk)                         ::r0_pip,r0_ple,width,l_pip
  real(kind=rk)                         ::sp1, sp2,sp3
  real(kind=rk)                         ::wall_thickness 
  character(len=90)                     :: geo_case
!-----------------------------------------------------------------------!
  mask=0.0_rk
  wall_thickness=plenum%wall_thickness
  r0_ple=0.5_rk*plenum%diameter_ple-plenum%wall_thickness
  width=R_domain
  l_pip = plenum%length_pip
  geo_case = plenum%name
  sp1 = 0.15_rk/0.8_rk*l_pip+plenum%wall_thickness
  sp2 = 0.1603_rk/0.8_rk*l_pip+plenum%wall_thickness
  sp3 = 0.317_rk/0.8_rk*l_pip+plenum%wall_thickness!*2.0_rk+wall_thickness
  !-------CONSTANT RADIUS PIPE CASE----------!
  
  
  
!-----------------------------------------------------------!
!-----------------____________________----------------------!
!----------------|CONSTANT RADIUS PIPE|---------------------!
!----------------|____________________|---------------------!
  select case(geo_case)
  case('const')
  if (x<= plenum%length_pip) then
    ! mask for r>R_domain-wall_thickness   (outer wall)
         !mask=mask+smoothstep(R_domain-funnel%wall_thickness-r,h)
        r0_pip=0.5_rk*plenum%diameter_pip
       mask=hard_bump(r,r0_pip,width)
  elseif (x >= plenum%length_pip) then

       mask=soft_bump(r,r0_ple,width,h)

  endif


!  draw_plenum_walls=mask


!-----------------------------------------------------------!
!-----------------____________------------------------------!
!----------------|SFC_PDC_CASE|-----------------------------!
!----------------|____________|-----------------------------!
  case('sfb_pdc')
  if (x<=sp1)then 
     r0_pip = 0.0403_rk
     mask=hard_bump(r,r0_pip,width)
  elseif (x>sp1 .and. x<=sp2) then
     r0_pip = 0.0403_rk - (x-sp1)/(sp1-sp2)*(0.02_rk-0.0403_rk)
     mask=soft_bump(r,r0_pip,width,h)
  elseif(x>sp2 .and. x<=sp3) then
     r0_pip = 0.02_rk + (x-sp2)/(sp2-sp3)*(0.02_rk-0.0403_rk)
     mask=soft_bump(r,r0_pip,width,h)
  elseif(x>sp3 .and. x <l_pip) then
     r0_pip = 0.0403_rk
     mask=soft_bump(r,r0_pip,width,h)
  end if
  
  case default
  call abort(8546533,"ERROR: No geometry for the plenum. "//plenum%name)
  
  
  end select
  if (x< plenum%wall_thickness) then
  mask = mask+1.0_rk
  endif
  
  !---------!
  
  
    ! wall in east
  if (x< plenum%wall_thickness) then
  mask = mask+1.0_rk
  endif
  ! wall in west
  if (x>domain_size(1)-plenum%wall_thickness) then
  mask=mask+1.0_rk
  endif
  
  if (mask>1.0_rk) then
  mask=1.0_rk
  endif
  
  draw_plenum_walls=mask
  
end function

!=======================================================================

           


function draw_stat_blades(x,y,r,plenum,y_p_stag)
implicit none
!-----------------------------------------------------------------------------------!
real(kind=rk), intent(in)               :: x,y,r
type(type_plenum),intent(in)            :: plenum
real(kind=rk)                           :: alpha,AoA_rad, length_blade,thickness_blade
real(kind=rk)                           :: eta, nu
real(kind=rk)                           :: blade_gap
real(kind=rk)                           :: x_m,y_m, nu_m,eta_m
real(kind=rk)                           ::x_p_stag, y_p_stag, nu_p_stag, eta_p_stag
real(kind=rk)                           ::draw_stat_blades,mask
!real(kind=rk),dimension(2)::A1,A2,A3,A4 
!====================================================================================

length_blade = plenum%L_blade! how long is the turbine blade (chord length)
thickness_blade = plenum%t_blade ! how thick is the turbine blade
blade_gap = plenum%gap_blade
x_p_stag = plenum%length_pip+blade_gap! where is the x location of the blade nose
alpha = plenum%AoA
AoA_rad = alpha*3.14159265359_rk/180.0_rk




! COORDINATE SYSTEM TRANSFORMATION
nu = cos(-AoA_rad)*x + sin(-AoA_rad)*y
eta = -sin(-AoA_rad)*x + cos(-AoA_rad)*y 
! STAGNATION POINT COORDINATES
nu_p_stag = cos(-AoA_rad)*x_p_stag + sin(-AoA_rad)*y_p_stag
eta_p_stag = -sin(-AoA_rad)*x_p_stag + cos(-AoA_rad)*y_p_stag
! center point coordinates
x_m = plenum%length_pip+ 0.5_rk*(domain_size(1)-plenum%length_pip)
y_m = y_p_stag
nu_m = cos(-AoA_rad)*x_m + sin(-AoA_rad)*y_m
eta_m =  -sin(-AoA_rad)*x_m+ cos(-AoA_rad)*y_m


!FIND CORNERS OF RECTANGLE in nu,eta coordinate system
!          A4---------------A3
!           |               |
!           |       x       |
!           |     (x_m,y_m) |
!          A1---------------A2

!A1(1)= nu_p_stag
!A1(2)= eta_p_stag-T_blade/2.0_rk
!A2(1)= nu_p_stag+L_blade
!A2(2)= eta_p_stag -T_blade/2.0_rk
!A3(1)= nu_p_stag+L_blade
!A3(2)= eta_p_stag + T_blade/2.0_rk
!A4(1)= nu_p_stag
!A4(2)= eta_p_stag + T_blade/2.0_rk
  if ( nu >(nu_m-length_blade/2.0_rk) .and. &
       nu < (nu_m+length_blade/2.0_rk) .and. &
       eta > (eta_m-thickness_blade/2.0_rk) .and. &
       eta < (eta_m+thickness_blade/2.0_rk)) then
       !write(*,*)" We re inside the stator blade"

       mask=1.0_rk
       
       if((y<plenum%wall_thickness) .or. (y>domain_size(2)-plenum%wall_thickness))  then
          call abort(90790,"Error: The blades do not fit the plenum domain. Reduce number of blades or blade length")
       endif
  endif

draw_stat_blades=mask

end function draw_stat_blades


function draw_plenum_rad_sponge(x, r, plenum, h)
implicit none
!------------------------------------------------------------------!
  real(kind=rk),    intent(in)          :: x, r, h
  type(type_plenum),intent(in)          ::plenum

  real(kind=rk)                         ::draw_plenum_rad_sponge,mask
  real(kind=rk)                         ::r0,width
!-------------------------------------------------------------------!
  mask=0.0_rk
  r0=0.5_rk*plenum%diameter_ple-plenum%wall_thickness
  width=plenum%wall_thickness
if (x>plenum%length_pip.and. &
   x< domain_size(1)-plenum%wall_thickness) then
   mask=soft_bump(r,r0,width,h)
endif
draw_plenum_rad_sponge=mask
end function



function draw_parabolic_sponge(x,y,r,plenum)
implicit none
!------------------------------------------------------------------------!
! DECLARE VARIABLES!
real(kind=rk), intent(in)              ::x,y,r
type(type_plenum), intent(in)          ::plenum

real(kind=rk)                          ::draw_parabolic_sponge, x_min,x_max
real(kind=rk)                          ::mask, mask1, mask2

!------------------------------------------------------------------------!


mask1 = 0.0_rk
mask2 = 0.0_rk
! quadratic power law determines strength of the sponge in radial plenum direction
    if (x>=plenum%length_pip .and. abs(r)>domain_size(2)/2.0_rk-plenum%wall_thickness) then
        x_min=domain_size(2)/2.0_rk-plenum%wall_thickness
        x_max=domain_size(2)/2.0_rk
        mask=((abs(r)-x_min)/(x_max-x_min))**2
! quadratic power law determines strength of the sponge in RIGHT BOUNDARY DOMAIN (OUTLET)
   elseif(x>domain_size(1)-plenum%wall_thickness .and. abs(r)<domain_size(2)/2.0_rk-plenum%wall_thickness) then 
        x_min=domain_size(1)-plenum%wall_thickness
        x_max=domain_size(1)
        mask=((x-x_min)/(x_max-x_min))**2
! multiplying power law for "bottom right corner" and "top right corner" of the domain (--> x^4 Law?)
    elseif(x>domain_size(1)-plenum%wall_thickness .and. abs(r)>=domain_size(2)/2.0_rk-plenum%wall_thickness) then
        x_min=domain_size(2)/2.0_rk-plenum%wall_thickness
        x_max=domain_size(2)/2.0_rk
        mask1=((abs(r)-x_min)/(x_max-x_min))**2

        x_min=domain_size(1)-plenum%wall_thickness
        x_max=domain_size(1)
        mask2=((x-x_min)/(x_max-x_min))**2
        mask = mask1*mask2
    endif 

    draw_parabolic_sponge=mask
end function




