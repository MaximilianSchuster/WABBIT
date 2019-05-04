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
    character(len=90)                     :: blade_position
    real(kind=rk)     :: x, y, r, h
    real(kind=rk)     :: chi,Ly_eff,Ly_div
    integer(kind=ik)  :: ix, iy,n,blade_no ! loop variables
    integer(kind=ik) ::NoB
    real(kind=rk)     ::y_p_stag,gap_blade,boxwidth,boxheight
  ! -----------------------------------------------------------------
NoB = plenum%n_blades
gap_blade = plenum%gap_blade
blade_position = plenum%blade_position
y_p_stag = 0.2_rk
h  = 2.5_rk*max(dx(1), dx(2))
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
            chi = draw_plenum_walls(x,y,r,plenum,h)            
            if (chi>0.0_rk) then                       ! default values on the walls
              mask_color(ix,iy)  = color_walls
              mask(ix,iy,UxF)    = mask(ix,iy,UxF) + chi
              mask(ix,iy,UyF)    = mask(ix,iy,UyF) + chi
            !  mask(ix,iy,pF)     = mask(ix,iy,pF)  + chi
            ! mask(ix,iy,rhoF)   = mask(ix,iy,rhoF)  + chi
            end if
          
       if (NoB>0) then
       ! do blade_no=1,NoB                        
       ! ONLY USE INSIDE A DEFINED BOX
       boxwidth=plenum%boxwidth
       boxheight=plenum%boxheight
       
       
       
         ! SET POSITION OF THE FIRST
         select case(blade_position)
         case('center')
       y_p_stag = 0.5_rk*domain_size(2)      
       if (x>=plenum%length_pip+gap_blade .and. x<=plenum%length_pip+gap_blade+boxwidth) then
           if (y>=y_p_stag .and. y<=y_p_stag+boxheight) then
                  chi = draw_poly_blades(x,y,r,plenum,y_p_stag)
                     if (chi>0.0_rk) then
                        mask_color(ix,iy)  = color_blades
                        mask(ix,iy,UxF)    = mask(ix,iy,UxF) + chi
                        mask(ix,iy,UyF)    = mask(ix,iy,UyF) + chi
                        mask(ix,iy,pF)     = mask(ix,iy,pF)  + chi
       !                 mask(ix,iy,rhoF) = mask(ix,iy,rhoF) + chi
                     endif
           endif
        endif
        do blade_no=0,NoB/2-1
 
                  y_p_stag = 0.5_rk*domain_size(2)+0.0763_rk*(blade_no-1)
               if (x>=plenum%length_pip+gap_blade .and. x<=plenum%length_pip+gap_blade+boxwidth) then
                  if (y>=y_p_stag .and. y<=y_p_stag+boxheight) then
                  chi = draw_poly_blades(x,y,r,plenum,y_p_stag)
                     if (chi>0.0_rk) then
                        mask_color(ix,iy)  = color_blades
                        mask(ix,iy,UxF)    = mask(ix,iy,UxF) + chi
                        mask(ix,iy,UyF)    = mask(ix,iy,UyF) + chi
                        mask(ix,iy,pF)     = mask(ix,iy,pF)  + chi
       !                 mask(ix,iy,rhoF) = mask(ix,iy,rhoF) + chi
                     endif
                     endif
           endif
          end do 
          
                do blade_no=0,NoB/2-1
 
                  y_p_stag = 0.5_rk*domain_size(2)-0.0763_rk*(blade_no-1)
               if (x>=plenum%length_pip+gap_blade .and. x<=plenum%length_pip+gap_blade+boxwidth) then
                  if (y>=y_p_stag .and. y<=y_p_stag+boxheight) then
                  chi = draw_poly_blades(x,y,r,plenum,y_p_stag)
                     if (chi>0.0_rk) then
                        mask_color(ix,iy)  = color_blades
                        mask(ix,iy,UxF)    = mask(ix,iy,UxF) + chi
                        mask(ix,iy,UyF)    = mask(ix,iy,UyF) + chi
                        mask(ix,iy,pF)     = mask(ix,iy,pF)  + chi
       !                 mask(ix,iy,rhoF) = mask(ix,iy,rhoF) + chi
                     endif
                     endif
           endif
          end do 
         
         case('channel')
         y_p_stag = 0.5_rk*domain_size(2)+0.043_rk    
            if (x>=plenum%length_pip+gap_blade .and. x<=plenum%length_pip+gap_blade+boxwidth) then
               if (y>=y_p_stag .and. y<=y_p_stag+boxheight) then
                  chi = draw_poly_blades(x,y,r,plenum,y_p_stag)
                     if (chi>0.0_rk) then
                        mask_color(ix,iy)  = color_blades
                        mask(ix,iy,UxF)    = mask(ix,iy,UxF) + chi
                        mask(ix,iy,UyF)    = mask(ix,iy,UyF) + chi
                        mask(ix,iy,pF)     = mask(ix,iy,pF)  + chi
       !                 mask(ix,iy,rhoF) = mask(ix,iy,rhoF) + chi
                     endif
               endif
            endif
            
            y_p_stag = 0.5_rk*domain_size(2)+0.043_rk-0.0763_rk      
            if (x>=plenum%length_pip+gap_blade .and. x<=plenum%length_pip+gap_blade+boxwidth) then
               if (y>=y_p_stag .and. y<=y_p_stag+boxheight) then
                  chi = draw_poly_blades(x,y,r,plenum,y_p_stag)
                     if (chi>0.0_rk) then
                        mask_color(ix,iy)  = color_blades
                        mask(ix,iy,UxF)    = mask(ix,iy,UxF) + chi
                        mask(ix,iy,UyF)    = mask(ix,iy,UyF) + chi
                        mask(ix,iy,pF)     = mask(ix,iy,pF)  + chi
       !                 mask(ix,iy,rhoF) = mask(ix,iy,rhoF) + chi
                     endif
               endif
            endif
    
         case('default')
         end select

                  

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
    integer(kind=ik), intent(in)               :: g          !< grid parameter
    integer(kind=ik), dimension(3), intent(in) :: Bs
    real(kind=rk), intent(in)                  :: x0(2), dx(2)   !< coordinates of block and block spacinf
    real(kind=rk), intent(inout)               :: mask(:,:,:)    !< mask function
    integer(kind=2), intent(inout), optional   :: mask_color(:,:)!< identifyers of mask parts (plates etc)
    ! -----------------------------------------------------------------
    real(kind=rk)     :: x, y, r, h, xstart, xEnd
    real(kind=rk)     :: chi,inibox_length
    real(kind=rk)     :: d_ple
    integer(kind=ik)  :: ix, iy,n ! loop variables
    integer(kind=ik)  :: set_inibox ! flag for initial conditions

  ! -----plenum%length_pip+plenum%wall_thickness------------------------------------------------------------
!INITIALISE PARAMETERS
!sp_thickness = plenum%sponge_thickness
d_ple = plenum%diameter_ple
    ! parameter for smoothing function (width)
    h  = 2.5_rk*max(dx(1), dx(2))

    ! smooth width in x and y direction
    do iy=g+1, Bs(2)+g
       y = dble(iy-(g+1)) * dx(2) + x0(2)
       r = abs(y-domain_size(2)*0.5_rk)
        do ix=g+1, Bs(1)+g
            x = dble(ix-(g+1)) * dx(1) + x0(1)    
           ! AUSKOMMENTIERT FÜR PRODUCTION RUN 1 
            ! WEST DOMAIN SPONGE
             xstart = domain_size(1)-3.0_rk*plenum%wall_thickness
             xEnd = domain_size(1)           
             chi = draw_double_sponge(x,y,r,plenum,h,xstart,xEnd)   
             if(chi>0.0_rk) then
                    mask_color(ix,iy) = color_west_sponge
                    mask(ix,iy,rhoF)  = mask(ix,iy,rhoF)+chi
                    mask(ix,iy, pF)   = mask(ix,iy,pF)+chi
                    !mask(ix,iy, UxF)  = mask(ix,iy,UxF)+chi
                    !mask(ix,iy, UyF) = mask(ix,iy, UyF)+chi
            endif       
           chi = draw_plenum_rad_sponge(x, r, plenum, h)
           if (chi>0.0_rk) then
                    mask_color(ix,iy) = color_rad_sponge
                    mask(ix,iy,rhoF)  = mask(ix,iy,rhoF)+chi
                    mask(ix,iy, pF)   = mask(ix,iy,pF)+chi
                    mask(ix,iy, UxF)  = mask(ix,iy,UxF)+chi
                    mask(ix,iy, UyF) = mask(ix,iy, UyF)+chi 
           endif
           ! AUSKOMMENTIERT FÜR FU 
           !EAST DOMAIN WALL
            if (x< plenum%wall_thickness) then
                chi = 1.0_rk
                if (chi>0.0_rk) then
                       mask_color(ix,iy) = color_east_wall
                       mask(ix,iy,UxF)   = mask(ix,iy,UxF)+chi
                       mask(ix,iy,UyF)   = mask(ix,iy,UyF)+chi
                      mask(ix,iy,pF)    = mask(ix,iy,pF) +chi
                endif
            endif
            ! SET A BOX SHAPED INITIAL CONDITION IN THE PIPE FOR A STATIONARY FLOW
            set_inibox = plenum%set_inibox
            if (set_inibox ==1) then
            ! ENDE AUSKOMMENTIEREN FU
               ! SPECIAL HARDCODE FU INLET
               inibox_length = plenum%inibox_length
               if (x>= 0.01 .and. x<=0.01+inibox_length .and. r <=0.5_rk*plenum%diameter_pip) then
                   chi = 1.0_rk
                   if (chi>0.0_rk) then
                       mask_color(ix,iy) = color_inlet
                       mask(ix,iy,UxF)   = mask(ix,iy,UxF)+chi
                       mask(ix,iy,rhoF)  = mask(ix,iy,rhoF)+chi
                       mask(ix,iy,pF)    = mask(ix,iy,pF) +chi
                   endif
               endif
            endif
            ! END SPECIAL HARDCODE FU INLET
           
        
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
real(kind=rk) :: rho_sponge, p_sponge,v_sponge,temperature
integer(kind=ik) :: ix,iy,n
real(kind=rk)    :: length_pip,diameter_pip,radius_pip
real(kind=rk)    ::u_in, rho_in, p_in
!real(kind=rk)    ::sponge_thickness
! INITIALIZE PARAMETERS
! -------------------------------------
phi_ref=0.0_rk
temperature = plenum%temperature
!--------------------------------------
rho_sponge = plenum%sponge_density
p_sponge = plenum%sponge_pressure
length_pip = plenum%length_pip
diameter_pip = plenum%diameter_pip
radius_pip = 0.5_rk*diameter_pip
u_in = plenum%inlet_velocity
p_in = plenum%inlet_pressure
rho_in=plenum%inlet_density
!--------------------------------------
 !set smoothing parameter
h  = 2.5_rk*max(dx(1), dx(2))
    do iy=g+1, Bs(2)+g
       y = dble(iy-(g+1)) * dx(2) + x0(2)
       r = abs(y-domain_size(2)*0.5_rk)
       do ix=g+1,Bs(1)+g
            x = dble(ix-(g+1)) * dx(1) + x0(1)

            rho = phi(ix,iy,rhoF)
            u   = phi(ix,iy,UxF)
            v   = phi(ix,iy,UyF)
            p   = phi(ix,iy,pF)
            C_inv=C_eta_inv

     !what happens in solid obstacles (walls) no slip wall is penalized by u=v=0, rho == arbitrary value and p == computed via ideal gas law
            ! LEFT WALL
               if (mask_color(ix,iy)==color_east_wall) then
               
                   Phi_ref(ix,iy,UxF)   = 0.0_rk ! no velocity in x
                   
                   Phi_ref(ix,iy,UyF)   = 0.0_rk ! no velocity in y
                   
                   Phi_ref(ix,iy,pF)    = rho*Rs*temperature ! 
                   
                   !C_inv = C_eta_inv
               end if

               if (mask_color(ix,iy)==color_west_sponge) then
               
                  !Phi_ref(ix,iy,UxF)    = 275.0_rk ! no velocity in x
                  
                  !Phi_ref(ix,iy,UyF)    = 0.0_rk ! no velocity in y
                  
                  Phi_ref(ix,iy, rhoF)  = plenum%sponge_density!
                  
                  Phi_ref(ix,iy,pF)     = plenum%sponge_pressure!
                  
                  !C_inv = C_eta_inv
               endif
               
               
               
               
               if (mask_color(ix,iy)==color_rad_sponge) then
               
                  Phi_ref(ix,iy,UxF)    = 0.0_rk ! no velocity in x
                  
                  Phi_ref(ix,iy,UyF)    = 0.0_rk ! no velocity in y
                  
                  Phi_ref(ix,iy, rhoF)  = plenum%sponge_density!
                  
                  Phi_ref(ix,iy,pF)     = plenum%sponge_pressure!
                  
                  !C_inv = C_eta_inv
               endif
               
               if (mask_color(ix,iy)==color_walls) then
               
                  
                  Phi_ref(ix,iy,UxF)   = 0.0_rk ! no velocity in x
                  
                  Phi_ref(ix,iy,UyF)   = 0.0_rk ! no velocity in y
                  
             !     Phi_ref(ix,iy,pF)    = rho*Rs*temperature
                  
                  !C_inv = C_eta_inv
                  
                  
               end if
              if (mask_color(ix,iy)==color_blades) then
               
                  
                  Phi_ref(ix,iy,UxF)   = 0.0_rk ! no velocity in x
                  
                  Phi_ref(ix,iy,UyF)   = 0.0_rk ! no velocity in y
                  
                  Phi_ref(ix,iy,pF)    = rho*Rs*temperature
                  
                  !C_inv = C_eta_inv
                  
                  
               end if
               
               ! FU HARDCODE
               if (mask_color(ix,iy)==color_inlet) then
               
                  
                 Phi_ref(ix,iy,UxF)   =  u_in! no velocity in x
                  
                 Phi_ref(ix,iy,rhoF)  = rho_in  ! no velocity in y
                  
                 Phi_ref(ix,iy,pF)    =  p_in
                  
                  !C_inv = C_eta_inv
                  
                  
               end if
               !ENDE FU HARDCODE
               mask(ix,iy,:) = C_inv*mask(ix,iy,:)
       end do
    end do

end subroutine compute_plenum_penal2D

!=======================================================================

!###SECTION FUNCTIONS

!=======================================================================

function draw_plenum_walls(x,y, r,plenum,h)
implicit none
!-----------------------------------------------------------------------!
  real(kind=rk)                         :: x,r, h
  type(type_plenum),intent(in)          ::plenum
  real(kind=rk)                         ::y
  real(kind=rk)                         ::mask, draw_plenum_walls
  real(kind=rk)                         ::r0_pip,r0_ple,width,l_pip
  real(kind=rk)                         ::sp1, sp2,sp3
  real(kind=rk)                         ::wall_thickness 
  character(len=90)                     :: geo_case
  real(kind=rk)                         ::geom_offset
!-----------------------------------------------------------------------!
  mask=0.0_rk
  geom_offset = 8.0_rk/1000.0_rk
  wall_thickness=plenum%wall_thickness
  r0_ple=0.5_rk*plenum%diameter_ple-plenum%wall_thickness
  width=R_domain
  l_pip = plenum%length_pip
  geo_case = plenum%name
  sp1 = 0.5_rk*0.3_rk+wall_thickness+geom_offset!/0.8_rk*l_pip+plenum%wall_thickness+geom_offset
  sp2 = 0.5_rk*0.3206_rk+wall_thickness+geom_offset!/0.8_rk*l_pip+plenum%wall_thickness+geom_offset
  sp3 = 0.5_rk*0.634_rk+wall_thickness+geom_offset!/0.8_rk*l_pip+plenum%wall_thickness+geom_offset!*2.0_rk+wall_thickness
  
!-----------------------------------------------------------!
!-----------------____________________----------------------!
!----------------|CONSTANT RADIUS PIPE|---------------------!
!----------------|____________________|---------------------!
  select case(geo_case)
  case('const')
  r0_pip=0.5_rk*plenum%diameter_pip
  if (x<=l_pip) then
  mask=soft_bump(r,r0_pip,width,h)
  endif

!-----------------------------------------------------------!
!-----------------____________------------------------------!
!----------------|SFC_PDC_CASE|-----------------------------!
!----------------|____________|-----------------------------!
  case('sfb_pdc')
  if (x<=sp1.and. x>plenum%wall_thickness)then 
     r0_pip = 0.0403_rk
     !mask=soft_bump(r,r0_pip,width,h)
     if (r>r0_pip) then
     mask = 1.0
     endif
     
  elseif (x>sp1 .and. x<=sp2) then
     r0_pip = 0.0403_rk - (x-sp1)/(sp1-sp2)*(0.02_rk-0.0403_rk)
     !mask=soft_bump(r,r0_pip,width,h)!(h needed for soft_bump
     if (r>r0_pip) then
     mask = 1.0
     endif
     
  elseif(x>sp2 .and. x<=sp3) then
     r0_pip = 0.02_rk + (x-sp2)/(sp2-sp3)*(0.02_rk-0.0403_rk)
     !mask=soft_bump(r,r0_pip,width,h)
     if (r>r0_pip) then
     mask = 1.0
     endif
     
  elseif(x>sp3 .and. x <=l_pip) then
     r0_pip = 0.0403_rk
     !mask=soft_bump(r,r0_pip,width,h)
     if (r>r0_pip) then
     mask = 1.0
     endif
  !elseif(x>l_pip) then
     !r0_pip = 0.5_rk*plenum%diameter_ple
!     mask=soft_bump(r,r0_pip,width,h)
   !  if (r>r0_pip) then
    ! mask = 1.0
     !endif
  endif
  
  case default
  call abort(8546533,"ERROR: DONT YOU FORGET ABOUT CHOSING A GEOMETRY. "//plenum%name)
  
  
  end select  
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
x_m = plenum%length_pip+ blade_gap+0.5*length_blade
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
       
       if((y<0) .or. (y>domain_size(2)))  then
          call abort(90790,"Error: The blades do not fit the plenum domain. Reduce number of blades or blade length")
       endif
  endif

draw_stat_blades=mask

end function draw_stat_blades

function draw_poly_blades(x,y,r,plenum,y_p_stag)
implicit none
!-----------------------------------------------------------------!
real(kind=rk), dimension(9)   :: low_coeffs,up_coeffs
real(kind=rk),intent(in)      ::x,y,r,y_p_stag
type(type_plenum),intent(in)  ::plenum
real(kind=rk)                 ::y_up,y_low,x_rel,y_rel,gap_blade
real(kind=rk)                 ::draw_poly_blades,mask
!-----------------------------------------------------------------!

gap_blade = plenum%gap_blade
x_rel = x-(plenum%length_pip+plenum%gap_blade)
y_rel = y-y_p_stag
up_coeffs = plenum%upper_coefficients
low_coeffs = plenum%lower_coefficients
!up_coeffs(1) =plenum%upper_coefficients(1)
!up_coeffs(2) = -168937742.016
!up_coeffs(3) =31724549.006
!up_coeffs(4) = -3110982.882
!up_coeffs(5) = 170019.450
!up_coeffs(6) = -5078.268
!up_coeffs(7) = 75.118
!up_coeffs(8) = -0.014
!up_coeffs(9) =0.002

!low_coeffs(1) =-30421126.879
!low_coeffs(2) = 18940577.023
!low_coeffs(3) =-4031939.655
!low_coeffs(4) = 404883.625
!low_coeffs(5) =-21172.482
!low_coeffs(6) =568.821
!low_coeffs(7) =-0.726
!low_coeffs(8) =-0.053
!low_coeffs(9) = 0.000
! FIND DEFINITION ZONE OF POLYNOMIALS [-0.005,0.05]
y_up = up_coeffs(1)*x_rel**8+up_coeffs(2)*x_rel**7+up_coeffs(3)*x_rel**6+up_coeffs(4)*x_rel**5+up_coeffs(5)*x_rel**4+up_coeffs(6)*x_rel**3+up_coeffs(7)*x_rel**2+up_coeffs(8)*x_rel**1+up_coeffs(9)
! HORNER SCHEMA n multiplications instead of 2*n-1
!y_up = (((((((((up_coeffs(1)*x_rel)*x_rel+up_coeffs(2))*x_rel+up_coeffs(3))*x_rel+up_coeffs(4))*x_rel+up_coeffs(5))*x_rel+up_coeffs(6))*x_rel+up_coeffs(7))*x_rel+up_coeffs(8))*x_rel+up_coeffs(9))
y_low = low_coeffs(1)*x_rel**8+low_coeffs(2)*x_rel**7+low_coeffs(3)*x_rel**6+low_coeffs(4)*x_rel**5+low_coeffs(5)*x_rel**4+low_coeffs(6)*x_rel**3+low_coeffs(7)*x_rel**2+low_coeffs(8)*x_rel**1+low_coeffs(9)
!HORNER SCHEMA
!y_low = (((((((((low_coeffs(1)*x_rel)*x_rel+low_coeffs(2))*x_rel+low_coeffs(3))*x_rel+low_coeffs(4))*x_rel+low_coeffs(5))*x_rel+low_coeffs(6))*x_rel+low_coeffs(7))*x_rel+low_coeffs(8))*x_rel+low_coeffs(9))

if (y_rel>=y_low .and. y_rel<=y_up) then
    mask = 1.0
endif

draw_poly_blades=mask



end function draw_poly_blades




function draw_plenum_rad_sponge(x, r, plenum, h)
implicit none
!------------------------------------------------------------------!
  real(kind=rk),    intent(in)          :: x, r, h
  type(type_plenum),intent(in)          ::plenum

  real(kind=rk)                         ::draw_plenum_rad_sponge,mask
  real(kind=rk)                         ::r0,width
!-------------------------------------------------------------------!
  mask=0.0_rk
  r0=0.5_rk*(plenum%diameter_ple-plenum%wall_thickness)
  width=plenum%wall_thickness
if (x>=plenum%length_pip.and. &
   x<= domain_size(1)-3.0_rk*plenum%wall_thickness) then
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


function draw_double_sponge(x,y,r,plenum,h,x0,xEnd)
! draws a      ____
!             /    \
!           _/      \_ -shaped sponge to break periodic BC
implicit none
!--------------------------------------------------------------------------!
real(kind=rk), intent(in)   ::x,y,r,h,x0,xEnd
type(type_plenum),intent(in)::plenum
real(kind=rk)               ::draw_double_sponge
real(kind=rk)               ::mask
real(kind=rk)               ::x_hi1,x_hi2
real(kind=rk)               ::x3rd
real(kind=rk)               ::mu_tmp, a, b,pi
!--------------------------------------------------------------------------!
pi = 3.14159265359_rk
x3rd = (xEnd-x0)/3.0_rk
x_hi1 = x0+x3rd
x_hi2 = x0+2.0_rk*x3rd
       !RAMP UP
       if (x>x0 .and. x<=x_hi1) then
          a = 2.0_rk*pi/(x_hi1-x0);
          !write(*,*) pi
          b = -a*x0;
          mu_tmp = a*x+b
          !write(*,*) mu_tmp
          mask = tanh(mu_tmp)
          !write(*,*) mask
        endif
        ! CONSTANT MIDDLE PART
        if (x>x_hi1 .and. x<=x_hi2) then
           mask = 1.0_rk
           !write(*,*) mask
        endif
        ! RAMP DOWN
        if (x>x_hi2 .and. x<=xEnd) then
            a = 2.0_rk*pi/(x_hi2-xEnd)
            b = -a*xEnd
            mu_tmp = a*x+b
            !write(*,*) mu_tmp
            mask = tanh(mu_tmp)
            !write(*,*) mask
        endif
draw_double_sponge = mask
end function

