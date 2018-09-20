  !-----------------------------------------------------------------------------
  !>\file
  !> Initial Conditions of the Navier Stokes module
  !-----------------------------------------------------------------------------
  subroutine INICOND_NStokes( time, u, g, x0, dx )

    use module_shock, only : set_shock_1D,moving_shockVals,standing_shockVals

    implicit none
    ! it may happen that some source terms have an explicit time-dependency
    ! therefore the general call has to pass time
    real(kind=rk), intent (in) :: time

    ! block data, containg the state vector. In general a 4D field (3 dims+components)
    ! in 2D, 3rd coindex is simply one. Note assumed-shape arrays
    real(kind=rk), intent(inout) :: u(1:,1:,1:,1:)

    ! as you are allowed to compute the RHS only in the interior of the field
    ! you also need to know where 'interior' starts: so we pass the number of ghost points
    integer, intent(in) :: g

    ! for each block, you'll need to know where it lies in physical space. The first
    ! non-ghost point has the coordinate x0, from then on its just cartesian with dx spacing
    real(kind=rk), intent(in) :: x0(1:3), dx(1:3)

    integer(kind=ik)          :: Bs,ix,iy,iz,dF
    real(kind=rk)             :: x,y_rel,tmp(1:3),b,mach,T_init,&
                                left(size(u,4)),right(size(u,4)),phi_init(size(u,4))
    real(kind=rk),allocatable:: mask(:,:,:) ! we dont save this datafield, since it is only called once

    ! compute the size of blocks
    Bs = size(u,1) - 2*g
    !------------------------------------------------------------------
    ! READ INITIAL CONDITION FROM FILE
    !------------------------------------------------------------------
    ! NStokes saves the temporary state in the "pure" format, i.e. rho, u, v, w, p
    ! This is why we have to convert them back to skew-symetric, when they have been
    ! red from a file
    if ( params_ns%read_from_files ) then
      ! convert (rho,u,v,p) to (sqrt(rho),sqrt(rho)u,sqrt(rho)v,p) if data was read from file
        call pack_statevector(u(:,:,:,:),'pure_variables')
        return
    else
        ! if not read form file, then reset the statevector
        u = 0.0_rk
    endif
    !------------------------------------------------------------------
    ! CASE SPECIFIC INITIIAL CONDITIONS
    !------------------------------------------------------------------
    ! generate case specific initial conditions:
    ! if the initial conditions are case specific they are set here
    if (set_inicond_case(params_ns,x0, dx, Bs, g, u)) return
    !------------------------------------------------------------------

    !------------------------------------------------------------------
    ! GENERIC INITIIAL CONDITIONS
    !------------------------------------------------------------------
    ! if set_inicond_case is .false. we provide some generic initial
    ! conditions here, which are available for every case study.
    phi_init(pF)    =params_ns%initial_pressure
    phi_init(rhoF)  =params_ns%initial_density
    phi_init(UxF)    =params_ns%initial_velocity(1)
    phi_init(UyF)    =params_ns%initial_velocity(2)
    if (params_ns%dim==3) then
      phi_init(UzF)    =params_ns%initial_velocity(3)
    end if
    T_init    =params_ns%initial_temp
    ! check if the input values of the initial conditions make sense
    if (phi_init(pF)<=0.0_rk .or. phi_init(rhoF) <=0.0) then
      call abort(6032,  "Error [module_navier_stokes.f90]: initial pressure and density "//&
      "must be larger then 0")
    endif
    ! The generic initial conditions are selected here
    select case( params_ns%inicond )
    case("shear_layer")
      if (params_ns%dim==2) then
        call inicond_shear_layer(  u, x0, dx ,Bs, g)
      else
        call abort(4832,"ERROR [navier_stokes_new.f90]: no 3d shear layer implemented")
      endif

    case ("zeros")
      ! add ambient pressure
      u( :, :, :, pF)   = params_ns%initial_pressure
      u( :, :, :, rhoF) = sqrt(phi_init(rhoF))
      u( :, :, :, UxF)  = 0.0_rk
      u( :, :, :, UyF)  = 0.0_rk
      if (params_ns%dim==3)  u( :, :, :, UzF) = 0.0_rk

    case ("standing-shock","moving-shock")
      ! following values are imposed and smoothed with tangens:
      ! ------------------------------------------
      !   rhoL    | rhoR                  | rho
      !   uL      | uR                    | uL
      !   pL      | pR                    | pL
      ! 0-----------------------------------------xLength
      !           x0_shock                0.95*Length
      ! reset the statevector values left and right of the shock
      right=0
      left =0
      ! compute the values from the initial conditions which have been read from params file
      if ( params_ns%inicond == "standing-shock" ) then
        ! in standing_shockVals we compute the shock values of the other side (i.e right) of the shock
        call standing_shockVals(phi_init(rhoF),phi_init(UxF),phi_init(pF),right(rhoF),right(UxF),right(pF),params_ns%gamma_)
      else
        ! compute the standing_shockVals of the left and right values
        call moving_shockVals(phi_init(rhoF),phi_init(UxF),phi_init(pF), right(rhoF),right(UxF),right(pF) &
                              ,params_ns%gamma_,params_ns%machnumber)
        params_ns%initial_velocity(1)=phi_init(UxF)
      end if

      left     =phi_init  ! the values to the left of the shock are taken from the ini file
      left(UyF)= 0        ! make sure that only the z and y component of the velocity is 0
      if (params_ns%dim==3) left(UzF)= 0        ! make sure that only the z and y component of the velocity is 0

      ! check for usefull inital values
      if ( right(rhoF)<0 .or. right(pF)<0 ) then
        call abort(3572,"ERROR: initial values are insufficient for simple-shock")
      end if
    ! the shock tube is only 1D, therefore we loop over the second component in 2D and
    ! 2cnd and 3rd component in 3D
    if (params_ns%dim==3) then
      do iz = 1, Bs+2*g
        do iy = 1, Bs+2*g
          call set_shock_1D(x0(1), dx(1), Bs, g, u(:,iy,iz,:), &
           left, right, params_ns%inicond_width,params_ns%domain_size(1))
        end do
      end do
    else
      do iy = 1, Bs+2*g
        call set_shock_1D(x0(1), dx(1), Bs, g, u(:,iy,1,:), left, right, &
                          params_ns%inicond_width, params_ns%domain_size(1))
      end do
    endif
    call pack_statevector(u,'pure_variables')

    case ("mask")
      if (.not. params_ns%penalization) call abort(110918,"SAY WHAAAT? can't hear you! Please switch on penalization for inicond mask!")
      if (.not. allocated(mask))        allocate(mask(size(u,1), size(u,2), size(u,3)))

      ! set velocity field u(x)=1 for x in mask
      ! u(x)=(1-mask(x))*u0 to make sure that flow is zero at mask values
      call get_mask(params_ns, x0, dx, Bs, g , mask)
      u( :, :, :, pF) = phi_init(rhoF)*params_ns%Rs*T_init
      u( :, :, :, rhoF) = sqrt(phi_init(rhoF))
      u( :, :, :, UxF) = ( 1 - mask ) * phi_init(UxF)*sqrt(phi_init(rhoF)) !flow in x
      u( :, :, :, UyF) = (1-mask)*phi_init(UyF)*sqrt(phi_init(rhoF)) !flow in y
      if (params_ns%dim==3) then
        u( :, :, :, UzF) = (1-mask)*phi_init(UyF)*sqrt(phi_init(rhoF)) !flow in z
      endif

    case ("pressure_blob")
        ! pressure component has a gaus function with maximum given by the
        ! initial value of the pressure.
        ! all other statevariables are set constant to the intial values
        call gauss_function( params_ns%inicond_width,params_ns%dim, &
                                Bs,g, params_ns%domain_size, u(:,:,:,pF), x0, dx, phi_init(pF))
        u( :, :, :, rhoF)= sqrt(phi_init(rhoF))
        if (params_ns%dim==3) then
          u( :, :, :, UxF) = phi_init(UxF)*sqrt(phi_init(rhoF))
          u( :, :, :, UyF) = phi_init(UyF)*sqrt(phi_init(rhoF))
          u( :, :, :, UzF) = phi_init(UzF)*sqrt(phi_init(rhoF))
        else
          u( :, :, :, UxF) = phi_init(UxF)*sqrt(phi_init(rhoF))
          u( :, :, :, UyF) = phi_init(UyF)*sqrt(phi_init(rhoF))
        endif
    case default
        call abort(7771,"the initial condition is unkown: "//trim(adjustl(params_ns%inicond)))
    end select

    !------------------------------------------------------------------
    !check if initial conditions are set properly
    do dF=1,params_ns%n_eqn
      if (  block_contains_NaN(u(:,:,:,dF)) ) then
        call abort(46924,"ERROR [module_navier_stokes]: NaN in inicond!"//&
        " Computer says NOOOOOOOOOOOOOOOO!")
      endif
    end do
    !------------------------------------------------------------------

    if (allocated(mask)) deallocate(mask)

  end subroutine INICOND_NStokes





  !> \brief initialize shear layer setup, setup for all datafields here,
  subroutine inicond_shear_layer(  u, x0, dx ,Bs, g)
    implicit none
      !------------------------------------------------------------
      !> actual block data (note this routine acts only on one block)
      real(kind=rk), intent(inout) :: u(:,:,:,:)
      !> spacing and origin of block
      real(kind=rk), intent(in)    :: x0(1:3),dx(1:3)
      ! grid
      integer(kind=ik),intent(in)  :: Bs, g
      !------------------------------------------------------------
      ! variable for shear layer position
      real(kind=rk)                           :: mux1, mux2, muy, x, y, sigma, w
      ! loop variables
      integer(kind=ik)                        :: ix, iy
      ! p0 value \todo get from ini file
      real(kind=rk)                           :: p0, rho0, u0


      rho0 = params_ns%initial_density
      p0   = params_ns%initial_pressure
      u0   = params_ns%initial_velocity(1)

      ! place layer
      mux1 = params_ns%domain_size(1)/2.0_rk - 0.25_rk
      mux2 = params_ns%domain_size(1)/2.0_rk + 0.25_rk

      muy  = 0.5_rk * params_ns%domain_size(2)

      ! boundary layer width
      sigma = params_ns%inicond_width * params_ns%domain_size(1)
      ! shear layer width
      w = params_ns%inicond_width

      if (size(u,3)>1) then
        call abort(110918, 'Shear layer is only available for 2D calculations')
      else
          ! 2D case
          ! create shear layer, Uy field
          do ix = 1,Bs+2*g
              do iy = 1,Bs+2*g

                  ! compute x,y coordinates from spacing and origin
                  x = dble(ix-(g+1)) * dx(1) + x0(1)
                  y = dble(iy-(g+1)) * dx(2) + x0(2)
                  ! ensure that coordinates are continued periodicaly
                  ! otherwise problems in redudant nodes.
                  ! shear layer, setup [1]
                  ! Uy
                  if ( x <= 0.5_rk*params_ns%domain_size(1) ) then
                      u(ix, iy, 1, UyF) = dtanh( w/params_ns%domain_size(1) * ( x - mux1 ) ) + u0
                      u(ix, iy, 1, rhoF) = ( rho0 + dtanh( w * ( x - mux1 ) ) ) / 2.0_rk + rho0
                  else
                      u(ix, iy, 1, UyF) = dtanh( w/params_ns%domain_size(1) * ( mux2 - x ) ) + u0
                      u(ix, iy, 1, rhoF) = ( rho0 + dtanh( w * ( mux2 - x ) ) ) / 2.0_rk + rho0
                  end if

                  ! Ux
                  !u(ix, iy, 1, UxF) = 0.05_rk * dsin( 8.0_rk * pi * ( y - muy  ) / params_ns%domain_size(2) )
                  u(ix, iy, 1, UxF) = 0.1_rk * dsin( 2.0_rk * pi * ( y - muy  ) )
              end do
          end do

          u(:, :, 1, pF)   = p0
          u(:, :, 1, rhoF) = dsqrt(u(:, :, 1, rhoF))
          u(:, :, 1, UxF)  = u(:, :, 1, UxF) * (u(:, :, 1, rhoF))
          u(:, :, 1, UyF)  = u(:, :, 1, UyF) * (u(:, :, 1, rhoF))
      endif

  end subroutine inicond_shear_layer




    !-------------------------------------------------------------------------------
    !> \brief Initialize gauss pulse for 2D/3D case: with \f$\mu=L/2\f$   \n
    !> \f{eqnarray*}{
    !> u(x\mid \mu ,\sigma ^{2})=A \exp \left\{-{\frac {\Vert\vec{x}-\vec{\mu} \Vert^{2}}{2\sigma ^{2}}}\right\}
    !> \f}
    !> if not defined by user input \f$A={\frac {1}{(\sqrt {2\pi} \sigma)^d}}\f$
    !> (\f$d\f$ number of dimensions)
    !------------------------------------------------------------------------------
    subroutine gauss_function(sigma, dimension, Bs, g, L, u, x0, dx, amplitude )

        implicit none
        real(kind=rk), intent(inout)        :: u(:,:,:)                      !< actual block data
        real(kind=rk), intent(in)          :: x0(1:3),dx(1:3),L(1:3),sigma  !<standard diviation/width of gauss pulse \f$\sigma\f$
        integer(kind=ik),intent(in)        :: Bs, g, dimension              !< note if the dimension is 0 the gaus blob is not normalized
        real(kind=rk), intent(in),optional :: amplitude  !< maximum of the gauss blob \f$A\f$ is optional


        ! auxiliary variable for gauss pulse
        real(kind=rk)                           :: x_cntr(3), x, z ,y,A
        ! loop variables
        integer(kind=ik)                        :: ix, iy, iz

        ! center point of the pressure blob
        x_cntr = 0.5_rk * L
        ! amplitude
        if ( present(amplitude) ) then
          A = amplitude
        else
          ! if not present normalize the gauss blob
          A = 1/(dsqrt(2*PI)*sigma)**dimension
        end if

        if (dimension==2) then
          ! 2D case
          ! create gauss pulse
          do ix = g+1,Bs+g
            do iy = g+1,Bs+g
              ! compute x,y coordinates from spacing and origin
              x = dble(ix-(g+1)) * dx(1) + x0(1)
              y = dble(iy-(g+1)) * dx(2) + x0(2)
              ! shift to new gauss blob center
              call continue_periodic(x,L(1))
              call continue_periodic(y,L(2))
              ! set actual inicond gauss blob
              u(ix,iy,1) = A*dexp( -( (x-x_cntr(1))**2 + (y-x_cntr(2))**2 ) / (2*sigma**2) )
            end do
          end do

        else
          ! 3D case
          ! create gauss pulse
          do ix = g+1,Bs+g
            do iy = g+1,Bs+g
              do iz = g+1,Bs+g

                x = dble(ix-(g+1)) * dx(1) + x0(1)
                y = dble(iy-(g+1)) * dx(2) + x0(2)
                z = dble(iz-(g+1)) * dx(3) + x0(3)
                ! shift to new gauss blob center
                call continue_periodic(x,L(1))
                call continue_periodic(y,L(2))
                call continue_periodic(z,L(3))
                ! set actual inicond gauss blob
                u(ix,iy,iz) = A*dexp( -( (x-x_cntr(1))**2 + (y-x_cntr(2))**2 +(z-x_cntr(3))**2 ) / (2*sigma)**2 )
              end do
            end do
          end do
      endif
    end subroutine gauss_function
