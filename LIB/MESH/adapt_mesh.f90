!> \file
!> \callgraph
!********************************************************************************************
! WABBIT
! ============================================================================================
!> \name     adapt_mesh.f90
!! \version  0.4
!! \author   msr, engels
!
!> \brief This routine performs the coarsing of the mesh, where possible. For the given mesh
!! we compute the details-coefficients on all blocks. If four sister blocks have maximum
!! details below the specified tolerance, (so they are insignificant), they are merged to
!! one coarser block one level below. This process is repeated until the grid does not change
!! anymore.
!!
!! As the grid changes, active lists and neighbor relations are updated, and load balancing
!! is applied.
!
!> \note The block thresholding is done with the restriction/prediction operators acting on the
!! entire block, INCLUDING GHOST NODES. Ghost node syncing is performed in threshold_block.
!
!> \note It is well possible to start with a very fine mesh and end up with only one active
!! block after this routine. You do *NOT* have to call it several times.
!
!> \details
!! input:    - params, light and heavy data \n
!! output:   - light and heavy data arrays
!!
!> = log ======================================================================================
!!
!! 10/11/16 - switch to v0.4
! ==========================================================================================
!********************************************************************************************
!> \image html adapt_mesh.svg width=400

subroutine adapt_mesh( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, &
    lgt_sortednumlist, hvy_active, hvy_n, indicator, com_lists, com_matrix, int_send_buffer,&
     int_receive_buffer, real_send_buffer, real_receive_buffer, hvy_synch )

!---------------------------------------------------------------------------------------------
! modules
    use module_indicators

!---------------------------------------------------------------------------------------------
! variables

    implicit none
    integer(kind=1), intent(inout)      :: hvy_synch(:, :, :, :)

    !> user defined parameter structure
    type (type_params), intent(in)      :: params
    !> light data array
    integer(kind=ik), intent(inout)     :: lgt_block(:, :)
    !> heavy data array
    real(kind=rk), intent(inout)        :: hvy_block(:, :, :, :, :)
    !> heavy data array - neighbor data
    integer(kind=ik), intent(inout)     :: hvy_neighbor(:,:)
    !> list of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_active(:)
    !> number of active blocks (light data)
    integer(kind=ik), intent(inout)     :: lgt_n
    !> sorted list of numerical treecodes, used for block finding
    integer(kind=tsize), intent(inout)   :: lgt_sortednumlist(:,:)
    !> list of active blocks (heavy data)
    integer(kind=ik), intent(inout)     :: hvy_active(:)
    !> number of active blocks (heavy data)
    integer(kind=ik), intent(inout)     :: hvy_n
    !> coarsening indicator
    character(len=*), intent(in)        :: indicator
    ! communication lists:
    integer(kind=ik), intent(inout)     :: com_lists(:, :, :, :)
    ! communications matrix:
    integer(kind=ik), intent(inout)     :: com_matrix(:,:,:)
    ! send/receive buffer, integer and real
    integer(kind=ik), intent(inout)      :: int_send_buffer(:,:), int_receive_buffer(:,:)
    real(kind=rk), intent(inout)         :: real_send_buffer(:,:), real_receive_buffer(:,:)

    ! loop variables
    integer(kind=ik)                    :: lgt_n_old, iteration, k, max_neighbors
    ! cpu time variables for running time calculation
    real(kind=rk)                       :: t0, t1, t_misc
    ! MPI error variable
    integer(kind=ik)                    :: ierr
    logical::test

!---------------------------------------------------------------------------------------------
! variables initialization

    ! start time
    t0 = MPI_Wtime()
    t1 = t0
    t_misc = 0.0_rk
    lgt_n_old = 0
    iteration = 0

    if ( params%threeD_case ) then
        max_neighbors = 56
    else
        max_neighbors = 12
    end if

!---------------------------------------------------------------------------------------------
! main body

    !> we iterate until the number of blocks is constant (note: as only coarseing
    !! is done here, no new blocks arise that could compromise the number of blocks -
    !! if it's constant, its because no more blocks are refined)
    do while ( lgt_n_old /= lgt_n )

        lgt_n_old = lgt_n

        !> (a) check where coarsening is possible
        ! ------------------------------------------------------------------------------------
        ! first: synchronize ghost nodes - thresholding on block with ghost nodes
        ! synchronize ghostnodes, grid has changed, not in the first one, but in later loops
        !call synchronize_ghosts( params, lgt_block, hvy_block, hvy_neighbor, hvy_active, hvy_n, com_lists(1:hvy_n*max_neighbors,:,:,:), com_matrix, .true., int_send_buffer, int_receive_buffer, real_send_buffer, real_receive_buffer )
        call check_redundant_nodes( params, lgt_block, hvy_block, hvy_synch, hvy_neighbor, hvy_active, &
hvy_n, int_send_buffer, int_receive_buffer, real_send_buffer, real_receive_buffer, test )
        ! calculate detail
        call coarsening_indicator( params, lgt_block, hvy_block, lgt_active, lgt_n, hvy_active, hvy_n, indicator, iteration)


        !> (b) check if block has reached maximal level, if so, remove refinement flags
        t0 = MPI_Wtime()
        call respect_min_max_treelevel( params, lgt_block, lgt_active, lgt_n )
        ! CPU timing (only in debug mode)
        call toc( params, "adapt_mesh (min/max)", MPI_Wtime()-t0 )


        !> (c) unmark blocks that cannot be coarsened due to gradedness
        t0 = MPI_Wtime()
        call ensure_gradedness( params, lgt_block, hvy_neighbor, lgt_active, lgt_n )
        ! CPU timing (only in debug mode)
        call toc( params, "adapt_mesh (gradedness)", MPI_Wtime()-t0 )


        !> (d) ensure completeness
        t0 = MPI_Wtime()
        call ensure_completeness( params, lgt_block, lgt_active, lgt_n, lgt_sortednumlist )
        ! CPU timing (only in debug mode)
        call toc( params, "adapt_mesh (completeness)", MPI_Wtime()-t0 )


        !> (e) adapt the mesh, i.e. actually merge blocks
        t0 = MPI_Wtime()
        call coarse_mesh( params, lgt_block, hvy_block, lgt_active, lgt_n, lgt_sortednumlist )
        ! CPU timing (only in debug mode)
        call toc( params, "adapt_mesh (coarse mesh)", MPI_Wtime()-t0 )


        ! the following calls are indeed required (threshold->ghosts->neighbors->active)
        ! update lists of active blocks (light and heavy data)
        ! update list of sorted nunmerical treecodes, used for finding blocks
        t0 = MPI_Wtime()
        call create_active_and_sorted_lists( params, lgt_block, lgt_active, lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true. )
        t_misc = t_misc + (MPI_Wtime() - t0)


        t0 = MPI_Wtime()
        ! update neighbor relations
        call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n )
        ! CPU timing (only in debug mode)
        call toc( params, "adapt_mesh (update neighbors)", MPI_Wtime()-t0 )

        iteration = iteration + 1
    end do

    !> At this point the coarsening is done. All blocks that can be coarsened are coarsened
    !! they may have passed several level also. Now, the distribution of blocks may no longer
    !! be balanced, so we have to balance load now
    call balance_load( params, lgt_block, hvy_block, hvy_neighbor, lgt_active, lgt_n, hvy_active, hvy_n )


    !> load balancing destroys the lists again, so we have to create them one last time to
    !! end on a valid mesh
    !! update lists of active blocks (light and heavy data)
    ! update list of sorted nunmerical treecodes, used for finding blocks
    t0 = MPI_wtime()
    call create_active_and_sorted_lists( params, lgt_block, lgt_active, lgt_n, hvy_active, hvy_n, lgt_sortednumlist, .true. )
    t_misc = t_misc + (MPI_Wtime() - t0)


    ! update neighbor relations
    t0 = MPI_Wtime()
    call update_neighbors( params, lgt_block, hvy_neighbor, lgt_active, lgt_n, lgt_sortednumlist, hvy_active, hvy_n )
    ! CPU timing (only in debug mode)
    call toc( params, "adapt_mesh (update neighbors) ", MPI_Wtime()-t0 )


    ! time remaining parts of this routine.
    call toc( params, "adapt_mesh (...)", t_misc )
    call toc( params, "adapt_mesh (TOTAL)", MPI_wtime()-t1)
end subroutine adapt_mesh
