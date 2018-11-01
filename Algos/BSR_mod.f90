module bsr

   use crossValid
   use coordDescent
   use utilities
   implicit none

   contains
   
   !================================================================================================!
   ! Block Separated Representation algorithm, where each subproblem is solved using LBFGS or GLMNET.
   ! *** Computation steps ***
   !  STEP 1: cross-validation procedure using a coarse lambda-grid, no warm starts 
   !          ==> returns 'lbd_star'
   !  STEP 2: cross-validation procedure using a fine lambda-grid, war starts 
   !          ==> returns 'lbd_opt'
   !  STEP 3: production run using lbd_opt
   !          ==> returns 'coeffs_mat'
   !
   ! Calculations are followed by validation error computation and diagnostics (timings & sparsity).
   !================================================================================================!
   subroutine BSR_resol(coeffs_mat, rank, d, basissize, N_train, varphi_train, Xdata_train, & 
    Ydata_train, nb_lbd_try_c, epsi_c, nb_lbd_try_f, beta1, beta2, kfold, iprint, opt_resol, int_resp, &
     lbd_grid_c, err_test_aver_c, stock_err_test_c, lbd_grid_f, err_test_aver_f, stock_err_test_f, &
      flag_cvg_out, lbd_star, lbd_opt, tstep1, tstep2, tstep3, tol_fobj_cv, tol_fobj_run)
   
   !inputs:
   integer, intent(in)                                 :: rank, d, basissize, N_train, nb_lbd_try_c, &
                                                          nb_lbd_try_f, iprint, kfold
   real *8, dimension(N_train,basissize,d), intent(in) :: varphi_train
   real *8, dimension(d, N_train), intent(in)          :: Xdata_train
   real *8, dimension(N_train), intent(in)             :: Ydata_train
   real *8, intent(in)                                 :: epsi_c, beta1, beta2, tol_fobj_cv, tol_fobj_run
   character(len=4), intent(in)                        :: int_resp, opt_resol
   
   !outputs:
   real *8, dimension(rank,d,basissize), intent(out)   :: coeffs_mat
   real *8, dimension(nb_lbd_try_c), intent(out)       :: lbd_grid_c, err_test_aver_c
   real *8, dimension(kfold,nb_lbd_try_c), intent(out) :: stock_err_test_c
   real *8, dimension(nb_lbd_try_f), intent(out)       :: lbd_grid_f, err_test_aver_f
   real *8, dimension(kfold,nb_lbd_try_f), intent(out) :: stock_err_test_f
   logical, intent(out)                                :: flag_cvg_out
   real *8, intent(out)                                :: lbd_star, lbd_opt, tstep1, tstep2, tstep3
   
   !local variables:
   integer                                             :: N, i, j, k, l, count_iter
   real *8, dimension(:), allocatable                  :: coeffs_init, vals
   real *8, dimension(:,:,:), allocatable              :: comp1d_train_init, coeffs_mat_init, &
                                                          comp1d_train
   real *8                                             :: lmax_lasso, lbd_max, lmin_lasso, lbd_min, &
                                                          lbd_step, t0, t1, lmax, lmin
   character(len=4)                                    :: opt_warm, opt_blockCD
   logical                                             :: flag_cvg
   
   ! Random initialization of unknown coefficients:
   N = rank*d*basissize
   allocate(coeffs_init(N))
   call random_seed()
   call random_number(coeffs_init)
   coeffs_init = 2.0D0*coeffs_init - 1.0D0 ! take values in the interval [-1,1]
   
!    do i=1,20
!      print*,'coeffs=', coeffs_init(i)
!    end do
!    STOP

   ! Recast initial coefficients in a 3d data structure:
   allocate(coeffs_mat_init(rank,d,basissize))
   do l = 1,rank
     do j = 1,d
       do k = 1,basissize
         coeffs_mat_init(l,j,k) = coeffs_init((l-1)*d*basissize+(j-1)*basissize+k) 
       end do
     end do
   end do  
   
   ! Evaluate 1d component functions of the model structure on each training point:
   allocate(comp1d_train_init(rank,d,N_train))
   do l = 1,rank
     do j = 1,d
       do i =1,N_train
         comp1d_train_init(l,j,i) = dot_product(coeffs_mat_init(l,j,:),varphi_train(i,:,j))
       end do
     end do
   end do
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! STEP 1: Cross-validation with NO warm starts over a coarse lambda-grid  !
   !         in order to get a rough estimation of penalization parameter    !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   ! --- Define a coarse lambda-grid in log10 scale
   allocate(vals(d))
   do j = 1,d
     vals(j) = ABS(dot_product(Xdata_train(j,:),Ydata_train(:)))/N_train
   end do
   lmax_lasso = maxval(vals) !classical lasso upper bound
   lbd_max = log10(lmax_lasso)
   lmin_lasso = epsi_c*lmax_lasso
   lbd_min = log10(lmin_lasso)
   deallocate(vals)
   lbd_step = (lbd_max-lbd_min)/float(nb_lbd_try_c-1)
   do i = 1, nb_lbd_try_c
     lbd_grid_c(i) = 10**(lbd_min + (i-1)*lbd_step)
!      print*, lbd_grid_c(i)
   end do

   print *, ' Coarse grid for lambda based on lasso:'
   print *, '  ... lbd_min =', 10**lbd_min
   print *, '  ... lbd_max =', 10**lbd_max
   print *, '  ... lbd tries =', nb_lbd_try_c
   print *, ''
  
   ! --- cross-validation procedure
   allocate(comp1d_train(rank,d,N_train))
   opt_warm = 'n' !option for warm starts ('y' or 'n')
   opt_blockCD = 'y' !option for block version ('y' or 'n')
   coeffs_mat(:,:,:) = coeffs_mat_init(:,:,:)
   comp1d_train(:,:,:) = comp1d_train_init(:,:,:)
   call CPU_TIME(t0)
   call crossvalid_kfold_SR(coeffs_mat, Ydata_train, N_train, rank, d, basissize, & 
         varphi_train, nb_lbd_try_c, lbd_grid_c, comp1d_train, stock_err_test_c, & 
          err_test_aver_c, iprint, kfold, int_resp, opt_blockCD, opt_warm, opt_resol, &
           lbd_star, tol_fobj_cv)
   call CPU_TIME(t1)
   tstep1 = t1-t0

   print*, ''
   print*, 'lbd_star', lbd_star
      
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! STEP 2: Cross-validation with warm starts over a fine lambda-grid       !
   !         in order to get an optimal value for the penalization parameter !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   ! --- Define a fine lambda-grid in log10 scale:
   lmax = beta2*lbd_star
   lbd_max = log10(lmax)
   lmin = beta1*lbd_star
   lbd_min = log10(lmin)
   lbd_step = (lbd_max-lbd_min)/float(nb_lbd_try_f-1)
   do i = 1, nb_lbd_try_f
     lbd_grid_f(i) = 10**(lbd_min + (i-1)*lbd_step)
!      print*,lbd_grid_f(i)
   end do
   
   print *, ' Refined grid for lambda:'
   print *, '  ... lbd_min =', 10**lbd_min
   print *, '  ... lbd_max =', 10**lbd_max
   print *, '  ... lbd tries =', nb_lbd_try_f
   print *, ''
!    print*, 'paused, type [enter] to continue'
!    read (*,*) 
   
   ! --- Preliminary run (using 'lmax') to generate initial coefficients 
   !     and corresponding 1-d component functions (needed to initialize
   !     the refined cross-validation procedure)
   coeffs_mat(:,:,:) = coeffs_mat_init(:,:,:)
   comp1d_train(:,:,:) = comp1d_train_init(:,:,:)
   select case (opt_resol)
     case('bfgs')
       call BSR_subpbs_bfgs(coeffs_mat, Ydata_train, N_train, rank, d, basissize, & 
             varphi_train, lmax, comp1d_train, iprint, count_iter, int_resp, flag_cvg, tol_fobj_cv)
     case('glmn')
       call BSR_subpbs_glmnet(coeffs_mat, Ydata_train, N_train, rank, d, basissize, & 
             varphi_train, lmax, comp1d_train, iprint, count_iter, int_resp, flag_cvg, tol_fobj_cv)
   end select  
   
   ! --- Refined cross-validation procedure
   call CPU_TIME(t0)
   opt_warm = 'y' !option for warm starts ('y' or 'n')
   call crossvalid_kfold_SR(coeffs_mat, Ydata_train, N_train, rank, d, basissize, & 
         varphi_train, nb_lbd_try_f, lbd_grid_f, comp1d_train, stock_err_test_f, & 
          err_test_aver_f, iprint, kfold, int_resp, opt_blockCD, opt_warm, opt_resol, &
           lbd_opt, tol_fobj_cv)
   call CPU_TIME(t1)
   tstep2 = t1-t0
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! STEP 3: Production run using lambda_opt                                 !
   ! Returns the set of coefficients of the model structure 'coeffs_mat'     !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   coeffs_mat(:,:,:) = coeffs_mat_init(:,:,:)
   comp1d_train(:,:,:) = comp1d_train_init(:,:,:)
   call CPU_TIME(t0)
   select case (opt_resol)
     case('bfgs')
       call BSR_subpbs_bfgs(coeffs_mat, Ydata_train, N_train, rank, d, basissize, & 
             varphi_train, lbd_opt, comp1d_train, iprint, count_iter, int_resp, flag_cvg, tol_fobj_run)
     case('glmn')
       call BSR_subpbs_glmnet(coeffs_mat, Ydata_train, N_train, rank, d, basissize, & 
             varphi_train, lbd_opt, comp1d_train, iprint, count_iter, int_resp, flag_cvg, tol_fobj_run)
   end select  
   call CPU_TIME(t1)
   tstep3 = t1-t0
   
   flag_cvg_out = flag_cvg
   
   deallocate(coeffs_init, coeffs_mat_init, comp1d_train, comp1d_train_init)

   end subroutine BSR_resol
   !==========================================================================================!

end module bsr
