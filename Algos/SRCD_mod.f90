module srcd

   use crossValid
   use coordDescent
   use utilities
   implicit none

   contains

   subroutine srcd_resol(coeffs_mat, rank, d, basissize, N_train, varphi_train, Xdata_train, & 
    Ydata_train, nb_lbd_try, epsi, kfold, iprint, opt_resol, int_resp, flag_cvg_out, lbd_grid, & 
     err_test_aver, stock_err_test, lbd_opt, tstep1, tstep2, bias_cvg, tol_fobj_cv, tol_fobj_run)

   !inputs:
   integer, intent(in)                                 :: rank, d, basissize, N_train, nb_lbd_try, &
                                                          kfold, iprint
   real *8, dimension(N_train,basissize,d), intent(in) :: varphi_train
   real *8, dimension(d, N_train), intent(in)          :: Xdata_train
   real *8, dimension(N_train), intent(in)             :: Ydata_train
   real *8, intent(in)                                 :: epsi, tol_fobj_cv, tol_fobj_run
   character(len=4), intent(in)                        :: int_resp, opt_resol
   
   !outputs:
   real *8, dimension(rank,d,basissize), intent(out)   :: coeffs_mat
   real *8, dimension(nb_lbd_try), intent(out)         :: lbd_grid, err_test_aver
   real *8, dimension(kfold,nb_lbd_try), intent(out)   :: stock_err_test
   logical, intent(out)                                :: flag_cvg_out
   real *8, intent(out)                                :: lbd_opt, tstep1, tstep2, bias_cvg
   
   !local variables:
   integer                                             :: N, i, j, k, l, count_iter
   real *8, dimension(:), allocatable                  :: coeffs_init, vals
   real *8, dimension(:,:,:), allocatable              :: comp1d_train_init, coeffs_mat_init, comp1d_train
   real *8                                             :: lmax_lasso, lbd_max, lmin_lasso, lbd_min, &
                                                          lbd_step, t0, t1, lmax, lmin
   logical                                             :: flag_cvg
   character(len=4)                                    :: opt_warm, opt_blockCD
                                                          
   ! Random initialization of unknown coefficients:
   N = rank*d*basissize
   allocate(coeffs_init(N))
   call random_seed()
   call random_number(coeffs_init)
   coeffs_init = 2.0D0*coeffs_init - 1.0D0 ! take values in the interval [-1,1]

!    DO i = 1,20
!      print*, 'coeffs=', coeffs_init(i)
!    ENDDO
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
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! STEP 1: Cross-validation to get a suitable value for the penalization parameter !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   ! --- Define a lambda-grid in log10 scale
   allocate(vals(d))
   do j = 1,d
     vals(j) = ABS(dot_product(Xdata_train(j,:),Ydata_train(:)))/N_train
   end do
   lmax_lasso = maxval(vals) !classical lasso upper bound
   lbd_max = log10(lmax_lasso)
   lmin_lasso = epsi*lmax_lasso
   lbd_min = log10(lmin_lasso)
   deallocate(vals)
   lbd_step = (lbd_max-lbd_min)/float(nb_lbd_try-1)
   do i = 1, nb_lbd_try
     lbd_grid(i) = 10**(lbd_min + (i-1)*lbd_step)
!      print*, lbd_grid(i)
   end do
  
   if (iprint == 1) then
     print *, ' Coarse grid for lambda based on lasso:'
     print *, '  ... lbd_min =', 10**lbd_min
     print *, '  ... lbd_max =', 10**lbd_max
     print *, '  ... lbd tries =', nb_lbd_try
     print *, ''
   endif
   
   allocate(comp1d_train(rank,d,N_train))
   opt_warm = 'n' !option for warm starts ('y' or 'n')
   opt_blockCD = 'n' !option for block version ('y' or 'n')
   coeffs_mat(:,:,:) = coeffs_mat_init(:,:,:)
   comp1d_train(:,:,:) = comp1d_train_init(:,:,:)
   call CPU_TIME(t0)
   call crossvalid_kfold_SR(coeffs_mat, Ydata_train, N_train, rank, d, basissize, & 
         varphi_train, nb_lbd_try, lbd_grid, comp1d_train, stock_err_test, & 
          err_test_aver, iprint, kfold, int_resp, opt_blockCD, opt_warm, & 
           opt_resol, lbd_opt, tol_fobj_cv)
   call CPU_TIME(t1)
   tstep1 = t1-t0

   if (iprint == 1) then
     print*, ''
     print*, 'lbd_opt', lbd_opt
   endif
   
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! STEP 2: Production run using lbd_opt                                    !
   ! Returns the set of coefficients of the model structure 'coeffs_mat'     !
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   coeffs_mat(:,:,:) = coeffs_mat_init(:,:,:)
   comp1d_train(:,:,:) = comp1d_train_init(:,:,:)
   call CPU_TIME(t0)
   select case (opt_resol)
     case('cd')
       call coordinate_descent_with_bias(coeffs_mat, Ydata_train, N_train, rank, d, basissize, &
             varphi_train, lbd_opt, comp1d_train, bias_cvg, iprint, count_iter, int_resp, flag_cvg, tol_fobj_run)   
     end select
   call CPU_TIME(t1)
   tstep2 = t1-t0
   
   flag_cvg_out = flag_cvg
   
   deallocate(coeffs_init, coeffs_mat_init, comp1d_train, comp1d_train_init)

   end subroutine srcd_resol

end module srcd