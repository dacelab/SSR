module crossValid

   use utilities
   use coordDescent
   
   implicit none

   contains
   
   !==========================================================================================!
   !
   ! Cross-validation procedure for LBFGS solver. Here 1/3 is used for training and 2/3 for
   ! validation (this should be the converse!!!)
   !
!    subroutine crossvalid_LBFGS(xinit, N, Ydata_train, N_train, rank, d, basissize, varphi_train, & 
!     nb_lbd_try, lbd_grid, factr, pgtol, nb_split, stock_err_test, err_test_aver)
!    
!    !inputs:
!    integer, intent(in)                                 :: N, N_train, rank, d, basissize, nb_lbd_try
!    real *8, dimension(nb_lbd_try), intent(in)          :: lbd_grid
!    real *8, dimension(N_train), intent(in)             :: Ydata_train
!    real *8, dimension(N_train,basissize,d), intent(in) :: varphi_train
!    real *8, intent(in)                                 :: factr, pgtol
!    integer, intent(in)                                 :: nb_split
!    real *8, dimension(N), intent(in)                   :: xinit
! 
!    !outputs:
!    real *8, dimension(nb_split,nb_lbd_try), intent(out) :: stock_err_test
!    real *8, dimension(nb_lbd_try), intent(out)          :: err_test_aver
!    
!    !local variables:
!    integer                                             :: isplit, ilbd, N_slice, i, j, l, k, cpt, N_test, optf
!    real *8                                             :: lbd_try, fobj, err_test, fobj_test, normD_ur_test, &
!                                                           normD_u_test, relativ_err_test
!    real *8, dimension(:), allocatable                  :: Ydata_train_wk, Ydata_test_wk, ur_test
!    real *8, dimension(:,:,:), allocatable              :: varphi_train_wk, varphi_test_wk, x_mat
!    real *8, dimension(N)                               :: x
!    
!    !Make sure N_train is a multiple of the number of splits (3):
!    if (mod(N_train,3) /= 0) then
!       print*, 'cross validation (LBFGS): N_train is not divisible by 3!!!'
!       STOP
!    else
!       N_slice = N_train/3
!       N_test = 2*N_slice
!       allocate(Ydata_train_wk(N_slice), varphi_train_wk(N_slice,basissize,d))
!       allocate(Ydata_test_wk(N_test), varphi_test_wk(N_test,basissize,d))
!    end if
!    
!    allocate(x_mat(rank,d,basissize), ur_test(N_test))
!    optf = 0 ! 0/1: don't print/print objective function values in LBFGS
!    
!    do isplit = 1,nb_split
!    
!       print*, 'isplit=', isplit
!       
!       ! Assign working data structure for training (one third of training dataset):
!       Ydata_train_wk = Ydata_train((isplit-1)*N_slice+1:isplit*N_slice)
!       varphi_train_wk = varphi_train((isplit-1)*N_slice+1:isplit*N_slice,:,:)
!       ! Assign working data structure for testing (two remaining thirds):
!       cpt = 0
!       do i = 1,nb_split
!         if (i /= isplit) then
!           cpt = cpt + 1
!           Ydata_test_wk((cpt-1)*N_slice+1:cpt*N_slice) = Ydata_train((i-1)*N_slice+1:i*N_slice)
!           varphi_test_wk((cpt-1)*N_slice+1:cpt*N_slice,:,:) = varphi_train((i-1)*N_slice+1:i*N_slice,:,:)
!         end if
!       end do
!       
!       normD_u_test = sqrt((sum(Ydata_test_wk**2))/N_test)
!       
!       do ilbd = 1,nb_lbd_try
!       
! 	  lbd_try = lbd_grid(ilbd)
! 	  print *, 'lbda=', lbd_try
! 	  x = xinit
!           call LBFGSDriver(x, N, Ydata_train_wk, N_slice, rank, d, basissize, varphi_train_wk, & 
!            lbd_try, fobj, factr, pgtol, optf)
! 
! 	  ! Recast coefficients in a convenient 3d data structure:
!           do l = 1,rank
!             do j = 1,d
!               do k = 1,basissize
!                 x_mat(l,j,k) = x((l-1)*d*basissize+(j-1)*basissize+k) 
!               end do
!             end do
!           end do  
! 	  
! 	  ! Compute testing error:
!           call err_normD(Ydata_test_wk, d, N_test, rank, basissize, varphi_test_wk, x_mat, & 
!            lbd_try, err_test, fobj_test, normD_ur_test, ur_test)
!  
!           relativ_err_test = err_test/normD_u_test
! 	  print*, 'relative err=', relativ_err_test
! 	  
! 	  stock_err_test(isplit,ilbd) = relativ_err_test
! 	  
!       end do !ilbd
!    
!    end do !isplit
!    
!    ! Compute the averaged curve:
!    do ilbd = 1,nb_lbd_try
!       err_test_aver(ilbd) = sum(stock_err_test(:,ilbd))/nb_split
!    end do
!    
!    ! Save testing errors/lambda-grid:
!    open (20, FILE='./Data/errors_LBFGS_cv.dat', STATUS='UNKNOWN', FORM='FORMATTED') 
!    do isplit = 1,nb_split
!      do ilbd = 1,nb_lbd_try
!        write (20, *) stock_err_test(isplit,ilbd)
!      end do
!    end do
!    close (20)
!    
!    open (21, FILE='./Data/errorAver_LBFGS_cv.dat', STATUS='UNKNOWN', FORM='FORMATTED') 
!    do ilbd = 1,nb_lbd_try
!      write (21, *) err_test_aver(ilbd)
!    end do
!    close (21)
!    
!    open (22, FILE='./Data/lbdgrid_LBFGS_cv.dat', STATUS='UNKNOWN', FORM='FORMATTED') 
!    do ilbd = 1,nb_lbd_try
!      write (22, *) lbd_grid(ilbd)
!    end do
!    close (22)
!    
!    !deallocations:
!    deallocate(Ydata_train_wk, varphi_train_wk, Ydata_test_wk, varphi_test_wk, x_mat, ur_test)
!    
!    end subroutine crossvalid_LBFGS
   !==========================================================================================!

   !==========================================================================================!
   !
   ! Cross-validation procedure for LBFGS solver using k-folds (typically, k=5 or k=10).
   ! (k-1) folds are used for training and one fold for testing.
   ! A pathwise algorithm with restarts is used.
   !
   subroutine crossvalid_kfold_LBFGS(xinit, N, Ydata_train, N_train, rank, d, basissize, varphi_train, &
    nb_lbd_try, lbd_grid, factr, pgtol, stock_err_test, err_test_aver, kfold, int_resp, opt_warm, lbd_opt)
   
   !inputs:
   integer, intent(in)                                 :: N, N_train, rank, d, basissize, nb_lbd_try, &
                                                          kfold
   real *8, dimension(nb_lbd_try), intent(in)          :: lbd_grid
   real *8, dimension(N_train), intent(in)             :: Ydata_train
   real *8, dimension(N_train,basissize,d), intent(in) :: varphi_train
   real *8, intent(in)                                 :: factr, pgtol
   real *8, dimension(N), intent(in)                   :: xinit
   character(len=10), intent(in)                       :: int_resp
   character(len=4), intent(in)                        :: opt_warm

   !outputs:
   real *8, dimension(kfold,nb_lbd_try), intent(out)   :: stock_err_test
   real *8, dimension(nb_lbd_try), intent(out)         :: err_test_aver
   real *8, intent(out)                                :: lbd_opt
   
   !local variables:
   integer                                             :: ifold, ilbd, N_slice, i, j, l, k, cpt, N_testCV, & 
                                                          N_trainCV, optf, indmin
   real *8                                             :: lbd_try, fobj, err_test, fobj_test, normD_ur_test, &
                                                          normD_u_test, relativ_err_test
   real *8, dimension(:), allocatable                  :: Ydata_train_wk, Ydata_test_wk, ur_test
   real *8, dimension(:,:,:), allocatable              :: varphi_train_wk, varphi_test_wk, x_mat
   real *8, dimension(N)                               :: x
   
   !Make sure N_train is a multiple of kfold:
   if (mod(N_train,kfold) /= 0) then
      print*, 'cross validation (LBFGS): N_train is not divisible by kfold!!!'
      STOP
   else
      N_slice = N_train/kfold
      N_trainCV = (kfold-1)*N_slice
      N_testCV = N_slice
      allocate(Ydata_train_wk(N_trainCV), varphi_train_wk(N_trainCV,basissize,d))
      allocate(Ydata_test_wk(N_testCV), varphi_test_wk(N_testCV,basissize,d))
   end if
   
   allocate(x_mat(rank,d,basissize), ur_test(N_testCV))
   optf = 0 ! 0/1: don't print/print objective function values in LBFGS
   
!    print*, 'N_slice=', N_slice
!    print*, 'N_trainCV=', N_trainCV
!    print*, 'N_testCV=', N_testCV
   
   do ifold = 1,kfold
   
      print*, 'ifold=', ifold
      
      ! Assign working data structure for training (using all folds except one):
      cpt = 0
      do i = 1,kfold
        if (i /= ifold) then
          cpt = cpt + 1
          Ydata_train_wk((cpt-1)*N_slice+1:cpt*N_slice) = Ydata_train((i-1)*N_slice+1:i*N_slice)
          varphi_train_wk((cpt-1)*N_slice+1:cpt*N_slice,:,:) = varphi_train((i-1)*N_slice+1:i*N_slice,:,:)
        end if
      end do
      
      ! Assign working data structure for testing (using the remaining fold):
      Ydata_test_wk = Ydata_train((ifold-1)*N_slice+1:ifold*N_slice)
      varphi_test_wk = varphi_train((ifold-1)*N_slice+1:ifold*N_slice,:,:)
            
      normD_u_test = sqrt((sum(Ydata_test_wk**2))/N_testCV)
      x = xinit
      
      do ilbd = 1,nb_lbd_try
      
          lbd_try = lbd_grid(nb_lbd_try-ilbd+1)
          print *, 'lbda=', lbd_try
          
          call LBFGSDriver(x, N, Ydata_train_wk, N_trainCV, rank, d, basissize, varphi_train_wk, &
           lbd_try, fobj, factr, pgtol, optf)
                     
	  ! Recast coefficients in a convenient 3d data structure:
          do l = 1,rank
            do j = 1,d
              do k = 1,basissize
                x_mat(l,j,k) = x((l-1)*d*basissize+(j-1)*basissize+k) 
              end do
            end do
          end do  

          ! Compute testing error:
          call err_normD(Ydata_test_wk, d, N_testCV, rank, basissize, varphi_test_wk, x_mat, & 
           lbd_try, err_test, fobj_test, normD_ur_test, ur_test, int_resp)
           
          relativ_err_test = err_test/normD_u_test
          print*, 'relative err=', relativ_err_test
          
          stock_err_test(ifold,ilbd) = relativ_err_test

          if (opt_warm == 'n') then
            !Use the same (random) initialization
            x = xinit
          elseif (opt_warm == 'y') then   
            !Use the output solution 'x' for the next lambda-try
          end if
  
      end do !ilbd
   
   end do !kfold
      
   ! Compute the averaged curve:
   do ilbd = 1,nb_lbd_try
      err_test_aver(ilbd) = sum(stock_err_test(:,ilbd))/kfold
   end do
   
   ! Get the optimal lambda parameter:
   indmin = minloc(err_test_aver,1)
   lbd_opt = lbd_grid(nb_lbd_try-indmin+1)    
       
   !deallocations:
   deallocate(Ydata_train_wk, varphi_train_wk, Ydata_test_wk, varphi_test_wk, x_mat, ur_test)
      
   end subroutine crossvalid_kfold_LBFGS
   !==========================================================================================!

   !==========================================================================================!
   !
   ! Cross-validation procedure for Coordinate Descent algorithm. Here 1/3 is used for training 
   ! and 2/3 for validation (this should be the converse!!!)
   !
!    subroutine crossvalid_CD(coeffs_init, Ydata_train, N_train, rank, d, basissize, varphi_train, &
!     opt_bias, nb_lbd_try, lbd_grid, comp1d_train, nb_split, stock_err_test, err_test_aver, iprint)
!    
!    !inputs:
!    integer, intent(in)                                  :: N_train, rank, d, basissize, nb_lbd_try, &
!                                                            nb_split, iprint
!    character(len=4)                                     :: opt_bias                                                          
!    real *8, dimension(nb_lbd_try), intent(in)           :: lbd_grid
!    real *8, dimension(N_train), intent(in)              :: Ydata_train
!    real *8, dimension(N_train,basissize,d), intent(in)  :: varphi_train
!    real *8, dimension(rank,d,N_train), intent(in)       :: comp1d_train   
! 
!    !outputs:
!    real *8, dimension(nb_split,nb_lbd_try), intent(out) :: stock_err_test
!    real *8, dimension(nb_lbd_try), intent(out)          :: err_test_aver
!    
!    !local variables:
!    integer                                             :: N_slice, N_test, isplit, ilbd, i, cpt, count_iter
!    real *8, dimension(:), allocatable                  :: Ydata_train_wk, Ydata_test_wk, ur_test
!    real *8, dimension(:,:,:), allocatable              :: varphi_train_wk, varphi_test_wk, comp1d_train_wk, &
!                                                           comp1d_train_copy
!    real *8                                             :: normD_u_test, lbd_try, bias, relativ_err_test, & 
!                                                           err_test, fobj_test, normD_ur_test
!    real *8, dimension(rank,d,basissize)                :: coeffs
!    
!    !Make sure N_train is a multiple of the number of splits (3):
!    if (mod(N_train,3) /= 0) then
!       print*, 'cross validation (CD): N_train is not divisible by 3!!!'
!       STOP
!    else
!       N_slice = N_train/3
!       N_test = 2*N_slice
!       allocate(Ydata_train_wk(N_slice))
!       allocate(varphi_train_wk(N_slice,basissize,d))
!       allocate(comp1d_train_wk(rank,d,N_slice), comp1d_train_copy(rank,d,N_slice))
!       allocate(Ydata_test_wk(N_test))
!       allocate(varphi_test_wk(N_test,basissize,d))
!       allocate(ur_test(N_test))
!    end if
!    
!    do isplit = 1,nb_split
!    
!       print*, 'isplit=', isplit
!       
!       ! Assign working data structure for training (one third of training dataset):
!       Ydata_train_wk = Ydata_train((isplit-1)*N_slice+1:isplit*N_slice)
!       varphi_train_wk = varphi_train((isplit-1)*N_slice+1:isplit*N_slice,:,:)
!       comp1d_train_wk = comp1d_train(:,:,(isplit-1)*N_slice+1:isplit*N_slice)
!       ! Assign working data structure for testing (two remaining thirds):
!       cpt = 0
!       do i = 1,nb_split
!         if (i /= isplit) then
!           cpt = cpt + 1
!           Ydata_test_wk((cpt-1)*N_slice+1:cpt*N_slice) = Ydata_train((i-1)*N_slice+1:i*N_slice)
!           varphi_test_wk((cpt-1)*N_slice+1:cpt*N_slice,:,:) = varphi_train((i-1)*N_slice+1:i*N_slice,:,:)
!         end if
!       end do
!       
!       normD_u_test = sqrt((sum(Ydata_test_wk**2))/N_test)
!    
!       do ilbd = 1,nb_lbd_try
!       
! 	  lbd_try = lbd_grid(ilbd)
! 	  print *, 'lbda=', lbd_try
!           
!           coeffs = coeffs_init
!           comp1d_train_copy = comp1d_train_wk
!           
!           if (opt_bias == 'n') then 
!              ! call CD solver:
!              call coordinate_descent(coeffs, Ydata_train_wk, N_slice, rank, d, basissize, &
!               varphi_train_wk, lbd_try, comp1d_train_copy, iprint)
!              ! Compute testing error:
!              call err_normD(Ydata_test_wk, d, N_test, rank, basissize, varphi_test_wk, &
!               coeffs, lbd_try, err_test, fobj_test, normD_ur_test, ur_test)
!           elseif (opt_bias == 'y') then
!              ! call CD solver:
!              call coordinate_descent_with_bias(coeffs, Ydata_train_wk, N_slice, rank, d, basissize, &
!               varphi_train_wk, lbd_try, comp1d_train_copy, bias, iprint, count_iter)
!              ! Compute testing error:
!              call err_normD_with_bias(Ydata_test_wk, d, N_test, rank, basissize, varphi_test_wk, &
!               coeffs, lbd_try, err_test, fobj_test, normD_ur_test, ur_test, bias)
!           end if
!           
!           relativ_err_test = err_test/normD_u_test
! 	  print*, 'relative err=', relativ_err_test
! 
! 	  stock_err_test(isplit,ilbd) = relativ_err_test
!    
!       end do !ilbd
!    
!    end do !isplit
!    
!    ! Compute the averaged curve:
!    do ilbd = 1,nb_lbd_try
!       err_test_aver(ilbd) = sum(stock_err_test(:,ilbd))/nb_split
!    end do
!    
!    ! Save testing errors/lambda-grid:
!    open (20, FILE='./Data/errors_CD_cv.dat', STATUS='UNKNOWN', FORM='FORMATTED') 
!    do isplit = 1,nb_split
!      do ilbd = 1,nb_lbd_try
!        write (20, *) stock_err_test(isplit,ilbd)
!      end do
!    end do
!    close (20)
!    
!    open (21, FILE='./Data/errorAver_CD_cv.dat', STATUS='UNKNOWN', FORM='FORMATTED') 
!    do ilbd = 1,nb_lbd_try
!      write (21, *) err_test_aver(ilbd)
!    end do
!    close (21)
!    
!    open (22, FILE='./Data/lbdgrid_CD_cv.dat', STATUS='UNKNOWN', FORM='FORMATTED') 
!    do ilbd = 1,nb_lbd_try
!      write (22, *) lbd_grid(ilbd)
!    end do
!    close (22)
!    
!    !deallocations:
!    deallocate(Ydata_train_wk, varphi_train_wk, Ydata_test_wk, varphi_test_wk, comp1d_train_wk, ur_test)
!    
!    end subroutine crossvalid_CD
   !==========================================================================================! 
   
   !==========================================================================================!
   !
   ! Cross-validation procedure using k-folds (typically, k=5 or k=10) in our separated rank
   ! (SR) algorithm to determine suitable value for penalization parameter.
   !  * (k-1) folds are used for training and one fold is used for testing.
   !  * A pathwise algorithm with or without warm starts is used.
   !  * A non-block or block version of SR can be used.
   !     -- non-block version based on CD updates.
   !     -- block version based on BFGS or GLMNET resolution.
   subroutine crossvalid_kfold_SR(coeffs_in, Ydata_train, N_train, rank, d, basissize, & 
    varphi_train, nb_lbd_try, lbd_grid, comp1d_train, stock_err_test, err_test_aver, & 
     iprint, kfold, int_resp, opt_blockCD, opt_warm, opt_resol, lbd_opt, tol_fobj_cv)
   
   !inputs:
   integer, intent(in)                                  :: N_train, rank, d, basissize, nb_lbd_try, &
                                                           kfold, iprint
   character(len=4)                                     :: opt_blockCD, opt_warm, opt_resol                                             
   real *8, dimension(nb_lbd_try), intent(in)           :: lbd_grid
   real *8, dimension(N_train), intent(in)              :: Ydata_train
   real *8, dimension(N_train,basissize,d), intent(in)  :: varphi_train
   real *8, dimension(rank,d,N_train), intent(in)       :: comp1d_train
   real *8, dimension(rank,d,basissize), intent(in)     :: coeffs_in
   character(len=4), intent(in)                         :: int_resp
   real *8, intent(in)                                  :: tol_fobj_cv
 
   !outputs:
   real *8, dimension(kfold,nb_lbd_try), intent(out)    :: stock_err_test
   real *8, dimension(nb_lbd_try), intent(out)          :: err_test_aver
   real *8, intent(out)                                 :: lbd_opt
   
   !local variables:
   integer                                              :: N_slice, N_testCV, N_trainCV, ifold, ilbd, i, j, l, &
                                                           cpt, count_iter, indmin
   real *8, dimension(:), allocatable                   :: Ydata_train_wk, Ydata_test_wk, ur_test
   real *8, dimension(:,:,:), allocatable               :: varphi_train_wk, varphi_test_wk, comp1d_train_wk, &
                                                           comp1d_train_copy
   real *8                                              :: normD_u_test, lbd_try, bias, relativ_err_test, & 
                                                           err_test, fobj_test, normD_ur_test, err_sum
   real *8, dimension(rank,d,basissize)                 :: coeffs
   logical                                              :: flag_cvg   
   
   !Make sure N_train is a multiple of kfold:
   if (mod(N_train,kfold) /= 0) then
      print*, 'cross validation (CD): N_train is not divisible by kfold!!!'
      STOP
   else
      N_slice = N_train/kfold
      N_trainCV = (kfold-1)*N_slice
      N_testCV = N_slice
      allocate(Ydata_train_wk(N_trainCV))
      allocate(varphi_train_wk(N_trainCV,basissize,d))
      allocate(comp1d_train_wk(rank,d,N_trainCV))
      allocate(comp1d_train_copy(rank,d,N_trainCV))
      allocate(Ydata_test_wk(N_testCV))
      allocate(varphi_test_wk(N_testCV,basissize,d))
      allocate(ur_test(N_testCV))
   end if
   
   do ifold = 1,kfold
      
      if (iprint == 1) then
        print*, 'ifold=', ifold
      endif 
      
      ! Assign working data structure for training (using all the folds except one):
      Ydata_train_wk(:) = 0.0D0
      varphi_train_wk(:,:,:) = 0.0D0
      comp1d_train_wk(:,:,:) = 0.0D0
      
      cpt = 0
      do i = 1,kfold
        if (i /= ifold) then
          cpt = cpt + 1
!           print*,'i=', i
!           print*,'istart (train)=', (i-1)*N_slice+1
          Ydata_train_wk((cpt-1)*N_slice+1:cpt*N_slice) = Ydata_train((i-1)*N_slice+1:i*N_slice)
          varphi_train_wk((cpt-1)*N_slice+1:cpt*N_slice,:,:) = varphi_train((i-1)*N_slice+1:i*N_slice,:,:)
          comp1d_train_wk(:,:,(cpt-1)*N_slice+1:cpt*N_slice) = comp1d_train(:,:,(i-1)*N_slice+1:i*N_slice)
        end if
      end do
      
      ! Assign working data structure for testing (using the remaining fold):
!       print*,'istart (test)=',(ifold-1)*N_slice+1
      Ydata_test_wk(:) = 0.0D0
      varphi_test_wk(:,:,:) = 0.0D0
      Ydata_test_wk = Ydata_train((ifold-1)*N_slice+1:ifold*N_slice)
      varphi_test_wk = varphi_train((ifold-1)*N_slice+1:ifold*N_slice,:,:)
      
      normD_u_test = sqrt((sum(Ydata_test_wk**2))/N_testCV)
      
      coeffs(:,:,:) = coeffs_in(:,:,:)
      comp1d_train_copy(:,:,:) = comp1d_train_wk(:,:,:)
          
      do ilbd = 1,nb_lbd_try
 
          lbd_try = lbd_grid(nb_lbd_try-ilbd+1)
          if (iprint == 1) then
            print *, 'lbda=', lbd_try
          endif 
          
          if (opt_blockCD == 'y') then
            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! BLOCK VERSION ALGORITHM !
            !!!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            select case (opt_resol)
              case('bfgs')
                ! No bias term in this formulation
                call BSR_subpbs_bfgs(coeffs, Ydata_train_wk, N_trainCV, rank, d, basissize, & 
                      varphi_train_wk, lbd_try, comp1d_train_copy, iprint, count_iter, int_resp, & 
                       flag_cvg, tol_fobj_cv)
              case('glmn')
                call BSR_subpbs_glmnet(coeffs, Ydata_train_wk, N_trainCV, rank, d, basissize, & 
                      varphi_train_wk, lbd_try, comp1d_train_copy, iprint, count_iter, int_resp, &
                       flag_cvg, tol_fobj_cv)    
            end select    
            ! Compute testing error (no bias term):
            if (flag_cvg .eqv. .true.) then
              call err_normD(Ydata_test_wk, d, N_testCV, rank, basissize, varphi_test_wk, &
                    coeffs, lbd_try, err_test, fobj_test, normD_ur_test, ur_test, int_resp)
            end if
            !
          elseif (opt_blockCD == 'n') then
            !
            !!!!!!!!!!!!!!!!!!!!!!!!!!
            ! NON-BLOCK SR ALGORITHM !
            !!!!!!!!!!!!!!!!!!!!!!!!!!
            !
            select case (opt_resol)
              case('cd')
                !SRCD algorithm (with bias term)
                call coordinate_descent_with_bias(coeffs, Ydata_train_wk, N_trainCV, rank, d, basissize, &
                       varphi_train_wk, lbd_try, comp1d_train_copy, bias, iprint, count_iter, int_resp, flag_cvg, tol_fobj_cv)
                ! Compute testing error:
                if (flag_cvg .eqv. .true.) then
                  call err_normD_with_bias(Ydata_test_wk, d, N_testCV, rank, basissize, varphi_test_wk, &
                        coeffs, lbd_try, err_test, fobj_test, normD_ur_test, ur_test, bias, int_resp)
                end if
            end select    
            !   
          end if
          
          if (flag_cvg .eqv. .true.) then
            relativ_err_test = err_test/normD_u_test
            if (iprint == 1) then
              print*, 'validation error=', relativ_err_test
              print*, 'flag_cvg=', flag_cvg
!             print*, 'nb of iters=', count_iter
            endif
          elseif (flag_cvg .eqv. .false.) then
            relativ_err_test = -666 
          end if

	  stock_err_test(ifold,ilbd) = relativ_err_test
	  
	  if (opt_warm == 'n') then
	    !Use the same initialization for each try of lambda:
	    coeffs(:,:,:) = coeffs_in(:,:,:)
            comp1d_train_copy(:,:,:) = comp1d_train_wk(:,:,:)
	  elseif (opt_warm == 'y') then
            !Use outputs 'coeffs' and 'comp1d_train_copy' as initial guesses for next try of lambda
          end if
   
      end do !ilbd
   
   end do !ifold
   
!    ! Compute the averaged curve (if SRCD converges for any lambda):
!    do ilbd = 1,nb_lbd_try
!       err_test_aver(ilbd) = sum(stock_err_test(:,ilbd))/kfold
!    end do

   ! Compute the averaged curve (discard the cases when SRCD diverges):
   do ilbd = 1,nb_lbd_try 
     err_sum = 0.0D0
     cpt = 0
     do ifold = 1,kfold
       if ((stock_err_test(ifold,ilbd)) .NE. -666) then
         err_sum = err_sum + stock_err_test(ifold,ilbd)
         cpt = cpt + 1
       end if
     end do
     err_test_aver(ilbd) = err_sum/cpt
   end do
   
   ! Get the optimal lambda parameter:
   indmin = minloc(err_test_aver,1)
   lbd_opt = lbd_grid(nb_lbd_try-indmin+1)
   
   !deallocations:
   deallocate(Ydata_train_wk, varphi_train_wk, Ydata_test_wk, varphi_test_wk, comp1d_train_wk, & 
    comp1d_train_copy, ur_test)  
   
   end subroutine crossvalid_kfold_SR
   !==========================================================================================! 
   
   
   
   !==========================================================================================! 
   ! Occasionally block_SRCD() diverges. For example, consider Friedman1 dataset with
   ! N_train=300, rank=4, 3rd-order Legendre polynomials. When running crossvalid_kfold_CD() with 
   ! kfold=5, opt_warm = 'n', nb_lbd_try_c = 20 and epsi = 1e-5 block_SRCD() blows up for ifold=2
   ! (while everything is working perfectly fine for the other folds!). Nothing is apparently wrong
   ! in the dataset, maybe linked to the convergence of SRBCD to a Nash point?
   !
   subroutine debug_CV(N_train, kfold, lbd, rank, d, basissize, iprint, Ydata_train, varphi_train, &
    comp1d_train, coeffs_mat, int_resp, tol_fobj_cv)
   
   !inputs:
   integer, intent(in)                                  :: N_train, kfold, rank, d, basissize, iprint
   real *8, intent(in)                                  :: lbd, tol_fobj_cv
   real *8, dimension(N_train), intent(in)              :: Ydata_train
   real *8, dimension(N_train,basissize,d), intent(in)  :: varphi_train
   real *8, dimension(rank,d,N_train), intent(in)       :: comp1d_train
   real *8, dimension(rank,d,basissize), intent(inout)  :: coeffs_mat
   character(len=10)                                    :: int_resp
  
   !local vars:
   integer                                :: N_slice, N_trainCV, ifold, cpt, i, j, count_iter
   real *8, dimension(:), allocatable     :: Ydata_train_wk
   real *8, dimension(:,:,:), allocatable :: varphi_train_wk, comp1d_train_wk
   real *8                                :: bias
   logical                                :: flag_cvg
   
   N_slice = N_train/kfold
   N_trainCV = (kfold-1)*N_slice
   
   print*, 'N_slice=', N_slice
   print*, 'N_trainCV=', N_trainCV
   
   allocate(Ydata_train_wk(N_trainCV))
   allocate(varphi_train_wk(N_trainCV,basissize,d))
   allocate(comp1d_train_wk(rank,d,N_trainCV))
   
   ifold = 1
   cpt = 0
   do i = 1,kfold
      if (i /= ifold) then
          cpt = cpt + 1
          Ydata_train_wk((cpt-1)*N_slice+1:cpt*N_slice) = Ydata_train((i-1)*N_slice+1:i*N_slice)
          varphi_train_wk((cpt-1)*N_slice+1:cpt*N_slice,:,:) = varphi_train((i-1)*N_slice+1:i*N_slice,:,:)
          comp1d_train_wk(:,:,(cpt-1)*N_slice+1:cpt*N_slice) = comp1d_train(:,:,(i-1)*N_slice+1:i*N_slice)
      end if
   end do
   
!    print*,'max(|Ydata_train_wk|)=', maxval(abs(Ydata_train_wk))
!    print*,'max(|varphi_train_wk|)=', maxval(abs(varphi_train_wk))
!    print*,'max(|comp1d_train_wk|)=', maxval(abs(comp1d_train_wk))
!    print *, 'paused, type [enter] to continue'
!    read (*,*) 
  
   ! Block SRCD:
!    call block_SRCD(coeffs_mat, Ydata_train_wk, N_trainCV, rank, d, basissize, &
!                     varphi_train_wk, lbd, comp1d_train_wk, bias, iprint, count_iter, int_resp, flag_cvg)
                    
   ! Block SRCD using LBFGS for each supbproblem:               
   call BSR_subpbs_bfgs(coeffs_mat, Ydata_train_wk, N_trainCV, rank, d, basissize, & 
                         varphi_train_wk, lbd, comp1d_train_wk, iprint, count_iter, int_resp, flag_cvg, tol_fobj_cv)
                         
   ! SRCD:                          
!    call coordinate_descent_with_bias(coeffs_mat, Ydata_train_wk, N_trainCV, rank, d, basissize, &
!          varphi_train_wk, lbd, comp1d_train_wk, bias, iprint, count_iter, int_resp)    
   
   deallocate(Ydata_train_wk, varphi_train_wk, comp1d_train_wk)
   
   end subroutine debug_CV

end module crossValid