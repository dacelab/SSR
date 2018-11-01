module sigmatune

   use coordDescent
   use utilities
   implicit none

   contains

   !==========================================================================================!
   !
   ! Cross-validation procedure (typically with 2 folds) for sigma-tuning while using RBF
   ! discretizations. This algorithm determines a rough estimate of a suitable value for the 
   ! shape parameter in RBFs.
   subroutine crossvalid_sigmatune(Xdata, Ydata, d, nobs, basissize, rank, opt_basis, kfold, & 
    opt_resol, int_resp, iprint, sigmabest)
   
   !inputs:
   integer, intent(in)                    :: d, nobs, basissize, rank, opt_basis, kfold, iprint
   real *8, dimension(d,nobs), intent(in) :: Xdata
   real *8, dimension(nobs), intent(in)   :: Ydata
   character(len=4), intent(in)           :: opt_resol, int_resp
   
   !outputs:
   real *8, intent(out)                   :: sigmabest
   
   !local variables:
   integer                                :: n1, n2, n3, i, j, k, l, ifold, isig, nsig_try, nslice, nrem, &
                                             N_train, N_test, iprint_srcd, nb_lbd_try, kfold_lbd, cpt, N, count_iter, indmin
   real *8, dimension(:), allocatable     :: sigmatries, Ydata_train, Ydata_test, coeffs_init, ur_test, rmse_aver
   real *8, dimension(:,:), allocatable   :: Xdata_train, Xdata_test, stock_rmse_test
   real *8, dimension(:,:,:), allocatable :: varphi, varphi_train, varphi_test, coeffs_mat, coeffs_mat_init, & 
                                             comp1d_train, comp1d_train_init
   real *8                                :: delta, epsi, tol_fobj_cv, tol_fobj_run, lbd_opt, & 
                                             tstep1, tstep2, bias, lbda, err_test, fobj_test, normD_ur_test, err_sum
   logical                                :: flag_cvg
   
   n1 = 5; n2 = 5; n3 = 5
   nsig_try = n1+n2+n3
   ! Determine the values of sigma to be attempted
   ALLOCATE(sigmatries(nsig_try))
   delta = (0.1D0-0.001D0)/dble(n1-1)
   sigmatries(1:n1) = (/(0.001D0 + dble(j-1)*delta,j=1,n1)/)
   
   delta = (0.9D0-0.2D0)/dble(n2-1)
   sigmatries(n1+1:n1+n2) = (/(0.2D0 + dble(j-1)*delta,j=1,n2)/)
   
   delta = (10.0D0-1.0D0)/dble(n3-1)
   sigmatries(n1+n2+1:n1+n2+n3) = (/(1D0 + dble(j-1)*delta,j=1,n3)/)
   
   !Smaller set of try values:
!    n1 = 5; n2 = 5
!    nsig_try = n1+n2
!    ALLOCATE(sigmatries(nsig_try))
!    delta = (0.9D0-0.5D0)/dble(n1-1)
!    sigmatries(1:n1) = (/(0.5D0 + dble(j-1)*delta,j=1,n1)/)
!    
!    delta = (5.0D0-1.0D0)/dble(n2-1)
!    sigmatries(n1+1:n1+n2) = (/(1.0D0 + dble(j-1)*delta,j=1,n2)/)
   
   do j=1,nsig_try
    print*, 'sigma=', sigmatries(j)
   enddo
   
   nslice = nobs/kfold
   
   !Allocate data structures for TRAINING:
   N_train = (kfold-1)*nslice
   allocate(Xdata_train(d,N_train)) 
   allocate(Ydata_train(N_train))
   allocate(varphi_train(N_train,basissize,d))
   
   !Allocate data structures for TESTING:
   if (MOD(nobs,kfold) == 0) then
     N_test = nslice
   else
     nrem = nobs-kfold*int(nobs/kfold)
     N_test = nslice+nrem
   endif
   allocate(Xdata_test(d,N_test)) 
   allocate(Ydata_test(N_test))
   allocate(varphi_test(N_test,basissize,d))
   
   !Define parameters for SRCD resolution:
   iprint_srcd = -1 !1: display messages / -1: no printing
   tol_fobj_run = 1e-9 
   
   !Other allocations:
   allocate(varphi(nobs,basissize,d))
   allocate(coeffs_mat(rank,d,basissize))
   allocate(comp1d_train_init(rank,d,N_train),comp1d_train(rank,d,N_train))
   allocate(ur_test(N_test))
   allocate(stock_rmse_test(kfold,nsig_try),rmse_aver(nsig_try))
   
   ! Random initialization of unknown coefficients:
   N = rank*d*basissize
   allocate(coeffs_init(N))
   call random_seed()
   call random_number(coeffs_init)
   coeffs_init = 2.0D0*coeffs_init - 1.0D0 ! take values in the interval [-1,1]
   
   ! Recast initial coefficients in a 3d data structure:
   allocate(coeffs_mat_init(rank,d,basissize))
   do l = 1,rank
     do j = 1,d
       do k = 1,basissize
         coeffs_mat_init(l,j,k) = coeffs_init((l-1)*d*basissize+(j-1)*basissize+k) 
       end do
     end do
   end do 
   
   lbda = 0.0D0 !fixed penalization parameter in SRCD resolutions
   
   !==================================!
   !     BEGINNING OF MAIN LOOP
   !==================================!
   
   do ifold = 1,kfold
     
     if (iprint == 1) then
       print*, 'ifold=', ifold
     endif
     
     !Assign data structures for TRAINING:
     Xdata_train(:,:) = -666; Ydata_train(:) = -666
     cpt = 0
     do j = 1,kfold
       if (j /= ifold) then
         cpt = cpt + 1
         Xdata_train(:,(cpt-1)*nslice+1:cpt*nslice) = Xdata(:,(j-1)*nslice+1:j*nslice)
         Ydata_train((cpt-1)*nslice+1:cpt*nslice) = Ydata((j-1)*nslice+1:j*nslice)
       endif
     enddo
     
     !Assign data structures for TESTING:
     Xdata_test(:,:) = -666; Ydata_test(:) = -666
     Xdata_test(:,1:nslice) = Xdata(:,(ifold-1)*nslice+1:ifold*nslice)
     Ydata_test(1:nslice) = Ydata((ifold-1)*nslice+1:ifold*nslice)
     if (MOD(nobs,kfold) /= 0) then 
       !add remaining data et the very end:
       Xdata_test(:,nslice+1:N_test) = Xdata(:,kfold*nslice+1:nobs)
       Ydata_test(nslice+1:N_test) = Ydata(kfold*nslice+1:nobs)
     endif
     
     do isig = 1,nsig_try
     
       if (iprint == 1) then
         print*, 'sigma=', sigmatries(isig)
       endif 
       
       !Evaluate RBF basis functions on the whole dataset:
       CALL eval_basisfunc_work(nobs, basissize, d, opt_basis, Xdata, varphi, sigmatries(isig))
  
       !Assign varphi_train:
       varphi_train(:,:,:) = -666; varphi_test(:,:,:) = -666
       cpt = 0
       do j = 1,kfold
         if (j /= ifold) then
           cpt = cpt + 1
           varphi_train((cpt-1)*nslice+1:cpt*nslice,:,:) = varphi((j-1)*nslice+1:j*nslice,:,:)
         endif
       enddo
       
       !Assign varphi_test:
       varphi_test(1:nslice,:,:) = varphi((ifold-1)*nslice+1:ifold*nslice,:,:)
       if (MOD(nobs,kfold) /= 0) then 
         !add remaining data et the very end:
         varphi_test(nslice+1:N_test,:,:) = varphi(kfold*nslice+1:nobs,:,:)
       endif
      
      !=====================================================================================================!
      ! Apply pathwise SRCD with CV for lambda:
      !
!       call srcd_resol(coeffs_mat, rank, d, basissize, N_train, varphi_train, Xdata_train, & 
!        Ydata_train, nb_lbd_try, epsi, kfold_lbd, iprint, opt_resol, int_resp, flag_cvg, lbd_grid, & 
!         err_test_aver, stock_err_test, lbd_opt, tstep1, tstep2, bias, tol_fobj_cv, tol_fobj_run)
      !
      !UNFORTUNATELY 'N_train' MUST BE DIVISIBLE BY 'kfold_lbd' IN 'srcd_resol'...
      !TO AVOID MODIFYING AGAIN THE PRESENT ROUTINE AND GET A QUICK IMPLEMENTATION WE DIRECTLY USE SRCD
      !WITH LAMBDA=0 (NO SPARSITY)
      !
      
      ! Initializations:
      do l = 1,rank
        do j = 1,d
          do i =1,N_train
            ! Evaluate 1d component functions of the model structure on each training point
            comp1d_train_init(l,j,i) = dot_product(coeffs_mat_init(l,j,:),varphi_train(i,:,j))
          end do
        end do
      end do 
      comp1d_train(:,:,:) = comp1d_train_init(:,:,:)
      coeffs_mat(:,:,:) = coeffs_mat_init(:,:,:)
      
      CALL coordinate_descent_with_bias(coeffs_mat, Ydata_train, N_train, rank, d, basissize, &
             varphi_train, lbda, comp1d_train, bias, iprint_srcd, count_iter, int_resp, flag_cvg, tol_fobj_run)   
     !=====================================================================================================!

     ! Compute validation error on testing dataset:
     if (flag_cvg .eqv. .true.) then     
       CALL err_normD_with_bias(Ydata_test, d, N_test, rank, basissize, varphi_test, &
             coeffs_mat, lbda, err_test, fobj_test, normD_ur_test, ur_test, bias, int_resp)
     end if
     
     if (flag_cvg .eqv. .true.) then
       stock_rmse_test(ifold,isig) = err_test
     elseif (flag_cvg .eqv. .false.) then
       stock_rmse_test(ifold,isig) = -666 
     end if
     
     if (iprint == 1) then
       print*, 'RMSE=', err_test
     endif
     
     enddo !isig
  
   enddo !ifold
   
   !==================================!
   !         END OF MAIN LOOP
   !==================================!
   
   ! Compute the averaged curve (discard the cases when SRCD diverges):
   do isig = 1,nsig_try
     err_sum = 0.0D0
     cpt = 0
     do ifold = 1,kfold
       if ((stock_rmse_test(ifold,isig)) .NE. -666) then
         err_sum = err_sum + stock_rmse_test(ifold,isig)
         cpt = cpt + 1
       end if
     end do
     rmse_aver(isig) = err_sum/cpt
     print*, 'averaged RMSE=', rmse_aver(isig)
   end do
   
   ! Get the optimal sigma value:
   indmin = minloc(rmse_aver,1)
   sigmabest = sigmatries(indmin)
   
   !deallocations:
   deallocate(sigmatries, varphi, Xdata_train, Ydata_train, varphi_train, Xdata_test, Ydata_test, varphi_test, &
    coeffs_mat, coeffs_init, coeffs_mat_init, comp1d_train_init, comp1d_train, ur_test, stock_rmse_test, rmse_aver)
   
   end subroutine  crossvalid_sigmatune

end module sigmatune