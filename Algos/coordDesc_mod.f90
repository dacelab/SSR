module coordDescent

   use utilities
   implicit none

   contains
   
   !==========================================================================================!
   !
   subroutine coordinate_descent(coeffs, Ydata_train, N_train, rank, d, basissize, & 
    varphi_train, lbda, comp1d_train, iprint, int_resp, flag_cvg)
   
   !inputs:
   integer, intent(in)                                 :: N_train, rank, d, basissize, iprint
   real *8, intent(in)                                 :: lbda
   real *8, dimension(N_train), intent(in)             :: Ydata_train
   real *8, dimension(N_train,basissize,d), intent(in) :: varphi_train
   character(len=10), intent(in)                       :: int_resp
   
   !outputs:
   real *8, dimension(rank,d,basissize), intent(inout) :: coeffs
   real *8, dimension(rank,d,N_train), intent(inout)   :: comp1d_train
   logical, intent(out)                                :: flag_cvg

   !local variables:
   real *8                                             :: scaling, Ytmp1, Ytmp2, prod_tmp, betahat, &
                                                          coeffnew, err, err_old, diff_err, tol_diff, &
                                                          tol_err, tol_fobj, normD_u, fobj, fobj_old, & 
                                                          diff_fobj, tmp, normD_ur, denom
   real *8, dimension(N_train)                         :: partproduct, Xinputs, Yresp, ur_train
   integer                                             :: i, j, k, l, m, jj, ll, count_iter

   normD_u = sqrt((sum(Ydata_train**2))/N_train)
   
   tol_diff = 1e-2 !tolerance for saturation of error
   tol_err = 1e-2  !tolerance for convergence of relative error
   tol_fobj = 1e-9 !tolerance for saturation of objective function
   err_old = 0.0D0
   
   fobj_old = 0.0D0
   
   count_iter = 0
   
   flag_cvg = .false. !test for saturation of the error norm ||u-ur||_D
   
!    do while (count_iter < 100) !FOR DEBUG
   
   wloop: do while (flag_cvg .eqv. .false.)
   
     !------------Big loop starts here------------
     do l = 1,rank
!        print *,'l=', l
       do j =1,d
!        print *,'j=', j
         do k = 1,basissize
           !
           ! Compute partial products:
           do i = 1,N_train
	     partproduct(i) = 1.0D0
	     do jj = 1,d
	       if (jj /= j) then
	         partproduct(i) = partproduct(i) * comp1d_train(l,jj,i)
	       end if
	     end do !jj
           end do !i
           !
           ! Compute "input data points" to be fitted:
           do i = 1,N_train
             Xinputs(i) = partproduct(i) * varphi_train(i,k,j) 
           end do !i
           !
           ! Compute "modified response":
           do i = 1,N_train
             Ytmp1 = 0.0D0
             do ll = 1,rank
               if (ll /= l) then
                 Ytmp1 = Ytmp1 + product(comp1d_train(ll,1:d,i))
               end if
             end do !ll
             !
             Ytmp2 = 0.0D0
             do m = 1,basissize
               if (m /= k) then
                 Ytmp2 = Ytmp2 + coeffs(l,j,m)*varphi_train(i,m,j)
               end if
             end do !m
             prod_tmp = 1.0D0
             do jj = 1,d
               if (jj /= j) then
                 prod_tmp = prod_tmp * comp1d_train(l,jj,i)
               end if
             end do !jj
             Ytmp2 = prod_tmp * Ytmp2
             !
             Yresp(i) = Ydata_train(i) - Ytmp1 - Ytmp2
           end do !i
           !
           ! Compute first input variable of soft-thresholding operator:
           betahat = (dot_product(Xinputs,Yresp))/N_train
           !
           ! Update coefficient by applying the soft-thresholding operator:
           call soft_threshold(betahat, lbda, coeffnew)
           denom = (dot_product(Xinputs,Xinputs))/N_train
           if (abs(denom) < 1e-16) then
             !'partproduct' is too small (can happen when 'lbda' is too large)
             if (iprint == 1) then
               print*, "Impossible to update coefficient (denom too small)!!!"
             end if  
             return
           else
             coeffs(l,j,k) = coeffnew/denom
           endif
           !
           count_iter = count_iter + 1
           !
         end do !k
         !
         ! Update 1d component function with new coefficients: 
         do i = 1,N_train
           comp1d_train(l,j,i) = dot_product(coeffs(l,j,1:basissize),varphi_train(i,1:basissize,j)) 
         end do !i
         !
         ! Compute error on the training dataset:
         call err_normD(Ydata_train, d, N_train, rank, basissize, varphi_train, coeffs, lbda, err, &
                        fobj, normD_ur, ur_train, int_resp) 
         
!          print*, '||u-ur||_D/||u||_D=', err/normD_u
!          diff_err = abs(err-err_old)
!          err_old = err
!          print *,'fobj=', fobj
         diff_fobj = fobj-fobj_old
         fobj_old = fobj
         !
         if (abs(diff_fobj) < tol_fobj*(normD_u**2)) then ! test on the difference between two successive objective values
!          if (diff_err < tol_diff*normD_u) then ! test on the difference between two successive errors
!           if (err < tol_err*normD_u) then ! WARNING: even if ||u-ur||_D always decreases it doesn't 
!                                           ! necessarily converge to 0 (especially when the rank is small)
	    flag_cvg = .true.
	    if (iprint == 1) then
              print*, 'Relative training error=', err/normD_u
	    end if
	    return
	 elseif (diff_fobj > 0.0D0 .and. count_iter > basissize) then
! 	    !Objective function is increasing (after first pass):
! 	    print*, 'Objective function is increasing'
! 	    test = .true.
! 	    print*, 'Relative training error=', err/normD_u
! 	    return
         end if
         if (fobj > 1e6) then
           print *,'fobj=', fobj 
           print*, 'SRCD: objective function diverges!!!' 
           exit wloop
         endif
         ! 
       end do !j  
     end do !l
   
     !-----------------End of Big loop-----------------
   
   end do wloop !while
   
   end subroutine coordinate_descent
   !==========================================================================================!

   !==========================================================================================!
   !
   subroutine coordinate_descent_with_bias(coeffs, Ydata_train, N_train, rank, d, basissize, & 
    varphi_train, lbda, comp1d_train, bias_cvg, iprint, count_iter, int_resp, flag_cvg, tol_fobj)
   
   !inputs:
   integer, intent(in)                                 :: N_train, rank, d, basissize, iprint
   real *8, intent(in)                                 :: lbda
   real *8, dimension(N_train), intent(in)             :: Ydata_train
   real *8, dimension(N_train,basissize,d), intent(in) :: varphi_train
   character(len=4), intent(in)                        :: int_resp
   real *8, intent(in)                                 :: tol_fobj
    
   !outputs:
   real *8, dimension(rank,d,basissize), intent(inout) :: coeffs
   real *8, dimension(rank,d,N_train), intent(inout)   :: comp1d_train
   real *8, intent(out)                                :: bias_cvg
   integer, intent(out)                                :: count_iter
   logical, intent(out)                                :: flag_cvg

   !local variables:
   real *8                                             :: scaling, Ytmp1, Ytmp2, prod_tmp, betahat, &
                                                          coeffnew, err, err_old, diff_err, tol_diff, &
                                                          tol_err, normD_u, fobj, fobj_old, & 
                                                          diff_fobj, tmp, normD_ur, denom, Ydata_mean, &
                                                          Xinputs_mean, Yresp_mean, bias, normGradObj
   real *8, dimension(d)                               :: normGradObj_block
   real *8, dimension(N_train)                         :: partproduct, Xinputs, Yresp, ur_train
   integer                                             :: i, j, k, l, m, jj, ll, N
            
   normD_u = sqrt((sum(Ydata_train**2))/N_train)
   Ydata_mean = sum(Ydata_train)/N_train
   N = rank*d*basissize !total number of inputs
   
   tol_diff = 1e-2 !tolerance for saturation of error
   tol_err = 1e-2  !tolerance for convergence of relative error
!    tol_fobj = 1e-9 !tolerance for saturation of objective function
   err_old = 0.0D0
   
   fobj_old = 0.0D0
   
   count_iter = 0
   
   flag_cvg = .false. !test for saturation of the error norm ||u-ur||_D
   
!    do while (count_iter < 100) !FOR DEBUG
   
   wloop: do while (flag_cvg .eqv. .false.)
      
     !------------Big loop starts here------------
     do l = 1,rank
       do j = 1,d
         do k = 1,basissize
           !
           ! Compute partial products of one-dimensional functions
           ! In SRCD paper, see coefficients a_{jl}^{(i)}
           do i = 1,N_train
	     partproduct(i) = 1.0D0
	     do jj = 1,d
	       if (jj /= j) then
	         partproduct(i) = partproduct(i) * comp1d_train(l,jj,i)
	       end if
	     end do !jj
           end do !i
           !
           ! Compute "input data points" to be fitted
           ! In SRCD paper, see coefficients z_i^{(jlk)}
           Xinputs_mean = 0.0D0
           do i = 1,N_train
             Xinputs(i) = partproduct(i) * varphi_train(i,k,j)
             Xinputs_mean = Xinputs_mean + Xinputs(i)
           end do !i
           Xinputs_mean = Xinputs_mean/N_train
           ! Update input data points with mean value:
           do i = 1,N_train
             Xinputs(i) = Xinputs(i) - Xinputs_mean
           end do
           !
           ! Compute "modified response"
           Yresp_mean = 0.0D0
           do i = 1,N_train
             ! In SRCD paper, see coefficients c_l^{(i)}
             Ytmp1 = 0.0D0
             do ll = 1,rank
               if (ll /= l) then
                 Ytmp1 = Ytmp1 + product(comp1d_train(ll,1:d,i)) 
               end if
             end do !ll
             !
             ! In SRCD paper, see coefficients b_{jlk}^{(i)}
             Ytmp2 = 0.0D0
             do m = 1,basissize
               if (m /= k) then
                 Ytmp2 = Ytmp2 + coeffs(l,j,m)*varphi_train(i,m,j)
               end if
             end do !m
             prod_tmp = 1.0D0
             do jj = 1,d
               if (jj /= j) then
                 prod_tmp = prod_tmp * comp1d_train(l,jj,i)
               end if
             end do !jj
             Ytmp2 = prod_tmp * Ytmp2 !In SRCD paper, see coefficients a_{jl}^{(i)}*b_{jlk}^{(i)}
             !
             Yresp(i) = Ydata_train(i) - Ydata_mean - Ytmp1 - Ytmp2
             Yresp_mean = Yresp_mean + Ytmp1 + Ytmp2
           end do !i
           Yresp_mean = Yresp_mean/N_train
           ! Update response with mean value:
           do i = 1,N_train
             Yresp(i) = Yresp(i) + Yresp_mean !In SRCD paper, see coefficients r_i^{(jlk)}
           end do
           !
           ! Compute first input variable of soft-thresholding operator:
           betahat = (dot_product(Xinputs,Yresp))/N_train
           !
           ! Update coefficient by applying the soft-thresholding operator:
           call soft_threshold(betahat, lbda, coeffnew)
           denom = (dot_product(Xinputs,Xinputs))/N_train
           
!            ! ==== SO FAR: STOP AS SOON AS 'denom' IS TOO SMALL ===
!            if (abs(denom) < 1e-16) then
!              !'partproduct' is too small (can happen when 'lbda' is too large)
!              if (iprint == 1) then
!                print*, "Impossible to update coefficient (denom too small)!!!"
!                print*, 'Number of iterations=', count_iter
!              end if
!              return
!            else
!              coeffs(l,j,k) = coeffnew/denom
!            endif
           ! ===
           
           ! === NOW: DON'T STOP, UPDATE COEFFICIENT WITH ZERO ===
           if (coeffnew == 0.0D0) then
             ! in this case 'denom' is likely to be extremely small
!              if (iprint == 1) then
!                print*, 'j,l,k=', j, l, k
!                print*, 'numer=', coeffnew
!                print*, 'denom=', denom
!              end if
             coeffs(l,j,k) = coeffnew
           else
             coeffs(l,j,k) = coeffnew/denom
           endif
           ! ===
           count_iter = count_iter + 1
           !
         end do !k
         !
         ! Update 1d component function with new coefficients: 
         do i = 1,N_train
           comp1d_train(l,j,i) = dot_product(coeffs(l,j,1:basissize),varphi_train(i,1:basissize,j)) 
         end do !i
         
         ! Compute error on the training dataset:
         call eval_bias(Ydata_train, N_train, rank, d, basissize, varphi_train, coeffs, bias)
!          print*, 'bias=', bias

         call err_normD_with_bias(Ydata_train, d, N_train, rank, basissize, varphi_train, coeffs, lbda, err, &
                        fobj, normD_ur, ur_train, bias, int_resp) 
!          print*, 'Training error -->', err/normD_u
!          diff_err = abs(err-err_old)
!          err_old = err

!          if (iprint == 1) then
!            print *,'fobj=', fobj
!          end if
         diff_fobj = fobj-fobj_old
         fobj_old = fobj
         !
         if (abs(diff_fobj) < tol_fobj*(normD_u**2)) then ! test on the difference between two successive objective values
!          if (abs(diff_fobj) < tol_fobj*(normD_u**2) .and. normGradObj < 0.05) then
!          if (diff_err < tol_diff*normD_u) then ! test on the difference between two successive errors
!           if (err < tol_err*normD_u) then ! WARNING: even if ||u-ur||_D always decreases it doesn't 
!                                           ! necessarily converge to 0 (especially when the rank is small)
            !
            !Compute final bias term:
            call eval_bias(Ydata_train, N_train, rank, d, basissize, varphi_train, coeffs, bias_cvg)   
	    !
	    !Compute the (normalized) L2-norm of the gradient of the objective functional:
	    call eval_NormGradObj(N, coeffs, normGradObj, normGradObj_block, N_train, rank, d, basissize, lbda, & 
             varphi_train, Ydata_train)
	    !
	    flag_cvg = .true.
	    if (iprint == 1) then
              print*, 'Relative training error =', err/normD_u
              print*, 'RMSE =', err
              print*, 'Number of updates =', count_iter
              print*, 'normalized ||nabla(obj)||_2 =', normGradObj
              if (normGradObj > 0.05) then
                print*, " ...Warning: 1st-order optimality conditions not satisfied!!!"
                print*, "             --> decrease 'tol_fobj'"
              end if
            end if
	    return
	 elseif (diff_fobj > 0.0D0 .and. count_iter > basissize) then
! 	    !Objective function is increasing (after first pass):
! 	    print*, 'Objective function is increasing'
! 	    !
! 	    !Compute final bias term:
! 	    call eval_bias(Ydata_train, N_train, rank, d, basissize, varphi_train, coeffs, bias_cvg)
! 	    test = .true.
! 	    print*, 'Relative training error=', err/normD_u
! 	    return
         end if
         if (fobj > 1e6) then
           print *,'fobj=', fobj 
           print*, 'SRCD: objective function diverges!!!' 
           exit wloop
         endif
         ! 
         ! Compute the (normalized) L2-norm of the gradient of the objective functional 
         ! (check to see if the first-order optimality conditions are satisfied):
!          call eval_NormGradObj(N, coeffs, normGradObj, normGradObj_block, N_train, rank, d, basissize, lbda, & 
!           varphi_train, Ydata_train)
!          print*, 'normalized ||nabla(obj)||_2=', normGradObj
         !
       end do !j  
     end do !l   
     !-----------------End of Big loop-----------------
   
   end do wloop !while
      
   end subroutine coordinate_descent_with_bias
   !==========================================================================================!
   
   !==========================================================================================!
   !
   ! block version of SRCD algorithm (d blocks in total). Each block of rank*basissize variables 
   ! is updated singly (i.e., one of its variable after the other) using SRCD updates.
   !
   subroutine block_SRCD(coeffs, Ydata_train, N_train, rank, d, basissize, & 
    varphi_train, lbda, comp1d_train, bias_cvg, iprint, count_iter, int_resp, flag_cvg)
   
   !inputs:
   integer, intent(in)                                 :: N_train, rank, d, basissize, iprint
   real *8, intent(in)                                 :: lbda
   real *8, dimension(N_train), intent(in)             :: Ydata_train
   real *8, dimension(N_train,basissize,d), intent(in) :: varphi_train
   character(len=4), intent(in)                        :: int_resp
   
   !outputs:
   real *8, dimension(rank,d,basissize), intent(inout) :: coeffs
   real *8, dimension(rank,d,N_train), intent(inout)   :: comp1d_train
   real *8, intent(out)                                :: bias_cvg
   integer, intent(out)                                :: count_iter
   logical, intent(out)                                :: flag_cvg

   !local variables:
   real *8                                             :: scaling, Ytmp1, Ytmp2, prod_tmp, betahat, &
                                                          coeffnew, err, err_old, diff_err, tol_diff, &
                                                          tol_err, tol_fobj, normD_u, fobj, fobj_old, & 
                                                          diff_fobj, tmp, normD_ur, denom, Ydata_mean, &
                                                          Xinputs_mean, Yresp_mean, bias, normGradObj
   real *8, dimension(d)                               :: normGradObj_block                                                     
   real *8, dimension(N_train)                         :: partproduct, Xinputs, Yresp, ur_train
   integer                                             :: i, j, k, l, m, jj, ll, N, count_iter_max
      
   normD_u = sqrt((sum(Ydata_train**2))/N_train)
   Ydata_mean = sum(Ydata_train)/N_train
   N = rank*d*basissize !total number of inputs
!    count_iter_max = 5000*N !total number of iterations
   
   tol_diff = 1e-2 !tolerance for saturation of error
   tol_err = 1e-2  !tolerance for convergence of relative error
   tol_fobj = 1e-9 !tolerance for saturation of objective function
   err_old = 0.0D0
   fobj_old = 0.0D0
   
   count_iter = 0
   
   flag_cvg = .false. !test for saturation of the error norm ||u-ur||_D
      
   wloop: do while (flag_cvg .eqv. .false.)
     !
     !------------Big loop starts here------------
     !
     do j = 1,d !LOOP OVER THE BLOCKS
       !
!        print*,'dim=', j
       ! Update each variable inside a block one after the other:
       do l = 1,rank
         do k = 1,basissize
           !
           ! Compute partial products of one-dimensional functions
           ! In SRCD paper, see coefficients a_{jl}^{(i)}
           do i = 1,N_train
	     partproduct(i) = 1.0D0
	     do jj = 1,d
	       if (jj /= j) then
	         partproduct(i) = partproduct(i) * comp1d_train(l,jj,i)
	       end if
	     end do !jj
           end do !i
           !
           ! Compute "input data points" to be fitted
           ! In SRCD paper, see coefficients z_i^{(jlk)}
           Xinputs_mean = 0.0D0
           do i = 1,N_train
             Xinputs(i) = partproduct(i) * varphi_train(i,k,j)
             Xinputs_mean = Xinputs_mean + Xinputs(i)
           end do !i
           Xinputs_mean = Xinputs_mean/N_train
           ! Update input data points with mean value:
           do i = 1,N_train
             Xinputs(i) = Xinputs(i) - Xinputs_mean
           end do
           !
           ! Compute "modified response"
           Yresp_mean = 0.0D0
           do i = 1,N_train
             ! In SRCD paper, see coefficients c_l^{(i)}
             Ytmp1 = 0.0D0
             do ll = 1,rank
               if (ll /= l) then
                 Ytmp1 = Ytmp1 + product(comp1d_train(ll,1:d,i)) 
               end if
             end do !ll
             !
             ! In SRCD paper, see coefficients b_{jlk}^{(i)}
             Ytmp2 = 0.0D0
             do m = 1,basissize
               if (m /= k) then
                 Ytmp2 = Ytmp2 + coeffs(l,j,m)*varphi_train(i,m,j)
               end if
             end do !m
             prod_tmp = 1.0D0
             do jj = 1,d
               if (jj /= j) then
                 prod_tmp = prod_tmp * comp1d_train(l,jj,i)
               end if
             end do !jj
             Ytmp2 = prod_tmp * Ytmp2 !In SRCD paper, see coefficients a_{jl}^{(i)}*b_{jlk}^{(i)}
             !
             Yresp(i) = Ydata_train(i) - Ydata_mean - Ytmp1 - Ytmp2
             Yresp_mean = Yresp_mean + Ytmp1 + Ytmp2
           end do !i
           Yresp_mean = Yresp_mean/N_train
           ! Update response with mean value:
           do i = 1,N_train
             Yresp(i) = Yresp(i) + Yresp_mean !In SRCD paper, see coefficients r_i^{(jlk)}
           end do
           !
           ! Compute first input variable of soft-thresholding operator:
           betahat = (dot_product(Xinputs,Yresp))/N_train
           !
           ! Update coefficient by applying the soft-thresholding operator:
           call soft_threshold(betahat, lbda, coeffnew)
           denom = (dot_product(Xinputs,Xinputs))/N_train
           
!            ! ==== SO FAR: STOP AS SOON AS 'denom' IS TOO SMALL ===
!            if (abs(denom) < 1e-16) then
!              !'partproduct' is too small (can happen when 'lbda' is too large)
!              if (iprint == 1) then
!                print*, "Impossible to update coefficient (denom too small)!!!"
!                print*, 'Number of iterations=', count_iter
! !                print *,'fobj=', fobj
!              end if
!              return
!            else
!              coeffs(l,j,k) = coeffnew/denom
!            endif
!            ! ===
           
           ! === NOW: DON'T STOP, UPDATE COEFFICIENT WITH ZERO ===
           if (coeffnew == 0.0D0) then
             ! in this case 'denom' can be extremely small
!              if (iprint == 1) then
!                print*, 'j,l,k=', j, l, k
!                print*, 'numer=', coeffnew
!                print*, 'denom=', denom
!              end if
             coeffs(l,j,k) = coeffnew
           else
             coeffs(l,j,k) = coeffnew/denom
           endif
           ! ===
           !
           count_iter = count_iter + 1
           !      
         end do !k
         !
       end do !l
       !
       ! Update 1d component functions in the j-th dimension (to mimic block update):
       do l = 1,rank
         do i = 1,N_train
           comp1d_train(l,j,i) = dot_product(coeffs(l,j,1:basissize),varphi_train(i,1:basissize,j)) 
         end do
       end do
       !
       ! Compute error on the training dataset:
       call eval_bias(Ydata_train, N_train, rank, d, basissize, varphi_train, coeffs, bias)
!          print*, 'bias=', bias

       call err_normD_with_bias(Ydata_train, d, N_train, rank, basissize, varphi_train, coeffs, lbda, err, &
        fobj, normD_ur, ur_train, bias, int_resp) 
!        print*, 'Training error -->', err/normD_u
!      diff_err = abs(err-err_old)
!      err_old = err        
       if (iprint == 1) then
         print *,'fobj=', fobj
       end if
       
       diff_fobj = fobj-fobj_old
       fobj_old = fobj
       !
       if (abs(diff_fobj) < tol_fobj*(normD_u**2)) then ! test on the difference between two successive objective values
!       if (diff_err < tol_diff*normD_u) then ! test on the difference between two successive errors
!       if (err < tol_err*normD_u) then ! WARNING: even if ||u-ur||_D always decreases it doesn't 
!                                       ! necessarily converge to 0 (especially when the rank is small)
        !
        !Compute final bias term:
        call eval_bias(Ydata_train, N_train, rank, d, basissize, varphi_train, coeffs, bias_cvg)
        
        !Compute the (normalized) L2-norm of the gradient of objective w.r.t each block of variables:
! 	call eval_NormGradObj(N, coeffs, normGradObj, normGradObj_block, N_train, rank, d, basissize, lbda, & 
!               varphi_train, Ydata_train)
        
        flag_cvg = .true.
        if (iprint == 1) then
          print*, 'Relative training error=', err/normD_u
          print*, 'Number of iterations=', count_iter
!           print "(a,10f6.3)", ' normalized ||nabla_j(obj)||_2 =', normGradObj_block
!           if (maxval(normGradObj_block) > 0.05) then
!             print*, " ...Warning: 1st-order optimality conditions not satisfied!!!"
!             print*, "             --> decrease 'tol_fobj'"
!           end if
        end if
        return
       elseif (diff_fobj > 0.0D0 .and. count_iter > basissize) then
! 	    !Objective function is increasing (after first pass):
! 	    print*, 'Objective function is increasing'
! 	    !
! 	    !Compute final bias term:
! 	    call eval_bias(Ydata_train, N_train, rank, d, basissize, varphi_train, coeffs, bias_cvg)
! 	    test = .true.
! 	    print*, 'Relative training error=', err/normD_u
! 	    return
       end if
       if (fobj > 1e6) then
         print *,'fobj=', fobj 
         print*, 'block_SRCD: objective function diverges!!!' 
         exit wloop
       endif
!        if (count_iter .eq. count_iter_max) then
!          print *,'Nb. of iterations=', count_iter
!          print *,'Max iters reached.'
!          print*, 'Relative training error=', err/normD_u
!          flag_cvg = .true. !to impose the computation of testing error in cross-validation
!                            !(even if the convergence is not necessarily reached...)
!          exit wloop
!        end if
       !
       end do !j
       !
       !-----------------End of Big loop-----------------
       !
   end do wloop !while
   
   end subroutine block_SRCD
   !==========================================================================================!
   
   !==========================================================================================!
   !
   ! BSR algorithm (d blocks in total). Each block of rank*basissize variables 
   ! is updated using a BFGS algorithm. To be tested to see if it performs better than CD updates
   ! and if it allows to get rid of the divergence observed while using block_SRCD().
   !
   ! Remark: There is NO bias term in this formulation.
   !
   subroutine BSR_subpbs_bfgs(coeffs, Ydata_train, N_train, rank, d, basissize, & 
    varphi_train, lbda, comp1d_train, iprint, count_iter, int_resp, flag_cvg, tol_fobj)
   
   !inputs:
   integer, intent(in)                                 :: N_train, rank, d, basissize, iprint
   real *8, intent(in)                                 :: lbda
   real *8, dimension(N_train), intent(in)             :: Ydata_train
   real *8, dimension(N_train,basissize,d), intent(in) :: varphi_train
   character(len=4), intent(in)                        :: int_resp
   real *8, intent(in)                                 :: tol_fobj
   
   !outputs:
   real *8, dimension(rank,d,basissize), intent(inout) :: coeffs
   real *8, dimension(rank,d,N_train), intent(inout)   :: comp1d_train
   integer, intent(out)                                :: count_iter
   logical, intent(out)                                :: flag_cvg
   
   !local variables:
   real *8                                             :: normD_u, Ydata_mean, tol_diff, tol_err, &
                                                          err, err_old, fobj_old, fobj, &
                                                          diff_fobj, factr, pgtol, normD_ur
   integer                                             :: i, j, k, l, N, optf
   real *8, dimension(rank*basissize)                  :: coeffs_block
   real *8, dimension(N_train)                         :: ur_train

   normD_u = sqrt((sum(Ydata_train**2))/N_train)
   Ydata_mean = sum(Ydata_train)/N_train
   N = rank*d*basissize !total number of inputs
   
   tol_diff = 1e-2 !tolerance for saturation of error
   tol_err = 1e-2  !tolerance for convergence of relative error
!    tol_fobj = 1e-9 !tolerance for saturation of objective function
   err_old = 0.0D0
   fobj_old = 0.0D0
   
   !BFGS parameters:
   optf = -1 !display values of objective function during optimization
   factr = 1.d+10
!    factr = 1.d+12 !low accuracy
!    factr = 1.d+7  !moderate accuracy
!    factr = 1.d+1  !extremely high accuracy

   pgtol = 1.0d-5
   
   count_iter = 0
   
   flag_cvg = .false. !test for saturation of the error norm ||u-ur||_D
   
   wloop: do while (flag_cvg .eqv. .false.)
     !
     !------------Big loop starts here------------
     !
     do j = 1,d !LOOP OVER THE BLOCKS (or dimension)
!        print*,'dim=', j 
       ! Extract the j-th block of variables in a 1d data structure:
       do l = 1,rank
         do k = 1,basissize
           coeffs_block((l-1)*basissize+k) = coeffs(l,j,k)
         end do
       end do  
       ! BFGS resolution to update the j-th block of variables:
       call LBFGSDriver_block(coeffs_block, rank*basissize, Ydata_train, N_train, & 
        rank, d, basissize, varphi_train, comp1d_train, j, lbda, fobj, factr, pgtol, optf)
       
       ! Update block of coefficients in 3d data structure:
       do l = 1,rank
         do k = 1,basissize
           coeffs(l,j,k) = coeffs_block((l-1)*basissize+k)
         end do
       end do 
       
       ! Update 1d component functions:
       do l = 1,rank
         do i = 1,N_train
           comp1d_train(l,j,i) = dot_product(coeffs(l,j,1:basissize),varphi_train(i,1:basissize,j)) 
         end do
       end do
       
       ! Compute error on the training dataset:
       call err_normD(Ydata_train, d, N_train, rank, basissize, varphi_train, coeffs, lbda, err, fobj, & 
             normD_ur, ur_train, int_resp)
!        print*, 'Training error -->', err/normD_u
!      diff_err = abs(err-err_old)
!      err_old = err        
       if (iprint == 1) then
!          print *,'fobj=', fobj
       end if
       
       diff_fobj = fobj-fobj_old
       fobj_old = fobj
       count_iter = count_iter + 1
       
       if (abs(diff_fobj) < tol_fobj*(normD_u**2)) then ! test on the difference between two successive objective values
!       if (diff_err < tol_diff*normD_u) then ! test on the difference between two successive errors
!       if (err < tol_err*normD_u) then ! WARNING: even if ||u-ur||_D always decreases it doesn't 
!                                       ! necessarily converge to 0 (especially when the rank is small)  
         flag_cvg = .true.
         if (iprint == 1) then
           print*, 'Relative training error=', err/normD_u
           print*, 'Number of iterations=', count_iter
         end if
         return
       end if
       if (fobj > 1e6) then
         print *,'fobj=', fobj 
         print*, 'BSR_subpbs_bfgs: objective function diverges!!!' 
         STOP
!          exit wloop
       endif
     
     end do !j
     !
     !-----------------End of Big loop-----------------
     !
   end do wloop !while
   
   end subroutine BSR_subpbs_bfgs
   !==========================================================================================!
   
   !==========================================================================================!
   !
   ! BSR algorithm (d blocks in total). Each block of rank*basissize variables 
   ! is updated using GLMNET software library (by Friedman, Hastie, Tibshirani and Simon).
   !
   subroutine BSR_subpbs_glmnet(coeffs, Ydata_train, N_train, rank, d, basissize, & 
    varphi_train, lbda, comp1d_train, iprint, count_iter, int_resp, flag_cvg, tol_fobj)
   
   !inputs:
   integer, intent(in)                                 :: N_train, rank, d, basissize, iprint
   real *8, intent(in)                                 :: lbda
   real *8, dimension(N_train), intent(in)             :: Ydata_train
   real *8, dimension(N_train,basissize,d), intent(in) :: varphi_train
   character(len=4), intent(in)                        :: int_resp
   real *8, intent(in)                                 :: tol_fobj
   
   !outputs:
   real *8, dimension(rank,d,basissize), intent(inout) :: coeffs
   real *8, dimension(rank,d,N_train), intent(inout)   :: comp1d_train
   integer, intent(out)                                :: count_iter
   logical, intent(out)                                :: flag_cvg
   
   !local variables:
   real *8                                             :: normD_u, Ydata_mean, tol_diff, tol_err, &
                                                          err, err_old, fobj_old, fobj, diff_fobj, &
                                                          normD_ur, partproduct                                      
   integer                                             :: i, ii, j, jj, k, l, lr, lb
   real *8, dimension(N_train)                         :: ur_train
   integer, dimension(:,:), allocatable                :: inds2d
   real *8, dimension(:,:), allocatable                :: Aj
   
   !local variables needed for elnet():
   real *4                                             :: parm, flmin, thr 
   integer                                             :: ka, no, ni, jd1, ne, nx, nlam, isd, intr, & 
                                                          maxit, lmu, nlp, jerr
   integer, dimension(:), allocatable                  :: jd, ia, nin
   real *4, dimension(:), allocatable                  :: y, w, vp, ulam, a0, rsq, alm
   real *4, dimension(:,:), allocatable                :: cl, ca, x
   real *8, dimension(:,:), allocatable                :: ca_dble
   
   normD_u = sqrt((sum(Ydata_train**2))/N_train)
   Ydata_mean = sum(Ydata_train)/N_train
   
   tol_diff = 1e-2 !tolerance for saturation of error
   tol_err = 1e-2  !tolerance for convergence of relative error
   err_old = 0.0D0
   fobj_old = 0.0D0
   
   !======================= GLMNET parameters =======================!
   ! --- input params:
   ! Remark: real parameters must be defined as real*4 when calling elnet()
   ka = 2 !*** NOT SURE (algo flag) ***
   parm = 1.0D0 ! parm=1 => lasso minimiz / parm=0 => ridge minimiz
   no = N_train !nb of observations
   ni = rank*basissize !nb of predictor variables
   allocate(w(no), x(no,ni), y(no), Aj(no,ni))
   w(:) = 1.0D0 !observation weights
   jd1 = 0
   allocate(jd(jd1+1))
   jd(1) = jd1
   allocate(vp(ni))
   vp(:) = 1.0D0 !*** NOT SURE (relative penalties) ***
   allocate(cl(2,ni))
   cl(1,:) = -100.0D0 !*** NOT SURE (lower bounds for input vars) ***
   cl(2,:) = 100.0D0 !*** NOT SURE (upper bounds for input vars) ***
   nx = ni !*** NOT SURE (max nb of vars allowed to enter all models) ***
   ne = nx-1 !*** NOT SURE (max nb of vars allowed to enter largest model) ***
   nlam = 1
   flmin = 1.0D0 !*** NOT SURE (user control of lambda values) ***
   allocate(ulam(nlam))
   ulam(1) = lbda
   thr = 1e-5
   isd = 1 !*** NOT SURE (standarization flag) ***
   intr = 0 !no bias (intercept) term
   maxit = 100000
   y(:) = Ydata_train(:)
   ! --- output params:
   lmu = 1
   allocate(a0(lmu), ca(nx,lmu), ca_dble(nx,lmu), ia(nx), nin(lmu), rsq(lmu), alm(lmu))
   !=================================================================!
   
   count_iter = 0
   
   flag_cvg = .false. !test for saturation of the error norm ||u-ur||_D
      
   wloop: do while (flag_cvg .eqv. .false.)
     !
     !------------Big loop starts here------------
     !
     do j = 1,d !LOOP OVER THE BLOCKS
!        print*,'dim=', j 
       
        ! assemble the data matrix Aj (real*8):
        do l = 1,rank !column-block index
          ! Compute partial products of one-dimensional functions (Cf. coefficients a_{jl}^{(i)})
          do i = 1,no !row index
            partproduct = 1.0D0
            do jj = 1,d
              if (jj /= j) then
                partproduct = partproduct * comp1d_train(l,jj,i)
              end if
            end do !jj
            Aj(i,(l-1)*basissize+1:l*basissize) = partproduct*varphi_train(i,:,j)
          end do !i
        end do !l
        x = real(Aj,4) !real*4 data matrix to be used in elnet()
        
!         do jj=1,no
!           print*,'row=', jj, 'x=', x(jj,:)
!         enddo
        
        ! call elnet solver (dense predictor matrix):
        ca(:,1) = 12345
        ia(:) = -666
        print*,'before call to elnet(): ulam=', ulam(1)
        call elnet(ka, parm, no, ni, x, y, w, jd, vp, cl, ne, nx, nlam, flmin, & 
              ulam, thr, isd, intr, maxit, lmu, a0, ca, ia, nin, rsq, alm, nlp, jerr)
        
        print*, ''
        print*, 'block nb=', j
        print*,'after call to elnet(): ca=', ca(:,1)
        print*,'after call to elnet(): ia=', ia(:)
        print*,'after call to elnet(): nin=', nin(1)
        print*,'after call to elnet(): rsq=', rsq
        print*,'after call to elnet(): jerr=', jerr
        
        if (jerr .ne. 0) then
          print*, ''
          print*,"Something's wrong in elnet()!!!"
          print*,'   jerr=', jerr
          STOP
        end if
        
        ! Make a copy in double precision for f90 consistency:
        ca_dble = real(ca,8) 

        ! Remap each index of compressed coefficients into a pair of 
        ! indices which can be used in our 3d data structure:
        allocate(inds2d(2,nin(1)))
        call remap_indices(nin(1), ia, basissize, inds2d)
                
        ! Update block of coefficients in 3d data structure:
        coeffs(:,j,:) = 0.0D0
        ! Insert the (nonzero) compressed coefficients:
        do ii = 1,nin(1)
          lr = inds2d(1,ii) !level rank index
          lb = inds2d(2,ii) !level basis index
          coeffs(lr,j,lb) = ca_dble(ii,1)
!           print*, 'ii=', ii, 'lr=', lr, 'lb=', lb, 'ca=', ca_dble(ii,1)
        end do !ii
        deallocate(inds2d)
        
        ! If no compressed coefficients in elnet() we'd have updates like:
!         do l = 1,rank
!           do k = 1,basissize
!             coeffs(l,j,k) = ca_dble((l-1)*basissize+k,1)
!           end do
!         end do 
       
        ! Update 1d component functions:
        do l = 1,rank
          do i = 1,N_train
            comp1d_train(l,j,i) = dot_product(coeffs(l,j,1:basissize),varphi_train(i,1:basissize,j)) 
          end do
        end do
   
        ! Compute error on the training dataset:
        call err_normD(Ydata_train, d, N_train, rank, basissize, varphi_train, coeffs, lbda, err, fobj, & 
              normD_ur, ur_train, int_resp)
        
        if (iprint == 1) then
          print *,'fobj=', fobj
        end if
        
        diff_fobj = fobj-fobj_old
        fobj_old = fobj
        count_iter = count_iter + 1
       
        ! test on the difference between 2 successive objective values
        if (abs(diff_fobj) < tol_fobj*(normD_u**2)) then
          flag_cvg = .true.
          if (iprint == 1) then
            print*, 'max|coeffs|=', maxval(abs(coeffs))
            print*, 'Relative training error=', err/normD_u
            print*, 'Number of iterations=', count_iter
          end if
          deallocate(x, y, Aj, w, jd, vp, cl, ulam, a0, ca, ca_dble, ia, nin, rsq, alm)
          return
        end if
        if (fobj > 1e6) then
          print *,'fobj=', fobj 
          print*, 'BSR_subpbs_glmnet: objective function diverges!!!' 
          STOP
!           exit wloop
        endif
   
     end do !j
     !
     !-----------------End of Big loop-----------------
     !
   end do wloop !while  
      
   deallocate(x, y, Aj, w, jd, vp, cl, ulam, a0, ca, ca_dble, ia, nin, rsq, alm)
   
   end subroutine BSR_subpbs_glmnet
   !==========================================================================================!
   
   !==========================================================================================!
   !
   ! Soft-thresholding operator
   !
   subroutine soft_threshold(x, lbda, y)
   
   !inputs:
   real *8, intent(in)    :: x, lbda
   !outputs:
   real *8, intent(out)   :: y
   
!    y = sign(1.0D0,x) * max(abs(x)-lbda,0.0D0)

   if (abs(x) .le. lbda) then
     y = 0.0D0
   else
     if (x > 0.0D0) then
       y = x - lbda
     elseif (x < 0.0D0) then
       y = x + lbda
     end if
   end if
   
   end subroutine soft_threshold
   !==========================================================================================!

   !==========================================================================================!
   ! 
   ! Evaluate the bias term
   !
   subroutine eval_bias(Ydata_train, N_train, rank, d, basissize, varphi_train, coeffs, bias)
   
   !inputs:
   integer, intent(in)                                 :: N_train, rank, d, basissize
   real *8, dimension(N_train), intent(in)             :: Ydata_train
   real *8, dimension(N_train,basissize,d), intent(in) :: varphi_train
   real *8, dimension(rank,d,basissize), intent(in)    :: coeffs
   
   !outputs:
   real *8, intent(out)                                :: bias
   
   !local variables:
   integer                                             :: i, j, l
   real *8                                             :: Ydata_mean, ur, ur_l
   real *8, dimension(N_train)                         :: urdata
   
   Ydata_mean = sum(Ydata_train)/N_train
   
   do i = 1,N_train
      ur = 0.0d0
      do l = 1,rank !rank loop
         ur_l = 1.0d0
         do j = 1,d ! loop on dimension
            ur_l = ur_l * dot_product(coeffs(l,j,:),varphi_train(i,:,j))
         end do !j
         ur = ur + ur_l !add contribution of rank level l
      end do !l
      urdata(i) = ur !store approximate value
   end do !i
   bias = Ydata_mean - sum(urdata)/N_train 
   
   end subroutine eval_bias
   !==========================================================================================!
   
   !==========================================================================================!
   !
   ! block version of SRCD algorithm (d blocks in total) using LASSO updates for each block of
   ! variables. Warning: THIS VERSION NEEDS TO BE DEBUGGED!!!
   !
   subroutine block_SRCD_lasso(coeffs, Ydata_train, N_train, rank, d, basissize, & 
    varphi_train, lbda, comp1d_train, bias_cvg, iprint, count_iter, int_resp)
   
   !inputs:
   integer, intent(in)                                 :: N_train, rank, d, basissize, iprint
   real *8, intent(in)                                 :: lbda
   real *8, dimension(N_train), intent(in)             :: Ydata_train
   real *8, dimension(N_train,basissize,d), intent(in) :: varphi_train
   character(len=10), intent(in)                       :: int_resp
   
   !outputs:
   real *8, dimension(rank,d,basissize), intent(inout) :: coeffs
   real *8, dimension(rank,d,N_train), intent(inout)   :: comp1d_train
   real *8, intent(out)                                :: bias_cvg
   integer, intent(out)                                :: count_iter

   !local variables:
   logical                                             :: test
   real *8                                             :: scaling, Ytmp1, Ytmp2, prod_tmp, betahat, &
                                                          coeffnew, err, err_old, diff_err, tol_diff, &
                                                          tol_err, tol_fobj, normD_u, fobj, fobj_old, & 
                                                          diff_fobj, tmp, normD_ur, denom, Ydata_mean, &
                                                          Xinputs_mean, Yresp_mean, bias, normGradObj
   real *8, dimension(d)                               :: normGradObj_block                                                     
   real *8, dimension(N_train)                         :: partproduct, Xinputs, Yresp, ur_train
   integer                                             :: i, j, k, l, m, jj, jjj, ll, N
   real *8, dimension(:,:,:), allocatable              :: matA, matA_shift
   real *8, dimension(:,:), allocatable                :: matA_mean, coeffs_block
   
   normD_u = sqrt((sum(Ydata_train**2))/N_train)
   Ydata_mean = sum(Ydata_train)/N_train
   N = rank*d*basissize !total number of inputs
   
   ! Recast variables by blocks:
   allocate(coeffs_block(d,rank*basissize))
   do j = 1,d
     do l = 1,rank
       coeffs_block(j,(l-1)*basissize+1:l*basissize) = coeffs(l,j,:)
     end do !l   
   end do !j
   
   ! Assemble matA:
   allocate(matA(d,N_train,rank*basissize), matA_shift(d,N_train,rank*basissize))
   matA(:,:,:) = 0.0D0
   matA_shift(:,:,:) = 0.0D0
   do j = 1,d
     ! assemble A^j (matrix of size N_train x rank*basissize)
     do l = 1,rank
       ! Compute partial products of one-dimensional functions (Cf. coefficients a_{jl}^{(i)})
       do i = 1,N_train
         partproduct(i) = 1.0D0
         do jj = 1,d
           if (jj /= j) then
             partproduct(i) = partproduct(i) * comp1d_train(l,jj,i)
           end if
         end do !jj
         ! assemble submatrices of size N_train x basissize
         matA(j,i,(l-1)*basissize+1:l*basissize) = partproduct(i)*varphi_train(i,:,j)
       end do !i
     end do !l
   end do !j
   ! compute the mean of the columns of each submatrix:
   allocate(matA_mean(d,rank*basissize))
   matA_mean(:,:) = 0.0D0
   do j = 1,d
     do m = 1,rank*basissize
       matA_mean(j,m) = sum(matA(j,:,m))/N_train
     end do !m
   end do !j
   
   tol_diff = 1e-2 !tolerance for saturation of error
   tol_err = 1e-2  !tolerance for convergence of relative error
   tol_fobj = 1e-9 !tolerance for saturation of objective function
   err_old = 0.0D0
   
   fobj_old = 0.0D0
   
   count_iter = 0
   
   test = .false. !test for saturation of the error norm ||u-ur||_D
      
   do while (test .eqv. .false.)
     !
     !------------Big loop starts here------------
     !
     do j = 1,d !LOOP OVER THE BLOCKS
       !
       ! Remove mean of column to submatrix column:
       do k = 1,rank*basissize
         matA_shift(j,:,k) = matA(j,:,k) - matA_mean(j,k)
       end do !k
       
       ! Update each variable of the j-th block:
       do m = 1,rank*basissize
         ! Compute the denominator:
         denom = sum((matA_shift(j,:,m))**2)/N_train
         ! Compute the partial residual:
         partproduct(:) = Ydata_train(:) - Ydata_mean - matmul(matA_shift(j,:,:),coeffs_block(j,:)) + &
          coeffs_block(j,m)*matA_shift(j,:,m)
         ! Compute the numerator:
         betahat = (dot_product(partproduct,matA_shift(j,:,m)))/N_train
         ! Apply the soft-thresholding operator:
         call soft_threshold(betahat, lbda, coeffnew)
         ! Update coefficient:
         if (coeffnew == 0.0D0) then
           coeffs_block(j,m) = coeffnew
         else
           coeffs_block(j,m) = coeffnew/denom
         endif
       end do !m
       
       ! Update 3d data structure coefficients:
       do l = 1,rank
         do k = 1,basissize
           coeffs(l,j,k) = coeffs_block(j,(l-1)*basissize+k) 
         end do
       end do
       
       ! Update 1d component functions in the j-th dimension:
       do l = 1,rank
         do i = 1,N_train
           comp1d_train(l,j,i) = dot_product(coeffs(l,j,1:basissize),varphi_train(i,1:basissize,j)) 
         end do
       end do
       
       ! Update submatrices A^jjj for jjj different from j:
       do jjj = 1,d
         if (jjj /= j) then
           do l = 1,rank
             ! Compute partial products of one-dimensional functions (Cf. coefficients a_{jl}^{(i)})
             do i = 1,N_train
               partproduct(i) = 1.0D0
               do jj = 1,d
                 if (jj /= jjj) then
                   partproduct(i) = partproduct(i) * comp1d_train(l,jj,i)
                 end if
               end do !jj
               ! assemble submatrices of size N_train x basissize
               matA(jjj,i,(l-1)*basissize+1:l*basissize) = partproduct(i)*varphi_train(i,:,jjj)
             end do !i
           end do !l
         end if   
       end do !jjj
       
       ! Compute error on the training dataset:
       call eval_bias(Ydata_train, N_train, rank, d, basissize, varphi_train, coeffs, bias)
       
       call err_normD_with_bias(Ydata_train, d, N_train, rank, basissize, varphi_train, coeffs, lbda, err, &
        fobj, normD_ur, ur_train, bias, int_resp)

       if (fobj > 1e9) then
          print *,'fobj=', fobj 
          print*, 'block_SRCD: objective function diverges!!!' 
          print*, 'STOP.'
          STOP
       endif
        
       if (iprint == 1) then
!         print *,'fobj=', fobj
       end if
       
       diff_fobj = fobj-fobj_old
       fobj_old = fobj
     
       if (abs(diff_fobj) < tol_fobj*(normD_u**2)) then ! test on the difference between two successive objective values
!       if (diff_err < tol_diff*normD_u) then ! test on the difference between two successive errors
!       if (err < tol_err*normD_u) then ! WARNING: even if ||u-ur||_D always decreases it doesn't 
!                                       ! necessarily converge to 0 (especially when the rank is small)
        !
        !Compute final bias term:
        call eval_bias(Ydata_train, N_train, rank, d, basissize, varphi_train, coeffs, bias_cvg)
        
        !Compute the (normalized) L2-norm of the gradient of objective w.r.t each block of variables:
! 	call eval_NormGradObj(N, coeffs, normGradObj, normGradObj_block, N_train, rank, d, basissize, lbda, & 
!               varphi_train, Ydata_train)
        
        test = .true.
        if (iprint == 1) then
          print*, 'Relative training error=', err/normD_u
          print*, 'Number of iterations=', count_iter
!           print "(a,10f6.3)", ' normalized ||nabla_j(obj)||_2 =', normGradObj_block
!           if (maxval(normGradObj_block) > 0.05) then
!             print*, " ...Warning: 1st-order optimality conditions not satisfied!!!"
!             print*, "             --> decrease 'tol_fobj'"
!           end if
        end if
        return
       elseif (diff_fobj > 0.0D0 .and. count_iter > basissize) then
! 	    !Objective function is increasing (after first pass):
! 	    print*, 'Objective function is increasing'
! 	    !
! 	    !Compute final bias term:
! 	    call eval_bias(Ydata_train, N_train, rank, d, basissize, varphi_train, coeffs, bias_cvg)
! 	    test = .true.
! 	    print*, 'Relative training error=', err/normD_u
! 	    return
       end if
     
     end do !j
     !
     !-----------------End of Big loop-----------------
     !
   end do !while
   
   deallocate(matA, matA_shift, matA_mean, coeffs_block)
   
   end subroutine block_SRCD_lasso
   !==========================================================================================!
   
   
end module coordDescent