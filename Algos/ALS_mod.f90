module ALS

   use utilities
   implicit none

   contains
   
   !==========================================================================================!
   !
   ! ALS scheme in each dimension along with Tykhonov regularization (the rank is fixed here).
   ! This version follows 2014 Validi's paper, where a criterion based on SVD 
   ! is used to decide whether or not a 1d least-squares problem needs to be regularized.
   ! Note: for a given dimension, search for related unknown coefficients of all rank levels. 
   !
   subroutine classicalALS(d, mpts, rank, Mbasis, sconsts, coeffs, varphi, Ydata, &
    normD_u, iterALS, int_resp)
   
   !inputs:
   integer, intent(in)                              :: d, mpts, rank, Mbasis
   real *8, dimension(mpts,Mbasis,d), intent(in)    :: varphi
   real *8, dimension(mpts), intent(in)             :: Ydata
   real *8, intent(in)                              :: normD_u
   character(len=4), intent(in)                     :: int_resp
   
   !Outputs:
   integer, intent(out)                             :: iterALS
   real *8, dimension(rank), intent(inout)          :: sconsts
   real *8, dimension(rank,d,Mbasis), intent(inout) :: coeffs
   
   !Local:
   integer                                          :: iter, i, j, k, l, m, info, l1, l2, & 
                                                       lwork_dgesvd, lbd_pts, ind_opt
   logical                                          :: test
   real *8                                          :: err_old, err, tmp, lbd_opt, lbd_min, &
                                                       lbd_max, lbd_step, lbd
   real *8                                          :: normD_uil
   real *8                                          :: tol, normD_ur, diff_err, tracehatmat
   real *8, dimension(mpts,rank*Mbasis)             :: Ai
   real *8, dimension(mpts,Mbasis)                  :: Ai_l
   real *8, dimension(rank*Mbasis,rank*Mbasis)      :: Bi, mat, mat_backup, invmat
   integer, dimension(rank*Mbasis)                  :: ipiv
   real *8, dimension(mpts,mpts)                    :: hatmat
   real *8, dimension(rank*Mbasis)                  :: rhs, rhs_work, ci, ci_lbd, rhs_ek
   real *8, dimension(Mbasis,Mbasis)                :: IdM
   real *8, dimension(:), allocatable               :: work_dgesvd
   real *8, dimension(mpts,mpts)                    :: U
   real *8, dimension(rank*Mbasis,rank*Mbasis)      :: VT
   real *8, dimension(:), allocatable               :: S
   real *8, dimension(:), allocatable               :: numer, denom, GCV
   real *8, dimension(mpts)                         :: ur_train
   
!    print *,' ...classicalALS: Be careful with the definition of search range in GCV!!!'
   
   test = .false. !test for saturation of the error norm ||u-ur||_D
   err_old = 0.0D0
   iter = 0
   tol = 1e-6
   
   IdM(:,:) = 0.0D0
   do i = 1,Mbasis
    IdM(i,i) = 1.0D0
   end do
   
   !Allocations for SVD:
   lwork_dgesvd = max(3*min(mpts, rank*Mbasis) + max(mpts,rank*Mbasis), 5*min(mpts,rank*Mbasis)-4)
   allocate(work_dgesvd(lwork_dgesvd))
   allocate(S(min(mpts,rank*Mbasis)))
   !Allocations for cross-validation (optimal regularization parameters):
   lbd_pts = 100
   allocate(numer(lbd_pts+1))
   allocate(denom(lbd_pts+1))
   allocate(GCV(lbd_pts+1))
      
   do while (test .eqv. .false.)
      
      do i = 1,d !loop over dimension
            
	  !Assemble matrix for least square minimization in ith-dimension:
          Ai(:,:) = 0.0D0
          do l = 1,rank
	      !Assemble each column block:
              Ai_l(:,:) = 0.0D0
	      do j = 1,mpts
		do k = 1,Mbasis
		  !compute \prod_m [u_m^l(x^{(j)}_m)] for m \neq i
		  tmp = 1.0D0
		  do m = 1,d
		    if (m /= i) then
		      tmp = tmp * dot_product(varphi(j,:,m),coeffs(l,m,:))
		    end if
		  end do !m
		  Ai_l(j,k) = sconsts(l)*varphi(j,k,i)*tmp
		end do !k
	      end do !j
	      Ai(:,(l-1)*Mbasis+1:l*Mbasis) = Ai_l
          end do !l
          
          !********************* TYKHONOV REGULARIZATION ***************************
          !---Assemble regularization matrix:
          Bi(:,:) = 0.0D0
          do l1 = 1,rank
	    do l2 = 1,rank
	      tmp = 1.0D0
              do k = 1,d
		if (k /= i) then
		  tmp = tmp * dot_product(coeffs(l1,k,:),coeffs(l2,k,:))
		end if
              end do !k
              Bi((l1-1)*Mbasis+1:l1*Mbasis,(l2-1)*Mbasis+1:l2*Mbasis) = &
               sconsts(l1) * sconsts(l2) * tmp * IdM(:,:)
            end do !l2
          end do !l1
          
          !---Selection of Tykhonov regularization parameter (using Generalized Cross-Validation):
	  ! Singular Value Decomposition of Ai:
	  call DGESVD('A', 'A', mpts, rank*Mbasis, real(Ai,8), mpts, S, U, mpts, VT, rank*Mbasis, &
	   work_dgesvd, lwork_dgesvd, info)
          lbd_min = 0.75D0*minval(S) 
          lbd_max = 5.0D0*maxval(S)
          lbd_step = (lbd_max-lbd_min)/float(lbd_pts)
          
          do j = 1,lbd_pts+1
	    lbd = lbd_min + (j-1)*lbd_step
	    !--compute numerator:
	    mat = matmul(transpose(Ai),Ai) + (lbd**2)*Bi
	    mat_backup = mat
	    rhs = matmul(transpose(Ai),Ydata)
	    
! 	    call DGETRF(rank*Mbasis, rank*Mbasis, mat, rank*Mbasis, ipiv, info) !LU factorization of mat 
! 	    if (info < 0) then
! 	      print*, 'algos_mod (DGETRF): Argument number', info, 'has an illegal value!!!'
!               stop
!             end if
!             rhs_work = rhs
! 	    call DGETRS('N', rank*Mbasis, 1, mat, rank*Mbasis, ipiv, rhs_work, rank*Mbasis, info) !linear system resolution
!             if (info < 0) then
! 	      print*, 'algos_mod (DGETRS): Argument number', info, 'has an illegal value!!!'
!               stop
!             end if
            
            rhs_work = rhs
            call DPOSV('U', rank*Mbasis, 1, mat, rank*Mbasis, rhs_work, rank*Mbasis, info)
	    if (info < 0) then
	      print*, 'algos_mod (DPOSV): Argument number', info, 'has an illegal value!!!'
	      stop
	    elseif (info > 0) then  
	      print*, 'algos_mod (DPOSV): Factorization could not be completed!!!'
              stop
            end if
	    ci_lbd = rhs_work
	    tmp = mpts*sum((matmul(Ai,ci_lbd)-Ydata)**2)
	    numer(j) = tmp
	    !-- compute inverse of mat:
	    invmat(:,:) = 0.0D0
	    do k = 1,rank*Mbasis
	      rhs_ek(:) = 0.0D0
	      rhs_ek(k) = 1.0D0
	      mat = mat_backup
	      call DPOSV('U', rank*Mbasis, 1, mat, rank*Mbasis, rhs_ek, rank*Mbasis, info)
	      if (info < 0) then
		print*, 'algos_mod (DPOSV): Argument number', info, 'has an illegal value!!!'
		stop
	      elseif (info > 0) then  
		print*, 'algos_mod (DPOSV): Factorization could not be completed!!!'
		stop
	      end if      
	      invmat(:,k) = rhs_ek
	    end do
	    !--compute denominator:
	    hatmat = matmul(Ai,matmul(invmat,transpose(Ai)))
	    call trace(hatmat, mpts, tracehatmat)
	    denom(j) = (mpts-tracehatmat)**2
	    GCV(j) = numer(j)/denom(j)
          end do !j
          
!           ! determine best value for lambda (see eq. (28) in Doostan & Iaccarino 2013):
          ind_opt = minloc(GCV, dim=1)
          lbd_opt = lbd_min + (ind_opt-1)*lbd_step
          
!           print*,'lbd_min=', lbd_min
!           print*,'lbd_opt=', lbd_opt
          
          if ( (lbd_opt > 0.99D0*lbd_min) .and. (lbd_opt < 1.01D0*lbd_min) ) then
	    !regularization is not needed:
	    lbd_opt = 0.0D0
	    print*, 'no regularization is needed'
          else
	    print*, ' ...Regularization: lbd_opt=', lbd_opt
          end if
          !*************************************************************************  
          
          !Skip regularization:
!           lbd_opt = 0.0D0
!           Bi(:,:) = 0.0D0
          
          !---Solve for PC coefficients in i-th dimension using (regularized) normal equations:
	  mat = matmul(transpose(Ai),Ai) + (lbd_opt**2)*Bi
	  rhs = matmul(transpose(Ai),Ydata)
	  
! 	  do j=1,rank*Mbasis
! 	    print*,'mat=', mat(j,:)
! 	  enddo
	  
          call DPOSV('U', rank*Mbasis, 1, mat, rank*Mbasis, rhs, rank*Mbasis, info)
          ci = rhs
          
          !---Update PC coefficients:
          do l = 1,rank
	    coeffs(l,i,1:Mbasis) = ci((l-1)*Mbasis+1:l*Mbasis)
          end do
          
          !---Normalizations (Cf algorithm 1 in Doostan-Iaccarino 2013):
          do l = 1,rank
	    call normD_1dfactor(coeffs, i, l, mpts, Mbasis, rank, d, varphi, normD_uil)
	    sconsts(l) = sconsts(l)*normD_uil
	    coeffs(l,i,1:Mbasis) = coeffs(l,i,1:Mbasis)/normD_uil
          end do
          iter = iter+1
                              
          !---Compute error ||u-ur||_D and ||ur||_D:        
          call err_normD_ALS(Ydata, sconsts, coeffs, rank, Mbasis, varphi, d, mpts, & 
                             err, normD_ur, ur_train, int_resp)
   
          print*, '||u-ur||_D/||u||_D=', err/normD_u
          
          diff_err = abs(err-err_old)
          err_old = err
          
          if (diff_err < tol*normD_u) then ! test on the difference between two successive errors
!           if (err < 1e-2*normD_u) then ! don't use this stopping criterion in classical ALS; even if ||u-ur||_D 
!                                        ! always decreases it doesn't necessarily converge to 0 (for small rank).
	    test = .true.
	    print*, '||u-ur||_D/||u||_D=', err/normD_u
	    iterALS = iter
	    return
          end if
 
      end do !i 
   
   end do !while
   
   deallocate(work_dgesvd, S, numer, denom, GCV)

   end subroutine classicalALS
   !==========================================================================================!

end module ALS
