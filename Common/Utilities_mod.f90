module utilities
   
   use lpp
   implicit none

   contains
   
   !==========================================================================================!
   !
   ! Evaluate functions on dataset of points to generate a vector of response values.
   !
   subroutine evalfunc(Xin, d, mpts, opt_func, Y)
   
   !inputs:
   integer, intent(in)                       :: d, mpts
   real *8, dimension(d,mpts), intent(in)    :: Xin
   character(len=10), intent(in)             :: opt_func
   !Outputs:
   real *8, dimension(mpts), intent(out)     :: Y 
   !Local:
   integer                                   :: i, j, k
   real *8, dimension(d)                     :: tmp, ci1, ci2, ci3, xj
   real *8, parameter                        :: pi = 4.0D0*atan(1.0D0)
   real *8, dimension(5)                     :: coeffs, xi
   real *8                                   :: omega, partsum
   
   select case(opt_func)
   
      case('Rastrigin')
	! Rastrigin function: 
	! Original function: 10*d + \sum_{i=1}^d (xi^2 - 10*cos(2*pi*xi)) with xi \in [-5.12,5.12]
	do j = 1,mpts
! 	  xj = 10.24D0*Xin(:,j) - 5.12D0 !rescale input variables in [-5.12,5.12]
	  tmp(:) = 0.0D0
	  do i = 1,d
	    tmp(i) = cos(2.0D0*pi*Xin(i,j)) 
	  end do
	  Y(j) = 10.0D0*d + sum(Xin(:,j)**2) - 10.0D0*sum(tmp)
	end do 
	
      case('Friedman1')
	 ! Friedman1 function:
	 do j = 1,mpts
	    Y(j) = 10.0D0*sin(pi*Xin(1,j)*Xin(2,j)) + 20.0D0*(Xin(3,j)-0.5D0)**2 + 10.0D0*Xin(4,j) + 5.0D0*Xin(5,j)
	 end do
	 
      case('corner')
         ! Corner peak function:
!          do i = 1,d
!            ci1(i) = (float(i) - 0.5D0)/float(d)
!            ci2(i) = 1.0D0/float(i**2)
!            ci3(i) = exp(i*log(10.0D0**(-8))/d)
!          end do
!          do j = 1,mpts
!             Y(j) = (d + 1 + dot_product(Xin(:,j),ci1))**(-(d+1))
!          end do
         coeffs = (/ 0.75127, 0.25510, 0.50596, 0.69908, 0.89090 /)
         do j = 1,mpts
            Y(j) = (1.0D0 + dot_product(Xin(:,j),coeffs))**(-(d+1))
         end do
         
      case ('polyn2')
         ! Second-order Polynomial
         coeffs = (/ 0.75127, 0.25510, 0.50596, 0.69908, 0.89090 /)
         do k = 1,mpts
            partsum = 1.0D0 + dot_product(Xin(:,k),coeffs)
            do i = 1,d-1
              do j = i+1,d
                partsum = partsum + coeffs(i)*coeffs(j)*Xin(i,k)*Xin(j,k)
              end do
            end do  
            Y(k) = partsum
         end do
         
      case ('oscil')
         ! Positive oscillatory function
         omega = 0.049654
         coeffs = (/ 0.90272, 0.94479, 0.49086, 0.48925, 0.33772 /)
         do j = 1,mpts
            Y(j) = 2.0D0 + cos(2*pi*omega + dot_product(Xin(:,j),coeffs))
         end do   
         
      case ('rational')
         ! Rational Polynomial
         coeffs = (/ 0.75127, 0.25510, 0.50596, 0.69908, 0.89090 /)
         do k = 1,mpts
            partsum = 1.0D0 + dot_product(Xin(:,k),coeffs)
            do i = 1,d-1
              do j = i+1,d
                partsum = partsum + coeffs(i)*coeffs(j)*Xin(i,k)*Xin(j,k)
              end do
            end do  
            Y(k) = 1.0D0/partsum
         end do	
         
       case ('Gaussian')
          ! Gaussian function
          coeffs = (/ 0.75127, 0.25510, 0.50596, 0.69908, 0.89090 /)
          xi = (/ -0.68477, 0.94118, 0.91433, -0.02924, 0.60056 /)
          do k = 1,mpts
            partsum = 0.0D0
            do j = 1,d
              partsum = partsum + (coeffs(j)**2)*(Xin(j,k)-xi(j))**2
            end do
            Y(k) = exp(-partsum)
          end do
          
       case ('contin')
         ! Continuous function
         coeffs = (/ 0.75127, 0.25510, 0.50596, 0.69908, 0.89090 /)
         xi = (/ -0.68477, 0.94118, 0.91433, -0.02924, 0.60056 /)
         do k = 1,mpts
            partsum = 0.0D0
            do j = 1,d
              partsum = partsum + coeffs(j)*abs(Xin(j,k)-xi(j))
            end do
            Y(k) = exp(-partsum)
        end do
        
       case ('Friedman2')
         ! Friedman2 function
         do j = 1,mpts
            Y(j) = sqrt(Xin(1,j)**2 + (Xin(2,j)*Xin(3,j)-1.0D0/(Xin(2,j)*Xin(4,j)))**2) 
         end do
         
       case ('Friedman3')
         ! Friedman3 function  
         do j = 1,mpts
            Y(j) = atan((Xin(2,j)*Xin(3,j)-1.0D0/(Xin(2,j)*Xin(4,j)))/Xin(1,j))
         end do
        
       case('curved')
         ! curved function from Dette & Pepelyshev 2010
         do j = 1,mpts
           Y(j) = 4.0D0*(Xin(1,j)-2.0D0+8.0D0*Xin(2,j)-8.0D0*Xin(2,j)**2)**2 &
                  + (3.0D0-4.0D0*Xin(2,j))**2 + 16.0D0*sqrt(Xin(3,j)+1.0D0)*(2.0D0*Xin(3,j)-1)**2
         end do
         
   end select   
   
   end subroutine evalfunc
   !==========================================================================================!

   !==========================================================================================!
   !
   !Computation of discrete semi-norm ||u-ur||_D.
   ! Check if post-processing of predicted values is needed or not.
   !
   subroutine err_normD(Ydata, d, mpts, rank, basissize, varphi, pccoeffs, lbda, err, fobj, & 
                        normD_ur, urdata, int_resp)
      
   !inputs:
   integer, intent(in)                               :: d, mpts, rank, basissize
   real *8, intent(in)                               :: lbda
   real *8, dimension(mpts), intent(in)              :: Ydata
   real *8, dimension(mpts,basissize,d), intent(in)  :: varphi
   real *8, dimension(rank,d,basissize), intent(in)  :: pccoeffs
   character(len=4), intent(in)                      :: int_resp
   
   !outputs:
   real *8, intent(out)                              :: err, fobj, normD_ur
   real *8, dimension(mpts), intent(out)             :: urdata
   
   !local variables:
   integer                                           :: i, j, l
   real *8                                           :: ur, ur_l       
      
   ! Compute approximate function at each data point
   do i = 1,mpts
      ur = 0.0d0
      do l = 1,rank !rank loop
         ur_l = 1.0d0
         do j = 1,d ! loop on dimension
            ur_l = ur_l * dot_product(pccoeffs(l,j,:),varphi(i,:,j))
         end do
         ur = ur + ur_l !add contribution of rank level l
      end do 
      urdata(i) = ur !store approximate value
   end do
   ! Post-processing of predicted values (only needed for ABALONE, REDWINE datasets):
   if (int_resp == 'y') then
      urdata = NINT(urdata) !Take the nearest integers (discrete output)
   end if
   
   err = sqrt((sum((Ydata-urdata)**2))/mpts) !discrete error norm:
   normD_ur = sqrt((sum((urdata)**2))/mpts)
  
   !objective function (with l1 penalization):
   fobj = err**2 + lbda * sum(abs(pccoeffs))
  
   end subroutine err_normD
  !==========================================================================================!

  
   !==========================================================================================!
   !
   !Computation of discrete semi-norm ||u-ur||_D with bias term in ur.
   !
   subroutine err_normD_with_bias(Ydata, d, mpts, rank, basissize, varphi, coeffs, lbda, &
                                  err, fobj, normD_ur, urdata, bias, int_resp)
      
   !inputs:
   integer, intent(in)                               :: d, mpts, rank, basissize
   real *8, intent(in)                               :: lbda, bias
   real *8, dimension(mpts), intent(in)              :: Ydata
   real *8, dimension(mpts,basissize,d), intent(in)  :: varphi
   real *8, dimension(rank,d,basissize), intent(in)  :: coeffs
   character(len=4), intent(in)                      :: int_resp
   
   !outputs:
   real *8, intent(out)                              :: err, fobj, normD_ur
   real *8, dimension(mpts), intent(out)             :: urdata
   
   !local variables:
   integer                                           :: i, j, l
   real *8                                           :: ur, ur_l       
      
   ! Compute approximate function at each data point
   do i = 1,mpts
      ur = 0.0d0
      do l = 1,rank !rank loop
         ur_l = 1.0d0
         do j = 1,d ! loop on dimension
            ur_l = ur_l * dot_product(coeffs(l,j,:),varphi(i,:,j))
         end do
         ur = ur + ur_l !add contribution of rank level l
      end do 
!       print*, 'ur=', ur
      urdata(i) = ur + bias !store approximate value
   end do
   ! Post-processing of predicted values (e.g., ABALONE, REDWINE datasets):
   if (int_resp == 'y') then
      urdata = NINT(urdata) !Take the nearest integers (discrete output)
   end if
  
   err = sqrt((sum((Ydata-urdata)**2))/mpts) !discrete error norm:
   normD_ur = sqrt((sum((urdata)**2))/mpts)
  
   !objective function (with l1 penalization):
   fobj = err**2 + lbda * sum(abs(coeffs))
   
!    print*, 'max(Ydata)=', maxval(Ydata)
!    print*, 'max(urdata)=', maxval(urdata)
!    print*, 'max(|coeffs|)=', maxval(coeffs)
  
   end subroutine err_normD_with_bias
  !==========================================================================================!
  
  !==========================================================================================!
  !
  !Computation of discrete semi-norm ||u-ur||_D in ALS algorithm. Slightly different version
  ! compared to 'err_normD' subroutine since here we have weights 'sconsts' for each level 'l'. 
  !
  subroutine err_normD_ALS(Ydata, sconsts, coeffs, rank, Mbasis, varphi, d, mpts, err, normD_ur, &
                           urdata, int_resp)
   
  !inputs:
  integer, intent(in)                            :: rank, Mbasis, d, mpts
  real *8, dimension(mpts), intent(in)           :: Ydata
  real *8, dimension(rank), intent(in)           :: sconsts
  real *8, dimension(rank,d,Mbasis), intent(in)  :: coeffs
  real *8, dimension(mpts,Mbasis,d), intent(in)  :: varphi
  character(len=4), intent(in)                   :: int_resp
  
  !outputs:
  real *8, intent(out)                           :: err, normD_ur
  real *8, dimension(mpts), intent(out)          :: urdata
  !local:
  integer                                        :: i, k, l
  real *8                                        :: ur, ur_l
    
  !Compute the approximate solution at each data point:
  do k = 1,mpts
     ur = 0.0D0
      do l = 1,rank
	 ur_l = 1.0D0
         do i = 1,d
           !PC expansion in dimension i:
	    ur_l = ur_l*dot_product(varphi(k,:,i),coeffs(l,i,:)) 
         end do 
         ur_l = sconsts(l)*ur_l
         ur = ur + ur_l !add contribution of rank level l
      end do !rank loop
      urdata(k) = ur !store approximate value in a vector
  end do !data point loop
  
  ! Post-processing of predicted values (only needed for ABALONE, REDWINE datasets):
  if (int_resp == 'y') then
     urdata = NINT(urdata) !Take the nearest integers (discrete output)
  end if
  
  err = sqrt((sum((Ydata-urdata)**2))/mpts) !discrete error norm
  normD_ur = sqrt((sum((urdata)**2))/mpts)
   
  end subroutine err_normD_ALS
  !==========================================================================================!     
  
  !==========================================================================================!   
  !
  ! Computation of ||ui^l||_D, where i refers to the i-th dimension and l denotes the rank level.
  ! Useful for ALS_mod.f90/classicalALS()
  !
  subroutine normD_1dfactor(coeffs, i, l, mpts, Mbasis, rank, d, varphi, normD_uil)

  !inputs:
  integer, intent(in)                            :: i, l, mpts, Mbasis, rank, d
  real *8, dimension(rank,d,Mbasis), intent(in)  :: coeffs
  real *8, dimension(mpts,Mbasis,d), intent(in)  :: varphi
  !outputs:
  real *8, intent(out)                           :: normD_uil
  !local:
  integer                                        :: k
  real *8, dimension(mpts)                       :: uil
   
  do k = 1,mpts !loop over data points
     uil(k) = dot_product(varphi(k,:,i),coeffs(l,i,:)) 
  end do

  normD_uil = sqrt((sum(uil**2))/mpts)

  end subroutine normD_1dfactor
  !==========================================================================================!   

  !==========================================================================================!   
  !
  ! Compute the trace of a matrix.
  ! 
  subroutine trace(A, n, traceA)
   
  !inputs:
  integer, intent(in)                 :: n
  real *8, dimension(n,n), intent(in) :: A
  !outputs:
  real *8, intent(out)                :: traceA
  !local
  integer                             :: i
  real *8                             :: tmp
   
  tmp = 0.0D0
  do i=1,n
     tmp = tmp + A(i,i)
  end do
  traceA = tmp
   
  end subroutine trace
  !==========================================================================================!   
  
  !==========================================================================================!   
  !
  ! Evaluate RBF at a given vector of points.
  ! 
  subroutine evalrbf(xpts, dimx, rbfcenter, rbfsigma, y)
   
  !inputs:
  integer, intent(in)                    :: dimx
  real *8, dimension(dimx), intent(in)   :: xpts
  real *8, intent(in)                    :: rbfcenter, rbfsigma
  !outputs:
  real *8, dimension(dimx), intent(out)  :: y
  
  y = exp(-(xpts-rbfcenter)**2/(rbfsigma**2))
  
  end subroutine evalrbf
  !==========================================================================================!   
   
  !==========================================================================================!   
  !
  ! Compute one-dimensional functions using Legendre basis functions.
  ! We assume here that the input variables belong to the SAME intervals [xmin,xmax]
  ! 
  subroutine compute_1dfactors_Legendre(rank, d, basissize, coeffs_mat, ngrid, xmin, xmax, &
   factors1d, xgrid)
  
  !inputs:
  integer, intent(in)                              :: rank, d, basissize, ngrid
  real *8, dimension(rank,d,basissize), intent(in) :: coeffs_mat
  real *8, intent(in)                              :: xmin, xmax
  
  !outputs:
  real *8, dimension(rank,d,ngrid), intent(out)    :: factors1d
  real *8, dimension(ngrid), intent(out)           :: xgrid
  
  !local variables:
  integer                                          :: i, j, k, l, n, order
  real *8, dimension(ngrid)                        :: PCbasis
  real *8, dimension(ngrid,basissize)              :: varphi_grid
  real *8                                          :: dx, normPC
  
  dx = (xmax-xmin)/(float(ngrid)-1.0D0)
  do i = 1,ngrid
    xgrid(i) = xmin + (i-1)*dx
  end do
 
  !Evaluate the Legendre polynomials on a regular grid:
  do k = 1,basissize
     order = k-1
     call lp_value (ngrid, order, xgrid, PCbasis) !one-dimensional Legendre polynomials
     normPC = sqrt(1.0D0/(2.0D0*order+1.0D0))
     varphi_grid(:,k) = PCbasis(:)/normPC
  end do 
    
  !DEBUG:
!   do i=1,ngrid
!      print*, 'PC1=', varphi_grid(i,1)
!   end do
!   do i=1,ngrid
!      print*, 'PC2=', varphi_grid(i,2)
!   end do
!   do i=1,ngrid
!      print*, 'PC3=', varphi_grid(i,3)
!   end do
!
!   do i = 1,basissize
!     print*,'coeffs=', coeffs_mat(1,1,i), coeffs_mat(1,2,i), coeffs_mat(1,3,i), coeffs_mat(1,4,i), coeffs_mat(1,5,i)
!   end do
!   
!   do i = 1,ngrid
!     print*,'varphi_grid=', varphi_grid(i,:)
!   end do  
  !END OF DEBUG  
  
  !Compute one-dimensional factors on regular grid points:
  do l = 1,rank
    do j = 1,d
      factors1d(l,j,:) = matmul(varphi_grid,coeffs_mat(l,j,1:basissize))
!       do n = 1,ngrid
!         factors1d(l,j,n) = dot_product(coeffs_mat(l,j,1:basissize),varphi_grid(n,1:basissize))
!       end do 
    end do
  end do
  
  open (50, FILE='./Data/coeffs.dat', STATUS='UNKNOWN', FORM='FORMATTED') 
  do l = 1,rank
   do j = 1,d
    do k = 1,basissize
      write (50, *) coeffs_mat(l,j,k)
    end do
   end do  
  end do
  close (50)
  
  open (51, FILE='./Data/varphi_grid.dat', STATUS='UNKNOWN', FORM='FORMATTED') 
  do n = 1,ngrid
   do k = 1,basissize
      write (51, *) varphi_grid(n,k)
   end do
  end do  
  close (51)
  
  end subroutine compute_1dfactors_Legendre
  
  !==========================================================================================!   
  !
  ! Compute one-dimensional functions using RBF basis functions.
  ! We assume here that the input variables belong to the SAME intervals [xmin,xmax]
  ! 
  subroutine compute_1dfactors_rbf(rank, d, basissize, coeffs_mat, ngrid, xmin, xmax, & 
   factors1d, xgrid, dir_data)
  
  !inputs:
  integer, intent(in)                              :: rank, d, basissize, ngrid
  real *8, dimension(rank,d,basissize), intent(in) :: coeffs_mat
  real *8, intent(in)                              :: xmin, xmax
  character(len=60), intent(in)                    :: dir_data

  !outputs:
  real *8, dimension(rank,d,ngrid), intent(out)    :: factors1d
  real *8, dimension(ngrid), intent(out)           :: xgrid
  
  !local variables:
  integer                                          :: i, j, k, l, n, order
  real *8, dimension(ngrid)                        :: RBFbasis
  real *8, dimension(ngrid,basissize)              :: varphi_grid
  real *8                                          :: dx, normPC
  real *8, dimension(basissize)                    :: rbfcenters, rbfsigmas
  character(len=80)                                :: file_data
  
  !load RBF data:
  file_data = trim(dir_data) // trim('rbfcenters.dat')
  open (20, FILE=file_data, STATUS='OLD', FORM='FORMATTED') 
  do i = 1,basissize
    read (20, *) rbfcenters(i)
  enddo
  close(20)
  file_data = trim(dir_data) // trim('rbfsigmas.dat')
  open (21, FILE=file_data, STATUS='OLD', FORM='FORMATTED') 
  do i = 1,basissize
    read (21, *) rbfsigmas(i)
  enddo
  close(21)
  
  dx = (xmax-xmin)/(float(ngrid)-1.0D0)
  do i = 1,ngrid
    xgrid(i) = xmin + (i-1)*dx
  end do
 
  !Evaluate the RBFs on a regular grid:
  do k = 1,basissize
     call evalrbf(xgrid, ngrid, rbfcenters(k), rbfsigmas(k), RBFbasis)
     varphi_grid(:,k) = RBFbasis
  end do 
  
  !Compute one-dimensional CD factors on regular grid points:
  do l = 1,rank
    do j = 1,d
      do n = 1,ngrid
        factors1d(l,j,n) = dot_product(coeffs_mat(l,j,:),varphi_grid(n,:))
      end do 
    end do
  end do
  
  end subroutine compute_1dfactors_rbf
   
  !==========================================================================================!
  !
  ! Compute the L^2-norm of the gradient of the objective functional used in SRCD.
  !
  ! Evaluation of the gradient of the objective functional: see the Fortran 77 routine used
  ! in LBFGS (/Algos/minimiz_routines.f)
  !
  subroutine eval_NormGradObj(N, coeffs, normGradObj, normGradObj_block, N_train, rank, d, & 
   basissize, lbda, varphi_train, Ydata_train)  
  
  !inputs:
  integer, intent(in)                                 :: N, N_train, rank, d, basissize
  real *8, intent(in)                                 :: lbda
  real *8, dimension(rank,d,basissize), intent(in)    :: coeffs
  real *8, dimension(N_train,basissize,d), intent(in) :: varphi_train
  real *8, dimension(N_train), intent(in)             :: Ydata_train
      
  !outputs:
  real *8, intent(out)                                :: normGradObj
  real *8, dimension(d), intent(out)                  :: normGradObj_block
  
  !Local variables:
  integer                                             :: i, j, k, l, jj, indg
  real *8, dimension(N_train,rank,d,basissize)        :: partial_urdata
  real *8, dimension(N_train)                         :: urdata
  real *8, dimension(N)                               :: gradF
  real *8                                             :: tmp, ur, ur_l, tmp_part, sgn, Ymean, &
                                                         urdata_mean, grad_block
          
  ! Remark: N represents the number of variables, that is N = rank*d*basissize
      
  Ymean = sum(Ydata_train)/N_train !averaged output
   
  ! Compute the approximate function at each training data point:
  do i = 1,N_train
    ur = 0.0d0
    do l = 1,rank !rank loop
      ur_l = 1.0d0
      do j = 1,d ! loop on dimension
        call dot_prod(coeffs(l,j,1:basissize), varphi_train(i,1:basissize,j), basissize, tmp)
        ur_l = ur_l * tmp
      end do
      ur = ur + ur_l !add contribution of rank level l
    end do 
    urdata(i) = ur !store \tilde{g} values
  end do  
  urdata_mean = sum(urdata)/N_train
  ! Remove the mean value to each of these approximations (coming from bias term):
  do i = 1,N_train
    urdata(i) = urdata(i) - urdata_mean
  end do
      
  ! Compute partial derivatives of \tilde{g} at each training point:
  do i = 1,N_train
    do l = 1,rank
      do j = 1,d
        do k = 1,basissize
          tmp_part = 1.0d0
          do jj = 1,d
            if (jj /= j) then
              call dot_prod(coeffs(l,jj,1:basissize), varphi_train(i,1:basissize,jj), basissize, tmp)
              tmp_part = tmp_part * tmp
            end if
          end do !jj
          partial_urdata(i,l,j,k) = tmp_part * varphi_train(i,k,j)
        end do !k
      end do !j
    end do !l
  end do !i
  ! Remove the derivative of the mean to each of these derivatives (coming from bias term):
  do i = 1,N_train
    do l = 1,rank
      do j = 1,d
        do k = 1,basissize
          partial_urdata(i,l,j,k) = partial_urdata(i,l,j,k) - (sum(partial_urdata(:,l,j,k)))/N_train
        end do
      end do
    end do
  end do
  
  ! Compute the gradient of the objective function
  do l = 1,rank
    do j = 1,d
      do k = 1,basissize
        tmp = 0.0d0
        do i = 1,N_train
          tmp = tmp + (Ydata_train(i) - Ymean - urdata(i)) * partial_urdata(i,l,j,k)  
        end do !i
        indg = (l-1)*d*basissize + (j-1)*basissize + k
        call absprime(coeffs(l,j,k), sgn)
        gradF(indg) = -tmp/N_train + lbda * sgn
      end do !k
    end do !j
  end do !l
  
  ! Compute the (normalized) L2-norm of the gradient of F:
  normGradObj = sqrt(sum(gradF**2)/N) 
  
  ! Compute the (normalized) L2-norm of the gradient of F w.r.t each block of variables:
  normGradObj_block(:) = 0.0d0
  do j = 1,d
    grad_block = 0.0d0
    do l = 1,rank
      do k = 1,basissize
        indg = (l-1)*d*basissize + (j-1)*basissize + k
        grad_block = grad_block + gradF(indg)**2
      end do
    end do
    normGradObj_block(j) = sqrt(grad_block/(rank*basissize))
  end do
  
  end subroutine eval_NormGradObj
  
  !==========================================================================================!
  !
  ! Computes the dot product of arrays a and b of n elements 
  !
  subroutine dot_prod(a, b, n, z)
  
  !inputs:
  integer, intent(in)               ::  n
  real *8, dimension(n), intent(in) ::  a, b
  
  !output:
  real *8, intent(out)              :: z
  
  !Local variables:
  integer                           ::  i
      
  z = 0.0d0
  do i = 1,n
    z = z + a(i) * b(i) 
  end do
  
  end subroutine dot_prod
  
  !==========================================================================================!
  !
  ! Computes the derivative of absolute value of a real, i.e. the sign of this real.
  !      
  subroutine absprime(x, sgn)
      
  !input:
  real *8, intent(in)  :: x
  
  !output:
  real *8, intent(out)  :: sgn
      
  if (x > 0.0d0) then
    sgn = 1.0d0
  else if (x < 0.0d0) then
    sgn = -1.0d0
  else if (x == 0.0d0) then
    sgn = 0.0d0
  end if
      
  end subroutine absprime
  
  !==========================================================================================!
  !
  !  Evaluates the basis functions on the training/testing datasets.
  !      
  subroutine eval_basisfunc_train_test(nobs, basissize, d, opt_basis, Xdata_work, & 
   N_train, N_test, varphi_train, varphi_test, dir_data)
  
  !inputs:
  integer, intent(in)                      :: nobs, basissize, d, opt_basis, N_train, N_test
  real *8, dimension(d,nobs), intent(in)   :: Xdata_work
  character(len=60), intent(in)            :: dir_data
   
  !outputs:
  real *8, dimension(N_train,basissize,d), intent(out) :: varphi_train
  real *8, dimension(N_test,basissize,d), intent(out)  :: varphi_test
  
  !local variables:
  integer                                  :: i, j, k, order(1)
  real *8, dimension(:,:,:), allocatable   :: varphi_work
  real *8, dimension(:), allocatable       :: xpts, PCbasis, RBFbasis, rbfcenters, rbfsigmas
  real *8                                  :: normPC, step_rbf, sigmaval, xmin, xmax
  character(len=80)                        :: file_data
  
  allocate(varphi_work(nobs,basissize,d))

  if (opt_basis == 0) then 
    !
    ! PC Legendre basis functions
    !
    allocate(xpts(d),PCbasis(d))
    do j = 1,nobs
       xpts = Xdata_work(:,j)
       do k = 1,basissize
         order(1) = k-1
         call lpp_value(1, d, order, xpts, PCbasis)
         normPC = sqrt(1.0D0/(2.0D0*order(1)+1.0D0))
         varphi_work(j,k,:) = PCbasis(:)/normPC
       end do 
    end do
    deallocate(xpts, PCbasis)
    !
  elseif (opt_basis == 1) then 
    !
    ! RBF basis functions
    !
    !Define RBF centers and (unique set of) shape parameters:
    allocate(rbfcenters(basissize),rbfsigmas(basissize))
    xmin = 0.0D0
    xmax = 1.0D0
    step_rbf = (xmax-xmin)/(float(basissize)+1.0D0)
    sigmaval = 0.9D0 !larger sigma (flatter RBFs) => should lead to more sparsity
    print*, 'RBF shape parameter=', sigmaval
    do i = 1,basissize
      rbfcenters(i) = xmin + i*step_rbf
      rbfsigmas(i) = sigmaval
    end do
    !
    !Define multiple RBFs (i.e., with different shape parameters)
    !(in this example take basissize=20):
!   step_rbf = (xmax-xmin)/(11.0D0)
!   do i = 1,10
!     rbfcenters(i) = xmin + i*step_rbf
!     rbfsigmas(i) = 0.1D0
!   end do
!   do i = 1,10
!     rbfcenters(i+10) = xmin + i*step_rbf
!     rbfsigmas(i+10) = 0.25D0
!   end do
    !
    allocate(xpts(d),RBFbasis(d))
    do j = 1,nobs
      xpts = Xdata_work(:,j)
      do k = 1,basissize
        call evalrbf(xpts, d, rbfcenters(k), rbfsigmas(k), RBFbasis)
        varphi_work(j,k,:) = RBFbasis(:)   
      end do
    end do
    deallocate(xpts,RBFbasis)
    
!     print*, 'rbfcenters=', rbfcenters(:)
!     print*, ''
!     print*, 'rbfsigmas=', rbfsigmas(:)
    
    !save RBF data:
    file_data = trim(dir_data) // trim('rbfcenters.dat')
    open (20, FILE=file_data, STATUS='UNKNOWN', FORM='FORMATTED') 
    do i = 1,basissize
      write (20, *) rbfcenters(i)
    enddo
    close(20)
    file_data = trim(dir_data) // trim('rbfsigmas.dat')
    open (21, FILE=file_data, STATUS='UNKNOWN', FORM='FORMATTED') 
    do i = 1,basissize
      write (21, *) rbfsigmas(i)
    enddo
    close(21)
    
    deallocate(rbfcenters, rbfsigmas)
    !
  end if
  
  ! Evaluate basis functions on training dataset:
  varphi_train = varphi_work(1:N_train,:,:)
  ! Evaluate basis functions on testing dataset:
  varphi_test = varphi_work(N_train+1:nobs,:,:)
  
  deallocate(varphi_work)
  
  end subroutine eval_basisfunc_train_test
  
  !==========================================================================================!
  !
  !  Evaluates the basis functions on the WHOLE dataset.
  !      
  subroutine eval_basisfunc_work(nobs, basissize, d, opt_basis, Xdata_work, varphi_work, sigmaval)
  
  !inputs:
  integer, intent(in)                      :: nobs, basissize, d, opt_basis
  real *8, dimension(d,nobs), intent(in)   :: Xdata_work
  real *8, intent(in)                      :: sigmaval
   
  !outputs:
  real *8, dimension(nobs,basissize,d), intent(out) :: varphi_work
  
  !local variables:
  integer                                  :: i, j, k, order(1)
  real *8, dimension(:), allocatable       :: xpts, PCbasis, RBFbasis, rbfcenters, rbfsigmas
  real *8                                  :: normPC, step_rbf, xmin, xmax
  
  if (opt_basis == 0) then 
    !
    ! PC Legendre basis functions
    !
    allocate(xpts(d),PCbasis(d))
    do j = 1,nobs
       xpts = Xdata_work(:,j)
       do k = 1,basissize
         order(1) = k-1
         call lpp_value(1, d, order, xpts, PCbasis)
         normPC = sqrt(1.0D0/(2.0D0*order(1)+1.0D0))
         varphi_work(j,k,:) = PCbasis(:)/normPC
       end do 
    end do
    deallocate(xpts, PCbasis)
    !
  elseif (opt_basis == 1) then 
    !
    ! RBF basis functions
    !
!     print*, 'RBF shape parameter=', sigmaval
    !Define RBF centers and (unique set of) shape parameters:
    allocate(rbfcenters(basissize),rbfsigmas(basissize))
    xmin = 0.0D0
    xmax = 1.0D0
    step_rbf = (xmax-xmin)/(float(basissize)+1.0D0)
    do i = 1,basissize
      rbfcenters(i) = xmin + i*step_rbf
      rbfsigmas(i) = sigmaval
    end do
    !
    !Define multiple RBFs (i.e., with different shape parameters)
    !(in this example take basissize=20):
!   step_rbf = (xmax-xmin)/(11.0D0)
!   do i = 1,10
!     rbfcenters(i) = xmin + i*step_rbf
!     rbfsigmas(i) = 0.1D0
!   end do
!   do i = 1,10
!     rbfcenters(i+10) = xmin + i*step_rbf
!     rbfsigmas(i+10) = 0.25D0
!   end do
    !
    allocate(xpts(d),RBFbasis(d))
    do j = 1,nobs
      xpts = Xdata_work(:,j)
      do k = 1,basissize
        call evalrbf(xpts, d, rbfcenters(k), rbfsigmas(k), RBFbasis)
        varphi_work(j,k,:) = RBFbasis(:)   
      end do
    end do
    deallocate(xpts,RBFbasis)
    
!     print*, 'rbfcenters=', rbfcenters(:)
!     print*, ''
!     print*, 'rbfsigmas=', rbfsigmas(:)
    
    deallocate(rbfcenters, rbfsigmas)
    !
  end if
  
  end subroutine eval_basisfunc_work
  
  !==========================================================================================!
  !
  !  Creates the dataset directory and sets up some parameters.
  !     nobs: total number of observations
  !     nfeat: number of features (inputs + single output)
  !     int_resp: continuous ('n') or discrete ('y') output
  !     dir_data: directory of dataset    
  subroutine setup_dataset(opt_data, nobs, nfeat, int_resp, dir_data)
  
  !inputs:
  character(len=10), intent(in)                        :: opt_data
  
  !outputs:
  integer, intent(out)                                 :: nobs, nfeat
  character(len=4), intent(out)                        :: int_resp
  character(len=60), intent(out)                       :: dir_data
    
  select case(opt_data)
  case('challenger')
    nobs = 23
    nfeat = 5
    int_resp = 'n'
    dir_data = './Data/UCI_datasets/challenger/'
  case('fertility')
    nobs = 100
    nfeat = 10
    int_resp = 'n'
    dir_data = './Data/UCI_datasets/fertility/'
  case('servo')
    nobs = 167
    nfeat = 5
    int_resp = 'n'
    dir_data = './Data/UCI_datasets/servo/'
  case('yacht')
    nobs = 308
    nfeat = 7
    int_resp = 'n'
    dir_data = './Data/UCI_datasets/yacht/'
  case('autompg')
    nobs = 392
    nfeat = 8
    int_resp = 'n'
    dir_data = './Data/UCI_datasets/autoMPG/'
  case('housing')
    nobs = 506
    nfeat = 14
    int_resp = 'n'
    dir_data = './Data/UCI_datasets/housing/'  
  case('forest')
    nobs = 517
    nfeat = 13
    int_resp = 'n'
    dir_data = './Data/UCI_datasets/forest/'
  case('stock')
    nobs = 536
    nfeat = 12
    int_resp = 'n'
    dir_data = './Data/UCI_datasets/stock/'
  case('pendulum')
    nobs = 630
    nfeat = 10
    int_resp = 'n'
    dir_data = './Data/UCI_datasets/pendulum/'
  case('energy')
    nobs = 768
    nfeat = 9
    int_resp = 'n'
    dir_data = './Data/UCI_datasets/energy/' 
  case('concrete')
    nobs = 1030
    nfeat = 9
    int_resp = 'n'
    dir_data = './Data/UCI_datasets/concrete/'
  case('solar')
    nobs = 1066
    nfeat = 10
    int_resp = 'n'
    dir_data = './Data/UCI_datasets/solar/'
  case('airfoil')
    nobs = 1503
    nfeat = 6
    int_resp = 'n'
    dir_data = './Data/UCI_datasets/airfoil/'
  case('wine')
    nobs = 1599
    nfeat = 12
    int_resp = 'n'
    dir_data = './Data/UCI_datasets/wine/'
  case('gas')
    nobs = 2565
    nfeat = 129
    int_resp = 'n'
    dir_data = './Data/UCI_datasets/gassensor/' 
  case('skillcraft')
    nobs = 3338
    nfeat = 20
    int_resp = 'n'
    dir_data = './Data/UCI_datasets/skillcraft/'
  case('sml')
    nobs = 4137
    nfeat = 27
    int_resp = 'n'
    dir_data = './Data/UCI_datasets/sml/'   
  case('parkinsons')
    nobs = 5875
    nfeat = 21
    int_resp = 'n'
    dir_data = './Data/UCI_datasets/parkinsons/'   
  case('pumadyn')
    nobs = 8192
    nfeat = 33
    int_resp = 'n'
    dir_data = './Data/UCI_datasets/pumadyn/' 
  case('pol')
    nobs = 15000
    nfeat = 27
    int_resp = 'n'
    dir_data = './Data/UCI_datasets/poltele/'   
  case('elevators')
    nobs = 16599
    nfeat = 19
    int_resp = 'n'
    dir_data = './Data/UCI_datasets/elevators/'  
  case('kin40k')
    nobs = 40000
    nfeat = 9
    int_resp = 'n'
    dir_data = './Data/UCI_datasets/kin40k/'  
  case('protein')
    nobs = 45730
    nfeat = 10
    int_resp = 'n'
    dir_data = './Data/UCI_datasets/protein/' 
  case('kegg')
    nobs = 48827
    nfeat = 21 !to be double-checked
    int_resp = 'n'
    dir_data = './Data/UCI_datasets/kegg/'
  case('slice')
    nobs = 53500
    nfeat = 386
    int_resp = 'n'
    dir_data = './Data/UCI_datasets/slice/'
  case('keggu')
    nobs = 63608
    nfeat = 28
    int_resp = 'n'
    dir_data = './Data/UCI_datasets/keggundir/'
  case default
    print*, 'opt_data is not correct!'
    print*, 'STOP.'
    STOP
  end select
  
  end subroutine setup_dataset
  
  !
  !  Creates dataset directory and set up some parameters (synthetic test-cases)
  !   dir_data: directory of dataset  
  !   d: input dimensionality
  !   int_resp: continuous ('n') or discrete ('y') output
  subroutine setup_dir_synth(opt_func, dir_data, d, int_resp)
  
  !inputs:
  character(len=10), intent(in)        :: opt_func
  
  !outputs:
  character(len=60), intent(out)       :: dir_data
  integer, intent(out)                 :: d
  character(len=4), intent(out)        :: int_resp
    
  select case(opt_func)
  case('Friedman1')
    dir_data = './Data/synthetic_datasets/Friedman1/'
    d = 5
    int_resp = 'n'
  case('Friedman2')
    dir_data = './Data/synthetic_datasets/Friedman2/'
    d = 4
    int_resp = 'n'
  case('Friedman3')
    dir_data = './Data/synthetic_datasets/Friedman3/'
    d = 4
    int_resp = 'n'
  case('Rastrigin')
    dir_data = './Data/synthetic_datasets/Rastrigin/'
    d = 5
    int_resp = 'n'
  case('corner')
    dir_data = './Data/synthetic_datasets/corner/'
    d = 5
    int_resp = 'n'
  case('oscil')
    dir_data = './Data/synthetic_datasets/oscil/'
    d = 5
    int_resp = 'n'  
  case('contin')
    dir_data = './Data/synthetic_datasets/contin/'
    d = 5
    int_resp = 'n'
  case('curved')
    dir_data = './Data/synthetic_datasets/curved/'
    d = 3
    int_resp = 'n'  
  case default
    print*, 'setup_dir_synth(): opt_func is not correct!'
    print*, 'STOP.'
    STOP
  end select
  
  end subroutine setup_dir_synth
  
  !
  ! Converts an integer to string.
  !
  character(len=20) function str(k)
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
  end function str
  
  !==========================================================================================!
  !
  ! Remaps the indices of compressed coefficients returned by glmnet/elnet() into pairs of
  ! indices (lr,lb) which can be used in the 3d data structure of BSR algorithm (see BSR_subpbs_glmnet()).
  ! lr: level rank indices (first row of inds_out)
  ! lb: level basis indices (second row of inds_out)
  subroutine remap_indices(nin, ia, basissize, inds_out)
  
  !inputs:
  integer, intent(in)                     :: nin, basissize
  integer, dimension(nin), intent(in)     :: ia
  
  !outputs:
  integer, dimension(2,nin), intent(out)  :: inds_out
  
  !local variables:
  integer                                 :: i, ind, lr, lb
  
  inds_out(:,:) = -111
  do i = 1,nin
    ind = ia(i)
    lb = MOD(ind,basissize)
    if (lb == 0) then
      inds_out(2,i) = basissize
    else
      inds_out(2,i) = lb
    end if
    lr = ceiling(real(ind)/real(basissize))
    inds_out(1,i) = lr
  end do
!   print*, 'inds_l=', inds_out(1,:)
!   print*, 'inds_k=', inds_out(2,:)
  
  end subroutine remap_indices
  
end module utilities
