!
! Numerical studies on real datasets coming from Machine Learning (UCI depository from GitLab)
! using SSR-CD algorithm. Here one single split of the data is used (simplified version of 
! 'main_SSRCD_real.split.f90').
!

PROGRAM main_driver

USE utilities
USE lpp
USE crossValid
USE coordDescent
USE srcd

IMPLICIT NONE

integer                                  :: nobs, nfeat, d, N_train, N_test, rank, opt_basis, &
                                            basissize, orderPol, nb_lbd_try, kfold, iprint, nzcoeffs
integer                                  :: i, j, k, l, ifold, ilbd
character(len=4)                         :: int_resp, opt_resol
character(len=10)                        :: opt_data
character(len=60)                        :: file_res, file_data, dir_data, data_name
real *8, dimension(:,:), allocatable     :: data_orig, Xdata_work, Xdata_train, Xdata_test, stock_err_test
real *8, dimension(:), allocatable       :: Ydata_work, Ydata_train, Ydata_test, xpts, PCbasis, ur_test, &
                                            lbd_grid, err_test_aver
real *8                                  :: xmin, xmax, epsi, err_test, fobj_test, &
                                            lbd_opt, normD_ur_test, normD_u_test, tol_spars, perc_spars, perc_zeros, &
                                            tstep1, tstep2, coeff_max, bias, tol_fobj_cv, tol_fobj_run
real *8, dimension(:,:,:), allocatable   :: varphi_train, varphi_test, varphi_work, coeffs_mat
logical                                  :: flag_cvg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Load input/output variables and define training/testing data structures !
! (for creating the .dat files see an example in create_dat_file.m)       !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! --- Define the number of features and the total number of observations
!      nobs: total number of observations
!      nfeat: number of features (inputs+output)
!      int_resp: continuous ('n') or discrete ('y') output
!      dir_data: directory of dataset
opt_data = 'challenger' 
call setup_dataset(opt_data, nobs, nfeat, int_resp, dir_data)

! --- Define dataset filename:
file_data = trim(dir_data) // trim('dataset.dat')

! --- Load data:
allocate(data_orig(nobs,nfeat))
open (10, FILE=file_data, STATUS='OLD', FORM='FORMATTED') 
read(10,*) ((data_orig(j,i),i=1,nfeat),j=1,nobs)   
close(10)

! --- Define input dimensionality:
d = nfeat-1

! --- Extract data to define X and Y in the right format:
allocate(Xdata_work(d,nobs))
allocate(Ydata_work(nobs))
do i = 1,d
  Xdata_work(i,:) = data_orig(:,i)
end do
Ydata_work(:) = data_orig(:,nfeat)
  
! --- Rescale the input variables between 0 and 1:
do i = 1,d
  xmin = minval(Xdata_work(i,:))
  xmax = maxval(Xdata_work(i,:))
  if (xmin < xmax) then
    Xdata_work(i,:) = 1.0D0/(xmax-xmin)*(Xdata_work(i,:)-xmin)
  elseif (xmin == xmax) then
     !special case when input vars are constant
     Xdata_work(i,:) = 0.0D0
  endif
end do
! Check the range of input variables:
! print*, 'min(X)=', minval(Xdata_work)
! print*, 'max(X)=', maxval(Xdata_work)

! --- Define the size of training/testing datasets:
N_train = ceiling(0.9*nobs) 
! N_train must also be a multiple of kfold:
kfold = 5 !number of folds in cross-validation
if (mod(N_train,kfold) /= 0) then
  N_train = kfold * floor(real(N_train)/real(kfold)) 
end if
N_test = nobs - N_train

allocate(Xdata_train(d, N_train))
allocate(Ydata_train(N_train))
Xdata_train = Xdata_work(:,1:N_train)
Ydata_train = Ydata_work(1:N_train)

allocate(Xdata_test(d, N_test))
allocate(Ydata_test(N_test))
Xdata_test = Xdata_work(:,N_train+1:nobs)
Ydata_test = Ydata_work(N_train+1:nobs)

deallocate(data_orig, Ydata_work)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define the rank, type of discretization and one-dimensional basis size !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

rank = 3 !low-rank decomposition rank

! --- Choice of discretization
!      0 --> Legendre PC basis functions 
!      1 --> RBF basis functions (with RBF centers on a regular grid)
opt_basis = 0

! --- Size of discretization basis
basissize = 2
if (opt_basis == 0) then
  orderPol = basissize-1 !PC order
end if

! --- Create file for writing some outputs (timings, errors, sparsity):
file_res = trim(dir_data) // trim('out_SRCD_r') // trim(str(rank)) // trim('_M') // trim(str(basissize)) // trim('.dat')
open (99, FILE=file_res, STATUS='UNKNOWN', FORM='FORMATTED')
print*, ''
print*, '========================================================'
! select case(opt_data)
!   case default
!     print*, 'opt_data is not correct!'
!     STOP
! end select
data_name = trim(opt_data) // trim(' dataset')
print*, data_name
write (99,*), data_name
write (99,*) ''
write (99,*) '*** Parameters ***'
write (99,*) 'Nb. of points in training dataset =', N_train
write (99,*) 'Nb. of points in testing dataset =', N_test
write (99,*) 'rank =', rank
write (99,*) 'basissize =', basissize
select case(opt_basis)
  case(0)
    write (99,*), 'discretization: Legendre PC'
  case(1)
    write (99,*), 'discretization: RBF'
  case default
    write (99,*), 'opt_basis is not correct!'  
end select
print*, 'dimension                            =', d
print*, 'rank                                 =', rank
print*, 'Number of points in training dataset =', N_train
print*, 'Number of points in testing dataset  =', N_test
print*, 'One-dimensional basis size           =', basissize
select case(opt_basis)
  case(0)
    print*, ' ...discretization: Legendre PC'
  case(1)
    print*, ' ...discretization: RBF'
  case default
    print*, 'opt_basis is not correct!'
    print*, 'STOP.'
    STOP
end select    
print*, '========================================================'
print*, ''

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Evaluate the basis functions on the training/testing datasets      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(varphi_train(N_train,basissize,d))
allocate(varphi_test(N_test,basissize,d))

call eval_basisfunc_train_test(nobs, basissize, d, opt_basis, Xdata_work, & 
 N_train, N_test, varphi_train, varphi_test, dir_data)  
      
deallocate(Xdata_work)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      Apply SRCD algorithm                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! --- Define parameters:
iprint = 1 !1: display messages / -1: no printing
nb_lbd_try = 20 !size of lambda-grid
epsi = 1e-5 !parameter for defining lower bound in lambda-grid
            ! Remark: if 'epsi' is not small enough SRCD can return zero coefficients,
            !         meaning the approximate model is a constant (equal to the bias term).
! Tolerances for saturation of objective function:
tol_fobj_cv = 1e-6 !in cross validation
tol_fobj_run = 1e-9 !in production run
! Choose the type of numerical resolution
opt_resol = 'cd' !

print*, 'nb_lbd_try=', nb_lbd_try
print*, 'epsi=', epsi
print*, 'tol_fobj_cv=', tol_fobj_cv
print*, 'tol_fobj_run=', tol_fobj_run
print*, 'num resol=', opt_resol

write (99,*) 'nb_lbd_try=', nb_lbd_try
write (99,*) 'epsi=', epsi
write (99,*) 'tol_fobj_cv=', tol_fobj_cv
write (99,*) 'tol_fobj_run=', tol_fobj_run
write (99,*) 'num resol=', opt_resol

print*, ''
print*, 'paused, type [enter] to continue'
read (*,*) 

allocate(coeffs_mat(rank,d,basissize))
allocate(lbd_grid(nb_lbd_try), err_test_aver(nb_lbd_try))
allocate(stock_err_test(kfold,nb_lbd_try))

call srcd_resol(coeffs_mat, rank, d, basissize, N_train, varphi_train, Xdata_train, & 
      Ydata_train, nb_lbd_try, epsi, kfold, iprint, opt_resol, int_resp, flag_cvg, lbd_grid, & 
       err_test_aver, stock_err_test, lbd_opt, tstep1, tstep2, bias, tol_fobj_cv, tol_fobj_run)

! --- Save errors and lambda-grids in local directory.
!     For plots use the routine 'plot_errors_SRCD.m'
!
! Validation errors for each fold and try of lambda
file_res = trim(dir_data) // trim('errors_cv_srcd.dat')
open (20, FILE=file_res, STATUS='UNKNOWN', FORM='FORMATTED') 
do ifold = 1,kfold
  do ilbd = 1,nb_lbd_try
    write (20, *) stock_err_test(ifold,ilbd) 
  end do
end do
close (20)

! Averaged validation errors for try of lambda
file_res = trim(dir_data) // trim('errorsAver_cv_srcd.dat')
open (21, FILE=file_res, STATUS='UNKNOWN', FORM='FORMATTED') 
do ilbd = 1,nb_lbd_try
  write (21, *) err_test_aver(ilbd)
end do
close (21)

! lambda-grid 
file_res = trim(dir_data) // trim('lbd_grid_srcd.dat')  
open (22, FILE=file_res, STATUS='UNKNOWN', FORM='FORMATTED') 
do ilbd = 1,nb_lbd_try
  write (22, *) lbd_grid(ilbd)
end do
close (22)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            Compute validation error on testing dataset                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(ur_test(N_test))
if (flag_cvg .eqv. .true.) then     
  call err_normD_with_bias(Ydata_test, d, N_test, rank, basissize, varphi_test, &
        coeffs_mat, lbd_opt, err_test, fobj_test, normD_ur_test, ur_test, bias, int_resp)    
elseif (flag_cvg .eqv. .false.) then
  print*, 'SRCD with lbd_opt diverges!!!'
  print*, 'STOP.'
  STOP
end if
  
! Compare a few predicted values against true ones:
do i = 1,N_test
  print*, 'Ytrue=', Ydata_test(i), 'Yhat=', ur_test(i)
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  Diagnostics (timings & sparsity)                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

normD_u_test = sqrt((sum(Ydata_test**2))/N_test) !norm of exact response on testing dataset

! Compute the percentage of almost zero coefficients (scaled by max value):
tol_spars = 1e-3
 coeff_max = maxval(abs(coeffs_mat))
nzcoeffs = count(abs(coeffs_mat)/coeff_max > tol_spars)
perc_spars = 100.0D0*(rank*d*basissize-nzcoeffs)/(rank*d*basissize)

nzcoeffs = count(abs(coeffs_mat)> 0.0D0)
perc_zeros = 100.0D0*(rank*d*basissize-nzcoeffs)/(rank*d*basissize)

write (99,*) ''
write (99,*) '*** Penalization parameter ***'
write (99,*) 'lbd_opt=', lbd_opt
write (99,*) ''
write (99,*) '*** Testing errors ***'
write (99,*) '||u-ur||_D / ||u||_D=', err_test/normD_u_test
write (99,*) 'RMSE=||u-ur||_D=', err_test
write (99,*) '||u||_D=', normD_u_test
write (99,*) '||ur||_D=', normD_ur_test
write (99,*) ''
write (99,*) '*** Timings ***'
write (99,*) ' step1=', tstep1
write (99,*) ' step2=', tstep2
write (99,*) ' total=', tstep1 + tstep2
write (99,*) ''
write (99,*) '*** Sparsity***'
write (99,*) ' tol_spars=', tol_spars
write (99,*) ' percentage of |coeffs/coeff_max| < tol_spars=', perc_spars
write (99,*) ' percentage of non-zeros=', perc_zeros
close (99)
print*, ''
print*, 'Check sparsity of model structure'
do l = 1,rank
  print*,'* level l=', l
  do j = 1,d
    print*,' - dimension=', j
    do k = 1,basissize
      print*,'     coeffs=', coeffs_mat(l,j,k)
    end do
  end do
end do
print*, 'max(|coeffs|)=', coeff_max
print*, 'bias=', bias
print*, ''
print*, '*** Penalization parameter ***'
print*, 'lbd_opt=', lbd_opt
print*, ''
print*, '*** Testing errors ***'
print*, '||u-ur||_D / ||u||_D=', err_test/normD_u_test
print*, 'RMSE=||u-ur||_D=', err_test
print*, '||u||_D=', normD_u_test
print*, '||ur||_D=', normD_ur_test
print*, ''
print*, '*** timings ***'
print*, ' step1=', tstep1
print*, ' step2=', tstep2
print*, ' total=', tstep1 + tstep2
print*, ''
print*, '*** sparsity***'
print*, ' tol_spars=', tol_spars
print "(a,10f5.2)", ' percentage of |coeffs/coeff_max| < tol_spars=', perc_spars
print "(a,10f5.2)", ' percentage of non-zeros=', perc_zeros
print*, ''

! Remark: When SRCD returns zero coefficients the approximate solution is a constant
!         equal to the bias term (equal to the mean of response):
if (coeff_max == 0.0D0) then
  print*, 'Approximate solution=bias=averaged response:'
  print*, 'mean(Y)=', sum(Ydata_train)/N_train !mean over training points
  print*, 'bias=', bias
endif 

!--- deallocations:
deallocate(Xdata_train, Xdata_test, Ydata_train, Ydata_test, varphi_train, varphi_test, & 
            coeffs_mat, lbd_grid, err_test_aver, stock_err_test, ur_test)
            
print*, ''
print*, 'DONE :-)'

END PROGRAM main_driver