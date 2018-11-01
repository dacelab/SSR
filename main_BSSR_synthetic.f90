!
! Block SSR algorithm applied to synthetic datasets (i.e. built from known functions).
! Subproblems are solved using LBFGS or GLMNET.
!

PROGRAM main_driver

USE sobol
USE utilities
USE lpp
USE crossValid
USE coordDescent
USE bsr

IMPLICIT NONE

integer                                  :: N_tot, d, N_train, N_test, rank, opt_basis, &
                                            basissize, orderPol, nb_lbd_try_c, nb_lbd_try_f, kfold, iprint, &
                                            nzcoeffs, ngrid
integer                                  :: i, j, k, l, m, ifold, ilbd
character(len=4)                         :: int_resp, opt_resol, opt_datagen, opt_fact1D
character(len=10)                        :: opt_func
character(len=60)                        :: file_data, dir_data, data_name
character(len=80)                        :: file_res
real *4, dimension(:,:), allocatable     :: Xdata_work
real *8, dimension(:,:), allocatable     :: data_orig, Xdata_work_dble, Xdata_train, Xdata_test, stock_err_test_c, &
                                            stock_err_test_f
real *8, dimension(:), allocatable       :: Ydata_work, Ydata_train, Ydata_test, xpts, PCbasis, &
                                            lbd_grid_c, lbd_grid_f, err_test_aver_c, err_test_aver_f, ur_test, xgrid
real *8                                  :: xmin, xmax, beta1, beta2, epsi_c, err_test, fobj_test, lbd_star, &
                                            lbd_opt, normD_ur_test, normD_u_test, tol_spars, perc_spars, perc_zeros, &
                                            tstep1, tstep2, tstep3, coeff_max, tol_fobj_cv, tol_fobj_run, ratioVar
real *8, dimension(:,:,:), allocatable   :: varphi_train, varphi_test, varphi_work, coeffs_mat, factors1d
logical                                  :: flag_cvg


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define input parameters and define training/testing data structures !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

opt_func = 'Friedman1' !choice of function
call setup_dir_synth(opt_func, dir_data, d, int_resp) !define data directory and input dimensionality

N_train = 1000  !number of training points (must be divisible by 5 or 10 for cross-validation)
N_test = 10000   !number of testing points
xmin = 0.0D0     !lower bound for variables
xmax = 1.0D0     !upper bound for variables
N_tot = N_train + N_test

opt_datagen = 'n' ! generating/loading option

if (opt_datagen == 'y') then 
  
  print*, 'Generate input/output datasets'
  
  allocate(Xdata_work(d, N_tot),Xdata_work_dble(d, N_tot))
  allocate(Ydata_work(N_tot))
  call i4_sobol_generate(d, N_tot, 0, Xdata_work) !generate values between 0 and 1
  
  ! ========================================= REMARKS =========================================
  ! 1. 'Xdata_work' must be defined as a real*4 when calling 'i4_sobol_generate'
  ! 2. In some cases the X dataset needs to be corrected manually. For example for Friedman2
  !    and Friedman3 the first generated data point (zero vector) is removed and replaced 
  !    since the output response is not defined.
  ! ============================================================================================
  
  Xdata_work_dble = real(Xdata_work,8) !make a copy in double precision
  Xdata_work_dble = (xmax-xmin)*Xdata_work_dble + xmin !rescale data between xmin and xmax
  call evalfunc(Xdata_work_dble, d, N_tot, opt_func, Ydata_work)
  deallocate(Xdata_work)

  ! Save X and Y:
  file_data = trim(dir_data) // trim('Xdataset.dat')
  open (10, FILE=file_data, STATUS='UNKNOWN', FORM='FORMATTED')
  do i = 1,d
    do j = 1,N_tot
      write (10, *) Xdata_work_dble(i,j)
    end do
  end do
  close (10)

  file_data = trim(dir_data) // trim('Ydataset.dat')
  open (11, FILE=file_data, STATUS='UNKNOWN', FORM='FORMATTED')
  do j = 1,N_tot
    write (11, *) Ydata_work(j)
  end do
  close (11)

elseif (opt_datagen == 'n') then
  
  print*, 'Load input/output datasets'

  allocate(Xdata_work_dble(d, N_tot))
  allocate(Ydata_work(N_tot))

  ! Load X and Y:
  file_data = trim(dir_data) // trim('Xdataset.dat')
  open (10, FILE=file_data, STATUS='OLD', FORM='FORMATTED')
  do i = 1,d
    do j = 1,N_tot
      read (10, *) Xdata_work_dble(i,j)
    end do
  end do
  close (10)

  file_data = trim(dir_data) // trim('Ydataset.dat')
  open (11, FILE=file_data, STATUS='OLD', FORM='FORMATTED')
  do j = 1,N_tot
    read (11, *) Ydata_work(j)
  end do
  close (11)

end if !opt_datagen

! Check the range of input variables:
! print*, 'min(X)=', minval(Xdata_work_dble)
! print*, 'max(X)=', maxval(Xdata_work_dble)
! print*, 'min(Y)=', minval(Ydata_work)
! print*, 'max(Y)=', maxval(Ydata_work)
! STOP

! Extract data for creating training/testing datasets:
allocate(Xdata_train(d, N_train))
allocate(Ydata_train(N_train))
allocate(Xdata_test(d, N_test))
allocate(Ydata_test(N_test))
Xdata_train = Xdata_work_dble(:,1:N_train)
Ydata_train = Ydata_work(1:N_train)
Xdata_test = Xdata_work_dble(:,N_train+1:N_tot)
Ydata_test = Ydata_work(N_train+1:N_tot)
deallocate(Ydata_work)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define the rank, type of discretization and one-dimensional basis size !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

rank = 4 !low-rank decomposition rank

! --- Choice of discretization
opt_basis = 0
! 0 --> Legendre PC basis functions 
! 1 --> RBF basis functions (with RBF centers on a regular grid)

! --- Size of discretization basis
basissize = 3
if (opt_basis == 0) then
  orderPol = basissize-1 !PC order
end if

! --- Create file for writing some outputs (timings, errors, sparsity):
select case(opt_basis)
  case(0)
    file_res = trim(dir_data) // trim('out_BSR_r') // trim(str(rank)) // trim('_M') // &
               trim(str(basissize)) // trim('_Legendre.dat')
  case(1)
    file_res = trim(dir_data) // trim('out_BSR_r') // trim(str(rank)) // trim('_M') // & 
               trim(str(basissize)) // trim('_RBF.dat')
end select
open (99, FILE=file_res, STATUS='UNKNOWN', FORM='FORMATTED')
print*, ''
print*, '========================================================'
! select case(opt_func)
!   case default
!     print*, 'opt_func is not correct!'
!     STOP
! end select
data_name = trim(opt_func) // trim(' dataset')
print*, 'Block-SSR algorithm'
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

call eval_basisfunc_train_test(N_tot, basissize, d, opt_basis, Xdata_work_dble, & 
 N_train, N_test, varphi_train, varphi_test, dir_data)
      
deallocate(Xdata_work_dble)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        Apply the Block Separated Represensation algorithm              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! --- Define parameters for cross-validation procedures
iprint = 1 !1: display messages / -1: no printing
kfold = 5 !number of folds in cross-validation
! Coarse lambda-grid in log10 scale (step 1):
nb_lbd_try_c = 20 
epsi_c = 1e-5
! Fine lambda-grid in log10 scale (step 2):
nb_lbd_try_f = 100 
beta1 = 1e-3
beta2 = 1.0D0
! Tolerances for saturation of objective function:
tol_fobj_cv = 1e-6 !in cross validation
tol_fobj_run = 1e-9 !in production run
! Choose the type of numerical resolution for each subproblem:
opt_resol = 'bfgs' ! 'bfgs' or 'glmn' (LBFGS or GLMNET)

print*, 'nb_lbd_try_c=', nb_lbd_try_c
print*, 'epsi_c=', epsi_c
print*, 'nb_lbd_try_f=', nb_lbd_try_f
print*, 'beta1=', beta1
print*, 'beta2=', beta2
print*, 'tol_fobj_cv=', tol_fobj_cv
print*, 'tol_fobj_run=', tol_fobj_run
print*, 'num resol=', opt_resol

write (99,*) 'nb_lbd_try_c=', nb_lbd_try_c
write (99,*) 'epsi_c=', epsi_c
write (99,*) 'nb_lbd_try_f=', nb_lbd_try_f
write (99,*) 'beta1=', beta1
write (99,*) 'beta2=', beta2
write (99,*) 'tol_fobj_cv=', tol_fobj_cv
write (99,*) 'tol_fobj_run=', tol_fobj_run
write (99,*) 'num resol=', opt_resol

print*, ''
print*, 'paused, type [enter] to continue'
read (*,*) 

allocate(coeffs_mat(rank,d,basissize))
allocate(lbd_grid_c(nb_lbd_try_c), err_test_aver_c(nb_lbd_try_c))
allocate(stock_err_test_c(kfold,nb_lbd_try_c))
allocate(lbd_grid_f(nb_lbd_try_f), err_test_aver_f(nb_lbd_try_f))
allocate(stock_err_test_f(kfold,nb_lbd_try_f))

call BSR_resol(coeffs_mat, rank, d, basissize, N_train, varphi_train, Xdata_train, & 
      Ydata_train, nb_lbd_try_c, epsi_c, nb_lbd_try_f, beta1, beta2, kfold, iprint, opt_resol, & 
       int_resp, lbd_grid_c, err_test_aver_c, stock_err_test_c, lbd_grid_f, err_test_aver_f, &
        stock_err_test_f, flag_cvg, lbd_star, lbd_opt, tstep1, tstep2, tstep3, tol_fobj_cv, tol_fobj_run)

! --- Save errors and lambda-grids in local directory.
!     (to display graphs run the MATLAB code 'plot_errors_BSR.m')
!
! Validation errors for each fold and try of lambda (STEP 1)
file_res = trim(dir_data) // trim('errors_cv_step1_bsr.dat')
open (20, FILE=file_res, STATUS='UNKNOWN', FORM='FORMATTED') 
do ifold = 1,kfold
  do ilbd = 1,nb_lbd_try_c
    write (20, *) stock_err_test_c(ifold,ilbd) 
  end do
end do
close (20)
!
! Averaged validation errors for try of lambda (STEP 1)
file_res = trim(dir_data) // trim('errorsAver_cv_step1_bsr.dat')
open (21, FILE=file_res, STATUS='UNKNOWN', FORM='FORMATTED') 
do ilbd = 1,nb_lbd_try_c
  write (21, *) err_test_aver_c(ilbd)
end do
close (21)
!
! Coarse lambda-grid (STEP 1)
file_res = trim(dir_data) // trim('lbd_grid_step1_bsr.dat')  
open (22, FILE=file_res, STATUS='UNKNOWN', FORM='FORMATTED') 
do ilbd = 1,nb_lbd_try_c
  write (22, *) lbd_grid_c(ilbd)
end do
close (22)
!
! Validation errors for each fold and try of lambda (STEP 2)
file_res = trim(dir_data) // trim('errors_cv_step2_bsr.dat')
open (30, FILE=file_res, STATUS='UNKNOWN', FORM='FORMATTED') 
do ifold = 1,kfold
  do ilbd = 1,nb_lbd_try_f
    write (30, *) stock_err_test_f(ifold,ilbd) 
  end do
end do
close (30)
!
! Averaged validation errors for try of lambda (STEP 2)
file_res = trim(dir_data) // trim('errorsAver_cv_step2_bsr.dat')
open (31, FILE=file_res, STATUS='UNKNOWN', FORM='FORMATTED') 
do ilbd = 1,nb_lbd_try_f
  write (31, *) err_test_aver_f(ilbd)
end do
close (31)
!
! Fine lambda-grid (STEP 2)
file_res = trim(dir_data) // trim('lbd_grid_step2_bsr.dat')  
open (32, FILE=file_res, STATUS='UNKNOWN', FORM='FORMATTED') 
do ilbd = 1,nb_lbd_try_f
  write (32, *) lbd_grid_f(ilbd)
end do
close (32)
        
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            Compute validation error on testing dataset                 !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(ur_test(N_test))
if (flag_cvg .eqv. .true.) then
  call err_normD(Ydata_test, d, N_test, rank, basissize, varphi_test, &
        coeffs_mat, lbd_opt, err_test, fobj_test, normD_ur_test, ur_test, int_resp)
elseif (flag_cvg .eqv. .false.) then
  print*, 'BSR with lbd_opt diverges!!!'
  print*, 'STOP.'
  STOP
end if
  
! Compare a few predicted values against true ones:
! do i = 1,50
!   print*, 'Ytrue=', Ydata_test(i), 'Yhat=', ur_test(i)
! end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute one-dimensional factors (optional: useful for vizualization)   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

opt_fact1D = 'y'

if (opt_fact1D == 'y') then

  ! Compute one-dimensional factors  g_j^l on regular grids:
  ngrid = 50
  allocate(factors1d(rank,d,ngrid), xgrid(ngrid))
  if (opt_basis == 0) then
     call compute_1dfactors_Legendre(rank, d, basissize, coeffs_mat, ngrid, xmin, xmax, &
      factors1d, xgrid)  
  elseif (opt_basis == 1) then
     call compute_1dfactors_rbf(rank, d, basissize, coeffs_mat, ngrid, xmin, xmax, & 
      factors1d, xgrid, dir_data)
  end if  
  
  ! Save data:
  file_res = trim(dir_data) // trim('factors1d.dat')  
  open (20, FILE=file_res, STATUS='UNKNOWN', FORM='FORMATTED') 
  do l = 1,rank
    do j = 1,d
      do m = 1,ngrid
        write (20, *) factors1d(l,j,m)
      end do
    end do  
  end do
  close (20)
   
  file_res = trim(dir_data) // trim('grid1d.dat')   
  open (21, FILE=file_res, STATUS='UNKNOWN', FORM='FORMATTED') 
  do i = 1,ngrid
    write (21, *) xgrid(i)
  end do
  close (21)
  deallocate(factors1d, xgrid)

end if
  
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
write (99,*) '*** Penalization parameters ***'
write (99,*) 'lbd_star=', lbd_star
write (99,*) 'lbd_opt=', lbd_opt
write (99,*) ''
write (99,*) '*** Testing errors ***'
write (99,*) '||u-ur||_D / ||u||_D=', err_test/normD_u_test
write (99,*) 'RMSE=||u-ur||_D=', err_test
write (99,*) '||u||_D=', normD_u_test
write (99,*) '||ur||_D=', normD_ur_test
write (99,*) ''
write (99,*) '*** timings ***'
write (99,*) ' step1=', tstep1
write (99,*) ' step2=', tstep2
write (99,*) ' step3=', tstep3
write (99,*) ' total=', tstep1 + tstep2 + tstep3
write (99,*) ''
write (99,*) '*** sparsity***'
write (99,*) ' tol_spars=', tol_spars
write (99,*) ' percentage of |coeffs/coeff_max| < tol_spars=', perc_spars
write (99,*) ' percentage of non-zeros=', perc_zeros
close (99)
! print*, ''
! print*, 'Check sparsity of model structure'
! do l = 1,rank
!   print*,'* level l=', l
!   do j = 1,d
!     print*,' - dimension=', j
!     do k = 1,basissize
!       print*,'     coeffs=', coeffs_mat(l,j,k)
!     end do
!   end do
! end do
print*, 'max(|coeffs|)=', coeff_max
print*, ''
print*, '*** Penalization parameters ***'
print*, 'lbd_star=', lbd_star
print*, 'lbd_opt=', lbd_opt
print*, ''
print*, '*** Testing errors ***'
print*, '||u-ur||_D / ||u||_D=', err_test/normD_u_test
print*, 'RMSE=||u-ur||_D=', err_test
print*, '||u||_D=', normD_u_test
print*, '||ur||_D=', normD_ur_test
print*, ''
print*, '*** Timings ***'
print*, ' step1=', tstep1
print*, ' step2=', tstep2
print*, ' step3=', tstep3
print*, ' total=', tstep1 + tstep2 + tstep3
print*, ''
print*, '*** Sparsity***'
print*, ' tol_spars=', tol_spars
print "(a,10f5.2)", ' percentage of |coeffs/coeff_max| < tol_spars=', perc_spars
print "(a,10f5.2)", ' percentage of non-zeros=', perc_zeros
print*, ''

! Remark: When SRCD returns zero coefficients the approximate solution is a constant
!         equal to the bias term (equal to the mean of response):
if (coeff_max == 0.0D0) then
  print*, 'Approximate solution = averaged response:'
  print*, 'mean(Y)=', sum(Ydata_train)/N_train !mean over training points
endif 

if (rank == 1 .and. opt_basis == 0 .and. basissize == 1) then 
  if (opt_func=='curved') then 
    print*, 'variance of response=', 52.01
    print*, 'estimated variance=', err_test**2
  end if
end if  
    
! --- deallocations:
deallocate(Xdata_train, Xdata_test, Ydata_train, Ydata_test, varphi_train, varphi_test, coeffs_mat, &
 lbd_grid_c, err_test_aver_c, stock_err_test_c, lbd_grid_f, err_test_aver_f, stock_err_test_f, ur_test)

print*, ''
print*, 'DONE :-)'

END PROGRAM main_driver