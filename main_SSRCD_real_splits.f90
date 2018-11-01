!
! Numerical studies on real datasets coming from Machine Learning (UCI depository from GitLab)
! using original SSR-CD algorithm. Several splits of the data is used for a fair comparison with
! Yang, Smola, Song and Wilson, 'A la Carte - learning fast kernels' (2015).
! This version was used for generating results in the Proc. Royal Society paper.
!

PROGRAM main_driver

USE utilities
USE lpp
USE crossValid
USE coordDescent
USE srcd
USE sigmatune

IMPLICIT NONE

integer                                  :: nobs, nfeat, d, N_train, N_test, rank, opt_basis, &
                                            basissize, orderPol, nb_lbd_try, kfold, iprint, nzcoeffs, &
                                            nsplit, nslice, isplit, cpt, kfold_sigma
integer                                  :: i, j, k, l, ifold, ilbd
integer, dimension(1)                    :: ind_min
character(len=4)                         :: int_resp, opt_resol, opt_sigmatune
character(len=10)                        :: opt_data
character(len=60)                        :: file_res, file_data, dir_data, data_name
real *8, dimension(:,:), allocatable     :: data_orig, Xdata, Xdata_wk, Xdata_train, Xdata_test, stock_err_test, &
                                            testing_norms
real *8, dimension(:), allocatable       :: Ydata, Ydata_wk, Ydata_train, Ydata_test, xpts, PCbasis, ur_test, &
                                            lbd_grid, err_test_aver, rmse, perc_spars
real *8                                  :: xmin, xmax, epsi, err_test, fobj_test, &
                                            lbd_opt, normD_ur_test, normD_u_test, tol_spars, t0, t1, tstep, &
                                            tstep1, tstep2, coeff_max, bias, tol_fobj_cv, tol_fobj_run, &
                                            rmse_mean, rmse_sdv, spars_mean, sigmaval, sigmabest
real *8, dimension(:,:,:), allocatable   :: varphi, varphi_train, varphi_test, varphi_wk, coeffs_mat, testing_vals
logical                                  :: flag_cvg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                       Load input/output variables                       !
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
allocate(Xdata(d,nobs))
allocate(Ydata(nobs))
do i = 1,d
  Xdata(i,:) = data_orig(:,i)
end do
Ydata(:) = data_orig(:,nfeat)
  
! --- Rescale the input variables between 0 and 1:
do i = 1,d
  xmin = minval(Xdata(i,:))
  xmax = maxval(Xdata(i,:))
  if (xmin < xmax) then
    Xdata(i,:) = 1.0D0/(xmax-xmin)*(Xdata(i,:)-xmin)
  elseif (xmin == xmax) then
     !special case when input vars are constant
     Xdata(i,:) = 0.0D0
  endif
end do

! Check the range of input variables (DEBUG):
! print*, 'min(X)=', minval(Xdata)
! print*, 'max(X)=', maxval(Xdata)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!             Define the numbers of training/testing points               !
! Nb of training points: roughly (nsplit-1) partitions of data            !
! Nb of testing points: roughly 1 slit of data                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

kfold = 5 !number of folds in cross-validation (in lambda)
nsplit = 10 !number of splits of data

nslice = int(nobs/nsplit)
! print*, 'nslice=', nslice

! Define the number of TRAINING points:
N_train = (nsplit-1)*nslice !nsplit-1 partitions of data
! N_train must also be a multiple of kfold:
if (mod(N_train,kfold) /= 0) then
  !slighly modify the nb of training points
  N_train = kfold * floor(real(N_train)/real(kfold)) 
end if

! Define the number of TESTING points:
N_test = nobs - N_train

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!               Evaluate the basis functions on the whole dataset         !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

opt_basis = 0  !0 --> Legendre polynomials / 1 --> RBFs
basissize = 2 !Size of discretization basis
rank = 3      !low-rank decomposition rank

if (opt_basis == 1) then
  opt_sigmatune = 'n'
  if (opt_sigmatune == 'y') then
    ! Preliminarty step: use cross-validation to tune RBF shape parameter 
    ! (large sigma means flat RBFs while small sigma means RBFs are peaked)
    kfold_sigma = 2
    opt_resol = 'cd'
    iprint = 1 !1: display messages / -1: no printing
    print*, '*** sigma-tuning for RBF discretization ***'
    print*, opt_data
    print*, 'r=', rank 
    print*, 'M=', basissize
    CALL CPU_TIME(t0)
    CALL crossvalid_sigmatune(Xdata, Ydata, d, nobs, basissize, rank, opt_basis, kfold_sigma, & 
          opt_resol, int_resp, iprint, sigmabest)
    CALL CPU_TIME(t1)
    print*, ' ...sigmabest=', sigmabest
    print*, ' ...Timing sigma-tuning=', t1-t0
    print*, ''
    STOP
  elseif (opt_sigmatune == 'n') then
    ! Reuse value obtained in a previous computation:
    sigmabest = 3.0D0
  endif  
else
  sigmabest = 1.0D0 !dummy variable
endif

! Allocate data structures:
allocate(Xdata_wk(d,nobs), Ydata_wk(nobs), varphi_wk(nobs,basissize,d)) !working data structures
allocate(Xdata_train(d,N_train), Ydata_train(N_train), varphi_train(N_train,basissize,d)) !data structures for training
allocate(Xdata_test(d,N_test), Ydata_test(N_test), varphi_test(N_test,basissize,d)) !data structures for testing

allocate(varphi(nobs,basissize,d))
CALL eval_basisfunc_work(nobs, basissize, d, opt_basis, Xdata, varphi, sigmabest)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         Apply pathwise SRCD algorithm for each split of data            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Define parameters for SRCD:
iprint = -1 !1: display messages / -1: no printing
nb_lbd_try = 20 !size of lambda-grid
epsi = 1e-5 !parameter for defining lower bound in lambda-grid
            ! Remark: if 'epsi' is not small enough SRCD can return zero coefficients,
            !         meaning the approximate model is a constant (equal to the bias term).
! Tolerances for saturation of objective function:
tol_fobj_cv = 1e-6  !in cross validation
tol_fobj_run = 1e-9 !in production run 
opt_resol = 'cd'    !type of numerical resolution
tol_spars = 1e-3    !tolerance for computing level of sparsity

data_name = trim(opt_data) // trim(' dataset')

print*, ''
print*, '========================================================'
print*, data_name
print*, 'dimension                            =', d
print*, 'rank                                 =', rank
print*, 'One-dimensional basis size           =', basissize
select case(opt_basis)
  case(0)
    print*, ' ...discretization: Legendre PC'
  case(1)
    print*, ' ...discretization: RBF with sigma   =', sigmabest
  case default
    print*, 'opt_basis is not correct!'
    print*, 'STOP.'
    STOP
end select
print*, 'Number of points in training dataset =', N_train
print*, 'Number of points in testing dataset  =', N_test
print*, 'Number of splits in data             =', nsplit
print*, '========================================================'
print*, ''

print*, ''
print*, 'paused, type [enter] to continue'
read (*,*) 

allocate(coeffs_mat(rank,d,basissize))
allocate(lbd_grid(nb_lbd_try), err_test_aver(nb_lbd_try))
allocate(stock_err_test(kfold,nb_lbd_try))
allocate(ur_test(N_test))
allocate(rmse(nsplit), perc_spars(nsplit))
allocate(testing_vals(nsplit,N_test,2))
allocate(testing_norms(nsplit,2))

!==============================!
!    Main loop over splits     !
!==============================!

call CPU_TIME(t0)

do isplit = 1,nsplit

  print*, 'isplit=', isplit
  
  ! --- Initialize data structures ---
  Xdata_wk(:,:) = -666; Ydata_wk(:) = -666; varphi_wk(:,:,:) = -666
  Xdata_train(:,:) = -666; Ydata_train(:) = -666; varphi_train(:,:,:) = -666
  Xdata_test(:,:) = -666; Ydata_test(:) = -666; varphi_test(:,:,:) = -666
  
  ! --- Reorganize data in temporary WORKING data structures ---
  cpt = 0
  ! fill in with (nsplit-1) slices:
  do i = 1,nsplit
    if (i /= isplit) then
      cpt = cpt + 1
      Xdata_wk(:,(cpt-1)*nslice+1:cpt*nslice) = Xdata(:,(i-1)*nslice+1:i*nslice)
      Ydata_wk((cpt-1)*nslice+1:cpt*nslice) = Ydata((i-1)*nslice+1:i*nslice)
      varphi_wk((cpt-1)*nslice+1:cpt*nslice,:,:) = varphi((i-1)*nslice+1:i*nslice,:,:)
    endif
  enddo
  ! add the last split at the end:
  Xdata_wk(:,(nsplit-1)*nslice+1:nsplit*nslice) = Xdata(:,(isplit-1)*nslice+1:isplit*nslice)
  Ydata_wk((nsplit-1)*nslice+1:nsplit*nslice) = Ydata((isplit-1)*nslice+1:isplit*nslice)
  varphi_wk((nsplit-1)*nslice+1:nsplit*nslice,:,:) = varphi((isplit-1)*nslice+1:isplit*nslice,:,:)
  if (MOD(nobs,nsplit*nslice) /=0) then 
    ! add the remaining data (if any) at the very end:
    Xdata_wk(:,nsplit*nslice+1:nobs) = Xdata(:,nsplit*nslice+1:nobs) 
    Ydata_wk(nsplit*nslice+1:nobs) = Ydata(nsplit*nslice+1:nobs) 
    varphi_wk(nsplit*nslice+1:nobs,:,:) = varphi(nsplit*nslice+1:nobs,:,:) 
  endif
  
  ! --- Assign data structures for TRAINING ---
  Xdata_train = Xdata_wk(:,1:N_train)
  Ydata_train = Ydata_wk(1:N_train)
  varphi_train = varphi_wk(1:N_train,:,:)
  
  ! --- Assign data structures for TESTING ---
  Xdata_test = Xdata_wk(:,N_train+1:nobs)
  Ydata_test = Ydata_wk(N_train+1:nobs)
  varphi_test = varphi_wk(N_train+1:nobs,:,:)
  
!   ! DEBUG: print testing input variables:
!   do j=1,N_test
!     print*, Xdata_test(:,j)
!   enddo  

  ! --- Apply SRCD algorithm ---
  call srcd_resol(coeffs_mat, rank, d, basissize, N_train, varphi_train, Xdata_train, & 
        Ydata_train, nb_lbd_try, epsi, kfold, iprint, opt_resol, int_resp, flag_cvg, lbd_grid, & 
         err_test_aver, stock_err_test, lbd_opt, tstep1, tstep2, bias, tol_fobj_cv, tol_fobj_run)

  ! --- Compute validation error on testing dataset ---
  if (flag_cvg .eqv. .true.) then     
  call err_normD_with_bias(Ydata_test, d, N_test, rank, basissize, varphi_test, &
        coeffs_mat, lbd_opt, err_test, fobj_test, normD_ur_test, ur_test, bias, int_resp)    
  elseif (flag_cvg .eqv. .false.) then
    print*, 'SRCD with lbd_opt diverges!!!'
    print*, 'STOP.'
    STOP
  endif
  rmse(isplit) = err_test
  ! Store predicted/true values:
  testing_vals(isplit,:,1) = Ydata_test
  testing_vals(isplit,:,2) = ur_test
  ! Store norms of predicted/true solution:
  normD_u_test = sqrt((sum(Ydata_test**2))/N_test)
  testing_norms(isplit,1) = normD_u_test
  testing_norms(isplit,2) = normD_ur_test
  
  ! --- Compute the percentage of almost zero coefficients (scaled by max value) ---
  coeff_max = maxval(abs(coeffs_mat))
  nzcoeffs = count(abs(coeffs_mat)/coeff_max > tol_spars)
  perc_spars(isplit) = 100.0D0*(rank*d*basissize-nzcoeffs)/(rank*d*basissize)
 
  print*, ' RMSE=', err_test
  print*, ' percentage of |coeffs/coeff_max| < tol_spars=', perc_spars(isplit)
  print*,' timing (step1)=', tstep1
  print*,' timing (step1)=', tstep2
  print*,' timing (total)=', tstep1 + tstep2
 
  print*, ''
  
enddo

call CPU_TIME(t1)
tstep = t1-t0

!==============================!
!   end of loop over splits    !
!==============================!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      Recompute best solution to compare true/predicted values          !
!      (not really needed, just for sanity check)                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Find the splitting of data which produces the best RMSE:
ind_min = MINLOC(rmse)
print*, 'best RMSE for split #', ind_min
! Compare predicted values with true ones:
do i = 1,N_test
  print*, 'Ytrue=', testing_vals(ind_min,i,1), 'Ypred=', testing_vals(ind_min,i,2)
end do
print*, '||u||_D=', testing_norms(ind_min,1)
print*, '||ur||_D=', testing_norms(ind_min,2)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                  Diagnostics (errors & sparsity)                      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! --- Compute statistics of errors ---
rmse_mean = SUM(rmse)/nsplit
rmse_sdv = SQRT((SUM((rmse-rmse_mean)**2))/(nsplit-1))

! --- Compute averaged level of sparsity ---
spars_mean = SUM(perc_spars)/nsplit

print*, ''
print*, '************************ Summary ************************'
print*, 'mean(RMSE)=', rmse_mean
print*, 'sdv(RMSE)=', rmse_sdv
print*, 'averaged sparsity (%)=', spars_mean
print*, 'overall timing (s)=', tstep
print*, '*********************************************************'

! --- Deallocate variables ---
deallocate(data_orig, Xdata, Ydata, Xdata_wk, Ydata_wk, Xdata_train, Ydata_train, Xdata_test, Ydata_test, &
 varphi, varphi_wk, varphi_train, varphi_test, coeffs_mat, lbd_grid, err_test_aver, stock_err_test, ur_test, &
  rmse, perc_spars, testing_vals, testing_norms)

print*, ''
print*, 'DONE :-)'

END PROGRAM main_driver