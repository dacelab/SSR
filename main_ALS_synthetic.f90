!
! Standard ALS algorithm (i.e. without sparsity) applied to synthetic datasets.
!

PROGRAM main_driver

USE sobol
USE utilities
USE lpp
USE ALS
USE srcd

IMPLICIT NONE

integer                                  :: N_tot, d, N_train, N_test, rank, opt_basis, &
                                            basissize, orderPol, i, j, k, l, iterALS
character(len=4)                         :: int_resp, opt_datagen
character(len=10)                        :: opt_func
character(len=60)                        :: file_data, dir_data, data_name
real *4, dimension(:,:), allocatable     :: Xdata_work
real *8, dimension(:,:), allocatable     :: Xdata_work_dble, Xdata_train, Xdata_test
real *8, dimension(:), allocatable       :: Ydata_work, Ydata_train, Ydata_test, coeffs, sconsts, & 
                                            ur_test
real *8                                  :: xmin, xmax, err_test, normD_ur_test, t0, t1, tcomp, &
                                            normD_u_test, normD_u_train
real *8, dimension(:,:,:), allocatable   :: varphi_train, varphi_test, coeffs_mat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define input parameters and define training/testing data structures !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

opt_func = 'Friedman3' !choice of function
call setup_dir_synth(opt_func, dir_data, d, int_resp) !define data directory and 
! input dimensionality (see Utilities_mod.f90)

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

rank = 2 !low-rank decomposition rank

! --- Choice of discretization
opt_basis = 0
! 0 --> Legendre PC basis functions 
! 1 --> RBF basis functions (with RBF centers on a regular grid)

! --- Size of discretization basis
basissize = 2
if (opt_basis == 0) then
  orderPol = basissize-1 !PC order
end if

print*, '========================================================'
! select case(opt_func)
!   case default
!     print*, 'opt_func is not correct!'
!     STOP
! end select
data_name = trim(opt_func) // trim(' dataset')
print*, 'Standard ALS algorithm'
print*, data_name
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
print*, 'paused, type [enter] to continue'
read (*,*) 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     Evaluate the basis functions on the training/testing datasets      !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(varphi_train(N_train,basissize,d))
allocate(varphi_test(N_test,basissize,d))

call eval_basisfunc_train_test(N_tot, basissize, d, opt_basis, Xdata_work_dble, & 
 N_train, N_test, varphi_train, varphi_test, dir_data)
      
deallocate(Xdata_work_dble)

normD_u_train = sqrt((sum(Ydata_train**2))/N_train)
normD_u_test = sqrt((sum(Ydata_test**2))/N_test)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                      Apply ALS algorithm                              !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Random initialization of unknown coefficients:
allocate(coeffs(rank*d*basissize))
call random_seed()
call random_number(coeffs)
 coeffs = 2.0D0*coeffs - 1.0D0 ! take values in the interval [-1,1]

! Recast initial coefficients in a 3d data structure:
allocate(coeffs_mat(rank,d,basissize))
do l = 1,rank
  do j = 1,d
    do k = 1,basissize
      coeffs_mat(l,j,k) = coeffs((l-1)*d*basissize+(j-1)*basissize+k) 
    end do
  end do
end do  

allocate(sconsts(rank))
sconsts(:) = 1.0D0 !normalization constants

call CPU_TIME(t0)
call classicalALS(d, N_train, rank, basissize, sconsts, coeffs_mat, varphi_train, & 
      Ydata_train, normD_u_train, iterALS, int_resp)
call CPU_TIME(t1)
tcomp = t1-t0
  
! Compute error on testing dataset:
allocate(ur_test(N_test))
call err_normD_ALS(Ydata_test, sconsts, coeffs_mat, rank, basissize, varphi_test, d, & 
      N_test, err_test, normD_ur_test, ur_test, int_resp)

! Compare a few predicted values against true ones:
! do i = 1,N_test
do i = 1,20
    print*, 'Y=', Ydata_test(i), 'Yhat=', ur_test(i)
end do
 
print*, ''
print*, '||u-ur||_D / ||u||_D=', err_test/normD_u_test
print*, 'RMSE=||u-ur||_D=', err_test 
print*, '||u||_D=', normD_u_test
print*, '||ur||_D=', normD_ur_test
print*, 'Run time=', tcomp
  
!When rank=1 and Legendre polynomials with degree 0 are used, the MSE
!estimates the variance in the realization of the data:
if (rank==1 .and. opt_basis==0 .and. basissize==1) then
    print*, 'Estimation of variance=', (err_test)**2
end if

!deallocations:
deallocate(Xdata_train, Ydata_train, Xdata_test, Ydata_test, varphi_train, varphi_test, & 
 coeffs, sconsts, coeffs_mat, ur_test)

!!! remove all pointless variables in declaration !!!

print*, ''
print*, 'DONE :-)'

END PROGRAM main_driver