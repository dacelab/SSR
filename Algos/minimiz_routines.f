!
! These routines are called by LBFGS main driver.
!

!
! Evaluation of objective function in the l1-penalized minimization problem
!
      subroutine evalFcn(N, X, F, N_train, rank, d, basissize, lbda, 
     +                   varphi_train, Ydata_train)
      
      implicit none
      integer          N, N_train, rank, d, basissize
      double precision lbda, X(N), F, varphi_train(N_train,basissize,d),
     +                 Ydata_train(N_train) 
      
      ! Local variables
      integer          i, j, k, l 
      double precision tmp, ur, ur_l, urdata(N_train), 
     +                 X3d(rank,d,basissize), errnormD2, normxl1
            
      ! Recast X using a 3d data structure
      do l = 1,rank
        do j = 1,d
          do k = 1,basissize
            X3d(l,j,k) = X((l-1)*d*basissize+(j-1)*basissize+k) 
          end do
        end do
      end do  
      
      ! Compute approximate function at each training data point
      do i = 1,N_train
        ur = 0.0d0
        do l = 1,rank !rank loop
          ur_l = 1.0d0
          do j = 1,d ! loop on dimension
            call dot_prod(X3d(l,j,1:basissize), 
     +           varphi_train(i,1:basissize,j), basissize, tmp)
            ur_l = ur_l * tmp
          end do
          ur = ur + ur_l !add contribution of rank level l
        end do 
        urdata(i) = ur !store approximate value
      end do  
      
      ! Compute the square of discrete semi-norm error
      errnormD2 = 0.0d0
      do i =1,N_train
         errnormD2 = errnormD2 + (Ydata_train(i)-urdata(i))**2
      end do
      errnormD2 = errnormD2/N_train
      
      ! Compute l1 norm of X
      normxl1 = 0.0d0
      do i = 1,N
        normxl1 = normxl1 + ABS(X(i))
      end do
      
      ! Final expression for objective function
      F = errnormD2 + lbda * normxl1
      
      RETURN
      END
      
!
! Evaluation of the gradient of objective function in the l1-penalized 
! minimization problem (no bias term here)
!
      subroutine evalGradFcn(N, X, gradF, N_train, rank, d, 
     + basissize, lbda, varphi_train, Ydata_train)  
      
      implicit none
      integer          N, N_train, rank, d, basissize
      double precision lbda, X(N), gradF(N), 
     +                 varphi_train(N_train,basissize,d),
     +                 Ydata_train(N_train)
      
      ! Local variables
      integer          i, j, k, l, jj, indg
      double precision X3d(rank,d,basissize), tmp, ur, ur_l, 
     +                 urdata(N_train), tmp_part, sgn,
     +                 partial_urdata(N_train,rank,d,basissize)
      
      ! Recast X using a 3d data structure
      do l = 1,rank
        do j = 1,d
          do k = 1,basissize
            X3d(l,j,k) = X((l-1)*d*basissize+(j-1)*basissize+k) 
          end do
        end do
      end do
      
      ! Compute approximate function at each training data point
      do i = 1,N_train
        ur = 0.0d0
        do l = 1,rank !rank loop
          ur_l = 1.0d0
          do j = 1,d ! loop on dimension
            call dot_prod(X3d(l,j,1:basissize), 
     +           varphi_train(i,1:basissize,j), basissize, tmp)
            ur_l = ur_l * tmp
          end do
          ur = ur + ur_l !add contribution of rank level l
        end do 
        urdata(i) = ur !store approximate value
      end do  
      
      ! Compute partial approximations at each training point
      do i = 1,N_train
       do l = 1,rank
        do j = 1,d
         do k = 1,basissize
           tmp_part = 1.0d0
           do jj = 1,d
            if (jj /= j) then
              call dot_prod(X3d(l,jj,1:basissize),
     +         varphi_train(i,1:basissize,jj), basissize, tmp)
              tmp_part = tmp_part * tmp
            end if
           end do !jj
           partial_urdata(i,l,j,k) = tmp_part * 
     +                               varphi_train(i,k,j)
         end do !k
        end do !j
       end do !l
      end do !i
      
      ! Compute the gradient of objective function
      do l = 1,rank
       do j = 1,d
        do k = 1,basissize
          tmp = 0.0d0
          do i = 1,N_train
            tmp = tmp + (Ydata_train(i)-urdata(i)) * 
     +       partial_urdata(i,l,j,k)  
          end do !i
          indg = (l-1)*d*basissize + (j-1)*basissize + k
          call absprime(X3d(l,j,k), sgn)
          gradF(indg) = -2.0d0*tmp/N_train + lbda * sgn
        end do !k
       end do !j
      end do !l
      
      RETURN
      END
      
!
! Computes the dot product of arrays a and b of n elements 
!
      subroutine dot_prod(a, b, n, z)
      
      implicit none
      integer          n
      double precision a(n), b(n), z
      ! Local variables
      integer          i
      
      z = 0.0d0
      do i = 1,n
          z = z + a(i) * b(i) 
      end do
      
      RETURN
      END 
      
!
! Computes the derivative of absolute value of a real, that is,
! the sign of this real.
!      
      subroutine absprime(x, sgn)
      
      implicit none
      double precision x, sgn
      
      if (x > 0.0d0) then
        sgn = 1.0d0
      else if (x < 0.0d0) then
        sgn = -1.0d0
      else if (x == 0.0d0) then
        sgn = 0.0d0
      end if
      
      RETURN
      END
      
!
! Computes the product y = A * x where A is a matrix 
! and x is a column vector.
!
      
      subroutine prod_mat_vec(A, nra, nca, x, y)
      
      implicit none
      integer          nra, nca
      double precision A(nra,nca), x(nca), y(nra)
      ! Local variables
      integer          i, j
      double precision sum
      
      do i = 1,nra
        sum = 0.0d0
        do j = 1,nca
          sum = sum + A(i,j) * x(j)
        end do  
        y(i) = sum
      end do
      
      RETURN
      END
      
!
! Computes the product of matrices C = A * B.
!     

      subroutine matmul77(A, nra, nca, B, nrb, ncb, C)
      
      integer          nra, nca, nrb, ncb
      double precision A(nra,nca), B(nrb,ncb), C(nra,ncb)
      ! Local variables
      integer          i, j, k
      double precision tmp
      
      if (nca .NE. nrb) then
        print*, 'Sizes of matrices are incompatible!!!'
        print*, 'STOP.'
        STOP
      end if
      
      C(:,:) = 0.0d0
      do i = 1,nra
        do j = 1,ncb
          tmp = 0.0d0
          do k = 1,nca
            tmp = tmp + A(i,k) * B(k,j) 
          end do
          C(i,j) = tmp
        end do
      end do  
      
      RETURN
      END
      
!
! Evaluation of objective function in the l1-penalized BLOCK minimization problem
! Note: N must be equal to rank * basissize.
!

      subroutine evalFcn_block(N, X, F, N_train, rank, d, basissize, 
     + lbda, varphi_train, comp1d_train, nbblock, Ydata_train)
      
      implicit none
      integer          N, N_train, rank, d, basissize, nbblock
      double precision lbda, X(N), F, varphi_train(N_train,basissize,d),
     +                 Ydata_train(N_train),
     +                 comp1d_train(rank,d,N_train)
     
      ! Local variables
      double precision matA(N_train,N), Yfit(N_train)
      double precision partproduct, normxl1, err2sq
      integer          i, k, l, jj
      
      ! Assemble block-structured matrix
      matA(:,:) = 0.0d0
      do l = 1,rank
        ! Compute partial products of one-dimensional functions 
        ! (Cf. coefficients a_{jl}^{(i)} in SRCD paper)
        do i = 1,N_train
          partproduct = 1.0d0
          do jj = 1,d
            if (jj /= nbblock) then
              partproduct = partproduct * comp1d_train(l,jj,i)
            end if
          end do !jj
          ! Assemble i-th row of submatrix
          matA(i,(l-1)*basissize+1:l*basissize) = partproduct 
     +                                      * varphi_train(i,:,nbblock)
        end do !i
      end do !l
      
      ! Compute squared l2 error norm 
      call prod_mat_vec(matA, N_train, N, X, Yfit)
      err2sq = 0.0d0
      do i =1,N_train
         err2sq = err2sq + (Ydata_train(i)-Yfit(i))**2
      end do
      err2sq = 0.5*err2sq/N_train
       
      ! Compute l1 norm of X
      normxl1 = 0.0d0
      do i = 1,N
        normxl1 = normxl1 + ABS(X(i))
      end do
      
      ! Final expression for objective function
      F = err2sq + lbda * normxl1
           
      RETURN
      END
      
!
! Evaluation of the gradient of objective function in the l1-penalized 
! BLOCK minimization problem.
!      

      subroutine evalGradFcn_block(N, X, gradF, N_train, rank, d, 
     + basissize, lbda, varphi_train, comp1d_train, nbblock, 
     + Ydata_train)
      
      implicit none
      integer          N, N_train, rank, d, basissize, nbblock
      double precision lbda, X(N), gradF(N),
     +                 varphi_train(N_train,basissize,d),
     +                 Ydata_train(N_train), 
     +                 comp1d_train(rank,d,N_train)
     
      ! Local variables
      double precision matA(N_train,N), matAT(N,N_train), 
     +                 Yfit(N_train), res(N_train)
      double precision partproduct, normxl1, err2sq, sgn
      integer          i, j, k, l, jj
      
      ! Assemble block-structured matrix
      matA(:,:) = 0.0d0
      do l = 1,rank
        ! Compute partial products of one-dimensional functions 
        ! (Cf. coefficients a_{jl}^{(i)} in SRCD paper)
        do i = 1,N_train
          partproduct = 1.0d0
          do jj = 1,d
            if (jj /= nbblock) then
              partproduct = partproduct * comp1d_train(l,jj,i)
            end if
          end do !jj
          ! Assemble i-th row of submatrix
          matA(i,(l-1)*basissize+1:l*basissize) = partproduct 
     +                                      * varphi_train(i,:,nbblock)
        end do !i
      end do !l
      
      ! Transpose matA
      matAT(:,:) = 0.0d0
      do i = 1,N
        do j = 1,N_train
          matAT(i,j) = matA(j,i)
        end do
      end do  
      
      ! Compute residual
      call prod_mat_vec(matA, N_train, N, X, Yfit)
      res = Ydata_train - Yfit
      
      ! Compute gradient of quadratic loss term
      call prod_mat_vec(matAT, N, N_train, res, gradF)
      gradF = -gradF/N_train
      
      ! Add the gradient of l1 penalization term
      do i = 1,N
        call absprime(X(i), sgn)
        gradF(i) = gradF(i) + lbda * sgn
      end do
      
      RETURN
      END