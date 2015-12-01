!================================================================================
!> @file main.f90
!! @brief Heat equation problem solver for a multiscale model
!! @author Popov A.M., Nikishin N.G., <mail@nikishinng.ru>
!! @copyright Department of the Automation for Scientific Research CMC MSU
!! @date 2015
!================================================================================
!#define MAGMA
!#define LAPACK
#define TDMA

program heateq
  implicit none
  double precision, parameter :: pi = 3.141592653589793238462643383279502884197
  
#ifdef MAGMA
  call magmaf_init()
#endif

  call test1d()

#ifdef MAGMA
  call magmaf_finalize()
#endif


contains
  !==============================================================================
  !============================= ONE-DIMENSIONAL  ===============================
  !==============================================================================
  subroutine test1d()
    double precision, dimension(:), allocatable :: U  
    double precision, dimension(:), allocatable :: B1
#ifdef TDMA
    double precision, dimension(:), allocatable :: d, du, du2, dl
#else
    double precision, dimension(:,:), allocatable :: A1
#endif

    ! Settings
    double precision :: maxt = 50.0d0
    double precision :: tau = 0.005d0
    integer :: N1 = 20000
    double precision :: Dx, q

!    character(len=32) :: filename
    double precision :: h, h2, coef, coef_a, coef_b, coef_c
    double precision :: err, maxerr
    integer :: maxstep, n, i 

    integer, dimension(:), allocatable :: ipiv
    integer :: lda,  ldb, nrhs, info


    real :: start, finish
    maxstep = int(maxt / tau) ! ps -> tau0
    
    !==============================================================================
    !============================= Allocation =====================================
    !==============================================================================
    allocate(U(0:N1))
    allocate(B1(0:N1))
    allocate(ipiv(N1+1))
    !==============================================================================
    !============================= End of allocation ==============================
    !==============================================================================

    !==============================================================================
    !============================= Initialization =================================
    !==============================================================================
    Dx = 1.0  
    q = 0.0  
    h = 1.0 / N1
    h2 = h*h
    
    ! Init condition
    do i = 0, N1
       U(i) = dsin(pi*i*h)
    end do
    call output1d(U,0,h)
    !==============================================================================
    !============================= End of initialization ==========================
    !==============================================================================

    !============================ Crank-Nicolson scheme ===========================
    lda = N1+1
    ldb = N1+1
    nrhs = 1

    coef = 0.5 * tau  / h2
    coef_a = coef * Dx
    coef_c = coef * Dx
    coef_b = coef_a + coef_c + 0.5 * tau * q

#ifdef TDMA
    allocate(dl(N1))
    allocate(d(N1+1))
    allocate(du(N1))
    allocate(du2(N1-1))

    dl = -coef_a
    d = 1.0 + coef_b
    du = -coef_c
    
    du(1) = 0.0
    d(1) = 1.0
    d(N1+1) = 1.0
    dl(N1+1) = 0.0

    call cpu_time(start)
    call dgttrf(N1+1, dl, d, du, du2, ipiv, info)
    call cpu_time(finish)
    call check_error(info)
    print '("[TDMA] Matrix factorization time = ",f8.4," s.")', finish-start
#else
    allocate(A1(0:N1, 0:N1))
#endif

#ifdef LAPACK
    A1 = 0.0
    A1(0, 0) = 1.0
    do i = 1, N1-1
       A1(i, i-1) =     - coef_a
       A1(i, i  ) = 1.0 + coef_b
       A1(i, i+1) =     - coef_c
    end do
    A1(N1, N1) = 1.0
    call cpu_time(start)
    call dgetrf(N1+1, N1+1, A1, lda, ipiv, info)
    call cpu_time(finish)
    print '("[LU] Matrix factorization time = ",f8.4," s.")', finish-start
    if (info /= 0) then
       print *, "Error: the ", -info, "-th argument had an illegal value"
       return
    end if
#endif 
    !==============================================================================
    !============================= Main evolution loop ============================
    !==============================================================================    
    do n = 1, maxstep
       call cpu_time(start)

#ifdef MAGMA
       A1 = 0.0
       A1(0, 0) = 1.0
       do i = 1, N1-1
          A1(i, i-1) =     - coef_a
          A1(i, i  ) = 1.0 + coef_b
          A1(i, i+1) =     - coef_c
       end do
       A1(N1, N1) = 1.0

       B1 = 0.0
       do i = 1, N1-1
          B1(i) =  coef_a * U(i-1) + coef_c * U(i+1) + &
               & (1.0 - coef_b) * U(i)
       end do

       call magmaf_dgesv(N1+1, nrhs, A1, lda, ipiv, B1, ldb, info)
       call check_error(info)
#endif 

#ifdef LAPACK
       B1 = 0.0
       do i = 1, N1-1
          B1(i) =  coef_a * U(i-1) + coef_c * U(i+1) + &
               & (1.0 - coef_b) * U(i)
       end do
!       call dgesv(N1+1, nrhs, A1, lda, ipiv, B1, ldb, info)
       call dgetrs('N', N1+1, nrhs, A1, lda, ipiv, B1, ldb, info)
       call check_error(info)
#endif 
#ifdef TDMA
       B1 = 0.0
       do i = 1, N1-1
          B1(i) =  coef_a * U(i-1) + coef_c * U(i+1) + &
               & (1.0 - coef_b) * U(i)
       end do

       call dgttrs('N', N1+1, nrhs, dl, d, du, du2, ipiv, B1, ldb, info)
       call check_error(info)
#endif

       U = B1

       call cpu_time(finish)
       
       
       ! Solution 
       maxerr = 0.0
       do i = 1, N1
          err = abs( U(i)  - dsin(pi*i*h) * dexp(- pi*pi*Dx*n*tau) )
          if ( err > maxerr) then 
             maxerr = err
          end if
       end do
       print '(i5,"| Max error = ",f10.8," | Time = ",f5.3," s.")', n, maxerr, finish-start
       if (mod(n,10) == 0) call output1d(U,n,h)
    end do
    !==============================================================================
    !========================== End of main evolution loop ========================
    !==============================================================================
    
    !==============================================================================
    !========================== DeAllocation ======================================
    !==============================================================================
    deallocate(U)
    deallocate(B1)
    deallocate(ipiv)
#ifdef TDMA
    deallocate(dl,d,du,du2)
#else
    deallocate(A1)
#endif
    !==============================================================================
    !========================== End of deallocation ===============================
    !==============================================================================

    !==============================================================================
    !============================= END ONE-DIMENSIONAL  ============================
    !==============================================================================
  end subroutine test1d

  subroutine check_error(info)
    integer, intent(in) :: info

    if (info == 0) return
    if (info < 0) then
       print *, "Error: the ", -info, "-th argument had an illegal value"
    else if (info > 0) then
       print *, "Error: U(", info, ",",info,") is exactly zero." 
    end if
    call exit()

  end subroutine check_error
  

!!$  
!!$  subroutine test3d()
!!$    !==============================================================================
!!$    !============================= THREE-DIMENSIONAL  ============================
!!$    !==============================================================================
!!$    !==================== First direction ===============================     
!!$    do i = 2, N1-1
!!$       A1(i, i-1) =     - coef * Dx
!!$       A1(i, i  ) = 1.0 + 2.0 * coef * Dx + tau * q
!!$       A1(i, i+1) =     - coef * Dx
!!$    end do
!!$    ! Boundaries
!!$    A1(1, 1) = 1.0
!!$    A1(N1, N1) = 1.0
!!$
!!$    do k = 1, N3
!!$       do j = 1, N2     
!!$          do i = 1, N1
!!$             B1(i, j + N2*(k-1)) =  coef * Dx * (U(i-1, j, k, n) + U(i+1, j, k, n)) + &
!!$                  & coef * Dy * (U(i, j-1, k, n) + U(i, j+1, k, n)) + &
!!$                  & coef * Dz * (U(i, j, k-1, n) + U(i, j, k+1, n)) + &
!!$                  & (1.0 - 2.0*coef*(Dx + 2.0*Dy + 2.0*Dz) + 2.5*tau*q) * U(i, j, k, n)
!!$          end do
!!$          ! Boundaries
!!$          B1(1, j + N2*(k-1)) = 0.0
!!$          B1(N1, j + N2*(k-1)) = 0.0
!!$       end do
!!$    end do
!!$    nrhs = N2*N3
!!$    !    call dgesv(N1, nrhs, A1, lda, ipiv, F1, ldb, info)
!!$    call magmaf_dgesv(N1, nrhs, A1, N1, ipiv, B1, N1, info)
!!$    !==================== Second direction ===============================
!!$
!!$    do j = 2, N2-1
!!$       A2(j+1, j) =     - coef * Dy
!!$       A2(j  , j) = 1.0 + 2.0 * coef * Dy + tau * q
!!$       A2(j-1, j) =     - coef * Dy
!!$    end do
!!$    do k = 1, N3
!!$       do i = 1, N1
!!$          do j = 1, N2
!!$             B2(j, i + N1*(k-1)) = B1(i, j, k) - &
!!$                  & coef * Dy * (U(i, j+1, k, n) + U(i, j-1, k, n)) + &
!!$                  & (2.0 * coef * Dy + 0.5 * tau * q) * U(i, j, k, n)
!!$          end do
!!$       end do
!!$    end do
!!$    nrhs = N1*N3
!!$    call magmaf_dgesv(N2, nrhs, A2, N2, ipiv, B2, N2, info)
!!$
!!$    !==================== Third direction ===============================
!!$    do k = 2, N3-1
!!$       A3(k+1, k) =     - coef * Dz
!!$       A3(k  , k) = 1.0 + 2.0 * coef * Dz + tau * q
!!$       A3(k-1, k) =     - coef * Dz
!!$    end do
!!$
!!$    do j = 1, N2
!!$       do i = 1, N1
!!$          do k = 1, N3
!!$             B3(k, i + N1*(j-1)) = B2(i, j, k) - &
!!$                  & coef * Dz * (U(i, j+1, k, n) + U(i, j-1, k, n)) + &
!!$                  & (2.0 * coef * Dz + 0.5 * tau * q) * U(i, j, k, n)
!!$          end do
!!$       end do
!!$    end do
!!$    !==============================================================================
!!$    !============================= THREE-DIMENSIONAL  ============================
!!$    !==============================================================================
!!$  end subroutine test3d
  
!!$  subroutine Set_D(Dx, Dy, Dz)
!!$    double precision, dimension(:), intent(out) :: Dx, Dy, Dz
!!$    integer :: n
!!$
!!$    do n = 1, size(Dx)
!!$       Dx(n) = 1.0
!!$       Dy(n) = 1.0
!!$       Dz(n) = 1.0
!!$    end do
!!$  end subroutine Set_D
!!$
!!$  subroutine Set_Q(q)
!!$    double precision, dimension(:), intent(out) :: q
!!$    integer :: n
!!$
!!$    do n = 1, size(q)
!!$       q(n) = 0.0
!!$    end do
!!$  end subroutine Set_Q
!!$
!!$  subroutine Set_U0(U)
!!$    double precision, dimension(:), intent(out) :: Dx, Dy, Dz
!!$    integer :: n
!!$
!!$    do n = 1, size(Dx)
!!$       Dx(n) = 1.0
!!$       Dy(n) = 1.0
!!$       Dz(n) = 1.0
!!$    end do
!!$  end subroutine Set_D


  subroutine output1d(array, n, h)
    double precision, dimension(:), intent(in) :: array
    double precision, intent(in) :: h
    integer, intent(in) :: n

    character(len=10) :: filename
    integer :: i

    write(filename,'(A1,I0.5,A4)') 'u',n,'.out'
    open(2, file=filename)
    do i = 1, size(array)
       write(2,'(F12.8, F12.8)') real(i*h), real(array(i))
    end do
    close(2)
  end subroutine output1d

!!$  subroutine XDRAW_output1d(filename, array, tau)
!!$    character(len=*), intent(in) :: filename
!!$    double precision, dimension(:), intent(in) :: array
!!$    double precision, intent(in) :: tau
!!$    integer :: step
!!$    
!!$    open(1, file=filename, form='unformatted')
!!$    t = 0.0
!!$    do step = 1, size(array)
!!$       write(1) real(t,4), real(array(step), 4)
!!$       t = t + tau
!!$    end do
!!$    close(1)
!!$  end subroutine XDRAW_output1d
!!$  
!!$  subroutine XDRAW_output2d(filename, array, tau)
!!$    character(len=*), intent(in) :: filename
!!$    double precision, dimension(:,:), intent(in) :: array
!!$    double precision, intent(in) :: tau
!!$    integer :: step
!!$    
!!$    open(2, file=filename, form='unformatted')
!!$    t = 0.0
!!$    do step = 1, size(array)
!!$       write(2) real(t,4), real(array(step,:), 4)
!!$       t = t + tau
!!$    end do
!!$    close(2)
!!$  end subroutine XDRAW_output2d


end program heateq
