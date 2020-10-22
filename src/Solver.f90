!
! =============================================================================
!
! Solver - LU Decomposition
! Computational Systems Design Laboratory (CSDL), https://csdlab.jbnu.ac.kr
! by Hyungmin Jun (hjun@jbnu.ac.kr)
!
! =============================================================================
!
! Copyright 2020 CSDL. All rights reserved.
!
! License - GPL version 3
! This program is free software: you can redistribute it and/or modify it under
! the terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or any later version.
! This program is distributed in the hope that it will be useful, but WITHOUT
! ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
! FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
! You should have received a copy of the GNU General Public License along with
! this program. If not, see <http://www.gnu.org/licenses/>.
!
! -----------------------------------------------------------------------------
!
module Solver

   public Solver_LU

contains

! --------------------------------------------------------------------------------

! Inverse matrix based on Doolittle LU factorization for Ax=b
subroutine Solver_LU(mat)
    double precision, intent(inout) :: mat(:,:)

    double precision, allocatable :: L(:,:), U(:,:), a(:,:), c(:,:)
    double precision, allocatable :: b(:), d(:), x(:)
    double precision :: coeff
    integer :: i, j, k, n

    n = ubound(mat, 1)

    allocate(L(n,n))
    allocate(U(n,n))
    allocate(a(n,n))
    allocate(c(n,n))
    allocate(b(n))
    allocate(d(n))
    allocate(x(n))

    a = mat

    ! Step 0: Initialization for matrices L and U and b
    L = 0.0d0
    U = 0.0d0
    b = 0.0d0

    ! Step 1: Forward elimination
    do k = 1, n - 1
        do i = k + 1, n
            coeff  = a(i,k) / a(k,k)
            L(i,k) = coeff
            do j = k + 1, n
                a(i,j) = a(i,j) - coeff * a(k,j)
            end do
        end do
    end do

    ! Step 2: Prepare L and U matrices
    ! L matrix is a matrix of the elimination coefficient
    do i = 1, n
        L(i,i) = 1.0d0
    end do

    ! U matrix is the upper triangular part of A
    do j = 1, n
        do i = 1, j
            U(i,j) = a(i,j)
        end do
    end do

    ! Step 3: Compute columns of the inverse matrix C
    do k = 1, n
        b(k) = 1.0d0
        d(1) = b(1)

        ! Step 3a: Solve Ld=b using the forward substitution
        do i = 2, n
            d(i) = b(i)
            do j = 1, i - 1
                d(i) = d(i) - L(i,j) * d(j)
            end do
        end do

        ! Step 3b: Solve Ux=d using the back substitution
        x(n) = d(n) / U(n,n)
        do i = n - 1, 1, -1
            x(i) = d(i)
            do j = n, i + 1, -1
                x(i) = x(i) - U(i,j) * x(j)
            end do
            x(i) = x(i) / u(i,i)
        end do

        ! Step 3c: Fill the solutions x(n) into column k of C
        do i = 1, n
            c(i,k) = x(i)
        end do
        b(k) = 0.0d0
    end do

    mat = c
end subroutine Solver_LU

! --------------------------------------------------------------------------------

end module Solver