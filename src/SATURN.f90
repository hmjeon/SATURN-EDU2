!
! =============================================================================
!
! SATURN Educational Ver. 2 - LU decomposition
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
program SATURN_EDU2

    use Data_Struct
    use DISP4
    use Solver

    implicit none

    call Main

contains

! --------------------------------------------------------------------------------

! Main subroutine
subroutine Main()

    ! Define variables
    type(NodeType), allocatable :: node(:)      ! Node
    type(ElemType), allocatable :: elem(:)      ! Element
    type(PropType) :: prop                      ! Property

    integer :: dof(3)                           ! # of DOFs (total, free, fixed)
    double precision, allocatable :: Kt(:,:)    ! Stiffness vector
    double precision, allocatable :: U(:)       ! Diplacement vector
    double precision, allocatable :: R(:)       ! Load vector
    integer :: n_eq                             ! Number of equations

    ! Open files
    open(unit=1, file="input.txt",      form="formatted")
    open(unit=2, file="SATURN_out.txt", form="formatted")
    open(unit=3, file="SATURN_res.txt", form="formatted")
    open(unit=4, file="SATURN_pos.txt", form="formatted")

    ! Read input file
    call Read_Input_File(node, elem, prop)

    ! # of DOFs and assign equation numbers
    n_eq = Set_DOF_Number(node, dof)

    ! Assemble index
    call Assemble_Index(node, elem, n_eq)

    ! Allocate memory
    allocate(Kt(n_eq, n_eq), U(dof(1)), R(dof(1)))

    ! Print information
    call Print_Information(node, elem, prop, dof, n_eq)

    ! Assemble global stiffness matrix
    write(0, "(a)"), " 1 - Assembling Stiffness and Load"
    call Assemble_Kt(node, elem, prop, Kt, n_eq)

    ! Assemble load vector
    call Assemble_Load(node, elem, prop, R)

    ! Solve linear system
    write(0, "(a)"), " 2 - Solving Linear System"
    call Solver_LU(Kt)
    U(1:n_eq) = matmul(Kt, R(1:n_eq))

    ! Calculate stress and print solutions
    write(0, "(a)"), " 3 - Printing Output Files"
    call Displacement_Stress(node, elem, prop, R, U)

    write(0, "(a)"), " 4 - Completed"
    write(0, "(a)")
    write(0, "(a, es17.10)"), " Strain energy = ", 0.5d0*dot_product(R, U)
    write(0, "(a)"), "   Ref. engrgy =  1.0760861791E-04"
    write(0, "(a)")

    ! Deallocate memory
    deallocate(node, elem, Kt, R, U)

    ! Close files
    close(unit=1); close(unit=2); close(unit=3); close(unit=4); close(unit=5)
end subroutine Main

! --------------------------------------------------------------------------------

! Read input file
subroutine Read_Input_File(node, elem, prop)
    type(NodeType), allocatable, intent(out) :: node(:)
    type(ElemType), allocatable, intent(out) :: elem(:)
    type(PropType), intent(out) :: prop

    character :: bufs
    integer :: bufi, i, node_n, elem_n

    ! Read nodal position vector
    read(1, *) bufs; read(1,*) node_n; read(1,*) bufs
    allocate(node(node_n))

    do i = 1, node_n
        read(1, *) bufi, node(i)%x(1:2), node(i)%bc(1:2), node(i)%pm(1:2)
    end do

    ! Read the element connectivity
    read(1, *) bufs; read(1,*) elem_n; read(1,*) bufs
    allocate(elem(elem_n))

    do i = 1, elem_n
        read(1, *) bufi, elem(i)%cn(:), elem(i)%q(:)
    end do

    ! Read properties
    read(1, *) bufs; read(1, *) prop.thick
    read(1, *) bufs; read(1, *) prop.young
    read(1, *) bufs; read(1, *) prop.poisson
end subroutine Read_Input_File

! --------------------------------------------------------------------------------

! # of DOFs and assign equation numbers to DOFs
! dof(1): # of total DOFs, dof(2): # of free DOFs, dof(3): # of fixed DOFs
function Set_DOF_Number(node, dof) result(n_eq)
    type(NodeType), intent(inout) :: node(:)
    integer, intent(out) :: dof(3)

    integer :: i, j, n_eq

    write(3, *) "EQUATION NUMBER"
    write(3, *) "---------------------"
    write(3, *) "    node   dof    eqn"

    dof(1) = size(node) * nDPN
    dof(2) = 0
    dof(3) = 0

    do i = 1, size(node)
        do j = 1, nDPN
            if (node(i)%bc(j) == 0) then
                dof(2) = dof(2) + 1
                node(i)%eq_n(j) = dof(2)
                write(3, "(3i7)") i, j, dof(2)
            else
                dof(3) = dof(3) + 1
                node(i)%eq_n(j) = dof(1) - dof(3) + 1
            end if
        end do
    end do
    write(3, *)

    n_eq = dof(2)
end function Set_DOF_Number

! --------------------------------------------------------------------------------

! Calculate assemble index
subroutine Assemble_Index(node, elem, n_eq)
    type(NodeType), intent(in) :: node(:)
    type(ElemType), intent(in) :: elem(:)
    integer, intent(in) :: n_eq

    integer, allocatable :: a_index(:)
    integer :: i, j, k, n_elem

     n_elem = size(elem)

    ! Allocate array
    allocate(a_index(nDPN*nNPE))

    do i = 1, n_elem

        ! Assemblage index
        do j = 1, nDPN
            do k = 1, nNPE
                a_index(nNPE*j+k-nNPE) = node(elem(i)%cn(k))%eq_n(j)
            end do
        end do
    end do

    ! Deallocate
    deallocate(a_index)
end subroutine Assemble_Index

! --------------------------------------------------------------------------------

! Assemble total stiffness matrix
subroutine Assemble_Kt(node, elem, prop, Kt, n_eq)
    type(NodeType), intent(in) :: node(:)
    type(ElemType), intent(in) :: elem(:)
    type(PropType), intent(in) :: prop
    double precision, intent(out) :: Kt(:,:)
    integer, intent(in) :: n_eq

    double precision, allocatable :: Ke(:,:)
    integer, allocatable :: a_index(:)
    double precision :: eNode(nNPE, nDIM)
    integer :: i, j, k

    allocate(Ke(nDPN*nNPE, nDPN*nNPE))
    allocate(a_index(nDPN*nNPE))

    do i = 1, n_eq
        do j = 1, n_eq
            Kt(i,j) = 0.0d0
        end do
    end do

    ! Assemble stiffness matrix
    do i = 1, size(elem)

        ! Nodal position of elements
        do j = 1, nNPE
            do k = 1, nDIM
                eNode(j, k) = node(elem(i)%cn(j))%x(k)
            end do
        end do

        ! Planestress element stiffness matrix
        Ke = Plane_Stiffness(prop.young, prop.poisson, prop.thick, eNode)

        ! Print all element stiffness
        call Print_Matrix(Ke)

        ! Assemblage index
        do j = 1, nDPN
            do k = 1, nNPE
                a_index(nNPE*j+k-nNPE) = node(elem(i)%cn(k))%eq_n(j)
            end do
        end do

        ! Assemble stiffness matrix
        do j = 1, nDPN * nNPE
            do k = 1, nDPN * nNPE
                if(a_index(j) <= n_eq .and. a_index(k) <= n_eq) then
                    Kt(a_index(j), a_index(k)) = Kt(a_index(j), a_index(k)) + Ke(j, k)
                end if
            end do
        end do
    end do

    ! Deallocate memory
    deallocate(Ke, a_index)
end subroutine Assemble_Kt

! --------------------------------------------------------------------------------

! assemble load vector
subroutine Assemble_Load(node, elem, prop, R)
    type(NodeType),   intent(in)  :: node(:)
    type(ElemType),   intent(in)  :: elem(:)
    type(PropType),   intent(in)  :: prop
    double precision, intent(out) :: R(:)

    double precision :: eNode(nNPE, nDIM)
    double precision :: NodalLoad(nDPN*nNPE)
    integer :: i, j, k

    R(:) = 0.0d0

    ! Assemble load vector for nodal load
    do i = 1, size(node)
        do j = 1, nDPN
            R(node(i)%eq_n(j)) = node(i)%pm(j)
        end do
    end do

    ! Assemble load vector for body force
    do i = 1, size(elem)

        ! Nodal position of elements
        do j = 1, nNPE
            do k = 1, nDIM
                eNode(j,k) = node(elem(i)%cn(j))%x(k)
            end do
        end do

        ! Calculate equivalent nodal load from body force
        NodalLoad = Plane_Load(eNode, elem(i)%q)

        ! Assemble load vector
        do j = 1, nDPN
            do k = 1, nNPE
                R(node(elem(i)%cn(k))%eq_n(j)) &
                    = R(node(elem(i)%cn(k))%eq_n(j)) + NodalLoad(nNPE*j+k-nNPE)
            end do
        end do
    end do
end subroutine Assemble_Load

! --------------------------------------------------------------------------------

! Calculate stress and print solutions
subroutine Displacement_Stress(node, elem, prop, R, U)
    type(NodeType),   intent(in) :: node(:)
    type(ElemType),   intent(in) :: elem(:)
    type(PropType),   intent(in) :: prop
    double precision, intent(in) :: R(:)
    double precision, intent(in) :: U(:)

    double precision :: eNode(nNPE,nDIM)                ! Nodal position of the element
    double precision :: displace(nNPE*nDPN)             ! Nodal displacement of the element
    double precision :: Stress(nNPE, 3)                 ! Sxx, Syy, Sxy according to # of Gauss points
    double precision :: scale_factor, max_pos, max_disp ! Scaling factor
    integer :: i, j, k

    ! Strain energy
    write(3, "(a17, E14.6)") "STRAIN ENERGY = ", 0.5d0 * dot_product(R, U)
    write(4, *) size(elem)

    ! Set scaling factor for the plotting
    max_pos  = 0.0d0
    max_disp = 0.0d0
    do i=1, size(node)
        if( max_disp < dabs(U(node(i)%eq_n(1))) ) then
            max_disp = dabs(U(node(i)%eq_n(1)))
            max_pos  = sqrt(node(i)%x(1)*node(i)%x(1)+node(i)%x(2)*node(i)%x(2))
        end if

        if( max_disp < dabs(U(node(i)%eq_n(2))) ) then
            max_disp = dabs(U(node(i)%eq_n(2)))
            max_pos  = sqrt(node(i)%x(1)*node(i)%x(1)+node(i)%x(2)*node(i)%x(2))
        end if
    end do

    ! 1.2 * max_pos = (scale_factor * max_disp + max_pos)
    scale_factor = (1.2d0 * max_pos - max_pos) / max_disp

    write(4, "(E14.6)"), scale_factor

    ! Print nodal displacement
    write(3, *)
    write(3, *) "DISPLACEMENT "
    write(3, *) "------------------------------"
    write(3, *) "  Node      Dx         Dy     "
   
    do i = 1, size(node)
        write(3, "(1x,i4,2x,2(1P,E11.3))") i, U(node(i)%eq_n(1)), U(node(i)%eq_n(2))
    end do
    write(3, *)

    do i = 1, size(elem)

        ! Nodal position of element
        do j = 1, nNPE
            do k = 1, nDIM
                eNode(j, k) = node( elem(i)%cn(j) )%x(k)
            end do
        end do
      
        ! Displacement vector of element 
        do j = 1, nDPN
            do k = 1, nNPE
                displace(nNPE*j+k-nNPE) = U(node(elem(i)%cn(k))%eq_n(j))
            end do
        end do

        ! Calculate stress of element
        Stress = Plane_Stress(prop.young, prop.poisson, eNode, displace)
      
        ! Print element stresses
        write(3, "(a21,i4)") " STRESS of ELEMENT : ", i
        write(3, *) "----------------------------------------------" 
        write(3, *) " Position       Sxx        Syy        Sxy     "
      
        do j = 1, nNPE
            write(3, "(1x,i4,5x, 3x,3(1P,E11.3))") j, Stress(j,:)
        end do
        write(3, *)

        ! Print deformed shape and stress for MATLAB post-processing
        write(4, "(1x,28(1P,E13.5))") eNode(1,:), displace(1), displace(5), Stress(1,:),&
                                      eNode(2,:), displace(2), displace(6), Stress(2,:),&
                                      eNode(3,:), displace(3), displace(7), Stress(3,:),&
                                      eNode(4,:), displace(4), displace(8), Stress(4,:)
    end do
end subroutine Displacement_Stress

! --------------------------------------------------------------------------------

! print matrix
subroutine Print_Matrix(M)
    double precision, intent(in) :: M(:,:)
    integer :: i

    write(2, *) "---------------------------"
    do i = 1, nNPE * nDPN
        write(2, "(8E12.4)" ) M(i,:)
    end do

    write(2, *)
end subroutine Print_Matrix

! --------------------------------------------------------------------------------

! Print information
subroutine Print_Information(node, elem, prop, dof, n_eq)
    type(NodeType), intent(in) :: node(:)
    type(ElemType), intent(in) :: elem(:)
    type(PropType), intent(in) :: prop
    integer, intent(in) :: dof(3)
    integer, intent(in) :: n_eq

    write(0, "(a)")
    write(0, "(a)"), " [SATURN - EDU Ver. 2 (LU decomposition)"
    write(0, "(a)")
    write(0, "(a)"), " ====================================================================="
    write(0, "(a, es10.3$)"), " Young's modulus: ", prop%Young
    write(0, "(a,   f6.3$)"), ", Poisson ratio: ", prop%Poisson
    write(0, "(a,   f6.3 )"), ", thick: ", prop%thick
    write(0, "(a, i5$)"), " # of elements: ", size(elem)
    write(0, "(a,  i6)"), ", # of nodes:  ", size(node)
    write(0, "(a, i6$)"), " # total DOFs: ", dof(1)
    write(0, "(a, i6$)"), ", # free DOFs: ", dof(2)
    write(0, "(a, i6 )"), ", # fixed DOFs: ", dof(3)
    write(0, "(a, i6$)")," # equations: ", n_eq
    write(0, "(a)")
    write(0, "(a)"), " ====================================================================="
    write(0, "(a)")
end subroutine Print_Information

! --------------------------------------------------------------------------------

end program SATURN_EDU2