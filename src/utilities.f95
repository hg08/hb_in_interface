
  REAL(KIND=8) FUNCTION distance2(r1,r2,boxsize)
      IMPLICIT NONE
      integer, parameter :: rk=8  
      real(kind=rk), DIMENSION(3), INTENT(IN) :: r1,r2
      real(kind=rk), DIMENSION(3), INTENT(IN) :: boxsize
      REAL(kind=rk) :: dx,dy,dz
      dx = r1(1) - r2(1)
      if (abs(dx) > boxsize(1)*0.5d0) then
          dx = boxsize(1) - dx
      endif
      dy = r1(2) - r2(2)
      if (abs(dy) > boxsize(2)*0.5d0) then
          dy = boxsize(2) - dy
      endif
      dz = r1(3) - r2(3)
      if (abs(dz) > boxsize(3)*0.5d0) then
          dz = boxsize(3) - dz
      endif
      distance2 = dx**2 + dy**2 + dz**2
  END FUNCTION distance2
 
  REAL(kind=8) FUNCTION diff_axis (u1, u2, h)
      IMPLICIT NONE
      ! u2 is used as origin
      integer,parameter :: rk=8  
      real(kind=rk),INTENT(IN) :: u1,u2
      real(kind=rk),INTENT(IN) :: h
      real(kind=rk) :: du 
      du = u1-u2
      if (abs(du) > 0.5d0*h) then
          diff_axis = h - du 
      else
          diff_axis= du
      endif
  END FUNCTION diff_axis
  
  REAL(kind=8) FUNCTION pm_adh(r1,r2,r3,boxsize)
      IMPLICIT NONE
      integer, parameter :: rk=8  
      real(kind=rk), DIMENSION(3), INTENT(IN) :: r1,r2,r3
      real(kind=rk), DIMENSION(3), INTENT(IN) :: boxsize 
      REAL(kind=rk) :: diff_axis
      pm_adh=diff_axis(r3(1),r2(1),boxsize(1))*     &
             diff_axis(r1(1),r2(1),boxsize(1))+     &
             diff_axis(r3(2),r2(2),boxsize(2))*     &
             diff_axis(r1(2),r2(2),boxsize(2))+     &
             diff_axis(r3(3),r2(3),boxsize(3))*     &
             diff_axis(r1(3),r2(3),boxsize(3))      ! pm: point multiplication. 
  END FUNCTION pm_adh

  REAL(kind=8) FUNCTION pm_ahd(r1,r2,r3,boxsize)
      IMPLICIT NONE
      integer, parameter :: rk=8  
      REAL(kind=rk), DIMENSION(3), INTENT(IN) :: r1,r2,r3
      REAL(kind=rk), DIMENSION(3), INTENT(IN) :: boxsize 
      REAL(kind=rk) :: diff_axis
      pm_ahd=diff_axis(r1(1),r3(1),boxsize(1))*     &
             diff_axis(r2(1),r3(1),boxsize(1))+     &
             diff_axis(r1(2),r3(2),boxsize(2))*     &
             diff_axis(r2(2),r3(2),boxsize(2))+     &
             diff_axis(r1(3),r3(3),boxsize(3))*     &
             diff_axis(r2(3),r3(3),boxsize(3))      ! pm: point multiplication. 
  END FUNCTION pm_ahd 
  
  INTEGER FUNCTION get_total_number_of_lines(file_)
      IMPLICIT NONE
      INTEGER :: nlines, io
      CHARACTER (len=*),INTENT(IN) :: file_

      !Initialization
      nlines = 0 
      
      OPEN (12, file = file_, status="old")
      DO
          READ(12,*,iostat=io)
          IF (io/=0) EXIT
          nlines = nlines + 1
      END DO
      CLOSE (12)

      get_total_number_of_lines = nlines
  END FUNCTION 

  FUNCTION get_number_of_oxygen(pos_filename, natoms) 
      implicit none
      !Parameters 
      INTEGER, PARAMETER :: rk=8

      !variables
      !=========
      character(LEN=200), INTENT(IN) :: pos_filename ! specific trajectory filename to analyz data
      INTEGER, INTENT(IN) :: natoms
      INTEGER :: get_number_of_oxygen 
      !Local 
      integer :: i,nmovie,iatom,imovie
      real(kind=rk), allocatable, dimension (:,:)           :: x,y,z
      character(LEN=3), allocatable, dimension (:) :: atom_type
      !==================
      !Read data in input
      !==================
      nmovie=1 ! Number of movie steps reqired.
      allocate(atom_type(natoms))
      allocate(x(natoms,nmovie))
      allocate(y(natoms,nmovie))
      allocate(z(natoms,nmovie))
      !=======================
      !read in trajectory file 
      !=======================
      open(10,file=trim(pos_filename))
      REWIND(10)     
      do imovie=1,nmovie  !Note that nmovie=1
          read(10,*)!Neglect data of this line
          read(10,*)                  
          i=0
          do iatom= 1,natoms
              read (10,*)atom_type(iatom),x(iatom,imovie),& 
                       y(iatom,imovie),z(iatom,imovie)
              if (trim(atom_type(iatom)) .eq. 'O') then
                  i=i+1
              endif
          enddo
      enddo
      get_number_of_oxygen = i
      close(10)
      deallocate(atom_type,x,y,z)
  END FUNCTION get_number_of_oxygen

  FUNCTION get_number_of_nitrogen(pos_filename, natoms) 
      implicit none
      !Parameters and variables
      !========================
      integer,parameter :: rk=8  
      character(LEN=200), INTENT(IN) :: pos_filename ! specific trajectory filename to analyz data
      INTEGER, INTENT(IN) :: natoms
      INTEGER :: get_number_of_nitrogen 
      !Local 
      integer :: i, nmovie, iatom, imovie
      real(kind=rk), allocatable, dimension (:,:)           :: x, y, z
      character(LEN=3), allocatable, dimension (:) :: atom_type
      !==================
      !Read data in input
      !==================
      nmovie=1 ! Number of movie steps reqired.
      allocate(atom_type(natoms))
      allocate(x(natoms,nmovie))
      allocate(y(natoms,nmovie))
      allocate(z(natoms,nmovie))
      !=======================
      !read in trajectory file 
      !=======================
      open(10,file=trim(pos_filename))
      REWIND(10)     
      do imovie=1,nmovie  !Note that nmovie=1
          read(10,*)!Neglect data of this line
          read(10,*)                  
          i=0
          do iatom= 1,natoms
              read (10,*)atom_type(iatom),x(iatom,imovie),& 
                       y(iatom,imovie),z(iatom,imovie)
              if (trim(atom_type(iatom)) .eq. 'N') then
                  i=i+1
              endif
          enddo
      enddo
      get_number_of_nitrogen = i
      close(10)
      deallocate(atom_type,x,y,z)
  END FUNCTION get_number_of_nitrogen

  FUNCTION get_number_of_lithium(pos_filename, natoms) 
      implicit none
      !Parameters and variables
      !========================
      integer,parameter :: rk=8  
      character(LEN=200), INTENT(IN) :: pos_filename ! specific trajectory filename to analyz data
      INTEGER, INTENT(IN) :: natoms
      INTEGER :: get_number_of_lithium 
      !Local 
      integer :: i, nmovie, iatom, imovie
      real(kind=rk), allocatable, dimension (:,:)           :: x, y, z
      character(LEN=3), allocatable, dimension (:) :: atom_type
      !==================
      !Read data in input
      !==================
      nmovie=1 ! Number of movie steps reqired.
      allocate(atom_type(natoms))
      allocate(x(natoms,nmovie))
      allocate(y(natoms,nmovie))
      allocate(z(natoms,nmovie))
      !=======================
      !read in trajectory file 
      !=======================
      open(10,file=trim(pos_filename))
      REWIND(10)     
      do imovie=1,nmovie  !Note that nmovie=1
          read(10,*)!Neglect data of this line
          read(10,*)                  
          i=0
          do iatom= 1,natoms
              read (10,*)atom_type(iatom),x(iatom,imovie),& 
                       y(iatom,imovie),z(iatom,imovie)
              if (trim(atom_type(iatom)) .eq. 'Li') then
                  i=i+1
              endif
          enddo
      enddo
      get_number_of_lithium = i
      close(10)
      deallocate(atom_type,x,y,z)
  END FUNCTION get_number_of_lithium

  FUNCTION hydrogen_ndx_list(ndx_oxygen_1, & 
           ndx_oxygen_2,pos_filename,natoms,boxsize) 
      implicit none
      ! 2020/07
      !========================
      !Parameters and variables
      !========================
      integer,parameter :: rk=8  
      INTEGER,INTENT(IN) :: ndx_oxygen_1, ndx_oxygen_2
      character(LEN=200),INTENT(IN) :: pos_filename        ! specific trajectory filename to analyz data
      INTEGER,INTENT(IN) :: natoms
      real(kind=rk),dimension(3),INTENT(IN) :: boxsize
      
      !Local 
      REAL(kind=rk) :: distance2
      INTEGER :: get_number_of_oxygen,get_number_of_hydrogen
      integer :: i,j,iatom,nmovie,& 
                    m1,m2,m3,i_H,&
                    i1,i2,ii,jj,i_OW, num
      INTEGER, DIMENSION(4) :: hydrogen_ndx_list
      real(kind=rk), allocatable,dimension (:,:) :: x,y,z
      character(LEN=3), allocatable, dimension (:) :: atom_type
      integer, allocatable, dimension (:) :: ndx_OW,ndx_H
      real(kind=rk), dimension(3) :: r_a, r_b
      real(kind=rk), parameter :: r_ohc=1.21d0   ! rOH (1.1**2)
      real(kind=rk) :: r ! distance between O and H.
      ! Initialization
      i1=0; i2=0
      !==================
      !Read data in input
      !==================
      nmovie=1 !Number of movie steps reqired.
      allocate(atom_type(natoms))
      allocate(x(natoms,nmovie))
      allocate(y(natoms,nmovie))
      allocate(z(natoms,nmovie))
      !=======================
      ! To obtain the # of H and O atoms
      i_H = get_number_of_hydrogen(pos_filename, natoms) ! By calling functions, we decouple different functions
      i_OW = get_number_of_oxygen(pos_filename, natoms)! By calling functions, we decouple different functions 
      !=======================================================
      ! Calculate the indices of O (H) atoms in water molecules
      !=======================================================
      allocate(ndx_OW(i_OW))    ! this should be put after i_OW is defined
      allocate(ndx_H(i_H))      ! this should be put after i_H is defined
      i=0; ii=0; jj=0
      OPEN(10,file=trim(pos_filename))
          REWIND(10)     
          read(10,*) !Neglect data of this line. Be careful, the two 'read' lines must be located before the 'do iatom=1,natoms' loop, otherwise, there will be an error.
          read(10,*)                  
          do iatom=1,natoms
              read (10,*)atom_type(iatom),x(iatom,nmovie),& 
                  y(iatom,nmovie),z(iatom,nmovie)
              if (trim(atom_type(iatom)) .eq. 'O')then
                     i=i+1      
                     ndx_OW(i)=iatom
              elseif(trim(atom_type(iatom)) .eq. 'H') then
                     ii=ii+1
                     ndx_H(ii)=iatom
              else
              endif
          enddo 
      CLOSE(10)
      deallocate(atom_type)
      !====================
      !Producing the H list      
      !====================
      m1=ndx_oxygen_1
      m2=ndx_oxygen_2
      num = 0 ! For counting
      ! Case one
      do j =1, 1 ! Consider one step
          do ii=1,i_H 
              m3=ndx_H(ii)
              ! I use pbc to consider the distance r
              r_a= (/x(m2,j),y(m2,j),z(m2,j)/) ! Coordinate of O2
              r_b = (/x(m3,j),y(m3,j),z(m3,j)/) ! Coordinate of H
              r=distance2(r_a,r_b,boxsize) ! squared distance between O2 and H
              if (r<r_ohc) then ! If the H is bound to O2
                  num = num + 1
                  hydrogen_ndx_list(num) = m3 ! Add the index of the H atom into the H list for the water pair (O1, O2)
              endif
          enddo
      enddo
      ! Case two
      do j =1, 1
          do ii=1,i_H
              m3=ndx_H(ii)
              ! I use pbc to consider the distance r
              r_a= (/x(m1,j),y(m1,j),z(m1,j)/) ! Coordinate of O1
              r_b = (/x(m3,j),y(m3,j),z(m3,j)/) ! Coordinate of H
              r=distance2(r_a,r_b,boxsize) ! squared distance between O2 and H
              if (r<r_ohc) then ! If the H is bound to O2
                  num = num + 1
                  hydrogen_ndx_list(num) = m3 ! Add the index of the H atom into the H list for the water pair (O1, O2)
              endif
          enddo
      enddo
      deallocate(ndx_OW,ndx_H,x,y,z)
      !=============================
  END FUNCTION hydrogen_ndx_list

  FUNCTION hydrogen_ndx_list_XO(ndx_X, & 
           ndx_oxygen_2,pos_filename, natoms,boxsize) 
      implicit none
      !============
      ! In this function, we use ndx_X replace ndx_oxygen_1. 
      ! ndx_X is not bonded with any Hydrogen atoms
      !========================
      !Parameters and variables
      !========================
      integer,parameter :: rk=8              
      INTEGER, INTENT(IN) :: ndx_X, ndx_oxygen_2
      character(LEN=200), INTENT(IN) :: pos_filename ! specific trajectory filename to analyz data
      INTEGER, INTENT(IN) :: natoms
      real(kind=rk), dimension(3), INTENT(IN) :: boxsize
      
      !Local 
      real(kind=rk), dimension(3) :: r_a, r_b
      real(kind=rk),parameter :: r_ohc=1.21d0  ! rOH (1.1**2)
      real(kind=rk) :: r ! distance between O and H.
      integer :: i,j,nmovie,iatom,& 
              m1,m2,m3,i_H,&
              i1,i2,ii,jj,i_O,num
      INTEGER, DIMENSION(2) :: hydrogen_ndx_list_XO ! the length is 2, because there are 2 Hydrogen are bonded to ndx_oxygen_2 but no Hydrogen is bonded to ndx_X
      real(kind=rk),allocatable,dimension (:,:) :: x,y,z
      character(LEN=3),allocatable,dimension (:) :: atom_type
      integer,allocatable,dimension (:) :: ndx_O,ndx_H
      real(kind=rk) :: distance2
      integer :: get_number_of_hydrogen, get_number_of_oxygen
      ! Initialization
      i1=0; i2=0
      m1=0; m2=0; m3=0
      
      !==================
      !Read data in input
      !==================
      nmovie=1!Number of movie steps reqired.
      allocate(atom_type(natoms))
      allocate(x(natoms,nmovie))
      allocate(y(natoms,nmovie))
      allocate(z(natoms,nmovie))
      !=======================
      ! To obtain the # of H and O atoms
      i_H = get_number_of_hydrogen(pos_filename, natoms) ! By calling functions, we decouple different functions
      i_O = get_number_of_oxygen(pos_filename, natoms)
      !=======================================================
      ! Calculate the indices of O (H) atoms in water molecules
      !=======================================================
      allocate(ndx_O(i_O)) ! this should be put after i_O is defined
      allocate(ndx_H(i_H)) ! this should be put after i_H is defined
      !allocate(ndx_I(i_I)) ! this should be put after i_I is defined
      i=0; ii=0; jj=0
      write(*,*) "pos_filenae:",pos_filename
      open(10,file=trim(pos_filename))
      REWIND(10)     
      read(10,*) !Neglect data of this line. Be careful, the two 'read' lines must be located before the 'do iatom=1,natoms' loop, otherwise, there will be an error.
      read(10,*)                  
      do iatom=1,natoms
          read (10,*)atom_type(iatom),x(iatom,nmovie),& 
                   y(iatom,nmovie),z(iatom,nmovie)
          if (trim(atom_type(iatom)) .eq. 'O') then
                 ndx_O(i+1)=iatom
                 i=i+1      
                 write(*,*) "iatom:(test O) ", iatom
          elseif(trim(atom_type(iatom)) .eq. 'H') then
                 ndx_H(ii+1)=iatom
                 ii=ii+1
                 write(*,*) "iatom:(test) ", iatom
          else
          endif
      enddo 
      close(10)
      deallocate(atom_type)
      !====================
      !Producing the H list      
      !====================
      m1=ndx_X
      m2=ndx_oxygen_2
      num = 0 ! For counting
      ! Case one
      do j =1, 1 ! Consider one step
          do ii=1,i_H 
              m3=ndx_H(ii)
              write(*,*) "i_H: ", i_H
              write(*,*) "m3: ", m3
              ! I use pbc to consider the distance r
              r_a= (/x(m2,j),y(m2,j),z(m2,j)/)  ! Coordinate of O2
              r_b = (/x(m3,j),y(m3,j),z(m3,j)/) ! Coordinate of H
              r=distance2(r_a,r_b,boxsize)  ! squared distance between O2 and H
              if (r<r_ohc) then  ! If the H is bound to O2
                  num = num + 1
                  hydrogen_ndx_list_XO(num) = m3  ! Add the index of the H atom into the H list for the water pair (O1, O2)
              endif
          enddo
      enddo
      ! Case two is canceled
      deallocate(ndx_O,ndx_H,x,y,z)
  END FUNCTION hydrogen_ndx_list_XO

  FUNCTION get_number_of_hydrogen(pos_filename, natoms) 
      implicit none
      !========================
      integer,parameter :: rk=8              
      character(LEN=200),INTENT(IN) :: pos_filename        ! specific trajectory filename to analyz data
      INTEGER,INTENT(IN) :: natoms
      INTEGER :: get_number_of_hydrogen
      !Local 
      integer :: i,nmovie,iatom,imovie
      real(kind=rk),allocatable,dimension (:,:) :: x,y,z
      character(LEN=3),allocatable,dimension (:) :: atom_type
      !==================
      !Read data in input
      !==================
      nmovie=1!Number of movie steps reqired.
      allocate(atom_type(natoms))
      allocate(x(natoms,nmovie))
      allocate(y(natoms,nmovie))
      allocate(z(natoms,nmovie))
      !=======================
      !read in trajectory file 
      !=======================
      open(10,file=trim(pos_filename))
      REWIND(10)     
      do imovie=1,nmovie  !Note that nmovie=1
          read(10,*)!Neglect data of this line
          read(10,*)                  
          i=0
          do iatom= 1,natoms
              read (10,*)atom_type(iatom),x(iatom,imovie),& 
                       y(iatom,imovie),z(iatom,imovie)
              if (trim(atom_type(iatom)) .eq. 'H') then
                  i=i+1
              endif
          enddo
      enddo
      get_number_of_hydrogen = i
      close(10)
      deallocate(atom_type,x,y,z)
  END FUNCTION get_number_of_hydrogen

  FUNCTION get_number_of_iodine(pos_filename, natoms) 
      implicit none
      integer,parameter :: rk=8              
      character(LEN=200),INTENT(IN) :: pos_filename        ! specific trajectory filename to analyz data
      INTEGER,INTENT(IN) :: natoms
      INTEGER :: get_number_of_iodine 
      !Local 
      integer :: i, nmovie, iatom, imovie
      real(KIND=rk),allocatable,dimension (:,:) :: x,y,z
      character(LEN=3),allocatable,dimension (:) :: atom_type
      !==================
      !Read data in input
      !==================
      nmovie=1!Number of movie steps reqired.
      allocate(atom_type(natoms))
      allocate(x(natoms,nmovie))
      allocate(y(natoms,nmovie))
      allocate(z(natoms,nmovie))
      !=======================
      !read in trajectory file 
      !=======================
      open(10,file=trim(pos_filename))
      REWIND(10)     
      do imovie=1,nmovie  !Note that nmovie=1
          read(10,*)!Neglect data of this line
          read(10,*)                  
          i=0
          do iatom= 1,natoms
              read (10,*)atom_type(iatom),x(iatom,imovie),& 
                       y(iatom,imovie),z(iatom,imovie)
              if (trim(atom_type(iatom)) .eq. 'I') then
                  i=i+1
              endif
          enddo
      enddo
      get_number_of_iodine = i
      close(10)
      deallocate(atom_type,x,y,z)
  END FUNCTION get_number_of_iodine

  INTEGER FUNCTION grid_index(x,y,divx,divy,n_divx)
      ! transfer the coordinates (x,y) to grid_index, which is an integer
      implicit none
      integer,parameter :: rk=8              
      REAL(kind=rk), INTENT(IN) :: x,y
      REAL(kind=rk), INTENT(IN) :: divx, divy
      INTEGER, INTENT(IN) :: n_divx
      INTEGER, DIMENSION(2) :: ind
      
      !Initialization
      ind = 0

      ind(1) = FLOOR(x/divx) 
      ind(2) = FLOOR(y/divy)
      
      grid_index = ind(2) * n_divx + ind(1)+1

  END FUNCTION grid_index

  !TODO: define get_number_of_oxygens_in_nitrate()
  !TODO: define get_number_of_oxygens_in_water()
  !TODO: define get_number_of_iodine() Done


  LOGICAL FUNCTION pair_in_surf1(surf1_mol1,z1,surf1_mol2,z2,thickness)
      implicit none
      integer, parameter :: rk=8
      LOGICAL :: mol1_in_surf1, mol2_in_surf1
      REAL(kind=rk), INTENT(IN) :: surf1_mol1,surf1_mol2
      REAL(kind=rk), INTENT(IN) :: thickness
      REAL(kind=rk), INTENT(IN) :: z1,z2
      mol1_in_surf1 = surf1_mol1 + & 
          thickness > z1
      mol2_in_surf1 = surf1_mol2 + &
          thickness > z2
      pair_in_surf1 = mol1_in_surf1 .and. mol2_in_surf1
  END FUNCTION pair_in_surf1    

  LOGICAL FUNCTION pair_in_surf2(surf2_mol1,z1,surf2_mol2,z2,thickness)
      implicit none
      integer, parameter :: rk=8
      LOGICAL :: mol1_in_surf2, mol2_in_surf2
      REAL(kind=rk), INTENT(IN) :: surf2_mol1,surf2_mol2
      REAL(kind=rk), INTENT(IN) :: thickness
      REAL(kind=rk), INTENT(IN) :: z1,z2
      mol1_in_surf2 = surf2_mol1 - & 
          thickness < z1
      mol2_in_surf2 = surf2_mol2 - &
          thickness < z2
      pair_in_surf2 = mol1_in_surf2 .and. mol2_in_surf2
  END FUNCTION pair_in_surf2    

  LOGICAL FUNCTION mol_in_surf1(surf1_mol,z1,thickness)
      !Check if a molecule is in the surf1 ( the lower layer)
      implicit none
      integer, parameter :: rk=8
      REAL(kind=rk), INTENT(IN) :: surf1_mol
      REAL(kind=rk), INTENT(IN) :: thickness
      REAL(kind=rk), INTENT(IN) :: z1
      mol_in_surf1 = surf1_mol + thickness > z1
  END FUNCTION mol_in_surf1    

  LOGICAL FUNCTION mol_in_surf2(surf2_mol,z2,thickness)
      !Check if a molecule is in the surf2 ( the upper layer)
      implicit none
      integer, parameter :: rk=8
      REAL(kind=rk), INTENT(IN) :: surf2_mol
      REAL(kind=rk), INTENT(IN) :: thickness
      REAL(kind=rk), INTENT(IN) :: z2
      mol_in_surf2 = surf2_mol - thickness < z2
  END FUNCTION mol_in_surf2    
  
  character(len=20) function str(k)
      ! Convert an integer to string.
      implicit none
      integer, intent(in) :: k
      write (str, *) k
      str = adjustl(str)
  end function str
  
  FUNCTION nth(k,n)
      ! To return a string containing the first N characters
      ! of the alphabet.
      IMPLICIT NONE
  
      ! Declaring calling parameters:
      CHARACTER(len=20), INTENT(IN) :: k ! String which is adjustl-ed
      INTEGER, INTENT(IN) :: n ! Length of string to return
      CHARACTER(len=n) nth ! Returned string
      
      ! Get string to return
      nth = k(1:n)
  END FUNCTION nth

  INTEGER FUNCTION sampling_number(nmo_start,nmo_end,ns)
      ! To calculate the total numbers of samples one WANT to include in their analysis.
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nmo_start, nmo_end  ! To get the total number of moves
      INTEGER, INTENT(IN) :: ns ! Get one sample from the trajectory every ns step.

      write(*,*) 'In function sampling_number: nmo_end = ', nmo_end
      ! no. of samples = INT({no. of moves}/ns)
      positive: IF (nmo_end <0 .OR. nmo_start < 0 .OR. ns <0) THEN
        write(*,*) 'Please enter non-negative values for the ns, starting step and ending step.'
      ELSE IF (nmo_end < nmo_start) THEN
        write(*,*) 'Please note that starting step shoud not larger than  ending step.'
      ELSE IF (ns ==0) THEN
        sampling_number = nmo_end-(nmo_start-1)
      ELSE 
        sampling_number = FLOOR(FLOAT(nmo_end-(nmo_start))/FLOAT(ns))
      END IF positive
  END FUNCTION sampling_number

  SUBROUTINE read_traj(indx,nmo_start,nmo_end,ns,nat,n_samples,sampled_movie,sampled_time,sampled_energy,atom_info)
      ! To read info from the trajectory file (format: ***.xyz)
      ! to READ data starting from a pattern-matched line.
      USE module_ihb, ONLY: atom
      IMPLICIT NONE

      integer, parameter :: rk=8

      CHARACTER(LEN=4) :: head_char
      INTEGER :: iatom, i_sample
      INTEGER,INTENT(IN) :: indx
      INTEGER,INTENT(IN) :: nat
      INTEGER,INTENT(IN) :: nmo_start, nmo_end  ! To get the total number of moves
      INTEGER,INTENT(IN) :: ns  ! Get one sample from the trajectory every ns step.
      INTEGER,INTENT(IN) :: n_samples  !n_samples = INT(nmo/ns)

      TYPE(atom),DIMENSION(nat,n_samples),INTENT(INOUT) :: atom_info
      INTEGER,DIMENSION(n_samples) :: sampled_movie
      REAL(kind=rk),DIMENSION(n_samples) :: sampled_time, sampled_energy
      INTEGER :: y
      
      i_sample = 1
      write(*,*) "In utilities.f95, read_traj(): New total time steps (n_samples):", n_samples
      DO WHILE (i_sample < n_samples+1) ! +1 means i_sample can take the value of n_samples 
          read(indx, '(1X,A4)') head_char  ! for some other format, one can use this format
          PRE_CHECK:IF (head_char=="i = ") THEN
              BACKSPACE(UNIT=indx) ! Because I am not able to read other lines with the format '(A4,I8)', and have not find any good way, so I try to read it in '(A4)' first 
              read(indx, '(1X,A4,I8)') head_char, y ! for some other format, one can use this format
              CHECK_HEAD:IF (head_char=="i = " .AND. (y>nmo_start-1 .and. y<nmo_end+1) .AND. MOD(y-nmo_start,ns) == 0) THEN
                  !-------------------------------------------------------------------------------------------------------
                  !We use y>nmo_start-1, because we want to consider the first step 'i=0'
                  !-------------------------------------------------------------------------------------------------------
                  WRITE(*,*)"read_traj():", head_char, y
                  BACKSPACE(UNIT=indx) ! Because we have to read the whole line with ' i = ' line.
                  read(indx,130) sampled_movie(i_sample), sampled_time(i_sample), sampled_energy(i_sample)
                  130 FORMAT (1X,4X,I8,9X,F12.3,6X,F20.10)
                  131 FORMAT (A4,3F20.10)
                  inner: do iatom= 1,nat
                    read (indx,131) atom_info(iatom, i_sample)%atom_name, atom_info(iatom,i_sample)%coord(1), & 
                      atom_info(iatom,i_sample)%coord(2), atom_info(iatom,i_sample)%coord(3)
                    if (atom_info(iatom, i_sample)%atom_name == "O") THEN
                        atom_info(iatom, i_sample)%mass = 16.00d0
                    elseif (atom_info(iatom, i_sample)%atom_name == "H") THEN
                        atom_info(iatom, i_sample)%mass = 1.00d0
                    elseif (atom_info(iatom, i_sample)%atom_name == "N") THEN
                        atom_info(iatom, i_sample)%mass = 14.00d0 
                    elseif (atom_info(iatom, i_sample)%atom_name == "Li") THEN
                        atom_info(iatom, i_sample)%mass = 6.94d0
                    elseif (atom_info(iatom, i_sample)%atom_name == "Na") THEN
                        atom_info(iatom, i_sample)%mass = 22.99d0
                    elseif (atom_info(iatom, i_sample)%atom_name == "K") THEN
                        atom_info(iatom, i_sample)%mass = 39.10d0
                    endif
                    !WRITE (*,131) & 
                    !atom_info(iatom, i_sample)%atom_name, atom_info(iatom,i_sample)%coord(1), &
                    !atom_info(iatom,i_sample)%coord(2), atom_info(iatom,i_sample)%coord(3)
                  enddo inner
                  i_sample = i_sample + 1 !The position is important. It must be located before ENDIF 
              ENDIF CHECK_HEAD
          ENDIF PRE_CHECK
      END DO
  END SUBROUTINE read_traj

  SUBROUTINE skip_lines(indx, i_input)
      ! To skip lines when read data from the input
      IMPLICIT NONE
      INTEGER :: i
      INTEGER,INTENT(IN) :: i_input,indx
      do i=1,i_input
         read(indx,*) !Neglect (nat+2)*(ns-1) lines
      enddo    
  END SUBROUTINE skip_lines
