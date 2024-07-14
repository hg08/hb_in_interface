  !REAL(KIND=8) FUNCTION distance2(r1,r2,boxsize)
  !    IMPLICIT NONE
  !    !NOTICE: This function is totally WRONG, I used this function before Mar. 24, 2021.!!
  !    INTEGER, parameter :: rk=8  
  !    REAL(KIND=rk), DIMENSION(3), INTENT(IN) :: r1,r2
  !    REAL(KIND=rk), DIMENSION(3), INTENT(IN) :: boxsize
  !    REAL(KIND=rk) :: dx,dy,dz
  !    dx = r1(1) - r2(1)
  !    if (abs(dx) > boxsize(1)*0.5d0) THEN
  !        dx = boxsize(1) - dx
  !    endif
  !    dy = r1(2) - r2(2)
  !    if (abs(dy) > boxsize(2)*0.5d0) THEN
  !        dy = boxsize(2) - dy
  !    endif
  !    dz = r1(3) - r2(3)
  !    if (abs(dz) > boxsize(3)*0.5d0) THEN
  !        dz = boxsize(3) - dz
  !    endif
  !    distance2 = dx**2 + dy**2 + dz**2
  !END FUNCTION distance2

  REAL(KIND=8) FUNCTION distance2(r1,r2,boxsize)
      ! Date: 2021-3-24
      ! This is the correct version for calculate distance between two particles with PBC considered.
      USE module_ihb, ONLY: rk
      IMPLICIT NONE

      REAL(KIND=rk), DIMENSION(3), INTENT(IN) :: r1,r2
      REAL(KIND=rk), DIMENSION(3), INTENT(IN) :: boxsize
      REAL(KIND=rk) :: dx,dy,dz
      REAL(KIND=rk) :: temp

      dx =0; dy=0; dz=0; temp=0
 
      dx = r1(1) - r2(1)
      temp = boxsize(1)*0.5d0
      if (dx > temp) THEN
          dx = boxsize(1) - dx
      else if (dx < - temp) THEN
          dx = boxsize(1) + dx
      endif
      dy = r1(2) - r2(2)
      temp = boxsize(2)*0.5d0
      if (dy > temp) THEN
          dy = boxsize(2) - dy
      else if (dy < - temp) THEN
          dy = boxsize(2) + dy
      endif
      dz = r1(3) - r2(3)
      temp = boxsize(3)*0.5d0
      if (dz > temp) THEN
          dz = boxsize(3) - dz
      else if (dz < - temp) THEN
          dz = boxsize(3) + dz
      endif
      distance2 = dx**2 + dy**2 + dz**2
  END FUNCTION distance2
 
  REAL(KIND=8) FUNCTION P2(x)
      IMPLICIT NONE
      ! P2: 2-order Legendre polynomial
      INTEGER,parameter :: rk=8  
      REAL(KIND=rk),INTENT(IN) :: x
      P2 = 0.5d0 * (3.d0 * x**2 -1.d0)
  END FUNCTION P2 

  REAL(KIND=8) FUNCTION diff_axis (u1, u2, h)
      IMPLICIT NONE
      ! Considering PBC
      ! u2 is used as origin
      INTEGER,parameter :: rk=8  
      REAL(KIND=rk),INTENT(IN) :: u1,u2
      REAL(KIND=rk),INTENT(IN) :: h
      REAL(KIND=rk) :: du 
      du = u1-u2
      if (abs(du) > 0.5d0*h) THEN
          diff_axis = h - du 
      else
          diff_axis= du
      endif
  END FUNCTION diff_axis
  
  REAL(KIND=8) FUNCTION pm_adh(r1,r2,r3,boxsize)
      IMPLICIT NONE
      INTEGER, parameter :: rk=8  
      REAL(KIND=rk), DIMENSION(3), INTENT(IN) :: r1,r2,r3
      REAL(KIND=rk), DIMENSION(3), INTENT(IN) :: boxsize 
      REAL(KIND=rk) :: diff_axis
      pm_adh=diff_axis(r3(1),r2(1),boxsize(1))*     &
             diff_axis(r1(1),r2(1),boxsize(1))+     &
             diff_axis(r3(2),r2(2),boxsize(2))*     &
             diff_axis(r1(2),r2(2),boxsize(2))+     &
             diff_axis(r3(3),r2(3),boxsize(3))*     &
             diff_axis(r1(3),r2(3),boxsize(3))      ! pm: point multiplication. 
  END FUNCTION pm_adh

  REAL(KIND=8) FUNCTION pm_ahd(r1,r2,r3,boxsize)
      IMPLICIT NONE
      INTEGER, parameter :: rk=8  
      REAL(KIND=rk), DIMENSION(3), INTENT(IN) :: r1,r2,r3
      REAL(KIND=rk), DIMENSION(3), INTENT(IN) :: boxsize 
      REAL(KIND=rk) :: diff_axis
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

  FUNCTION get_ndx_of_oxygen(pos_filename, natoms) 
      implicit none
      !Parameters 
      INTEGER, PARAMETER :: rk=8

      !variables
      !=========
      character(LEN=200), INTENT(IN) :: pos_filename ! specific trajectory filename to analyz data
      INTEGER, INTENT(IN) :: natoms
      INTEGER, DIMENSION(natoms) :: get_ndx_of_oxygen 
      !Local 
      INTEGER :: i,nmovie,iatom,imovie
      REAL(KIND=rk), ALLOCATABLE, DIMENSION (:,:)  :: x,y,z
      character(LEN=3), ALLOCATABLE, DIMENSION (:) :: atom_type
      
      !Initialization
      get_ndx_of_oxygen = 0
      !==================
      !Read data in input
      !==================
      nmovie=1 ! Number of movie steps reqired.
      ALLOCATE(atom_type(natoms))
      ALLOCATE(x(natoms,nmovie))
      ALLOCATE(y(natoms,nmovie))
      ALLOCATE(z(natoms,nmovie))
      !=======================
      !read in trajectory file 
      !=======================
      open(10,file=trim(pos_filename))
      REWIND(10)     
      do imovie=1,nmovie  !Note that nmovie=1
          READ(10,*)!Neglect data of this line
          READ(10,*)                  
          i=0
          do iatom= 1,natoms
              read (10,*)atom_type(iatom),x(iatom,imovie),& 
                       y(iatom,imovie),z(iatom,imovie)
              if (trim(atom_type(iatom)) .eq. 'O') THEN
                  i=i+1
                  get_ndx_of_oxygen(i) = iatom ! ONCE, I write imovie here. This mistake leads to obvious errors in the results
              endif
          enddo
      enddo
      close(10)
      deALLOCATE(atom_type,x,y,z)
  END FUNCTION get_ndx_of_oxygen

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
      INTEGER :: i,nmovie,iatom,imovie
      REAL(KIND=rk), ALLOCATABLE, DIMENSION (:,:)           :: x,y,z
      character(LEN=3), ALLOCATABLE, DIMENSION (:) :: atom_type
      !==================
      !Read data in input
      !==================
      nmovie=1 ! Number of movie steps reqired.
      ALLOCATE(atom_type(natoms))
      ALLOCATE(x(natoms,nmovie))
      ALLOCATE(y(natoms,nmovie))
      ALLOCATE(z(natoms,nmovie))
      !=======================
      !read in trajectory file 
      !=======================
      open(10,file=trim(pos_filename))
      REWIND(10)     
      do imovie=1,nmovie  ! Note that nmovie=1
          READ(10,*) ! Neglect data of this line
          READ(10,*)                  
          i=0
          do iatom= 1, natoms
              read (10,*) atom_type(iatom),x(iatom,imovie),& 
                       y(iatom,imovie),z(iatom,imovie)
              if (trim(atom_type(iatom)) .eq. 'O') THEN
                  i=i+1
              endif
          enddo
      enddo
      get_number_of_oxygen = i
      close(10)
      deALLOCATE(atom_type,x,y,z)
  END FUNCTION get_number_of_oxygen

  FUNCTION get_number_of_nitrogen(pos_filename, natoms) 
      implicit none
      !Parameters and variables
      !========================
      INTEGER,parameter :: rk=8  
      character(LEN=200), INTENT(IN) :: pos_filename ! specific trajectory filename to analyz data
      INTEGER, INTENT(IN) :: natoms
      INTEGER :: get_number_of_nitrogen 
      !Local 
      INTEGER :: i, nmovie, iatom, imovie
      REAL(KIND=rk), ALLOCATABLE, DIMENSION (:,:)           :: x, y, z
      character(LEN=3), ALLOCATABLE, DIMENSION (:) :: atom_type
      !==================
      !Read data in input
      !==================
      nmovie=1 ! Number of movie steps reqired.
      ALLOCATE(atom_type(natoms))
      ALLOCATE(x(natoms,nmovie))
      ALLOCATE(y(natoms,nmovie))
      ALLOCATE(z(natoms,nmovie))
      !=======================
      !read in trajectory file 
      !=======================
      open(10,file=trim(pos_filename))
      REWIND(10)     
      do imovie=1,nmovie  !Note that nmovie=1
          READ(10,*)!Neglect data of this line
          READ(10,*)                  
          i=0
          do iatom= 1,natoms
              read (10,*)atom_type(iatom),x(iatom,imovie),& 
                       y(iatom,imovie),z(iatom,imovie)
              if (trim(atom_type(iatom)) .eq. 'N') THEN
                  i=i+1
              endif
          enddo
      enddo
      get_number_of_nitrogen = i
      close(10)
      deALLOCATE(atom_type,x,y,z)
  END FUNCTION get_number_of_nitrogen

  FUNCTION get_number_of_lithium(pos_filename, natoms) 
      implicit none
      !Parameters and variables
      !========================
      INTEGER,parameter :: rk=8  
      character(LEN=200), INTENT(IN) :: pos_filename ! specific trajectory filename to analyz data
      INTEGER, INTENT(IN) :: natoms
      INTEGER :: get_number_of_lithium 
      !Local 
      INTEGER :: i, nmovie, iatom, imovie
      REAL(KIND=rk), ALLOCATABLE, DIMENSION (:,:)           :: x, y, z
      character(LEN=3), ALLOCATABLE, DIMENSION (:) :: atom_type
      !==================
      !Read data in input
      !==================
      nmovie=1 ! Number of movie steps reqired.
      ALLOCATE(atom_type(natoms))
      ALLOCATE(x(natoms,nmovie))
      ALLOCATE(y(natoms,nmovie))
      ALLOCATE(z(natoms,nmovie))
      !=======================
      !read in trajectory file 
      !=======================
      open(10,file=trim(pos_filename))
      REWIND(10)     
      do imovie=1,nmovie  !Note that nmovie=1
          READ(10,*)!Neglect data of this line
          READ(10,*)                  
          i=0
          do iatom= 1,natoms
              read (10,*)atom_type(iatom),x(iatom,imovie),& 
                       y(iatom,imovie),z(iatom,imovie)
              if (trim(atom_type(iatom)) .eq. 'Li') THEN
                  i=i+1
              endif
          enddo
      enddo
      get_number_of_lithium = i
      close(10)
      deALLOCATE(atom_type,x,y,z)
  END FUNCTION get_number_of_lithium

  FUNCTION hydrogen_ndx_list(ndx_oxygen_1, & 
           ndx_oxygen_2,pos_filename,natoms,boxsize) 
      implicit none
      ! 2020/07
      !========================
      !Parameters and variables
      !========================
      INTEGER,parameter :: rk=8  
      INTEGER,INTENT(IN) :: ndx_oxygen_1, ndx_oxygen_2
      character(LEN=200),INTENT(IN) :: pos_filename        ! specific trajectory filename to analyz data
      INTEGER,INTENT(IN) :: natoms
      REAL(KIND=rk),DIMENSION(3),INTENT(IN) :: boxsize
      
      !Local 
      REAL(KIND=rk) :: distance2
      INTEGER :: get_number_of_oxygen,get_number_of_hydrogen
      INTEGER :: i,j,iatom,nmovie,& 
                    m1,m2,m3,i_H,&
                    i1,i2,ii,jj,i_OW, num
      INTEGER, DIMENSION(4) :: hydrogen_ndx_list
      REAL(KIND=rk), ALLOCATABLE,DIMENSION (:,:) :: x,y,z
      character(LEN=3), ALLOCATABLE, DIMENSION (:) :: atom_type
      INTEGER, ALLOCATABLE, DIMENSION (:) :: ndx_OW,ndx_H
      REAL(KIND=rk), DIMENSION(3) :: r_a, r_b
      REAL(KIND=rk), parameter :: r_ohc=1.21d0   ! rOH (1.1**2)
      REAL(KIND=rk) :: r ! distance between O and H.
      ! Initialization
      i1=0; i2=0
      !==================
      !Read data in input
      !==================
      nmovie=1 !Number of movie steps reqired.
      ALLOCATE(atom_type(natoms))
      ALLOCATE(x(natoms,nmovie))
      ALLOCATE(y(natoms,nmovie))
      ALLOCATE(z(natoms,nmovie))
      !=======================
      ! To obtain the # of H and O atoms
      i_H = get_number_of_hydrogen(pos_filename, natoms) ! By calling functions, we decouple different functions
      i_OW = get_number_of_oxygen(pos_filename, natoms)! By calling functions, we decouple different functions 
      !=======================================================
      ! Calculate the indices of O (H) atoms in water molecules
      !=======================================================
      ALLOCATE(ndx_OW(i_OW))    ! this should be put after i_OW is defined
      ALLOCATE(ndx_H(i_H))      ! this should be put after i_H is defined
      i=0; ii=0; jj=0
      OPEN(10,file=trim(pos_filename))
          REWIND(10)     
          READ(10,*) !Neglect data of this line. Be careful, the two 'read' lines must be located before the 'do iatom=1,natoms' loop, otherwise, there will be an error.
          READ(10,*)                  
          do iatom=1,natoms
              read (10,*)atom_type(iatom),x(iatom,nmovie),& 
                  y(iatom,nmovie),z(iatom,nmovie)
              if (trim(atom_type(iatom)) .eq. 'O')THEN
                     i=i+1      
                     ndx_OW(i)=iatom
              ELSEIF(trim(atom_type(iatom)) .eq. 'H') THEN
                     ii=ii+1
                     ndx_H(ii)=iatom
              else
              endif
          enddo 
      CLOSE(10)
      deALLOCATE(atom_type)
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
              if (r<r_ohc) THEN ! If the H is bound to O2
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
              if (r<r_ohc) THEN ! If the H is bound to O2
                  num = num + 1
                  hydrogen_ndx_list(num) = m3 ! Add the index of the H atom into the H list for the water pair (O1, O2)
              endif
          enddo
      enddo
      deALLOCATE(ndx_OW,ndx_H,x,y,z)
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
      INTEGER,parameter :: rk=8              
      INTEGER, INTENT(IN) :: ndx_X, ndx_oxygen_2
      character(LEN=200), INTENT(IN) :: pos_filename ! specific trajectory filename to analyz data
      INTEGER, INTENT(IN) :: natoms
      REAL(KIND=rk), DIMENSION(3), INTENT(IN) :: boxsize
      
      !Local 
      REAL(KIND=rk), DIMENSION(3) :: r_a, r_b
      REAL(KIND=rk),parameter :: r_ohc=1.21d0  ! rOH (1.1**2)
      REAL(KIND=rk) :: r ! distance between O and H.
      INTEGER :: i,j,nmovie,iatom,& 
              m1,m2,m3,i_H,&
              i1,i2,ii,jj,i_O,num
      INTEGER, DIMENSION(2) :: hydrogen_ndx_list_XO ! the length is 2, because there are 2 Hydrogen are bonded to ndx_oxygen_2 but no Hydrogen is bonded to ndx_X
      REAL(KIND=rk),ALLOCATABLE,DIMENSION (:,:) :: x,y,z
      character(LEN=3),ALLOCATABLE,DIMENSION (:) :: atom_type
      INTEGER,ALLOCATABLE,DIMENSION (:) :: ndx_O,ndx_H
      REAL(KIND=rk) :: distance2
      INTEGER :: get_number_of_hydrogen, get_number_of_oxygen
      ! Initialization
      i1=0; i2=0
      m1=0; m2=0; m3=0
      
      !==================
      !Read data in input
      !==================
      nmovie=1!Number of movie steps reqired.
      ALLOCATE(atom_type(natoms))
      ALLOCATE(x(natoms,nmovie))
      ALLOCATE(y(natoms,nmovie))
      ALLOCATE(z(natoms,nmovie))
      !=======================
      ! To obtain the # of H and O atoms
      i_H = get_number_of_hydrogen(pos_filename, natoms) ! By calling functions, we decouple different functions
      i_O = get_number_of_oxygen(pos_filename, natoms)
      !=======================================================
      ! Calculate the indices of O (H) atoms in water molecules
      !=======================================================
      ALLOCATE(ndx_O(i_O)) ! this should be put after i_O is defined
      ALLOCATE(ndx_H(i_H)) ! this should be put after i_H is defined
      !ALLOCATE(ndx_I(i_I)) ! this should be put after i_I is defined
      i=0; ii=0; jj=0
      open(10,file=trim(pos_filename))
      REWIND(10)     
      READ(10,*) !Neglect data of this line. Be careful, the two 'read' lines must be located before the 'do iatom=1,natoms' loop, otherwise, there will be an error.
      READ(10,*)                  
      do iatom=1,natoms
          read (10,*)atom_type(iatom),x(iatom,nmovie),& 
                   y(iatom,nmovie),z(iatom,nmovie)
          if (trim(atom_type(iatom)) .eq. 'O') THEN
                 ndx_O(i+1)=iatom
                 i=i+1      
                 !write(*,*) "iatom:(test O) ", iatom
          ELSEIF(trim(atom_type(iatom)) .eq. 'H') THEN
                 ndx_H(ii+1)=iatom
                 ii=ii+1
                 !write(*,*) "iatom:(test) ", iatom
          else
          endif
      enddo 
      close(10)
      deALLOCATE(atom_type)
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
              !write(*,*) "i_H: ", i_H
              !write(*,*) "m3: ", m3
              ! I use pbc to consider the distance r
              r_a= (/x(m2,j),y(m2,j),z(m2,j)/)  ! Coordinate of O2
              r_b = (/x(m3,j),y(m3,j),z(m3,j)/) ! Coordinate of H
              r=distance2(r_a,r_b,boxsize)  ! squared distance between O2 and H
              if (r<r_ohc) THEN  ! If the H is bound to O2
                  num = num + 1
                  hydrogen_ndx_list_XO(num) = m3  ! Add the index of the H atom into the H list for the water pair (O1, O2)
              endif
          enddo
      enddo
      ! Case two is canceled
      deALLOCATE(ndx_O,ndx_H,x,y,z)
  END FUNCTION hydrogen_ndx_list_XO

  FUNCTION get_number_of_hydrogen(pos_filename, natoms) 
      implicit none
      !========================
      INTEGER,parameter :: rk=8              
      character(LEN=200),INTENT(IN) :: pos_filename        ! specific trajectory filename to analyz data
      INTEGER,INTENT(IN) :: natoms
      INTEGER :: get_number_of_hydrogen
      !Local 
      INTEGER :: i,nmovie,iatom,imovie
      REAL(KIND=rk),ALLOCATABLE,DIMENSION (:,:) :: x,y,z
      character(LEN=3),ALLOCATABLE,DIMENSION (:) :: atom_type
      !==================
      !Read data in input
      !==================
      nmovie=1!Number of movie steps reqired.
      ALLOCATE(atom_type(natoms))
      ALLOCATE(x(natoms,nmovie))
      ALLOCATE(y(natoms,nmovie))
      ALLOCATE(z(natoms,nmovie))
      !=======================
      !read in trajectory file 
      !=======================
      open(10,file=trim(pos_filename))
      REWIND(10)     
      do imovie=1,nmovie  !Note that nmovie=1
          READ(10,*)!Neglect data of this line
          READ(10,*)                  
          i=0
          do iatom= 1,natoms
              read (10,*)atom_type(iatom),x(iatom,imovie),& 
                       y(iatom,imovie),z(iatom,imovie)
              if (trim(atom_type(iatom)) .eq. 'H') THEN
                  i=i+1
              endif
          enddo
      enddo
      get_number_of_hydrogen = i
      close(10)
      deALLOCATE(atom_type,x,y,z)
  END FUNCTION get_number_of_hydrogen

  FUNCTION get_number_of_iodine(pos_filename, natoms) 
      implicit none
      INTEGER,parameter :: rk=8              
      character(LEN=200),INTENT(IN) :: pos_filename        ! specific trajectory filename to analyz data
      INTEGER,INTENT(IN) :: natoms
      INTEGER :: get_number_of_iodine 
      !Local 
      INTEGER :: i, nmovie, iatom, imovie
      REAL(KIND=rk),ALLOCATABLE,DIMENSION (:,:) :: x,y,z
      character(LEN=3),ALLOCATABLE,DIMENSION (:) :: atom_type
      !==================
      !Read data in input
      !==================
      nmovie=1!Number of movie steps reqired.
      ALLOCATE(atom_type(natoms))
      ALLOCATE(x(natoms,nmovie))
      ALLOCATE(y(natoms,nmovie))
      ALLOCATE(z(natoms,nmovie))
      !=======================
      !read in trajectory file 
      !=======================
      open(10,file=trim(pos_filename))
      REWIND(10)     
      do imovie=1,nmovie  !Note that nmovie=1
          READ(10,*)!Neglect data of this line
          READ(10,*)                  
          i=0
          do iatom= 1,natoms
              read (10,*)atom_type(iatom),x(iatom,imovie),& 
                       y(iatom,imovie),z(iatom,imovie)
              if (trim(atom_type(iatom)) .eq. 'I') THEN
                  i=i+1
              endif
          enddo
      enddo
      get_number_of_iodine = i
      close(10)
      deALLOCATE(atom_type,x,y,z)
  END FUNCTION get_number_of_iodine

  INTEGER FUNCTION grid_index(x,y,divx,divy,nb_divx,nb_divy)
      ! transfer the coordinates (x,y) to grid_index, which is an INTEGER
      implicit none
      INTEGER,parameter :: rk=8              
      REAL(KIND=rk), INTENT(IN) :: x,y
      REAL(KIND=rk), INTENT(IN) :: divx, divy
      INTEGER, INTENT(IN) :: nb_divx, nb_divy
      INTEGER, DIMENSION(2) :: ind
      
      !Initialization
      ind = 0

      ind(1) = FLOOR(x/divx) 
      ind(2) = FLOOR(y/divy)

      !Correction for avoiding 0 index
      IF (ind(1) == 0) THEN
          ind(1) = ind(1) + 1
      ENDIF
      IF (ind(2) == 0) THEN
          ind(2) = ind(2) + 1
      ENDIF

      grid_index = (ind(2)-1) * nb_divx + ind(1) 
      !grid_index = ind(1) * n_divy + ind(2)
      !WRITE (*,*) "NB_DIVX = ", nb_divx
      !WRITE (*,*) "nb_divy = ", nb_divy
      !WRITE (*,*) "index(1) = ", ind(1)
      !WRITE (*,*) "index(2) = ", ind(2)
      !WRITE (*,*) "grid_index = ", grid_index
  END FUNCTION grid_index

  !TODO: define get_number_of_oxygens_in_nitrate()
  !TODO: define get_number_of_oxygens_in_water()
  !TODO: define get_number_of_iodine() Done

  LOGICAL FUNCTION oh_in_surf1(surf1_mol1,z1,thickness)
      implicit none
      INTEGER, parameter :: rk=8
      LOGICAL :: mol1_in_surf1
      REAL(KIND=rk), INTENT(IN) :: surf1_mol1
      REAL(KIND=rk), INTENT(IN) :: thickness
      REAL(KIND=rk), INTENT(IN) :: z1
      mol1_in_surf1 = surf1_mol1 + & 
          thickness > z1
      oh_in_surf1 = mol1_in_surf1 
  END FUNCTION oh_in_surf1    

  LOGICAL FUNCTION oh_in_surf2(surf2_mol1,z1,thickness)
      implicit none
      INTEGER, parameter :: rk=8
      LOGICAL :: mol1_in_surf2
      REAL(KIND=rk), INTENT(IN) :: surf2_mol1
      REAL(KIND=rk), INTENT(IN) :: thickness
      REAL(KIND=rk), INTENT(IN) :: z1
      mol1_in_surf2 = surf2_mol1 - & 
          thickness < z1
      oh_in_surf2 = mol1_in_surf2 
  END FUNCTION oh_in_surf2    

  LOGICAL FUNCTION pair_in_surf1(surf1_mol1,z1,surf1_mol2,z2,thickness)
      implicit none
      INTEGER, parameter :: rk=8
      LOGICAL :: mol1_in_surf1, mol2_in_surf1
      REAL(KIND=rk), INTENT(IN) :: surf1_mol1,surf1_mol2
      REAL(KIND=rk), INTENT(IN) :: thickness
      REAL(KIND=rk), INTENT(IN) :: z1,z2
      mol1_in_surf1 = surf1_mol1 + & 
          thickness > z1
      mol2_in_surf1 = surf1_mol2 + &
          thickness > z2
      pair_in_surf1 = mol1_in_surf1 .and. mol2_in_surf1
  END FUNCTION pair_in_surf1    

  LOGICAL FUNCTION pair_in_surf2(surf2_mol1,z1,surf2_mol2,z2,thickness)
      implicit none
      INTEGER, parameter :: rk=8
      LOGICAL :: mol1_in_surf2, mol2_in_surf2
      REAL(KIND=rk), INTENT(IN) :: surf2_mol1,surf2_mol2
      REAL(KIND=rk), INTENT(IN) :: thickness
      REAL(KIND=rk), INTENT(IN) :: z1,z2
      mol1_in_surf2 = surf2_mol1 - & 
          thickness < z1
      mol2_in_surf2 = surf2_mol2 - &
          thickness < z2
      pair_in_surf2 = mol1_in_surf2 .and. mol2_in_surf2
  END FUNCTION pair_in_surf2    

  LOGICAL FUNCTION mol_in_surf1(surf1_mol,z1,thickness)
      !Check if a molecule is in the surf1 ( the lower layer)
      implicit none
      INTEGER, parameter :: rk=8
      REAL(KIND=rk), INTENT(IN) :: surf1_mol
      REAL(KIND=rk), INTENT(IN) :: thickness
      REAL(KIND=rk), INTENT(IN) :: z1
      mol_in_surf1 = surf1_mol + thickness > z1
  END FUNCTION mol_in_surf1    

  LOGICAL FUNCTION mol_in_surf2(surf2_mol,z2,thickness)
      !Check if a molecule is in the surf2 ( the upper layer)
      implicit none
      INTEGER, parameter :: rk=8
      REAL(KIND=rk), INTENT(IN) :: surf2_mol
      REAL(KIND=rk), INTENT(IN) :: thickness
      REAL(KIND=rk), INTENT(IN) :: z2
      mol_in_surf2 = surf2_mol - thickness < z2
  END FUNCTION mol_in_surf2    
  
  character(len=20) function int2str(k)
      !  "Convert an INTEGER to string."
      INTEGER, intent(in) :: k
      write (int2str, *) k
      int2str = adjustl(int2str)
  end function int2str

  character(len=20) function str(k)
      !  "Convert an INTEGER/REAL to string."
      REAL (KIND=8), intent(in) :: k
      write (str, '(F3.1)') k
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
      ELSE IF (ns == 0) THEN
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

      INTEGER, parameter :: rk=8

      CHARACTER(LEN=4) :: head_char
      INTEGER :: iatom, i_sample
      INTEGER,INTENT(IN) :: indx
      INTEGER,INTENT(IN) :: nat
      INTEGER,INTENT(IN) :: nmo_start, nmo_end  ! To get the total number of moves
      INTEGER,INTENT(IN) :: ns  ! Get one sample from the trajectory every ns step.
      INTEGER,INTENT(IN) :: n_samples  !n_samples = INT(nmo/ns)

      TYPE(atom),DIMENSION(nat,n_samples),INTENT(INOUT) :: atom_info
      INTEGER,DIMENSION(n_samples) :: sampled_movie
      REAL(KIND=rk),DIMENSION(n_samples) :: sampled_time, sampled_energy
      INTEGER :: y
      
      i_sample = 1
      write(*,*) "In utilities.f95, read_traj(): New total time steps (n_samples):", n_samples
      DO WHILE (i_sample < n_samples+1) ! +1 means i_sample can take the value of n_samples 
          READ(indx, '(1X,A4)') head_char  ! for some other format, one can use this format
          PRE_CHECK:IF (head_char=="i = ") THEN
              BACKSPACE(UNIT=indx) ! Because I am not able to read other lines with the format '(A4,I8)', and have not find any good way, so I try to read it in '(A4)' first 
              READ(indx, '(1X,A4,I8)') head_char, y ! for some other format, one can use this format
              CHECK_HEAD:IF (head_char=="i = " .AND. (y>nmo_start-1 .and. y<nmo_end+1) .AND. MOD(y-nmo_start,ns) == 0) THEN
                  !-------------------------------------------------------------------------------------------------------
                  !We use y>nmo_start-1, because we want to consider the first step 'i=0'
                  !-------------------------------------------------------------------------------------------------------
                  !WRITE(*,*)"read_traj():", head_char, y
                  BACKSPACE(UNIT=indx) ! Because we have to read the whole line with ' i = ' line.
                  READ(indx,130) sampled_movie(i_sample), sampled_time(i_sample), sampled_energy(i_sample)
                  130 FORMAT (1X,4X,I8,9X,F12.3,6X,F20.10)
                  131 FORMAT (A4,3F20.10)
                  inner: do iatom= 1,nat
                    read (indx,131) atom_info(iatom, i_sample)%atom_name, atom_info(iatom,i_sample)%coord(1), & 
                      atom_info(iatom,i_sample)%coord(2), atom_info(iatom,i_sample)%coord(3)
                    if (atom_info(iatom, i_sample)%atom_name == "O") THEN
                        atom_info(iatom, i_sample)%mass = 16.00d0
                    ELSEIF (atom_info(iatom, i_sample)%atom_name == "H") THEN
                        atom_info(iatom, i_sample)%mass = 1.00d0
                    ELSEIF (atom_info(iatom, i_sample)%atom_name == "N") THEN
                        atom_info(iatom, i_sample)%mass = 14.00d0 
                    ELSEIF (atom_info(iatom, i_sample)%atom_name == "Li") THEN
                        atom_info(iatom, i_sample)%mass = 6.94d0
                    ELSEIF (atom_info(iatom, i_sample)%atom_name == "Na") THEN
                        atom_info(iatom, i_sample)%mass = 22.99d0
                    ELSEIF (atom_info(iatom, i_sample)%atom_name == "K") THEN
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

  SUBROUTINE read_traj_sphere(indx,nmo_start,nmo_end,ns,nat,n_samples,sampled_movie,sampled_time,sampled_energy,sphere_info)
      ! To read info from the trajectory file (format: ***.xyz)
      ! to READ data starting from a pattern-matched line.
      USE module_ihb, ONLY: sphere
      IMPLICIT NONE

      INTEGER, parameter :: rk = 8

      CHARACTER(LEN=4) :: head_char
      INTEGER :: iatom, i_sample
      INTEGER,INTENT(IN) :: indx
      INTEGER,INTENT(IN) :: nat
      INTEGER,INTENT(IN) :: nmo_start, nmo_end  ! To get the total number of moves
      INTEGER,INTENT(IN) :: ns  ! Get one sample from the trajectory every ns step.
      INTEGER,INTENT(IN) :: n_samples  !n_samples = INT(nmo/ns)

      TYPE(sphere),DIMENSION(nat,n_samples),INTENT(INOUT) :: sphere_info
      INTEGER,DIMENSION(n_samples) :: sampled_movie
      REAL(KIND=rk),DIMENSION(n_samples) :: sampled_time, sampled_energy
      INTEGER :: y
      
      i_sample = 1
      write(*,*) "In utilities.f95, read_traj_sphere(): New total time steps (n_samples):", n_samples
      DO WHILE (i_sample < n_samples+1) ! +1 means i_sample can take the value of n_samples 
          READ(indx, '(1X,A4)') head_char  ! for some other format, one can use this format
          PRE_CHECK:IF (head_char=="i = ") THEN
              BACKSPACE(UNIT=indx) ! Because I am not able to read other lines with the format '(A4,I8)', and have not find any good way, so I try to read it in '(A4)' first 
              READ(indx, '(1X,A4,I8)') head_char, y ! for some other format, one can use this format
              CHECK_HEAD:IF (head_char=="i = " .AND. (y>nmo_start-1 .and. y<nmo_end+1) .AND. MOD(y-nmo_start,ns) == 0) THEN
                  !-------------------------------------------------------------------------------------------------------
                  !We use y>nmo_start-1, because we want to consider the first step 'i=0'
                  !-------------------------------------------------------------------------------------------------------
                  !WRITE(*,*)"read_traj():", head_char, y
                  BACKSPACE(UNIT=indx) ! Because we have to read the whole line with ' i = ' line.
                  READ(indx,130) sampled_movie(i_sample), sampled_time(i_sample), sampled_energy(i_sample)
                  130 FORMAT (1X,4X,I8,9X,F12.3,6X,F20.10)
                  131 FORMAT (A4,3F20.10)
                  inner: do iatom= 1,nat
                    read (indx,131) sphere_info(iatom, i_sample)%name, sphere_info(iatom,i_sample)%coord(1), & 
                      sphere_info(iatom,i_sample)%coord(2), sphere_info(iatom,i_sample)%coord(3)
                    sphere_info(iatom, i_sample)%id = iatom
                    if (sphere_info(iatom, i_sample)%name == "O") THEN
                        sphere_info(iatom, i_sample)%mass = 16.00d0
                        sphere_info(iatom, i_sample)%radius = 0.48d0
                    ELSEIF (sphere_info(iatom, i_sample)%name == "H") THEN
                        sphere_info(iatom, i_sample)%mass = 1.00d0
                        sphere_info(iatom, i_sample)%radius = 0.53d0
                    ELSEIF (sphere_info(iatom, i_sample)%name == "N") THEN
                        sphere_info(iatom, i_sample)%mass = 14.00d0 
                        sphere_info(iatom, i_sample)%radius = 0.56d0
                    ELSEIF (sphere_info(iatom, i_sample)%name == "Li") THEN
                        sphere_info(iatom, i_sample)%mass = 6.94d0
                        sphere_info(iatom, i_sample)%radius = 1.67d0
                    ELSEIF (sphere_info(iatom, i_sample)%name == "Na") THEN
                        sphere_info(iatom, i_sample)%mass = 22.99d0
                        sphere_info(iatom, i_sample)%radius = 1.90d0
                    ELSEIF (sphere_info(iatom, i_sample)%name == "K") THEN
                        sphere_info(iatom, i_sample)%mass = 39.10d0
                        sphere_info(iatom, i_sample)%radius = 2.43d0
                    ELSEIF (sphere_info(iatom, i_sample)%name == "I") THEN
                        sphere_info(iatom, i_sample)%mass = 126.9045d0
                        sphere_info(iatom, i_sample)%radius = 1.40d0
                    endif
                    !WRITE (*,131) & 
                    !sphere_info(iatom, i_sample)%atom_name, sphere_info(iatom,i_sample)%coord(1), &
                    !sphere_info(iatom,i_sample)%coord(2), sphere_info(iatom,i_sample)%coord(3)
                  enddo inner
                  i_sample = i_sample + 1 !The position is important. It must be located before ENDIF 
              ENDIF CHECK_HEAD
          ENDIF PRE_CHECK
      END DO
  END SUBROUTINE read_traj_sphere

  SUBROUTINE read_surf_coord(indx,n_samples,n_grid,surf_info_fortran)
    IMPLICIT NONE

    INTEGER :: ierror, i_sample, i_grid
    INTEGER :: skip_num ! The number of steps skipped before reading
    INTEGER,INTENT(IN) :: indx
    INTEGER,INTENT(IN) :: n_grid ! n_grid = nb_divx * nb_divy; KEEP
    INTEGER,INTENT(IN) :: n_samples ! n_samples = INT(nmo/ns); KEEP
    REAL(KIND=8),DIMENSION(2,n_samples,n_grid),INTENT(INOUT) :: surf_info_fortran
    
    ierror =0; i_sample=0; i_grid=0; skip_num=0
  
    131 FORMAT (11X,2F13.6)
      outer: DO i_sample = 1, n_samples 
          READ(indx, '(A4,I5)',IOSTAT=ierror) 
          inner: do i_grid = 1, n_grid
              read (indx,131,IOSTAT=ierror) surf_info_fortran(1,i_sample,i_grid), & 
                  surf_info_fortran(2,i_sample,i_grid)
          enddo inner
      END DO outer
  END SUBROUTINE read_surf_coord

  SUBROUTINE skip_lines(indx, i_input)
      ! To skip lines when read data from the input
      IMPLICIT NONE
      INTEGER :: i
      INTEGER,INTENT(IN) :: i_input,indx
      do i=1,i_input
         READ(indx,*) !Neglect (nat+2)*(ns-1) lines
      enddo    
  END SUBROUTINE skip_lines

  FUNCTION hydrogen_ndx_list_for_oxygen(ndx_oxygen_1, & 
           pos_filename,natoms,boxsize) 
      implicit none
      ! 2020/07
      !========================
      !Parameters and variables
      !========================
      INTEGER,parameter :: rk=8  
      INTEGER,INTENT(IN) :: ndx_oxygen_1
      character(LEN=200),INTENT(IN) :: pos_filename        ! specific trajectory filename to analyz data
      INTEGER,INTENT(IN) :: natoms
      REAL(KIND=rk),DIMENSION(3),INTENT(IN) :: boxsize
      
      !Local 
      REAL(KIND=rk) :: distance2
      INTEGER :: get_number_of_oxygen,get_number_of_hydrogen
      INTEGER :: i,j,iatom,nmovie,& 
                    m1,m3,i_H,&
                    i1,ii,jj,i_OW, num
      INTEGER, DIMENSION(2) :: hydrogen_ndx_list_for_oxygen
      REAL(KIND=rk), ALLOCATABLE,DIMENSION (:,:) :: x,y,z
      character(LEN=3), ALLOCATABLE, DIMENSION (:) :: atom_type
      INTEGER, ALLOCATABLE, DIMENSION (:) :: ndx_OW,ndx_H
      REAL(KIND=rk), DIMENSION(3) :: r_a, r_b
      REAL(KIND=rk), parameter :: r_ohc=1.21d0   ! rOH (1.1**2)
      REAL(KIND=rk) :: r ! distance between O and H.
      ! Initialization
      i1=0
      !==================
      !Read data in input
      !==================
      nmovie=1 !Number of movie steps reqired.
      ALLOCATE(atom_type(natoms))
      ALLOCATE(x(natoms,nmovie))
      ALLOCATE(y(natoms,nmovie))
      ALLOCATE(z(natoms,nmovie))
      !=======================
      ! To obtain the # of H and O atoms
      i_H = get_number_of_hydrogen(pos_filename, natoms) ! By calling functions, we decouple different functions
      i_OW = get_number_of_oxygen(pos_filename, natoms)! By calling functions, we decouple different functions 
      !=======================================================
      ! Calculate the indices of O (H) atoms in water molecules
      !=======================================================
      ALLOCATE(ndx_OW(i_OW))    ! this should be put after i_OW is defined
      ALLOCATE(ndx_H(i_H))      ! this should be put after i_H is defined
      i=0; ii=0; jj=0
      OPEN(10,file=trim(pos_filename))
          REWIND(10)     
          READ(10,*) !Neglect data of this line. Be careful, the two 'read' lines must be located before the 'do iatom=1,natoms' loop, otherwise, there will be an error.
          READ(10,*)                  
          do iatom=1,natoms
              read (10,*)atom_type(iatom),x(iatom,nmovie),& 
                  y(iatom,nmovie),z(iatom,nmovie)
              if (trim(atom_type(iatom)) .eq. 'O')THEN
                     i=i+1      
                     ndx_OW(i)=iatom
              ELSEIF(trim(atom_type(iatom)) .eq. 'H') THEN
                     ii=ii+1
                     ndx_H(ii)=iatom
              else
              endif
          enddo 
      CLOSE(10)
      DEALLOCATE(atom_type)
      !====================
      !Producing the H list      
      !====================
      m1 = ndx_oxygen_1
      !WRITE(*,*) "m1=",m1
      !m2 = ndx_oxygen_2
      num = 0 ! For counting
      ! 
      do j = 1, 1
          do ii = 1, i_H
              m3=ndx_H(ii)
              ! I use pbc to consider the distance r
              r_a = (/x(m1,j),y(m1,j),z(m1,j)/) ! Coordinate of O
              r_b = (/x(m3,j),y(m3,j),z(m3,j)/) ! Coordinate of H
              r = distance2(r_a,r_b,boxsize) ! squared distance between O and H
              if (r < r_ohc) THEN ! If the H is bound to the O atom
                  num = num + 1
                  hydrogen_ndx_list_for_oxygen(num) = m3 ! Add the index of the H atom into the H list for the water pair (O1, O2)
              endif
          enddo
      enddo
      DEALLOCATE(ndx_OW,ndx_H,x,y,z)
      !=============================
  END FUNCTION hydrogen_ndx_list_for_oxygen
   
  SUBROUTINE make_neighbors(nat,n_samples,boxsize,neighbors_number,neighbors_indices,sphere_info)
    USE module_ihb, only: rk, sphere, num_coord_max, r_OH_LU2, r_ON_LU2, distance2
    IMPLICIT NONE

    INTEGER :: i, j
    INTEGER,INTENT(IN) :: nat, n_samples
    REAL(KIND=rk), DIMENSION(3), INTENT(IN) :: boxsize
    TYPE(sphere),DIMENSION(nat,n_samples),INTENT(INOUT) :: sphere_info
    INTEGER,DIMENSION(nat),INTENT(INOUT) :: neighbors_number
    INTEGER,DIMENSION(num_coord_max,nat),INTENT(INOUT) :: neighbors_indices
    !REAL(KIND=rk) :: dist2
    
    INTEGER,PARAMETER :: start_step = 1
    REAL(KIND=rk) :: r_LU2 ! LU: Least Upper bound; 2: squared value 

    !Initialize arrays and variables
    neighbors_number = 0
    neighbors_indices = 0
    r_LU2 = MAX(r_OH_LU2,r_ON_LU2)

    Outer: DO i = 1, nat
      Inner: DO j = 1, nat
        IF (j==i) CYCLE
        IF ( distance2(sphere_info(i,start_step)%coord, sphere_info(j,start_step)%coord, boxsize)<   &
           r_LU2                                                                                     &
           ) THEN 
           neighbors_number(i) = neighbors_number(i) + 1
           neighbors_indices(neighbors_number(i),i) = j
           !WRITE(*,*) "Atom Name:",sphere_info(i,start_step)%name
           !WRITE(*,*) "Atom id:",sphere_info(i,start_step)%id
           !WRITE(*,*) "distance:",SQRT(distance2(sphere_info(i,start_step)%coord, sphere_info(j,start_step)%coord, boxsize))
        END IF
      END DO Inner
    END DO Outer
  END SUBROUTINE make_neighbors

  SUBROUTINE obtain_indices(atom_name,indices,nat,n_samples,sphere_info)
    !===================================================================
    ! This subroutine does one thing for an chemical element (X, a string denoted as atom_name):
    ! obtain the indices of atom X and store it in the array indices
    ! indices: the output array
    !===========================
    
    USE module_ihb, only: sphere, num_coord_max
    IMPLICIT NONE
 
    CHARACTER(len=*),INTENT(IN) :: atom_name
    INTEGER,INTENT(IN) :: nat,n_samples
    TYPE(sphere),DIMENSION(nat,n_samples),INTENT(INOUT) :: sphere_info
    INTEGER,DIMENSION(nat), INTENT(INOUT) :: indices
    
    INTEGER,PARAMETER :: start_step=1 
    INTEGER :: i, n_count
    
    !WRITE(*,*) "obtain_indices:  atom_name = ", atom_name
    ! Initialize indices
    indices = 0
    ! Generate the array of index of Li,I.
    n_count = 0
    DO i = 1, nat
      IF ( trim(sphere_info(i,start_step)%name) == trim(atom_name) ) THEN
        n_count = n_count + 1
        indices(n_count) = i
        !WRITE(*,*)"indices(n_count): ", indices(n_count)
      ENDIF
    END DO
  END SUBROUTINE obtain_indices

  SUBROUTINE make_wat_neighbors_for_X_shell(atom_name,step_t,nat,n_samples,boxsize, &
             wat_neighbors_number,wat_neighbors_indices,sphere_info)
    !===================================================================
    ! this subroutine do one thing:
    ! For an atom (index = i_ion), at time t (denoted by step_t), 
    ! 1: to find the number (N) of its neighboring water molecules and store N in wat_neighbors_number(i_ion); 
    ! 2: to add the index of O of this water molecule into the wat_neighbors_indices(ind, i_ion)
    !====================================================================
    
    USE module_ihb, only: rk,sphere,num_coord_mol_max,r_OH_LU2,r_NN_OW_LU2,r_ON_OW_LU2,r_Li_OW_LU2,&
      r_Na_OW_LU2,r_K_OW_LU2,r_I_OW_LU2,r_OW_OW_LU2, distance2
    IMPLICIT NONE

    INTEGER :: i, j, k
    !INTEGER :: n_count  ! To count the number of neighbors of X-ion/atom/etc.
    CHARACTER(len=*),INTENT(IN) :: atom_name
    INTEGER,INTENT(IN) :: nat,n_samples
    REAL(KIND=rk),DIMENSION(3),INTENT(IN) :: boxsize
    TYPE(sphere),DIMENSION(nat,n_samples),INTENT(INOUT) :: sphere_info
    INTEGER,DIMENSION(nat),INTENT(INOUT) :: wat_neighbors_number
    INTEGER,DIMENSION(num_coord_mol_max,nat),INTENT(INOUT) :: wat_neighbors_indices
     
    ! Local variable
    INTEGER,DIMENSION(nat) :: indices
    !REAL(KIND=rk) :: dist2
    
    INTEGER,INTENT(IN):: step_t 
    REAL(KIND=rk) :: r_LU2 ! LU: Least Upper bound; 2: squared value 

    !Initialize arrays and variables
    i=0;j=0;k=0
    wat_neighbors_number = 0
    wat_neighbors_indices = 0

    !================    
    ! Determine r_LU2
    !================    
    IF (atom_name =="Li") THEN
      r_LU2 = r_Li_OW_LU2
    ELSE IF (atom_name =="Na") THEN
      r_LU2 = r_Na_OW_LU2
    ELSE IF (atom_name =="K") THEN
      r_LU2 = r_K_OW_LU2
    ELSE IF (atom_name =="I") THEN
      r_LU2 = r_I_OW_LU2
    ELSE IF (atom_name =="N") THEN
      r_LU2 = r_NN_OW_LU2
    END IF

    ! Generate the array of index of Li, or I, etc.
    WRITE(*,*) "make_wat_neighbors_for_X_shell:  atom_name = ", atom_name
    CALL obtain_indices(atom_name,indices,nat,n_samples,sphere_info)
    !Test
    !WRITE(*,*) "indices of atom_name: ", indices
    Outer: DO i = 1, nat
      IF(indices(i) == 0) EXIT
      k = indices(i)
      !WRITE(*,*) "i=",i, "; indices(i)=", k
      Inner: DO j = 1, nat
      IF (j== k ) CYCLE
      IF (TRIM(sphere_info(j,step_t)%name)=="O" .AND. &
          distance2(sphere_info(k,step_t)%coord, sphere_info(j,step_t)%coord, boxsize)< r_LU2 &
         ) THEN 
         wat_neighbors_number(i) = wat_neighbors_number(i) + 1
         !If the number of neighbors are too large, we neglect the extral neighbors (Here, they are water molecules.)
         IF (wat_neighbors_number(i) > num_coord_mol_max) EXIT
         wat_neighbors_indices(wat_neighbors_number(i),i) = j ! I have to separate the O, and H indies.
         !WRITE(*,*) "ATOM Name (k):", sphere_info(k,step_t)%name
         !WRITE(*,*) "ATOM Name (j):", sphere_info(j,step_t)%name
         !WRITE(*,*) "Atom id (j):", sphere_info(j,step_t)%id
         !WRITE(*,*) "distance (k,j):", SQRT(distance2(sphere_info(k,step_t)%coord, sphere_info(j,step_t)%coord, boxsize))
      END IF
      END DO Inner
      !WRITE(*,*) "wat_neighbors_number = ", wat_neighbors_number(i)
      !WRITE(*,*) "wat_neighbors_indices:", wat_neighbors_indices(:,i)
    END DO Outer
  END SUBROUTINE make_wat_neighbors_for_X_shell

  CHARACTER(LEN=200) FUNCTION name_list_OH_at_shell_interface(filename,core_atom_str)
    IMPLICIT NONE

    CHARACTER(len=*),INTENT(IN) :: filename 
    character(LEN=*),INTENT(IN) :: core_atom_str        ! For example, core_atom_str can be "I","Na", etc.

    name_list_OH_at_shell_interface=trim(filename)//"_"//core_atom_str//'_OH_at_shell_interface_list.dat'
  END FUNCTION

  CHARACTER(LEN=200) FUNCTION name_c2_at_shell_interface(filename,core_atom_str)
    IMPLICIT NONE

    CHARACTER(len=*),INTENT(IN) :: filename 
    character(LEN=*),INTENT(IN) :: core_atom_str        ! For example, core_atom_str can be "I","Na", etc.

    name_c2_at_shell_interface=trim(filename)//"_"//core_atom_str//'_c2_at_shell_interface.dat'
  END FUNCTION
