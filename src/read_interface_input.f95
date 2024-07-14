SUBROUTINE read_interface_input(boxsize,delta_t0,filename,pos_filename,nmo_start,nmo_end,nat,ns,&
               criterion,surf_filename,thickness) 
  !2018/12/27
  ! 0) Purpose:     
  ! 1a) The subroutine read_interface_input.f95 read the PARAMETERs from the input file. 
  ! One important input file is: the surface trajectory file. (And its name should be given to argument surf_filename.)
  ! The surface trajectory file is obtained from the program XXX 
  ! 0a) Date    Programmer     version
  ! 20181225    Gang Huang     3c
  ! Declaration
  IMPLICIT NONE

  !==========
  !PARAMETERs
  !==========
  INTEGER,PARAMETER :: rk=8              

  !=========
  !Variables
  !=========
  REAL(KIND=rk),DIMENSION(3),INTENT(INOUT) :: boxsize  
  INTEGER,INTENT(INOUT) :: criterion
  REAL(KIND=rk),INTENT(INOUT) :: delta_t0  ! For reading data
  CHARACTER(LEN=200), INTENT(INOUT) :: filename
  INTEGER :: i,iatom,ierr,imovie
  CHARACTER(LEN=50) :: line 
  INTEGER,INTENT(INOUT) :: nat ! number of atoms
  INTEGER,INTENT(INOUT) :: nmo_start 
  INTEGER,INTENT(INOUT) :: nmo_end
  INTEGER,INTENT(INOUT) :: ns
  CHARACTER(LEN=200), INTENT(INOUT) :: pos_filename
  CHARACTER(LEN=200), INTENT(INOUT) :: surf_filename
  REAL(KIND=rk), INTENT(INOUT) :: thickness ! the thickness of the instantaneous interfaces

  !==============
  !Initialization
  !==============
  iatom=0; imovie=0; i=0

  !==================
  !read data in input
  !==================
  WRITE(6,*)'What is the size of box (a,b,c):'
  READ(*,'(a)') line
  if (len_trim(line)==0) then
      boxsize(1) = 5.32d0
      boxsize(2) = 5.32d0
      boxsize(3) = 5.32d0
  else
    READ(line, *, iostat = ierr) boxsize(1), boxsize(2), boxsize(3)
    WRITE(*,*)"boxsize:", boxsize(1), boxsize(2), boxsize(3)
  endif

  WRITE(6,*)'What is the time step in the traj. file (ps): (Default: 0.d0005)'
  READ(*,'(a)') line
  if (len_trim(line)==0) then
    delta_t0 = 0.0005d0
  else
    READ(line, *, iostat = ierr) delta_t0
    WRITE(*,*) delta_t0
  endif

  WRITE(6,*)'What is the name of the system:'
  READ(*,'(a)') line
  if (len_trim(line)==0) then
    filename = "8w_bulk_pbc"
  else
    READ(line, *, iostat = ierr) filename
    WRITE(*,*) "filename: ",filename
  endif

  WRITE(6,*)'What is the name of the trajecotry file:'
  READ(*,'(a)') line
  if (len_trim(line)==0) then
    pos_filename = "traj_pos_8w.xyz"
  else
    READ(line, *, iostat = ierr) pos_filename
    WRITE(*,*) "pos_filename: ",pos_filename
  endif

  WRITE(6,*)'What is the initial step of the trajecotry:'
  READ(*,'(a)') line
  if (len_trim(line)==0) then
    nmo_start = 0
  else
    READ(line, *, iostat = ierr) nmo_start
    WRITE(*,*) "nmo_start: ", nmo_start
  endif

  WRITE(6,*)'What is the end step of the trajecotry:'
  READ(*,'(a)') line
  if (len_trim(line)==0) then
    nmo_end = 0
  else
    READ(line, *, iostat = ierr) nmo_end
    WRITE(*,*) "nmo_end: ", nmo_end
  endif

  WRITE(6,*)'What is the total number of atoms in the system:'
  READ(*,'(a)') line
  if (len_trim(line)==0) then
    nat = 18
  else
    READ(line, *, iostat = ierr) nat
    WRITE(*,*) "nat:", nat
  endif

  WRITE(6,*)'What is the time step for calculating CORRELATION:' 
  ! [ns*0.d0005] ps is the new time step for calculating correl func.
  READ(*,'(a)') line
  if (len_trim(line)==0) then
    ns = 40
  else
    READ(line, *, iostat = ierr) ns
    WRITE(*,*) "ns: ", ns
  endif

  WRITE(6,*)'What is the criterion of HB:' 
  READ(*,'(a)') line
  if (len_trim(line)==0) then
    criterion = 1
  else
    READ(line, *, iostat = ierr) criterion
    WRITE(*,*) "HB criterion: ", criterion
    WRITE(*,*) "criterion 1 denote ADH and 2 denotes AHD definition. "
  endif  

  WRITE(6,*)'What is the name of the surface trajectory file:' 
  READ(*,'(a)') line
  if (len_trim(line)==0) then
    surf_filename = "surf_traj.dat"
  else
    READ(line, *, iostat = ierr) surf_filename
    WRITE(*,*) "surf traj. name: ", surf_filename
  endif  

  WRITE(6,*)'What is the thickness of the interface you want to define: (Angstrom)' 
  READ(*,'(a)') line
  if (len_trim(line)==0) then
    thickness = 3.d0
  else
    READ(line, *, iostat = ierr) thickness
    WRITE(*,*) "READ thickness: ", thickness
  endif  

END SUBROUTINE read_interface_input 
