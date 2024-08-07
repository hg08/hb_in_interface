PROGRAM main_interface
  ! 0) Purpose:     
  !    Determine the HB dynamics of interface layer of neat water, by using Sulpizi's selection algorithm.
  !
  ! 0a) Date     Programmer     version
  !  ======      ==========     =======
  !  20/9/23   Gang Huang     Original code
  !============
  ! Declaration
  !============
  USE module_ihb
  IMPLICIT NONE

  !==========
  !parameters
  !==========
  !integer, parameter :: rk=8              

  ! The array atom_info can be shared by subroutines  
  TYPE(atom), ALLOCATABLE, DIMENSION(:,:) :: atom_info
  !To declear data for calculating the total time cosummed.
  integer :: begin_time,end_time,rat
  !To declare data to share between routines.
  character(LEN=200) :: sampled_pos_filename
  INTEGER, ALLOCATABLE, DIMENSION(:) :: sampled_movie
  REAL(KIND=rk), ALLOCATABLE, DIMENSION(:) :: sampled_time, sampled_energy
  ! For griding
  !REAL(KIND=rk), PARAMETER :: whish_size=0.5d0! Angstrom
  INTEGER :: nb_divx, nb_divy, nb_divz, n_grid 
  REAL(KIND=rk) :: divx, divy, divz

  REAL(KIND=rk) :: thickness ! the thickness of the instantaneous interfaces
  REAL(KIND=rk), ALLOCATABLE, DIMENSION(:,:,:) :: surf_info 

  CHARACTER(LEN=2) :: atom_type
  REAL(KIND=rk), dimension(3) :: boxsize
  integer :: criterion
  REAL(kind=rk) :: delta_t0 ! For reading data
  character(LEN=200) :: filename
  CHARACTER(LEN=2) :: guest_atom
  CHARACTER(LEN=2) :: host_atom
  INTEGER :: i, iatom, imovie
  !To save the indices of the molecules for generating list file, we define an array for each time point (jj, in this code)
  INTEGER, ALLOCATABLE, DIMENSION(:,:) :: indx_array
  character(LEN=200) :: list_oxygen_pairs
  INTEGER :: nat ! number of atoms
  INTEGER :: nmo_start ! starting step index
  INTEGER :: nmo_end ! end step index
  INTEGER :: ns ! Get one sample from the trajectory every ns step.
  INTEGER :: ns_2nd
  INTEGER :: n_samples ! n_samples = INT(nmo/ns)
  INTEGER :: nwat ! number of water molecules
  character(LEN=200) :: surf_filename
  character(LEN=200) :: pos_filename
  
  !==============
  !Initialization
  !==============
  iatom=0;imovie=0;i=0
  atom_type=''
  guest_atom="H"
  host_atom="O"
  boxsize=(/0.d0,0.d0,0.d0/)
  nat=0
  nwat=nat/3
  criterion=1 ! 1 means ADH criterion of H-Bbond definition
  filename=""; pos_filename=""; sampled_pos_filename=""
  surf_filename=""
  list_oxygen_pairs=""
  thickness=1.d0
  ns_2nd = 1

  call system_clock(begin_time,rat) !To get the starting time
  
  !============================================
  ! To read the required controlling parameters
  !============================================
  CALL read_interface_input(boxsize,delta_t0,filename,pos_filename,nmo_start,nmo_end,nat,ns,&
          criterion,surf_filename,thickness) 

  ! Obatin n_samples
  n_samples = sampling_number(nmo_start, nmo_end,ns)
  allocate(sampled_movie(n_samples))
  allocate(sampled_time(n_samples))
  allocate(sampled_energy(n_samples))
  allocate(atom_info(nat,n_samples))
 
  !========================
  ! Sampling the trajectory
  !========================
  !CASE1: If one does not need to recenter, one can just call sample_format2()
  !CALL sample_format2(pos_filename,nmo_start,nmo_end,nat,ns,n_samples)
  !CASE2: If one have to recenter, one call sample_and_recenter_format2() instead.
  WRITE(*,*) "sample and recenter_format2 starting..."
  CALL sample_and_recenter_format2(pos_filename,nmo_start,nmo_end,nat,ns,n_samples,boxsize,&
       sampled_pos_filename,sampled_movie,sampled_time,sampled_energy, &
       nb_divx,nb_divy,nb_divz,n_grid,divx,divy,divz,whish_size,atom_info)

  ! After running the sample() or sample_format2() subroutine, therefore we have a new sampled trajectory file (stored in atom_info), 
  ! which is generally shorter than the original one.
  
  !====================
  !read surf trajectory
  !====================
  allocate(surf_info(2,n_grid,n_samples))
  surf_info = 0
  !WRITE(*,*) "SHAPE(surf_info)= ", SHAPE(surf_info)
  ns_2nd = 1 ! sample freq is 1, ie., all data are sampled
  WRITE(*,*) "read_surf_traj is starting..."
  CALL read_surf_traj(surf_filename,nmo_start,nmo_end,ns_2nd,n_grid,n_samples,surf_info)
  
  ! Use array instead of linked list, it may be faster. 
  allocate(indx_array(n_samples,nat))
  indx_array = 0

  WRITE(*,*) "molecules_in_interface() is starting..."
  CALL molecules_in_interface(n_samples,nat,indx_array,atom_info,&
     n_grid,divx,divy,divz,nb_divx,nb_divy,nb_divz,thickness,surf_info)

  !To determine the indices of Oxygens' pairs that located in one of the interfaces.
  WRITE(*,*) "ghbond_interface() is starting..."
  CALL ghbond_interface(filename,list_oxygen_pairs,n_samples,nat,indx_array)

  !Calculate n_HB(t),k(t),etc for pure water system. If the format of data is different, one may use another funtion, eg ghbacf_n_pbc_format2().
  WRITE(*,*) "c is starting..."
  CALL ghbacf_interface_c_pbc_format2(boxsize,delta_t0,filename,pos_filename,list_oxygen_pairs, &
                    n_samples,nat,ns,criterion,atom_info,n_grid,divx,divy,divz,nb_divx,nb_divy,&
                    nb_divz,thickness,surf_info)
  WRITE(*,*) "n is starting..."
  CALL ghbacf_interface_n_pbc_format2(boxsize,delta_t0,filename,pos_filename,list_oxygen_pairs, &
                    n_samples,nat,ns,criterion,atom_info,n_grid,divx,divy,divz,nb_divx,nb_divy,&
                    nb_divz,thickness,surf_info)
  WRITE(*,*) "k is starting..."
  CALL ghbacf_interface_k_pbc_format2(boxsize,delta_t0,filename,pos_filename,list_oxygen_pairs, &
                    n_samples,nat,ns,criterion,atom_info,n_grid,divx,divy,divz,nb_divx,nb_divy,&
                    nb_divz,thickness,surf_info)

  call system_clock(end_time,rat)
  WRITE(6, *)"elapsed time: ", real(end_time-begin_time)/real(rat) 

END PROGRAM main_interface 
