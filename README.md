## Top-to-Down 

Date                          Description
2020-7-11                 Steps to do

Language                  Fortran 95

### 0. Describe the problem

I want to calculate the correlation function $<h(0)h(t)>/<h>$, where $h=1$ when a **pair of water molecules** are H-bonded, and 0 otherwise. 
Note that there are 4 bonding forms for a pair of water molecules. 

some parameters:

```fortran
R_{oh} = 1.1  !Angstrom
```

### 1 Define the input and output
​    input: 

​           the trajectory file `pos.xyz`.
​    output: 

​           $h(n_{pair},t)$:   

​           C(t) for all the sampled time points. It is an 1-D array with the shape (STEPS,),
​           where STEPS is the number of sampled time points.

### 2 Design the algorithm which should be implemented in this program.

(Done) First, one have to read the trajectory file :

`1a_sample.f95`

`1b_label.f95`

Have to define a new type called `water_pair`, and define a 2-D array  `water_pair_info(num_water_pairs, imovie)`

```fortran
USE water_pair_module

! Purpose:
! To define the derived data type for a water pair
IMPLICIT NONE
TYPE :: water_pair
  INTEGER :: water_pair_id
  INTEGER :: head_id ! For example, 2
  INTEGER :: tail_id ! For example, 1
  TYPE (water) :: head ! water molecule
  TYPE (water) :: tail  !water molecule
  LOGICAL :: h_state 
  REAL :: h ! Hydrogen bond population operator 1 or 0 
  REAL :: hb_length 
END TYPE water_pair

!The array ```water_pair_info``` can be shared by subroutines  
TYPE(water_pair), ALLOCATABLE, DIMENSION(:,:) :: water_pair_info
END MODULE water_pair_module  
```

Naturally, we need  two types `water` and `atom`. The type `atom` is defined in `atom_module.f95` as follows.

```fortran
MODULE atom_module
! Purpose:
!   To define the derived data type for atom
IMPLICIT NONE
TYPE :: atom
  CHARACTER(LEN=2) :: atom_name
  INTEGER :: atom_id 
  INTEGER :: host_id  ! For O atom in water, host_id = atom_id
  REAL :: mass
  REAL, DIMENSION(3) :: coord 
END TYPE atom

! The array atom_info can be shared by subroutines  
TYPE(atom), ALLOCATABLE, DIMENSION(:,:) :: atom_info
END MODULE atom_module
```

The type `water_module`.

```fortran
MODULE water_module
! Purpose:
! To define the derived data type for water
IMPLICIT NONE
TYPE :: water
  INTEGER :: molecule_id
  TYPE (atom):: oxygen
  TYPE (atom):: hydrogen1
  TYPE (atom):: hydrogen2
END TYPE water

! The array water_info can be shared by subroutines  
TYPE(water), ALLOCATABLE, DIMENSION(:,:) :: water_info
END MODULE water_module
```

### 2. Create a list file for all O-O pairs



#### 2a For any atom, find out which molecule it belong to  

​      2a1) host molecule's symbol: `symb_host (indics of molecules must be defined)`, define a molecular ID number

​     2a2) For any atom, find out its host (group), we have to define an attribute for atom 

​     2a2a) PBC is considered before define the distance between the host atom and the guest atom. The PBC can be considered in this way

```fortran
if (dx = x(host) - x > boxsize(1)/2 ) then  
    dx = dx - boxsize(1) 
```

If `distance2()` betweeen H and O atom is less than `R_oh`,  then `H.host_id = O.atom_id`.

#### 2b For any H atom, find out which water molecule is it belong to

The same to 2a.



Calculate O-O distance by  `water_pair_info(k,imovie).head%oxygen%coord(1)` etc.,

I need a function called `distance2(r1,r2,size) `, in which I have to consider the PBC. Here `r1` is a 1D array with 3 elements 

```fortran
(/water_pair_info(k,imovie).head%oxygen%coord(1), water_pair_info(k,imovie).head%oxygen%coord(2), water_pair_info(k,imovie).head%oxygen%coord(3)/) 
```

And `r2` is

```fortran
(/water_pair_info(k,imovie).tail%oxygen%coord(1), water_pair_info(k,imovie).tail%oxygen%coord(2), water_pair_info(k,imovie).tail%oxygen%coord(3)/)
```

`boxsize` is a 1-D array with 3 elements `boxsize = (/a,b,c/)`.

The function can be defined as follows.

```fortran
! This function calculate the squared distance between r1,and r2.
REAL FUNCTION distance2(r1,r2,boxsize)
          real(kind=rk), DIMENSION(1:3), INTENT(IN) :: r1, r2, boxsize
          REAL(kind=rk) :: dx,dy,dz
          dx = r1(1) - r2(1)
          if (abs(dx) > boxsize(1)*0.5) then
              dx = boxsize(1) - dx
          endif
          dy = r1(2) - r2(2)
          if (abs(dy) > boxsize(2)*0.5) then
              dy = boxsize(2) - dy
          endif
          dz = r1(3) - r2(3)
          if (abs(dz) > boxsize(3)*0.5) then
              dz = boxsize(3) - dz
          endif
          distance2 = dx**2 + dy**2 + dz**2
      END FUNCTION distance2
```



We also need a function to calculate the cosine of the angle formed by r1, r2, r3, in which PBC also has to be considered.  

```fortran
! In modules_tools.f95
      REAL FUNCTION diff_axis(u1,u2,h)
          ! u2 is used as origin
          real(kind=rk),INTENT(IN) :: u1, u2
          real(kind=rk),INTENT(IN) :: h
          real(kind=rk) :: du 
          du = u1-u2
          if (abs(du) > 0.5*h) then
              du = h - du 
          endif
          diff_axis = du
      END FUNCTION diff_axis
      
      REAL FUNCTION pmADH(r1,r2,r3,boxsize)
          real(kind=rk),DIMENSION(3), INTENT(IN) :: r1,r2,r3
          real(kind=rk),DIMENSION(3), INTENT(IN) :: boxsize 
          pmADH= diff_axis(r3(1),r2(1),boxsize(1))*           &
                 diff_axis(r1(1),r2(1),boxsize(1))+      &
                 diff_axis(r3(2),r2(2),boxsize(2))*      &
                 diff_axis(r1(2),r2(2),boxsize(2))+      &
                 diff_axis(r3(3),r2(3),boxsize(3))*      &
                 diff_axis(r1(3),r2(3),boxsize(3))       ! pm: point multiplication. 
      END FUNCTION pmADH

      REAL FUNCTION pmAHD(r1,r2,r3,boxsize)
          real(kind=rk),DIMENSION(3), INTENT(IN) :: r1,r2,r3
          real(kind=rk),DIMENSION(3), INTENT(IN) :: boxsize 
          pmAHD= diff_axis(r1(1),r3(1),boxsize(1))*           &
                 diff_axis(r2(1),r3(1),boxsize(1))+      &
                 diff_axis(r1(2),r3(2),boxsize(2))*      &
                 diff_axis(r2(2),r3(2),boxsize(2))+      &
                 diff_axis(r1(3),r3(3),boxsize(3))*      &
                 diff_axis(r2(3),r3(3),boxsize(3)) ! pm: point multiplication. 
      END FUNCTION pmAHD 
      
      ! Then one can obtain cosPhi as follows.
      r21 = distance2(r1,r2,boxsize)
      r31 = distance2(r1,r3,boxsize)
      cosPhi = pmADN(r1,r2,r3,boxsize) /SQRT(r21 * r23) 
```

 We need to use a DO LOOP to check all the 4 forms of HB link between the two water molecules in the $k-$th water pair. At least one of them satisfies the geometry criterion of HB, we can set  (`water_pair_info(k,imovie).h = 1. `  The condition is 

```fortran
distance2(r1,r2,boxsize) < r21C .AND.
pm213(r2,r1,r3a,boxsize) < cosPhiC .OR.  
pm213(r2,r1,r3b,boxsize) < cosPhiC .OR. 
pm213(r1,r2,r3c,boxsize) < cosPhiC .OR.
pm213(r1,r2,r3d,boxsize) < cosPhiC 
```

where we have named the four Hydrogen atoms as a,b,c and d.

Consider the $k$-th water pair, check if it satisfy the geometry criterion of HB, if it does, $h=1$,  (`water_pair_info(k,imovie).h =1. `)  otherwise, $h=0$. (`water_pair_info(k,imovie).h =0. `)

2b1) PBC is also considered before define the distance between two atoms (eg, Donar and Acceptor).

 

This step is the same to the 2a2a-th step. 



 ### 3 Calculate the distance between the water molecule pairs in the object water-pair.
​    input: 

​          cosphiC
​          r21C 
​          r31C

Now that We have $h(k,t )$ stored in a 2-D array. The size is: $n_p \times n_{sample}$.

Here $n_p$ is the the number of pairs of water molecules.

### 4 Calculate the correlation function $C(t)$:

​     output: h(t) an 1-D array with nmo elements.
​            $$C(t) = <h(0)h(t)> / <h>    $$   

### 4A Calculate the correlation function $S(t)$: 

output: h(t) an 1-D array with nmo elements.
            $$S(t) = <h(0)H(t)> / <h> $$   

First, we can write a subroutine `ghbacf_s_pbc.f95` for pure water, and `ghbacf_s_pbc_format2.f95` for LiI solution system, etc.



### 5 Translate the algorithm to Fortran codes

### 6 Test the codes.

### 7 Finish up.



## Tutorial

#### (1) General HB

The input file: Here is an example.

```input_main```:

```shell
15.6404 15.6404 15.6404
0.001
128w
128w-pos_2delta_t.xyz
0
60000
384
1000
1
```

To compile this program, run

```
./main_compile.sh
```

To run this program:

```shell
./main< input_main_128w  > hbacf_pair_water.out
```

Or, one can compile and run:

```
./main_compile_and_run.sh
```

#### (2) Interfacial HB

If one want to calculate correlation functions about interfacial HBs, one has to use an additional file: `surf_filename`,  i.e., the input file looks like:

```shell
15.6404 15.6404 31.2808 
0.0005
128w_itp
traj_pos_128w_itp.xyz
0
80000
384
80
1
surf_traj.dat
```

Note:  In the current version, the `surf_traj.dat` is obtained from the Python code `1_chandler_fast.py`. We 'd better translate it into Fortran later.



If one want to calculate correlation functions about interfacial HBs,  the command to run is similar. For example, first, compile the main code by the script

```shell
# main_interface_n_compile.sh
gfortran -Wall -fcheck=all module_tools.f95 types.f95 parameter_shared.f95 atom_module.f95 surf_module.f95 water_module.f95 water_pair_module.f95 count_time.f95 traj_format2.f95 surf_traj.f95 read_interface_input.f95 1b_sample_and_recenter_format2.f95 ghbond.f95 ghbacf_interface_n_pbc.f95 ghbacf_k_pbc_format2.f95 main_interface_n.f95 -o main_interface_n
 
rm *.mod
```

with the command

```shell
./main_interface_n_compile.sh
```

Then, one can run the main code by 

```shell
./main_interface_n < input_main_sample_128w_itp
```

### 2. Calculate $C_{HB}(t)$ of I$^-$-water hydrogen bonds

(1). Edit the input file: `input_main_lii_new `

```shell
16.0963 16.0963 16.0963
0.0005
LiI_new
traj_lii_new.xyz
73001
146300
376
40
1
```

(2). Then run the `main_lii` function

```shell 
./main_lii < input_main_lii_new
```

This function can return all water-water hydrogen bond correlation function $C_{HB}(t)$ of water-water hydrogen bonds and  I$^-$-water hydrogen bonds.

### 3. Least Square Fits 

See `/home/gang/Github/hbacf/__hbacf_continuous/correlations/least_square_fit/128w_itp_pure_ihb`
