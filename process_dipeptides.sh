source /cluster/tufts/ysl8/jovan/gromacs_m4/gromacs/bin/GMXRC
export PLUMED_KERNEL=/cluster/tufts/ysl8/jovan/gromacs_gpu/plumed/lib/libplumedKernel.so

module load gcc/4.9.2
module load cuda/8.0.44
module load openmpi/1.8.2

/cluster/tufts/ysl8/jovan/1211_Jennifer/s1/

gmx_mpi trjconv -f prod.xtc -s start.tpr -o prod_all.xtc -pbc mol
gmx_mpi  angle -f prod_all.xtc -n phipsichi.ndx -or prod_angle.trr -type dihedral
gmx_mpi  angle -f prod_all.xtc -n phipsichi.ndx -type dihedral -all -ov prod_angle.xvg

/cluster/tufts/ysl8/jovan/1211_Jennifer/s2/
gmx_mpi trjconv -f prod.xtc -s start.tpr -o prod_all.xtc -pbc mol
gmx_mpi  angle -f prod_all.xtc -n phipsichi.ndx -or prod_angle.trr -type dihedral
gmx_mpi  angle -f prod_all.xtc -n phipsichi.ndx -type dihedral -all -ov prod_angle.xvg


/cluster/tufts/ysl8/jovan/1211_Jennifer/s1/bemeta/

gmx_mpi trjconv -f prod1.xtc -s start1.tpr -o prod1_50_100.xtc -b 50001 -e 100000 -pbc mol
gmx_mpi trjconv -f prod2.xtc -s start2.tpr -o prod2_50_100.xtc -b 50001 -e 100000 -pbc mol
gmx_mpi trjconv -f prod3.xtc -s start3.tpr -o prod3_50_100.xtc -b 50001 -e 100000 -pbc mol
gmx_mpi trjconv -f prod4.xtc -s start4.tpr -o prod4_50_100.xtc -b 50001 -e 100000 -pbc mol
gmx_mpi trjconv -f prod5.xtc -s start5.tpr -o prod5_50_100.xtc -b 50001 -e 100000 -pbc mol

gmx_mpi trjcat -f prod?_50_100.xtc -cat -nosort -o s1Leu_all.xtc

gmx_mpi  angle -f s1Leu_all.xtc -n phipsichi.ndx -or s1Leu_all.trr -type dihedral
gmx_mpi  angle -f s1Leu_all.xtc -n phipsichi.ndx -type dihedral -all -ov s1Leu_all.xvg


/cluster/tufts/ysl8/jovan/1211_Jennifer/s1/bemeta/

gmx_mpi trjconv -f prod1.xtc -s start1.tpr -o prod1_50_100.xtc -b 50001 -e 100000 -pbc mol
gmx_mpi trjconv -f prod2.xtc -s start2.tpr -o prod2_50_100.xtc -b 50001 -e 100000 -pbc mol
gmx_mpi trjconv -f prod3.xtc -s start3.tpr -o prod3_50_100.xtc -b 50001 -e 100000 -pbc mol
gmx_mpi trjconv -f prod4.xtc -s start4.tpr -o prod4_50_100.xtc -b 50001 -e 100000 -pbc mol
gmx_mpi trjconv -f prod5.xtc -s start5.tpr -o prod5_50_100.xtc -b 50001 -e 100000 -pbc mol

gmx_mpi trjcat -f prod?_50_100.xtc -cat -nosort -o s2Leu_all.xtc

gmx_mpi  angle -f s2Leu_all.xtc -n phipsichi.ndx -or s2Leu_all.trr -type dihedral
gmx_mpi  angle -f s2Leu_all.xtc -n phipsichi.ndx -type dihedral -all -ov s2Leu_all.xvg