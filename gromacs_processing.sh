source /cluster/tufts/ysl8/jovan/gromacs_m4/gromacs/bin/GMXRC
export PLUMED_KERNEL=/cluster/tufts/ysl8/jovan/gromacs_gpu/plumed/lib/libplumedKernel.so

module load gcc/4.9.2
module load cuda/8.0.44
module load openmpi/1.8.2

/cluster/tufts/ysl8/jovan/GNSRV_fixed/s1/bemeta/
/cluster/tufts/ysl8/jovan/GNSRV_fixed/s2/bemeta/
gmx_mpi trjconv -f prod10.xtc -s start10.tpr -o prod10_50_100.xtc -b 50001 -e 100000 -pbc mol
gmx_mpi trjconv -f prod11.xtc -s start11.tpr -o prod11_50_100.xtc -b 50001 -e 100000 -pbc mol
gmx_mpi trjconv -f prod12.xtc -s start12.tpr -o prod12_50_100.xtc -b 50001 -e 100000 -pbc mol
gmx_mpi trjconv -f prod13.xtc -s start13.tpr -o prod13_50_100.xtc -b 50001 -e 100000 -pbc mol
gmx_mpi trjconv -f prod14.xtc -s start14.tpr -o prod14_50_100.xtc -b 50001 -e 100000 -pbc mol

gmx_mpi trjcat -f prod1?_50_100_ns.xtc -cat -nosort -o s1cGNSRV_all.xtc


gmx_mpi  angle -f s1cGNSRV_all.xtc -n phipsichi.ndx -or GNSRV_s1_angle.trr -type dihedral

gmx_mpi  angle -f s1cGNSRV_all.xtc -n phipsichi.ndx -type dihedral -all -ov GNSRV_s1_angle.xvg


## Cluster analysis

copied from: /cluster/tufts/ysl8/jovan/GNSRV_fixed/dPCA/50_100_ns/cluster_traj/s1cGNSRV/

gmx_mpi angle -f cluster1.xtc -n phipsichi.ndx -type dihedral -all -ov s1_cluster1.xvg
gmx_mpi angle -f cluster2.xtc -n phipsichi.ndx -type dihedral -all -ov s1_cluster2.xvg
gmx_mpi angle -f cluster3.xtc -n phipsichi.ndx -type dihedral -all -ov s1_cluster3.xvg
gmx_mpi angle -f cluster4.xtc -n phipsichi.ndx -type dihedral -all -ov s1_cluster4.xvg
gmx_mpi angle -f cluster5.xtc -n phipsichi.ndx -type dihedral -all -ov s1_cluster5.xvg
gmx_mpi angle -f cluster6.xtc -n phipsichi.ndx -type dihedral -all -ov s1_cluster6.xvg
gmx_mpi angle -f cluster7.xtc -n phipsichi.ndx -type dihedral -all -ov s1_cluster7.xvg


# Hbond analysis
gmx_mpi hbond -noda -f cluster1.xtc -s start10.tpr -n hbond_G1_index.ndx -num hbond_c1_G1.xvg
gmx_mpi hbond -noda -f cluster1.xtc -s start10.tpr -n hbond_N2_index.ndx -num hbond_c1_N2.xvg
gmx_mpi hbond -noda -f cluster1.xtc -s start10.tpr -n hbond_S3_index.ndx -num hbond_c1_S3.xvg
gmx_mpi hbond -noda -f cluster1.xtc -s start10.tpr -n hbond_R4_index.ndx -num hbond_c1_R4.xvg
gmx_mpi hbond -noda -f cluster1.xtc -s start10.tpr -n hbond_V5_index.ndx -num hbond_c1_V5.xvg
gmx_mpi hbond -noda -f cluster1.xtc -s start10.tpr -n hbond_S3S3_index.ndx -num hbond_c3_S3S3.xvg

gmx_mpi hbond -noda -f cluster2.xtc -s start10.tpr -n hbond_G1_index.ndx -num hbond_c2_G1.xvg
gmx_mpi hbond -noda -f cluster2.xtc -s start10.tpr -n hbond_N2_index.ndx -num hbond_c2_N2.xvg
gmx_mpi hbond -noda -f cluster2.xtc -s start10.tpr -n hbond_S3_index.ndx -num hbond_c2_S3.xvg
gmx_mpi hbond -noda -f cluster2.xtc -s start10.tpr -n hbond_R4_index.ndx -num hbond_c2_R4.xvg
gmx_mpi hbond -noda -f cluster2.xtc -s start10.tpr -n hbond_V5_index.ndx -num hbond_c2_V5.xvg

gmx_mpi hbond -noda -f cluster3.xtc -s start10.tpr -n hbond_G1_index.ndx -num hbond_c3_G1.xvg
gmx_mpi hbond -noda -f cluster3.xtc -s start10.tpr -n hbond_N2_index.ndx -num hbond_c3_N2.xvg
gmx_mpi hbond -noda -f cluster3.xtc -s start10.tpr -n hbond_S3_index.ndx -num hbond_c3_S3.xvg
gmx_mpi hbond -noda -f cluster3.xtc -s start10.tpr -n hbond_R4_index.ndx -num hbond_c3_R4.xvg
gmx_mpi hbond -noda -f cluster3.xtc -s start10.tpr -n hbond_V5_index.ndx -num hbond_c3_V5.xvg

# gmx_mpi hbond -f cluster1.xtc -s start10.tpr -num hbond_c1_G1.xvg -contact
# gmx_mpi hbond -noda -f cluster1.xtc -s start10.tpr -num hbond_c1_N2.xvg -contact 
# gmx_mpi hbond -noda -f cluster1.xtc -s start10.tpr -num hbond_c1_S3.xvg -contact
# gmx_mpi hbond -noda -f cluster1.xtc -s start10.tpr -num hbond_c1_R4.xvg -contact
# gmx_mpi hbond -noda -f cluster1.xtc -s start10.tpr -num hbond_c1_V5.xvg -contact


# gmx_mpi hbond -f cluster2.xtc -s start10.tpr -n hbond_R4_index.ndx -num hbond_c2_R4.xvg -contact -noda
