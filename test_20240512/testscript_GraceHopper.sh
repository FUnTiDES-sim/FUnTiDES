DIR=/shared/data1/Users/j0535952/work2024/proxys_2024/SEM_Proxy/proxys_main/proxys_main_base/install/bin

declare -a arr=( 
"sem_Kokkos_H100GPU" 
"sem_Raja_H100GPU" 
"fd_Kokkos_H100GPU" 
"fd_Raja_H100GPU" 
"sem_OMP_GRACE" 
"sem_SEQUENTIAL_GRACE"
"fd_OMP_GRACE" 
"fd_SEQUENTIAL_GRACE" 
)

export OMP_PROC_BIND=spread; export OMP_PLACES=threads
export OMP_NUM_THREADS=72

for app in "${arr[@]}"
do 
  exec=${DIR}/${app}.exe
  echo $exec
  output=result_GraceHopper_${app}_01.txt
  $exec 2>&1 | tee ${output}; 
  output=result_GraceHopper_${app}_02.txt
  $exec 2>&1 | tee ${output}; 
  output=result_GraceHopper_${app}_03.txt
  $exec 2>&1 | tee ${output}; 
  #grep computeOneStep ${output} | awk '{print $9}' >>  ${output}
done
