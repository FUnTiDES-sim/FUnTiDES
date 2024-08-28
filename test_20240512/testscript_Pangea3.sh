DIR=/appli_RD/JIEMENG/SEMproxy/proxy_20240430/proxys/install/bin

declare -a arr=( 
"sem_Kokkos_V100GPU" 
"sem_Raja_V100GPU" 
"fd_Kokkos_V100GPU" 
"fd_Raja_V100GPU" 
"sem_OMP_Power9" 
"sem_SEQUENTIAL_Power9"
"fd_OMP_Power9" 
"fd_SEQUENTIAL_Power9" 
)

export OMP_PROC_BIND=spread; export OMP_PLACES=threads
export OMP_NUM_THREADS=36

for app in "${arr[@]}"
do 
  exec=${DIR}/${app}.exe
  echo $exec
  output=result_Pangea3_${app}_01.txt
  $exec 2>&1 | tee ${output}; 
  output=result_Pangea3_${app}_02.txt
  $exec 2>&1 | tee ${output}; 
  output=result_Pangea3_${app}_03.txt
  $exec 2>&1 | tee ${output}; 
  #grep computeOneStep ${output} | awk '{print $9}' >>  ${output}
done
