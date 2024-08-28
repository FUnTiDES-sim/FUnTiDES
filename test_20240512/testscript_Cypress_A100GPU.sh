DIR=/shared/data1/Users/j0535952/work2024/proxys_2024/SEM_Proxy/proxys_main/proxys_main_base/install/bin

declare -a arr=( 
"sem_Kokkos_A100GPU" 
"sem_Raja_A100GPU" 
"fd_Kokkos_A100GPU" 
"fd_Raja_A100GPU" 
"sem_OMP_EPYCROME" 
"sem_SEQUENTIAL_EPYCROME"
"fd_OMP_EPYCROME" 
"fd_SEQUENTIAL_EPYCROME" 
)

export OMP_PROC_BIND=spread; export OMP_PLACES=threads
export OMP_NUM_THREADS=16

for app in "${arr[@]}"
do 
  exec=${DIR}/${app}.exe
  echo $exec
  output=result_Cypress_${app}_01.txt
  $exec 2>&1 | tee ${output}; 
  output=result_Cypress_${app}_02.txt
  $exec 2>&1 | tee ${output}; 
  output=result_Cypress_${app}_03.txt
  $exec 2>&1 | tee ${output}; 
  #grep computeOneStep ${output} | awk '{print $9}' >>  ${output}
done
