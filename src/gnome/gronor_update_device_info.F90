subroutine gronor_update_device_info()

  use cidist

#ifdef _OPENACC
  use openacc
#ifdef CUDA
  use cuda_functions
#endif
#endif
  
#ifdef CUDA
  type(c_ptr) :: cpfre, cptot
#endif
  
#ifdef CUDA
  cpfre=c_loc(memfre)
  cptot=c_loc(memtot)
  istat=cudaMemGetInfo(cpfre,cptot)
#endif

  memavail=memfre

  return
  
end subroutine gronor_update_device_info
