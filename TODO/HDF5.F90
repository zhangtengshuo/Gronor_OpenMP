program test_hdf5_types
  use hdf5
  implicit none
  print *, "HID_T kind:", HID_T
  print *, "HSIZE_T kind:", HSIZE_T
  print *, "SIZE_T kind:", SIZE_T
end program
