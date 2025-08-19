#ifdef DEBUG_HDF5
module debug_hdf5
  use hdf5
  implicit none
  integer(HID_T) :: dbg_file = 0
  integer :: dbg_level = 0
  integer :: dbg_rank  = 0
contains
  subroutine dbg_open(filename, level, rank_in)
    character(len=*), intent(in) :: filename
    integer,          intent(in) :: level, rank_in
    integer :: ierr
    dbg_level = level
    dbg_rank  = rank_in
    call h5open_f(ierr)
    call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, dbg_file, ierr)
  end subroutine dbg_open

  subroutine dbg_close()
    integer :: ierr
    if (dbg_file > 0) then
      call h5fclose_f(dbg_file, ierr)
      dbg_file = 0
    end if
    call h5close_f(ierr)
  end subroutine dbg_close

  subroutine dbg_write_scalar(group, name, value, step)
    character(len=*), intent(in) :: group, name
    real(8),          intent(in) :: value
    integer,          intent(in), optional :: step
    character(len=256) :: dset_name
    integer(HID_T)     :: grp, space, dset, attr_space, attr
    integer            :: ierr
    logical(kind=4)    :: exists
    integer(HSIZE_T)   :: dims(1), adims(1)

    if (dbg_file <= 0) return
    if (present(step)) then
      write(dset_name, '(a,"_step",i4.4)') trim(name), step
    else
      dset_name = trim(name)
    end if

    call h5gopen_f(dbg_file, trim(group), grp, ierr)
    if (ierr < 0) then
      call h5gcreate_f(dbg_file, trim(group), grp, ierr)
    end if

    dims(1) = int(1, kind=HSIZE_T)
    call h5screate_simple_f(1, dims, space, ierr)
    call h5lexists_f(grp, trim(dset_name), exists, ierr)
    if (.not. exists) then
      call h5dcreate_f(grp, trim(dset_name), H5T_NATIVE_DOUBLE, space, dset, ierr)
    else
      call h5dopen_f(grp, trim(dset_name), dset, ierr)
    end if
    call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, value, dims, ierr)
    call h5sclose_f(space, ierr)

    ! rank attribute
    adims(1) = int(1, kind=HSIZE_T)
    call h5screate_simple_f(1, adims, attr_space, ierr)
    call h5acreate_f(dset, "rank", H5T_NATIVE_INTEGER, attr_space, attr, ierr)
    call h5awrite_f(attr, H5T_NATIVE_INTEGER, dbg_rank, adims, ierr)
    call h5aclose_f(attr, ierr)
    call h5sclose_f(attr_space, ierr)

    call h5dclose_f(dset, ierr)
    call h5gclose_f(grp, ierr)
  end subroutine dbg_write_scalar

  subroutine dbg_write_array(group, name, arr, step)
    character(len=*), intent(in) :: group, name
    real(8),          intent(in) :: arr(:)
    integer,          intent(in), optional :: step
    character(len=256) :: dset_name
    integer(HID_T)     :: grp, space, dset, attr_space, attr
    integer            :: ierr
    logical(kind=4)    :: exists
    integer(HSIZE_T)   :: dims(1), adims(1)

    if (dbg_file <= 0) return
    if (present(step)) then
      write(dset_name, '(a,"_step",i4.4)') trim(name), step
    else
      dset_name = trim(name)
    end if

    call h5gopen_f(dbg_file, trim(group), grp, ierr)
    if (ierr < 0) then
      call h5gcreate_f(dbg_file, trim(group), grp, ierr)
    end if

    dims(1) = int(size(arr), HSIZE_T)
    call h5screate_simple_f(1, dims, space, ierr)
    call h5lexists_f(grp, trim(dset_name), exists, ierr)
    if (.not. exists) then
      call h5dcreate_f(grp, trim(dset_name), H5T_NATIVE_DOUBLE, space, dset, ierr)
    else
      call h5dopen_f(grp, trim(dset_name), dset, ierr)
    end if
    call h5dwrite_f(dset, H5T_NATIVE_DOUBLE, arr, dims, ierr)
    call h5sclose_f(space, ierr)

    adims(1) = int(1, kind=HSIZE_T)
    call h5screate_simple_f(1, adims, attr_space, ierr)
    call h5acreate_f(dset, "rank", H5T_NATIVE_INTEGER, attr_space, attr, ierr)
    call h5awrite_f(attr, H5T_NATIVE_INTEGER, dbg_rank, adims, ierr)
    call h5aclose_f(attr, ierr)
    call h5sclose_f(attr_space, ierr)

    call h5dclose_f(dset, ierr)
    call h5gclose_f(grp, ierr)
  end subroutine dbg_write_array

  subroutine dbg_write_int_scalar(group, name, value, step)
    character(len=*), intent(in) :: group, name
    integer,          intent(in) :: value
    integer,          intent(in), optional :: step
    character(len=256) :: dset_name
    integer(HID_T)     :: grp, space, dset, attr_space, attr
    integer            :: ierr
    logical(kind=4)    :: exists
    integer(HSIZE_T)   :: dims(1), adims(1)

    if (dbg_file <= 0) return
    if (present(step)) then
      write(dset_name, '(a,"_step",i4.4)') trim(name), step
    else
      dset_name = trim(name)
    end if

    call h5gopen_f(dbg_file, trim(group), grp, ierr)
    if (ierr < 0) then
      call h5gcreate_f(dbg_file, trim(group), grp, ierr)
    end if

    dims(1) = int(1, kind=HSIZE_T)
    call h5screate_simple_f(1, dims, space, ierr)
    call h5lexists_f(grp, trim(dset_name), exists, ierr)
    if (.not. exists) then
      call h5dcreate_f(grp, trim(dset_name), H5T_NATIVE_INTEGER, space, dset, ierr)
    else
      call h5dopen_f(grp, trim(dset_name), dset, ierr)
    end if
    call h5dwrite_f(dset, H5T_NATIVE_INTEGER, value, dims, ierr)
    call h5sclose_f(space, ierr)

    adims(1) = int(1, kind=HSIZE_T)
    call h5screate_simple_f(1, adims, attr_space, ierr)
    call h5acreate_f(dset, "rank", H5T_NATIVE_INTEGER, attr_space, attr, ierr)
    call h5awrite_f(attr, H5T_NATIVE_INTEGER, dbg_rank, adims, ierr)
    call h5aclose_f(attr, ierr)
    call h5sclose_f(attr_space, ierr)

    call h5dclose_f(dset, ierr)
    call h5gclose_f(grp, ierr)
  end subroutine dbg_write_int_scalar

  subroutine dbg_write_int_array(group, name, arr, step)
    character(len=*), intent(in) :: group, name
    integer,          intent(in) :: arr(:)
    integer,          intent(in), optional :: step
    character(len=256) :: dset_name
    integer(HID_T)     :: grp, space, dset, attr_space, attr
    integer            :: ierr
    logical(kind=4)    :: exists
    integer(HSIZE_T)   :: dims(1), adims(1)

    if (dbg_file <= 0) return
    if (present(step)) then
      write(dset_name, '(a,"_step",i4.4)') trim(name), step
    else
      dset_name = trim(name)
    end if

    call h5gopen_f(dbg_file, trim(group), grp, ierr)
    if (ierr < 0) then
      call h5gcreate_f(dbg_file, trim(group), grp, ierr)
    end if

    dims(1) = int(size(arr), HSIZE_T)
    call h5screate_simple_f(1, dims, space, ierr)
    call h5lexists_f(grp, trim(dset_name), exists, ierr)
    if (.not. exists) then
      call h5dcreate_f(grp, trim(dset_name), H5T_NATIVE_INTEGER, space, dset, ierr)
    else
      call h5dopen_f(grp, trim(dset_name), dset, ierr)
    end if
    call h5dwrite_f(dset, H5T_NATIVE_INTEGER, arr, dims, ierr)
    call h5sclose_f(space, ierr)

    adims(1) = int(1, kind=HSIZE_T)
    call h5screate_simple_f(1, adims, attr_space, ierr)
    call h5acreate_f(dset, "rank", H5T_NATIVE_INTEGER, attr_space, attr, ierr)
    call h5awrite_f(attr, H5T_NATIVE_INTEGER, dbg_rank, adims, ierr)
    call h5aclose_f(attr, ierr)
    call h5sclose_f(attr_space, ierr)

    call h5dclose_f(dset, ierr)
    call h5gclose_f(grp, ierr)
  end subroutine dbg_write_int_array

  subroutine dbg_log_msg(group, msg, step)
    character(len=*), intent(in) :: group, msg
    integer,          intent(in), optional :: step
    character(len=256) :: dset_name
    integer(HID_T)     :: grp, space, dset, dtype, attr_space, attr
    integer            :: ierr
    logical(kind=4)    :: exists
    integer(HSIZE_T)   :: dims(1), adims(1)
    integer(SIZE_T)    :: len

    if (dbg_file <= 0) return
    if (present(step)) then
      write(dset_name, '(a,"_step",i4.4)') 'log', step
    else
      dset_name = 'log'
    end if

    call h5gopen_f(dbg_file, trim(group), grp, ierr)
    if (ierr < 0) then
      call h5gcreate_f(dbg_file, trim(group), grp, ierr)
    end if

    len = int(len_trim(msg), SIZE_T)
    dims(1) = int(1, kind=HSIZE_T)
    call h5screate_simple_f(1, dims, space, ierr)
    call h5tcopy_f(H5T_C_S1, dtype, ierr)
    call h5tset_size_f(dtype, len, ierr)
    call h5lexists_f(grp, trim(dset_name), exists, ierr)
    if (.not. exists) then
      call h5dcreate_f(grp, trim(dset_name), dtype, space, dset, ierr)
    else
      call h5dopen_f(grp, trim(dset_name), dset, ierr)
    end if
    call h5dwrite_f(dset, dtype, msg(1:int(len)), dims, ierr)

    ! rank attribute
    adims(1) = int(1, kind=HSIZE_T)
    call h5screate_simple_f(1, adims, attr_space, ierr)
    call h5acreate_f(dset, "rank", H5T_NATIVE_INTEGER, attr_space, attr, ierr)
    call h5awrite_f(attr, H5T_NATIVE_INTEGER, dbg_rank, adims, ierr)
    call h5aclose_f(attr, ierr)
    call h5sclose_f(attr_space, ierr)

    call h5tclose_f(dtype, ierr)
    call h5sclose_f(space, ierr)
    call h5dclose_f(dset, ierr)
    call h5gclose_f(grp, ierr)
  end subroutine dbg_log_msg

end module debug_hdf5

#else
module debug_hdf5
end module debug_hdf5
#endif
