module io_hdf5
    use hdf5
    use parameters, only: dp
    implicit none
    private

    public :: h5_init_run, h5_write_snapshot, h5_close_run

    integer(HID_T) :: file_id = -1

    integer(HID_T) :: dset_u = -1, dset_v = -1, dset_t = -1
    integer(HID_T) :: sid_u  = -1, sid_v  = -1, sid_t  = -1     ! dataspace ids (NO usar space_t)
    integer(HID_T) :: dcpl_u = -1, dcpl_v = -1, dcpl_t = -1     ! dataset create property lists

    integer :: nxp = 0, my = 0
    integer :: nsnap_written = 0

contains

    subroutine h5_init_run(filename, Nx_phys, M, hx, hy, dt, vl, Nfull, Mfull)
        character(len=*), intent(in) :: filename
        integer, intent(in) :: Nx_phys, M, Nfull, Mfull
        real(dp), intent(in) :: hx, hy, dt, vl

        integer :: hdferr
        integer(HID_T) :: aspace, attr_id
        integer(HSIZE_T), dimension(3) :: dims3, maxdims3, chunk3
        integer(HSIZE_T), dimension(1) :: dims1, maxdims1, chunk1, adims

        nxp = Nx_phys
        my  = M
        nsnap_written = 0

        call h5open_f(hdferr)

        call h5fcreate_f(trim(filename), H5F_ACC_TRUNC_F, file_id, hdferr)
        if (hdferr /= 0) stop "HDF5: error creating file"

        ! -------------------------
        ! Dataset u, v: (Nx_phys, M, nsnap) con nsnap extendible
        ! -------------------------
        dims3    = (/ int(nxp,HSIZE_T), int(my,HSIZE_T), 0_HSIZE_T /)  ! nsnap=0 al inicio
        maxdims3 = (/ int(nxp,HSIZE_T), int(my,HSIZE_T), H5S_UNLIMITED_F /)
        chunk3   = (/ int(nxp,HSIZE_T), int(my,HSIZE_T), 1_HSIZE_T /)

        call h5screate_simple_f(3, dims3, sid_u, hdferr, maxdims3)
        call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_u, hdferr)
        call h5pset_chunk_f(dcpl_u, 3, chunk3, hdferr)
        call h5dcreate_f(file_id, "u", H5T_NATIVE_DOUBLE, sid_u, dset_u, hdferr, dcpl_id=dcpl_u)

        call h5screate_simple_f(3, dims3, sid_v, hdferr, maxdims3)
        call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_v, hdferr)
        call h5pset_chunk_f(dcpl_v, 3, chunk3, hdferr)
        call h5dcreate_f(file_id, "v", H5T_NATIVE_DOUBLE, sid_v, dset_v, hdferr, dcpl_id=dcpl_v)

        ! -------------------------
        ! Dataset t: (nsnap) extendible
        ! -------------------------
        dims1    = (/ 0_HSIZE_T /)
        maxdims1 = (/ H5S_UNLIMITED_F /)
        chunk1   = (/ 256_HSIZE_T /)

        call h5screate_simple_f(1, dims1, sid_t, hdferr, maxdims1)
        call h5pcreate_f(H5P_DATASET_CREATE_F, dcpl_t, hdferr)
        call h5pset_chunk_f(dcpl_t, 1, chunk1, hdferr)
        call h5dcreate_f(file_id, "t", H5T_NATIVE_DOUBLE, sid_t, dset_t, hdferr, dcpl_id=dcpl_t)

        ! -------------------------
        ! Atributos (params)
        ! -------------------------
        adims = (/ 1_HSIZE_T /)
        call h5screate_simple_f(1, adims, aspace, hdferr)

        call h5acreate_f(file_id, "hx", H5T_NATIVE_DOUBLE, aspace, attr_id, hdferr)
        call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, hx, adims, hdferr)
        call h5aclose_f(attr_id, hdferr)

        call h5acreate_f(file_id, "hy", H5T_NATIVE_DOUBLE, aspace, attr_id, hdferr)
        call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, hy, adims, hdferr)
        call h5aclose_f(attr_id, hdferr)

        call h5acreate_f(file_id, "dt", H5T_NATIVE_DOUBLE, aspace, attr_id, hdferr)
        call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, dt, adims, hdferr)
        call h5aclose_f(attr_id, hdferr)

        call h5acreate_f(file_id, "vl", H5T_NATIVE_DOUBLE, aspace, attr_id, hdferr)
        call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, vl, adims, hdferr)
        call h5aclose_f(attr_id, hdferr)

        call h5sclose_f(aspace, hdferr)

        call write_attr_int(file_id, "N", Nfull)
        call write_attr_int(file_id, "M", Mfull)
        call write_attr_int(file_id, "Nx_phys", Nx_phys)
        call write_attr_int(file_id, "M_phys", M)

    end subroutine h5_init_run


    subroutine h5_write_snapshot(t_phys, u2d, v2d)
        real(dp), intent(in) :: t_phys
        real(dp), intent(in), contiguous :: u2d(:,:), v2d(:,:)

        integer :: hdferr
        integer(HSIZE_T), dimension(3) :: newdims3, start3, count3
        integer(HSIZE_T), dimension(1) :: newdims1, start1, count1
        integer(HID_T) :: mem3, mem1, fspace

        ! Buffer 3D seguro (evita mismatch de rango al escribir hyperslab 3D)
        real(dp) :: u3(nxp, my, 1), v3(nxp, my, 1)

        if (file_id < 0) stop "HDF5: file not initialized"

        nsnap_written = nsnap_written + 1

        u3(:,:,1) = u2d(:,:)
        v3(:,:,1) = v2d(:,:)

        ! ---- Extender u, v a (nxp,my,nsnap)
        newdims3 = (/ int(nxp,HSIZE_T), int(my,HSIZE_T), int(nsnap_written,HSIZE_T) /)
        call h5dset_extent_f(dset_u, newdims3, hdferr)
        call h5dset_extent_f(dset_v, newdims3, hdferr)

        start3 = (/ 0_HSIZE_T, 0_HSIZE_T, int(nsnap_written-1,HSIZE_T) /)
        count3 = (/ int(nxp,HSIZE_T), int(my,HSIZE_T), 1_HSIZE_T /)

        call h5screate_simple_f(3, count3, mem3, hdferr)

        call h5dget_space_f(dset_u, fspace, hdferr)
        call h5sselect_hyperslab_f(fspace, H5S_SELECT_SET_F, start3, count3, hdferr)
        call h5dwrite_f(dset_u, H5T_NATIVE_DOUBLE, u3, count3, hdferr, file_space_id=fspace, mem_space_id=mem3)
        call h5sclose_f(fspace, hdferr)

        call h5dget_space_f(dset_v, fspace, hdferr)
        call h5sselect_hyperslab_f(fspace, H5S_SELECT_SET_F, start3, count3, hdferr)
        call h5dwrite_f(dset_v, H5T_NATIVE_DOUBLE, v3, count3, hdferr, file_space_id=fspace, mem_space_id=mem3)
        call h5sclose_f(fspace, hdferr)

        call h5sclose_f(mem3, hdferr)

        ! ---- Extender t a (nsnap)
        newdims1 = (/ int(nsnap_written,HSIZE_T) /)
        call h5dset_extent_f(dset_t, newdims1, hdferr)

        start1 = (/ int(nsnap_written-1,HSIZE_T) /)
        count1 = (/ 1_HSIZE_T /)

        call h5screate_simple_f(1, count1, mem1, hdferr)
        call h5dget_space_f(dset_t, fspace, hdferr)
        call h5sselect_hyperslab_f(fspace, H5S_SELECT_SET_F, start1, count1, hdferr)
        call h5dwrite_f(dset_t, H5T_NATIVE_DOUBLE, t_phys, count1, hdferr, file_space_id=fspace, mem_space_id=mem1)
        call h5sclose_f(fspace, hdferr)
        call h5sclose_f(mem1, hdferr)

    end subroutine h5_write_snapshot


    subroutine h5_close_run()
        integer :: hdferr

        if (dset_u >= 0) call h5dclose_f(dset_u, hdferr)
        if (dset_v >= 0) call h5dclose_f(dset_v, hdferr)
        if (dset_t >= 0) call h5dclose_f(dset_t, hdferr)

        if (sid_u  >= 0) call h5sclose_f(sid_u, hdferr)
        if (sid_v  >= 0) call h5sclose_f(sid_v, hdferr)
        if (sid_t  >= 0) call h5sclose_f(sid_t, hdferr)

        if (dcpl_u >= 0) call h5pclose_f(dcpl_u, hdferr)
        if (dcpl_v >= 0) call h5pclose_f(dcpl_v, hdferr)
        if (dcpl_t >= 0) call h5pclose_f(dcpl_t, hdferr)

        if (file_id >= 0) call h5fclose_f(file_id, hdferr)

        call h5close_f(hdferr)

        file_id = -1
        dset_u = -1; dset_v = -1; dset_t = -1
        sid_u  = -1; sid_v  = -1; sid_t  = -1
        dcpl_u = -1; dcpl_v = -1; dcpl_t = -1
        nsnap_written = 0
        nxp = 0; my = 0

    end subroutine h5_close_run


    subroutine write_attr_int(obj_id, name, val)
        integer(HID_T), intent(in) :: obj_id
        character(len=*), intent(in) :: name
        integer, intent(in) :: val

        integer :: hdferr
        integer(HID_T) :: aspace, attr_id
        integer(HSIZE_T), dimension(1) :: adims

        adims = (/ 1_HSIZE_T /)
        call h5screate_simple_f(1, adims, aspace, hdferr)

        call h5acreate_f(obj_id, trim(name), H5T_NATIVE_INTEGER, aspace, attr_id, hdferr)
        call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, val, adims, hdferr)
        call h5aclose_f(attr_id, hdferr)

        call h5sclose_f(aspace, hdferr)
    end subroutine write_attr_int

end module io_hdf5
