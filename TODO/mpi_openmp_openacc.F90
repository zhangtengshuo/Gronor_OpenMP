module data_module
    use iso_c_binding
    use omp_lib
    implicit none

    ! --- MPI Tags ---
    integer, parameter :: TAG_TASK_REQUEST    = 100 ! Worker -> Master
    integer, parameter :: TAG_TASK_DESC       = 101 ! Master -> Worker [start_row, start_col]
    integer, parameter :: TAG_NO_MORE_TASKS   = 102 ! Master -> Worker
    integer, parameter :: TAG_TASK_DONE       = 103 ! Worker -> Master [start_row, start_col]
    integer, parameter :: TAG_GLOBAL_DATA     = 104 ! Master -> Worker (for Bcast)

    ! --- 问题尺寸 ---
    integer, parameter :: GLOBAL_N = 16000 
    integer, parameter :: GLOBAL_M = 16000
    integer, parameter :: SUB_N    = 200
    integer, parameter :: SUB_M    = 200
    integer, parameter :: NUM_TASKS = ((GLOBAL_N - 1) / SUB_N + 1) * ((GLOBAL_M - 1) / SUB_M + 1) ! ~100 tasks

    ! --- 模块变量 ---
    real(8), pointer, dimension(:,:), public :: global_matrix => null()
    logical, public :: data_on_host = .false.

contains

    ! 在主机上分配全局缓冲区
    subroutine allocate_global_buffer()
        if (associated(global_matrix)) then
            deallocate(global_matrix)
        end if
        allocate(global_matrix(GLOBAL_N, GLOBAL_M))
        data_on_host = .true.
    end subroutine allocate_global_buffer

    ! 释放全局缓冲区
    subroutine deallocate_global_buffer()
        if (associated(global_matrix)) then
            deallocate(global_matrix)
            global_matrix => null()
            data_on_host = .false.
        end if
    end subroutine deallocate_global_buffer

    ! Master Rank 用来初始化全局大矩阵
    subroutine initialize_global_matrix()
        integer :: i, j
        real(8) :: random_val
        integer, allocatable :: seed(:)
        integer :: seed_size

        if (.not. associated(global_matrix)) then
            call allocate_global_buffer()
        end if

        ! 设置随机种子
        call random_seed(size=seed_size)
        allocate(seed(seed_size))
        do i = 1, min(seed_size, 8)
            seed(i) = 123456789 + i * 1000
        end do
        if (seed_size > 8) then
            seed(9:) = 1
        end if
        call random_seed(put=seed)
        deallocate(seed)

        ! 初始化随机数据
        write(*,*) 'Master: Initializing global matrix (', GLOBAL_N, 'x', GLOBAL_M, ')...'
        !$omp parallel do private(i, j, random_val) collapse(2)
        do i = 1, GLOBAL_N
            do j = 1, GLOBAL_M
                call random_number(random_val)
                global_matrix(i, j) = random_val * 100.0_8
            end do
        end do
        !$omp end parallel do
        write(*,*) 'Master: Global matrix initialization completed.'
    end subroutine initialize_global_matrix

end module data_module

! ======================================================================

program mpi_openmp_openacc_svd_hybrid_main
    use data_module
    use mpi
    use omp_lib
    use openacc
    implicit none

    ! --- 变量声明部分 (必须在所有可执行语句之前) ---
    integer :: rank, size, ierr
    integer :: master_rank
    integer :: gpu_id
    integer :: provided_thread_level
    integer :: tasks_completed = 0
    integer :: current_task_index = 0
    integer :: total_tasks
    integer :: task_desc(2) ! [start_row, start_col]
    integer :: status(MPI_STATUS_SIZE)
    logical :: all_done = .false.
    logical :: flag
    real(8) :: start_time, end_time

    ! --- 1. 初始化MPI (线程安全) ---
    call MPI_Init_thread(MPI_THREAD_MULTIPLE, provided_thread_level, ierr)
    if (ierr /= MPI_SUCCESS) then
        write(*,*) 'MPI_Init_thread failed with error: ', ierr
        stop
    end if
    if (provided_thread_level /= MPI_THREAD_MULTIPLE) then
        write(*,*) 'Warning: MPI_THREAD_MULTIPLE not provided. Got level: ', provided_thread_level
    end if

    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, size, ierr)

    if (size /= 9) then
        if (rank == size - 1) then
            write(*,'(A,I0,A)') 'Error: This program requires exactly 9 MPI ranks, but ', size, ' were provided.'
        end if
        call MPI_Finalize(ierr)
        stop
    end if

    ! --- 2. 设置 Rank/GPU 对应关系 ---
    master_rank = size - 1
    if (rank < master_rank) then
        gpu_id = rank
        call acc_set_device_num(gpu_id, acc_device_nvidia)
        write(*,'(A,I0,A,I0,A)') 'Rank ', rank, ' assigned to GPU ', gpu_id
    else
        gpu_id = 0 ! Master uses GPU 0 if needed
        write(*,'(A,I0,A)') 'Rank ', rank, ' (Master) initialized'
    end if

    ! --- 3. 同步所有进程 ---
    call MPI_Barrier(MPI_COMM_WORLD, ierr)

    ! --- 4. Master rank 生成数据并广播 ---
    if (rank == master_rank) then
        call initialize_global_matrix()
    end if

    ! 所有 Rank (包括 Master) 都需要主机缓冲区来接收/持有数据
    if (.not. associated(global_matrix)) then
        call allocate_global_buffer()
    end if

    ! 广播完整数据
    write(*,'(A,I0,A)') 'Rank ', rank, ' receiving global matrix data...'
    call MPI_Bcast(global_matrix, GLOBAL_N * GLOBAL_M, MPI_DOUBLE_PRECISION, master_rank, MPI_COMM_WORLD, ierr)
    if (ierr /= MPI_SUCCESS) then
        write(*,'(A,I0,A,I0)') 'Rank ', rank, ': MPI_Bcast failed, ierr=', ierr
    else
        write(*,'(A,I0,A)') 'Rank ', rank, ' received global matrix data.'
    end if
    data_on_host = .true.

    ! --- 5. Worker ranks 上传数据到 GPU 并执行任务 ---
    if (rank < master_rank) then
        write(*,'(A,I0,A)') 'Rank ', rank, ' allocating device memory and uploading global matrix...'

        ! --- 关键步骤: 主线程将完整数据上传到 GPU ---
        !$acc data copyin(global_matrix(1:GLOBAL_N, 1:GLOBAL_M))
        ! 数据现在在设备上，整个 !$acc data ... !$acc end data 区域内都可用

        write(*,'(A,I0,A)') 'Rank ', rank, ' global matrix uploaded to GPU. Starting worker loop...'

        ! 设置OpenMP线程数
        call omp_set_num_threads(4) ! 例如4个线程

        ! 启动工作循环 (在 OpenMP 主线程中)
        call worker_main_loop(rank, master_rank)

        ! --- 关键步骤结束: 当此区域结束时，设备上的 global_matrix 副本才会被释放 ---
        !$acc end data

        write(*,'(A,I0,A)') 'Rank ', rank, ' finished all tasks.'
    end if

    ! --- 6. Master rank 管理任务 ---
    if (rank == master_rank) then
        total_tasks = NUM_TASKS
        write(*,'(A,I0,A)') 'Master: Total tasks to dispatch: ', total_tasks
        start_time = MPI_Wtime()

        do while (.not. all_done)
            ! 等待任何 worker 的任务请求
            call MPI_Recv(task_desc, 2, MPI_INTEGER, MPI_ANY_SOURCE, TAG_TASK_REQUEST, MPI_COMM_WORLD, status, ierr)
            if (ierr /= MPI_SUCCESS) then
                write(*,*) 'Master: MPI_Recv TASK_REQUEST failed, ierr=', ierr
                exit
            end if

            if (current_task_index < total_tasks) then
                current_task_index = current_task_index + 1
                ! 计算任务区域 (简化: 按行优先)
                task_desc(1) = mod((current_task_index - 1), (GLOBAL_N - 1) / SUB_N + 1) * SUB_N + 1
                task_desc(2) = ((current_task_index - 1) / ((GLOBAL_N - 1) / SUB_N + 1)) * SUB_M + 1
                task_desc(1) = min(task_desc(1), GLOBAL_N) ! 边界保护
                task_desc(2) = min(task_desc(2), GLOBAL_M)

                call MPI_Send(task_desc, 2, MPI_INTEGER, status(MPI_SOURCE), TAG_TASK_DESC, MPI_COMM_WORLD, ierr)
            else
                ! 没有更多任务，发送结束信号
                task_desc = -1 ! 标记为无任务
                call MPI_Send(task_desc, 2, MPI_INTEGER, status(MPI_SOURCE), TAG_NO_MORE_TASKS, MPI_COMM_WORLD, ierr)
            end if

            ! 检查任务完成消息
            do
                call MPI_Iprobe(MPI_ANY_SOURCE, TAG_TASK_DONE, MPI_COMM_WORLD, flag, status, ierr)
                if (ierr == MPI_SUCCESS .and. flag) then
                    call MPI_Recv(task_desc, 2, MPI_INTEGER, status(MPI_SOURCE), TAG_TASK_DONE, MPI_COMM_WORLD, status, ierr)
                    if (ierr == MPI_SUCCESS) then
                        tasks_completed = tasks_completed + 1
                        ! write(*,'(A,I0,A,2I6)') 'Master: Task completed by Rank ', status(MPI_SOURCE), ' [', task_desc(1:2), ']'
                        if (tasks_completed >= total_tasks) then
                            all_done = .true.
                            end_time = MPI_Wtime()
                            write(*,'(A,I0,A,F10.2,A)') 'Master: All ', tasks_completed, ' tasks completed in ', end_time - start_time, ' seconds.'
                            exit
                        end if
                    end if
                else
                    exit
                end if
            end do
        end do

        ! 清理
        call deallocate_global_buffer()
        write(*,'(A)') 'Master rank: Program execution completed successfully.'
    end if

    ! --- 7. 同步并结束 ---
    call MPI_Barrier(MPI_COMM_WORLD, ierr)
    call MPI_Finalize(ierr)

end program mpi_openmp_openacc_svd_hybrid_main

! ======================================================================
!                Worker 实现部分
! ======================================================================

! Worker 的主线程执行工作循环
subroutine worker_main_loop(my_rank, master_rank)
    use data_module
    use mpi
    use omp_lib
    use openacc
    implicit none

    ! --- 变量声明部分 ---
    integer, intent(in) :: my_rank, master_rank
    integer :: thread_id
    integer :: omp_rank
    integer :: ierr ! 添加缺失的 ierr 声明

    ! 在 OpenMP 并行区域外启动并行区域
    !$omp parallel private(thread_id, omp_rank, ierr) ! 将 ierr 添加到 private 列表
    thread_id = omp_get_thread_num()
    ! 在 OpenMP 线程中安全地获取 MPI rank
    call MPI_Comm_rank(MPI_COMM_WORLD, omp_rank, ierr)

    ! 每个 OpenMP 线程执行自己的工作循环
    call worker_thread_loop(omp_rank, thread_id, master_rank)

    !$omp end parallel

end subroutine worker_main_loop

! 每个 OpenMP 线程的工作循环
subroutine worker_thread_loop(my_rank, thread_id, master_rank)
    use data_module
    use mpi
    use openacc
    implicit none

    ! --- 变量声明部分 ---
    integer, intent(in) :: my_rank, thread_id, master_rank
    integer :: task_desc(2) ! [start_row, start_col]
    integer :: status(MPI_STATUS_SIZE)
    integer :: ierr ! 添加缺失的 ierr 声明

    do
        ! 1. 请求任务
        ! write(*,'(A,I0,A,I0,A)') 'Rank ', my_rank, ' Thread ', thread_id, ': Requesting task...'
        call MPI_Send(task_desc, 2, MPI_INTEGER, master_rank, TAG_TASK_REQUEST, MPI_COMM_WORLD, ierr)

        ! 2. 接收任务
        call MPI_Recv(task_desc, 2, MPI_INTEGER, master_rank, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)

        if (status(MPI_TAG) == TAG_NO_MORE_TASKS) then
            ! write(*,'(A,I0,A,I0,A)') 'Rank ', my_rank, ' Thread ', thread_id, ': No more tasks. Exiting.'
            exit
        else if (status(MPI_TAG) == TAG_TASK_DESC) then
            ! 3. 执行任务
            ! write(*,'(A,I0,A,I0,A,2I6)') 'Rank ', my_rank, ' Thread ', thread_id, ': Received task [', task_desc(1:2), ']'

            call execute_single_task(my_rank, thread_id, task_desc(1), task_desc(2))

            ! 4. 报告任务完成
            ! write(*,'(A,I0,A,I0,A,2I6,A)') 'Rank ', my_rank, ' Thread ', thread_id, ': Completed task [', task_desc(1:2), ']. Reporting...'
            call MPI_Send(task_desc, 2, MPI_INTEGER, master_rank, TAG_TASK_DONE, MPI_COMM_WORLD, ierr)
        end if
    end do

end subroutine worker_thread_loop

! 执行单个任务：CPU SVD -> GPU OpenACC 计算
! 注意：此子程序在 OpenMP 线程内调用，可以访问外层 !$acc data region 中的 global_matrix
subroutine execute_single_task(my_rank, thread_id, start_row, start_col)
    use data_module
    use mpi
    use openacc
    implicit none

    ! --- 变量声明部分 (所有声明必须在可执行语句之前) ---
    integer, intent(in) :: my_rank, thread_id, start_row, start_col
    integer, parameter :: local_n = SUB_N, local_m = SUB_M
    real(8), dimension(local_n, local_m) :: local_sub_matrix_svd
    real(8), dimension(local_n, local_n) :: U ! U matrix from SVD
    real(8), dimension(min(local_n, local_m)) :: S ! Singular values
    real(8), dimension(local_m, local_m) :: VT ! V transpose
    real(8), dimension(local_n, local_n) :: gpu_result ! GPU计算结果
    integer :: ierr, i, j, k
    integer :: actual_rows, actual_cols
    integer :: lwork
    real(8), allocatable :: work(:)
    real(8), parameter :: zero_threshold = 1.0e-10
    integer :: info
    real(8) :: svd_norm, gpu_norm
    real(8) :: combined_result ! 将 combined_result 的声明移到开头

    ! --- 1. 确定实际区域大小 ---
    actual_rows = min(local_n, GLOBAL_N - start_row + 1)
    actual_cols = min(local_m, GLOBAL_M - start_col + 1)

    ! --- 2. 在 CPU 上执行 SVD ---
    ! 从 Worker Rank 持有的完整 global_matrix 中拷贝子矩阵
    local_sub_matrix_svd(1:actual_rows, 1:actual_cols) = &
        global_matrix(start_row : start_row + actual_rows - 1, &
                      start_col : start_col + actual_cols - 1)

    lwork = -1
    allocate(work(1))
    call dgesvd('A', 'N', actual_rows, actual_cols, local_sub_matrix_svd, local_n, S, U, local_n, VT, local_m, work, lwork, info)
    if (info == 0) then
        lwork = int(work(1))
        deallocate(work)
        allocate(work(lwork))
        call dgesvd('A', 'N', actual_rows, actual_cols, local_sub_matrix_svd, local_n, S, U, local_n, VT, local_m, work, lwork, info)
        deallocate(work)
    else
        write(*,*) 'Rank ', my_rank, ' Thread ', thread_id, ': SVD query failed, info=', info
        deallocate(work)
        return
    end if

    if (info /= 0) then
        write(*,*) 'Rank ', my_rank, ' Thread ', thread_id, ': SVD failed, info=', info
        return
    end if
    ! 计算 SVD 结果的范数作为示例
    svd_norm = 0.0_8
    do i = 1, min(actual_rows, actual_cols)
        svd_norm = svd_norm + S(i) * S(i)
    end do
    svd_norm = sqrt(svd_norm)
    ! write(*,'(A,I0,A,I0,A,2I6,A,F10.4)') 'Rank ', my_rank, ' Thread ', thread_id, ': SVD norm for task [', start_row, start_col, '] is ', svd_norm

    ! --- 3. 在 GPU 上使用 OpenACC 执行计算 (使用外层已上传的 global_matrix) ---
    ! 示例：计算 global_matrix 区域的元素平方和
    gpu_norm = 0.0_8 ! 初始化归约变量
!    !$acc parallel loop gang vector collapse(2) reduction(+:gpu_norm) private(i,j)
!    do j = 1, actual_cols ! 注意OpenACC常用列优先
!        do i = 1, actual_rows
!            ! global_matrix 在外层 !$acc data copyin 中，此处可直接访问其设备副本
!            gpu_norm = gpu_norm + global_matrix(start_row + i - 1, start_col + j - 1) * &
!                                 global_matrix(start_row + i - 1, start_col + j - 1)
!        end do
!    end do
!    !$acc end parallel loop

    !$acc parallel loop gang vector collapse(2) reduction(+:gpu_norm) private(i,j,k)
    do j = 1, actual_cols
        do i = 1, actual_rows
            !$acc loop seq
            do k = 1, 100000
                gpu_norm = gpu_norm + global_matrix(start_row + i - 1, start_col + j - 1) * &
                                        global_matrix(start_row + i - 1, start_col + j - 1)
            end do
        end do
    end do
    !$acc end parallel loop

    gpu_norm = sqrt(gpu_norm)
    ! write(*,'(A,I0,A,I0,A,2I6,A,F10.4)') 'Rank ', my_rank, ' Thread ', thread_id, ': GPU norm for task [', start_row, start_col, '] is ', gpu_norm

    ! --- 4. (可选) 结合 SVD 和 GPU 计算结果 ---
    ! 例如，将两个结果相加或比较
    combined_result = svd_norm + gpu_norm
    ! write(*,'(A,I0,A,I0,A,2I6,A,F10.4)') 'Rank ', my_rank, ' Thread ', thread_id, ': Combined result for task [', start_row, start_col, '] is ', combined_result

end subroutine execute_single_task