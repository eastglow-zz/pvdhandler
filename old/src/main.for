      program main
      
      use mpi 
      use cla 
      use grid_module
      use boundary_conditions
      use pvoutputs
        
      implicit none 

      
      integer(4) :: i,j,k
      character(len=100) :: option_value
      logical :: opt_found
      real(8) :: timestepval
      integer(4) :: iter_num 
      type(grid) :: grid_main 
      real(8),allocatable :: myarray(:,:,:)
      integer(4) :: istatus(MPI_STATUS_SIZE), ierr
      integer(4) :: nprocs, myrank 
      character(len=4) :: myrankstr 

      CALL MPI_INIT(ierr)
      CALL MPI_COMM_SIZE(MPI_COMM_WORLD,nprocs,ierr)
      CALL MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ierr)

      write(myrankstr,'(I4.4)') myrank
      write(*,*)'myrankstr = ',myrankstr 

      call initialize_grid_from_file(grid_main, 
     &                  './exampledata/grid_main/grid_main_'//myrankstr)

      call grid_verbose(grid_main)

      allocate(myarray(1-1:grid_main%iml+1,1,1-1:grid_main%kml+1))


      call read_iternum_from_pvd(iter_num, './exampledata/','SI',2.8d0)

      write(*,*) 'inter_num from pvd @t=2.8 s',iter_num

      call mpi_barrier(MPI_COMM_WORLD,ierr)
      call read_vtr_r8_2d(grid_main, 1, myarray, './exampledata/','SI',
     &                                                         iter_num)
      call mpi_barrier(MPI_COMM_WORLD,ierr)

      ! do k=1,grid_main%kml 
      ! do j=1,1
      ! do i=1,grid_main%iml 
      !   write(*,*) myarray(i,j,k)
      ! enddo
      ! enddo
      ! enddo

      deallocate(myarray)

      call MPI_FINALIZE(ierr)

      end program main 