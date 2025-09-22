 
      program main_mini

      use grid_module

      implicit none 

      integer(4) :: i,j,k
      character(len=100) :: option_value
      logical :: opt_found
      real(8) :: timestepval
      integer(4) :: iter_num 
      character(300) :: path 
      type(grid) :: grid_main 
      real(8),allocatable :: myarray(:,:,:,:)
      integer(4) :: nprocs, pid 
      character(len=4) :: pid2s 

      path='./exampledata/fesi_eqx_tr_v_dsfar/outputs/'

      pid = 0 
      write(pid2s,'(i4.4)') pid
      write(*,*) trim(path)//'grid_main/grid_main_'//pid2s
      call initialize_grid_from_file(grid_main, 
     &                       trim(path)//'grid_main/grid_main_'//pid2s)

      call grid_verbose(grid_main)

      allocate(myarray(1-1:grid_main%iml+1,1,1-1:grid_main%kml+1,
     &                                                      0:nprocs-1))


      call read_iternum_from_pvd(iter_num, trim(path),'SI',2.8d0)

      write(*,*) 'inter_num from pvd @t=2.8 s',iter_num

      call read_vtr_r8_2d(grid_main, 1, myarray(:,:,:,pid), trim(path),
     &                                                    'SI',iter_num)

      do k=1,grid_main%kml 
      do j=1,1
      do i=1,grid_main%iml 
        write(*,*) myarray(i,j,k,0)
      enddo
      enddo
      enddo

      contains 

      ! ----------------------------------------------------------------

      subroutine get_option_value(option_string,key,option_value,found)
      implicit none
  
      ! Input arguments
      character(len=*), intent(in) :: option_string
      character(len=*), intent(in) :: key
  
      ! Output arguments
      character(len=*), intent(out) :: option_value
      logical, intent(out) :: found
  
      ! Local variables
      integer :: pos_key, pos_equals, pos_quote1, pos_quote2, len_key, 
     &                                                 len_option_string
  
      len_option_string = len_trim(option_string)
      len_key = len_trim(key)
      option_value = ''  ! Initialize the value to an empty string
      found = .false.  ! Initialize found to false
  
      ! Search for the key in the option string
      pos_key = index(option_string, trim(key) // '=')
  
      if (pos_key /= 0) then
        pos_equals = pos_key + len_key  ! Position after the '='
        ! Find the positions of the quotes around the value
        pos_quote1 = index(option_string(pos_equals+1:),'"')+pos_equals
        pos_quote2 = index(option_string(pos_quote1+1:),'"')+pos_quote1

        if (pos_quote1 > 0 .and. pos_quote2 > pos_quote1) then
            ! Extract the value between the quotes
            option_value = option_string(pos_quote1+1:pos_quote2-1)
            found = .true.
        end if
      end if

      return 
      end subroutine get_option_value


      ! ----------------------------------------------------------------
      subroutine extract_integer_from_pvtrpath(str, num)  
      implicit none
      character(len=*), intent(in) :: str
      integer, intent(out) :: num
      integer :: start_pos, end_pos
      character(len=100) :: number_str
  
      ! Find the position of the first underscore
      start_pos = index(str, '_')
  
      ! Find the position of the first period after the number
      end_pos = index(str, '.pvtr')
  
      ! Extract the substring that contains the number
      if (start_pos > 0 .and. end_pos > start_pos) then
          number_str = str(start_pos+1:end_pos-1)
          ! Convert the substring to an integer
          read(number_str, *) num
      else
          print *, "Error: Could not find the number in the string."
          num = -1
      end if

      return 
      end subroutine extract_integer_from_pvtrpath
      ! ----------------------------------------------------------------


      subroutine read_iternum_from_pvd(timestep,path,fnamebase,time)
      implicit none
      integer(4),parameter :: pvdheaderlines = 3
      integer(4),parameter :: pvdfooterlines = 2
      integer(4),intent(inout) :: timestep
      character(len=*),intent(in) :: path 
      character(len=*),intent(in) :: fnamebase
      real(8),intent(in) :: time

      character(len=100) :: fnamepvtr, fnamepvd
      character(len=8) :: ts2s
      character(len=20) :: time_str
      logical :: pvdexist
      integer(4) :: nlines, nbody, line
      character(len=300) :: linebuffer
      character(len=100) :: option_value
      logical :: opt_found
      real(8) :: timestepval 
      real(8) :: timestep_picked
      integer(4) :: iter_picked

      fnamepvd=trim(path)//trim(fnamebase)//'.pvd'

      write(*,*)trim(fnamepvd)

      ! Check if the file exists or not
      inquire(file=trim(fnamepvd),exist=pvdexist)
      if (pvdexist .eqv. .false.)then  ! .pvd doesn't exist. Give an error message and abort
         write(*,*)'Error, M801, read_iternum_from_pvd()'
         write(*,*)'  .pvd file cannot be found.'
         write(*,*)'  File name requested to read: '
         write(*,*)'  ', trim(fnamepvd)
         write(*,*)'  aborting...'
         call abort()
      else ! .pvd already exists. Update the lines.
        ! Find the # of lines in the pvd file
        call get_file_lines(nlines, trim(fnamepvd))
        nbody = nlines - pvdheaderlines - pvdfooterlines

        open(unit=900,file=trim(fnamepvd),status='old', action='read')

        ! Scan through the lines over the data
        iter_picked = -1
        do line = 1, pvdheaderlines + nbody
          read(900,'(A)') linebuffer 
          if(line.gt.pvdheaderlines) then 
            call get_option_value(trim(linebuffer),'timestep', 
     &                                           option_value,opt_found)
            if(opt_found.eqv..false.) then
              write(*,*)'timestep option not found.'
              write(*,*)'linebuffer: ',trim(linebuffer)
            endif
            read(option_value,*) timestepval

            if(timestepval.le.time) then 
              call get_option_value(trim(linebuffer),'file', 
     &                                           option_value,opt_found)
              fnamepvtr=option_value 
              timestep_picked = timestepval 
              call extract_integer_from_pvtrpath(fnamepvtr, iter_picked)
            endif

          endif !if(line.gt.pvdheaderlines) then 
        enddo ! do line = 1, pvdheaderlines + nbody
      
        close(900)

      endif ! if (pdvexist .eqv. .false.)then ... else

      if(iter_picked .eq. -1) then
        write(*,*)'Error, M801, read_iternum_from_pvd()'
        write(*,*)'  Matching data not found.'
        write(*,*)'  Input file value: ',time
        write(*,*)'  Aborting...'
        call abort()
      else 
        timestep = iter_picked
      endif

      return
      end subroutine read_iternum_from_pvd

      ! ----------------------------------------------------------------

      subroutine get_file_lines(nlines, fname)
      implicit none
      integer(4),intent(inout) :: nlines
      character(len=*),intent(in) :: fname

      integer(4):: io

      nlines = 0
      open(100, FILE=fname, status='old',action='read')
      DO
        READ(100,*,iostat=io)
        IF (io/=0) exit
        nlines = nlines + 1
      END DO
      close(100)
      return
      end subroutine get_file_lines


      ! ----------------------------------------------------------------

      subroutine parse_extents_6i(line, integers)
      implicit none
      character(len=*), intent(in) :: line
      integer, intent(out) :: integers(6)

      integer(4),parameter :: expected_width = 8
  
      integer :: istart, iend, ios
      character(len=100) :: temp_string
      character(len=8) :: padded_string

      character(len=expected_width*6) :: str_to_read=' '

      integer(4) :: len_input
      integer(4) :: missing_spaces

      integer(4) :: i
  
      ! Initialize the integers to zero
      integers = 0
  
      ! Find the start and end of the substring enclosed in quotes
      istart = index(line, '"') + 1
      iend = index(line(istart:), '"') + istart - 2
  
      ! Extract the substring containing the numbers
      temp_string = line(istart:iend)

      len_input = len(trim(temp_string))

      missing_spaces = expected_width*6 - len_input

      str_to_read(missing_spaces+1:len(str_to_read))=trim(temp_string)

      ! Read the numbers from the extracted substring
      read(str_to_read, '(6I8)', iostat=ios) integers
  
      ! Check for reading errors
      if (ios /= 0) then
        print *, "Error reading integers from the string:", temp_string
      end if
  
      return
      end subroutine parse_extents_6i
    
      ! ----------------------------------------------------------------

      subroutine get_extents_from_vtrheader(fnamevtr,wextent,pextent)
      implicit none 
      character(len=*),intent(in) :: fnamevtr 
      integer(4),intent(inout) :: wextent(6)
      integer(4),intent(inout) :: pextent(6)

      logical :: fexist
      character(len=300) :: wextentstr
      character(len=300) :: pextentstr

      inquire(file=trim(fnamevtr),exist=fexist)
      if (fexist .eqv. .false.) then
        write(*,*)'Error, M801, get_extents_from_vtrheader()'
        write(*,*)'  .vtr file cannot be found.'
        write(*,*)'  File name requested to read: '
        write(*,*)'  ', trim(fnamevtr)
        write(*,*)'  aborting...'
        call abort()
      else
        open(1001,FILE=trim(fnamevtr),status='old',action='read')
        read(1001,*) ! Skipping a line
        read(1001,*) ! Skipping another line 
        read(1001,'(A)') wextentstr
        read(1001,'(A)') pextentstr
        close(1001)

        call parse_extents_6i(wextentstr, wextent)
        call parse_extents_6i(pextentstr, pextent)

      endif  !if (fexist .eqv. .false.) then ... else 

      return 
      end subroutine get_extents_from_vtrheader

      ! ----------------------------------------------------------------

      subroutine read_vtr_r8_2d(gridobj, ng, a, path, fnamebase, iter)
      type(grid),intent(in) :: gridobj 
      integer(4),intent(in) :: ng
      real(8),intent(inout) :: 
     &                      a(1-ng:gridobj%iml+ng,1,1-ng:gridobj%kml+ng) 
      character(len=*),intent(in) :: path 
      character(len=*),intent(in) :: fnamebase 
      integer(4),intent(in) :: iter 

      character(len=300) :: fnamevtr 
      character(len=100) :: subdir 
      character(len=8) :: ts2s, pid2s
      logical :: fexist
      integer(4) :: wextent(6)
      integer(4) :: pextent(6)
      integer(4) :: numdat_x, numdat_y, numdat_z 
      integer(4) :: numgrid_x, numgrid_y, numgrid_z
      integer(4) :: l
      integer(4) :: i,j,k 
      integer(4) :: i0, j0, k0

      subdir = trim(fnamebase)//'/'

      write(ts2s,"(i8.8)") int(iter)
      write(pid2s,"(i4.4)") int(gridobj%pid)

      fnamevtr=trim(path)//trim(subdir)//trim(fnamebase)//
     & '_'//trim(ts2s)//'_'//trim(pid2s)//'.vtr'

      call get_extents_from_vtrheader(trim(fnamevtr), wextent, pextent)
      numdat_x = pextent(2)-pextent(1)+1
      numdat_y = pextent(4)-pextent(3)+1
      numdat_z = pextent(6)-pextent(5)+1

      numgrid_x = gridobj%iml 
      numgrid_y = gridobj%jml 
      numgrid_z = gridobj%kml
      if(gridobj%pidx.ge.1) numgrid_x = gridobj%iml + 1 
      if(gridobj%pidy.ge.1) numgrid_y = gridobj%jml + 1 
      if(gridobj%pidz.ge.1) numgrid_z = gridobj%kml + 1 

      ! Sanitation check
      if(numgrid_x.ne.numdat_x .or. numgrid_y.ne.numdat_y .or. 
     &                                      numgrid_z.ne.numdat_z ) then
        write(*,*)'Error, M801, read_vtr_r8_2d()'
        write(*,*)'  # of data and # of grid do not match.'
        write(*,*)'  pid: ',gridobj%pid 
        write(*,*)'  pidx, pidy, pidz: ',
     &                            gridobj%pidx,gridobj%pidy,gridobj%pidz
        write(*,*)'  numgrid_x, numgrid_y, numgrid_z: ', 
     &                                   numgrid_x, numgrid_y, numgrid_z 
        write(*,*)'  numdat_x, numdat_y, numdat_z: ',
     &                                      numdat_x, numdat_y, numdat_z
        write(*,*)'  Aborting...'
        call abort()
      endif

      inquire(file=trim(fnamevtr),exist=fexist)
      if (fexist .eqv. .false.) then
        write(*,*)'Error, M801, read_vtr_r8_2d()'
        write(*,*)'  .vtr file cannot be found.'
        write(*,*)'  File name requested to read: '
        write(*,*)'  ', trim(fnamevtr)
        write(*,*)'  Aborting...'
        call abort()
      else
        ! File exists. Start reading the vtr file. 
        open(1001,FILE=trim(fnamevtr),status='old',action='read')
        read(1001,*) ! Skipping <?xml version="1.0"?>
        read(1001,*) ! Skipping <VTKFile type="RectilinearGrid" version="0.1" byte_order="LittleEndian">
        read(1001,*) ! Skipping Whole Extent line which we already know
        read(1001,*) ! Skipping Piece Extent line which we already know 
        read(1001,*) ! Skipping <Coordinate> line
        read(1001,*) ! Skipping <DataArray Name="x" type="Float32" format="ascii">
        do l=1,numdat_x 
          read(1001,*) ! Now x coordinates comes (all skip)
        enddo
        read(1001,*) ! Skipping </DataArray>
        read(1001,*) ! Skipping <DataArray Name="y" type="Float32" format="ascii">
        do l=1,numdat_y 
          read(1001,*) ! Now y coordinates comes (all skip)
        enddo 
        read(1001,*) ! Skipping </DataArray>
        read(1001,*) ! Skipping <DataArray Name="z" type="Float32" format="ascii">
        do l=1,numdat_z 
          read(1001,*) ! Now z coordinates comes (all skip)
        enddo 
        read(1001,*) ! Skipping </DataArray>
        read(1001,*) ! Skipping </Coordinates>
        read(1001,*) ! Skipping <PointData Scalars="SI">
        read(1001,*) ! Skipping <DataArray type="Float32" Name="SI" format="ascii">
        ! Now actual data reading starts here 
        i0=1
        j0=1
        k0=1
        if(gridobj%pidx.ge.1) i0=0
        if(gridobj%pidy.ge.1) j0=0
        if(gridobj%pidz.ge.1) k0=0

        do k=k0,gridobj%kml
        do j=1,1
        do i=i0,gridobj%iml
          read(1001,*) a(i,j,k) 
        enddo
        enddo
        enddo

        read(1001,*) ! Skipping </DataArray>
        read(1001,*) ! Skipping </PointData>
        read(1001,*) ! Skipping </Piece>
        read(1001,*) ! Skipping </RectilinearGrid>
        read(1001,*) ! Skipping </VTKFile>

        close(1001)
      endif !if (fexist .eqv. .false.) then ... else
        
      return 
      end subroutine read_vtr_r8_2d


      ! ----------------------------------------------------------------

      end program main_mini 