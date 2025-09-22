      module pvoutputs

        use mpi
        use grid_module 
        use boundary_conditions
        use cla
        implicit none
        character(len=300) :: PVOPATH='./outputs/'
        character(len=300) :: CHECKPOINTPATH = './outputs/checkpoint/'
        integer(4),parameter :: pvdheaderlines = 3
        integer(4),parameter :: pvdfooterlines = 2

        real(8) :: x0_o, y0_o, z0_o  ! origin shift
        real(8) :: dx_o, dy_o, dz_o  ! grid sizes
        integer(4) :: pv_iml, pv_jml, pv_kml  ! # grid, local
        integer(4) :: pv_img, pv_jmg, pv_kmg  ! # grid, global
        integer(4) :: nprocs_o ! # of processors

        character(len=7), parameter :: fmt_r7  = '(G0.7)'
        character(len=7), parameter :: fmt_r16 = '(G0.16)'
        character(len=7) :: fmt_r = fmt_r7  ! Default, 7 significant digits for real valued data

        type :: pvhandle
          real(8) :: x0,y0,z0
          real(8) :: dx,dy,dz
          integer(4) :: iml,jml,kml
          integer(4) :: img,jmg,kmg
          integer(4) :: nprocs
        end type

        type :: pvddata
          real(8) :: timestep
          integer(4) :: group 
          integer(4) :: part 
          character(len=300) :: file 
          character(len=30) :: name 
        end type

        !---------------------------------------------------------------
        ! Included subroutine list
        !  get_file_lines
        !  output_i4_raw
        !  output_pvd
        !  output_i4_pvtr
        !  output_i4_pvtr_2d
        !  output_r8_raw
        !  output_r8_raw_2d
        !  output_r8_pvtr
        !  output_r8_pvtr_2d
        !  checkpoint_i4_pvtr
        !  checkpoint_i4_pvtr_2d
        !  checkpoint_r8_pvtr
        !  checkpoint_r8_pvtr_2d
      contains
      subroutine init_pvoutputs(xo,yo,zo,dx,dy,dz,imloc,jmloc,kmloc,
     &                                          img, jmg, kmg, numprocs)
      implicit none
      real(8),intent(in) :: xo, yo, zo
      real(8),intent(in) :: dx, dy, dz
      integer(4),intent(in) :: imloc, jmloc, kmloc
      integer(4),intent(in) :: img, jmg, kmg
      integer(4), intent(in) :: numprocs

      x0_o = xo
      y0_o = yo
      z0_o = zo
      dx_o = dx
      dy_o = dy
      dz_o = dz
      nprocs_o = numprocs
      pv_iml = imloc
      pv_jml = jmloc
      pv_kml = kmloc
      pv_img = img
      pv_jmg = jmg
      pv_kmg = kmg

      return
      end subroutine init_pvoutputs

      ! ----------------------------------------------------------------

      subroutine init_pvhandle(handle,xo,yo,zo,dx,dy,dz,
     &                       imloc,jmloc,kmloc, img, jmg, kmg, numprocs)
      implicit none
      type(pvhandle),intent(inout) :: handle
      real(8),intent(in) :: xo, yo, zo
      real(8),intent(in) :: dx, dy, dz
      integer(4),intent(in) :: imloc, jmloc, kmloc
      integer(4),intent(in) :: img, jmg, kmg
      integer(4), intent(in) :: numprocs

      handle%x0 = xo
      handle%y0 = yo
      handle%z0 = zo
      handle%dx = dx
      handle%dy = dy
      handle%dz = dz
      handle%nprocs = numprocs
      handle%iml = imloc
      handle%jml = jmloc
      handle%kml = kmloc
      handle%img = img
      handle%jmg = jmg
      handle%kmg = kmg

      return
      end subroutine init_pvhandle

      !-----------------------------------------------------------------

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

      !-----------------------------------------------------------------

      subroutine output_i4_raw(a, ng, fnamebase, timestep, 
     & gpid, description)
      implicit none
      integer(4),intent(in) :: 
     &                   a(1-ng:pv_iml+ng,1-ng:pv_jml+ng,1-ng:pv_kml+ng)
      integer(4),intent(in) :: ng
      character(len=*),intent(in) :: fnamebase
      integer(4),intent(in) :: timestep, gpid
      character(len=*),intent(in) :: description

      integer(4) :: i,j,k
      character(len=300) :: fnameout
      character(len=8) :: ts2s, pid2s
      integer(4) :: pidx, pidy, pidz
      integer(4) :: il, iu, jl, ju, kl, ku

      call get_pidxyz(pidx, pidy, pidz, gpid)
      call get_ijkrange(il, iu, jl, ju, kl, ku, gpid)

      write(ts2s,"(i8.8)") int(timestep)
      write(pid2s,"(i4.4)") int(gpid)
      fnameout = trim(PVOPATH)//'/'//trim(fnamebase)//'_'//trim(ts2s)//
     & '_'//trim(pid2s)//'.raw'

      open(1000,file=trim(fnameout),status='unknown',action='write')
      ! Header output
      write(1000,*)description
      write(1000,*)
     & 'GlobalPID, PIDx, PIDy, PIDz, il, iu, jl, ju, kl, ku, SortOrder'
      write(1000,*)
     & gpid, pidx, pidy, pidz, il, iu, jl, ju, kl, ku, 'fortran'

      ! Data output - fortran order
      do k=1,pv_kml
      do j=1,pv_jml
      do i=1,pv_iml
            write(1000,*) a(i,j,k)
      enddo
      enddo
      enddo

      close(1000)

      end subroutine output_i4_raw

      !-----------------------------------------------------------------

      subroutine output_pvd(path,fnamebase,subdir,timestep,time,
     &                                                      description)
      implicit none
      character(len=*),intent(in) :: path 
      character(len=*),intent(in) :: fnamebase
      character(len=*),intent(in) :: subdir
      integer(4),intent(in) :: timestep
      real(8),intent(in) :: time
      character(len=*),intent(in) :: description

      character(len=300) :: fnamepvtr, fnamepvd
      character(len=8) :: ts2s
      character(len=20) :: time_str
      logical :: pvdexist
      integer(4) :: nlines, nbody, line
      write(ts2s,"(i8.8)") int(timestep)
      write(time_str,"(E20.14)") time
      fnamepvtr = trim(subdir)//trim(fnamebase)//'_'//
     & trim(ts2s)//'.pvtr'
      fnamepvd=trim(path)//'/'//trim(fnamebase)//'.pvd'

      ! Check if the file exists or not
      inquire(file=trim(fnamepvd),exist=pvdexist)
      if (pvdexist .eqv. .false.)then  ! .pvd doesn't exist. Make a new one.

        open(unit=900,file=trim(fnamepvd),
     &      status='unknown', action='write')

        ! Write .pvd header
        write(900,'(A)')'<?xml version="1.0"?>'
        write(900,'(A)')'<VTKFile type="Collection" version="0.1" '//
     &    'byte_order="LittleEndian">'
        write(900,'(A)')'<Collection>'

        ! Write the pvd data
        write(900,'(A)')'<DataSet timestep="'//trim(time_str)//
     &    '" group="" part="0" file="'//trim(fnamepvtr)//'" name="'//
     &    trim(description)//'"/>'

        ! Write .pvd footer
        write(900,'(A)')'</Collection>'
        write(900,'(A)')'</VTKFile>'
        close(900)

      else ! .pvd already exists. Update the lines.
        ! Find the # of lines in the pvd file
        call get_file_lines(nlines, trim(fnamepvd))
        nbody = nlines - pvdheaderlines - pvdfooterlines

        open(unit=900,file=trim(fnamepvd),
     &      status='old', action='readwrite')

        ! Skip the lines over prior data
        do line = 1, pvdheaderlines + nbody
          read(900,*)
        enddo

        ! Write the pvd data
        write(900,'(A)')'<DataSet timestep="'//trim(time_str)//
     &    '" group="" part="0" file="'//trim(fnamepvtr)//'" name="'//
     &    trim(description)//'"/>'

        ! Write .pvd footer
        write(900,'(A)')'</Collection>'
        write(900,'(A)')'</VTKFile>'
        close(900)

      endif ! if (pdvexist .eqv. .false.)then ... else

      return
      end subroutine output_pvd

      !-----------------------------------------------------------------

      subroutine output_i4_pvtr(a, ng, fnamebase,timestep,timeval,gpid,
     & description)
      ! Output of pvtr requires at least one ghost layer (ng>=1)

      implicit none
      integer(4),intent(in) :: 
     &                   a(1-ng:pv_iml+ng,1-ng:pv_jml+ng,1-ng:pv_kml+ng)
      integer(4),intent(in) :: ng
      character(len=*),intent(in) :: fnamebase
      integer(4),intent(in) :: timestep
      real(8),intent(in) :: timeval
      integer(4),intent(in) :: gpid
      character(len=*),intent(in) :: description

      integer(4) :: i,j,k,nn
      character(len=300) :: fnameg, fnamel, subdir
      character(len=8) :: ts2s, pid2s
      integer(4) :: pidx, pidy, pidz
      integer(4) :: il, iu, jl, ju, kl, ku
      integer(4) :: itmp, jtmp, ktmp
      character(len=300) :: str_extent_g, str_extent_l
      integer(4) :: ierr

      call get_pidxyz(pidx, pidy, pidz, gpid)
      call get_ijkrange(il, iu, jl, ju, kl, ku, gpid)

      subdir = trim(fnamebase)//'/'

      ! Make a subdir
      if(gpid.eq.0)then
        call execute_command_line(
     &                         'mkdir -p '//trim(PVOPATH)//trim(subdir))
      endif
      call mpi_barrier(MPI_COMM_WORLD, ierr)

      write(ts2s,"(i8.8)") int(timestep)
      write(pid2s,"(i4.4)") int(gpid)
      write(str_extent_g,"(6(i7,1x))")
     &                            0,pv_img-1,0,pv_jmg-1,0,pv_kmg-1
      fnameg = trim(PVOPATH)//trim(subdir)//trim(fnamebase)//'_'//
     & trim(ts2s)//'.pvtr'

      ! Output of pvtr file (meta file to combine each local data(vtr))
      if(gpid.eq.0)then
        call output_pvd(PVOPATH, fnamebase,subdir,timestep,timeval,
     &                                                      description)
        open(1001,FILE=trim(fnameg),status='unknown',action='write')
        write(1001,'(A)')'<?xml version="1.0"?>'
        write(1001,'(A)')'<VTKFile type="PRectilinearGrid" '//
     &   'version="0.1" byte_order="LittleEndian">'
        write(1001,'(A)')'<PRectilinearGrid WholeExtent="'//
     &   trim(adjustl(str_extent_g))//'">'

        write(1001,'(A)')'<PCoordinates>'
        write(1001,'(A)')
     &   '<DataArray Name="x" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="y" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="z" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PCoordinates>'

        write(1001,'(A)')'<PPointData Scalars="'//
     &   trim(fnamebase)//'">'
        write(1001,'(A)')'<DataArray type="Int32" Name="'//
     &   trim(fnamebase)//'" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PPointData>'

        do nn=0,nprocs_o-1
          write(pid2s,"(i4.4)") int(nn)
          fnamel=trim(fnamebase)//'_'//trim(ts2s)//'_'//
     &     trim(pid2s)//'.vtr'
          call get_ijkrange(il, iu, jl, ju, kl, ku, nn)
          call get_pidxyz(pidx, pidy, pidz, nn)

          itmp=il-1
          jtmp=jl-1
          ktmp=kl-1
          if(pidx.ge.1)itmp=il-2
          if(pidy.ge.1)jtmp=jl-2
          if(pidz.ge.1)ktmp=kl-2
          write(str_extent_l,"(6(i7,1x))")
     &      itmp,iu-1,jtmp,ju-1,ktmp,ku-1

          write(1001,'(A)')
     &     '<Piece Extent="'//trim(adjustl(str_extent_l))//
     &                '" Source="'//trim(fnamel)//'">'
          write(1001,'(A)')'</Piece>'
        enddo !do nn=0,nprocs_o-1

        write(1001,'(A)')'</PRectilinearGrid>'
        write(1001,'(A)')'</VTKFile>'
        close(1001)      
      endif !if(gpid.eq.0)then


      ! Output of vtr file (actual local data)
      write(pid2s,"(i4.4)") int(gpid)
      call get_ijkrange(il, iu, jl, ju, kl, ku, gpid)
      call get_pidxyz(pidx, pidy, pidz, gpid)

      itmp=il-1
      jtmp=jl-1
      ktmp=kl-1
      if(pidx.ge.1)itmp=il-2
      if(pidy.ge.1)jtmp=jl-2
      if(pidz.ge.1)ktmp=kl-2
      write(str_extent_l,"(6(i7,1x))") itmp,iu-1,jtmp,ju-1,ktmp,ku-1
      fnamel=trim(PVOPATH)//trim(subdir)//trim(fnamebase)//
     & '_'//trim(ts2s)//'_'//trim(pid2s)//'.vtr'
      open(1002,FILE=trim(fnamel),status='unknown',action='write')
      write(1002,'(A)')'<?xml version="1.0"?>'
      write(1002,'(A)')'<VTKFile type="RectilinearGrid" '//
     & 'version="0.1" byte_order="LittleEndian">'
      write(1002,'(A)')'<RectilinearGrid WholeExtent="'//
     & trim(adjustl(str_extent_g))//'">'
      write(1002,'(A)')'<Piece Extent="'//
     & trim(adjustl(str_extent_l))//'">'
      write(1002,'(A)')'<Coordinates>'
      write(1002,'(A)')
     & '<DataArray Name="x" type="Float32" format="ascii">'
      do i=itmp,iu-1
            write(1002,fmt_r) x0_o + i*dx_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="y" type="Float32" format="ascii">'
      do j=jtmp,ju-1
            write(1002,fmt_r) y0_o + j*dy_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="z" type="Float32" format="ascii">'
      do k=ktmp,ku-1
            write(1002,fmt_r) z0_o + k*dz_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</Coordinates>'

      write(1002,'(A)')'<PointData Scalars="'//trim(fnamebase)//'">'
      write(1002,'(A)')'<DataArray type="Int32" Name="'//
     & trim(fnamebase)//'" format="ascii">'
      itmp=1-1
      jtmp=1-1
      ktmp=1-1
      if(pidx.ge.1)itmp=1-2
      if(pidy.ge.1)jtmp=1-2
      if(pidz.ge.1)ktmp=1-2
      do k=ktmp,pv_kml-1
      do j=jtmp,pv_jml-1
      do i=itmp,pv_iml-1
            write(1002,'(I0)')a(i+1,j+1,k+1)
      enddo
      enddo
      enddo

      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</PointData>'
      write(1002,'(A)')'</Piece>'
      write(1002,'(A)')'</RectilinearGrid>'
      write(1002,'(A)')'</VTKFile>'
      close(1002)

      end subroutine output_i4_pvtr

      ! ----------------------------------------------------------------

      subroutine output_i4_pvtr_with_grid(gridobj,iml,jml,kml,ng,a,
     &                      fnamebase,timestep,timeval,gpid,description)
      ! Output of pvtr requires at least one ghost layer (ng>=1)

      implicit none
      type(grid),intent(in) :: gridobj
      integer(4),intent(in) :: iml,jml,kml
      integer(4),intent(in) :: ng
      integer(4),intent(in) :: a(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)
      character(len=*),intent(in) :: fnamebase
      integer(4),intent(in) :: timestep
      real(8),intent(in) :: timeval
      integer(4),intent(in) :: gpid
      character(len=*),intent(in) :: description

      integer(4) :: i,j,k,nn
      character(len=300) :: fnameg, fnamel, subdir
      character(len=8) :: ts2s, pid2s
      integer(4) :: pidx, pidy, pidz
      integer(4) :: il, iu, jl, ju, kl, ku
      integer(4) :: itmp, jtmp, ktmp
      character(len=300) :: str_extent_g, str_extent_l
      integer(4) :: ierr
      integer(4) :: img, jmg, kmg
      real(8) :: x0, y0, z0
      real(8) :: dxl, dyl, dzl 

      if(gridobj%is_initialized .eqv. .false.)then
        write(*,*)"Error, M801 pvoutputs, output_i4_pvtr_with_grid()"
        write(*,*)"  Uninitialized grid object delivered."
        write(*,*)"  Aborting..."
        call abort()
      endif

      x0 = gridobj%x0 
      y0 = gridobj%y0 
      z0 = gridobj%z0 
      dxl = gridobj%dxl 
      dyl = gridobj%dyl 
      dzl = gridobj%dzl 

      img = gridobj%img 
      jmg = gridobj%jmg 
      kmg = gridobj%kmg 

      pidx = gridobj%pidx
      pidy = gridobj%pidy
      pidz = gridobj%pidz
      il = gridobj%il 
      iu = gridobj%iu 
      jl = gridobj%jl 
      ju = gridobj%ju 
      kl = gridobj%kl 
      ku = gridobj%ku 

      subdir = trim(fnamebase)//'/'

      ! Make a subdir
      if(gpid.eq.0)then
        call execute_command_line(
     &                         'mkdir -p '//trim(PVOPATH)//trim(subdir))
      endif
      call mpi_barrier(MPI_COMM_WORLD, ierr)

      write(ts2s,"(i8.8)") int(timestep)
      write(pid2s,"(i4.4)") int(gpid)
      write(str_extent_g,"(6(i7,1x))")
     &                            0,img-1,0,jmg-1,0,kmg-1
      fnameg = trim(PVOPATH)//trim(subdir)//trim(fnamebase)//'_'//
     & trim(ts2s)//'.pvtr'

      ! Output of pvtr file (meta file to combine each local data(vtr))
      if(gpid.eq.0)then
        call output_pvd(PVOPATH, fnamebase,subdir,timestep,timeval,
     &                                                      description)
        open(1001,FILE=trim(fnameg),status='unknown',action='write')
        write(1001,'(A)')'<?xml version="1.0"?>'
        write(1001,'(A)')'<VTKFile type="PRectilinearGrid" '//
     &   'version="0.1" byte_order="LittleEndian">'
        write(1001,'(A)')'<PRectilinearGrid WholeExtent="'//
     &   trim(adjustl(str_extent_g))//'">'

        write(1001,'(A)')'<PCoordinates>'
        write(1001,'(A)')
     &   '<DataArray Name="x" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="y" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="z" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PCoordinates>'

        write(1001,'(A)')'<PPointData Scalars="'//
     &   trim(fnamebase)//'">'
        write(1001,'(A)')'<DataArray type="Int32" Name="'//
     &   trim(fnamebase)//'" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PPointData>'

        do nn=0,nprocs_o-1
          write(pid2s,"(i4.4)") int(nn)
          fnamel=trim(fnamebase)//'_'//trim(ts2s)//'_'//
     &     trim(pid2s)//'.vtr'
          call grid_get_ijkrange(il, iu, jl, ju, kl, ku, gridobj, nn)
          call grid_get_pidxyz(pidx, pidy, pidz, nn, gridobj)

          itmp=il-1
          jtmp=jl-1
          ktmp=kl-1
          if(pidx.ge.1)itmp=il-2
          if(pidy.ge.1)jtmp=jl-2
          if(pidz.ge.1)ktmp=kl-2
          write(str_extent_l,"(6(i7,1x))")
     &      itmp,iu-1,jtmp,ju-1,ktmp,ku-1

          write(1001,'(A)')
     &     '<Piece Extent="'//trim(adjustl(str_extent_l))//
     &                '" Source="'//trim(fnamel)//'">'
          write(1001,'(A)')'</Piece>'
        enddo !do nn=0,nprocs_o-1

        write(1001,'(A)')'</PRectilinearGrid>'
        write(1001,'(A)')'</VTKFile>'
        close(1001)      
      endif !if(gpid.eq.0)then


      ! Output of vtr file (actual local data)
      write(pid2s,"(i4.4)") int(gpid)
      call grid_get_ijkrange(il, iu, jl, ju, kl, ku, gridobj, gpid)
      call grid_get_pidxyz(pidx, pidy, pidz, gpid, gridobj)

      itmp=il-1
      jtmp=jl-1
      ktmp=kl-1
      if(pidx.ge.1)itmp=il-2
      if(pidy.ge.1)jtmp=jl-2
      if(pidz.ge.1)ktmp=kl-2
      write(str_extent_l,"(6(i7,1x))") itmp,iu-1,jtmp,ju-1,ktmp,ku-1
      fnamel=trim(PVOPATH)//trim(subdir)//trim(fnamebase)//
     & '_'//trim(ts2s)//'_'//trim(pid2s)//'.vtr'
      open(1002,FILE=trim(fnamel),status='unknown',action='write')
      write(1002,'(A)')'<?xml version="1.0"?>'
      write(1002,'(A)')'<VTKFile type="RectilinearGrid" '//
     & 'version="0.1" byte_order="LittleEndian">'
      write(1002,'(A)')'<RectilinearGrid WholeExtent="'//
     & trim(adjustl(str_extent_g))//'">'
      write(1002,'(A)')'<Piece Extent="'//
     & trim(adjustl(str_extent_l))//'">'
      write(1002,'(A)')'<Coordinates>'
      write(1002,'(A)')
     & '<DataArray Name="x" type="Float32" format="ascii">'
      do i=itmp,iu-1
            write(1002,fmt_r) x0 + i*dxl
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="y" type="Float32" format="ascii">'
      do j=jtmp,ju-1
            write(1002,fmt_r) y0 + j*dyl
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="z" type="Float32" format="ascii">'
      do k=ktmp,ku-1
            write(1002,fmt_r) z0 + k*dzl
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</Coordinates>'

      write(1002,'(A)')'<PointData Scalars="'//trim(fnamebase)//'">'
      write(1002,'(A)')'<DataArray type="Int32" Name="'//
     & trim(fnamebase)//'" format="ascii">'
      itmp=1-1
      jtmp=1-1
      ktmp=1-1
      if(pidx.ge.1)itmp=1-2
      if(pidy.ge.1)jtmp=1-2
      if(pidz.ge.1)ktmp=1-2
      do k=ktmp,kml-1
      do j=jtmp,jml-1
      do i=itmp,iml-1
            write(1002,'(I0)')a(i+1,j+1,k+1)
      enddo
      enddo
      enddo

      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</PointData>'
      write(1002,'(A)')'</Piece>'
      write(1002,'(A)')'</RectilinearGrid>'
      write(1002,'(A)')'</VTKFile>'
      close(1002)

      end subroutine output_i4_pvtr_with_grid

      !-----------------------------------------------------------------

      subroutine output_i4_pvtr_2d(a,ng,fnamebase,timestep,timeval,gpid,
     & description)
      ! Output of pvtr requires at least one ghost layer (ng>=1)

      implicit none
      integer(4),intent(in) :: a(1-ng:pv_iml+ng,1:1,1-ng:pv_kml+ng)
      integer(4),intent(in) :: ng
      character(len=*),intent(in) :: fnamebase
      integer(4),intent(in) :: timestep
      real(8),intent(in) :: timeval
      integer(4),intent(in) :: gpid
      character(len=*),intent(in) :: description

      integer(4) :: i,j,k,nn
      character(len=300) :: fnameg, fnamel, subdir
      character(len=8) :: ts2s, pid2s
      integer(4) :: pidx, pidy, pidz
      integer(4) :: il, iu, jl, ju, kl, ku
      integer(4) :: itmp, jtmp, ktmp
      character(len=300) :: str_extent_g, str_extent_l
      integer(4) :: ierr

      call get_pidxyz(pidx, pidy, pidz, gpid)
      call get_ijkrange(il, iu, jl, ju, kl, ku, gpid)

      subdir = trim(fnamebase)//'/'

      ! Make a subdir
      if(gpid.eq.0)then
        call execute_command_line(
     &                         'mkdir -p '//trim(PVOPATH)//trim(subdir))
      endif
      call mpi_barrier(MPI_COMM_WORLD, ierr)

      write(ts2s,"(i8.8)") int(timestep)
      write(pid2s,"(i4.4)") int(gpid)
      write(str_extent_g,"(6(i7,1x))")
     &                            0,pv_img-1,0,pv_jmg-1,0,pv_kmg-1
      fnameg = trim(PVOPATH)//trim(subdir)//trim(fnamebase)//'_'//
     & trim(ts2s)//'.pvtr'

      ! Output of pvtr file (meta file to combine each local data(vtr))
      if(gpid.eq.0)then
        call output_pvd(PVOPATH, fnamebase,subdir,timestep,timeval,
     &                                                      description)
        open(1001,FILE=trim(fnameg),status='unknown',action='write')
        write(1001,'(A)')'<?xml version="1.0"?>'
        write(1001,'(A)')'<VTKFile type="PRectilinearGrid" '//
     &   'version="0.1" byte_order="LittleEndian">'
        write(1001,'(A)')'<PRectilinearGrid WholeExtent="'//
     &   trim(adjustl(str_extent_g))//'">'

        write(1001,'(A)')'<PCoordinates>'
        write(1001,'(A)')
     &   '<DataArray Name="x" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="y" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="z" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PCoordinates>'

        write(1001,'(A)')'<PPointData Scalars="'//
     &   trim(fnamebase)//'">'
        write(1001,'(A)')'<DataArray type="Int32" Name="'//
     &   trim(fnamebase)//'" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PPointData>'

        do nn=0,nprocs_o-1
          write(pid2s,"(i4.4)") int(nn)
          fnamel=trim(fnamebase)//'_'//trim(ts2s)//'_'//
     &     trim(pid2s)//'.vtr'
          call get_ijkrange(il, iu, jl, ju, kl, ku, nn)
          call get_pidxyz(pidx, pidy, pidz, nn)

          itmp=il-1
          jtmp=jl-1
          ktmp=kl-1
          if(pidx.ge.1)itmp=il-2
          if(pidy.ge.1)jtmp=jl-2
          if(pidz.ge.1)ktmp=kl-2
          write(str_extent_l,"(6(i7,1x))")
     &      itmp,iu-1,jtmp,ju-1,ktmp,ku-1

          write(1001,'(A)')
     &     '<Piece Extent="'//trim(adjustl(str_extent_l))//
     &                '" Source="'//trim(fnamel)//'">'
          write(1001,'(A)')'</Piece>'
        enddo !do nn=0,nprocs_o-1

        write(1001,'(A)')'</PRectilinearGrid>'
        write(1001,'(A)')'</VTKFile>'
        close(1001)      
      endif !if(gpid.eq.0)then


      ! Output of vtr file (actual local data)
      write(pid2s,"(i4.4)") int(gpid)
      call get_ijkrange(il, iu, jl, ju, kl, ku, gpid)
      call get_pidxyz(pidx, pidy, pidz, gpid)

      itmp=il-1
      jtmp=jl-1
      ktmp=kl-1
      if(pidx.ge.1)itmp=il-2
      if(pidy.ge.1)jtmp=jl-2
      if(pidz.ge.1)ktmp=kl-2
      write(str_extent_l,"(6(i7,1x))") itmp,iu-1,jtmp,ju-1,ktmp,ku-1
      fnamel=trim(PVOPATH)//trim(subdir)//trim(fnamebase)//
     & '_'//trim(ts2s)//'_'//trim(pid2s)//'.vtr'
      open(1002,FILE=trim(fnamel),status='unknown',action='write')
      write(1002,'(A)')'<?xml version="1.0"?>'
      write(1002,'(A)')'<VTKFile type="RectilinearGrid" '//
     & 'version="0.1" byte_order="LittleEndian">'
      write(1002,'(A)')'<RectilinearGrid WholeExtent="'//
     & trim(adjustl(str_extent_g))//'">'
      write(1002,'(A)')'<Piece Extent="'//
     & trim(adjustl(str_extent_l))//'">'
      write(1002,'(A)')'<Coordinates>'
      write(1002,'(A)')
     & '<DataArray Name="x" type="Float32" format="ascii">'
      do i=itmp,iu-1
            write(1002,fmt_r) x0_o + i*dx_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="y" type="Float32" format="ascii">'
      do j=jtmp,ju-1
            write(1002,fmt_r) y0_o + j*dy_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="z" type="Float32" format="ascii">'
      do k=ktmp,ku-1
            write(1002,fmt_r) z0_o + k*dz_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</Coordinates>'

      write(1002,'(A)')'<PointData Scalars="'//trim(fnamebase)//'">'
      write(1002,'(A)')'<DataArray type="Int32" Name="'//
     & trim(fnamebase)//'" format="ascii">'
      itmp=1-1
      jtmp=1-1
      ktmp=1-1
      if(pidx.ge.1)itmp=1-2
      if(pidy.ge.1)jtmp=1-2
      if(pidz.ge.1)ktmp=1-2
      do k=ktmp,pv_kml-1
      do j=1,1
      do i=itmp,pv_iml-1
            write(1002,'(I0)')a(i+1,j,k+1)
      enddo
      enddo
      enddo

      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</PointData>'
      write(1002,'(A)')'</Piece>'
      write(1002,'(A)')'</RectilinearGrid>'
      write(1002,'(A)')'</VTKFile>'
      close(1002)

      end subroutine output_i4_pvtr_2d

      ! ----------------------------------------------------------------

      subroutine output_i4_pvtr_2d_with_grid(gridobj,iml,jml,kml,ng,a,
     &                      fnamebase,timestep,timeval,gpid,description)
      ! Output of pvtr requires at least one ghost layer (ng>=1)

      implicit none
      type(grid),intent(in) :: gridobj
      integer(4),intent(in) :: iml,jml,kml
      integer(4),intent(in) :: ng
      integer(4),intent(in) :: a(1-ng:iml+ng,1:1,1-ng:kml+ng)
      character(len=*),intent(in) :: fnamebase
      integer(4),intent(in) :: timestep
      real(8),intent(in) :: timeval
      integer(4),intent(in) :: gpid
      character(len=*),intent(in) :: description

      integer(4) :: i,j,k,nn
      character(len=300) :: fnameg, fnamel, subdir
      character(len=8) :: ts2s, pid2s
      integer(4) :: pidx, pidy, pidz
      integer(4) :: il, iu, jl, ju, kl, ku
      integer(4) :: itmp, jtmp, ktmp
      character(len=300) :: str_extent_g, str_extent_l
      integer(4) :: ierr
      integer(4) :: img, jmg, kmg
      real(8) :: x0, y0, z0
      real(8) :: dxl, dyl, dzl 

      if(gridobj%is_initialized .eqv. .false.)then
        write(*,*)"Error, M801 pvoutputs, output_i4_pvtr_with_grid()"
        write(*,*)"  Uninitialized grid object delivered."
        write(*,*)"  Aborting..."
        call abort()
      endif

      x0 = gridobj%x0 
      y0 = gridobj%y0 
      z0 = gridobj%z0 
      dxl = gridobj%dxl 
      dyl = gridobj%dyl 
      dzl = gridobj%dzl 

      img = gridobj%img 
      jmg = gridobj%jmg 
      kmg = gridobj%kmg 

      pidx = gridobj%pidx
      pidy = gridobj%pidy
      pidz = gridobj%pidz
      il = gridobj%il 
      iu = gridobj%iu 
      jl = gridobj%jl 
      ju = gridobj%ju 
      kl = gridobj%kl 
      ku = gridobj%ku 

      subdir = trim(fnamebase)//'/'

      ! Make a subdir
      if(gpid.eq.0)then
        call execute_command_line(
     &                         'mkdir -p '//trim(PVOPATH)//trim(subdir))
      endif
      call mpi_barrier(MPI_COMM_WORLD, ierr)

      write(ts2s,"(i8.8)") int(timestep)
      write(pid2s,"(i4.4)") int(gpid)
      write(str_extent_g,"(6(i7,1x))")
     &                            0,img-1,0,jmg-1,0,kmg-1
      fnameg = trim(PVOPATH)//trim(subdir)//trim(fnamebase)//'_'//
     & trim(ts2s)//'.pvtr'

      ! Output of pvtr file (meta file to combine each local data(vtr))
      if(gpid.eq.0)then
        call output_pvd(PVOPATH, fnamebase,subdir,timestep,timeval,
     &                                                      description)
        open(1001,FILE=trim(fnameg),status='unknown',action='write')
        write(1001,'(A)')'<?xml version="1.0"?>'
        write(1001,'(A)')'<VTKFile type="PRectilinearGrid" '//
     &   'version="0.1" byte_order="LittleEndian">'
        write(1001,'(A)')'<PRectilinearGrid WholeExtent="'//
     &   trim(adjustl(str_extent_g))//'">'

        write(1001,'(A)')'<PCoordinates>'
        write(1001,'(A)')
     &   '<DataArray Name="x" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="y" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="z" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PCoordinates>'

        write(1001,'(A)')'<PPointData Scalars="'//
     &   trim(fnamebase)//'">'
        write(1001,'(A)')'<DataArray type="Int32" Name="'//
     &   trim(fnamebase)//'" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PPointData>'

        do nn=0,nprocs_o-1
          write(pid2s,"(i4.4)") int(nn)
          fnamel=trim(fnamebase)//'_'//trim(ts2s)//'_'//
     &     trim(pid2s)//'.vtr'
          call grid_get_ijkrange(il, iu, jl, ju, kl, ku, gridobj, nn)
          call grid_get_pidxyz(pidx, pidy, pidz, nn, gridobj)

          itmp=il-1
          jtmp=jl-1
          ktmp=kl-1
          if(pidx.ge.1)itmp=il-2
          if(pidy.ge.1)jtmp=jl-2
          if(pidz.ge.1)ktmp=kl-2
          write(str_extent_l,"(6(i7,1x))")
     &      itmp,iu-1,jtmp,ju-1,ktmp,ku-1

          write(1001,'(A)')
     &     '<Piece Extent="'//trim(adjustl(str_extent_l))//
     &                '" Source="'//trim(fnamel)//'">'
          write(1001,'(A)')'</Piece>'
        enddo !do nn=0,nprocs_o-1

        write(1001,'(A)')'</PRectilinearGrid>'
        write(1001,'(A)')'</VTKFile>'
        close(1001)      
      endif !if(gpid.eq.0)then


      ! Output of vtr file (actual local data)
      write(pid2s,"(i4.4)") int(gpid)
      call grid_get_ijkrange(il, iu, jl, ju, kl, ku, gridobj, gpid)
      call grid_get_pidxyz(pidx, pidy, pidz, gpid, gridobj)

      itmp=il-1
      jtmp=jl-1
      ktmp=kl-1
      if(pidx.ge.1)itmp=il-2
      if(pidy.ge.1)jtmp=jl-2
      if(pidz.ge.1)ktmp=kl-2
      write(str_extent_l,"(6(i7,1x))") itmp,iu-1,jtmp,ju-1,ktmp,ku-1
      fnamel=trim(PVOPATH)//trim(subdir)//trim(fnamebase)//
     & '_'//trim(ts2s)//'_'//trim(pid2s)//'.vtr'
      open(1002,FILE=trim(fnamel),status='unknown',action='write')
      write(1002,'(A)')'<?xml version="1.0"?>'
      write(1002,'(A)')'<VTKFile type="RectilinearGrid" '//
     & 'version="0.1" byte_order="LittleEndian">'
      write(1002,'(A)')'<RectilinearGrid WholeExtent="'//
     & trim(adjustl(str_extent_g))//'">'
      write(1002,'(A)')'<Piece Extent="'//
     & trim(adjustl(str_extent_l))//'">'
      write(1002,'(A)')'<Coordinates>'
      write(1002,'(A)')
     & '<DataArray Name="x" type="Float32" format="ascii">'
      do i=itmp,iu-1
            write(1002,fmt_r) x0 + i*dxl
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="y" type="Float32" format="ascii">'
      do j=jtmp,ju-1
            write(1002,fmt_r) y0 + j*dyl
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="z" type="Float32" format="ascii">'
      do k=ktmp,ku-1
            write(1002,fmt_r) z0 + k*dzl
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</Coordinates>'

      write(1002,'(A)')'<PointData Scalars="'//trim(fnamebase)//'">'
      write(1002,'(A)')'<DataArray type="Int32" Name="'//
     & trim(fnamebase)//'" format="ascii">'
      itmp=1-1
      jtmp=1-1
      ktmp=1-1
      if(pidx.ge.1)itmp=1-2
      if(pidy.ge.1)jtmp=1-2
      if(pidz.ge.1)ktmp=1-2
      do k=ktmp,kml-1
      do j=1,1
      do i=itmp,iml-1
            write(1002,'(I0)')a(i+1,j,k+1)
      enddo
      enddo
      enddo

      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</PointData>'
      write(1002,'(A)')'</Piece>'
      write(1002,'(A)')'</RectilinearGrid>'
      write(1002,'(A)')'</VTKFile>'
      close(1002)

      end subroutine output_i4_pvtr_2d_with_grid

      ! ----------------------------------------------------------------

      subroutine output_i4_sub_pvtr_2d_with_grid(gridobj,iml,jml,kml,ng,
     &                    imin,imax,jmin,jmax,kmin,kmax,
     &                    a,fnamebase,timestep,timeval,gpid,description)
      ! Output of pvtr requires at least one ghost layer (ng>=1)

      implicit none
      type(grid),intent(in) :: gridobj
      integer(4),intent(in) :: iml,jml,kml
      integer(4),intent(in) :: ng
      integer(4),intent(in) :: imin,imax,jmin,jmax,kmin,kmax 
      integer(4),intent(in) :: a(1-ng:iml+ng,1:1,1-ng:kml+ng)
      character(len=*),intent(in) :: fnamebase
      integer(4),intent(in) :: timestep
      real(8),intent(in) :: timeval
      integer(4),intent(in) :: gpid
      character(len=*),intent(in) :: description

      integer(4) :: i,j,k,nn
      character(len=300) :: fnameg, fnamel, subdir
      character(len=8) :: ts2s, pid2s
      integer(4) :: pidx, pidy, pidz
      integer(4) :: il, iu, jl, ju, kl, ku
      integer(4) :: itmp, jtmp, ktmp
      character(len=300) :: str_extent_g, str_extent_l
      integer(4) :: ierr
      integer(4) :: img, jmg, kmg
      real(8) :: x0, y0, z0
      real(8) :: dxl, dyl, dzl 
      integer(4) :: ilg, iug, jlg, jug, klg, kug

      if(gridobj%is_initialized .eqv. .false.)then
        write(*,*)
     &           "Error, M801 pvoutputs, output_i4_sub,pvtr_with_grid()"
        write(*,*)"  Uninitialized grid object delivered."
        write(*,*)"  Aborting..."
        call abort()
      endif

      x0 = gridobj%x0 
      y0 = gridobj%y0 
      z0 = gridobj%z0 
      dxl = gridobj%dxl 
      dyl = gridobj%dyl 
      dzl = gridobj%dzl 

      ! img = gridobj%img 
      ! jmg = gridobj%jmg 
      ! kmg = gridobj%kmg 

      pidx = gridobj%pidx
      pidy = gridobj%pidy
      pidz = gridobj%pidz
      il = gridobj%il 
      iu = gridobj%iu 
      jl = gridobj%jl 
      ju = gridobj%ju 
      kl = gridobj%kl 
      ku = gridobj%ku 

      img = imax-imin+1
      jmg = jmax-jmin+1
      kmg = kmax-kmin+1

      if(gridobj%npx.ne.1.and.img.ne.gridobj%img .or. 
     &   gridobj%npy.ne.1.and.jmg.ne.gridobj%jmg .or.
     &   gridobj%npz.ne.1.and.kmg.ne.gridobj%kmg )then
        write(*,*)
     &           "Error, M801 pvoutputs, output_i4_sub,pvtr_with_grid()"
        write(*,*)"  Invalid NPX, NPY, or NPZ."
        write(*,*)"  To make output of the shrunken array (< *MORI), "
        write(*,*)"  NP* should be only 1 along the shrunken direction."
        write(*,*)"  Aborting..."
        call abort()
      endif

      subdir = trim(fnamebase)//'/'

      ilg = 0
      iug = img-1
      jlg = 0
      jug = jmg-1
      klg = 0
      kug = kmg-1

      if(img.lt.gridobj%img) then
        ilg = imin-1
        iug = imax-1
      endif

      if(jmg.lt.gridobj%jmg) then
        jlg = jmin-1
        jug = jmax-1
      endif

      if(kmg.lt.gridobj%kmg) then
        klg = kmin-1
        kug = kmax-1
      endif

      ! Make a subdir
      if(gpid.eq.0)then
        call execute_command_line(
     &                         'mkdir -p '//trim(PVOPATH)//trim(subdir))
      endif
      call mpi_barrier(MPI_COMM_WORLD, ierr)

      write(ts2s,"(i8.8)") int(timestep)
      write(pid2s,"(i4.4)") int(gpid)
      write(str_extent_g,"(6(i7,1x))") ilg,iug,jlg,jug,klg,kug
      fnameg = trim(PVOPATH)//trim(subdir)//trim(fnamebase)//'_'//
     & trim(ts2s)//'.pvtr'

      ! Output of pvtr file (meta file to combine each local data(vtr))
      if(gpid.eq.0)then
        call output_pvd(PVOPATH, fnamebase,subdir,timestep,timeval,
     &                                                      description)
        open(1001,FILE=trim(fnameg),status='unknown',action='write')
        write(1001,'(A)')'<?xml version="1.0"?>'
        write(1001,'(A)')'<VTKFile type="PRectilinearGrid" '//
     &   'version="0.1" byte_order="LittleEndian">'
        write(1001,'(A)')'<PRectilinearGrid WholeExtent="'//
     &   trim(adjustl(str_extent_g))//'">'

        write(1001,'(A)')'<PCoordinates>'
        write(1001,'(A)')
     &   '<DataArray Name="x" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="y" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="z" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PCoordinates>'

        write(1001,'(A)')'<PPointData Scalars="'//
     &   trim(fnamebase)//'">'
        write(1001,'(A)')'<DataArray type="Int32" Name="'//
     &   trim(fnamebase)//'" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PPointData>'

        do nn=0,nprocs_o-1
          write(pid2s,"(i4.4)") int(nn)
          fnamel=trim(fnamebase)//'_'//trim(ts2s)//'_'//
     &     trim(pid2s)//'.vtr'
          call grid_get_ijkrange(il, iu, jl, ju, kl, ku, gridobj, nn)
          call grid_get_pidxyz(pidx, pidy, pidz, nn, gridobj)

          if(img.lt.gridobj%img) then
            il = imin
            iu = imax 
          endif

          if(jmg.lt.gridobj%jmg) then
            jl = jmin
            ju = jmax 
          endif

          if(kmg.lt.gridobj%kmg) then
            kl = kmin
            ku = kmax 
          endif

          itmp=il-1
          jtmp=jl-1
          ktmp=kl-1
          if(pidx.ge.1)itmp=il-2
          if(pidy.ge.1)jtmp=jl-2
          if(pidz.ge.1)ktmp=kl-2
          write(str_extent_l,"(6(i7,1x))")
     &      itmp,iu-1,jtmp,ju-1,ktmp,ku-1

          write(1001,'(A)')
     &     '<Piece Extent="'//trim(adjustl(str_extent_l))//
     &                '" Source="'//trim(fnamel)//'">'
          write(1001,'(A)')'</Piece>'
        enddo !do nn=0,nprocs_o-1

        write(1001,'(A)')'</PRectilinearGrid>'
        write(1001,'(A)')'</VTKFile>'
        close(1001)      
      endif !if(gpid.eq.0)then


      ! Output of vtr file (actual local data)
      write(pid2s,"(i4.4)") int(gpid)
      call grid_get_ijkrange(il, iu, jl, ju, kl, ku, gridobj, gpid)
      call grid_get_pidxyz(pidx, pidy, pidz, gpid, gridobj)

      if(img.lt.gridobj%img) then
        il = imin
        iu = imax 
      endif

      if(jmg.lt.gridobj%jmg) then
        jl = jmin
        ju = jmax 
      endif

      if(kmg.lt.gridobj%kmg) then
        kl = kmin
        ku = kmax 
      endif

      itmp=il-1
      jtmp=jl-1
      ktmp=kl-1
      if(pidx.ge.1)itmp=il-2
      if(pidy.ge.1)jtmp=jl-2
      if(pidz.ge.1)ktmp=kl-2
      write(str_extent_l,"(6(i7,1x))") itmp,iu-1,jtmp,ju-1,ktmp,ku-1
      fnamel=trim(PVOPATH)//trim(subdir)//trim(fnamebase)//
     & '_'//trim(ts2s)//'_'//trim(pid2s)//'.vtr'
      open(1002,FILE=trim(fnamel),status='unknown',action='write')
      write(1002,'(A)')'<?xml version="1.0"?>'
      write(1002,'(A)')'<VTKFile type="RectilinearGrid" '//
     & 'version="0.1" byte_order="LittleEndian">'
      write(1002,'(A)')'<RectilinearGrid WholeExtent="'//
     & trim(adjustl(str_extent_g))//'">'
      write(1002,'(A)')'<Piece Extent="'//
     & trim(adjustl(str_extent_l))//'">'
      write(1002,'(A)')'<Coordinates>'
      write(1002,'(A)')
     & '<DataArray Name="x" type="Float32" format="ascii">'
      do i=itmp,iu-1
            write(1002,fmt_r) x0 + i*dxl
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="y" type="Float32" format="ascii">'
      do j=jtmp,ju-1
            write(1002,fmt_r) y0 + j*dyl
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="z" type="Float32" format="ascii">'
      do k=ktmp,ku-1
            write(1002,fmt_r) z0 + k*dzl
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</Coordinates>'

      write(1002,'(A)')'<PointData Scalars="'//trim(fnamebase)//'">'
      write(1002,'(A)')'<DataArray type="Int32" Name="'//
     & trim(fnamebase)//'" format="ascii">'
      itmp=1-1
      jtmp=1-1
      ktmp=1-1
      if(pidx.ge.1)itmp=1-2
      if(pidy.ge.1)jtmp=1-2
      if(pidz.ge.1)ktmp=1-2

      il=itmp 
      iu=iml-1
      jl=1
      ju=1
      kl=ktmp
      ku=kml-1
      if(img.lt.gridobj%img) then
        il = imin-1
        iu = imax-1
      endif
      if(jmg.lt.gridobj%jmg) then
        jl = jmin-1
        ju = jmax-1
      endif
      if(kmg.lt.gridobj%kmg) then
        kl = kmin-1
        ku = kmax-1
      endif

      do k=kl,ku
      do j=1,1
      do i=il,iu
            write(1002,'(I0)')a(i+1,j,k+1)
      enddo
      enddo
      enddo

      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</PointData>'
      write(1002,'(A)')'</Piece>'
      write(1002,'(A)')'</RectilinearGrid>'
      write(1002,'(A)')'</VTKFile>'
      close(1002)

      end subroutine output_i4_sub_pvtr_2d_with_grid

      ! ----------------------------------------------------------------

      subroutine output_r8_raw(a, ng, fnamebase, timestep, 
     & gpid, description)

      implicit none
      real(8),intent(in) :: 
     &                   a(1-ng:pv_iml+ng,1-ng:pv_jml+ng,1-ng:pv_kml+ng)
      integer(4),intent(in) :: ng
      character(len=*),intent(in) :: fnamebase
      integer(4),intent(in) :: timestep, gpid
      character(len=*),intent(in) :: description

      integer(4) :: i,j,k
      character(len=300) :: fnameout, subdir
      character(len=8) :: ts2s, pid2s
      integer(4) :: pidx, pidy, pidz
      integer(4) :: il, iu, jl, ju, kl, ku
      integer(4) :: ierr

      call get_pidxyz(pidx, pidy, pidz, gpid)
      call get_ijkrange(il, iu, jl, ju, kl, ku, gpid)

      subdir = trim(fnamebase)//'/'

      ! Make a subdir
      if(gpid.eq.0)then
        call execute_command_line(
     &                         'mkdir -p '//trim(PVOPATH)//trim(subdir))
      endif
      call mpi_barrier(MPI_COMM_WORLD, ierr)

      write(ts2s,"(i8.8)") int(timestep)
      write(pid2s,"(i4.4)") int(gpid)
      fnameout = trim(PVOPATH)//trim(subdir)//trim(fnamebase)//
     & '_'//trim(ts2s)//'_'//trim(pid2s)//'.raw'

      open(1000,file=trim(fnameout),status='unknown',action='write')
      ! Header output
      write(1000,*)description
      write(1000,*)
     & 'GlobalPID, PIDx, PIDy, PIDz, il, iu, jl, ju, kl, ku, SortOrder'
      write(1000,*)
     & gpid, pidx, pidy, pidz, il, iu, jl, ju, kl, ku, 'fortran'

      ! Data output - fortran order
      do k=1,pv_kml
      do j=1,pv_jml
      do i=1,pv_iml
            write(1000,*) a(i,j,k)
      enddo
      enddo
      enddo

      close(1000)

      end subroutine output_r8_raw

      !-----------------------------------------------------------------

      subroutine output_r8_raw_2d(a, ng, fnamebase, timestep, 
     & gpid, description)

      implicit none
      real(8),intent(in) :: a(1-ng:pv_iml+ng,1:1,1-ng:pv_kml+ng)
      integer(4),intent(in) :: ng
      character(len=*),intent(in) :: fnamebase
      integer(4),intent(in) :: timestep, gpid
      character(len=*),intent(in) :: description

      integer(4) :: i,j,k
      character(len=300) :: fnameout, subdir
      character(len=8) :: ts2s, pid2s
      integer(4) :: pidx, pidy, pidz
      integer(4) :: il, iu, jl, ju, kl, ku
      integer(4) :: ierr

      call get_pidxyz(pidx, pidy, pidz, gpid)
      call get_ijkrange(il, iu, jl, ju, kl, ku, gpid)

      subdir = trim(fnamebase)//'/'

      ! Make a subdir
      if(gpid.eq.0)then
        call execute_command_line(
     &                         'mkdir -p '//trim(PVOPATH)//trim(subdir))
      endif
      call mpi_barrier(MPI_COMM_WORLD, ierr)

      write(ts2s,"(i8.8)") int(timestep)
      write(pid2s,"(i4.4)") int(gpid)
      fnameout = trim(PVOPATH)//trim(subdir)//trim(fnamebase)//
     & '_'//trim(ts2s)//'_'//trim(pid2s)//'.raw'

      open(1000,file=trim(fnameout),status='unknown',action='write')
      ! Header output
      write(1000,*)description
      write(1000,*)
     & 'GlobalPID, PIDx, PIDy, PIDz, il, iu, jl, ju, kl, ku, SortOrder'
      write(1000,*)
     & gpid, pidx, pidy, pidz, il, iu, jl, ju, kl, ku, 'fortran'

      ! Data output - fortran order
      do k=1,pv_kml
      do j=1,pv_jml
      do i=1,pv_iml
            write(1000,*) a(i,j,k)
      enddo
      enddo
      enddo

      close(1000)

      end subroutine output_r8_raw_2d

      !-----------------------------------------------------------------

      subroutine output_r8_pvtr(a,ng,fnamebase,timestep,timeval,gpid,
     & description)
      ! Output of pvtr requires at least one ghost layer (ng>=1)

      implicit none
      real(8),intent(in) :: 
     &                   a(1-ng:pv_iml+ng,1-ng:pv_jml+ng,1-ng:pv_kml+ng)
      integer(4),intent(in) :: ng
      character(len=*),intent(in) :: fnamebase
      integer(4),intent(in) :: timestep
      real(8),intent(in) :: timeval
      integer(4),intent(in) :: gpid
      character(len=*),intent(in) :: description

      integer(4) :: i,j,k,nn
      character(len=300) :: fnameg, fnamel, subdir
      character(len=8) :: ts2s, pid2s
      integer(4) :: pidx, pidy, pidz
      integer(4) :: il, iu, jl, ju, kl, ku
      integer(4) :: itmp, jtmp, ktmp
      character(len=300) :: str_extent_g, str_extent_l
      integer(4) :: ierr

      call get_pidxyz(pidx, pidy, pidz, gpid)
      call get_ijkrange(il, iu, jl, ju, kl, ku, gpid)

      subdir = trim(fnamebase)//'/'

      ! Make a subdir
      if(gpid.eq.0)then
        call execute_command_line(
     &                         'mkdir -p '//trim(PVOPATH)//trim(subdir))
      endif
      call mpi_barrier(MPI_COMM_WORLD, ierr)

      write(ts2s,"(i8.8)") int(timestep)
      write(pid2s,"(i4.4)") int(gpid)
      write(str_extent_g,"(6(i7,1x))")
     &                            0,pv_img-1,0,pv_jmg-1,0,pv_kmg-1
      fnameg = trim(PVOPATH)//trim(subdir)//trim(fnamebase)//'_'//
     & trim(ts2s)//'.pvtr'

      ! Output of pvtr file (meta file to combine each local data(vtr))
      if(gpid.eq.0)then
        call output_pvd(PVOPATH, fnamebase,subdir,timestep,timeval,
     &                                                      description)
        open(1001,FILE=trim(fnameg),status='unknown',action='write')
        write(1001,'(A)')'<?xml version="1.0"?>'
        write(1001,'(A)')'<VTKFile type="PRectilinearGrid" '//
     &   'version="0.1" byte_order="LittleEndian">'
        write(1001,'(A)')'<PRectilinearGrid WholeExtent="'//
     &   trim(adjustl(str_extent_g))//'">'

        write(1001,'(A)')'<PCoordinates>'
        write(1001,'(A)')
     &   '<DataArray Name="x" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="y" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="z" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PCoordinates>'

        write(1001,'(A)')'<PPointData Scalars="'//
     &   trim(fnamebase)//'">'
        write(1001,'(A)')'<DataArray type="Float32" Name="'//
     &   trim(fnamebase)//'" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PPointData>'

        do nn=0,nprocs_o-1
          write(pid2s,"(i4.4)") int(nn)
          fnamel=trim(fnamebase)//'_'//trim(ts2s)//'_'//
     &     trim(pid2s)//'.vtr'
          call get_ijkrange(il, iu, jl, ju, kl, ku, nn)
          call get_pidxyz(pidx, pidy, pidz, nn)

          itmp=il-1
          jtmp=jl-1
          ktmp=kl-1
          if(pidx.ge.1)itmp=il-2
          if(pidy.ge.1)jtmp=jl-2
          if(pidz.ge.1)ktmp=kl-2
          write(str_extent_l,"(6(i7,1x))")
     &      itmp,iu-1,jtmp,ju-1,ktmp,ku-1

          write(1001,'(A)')
     &     '<Piece Extent="'//trim(adjustl(str_extent_l))//
     &                '" Source="'//trim(fnamel)//'">'
          write(1001,'(A)')'</Piece>'
        enddo !do nn=0,nprocs_o-1

        write(1001,'(A)')'</PRectilinearGrid>'
        write(1001,'(A)')'</VTKFile>'
        close(1001)      
      endif !if(gpid.eq.0)then


      ! Output of vtr file (actual local data)
      write(pid2s,"(i4.4)") int(gpid)
      call get_ijkrange(il, iu, jl, ju, kl, ku, gpid)
      call get_pidxyz(pidx, pidy, pidz, gpid)

      itmp=il-1
      jtmp=jl-1
      ktmp=kl-1
      if(pidx.ge.1)itmp=il-2
      if(pidy.ge.1)jtmp=jl-2
      if(pidz.ge.1)ktmp=kl-2
      write(str_extent_l,"(6(i7,1x))") itmp,iu-1,jtmp,ju-1,ktmp,ku-1
      fnamel=trim(PVOPATH)//trim(subdir)//trim(fnamebase)//
     & '_'//trim(ts2s)//'_'//trim(pid2s)//'.vtr'
      open(1002,FILE=trim(fnamel),status='unknown',action='write')
      write(1002,'(A)')'<?xml version="1.0"?>'
      write(1002,'(A)')'<VTKFile type="RectilinearGrid" '//
     & 'version="0.1" byte_order="LittleEndian">'
      write(1002,'(A)')'<RectilinearGrid WholeExtent="'//
     & trim(adjustl(str_extent_g))//'">'
      write(1002,'(A)')'<Piece Extent="'//
     & trim(adjustl(str_extent_l))//'">'
      write(1002,'(A)')'<Coordinates>'
      write(1002,'(A)')
     & '<DataArray Name="x" type="Float32" format="ascii">'
      do i=itmp,iu-1
            write(1002,fmt_r) x0_o + i*dx_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="y" type="Float32" format="ascii">'
      do j=jtmp,ju-1
            write(1002,fmt_r) y0_o + j*dy_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="z" type="Float32" format="ascii">'
      do k=ktmp,ku-1
            write(1002,fmt_r) z0_o + k*dz_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</Coordinates>'

      write(1002,'(A)')'<PointData Scalars="'//trim(fnamebase)//'">'
      write(1002,'(A)')'<DataArray type="Float32" Name="'//
     & trim(fnamebase)//'" format="ascii">'
      itmp=1-1
      jtmp=1-1
      ktmp=1-1
      if(pidx.ge.1)itmp=1-2
      if(pidy.ge.1)jtmp=1-2
      if(pidz.ge.1)ktmp=1-2
      do k=ktmp,pv_kml-1
      do j=jtmp,pv_jml-1
      do i=itmp,pv_iml-1
            write(1002,fmt_r)a(i+1,j+1,k+1)
      enddo
      enddo
      enddo

      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</PointData>'
      write(1002,'(A)')'</Piece>'
      write(1002,'(A)')'</RectilinearGrid>'
      write(1002,'(A)')'</VTKFile>'
      close(1002)

      end subroutine output_r8_pvtr

      ! ----------------------------------------------------------------

      subroutine output_r8_pvtr_with_grid(gridobj,iml,jml,kml,ng,a,
     &                      fnamebase,timestep,timeval,gpid,description)
      ! Output of pvtr requires at least one ghost layer (ng>=1)

      implicit none
      type(grid),intent(in) :: gridobj
      integer(4),intent(in) :: iml,jml,kml
      integer(4),intent(in) :: ng
      real(8),intent(in) :: a(1-ng:iml+ng,1-ng:jml+ng,1-ng:kml+ng)
      character(len=*),intent(in) :: fnamebase
      integer(4),intent(in) :: timestep
      real(8),intent(in) :: timeval
      integer(4),intent(in) :: gpid
      character(len=*),intent(in) :: description

      integer(4) :: i,j,k,nn
      character(len=300) :: fnameg, fnamel, subdir
      character(len=8) :: ts2s, pid2s
      integer(4) :: pidx, pidy, pidz
      integer(4) :: il, iu, jl, ju, kl, ku
      integer(4) :: itmp, jtmp, ktmp
      character(len=300) :: str_extent_g, str_extent_l
      integer(4) :: ierr
      integer(4) :: img, jmg, kmg
      real(8) :: x0, y0, z0
      real(8) :: dxl, dyl, dzl 

      if(gridobj%is_initialized .eqv. .false.)then
        write(*,*)"Error, M801 pvoutputs, output_i4_pvtr_with_grid()"
        write(*,*)"  Uninitialized grid object delivered."
        write(*,*)"  Aborting..."
        call abort()
      endif

      x0 = gridobj%x0 
      y0 = gridobj%y0 
      z0 = gridobj%z0 
      dxl = gridobj%dxl 
      dyl = gridobj%dyl 
      dzl = gridobj%dzl 

      img = gridobj%img 
      jmg = gridobj%jmg 
      kmg = gridobj%kmg 

      pidx = gridobj%pidx
      pidy = gridobj%pidy
      pidz = gridobj%pidz
      il = gridobj%il 
      iu = gridobj%iu 
      jl = gridobj%jl 
      ju = gridobj%ju 
      kl = gridobj%kl 
      ku = gridobj%ku 

      subdir = trim(fnamebase)//'/'

      ! Make a subdir
      if(gpid.eq.0)then
        call execute_command_line(
     &                         'mkdir -p '//trim(PVOPATH)//trim(subdir))
      endif
      call mpi_barrier(MPI_COMM_WORLD, ierr)

      write(ts2s,"(i8.8)") int(timestep)
      write(pid2s,"(i4.4)") int(gpid)
      write(str_extent_g,"(6(i7,1x))")
     &                            0,img-1,0,jmg-1,0,kmg-1
      fnameg = trim(PVOPATH)//trim(subdir)//trim(fnamebase)//'_'//
     & trim(ts2s)//'.pvtr'

      ! Output of pvtr file (meta file to combine each local data(vtr))
      if(gpid.eq.0)then
        call output_pvd(PVOPATH, fnamebase,subdir,timestep,timeval,
     &                                                      description)
        open(1001,FILE=trim(fnameg),status='unknown',action='write')
        write(1001,'(A)')'<?xml version="1.0"?>'
        write(1001,'(A)')'<VTKFile type="PRectilinearGrid" '//
     &   'version="0.1" byte_order="LittleEndian">'
        write(1001,'(A)')'<PRectilinearGrid WholeExtent="'//
     &   trim(adjustl(str_extent_g))//'">'

        write(1001,'(A)')'<PCoordinates>'
        write(1001,'(A)')
     &   '<DataArray Name="x" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="y" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="z" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PCoordinates>'

        write(1001,'(A)')'<PPointData Scalars="'//
     &   trim(fnamebase)//'">'
        write(1001,'(A)')'<DataArray type="Float32" Name="'//
     &   trim(fnamebase)//'" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PPointData>'

        do nn=0,nprocs_o-1
          write(pid2s,"(i4.4)") int(nn)
          fnamel=trim(fnamebase)//'_'//trim(ts2s)//'_'//
     &     trim(pid2s)//'.vtr'
          call grid_get_ijkrange(il, iu, jl, ju, kl, ku, gridobj, nn)
          call grid_get_pidxyz(pidx, pidy, pidz, nn, gridobj)

          itmp=il-1
          jtmp=jl-1
          ktmp=kl-1
          if(pidx.ge.1)itmp=il-2
          if(pidy.ge.1)jtmp=jl-2
          if(pidz.ge.1)ktmp=kl-2
          write(str_extent_l,"(6(i7,1x))")
     &      itmp,iu-1,jtmp,ju-1,ktmp,ku-1

          write(1001,'(A)')
     &     '<Piece Extent="'//trim(adjustl(str_extent_l))//
     &                '" Source="'//trim(fnamel)//'">'
          write(1001,'(A)')'</Piece>'
        enddo !do nn=0,nprocs_o-1

        write(1001,'(A)')'</PRectilinearGrid>'
        write(1001,'(A)')'</VTKFile>'
        close(1001)      
      endif !if(gpid.eq.0)then


      ! Output of vtr file (actual local data)
      write(pid2s,"(i4.4)") int(gpid)
      call grid_get_ijkrange(il, iu, jl, ju, kl, ku, gridobj, gpid)
      call grid_get_pidxyz(pidx, pidy, pidz, gpid, gridobj)

      itmp=il-1
      jtmp=jl-1
      ktmp=kl-1
      if(pidx.ge.1)itmp=il-2
      if(pidy.ge.1)jtmp=jl-2
      if(pidz.ge.1)ktmp=kl-2
      write(str_extent_l,"(6(i7,1x))") itmp,iu-1,jtmp,ju-1,ktmp,ku-1
      fnamel=trim(PVOPATH)//trim(subdir)//trim(fnamebase)//
     & '_'//trim(ts2s)//'_'//trim(pid2s)//'.vtr'
      open(1002,FILE=trim(fnamel),status='unknown',action='write')
      write(1002,'(A)')'<?xml version="1.0"?>'
      write(1002,'(A)')'<VTKFile type="RectilinearGrid" '//
     & 'version="0.1" byte_order="LittleEndian">'
      write(1002,'(A)')'<RectilinearGrid WholeExtent="'//
     & trim(adjustl(str_extent_g))//'">'
      write(1002,'(A)')'<Piece Extent="'//
     & trim(adjustl(str_extent_l))//'">'
      write(1002,'(A)')'<Coordinates>'
      write(1002,'(A)')
     & '<DataArray Name="x" type="Float32" format="ascii">'
      do i=itmp,iu-1
            write(1002,fmt_r) x0 + i*dxl
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="y" type="Float32" format="ascii">'
      do j=jtmp,ju-1
            write(1002,fmt_r) y0 + j*dyl
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="z" type="Float32" format="ascii">'
      do k=ktmp,ku-1
            write(1002,fmt_r) z0 + k*dzl
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</Coordinates>'

      write(1002,'(A)')'<PointData Scalars="'//trim(fnamebase)//'">'
      write(1002,'(A)')'<DataArray type="Float32" Name="'//
     & trim(fnamebase)//'" format="ascii">'
      itmp=1-1
      jtmp=1-1
      ktmp=1-1
      if(pidx.ge.1)itmp=1-2
      if(pidy.ge.1)jtmp=1-2
      if(pidz.ge.1)ktmp=1-2
      do k=ktmp,kml-1
      do j=jtmp,jml-1
      do i=itmp,iml-1
            write(1002,fmt_r)a(i+1,j+1,k+1)
      enddo
      enddo
      enddo

      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</PointData>'
      write(1002,'(A)')'</Piece>'
      write(1002,'(A)')'</RectilinearGrid>'
      write(1002,'(A)')'</VTKFile>'
      close(1002)

      end subroutine output_r8_pvtr_with_grid

      !-----------------------------------------------------------------

      subroutine output_r8_pvtr_2d(a,ng,fnamebase,timestep,timeval,gpid,
     & description)
      ! Output of pvtr requires at least one ghost layer (ng>=1)
      use boundary_conditions
      implicit none
      real(8),intent(in) :: a(1-ng:pv_iml+ng,1:1,1-ng:pv_kml+ng)
      integer(4),intent(in) :: ng
      character(len=*),intent(in) :: fnamebase
      integer(4),intent(in) :: timestep
      real(8),intent(in) :: timeval
      integer(4),intent(in) :: gpid
      character(len=*),intent(in) :: description

      integer(4) :: i,j,k,nn
      character(len=300) :: fnameg, fnamel, subdir
      character(len=8) :: ts2s, pid2s
      integer(4) :: pidx, pidy, pidz
      integer(4) :: il, iu, jl, ju, kl, ku
      integer(4) :: itmp, jtmp, ktmp
      character(len=300) :: str_extent_g, str_extent_l
      integer(4) :: ierr

      call get_pidxyz(pidx, pidy, pidz, gpid)
      call get_ijkrange(il, iu, jl, ju, kl, ku, gpid)

      subdir = trim(fnamebase)//'/'

      ! Make a subdir
      if(gpid.eq.0)then
        call execute_command_line(
     &                         'mkdir -p '//trim(PVOPATH)//trim(subdir))
      endif
      call mpi_barrier(MPI_COMM_WORLD, ierr)

      write(ts2s,"(i8.8)") int(timestep)
      write(pid2s,"(i4.4)") int(gpid)
      write(str_extent_g,"(6(i7,1x))")
     &                            0,pv_img-1,0,pv_jmg-1,0,pv_kmg-1
      fnameg = trim(PVOPATH)//trim(subdir)//trim(fnamebase)//'_'//
     & trim(ts2s)//'.pvtr'

      ! Output of pvtr file (meta file to combine each local data(vtr))
      if(gpid.eq.0)then
        call output_pvd(PVOPATH, fnamebase,subdir,timestep,timeval,
     &                                                      description)
        open(1001,FILE=trim(fnameg),status='unknown',action='write')
        write(1001,'(A)')'<?xml version="1.0"?>'
        write(1001,'(A)')'<VTKFile type="PRectilinearGrid" '//
     &   'version="0.1" byte_order="LittleEndian">'
        write(1001,'(A)')'<PRectilinearGrid WholeExtent="'//
     &   trim(adjustl(str_extent_g))//'">'

        write(1001,'(A)')'<PCoordinates>'
        write(1001,'(A)')
     &   '<DataArray Name="x" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="y" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="z" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PCoordinates>'

        write(1001,'(A)')'<PPointData Scalars="'//
     &   trim(fnamebase)//'">'
        write(1001,'(A)')'<DataArray type="Float32" Name="'//
     &   trim(fnamebase)//'" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PPointData>'

        do nn=0,nprocs_o-1
          write(pid2s,"(i4.4)") int(nn)
          fnamel=trim(fnamebase)//'_'//trim(ts2s)//'_'//
     &     trim(pid2s)//'.vtr'
          call get_ijkrange(il, iu, jl, ju, kl, ku, nn)
          call get_pidxyz(pidx, pidy, pidz, nn)

          itmp=il-1
          jtmp=jl-1
          ktmp=kl-1
          if(pidx.ge.1)itmp=il-2
          if(pidy.ge.1)jtmp=jl-2
          if(pidz.ge.1)ktmp=kl-2
          write(str_extent_l,"(6(i7,1x))")
     &      itmp,iu-1,jtmp,ju-1,ktmp,ku-1

          write(1001,'(A)')
     &     '<Piece Extent="'//trim(adjustl(str_extent_l))//
     &                '" Source="'//trim(fnamel)//'">'
          write(1001,'(A)')'</Piece>'
        enddo !do nn=0,nprocs_o-1

        write(1001,'(A)')'</PRectilinearGrid>'
        write(1001,'(A)')'</VTKFile>'
        close(1001)      
      endif !if(gpid.eq.0)then


      ! Output of vtr file (actual local data)
      write(pid2s,"(i4.4)") int(gpid)
      call get_ijkrange(il, iu, jl, ju, kl, ku, gpid)
      call get_pidxyz(pidx, pidy, pidz, gpid)

      itmp=il-1
      jtmp=jl-1
      ktmp=kl-1
      if(pidx.ge.1)itmp=il-2
      if(pidy.ge.1)jtmp=jl-2
      if(pidz.ge.1)ktmp=kl-2
      write(str_extent_l,"(6(i7,1x))") itmp,iu-1,jtmp,ju-1,ktmp,ku-1
      fnamel=trim(PVOPATH)//trim(subdir)//trim(fnamebase)//
     & '_'//trim(ts2s)//'_'//trim(pid2s)//'.vtr'
      open(1002,FILE=trim(fnamel),status='unknown',action='write')
      write(1002,'(A)')'<?xml version="1.0"?>'
      write(1002,'(A)')'<VTKFile type="RectilinearGrid" '//
     & 'version="0.1" byte_order="LittleEndian">'
      write(1002,'(A)')'<RectilinearGrid WholeExtent="'//
     & trim(adjustl(str_extent_g))//'">'
      write(1002,'(A)')'<Piece Extent="'//
     & trim(adjustl(str_extent_l))//'">'
      write(1002,'(A)')'<Coordinates>'
      write(1002,'(A)')
     & '<DataArray Name="x" type="Float32" format="ascii">'
      do i=itmp,iu-1
            write(1002,fmt_r) x0_o + i*dx_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="y" type="Float32" format="ascii">'
      do j=jtmp,ju-1
            write(1002,fmt_r) y0_o + j*dy_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="z" type="Float32" format="ascii">'
      do k=ktmp,ku-1
            write(1002,fmt_r) z0_o + k*dz_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</Coordinates>'

      write(1002,'(A)')'<PointData Scalars="'//trim(fnamebase)//'">'
      write(1002,'(A)')'<DataArray type="Float32" Name="'//
     & trim(fnamebase)//'" format="ascii">'
      itmp=1-1
      jtmp=1-1
      ktmp=1-1
      if(pidx.ge.1)itmp=1-2
      if(pidy.ge.1)jtmp=1-2
      if(pidz.ge.1)ktmp=1-2
      do k=ktmp,pv_kml-1
      do j=1,1
      do i=itmp,pv_iml-1
            write(1002,fmt_r)a(i+1,j,k+1)
      enddo
      enddo
      enddo

      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</PointData>'
      write(1002,'(A)')'</Piece>'
      write(1002,'(A)')'</RectilinearGrid>'
      write(1002,'(A)')'</VTKFile>'
      close(1002)

      end subroutine output_r8_pvtr_2d

      !-----------------------------------------------------------------


      subroutine output_r8_pvtr_2d_with_grid(gridobj,iml,jml,kml,ng,a,
     &                      fnamebase,timestep,timeval,gpid,description)
      ! Output of pvtr requires at least one ghost layer (ng>=1)

      implicit none
      type(grid),intent(in) :: gridobj
      integer(4),intent(in) :: iml,jml,kml
      integer(4),intent(in) :: ng
      real(8),intent(in) :: a(1-ng:iml+ng,1:1,1-ng:kml+ng)
      character(len=*),intent(in) :: fnamebase
      integer(4),intent(in) :: timestep
      real(8),intent(in) :: timeval
      integer(4),intent(in) :: gpid
      character(len=*),intent(in) :: description

      integer(4) :: i,j,k,nn
      character(len=300) :: fnameg, fnamel, subdir
      character(len=8) :: ts2s, pid2s
      integer(4) :: pidx, pidy, pidz
      integer(4) :: il, iu, jl, ju, kl, ku
      integer(4) :: itmp, jtmp, ktmp
      character(len=300) :: str_extent_g, str_extent_l
      integer(4) :: ierr
      integer(4) :: img, jmg, kmg
      real(8) :: x0, y0, z0
      real(8) :: dxl, dyl, dzl 

      if(gridobj%is_initialized .eqv. .false.)then
        write(*,*)"Error, M801 pvoutputs, output_i4_pvtr_with_grid()"
        write(*,*)"  Uninitialized grid object delivered."
        write(*,*)"  Aborting..."
        call abort()
      endif

      x0 = gridobj%x0 
      y0 = gridobj%y0 
      z0 = gridobj%z0 
      dxl = gridobj%dxl 
      dyl = gridobj%dyl 
      dzl = gridobj%dzl 

      img = gridobj%img 
      jmg = gridobj%jmg 
      kmg = gridobj%kmg 

      pidx = gridobj%pidx
      pidy = gridobj%pidy
      pidz = gridobj%pidz
      il = gridobj%il 
      iu = gridobj%iu 
      jl = gridobj%jl 
      ju = gridobj%ju 
      kl = gridobj%kl 
      ku = gridobj%ku 

      subdir = trim(fnamebase)//'/'

      ! Make a subdir
      if(gpid.eq.0)then
        call execute_command_line(
     &                         'mkdir -p '//trim(PVOPATH)//trim(subdir))
      endif
      call mpi_barrier(MPI_COMM_WORLD, ierr)

      write(ts2s,"(i8.8)") int(timestep)
      write(pid2s,"(i4.4)") int(gpid)
      write(str_extent_g,"(6(i7,1x))")
     &                            0,img-1,0,jmg-1,0,kmg-1
      fnameg = trim(PVOPATH)//trim(subdir)//trim(fnamebase)//'_'//
     & trim(ts2s)//'.pvtr'

      ! Output of pvtr file (meta file to combine each local data(vtr))
      if(gpid.eq.0)then
        call output_pvd(PVOPATH, fnamebase,subdir,timestep,timeval,
     &                                                      description)
        open(1001,FILE=trim(fnameg),status='unknown',action='write')
        write(1001,'(A)')'<?xml version="1.0"?>'
        write(1001,'(A)')'<VTKFile type="PRectilinearGrid" '//
     &   'version="0.1" byte_order="LittleEndian">'
        write(1001,'(A)')'<PRectilinearGrid WholeExtent="'//
     &   trim(adjustl(str_extent_g))//'">'

        write(1001,'(A)')'<PCoordinates>'
        write(1001,'(A)')
     &   '<DataArray Name="x" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="y" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="z" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PCoordinates>'

        write(1001,'(A)')'<PPointData Scalars="'//
     &   trim(fnamebase)//'">'
        write(1001,'(A)')'<DataArray type="Float32" Name="'//
     &   trim(fnamebase)//'" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PPointData>'

        do nn=0,nprocs_o-1
          write(pid2s,"(i4.4)") int(nn)
          fnamel=trim(fnamebase)//'_'//trim(ts2s)//'_'//
     &     trim(pid2s)//'.vtr'
          call grid_get_ijkrange(il, iu, jl, ju, kl, ku, gridobj, nn)
          call grid_get_pidxyz(pidx, pidy, pidz, nn, gridobj)

          itmp=il-1
          jtmp=jl-1
          ktmp=kl-1
          if(pidx.ge.1)itmp=il-2
          if(pidy.ge.1)jtmp=jl-2
          if(pidz.ge.1)ktmp=kl-2
          write(str_extent_l,"(6(i7,1x))")
     &      itmp,iu-1,jtmp,ju-1,ktmp,ku-1

          write(1001,'(A)')
     &     '<Piece Extent="'//trim(adjustl(str_extent_l))//
     &                '" Source="'//trim(fnamel)//'">'
          write(1001,'(A)')'</Piece>'
        enddo !do nn=0,nprocs_o-1

        write(1001,'(A)')'</PRectilinearGrid>'
        write(1001,'(A)')'</VTKFile>'
        close(1001)      
      endif !if(gpid.eq.0)then


      ! Output of vtr file (actual local data)
      write(pid2s,"(i4.4)") int(gpid)
      call grid_get_ijkrange(il, iu, jl, ju, kl, ku, gridobj, gpid)
      call grid_get_pidxyz(pidx, pidy, pidz, gpid, gridobj)

      itmp=il-1
      jtmp=jl-1
      ktmp=kl-1
      if(pidx.ge.1)itmp=il-2
      if(pidy.ge.1)jtmp=jl-2
      if(pidz.ge.1)ktmp=kl-2
      write(str_extent_l,"(6(i7,1x))") itmp,iu-1,jtmp,ju-1,ktmp,ku-1
      fnamel=trim(PVOPATH)//trim(subdir)//trim(fnamebase)//
     & '_'//trim(ts2s)//'_'//trim(pid2s)//'.vtr'
      open(1002,FILE=trim(fnamel),status='unknown',action='write')
      write(1002,'(A)')'<?xml version="1.0"?>'
      write(1002,'(A)')'<VTKFile type="RectilinearGrid" '//
     & 'version="0.1" byte_order="LittleEndian">'
      write(1002,'(A)')'<RectilinearGrid WholeExtent="'//
     & trim(adjustl(str_extent_g))//'">'
      write(1002,'(A)')'<Piece Extent="'//
     & trim(adjustl(str_extent_l))//'">'
      write(1002,'(A)')'<Coordinates>'
      write(1002,'(A)')
     & '<DataArray Name="x" type="Float32" format="ascii">'
      do i=itmp,iu-1
            write(1002,fmt_r) x0 + i*dxl
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="y" type="Float32" format="ascii">'
      do j=jtmp,ju-1
            write(1002,fmt_r) y0 + j*dyl
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="z" type="Float32" format="ascii">'
      do k=ktmp,ku-1
            write(1002,fmt_r) z0 + k*dzl
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</Coordinates>'

      write(1002,'(A)')'<PointData Scalars="'//trim(fnamebase)//'">'
      write(1002,'(A)')'<DataArray type="Float32" Name="'//
     & trim(fnamebase)//'" format="ascii">'
      itmp=1-1
      jtmp=1-1
      ktmp=1-1
      if(pidx.ge.1)itmp=1-2
      if(pidy.ge.1)jtmp=1-2
      if(pidz.ge.1)ktmp=1-2
      do k=ktmp,gridobj%kml-1
      do j=1,1
      do i=itmp,gridobj%iml-1
            write(1002,fmt_r)a(i+1,j,k+1)
      enddo
      enddo
      enddo

      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</PointData>'
      write(1002,'(A)')'</Piece>'
      write(1002,'(A)')'</RectilinearGrid>'
      write(1002,'(A)')'</VTKFile>'
      close(1002)

      end subroutine output_r8_pvtr_2d_with_grid

      ! ----------------------------------------------------------------

      subroutine output_r8_sub_pvtr_2d_with_grid(gridobj,iml,jml,kml,ng,
     &                    imin,imax,jmin,jmax,kmin,kmax,
     &                    a,fnamebase,timestep,timeval,gpid,description)
      ! Output of pvtr requires at least one ghost layer (ng>=1)

      implicit none
      type(grid),intent(in) :: gridobj
      integer(4),intent(in) :: iml,jml,kml
      integer(4),intent(in) :: ng
      integer(4),intent(in) :: imin,imax,jmin,jmax,kmin,kmax
      real(8),intent(in) :: a(1-ng:iml+ng,1:1,1-ng:kml+ng)
      character(len=*),intent(in) :: fnamebase
      integer(4),intent(in) :: timestep
      real(8),intent(in) :: timeval
      integer(4),intent(in) :: gpid
      character(len=*),intent(in) :: description

      integer(4) :: i,j,k,nn
      character(len=300) :: fnameg, fnamel, subdir
      character(len=8) :: ts2s, pid2s
      integer(4) :: pidx, pidy, pidz
      integer(4) :: il, iu, jl, ju, kl, ku
      integer(4) :: itmp, jtmp, ktmp
      character(len=300) :: str_extent_g, str_extent_l
      integer(4) :: ierr
      integer(4) :: img, jmg, kmg
      real(8) :: x0, y0, z0
      real(8) :: dxl, dyl, dzl 
      integer(4) :: ilg, iug, jlg, jug, klg, kug

      if(gridobj%is_initialized .eqv. .false.)then
        write(*,*)"Error, M801 pvoutputs, output_r8_pvtr_with_grid()"
        write(*,*)"  Uninitialized grid object delivered."
        write(*,*)"  Aborting..."
        call abort()
      endif

      x0 = gridobj%x0 
      y0 = gridobj%y0 
      z0 = gridobj%z0 
      dxl = gridobj%dxl 
      dyl = gridobj%dyl 
      dzl = gridobj%dzl 

      ! img = gridobj%img 
      ! jmg = gridobj%jmg 
      ! kmg = gridobj%kmg 

      pidx = gridobj%pidx
      pidy = gridobj%pidy
      pidz = gridobj%pidz
      il = gridobj%il 
      iu = gridobj%iu 
      jl = gridobj%jl 
      ju = gridobj%ju 
      kl = gridobj%kl 
      ku = gridobj%ku 

      img = imax-imin+1
      jmg = jmax-jmin+1
      kmg = kmax-kmin+1

      if(gridobj%npx.ne.1.and.img.ne.gridobj%img .or. 
     &   gridobj%npy.ne.1.and.jmg.ne.gridobj%jmg .or.
     &   gridobj%npz.ne.1.and.kmg.ne.gridobj%kmg )then
        write(*,*)
     &           "Error, M801 pvoutputs, output_r8_sub,pvtr_with_grid()"
        write(*,*)"  Invalid NPX, NPY, or NPZ."
        write(*,*)"  To make output of the shrunken array (< *MORI), "
        write(*,*)"  NP* should be only 1 along the shrunken direction."
        write(*,*)"  Aborting..."
        call abort()
      endif

      ilg = 0
      iug = img-1
      jlg = 0
      jug = jmg-1
      klg = 0
      kug = kmg-1

      if(img.lt.gridobj%img) then
        ilg = imin-1
        iug = imax-1
      endif

      if(jmg.lt.gridobj%jmg) then
        jlg = jmin-1
        jug = jmax-1
      endif

      if(kmg.lt.gridobj%kmg) then
        klg = kmin-1
        kug = kmax-1
      endif

      subdir = trim(fnamebase)//'/'

      ! Make a subdir
      if(gpid.eq.0)then
        call execute_command_line(
     &                         'mkdir -p '//trim(PVOPATH)//trim(subdir))
      endif
      call mpi_barrier(MPI_COMM_WORLD, ierr)

      write(ts2s,"(i8.8)") int(timestep)
      write(pid2s,"(i4.4)") int(gpid)
      write(str_extent_g,"(6(i7,1x))") ilg,iug,jlg,jug,klg,kug
      fnameg = trim(PVOPATH)//trim(subdir)//trim(fnamebase)//'_'//
     & trim(ts2s)//'.pvtr'

      ! Output of pvtr file (meta file to combine each local data(vtr))
      if(gpid.eq.0)then
        call output_pvd(PVOPATH, fnamebase,subdir,timestep,timeval,
     &                                                      description)
        open(1001,FILE=trim(fnameg),status='unknown',action='write')
        write(1001,'(A)')'<?xml version="1.0"?>'
        write(1001,'(A)')'<VTKFile type="PRectilinearGrid" '//
     &   'version="0.1" byte_order="LittleEndian">'
        write(1001,'(A)')'<PRectilinearGrid WholeExtent="'//
     &   trim(adjustl(str_extent_g))//'">'

        write(1001,'(A)')'<PCoordinates>'
        write(1001,'(A)')
     &   '<DataArray Name="x" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="y" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="z" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PCoordinates>'

        write(1001,'(A)')'<PPointData Scalars="'//
     &   trim(fnamebase)//'">'
        write(1001,'(A)')'<DataArray type="Float32" Name="'//
     &   trim(fnamebase)//'" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PPointData>'

        do nn=0,nprocs_o-1
          write(pid2s,"(i4.4)") int(nn)
          fnamel=trim(fnamebase)//'_'//trim(ts2s)//'_'//
     &     trim(pid2s)//'.vtr'
          call grid_get_ijkrange(il, iu, jl, ju, kl, ku, gridobj, nn)
          call grid_get_pidxyz(pidx, pidy, pidz, nn, gridobj)

          if(img.lt.gridobj%img) then
            il = imin
            iu = imax 
          endif

          if(jmg.lt.gridobj%jmg) then
            jl = jmin
            ju = jmax 
          endif

          if(kmg.lt.gridobj%kmg) then
            kl = kmin
            ku = kmax 
          endif

          itmp=il-1
          jtmp=jl-1
          ktmp=kl-1
          if(pidx.ge.1)itmp=il-2
          if(pidy.ge.1)jtmp=jl-2
          if(pidz.ge.1)ktmp=kl-2
          write(str_extent_l,"(6(i7,1x))")
     &      itmp,iu-1,jtmp,ju-1,ktmp,ku-1

          write(1001,'(A)')
     &     '<Piece Extent="'//trim(adjustl(str_extent_l))//
     &                '" Source="'//trim(fnamel)//'">'
          write(1001,'(A)')'</Piece>'
        enddo !do nn=0,nprocs_o-1

        write(1001,'(A)')'</PRectilinearGrid>'
        write(1001,'(A)')'</VTKFile>'
        close(1001)      
      endif !if(gpid.eq.0)then


      ! Output of vtr file (actual local data)
      write(pid2s,"(i4.4)") int(gpid)
      call grid_get_ijkrange(il, iu, jl, ju, kl, ku, gridobj, gpid)
      call grid_get_pidxyz(pidx, pidy, pidz, gpid, gridobj)

      if(img.lt.gridobj%img) then
        il = imin
        iu = imax 
      endif

      if(jmg.lt.gridobj%jmg) then
        jl = jmin
        ju = jmax 
      endif

      if(kmg.lt.gridobj%kmg) then
        kl = kmin
        ku = kmax 
      endif

      itmp=il-1
      jtmp=jl-1
      ktmp=kl-1
      if(pidx.ge.1)itmp=il-2
      if(pidy.ge.1)jtmp=jl-2
      if(pidz.ge.1)ktmp=kl-2
      write(str_extent_l,"(6(i7,1x))") itmp,iu-1,jtmp,ju-1,ktmp,ku-1
      fnamel=trim(PVOPATH)//trim(subdir)//trim(fnamebase)//
     & '_'//trim(ts2s)//'_'//trim(pid2s)//'.vtr'
      open(1002,FILE=trim(fnamel),status='unknown',action='write')
      write(1002,'(A)')'<?xml version="1.0"?>'
      write(1002,'(A)')'<VTKFile type="RectilinearGrid" '//
     & 'version="0.1" byte_order="LittleEndian">'
      write(1002,'(A)')'<RectilinearGrid WholeExtent="'//
     & trim(adjustl(str_extent_g))//'">'
      write(1002,'(A)')'<Piece Extent="'//
     & trim(adjustl(str_extent_l))//'">'
      write(1002,'(A)')'<Coordinates>'
      write(1002,'(A)')
     & '<DataArray Name="x" type="Float32" format="ascii">'
      do i=itmp,iu-1
            write(1002,fmt_r) x0 + i*dxl
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="y" type="Float32" format="ascii">'
      do j=jtmp,ju-1
            write(1002,fmt_r) y0 + j*dyl
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="z" type="Float32" format="ascii">'
      do k=ktmp,ku-1
            write(1002,fmt_r) z0 + k*dzl
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</Coordinates>'

      write(1002,'(A)')'<PointData Scalars="'//trim(fnamebase)//'">'
      write(1002,'(A)')'<DataArray type="Float32" Name="'//
     & trim(fnamebase)//'" format="ascii">'
      itmp=1-1
      jtmp=1-1
      ktmp=1-1
      if(pidx.ge.1)itmp=1-2
      if(pidy.ge.1)jtmp=1-2
      if(pidz.ge.1)ktmp=1-2

      il=itmp 
      iu=iml-1
      jl=1
      ju=1
      kl=ktmp
      ku=kml-1
      if(img.lt.gridobj%img) then
        il = imin-1
        iu = imax-1
      endif
      if(jmg.lt.gridobj%jmg) then
        jl = jmin-1
        ju = jmax-1
      endif
      if(kmg.lt.gridobj%kmg) then
        kl = kmin-1
        ku = kmax-1
      endif

      do k=kl,ku
      do j=1,1
      do i=il,iu
            write(1002,fmt_r)a(i+1,j,k+1)
      enddo
      enddo
      enddo

      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</PointData>'
      write(1002,'(A)')'</Piece>'
      write(1002,'(A)')'</RectilinearGrid>'
      write(1002,'(A)')'</VTKFile>'
      close(1002)

      end subroutine output_r8_sub_pvtr_2d_with_grid

      ! ----------------------------------------------------------------

      subroutine checkpoint_i4_pvtr(a,ng,fnamebase,timestep,timeval,
     & gpid, description)
      ! Output of pvtr requires at least one ghost layer (ng>=1)

      implicit none
      integer(4),intent(in) :: 
     &                   a(1-ng:pv_iml+ng,1-ng:pv_jml+ng,1-ng:pv_kml+ng)
      integer(4),intent(in) :: ng
      character(len=*),intent(in) :: fnamebase
      integer(4),intent(in) :: timestep
      real(8),intent(in) :: timeval
      integer(4),intent(in) :: gpid
      character(len=*),intent(in) :: description

      integer(4) :: i,j,k,nn
      character(len=300) :: fnameg, fnamel, subdir
      character(len=8) :: ts2s, pid2s
      integer(4) :: pidx, pidy, pidz
      integer(4) :: il, iu, jl, ju, kl, ku
      integer(4) :: itmp, jtmp, ktmp
      character(len=300) :: str_extent_g, str_extent_l
      integer(4) :: ierr

      call get_pidxyz(pidx, pidy, pidz, gpid)
      call get_ijkrange(il, iu, jl, ju, kl, ku, gpid)

      subdir = trim(fnamebase)//'/'

      ! Make a subdir
      if(gpid.eq.0)then
      call execute_command_line(
     & 'mkdir -p '//trim(CHECKPOINTPATH))
      call execute_command_line(
     & 'mkdir -p '//trim(CHECKPOINTPATH)//trim(subdir))
      endif
      call mpi_barrier(MPI_COMM_WORLD, ierr)

      write(ts2s,"(i8.8)") int(timestep)
      write(pid2s,"(i4.4)") int(gpid)
      write(str_extent_g,"(6(i7,1x))")
     &                            0,pv_img-1,0,pv_jmg-1,0,pv_kmg-1
      fnameg = trim(CHECKPOINTPATH)//trim(subdir)//
     & trim(fnamebase)//'_'//trim(ts2s)//'.pvtr'

      ! Output of pvtr file (meta file to combine each local data(vtr))
      if(gpid.eq.0)then
        call output_pvd(CHECKPOINTPATH,fnamebase,subdir,timestep,
     &                                              timeval,description)
        open(1001,FILE=trim(fnameg),status='unknown',action='write')
        write(1001,'(A)')'<?xml version="1.0"?>'
        write(1001,'(A)')'<VTKFile type="PRectilinearGrid" '//
     &   'version="0.1" byte_order="LittleEndian">'
        write(1001,'(A)')'<PRectilinearGrid WholeExtent="'//
     &   trim(adjustl(str_extent_g))//'">'

        write(1001,'(A)')'<PCoordinates>'
        write(1001,'(A)')
     &   '<DataArray Name="x" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="y" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="z" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PCoordinates>'

        write(1001,'(A)')'<PPointData Scalars="'//
     &   trim(fnamebase)//'">'
        write(1001,'(A)')'<DataArray type="Int32" Name="'//
     &   trim(fnamebase)//'" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PPointData>'

        do nn=0,nprocs_o-1
          write(pid2s,"(i4.4)") int(nn)
          fnamel=trim(fnamebase)//'_'//trim(ts2s)//'_'//
     &     trim(pid2s)//'.vtr'
          call get_ijkrange(il, iu, jl, ju, kl, ku, nn)
          call get_pidxyz(pidx, pidy, pidz, nn)

          itmp=il-1
          jtmp=jl-1
          ktmp=kl-1
          if(pidx.ge.1)itmp=il-2
          if(pidy.ge.1)jtmp=jl-2
          if(pidz.ge.1)ktmp=kl-2
          write(str_extent_l,"(6(i7,1x))")
     &      itmp,iu-1,jtmp,ju-1,ktmp,ku-1

          write(1001,'(A)')
     &     '<Piece Extent="'//trim(adjustl(str_extent_l))//
     &                '" Source="'//trim(fnamel)//'">'
          write(1001,'(A)')'</Piece>'
        enddo !do nn=0,nprocs_o-1

        write(1001,'(A)')'</PRectilinearGrid>'
        write(1001,'(A)')'</VTKFile>'
        close(1001)      
      endif !if(gpid.eq.0)then


      ! Output of vtr file (actual local data)
      write(pid2s,"(i4.4)") int(gpid)
      call get_ijkrange(il, iu, jl, ju, kl, ku, gpid)
      call get_pidxyz(pidx, pidy, pidz, gpid)

      itmp=il-1
      jtmp=jl-1
      ktmp=kl-1
      if(pidx.ge.1)itmp=il-2
      if(pidy.ge.1)jtmp=jl-2
      if(pidz.ge.1)ktmp=kl-2
      write(str_extent_l,"(6(i7,1x))") itmp,iu-1,jtmp,ju-1,ktmp,ku-1
      fnamel=trim(CHECKPOINTPATH)//trim(subdir)//
     & trim(fnamebase)//'_'//trim(ts2s)//'_'//trim(pid2s)//'.vtr'
      open(1002,FILE=trim(fnamel),status='unknown',action='write')
      write(1002,'(A)')'<?xml version="1.0"?>'
      write(1002,'(A)')'<VTKFile type="RectilinearGrid" '//
     & 'version="0.1" byte_order="LittleEndian">'
      write(1002,'(A)')'<RectilinearGrid WholeExtent="'//
     & trim(adjustl(str_extent_g))//'">'
      write(1002,'(A)')'<Piece Extent="'//
     & trim(adjustl(str_extent_l))//'">'
      write(1002,'(A)')'<Coordinates>'
      write(1002,'(A)')
     & '<DataArray Name="x" type="Float32" format="ascii">'
      do i=itmp,iu-1
            write(1002,fmt_r) x0_o + i*dx_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="y" type="Float32" format="ascii">'
      do j=jtmp,ju-1
            write(1002,fmt_r) y0_o + j*dy_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="z" type="Float32" format="ascii">'
      do k=ktmp,ku-1
            write(1002,fmt_r) z0_o + k*dz_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</Coordinates>'

      write(1002,'(A)')'<PointData Scalars="'//trim(fnamebase)//'">'
      write(1002,'(A)')'<DataArray type="Int32" Name="'//
     & trim(fnamebase)//'" format="ascii">'
      itmp=1-1
      jtmp=1-1
      ktmp=1-1
      if(pidx.ge.1)itmp=1-2
      if(pidy.ge.1)jtmp=1-2
      if(pidz.ge.1)ktmp=1-2
      do k=ktmp,pv_kml-1
      do j=jtmp,pv_jml-1
      do i=itmp,pv_iml-1
            write(1002,'(I0)')a(i+1,j+1,k+1)
      enddo
      enddo
      enddo

      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</PointData>'
      write(1002,'(A)')'</Piece>'
      write(1002,'(A)')'</RectilinearGrid>'
      write(1002,'(A)')'</VTKFile>'
      close(1002)

      end subroutine checkpoint_i4_pvtr

      !-----------------------------------------------------------------

      subroutine checkpoint_i4_pvtr_2d(a,ng,fnamebase,timestep,timeval,
     & gpid, description)
      ! Output of pvtr requires at least one ghost layer (ng>=1)

      implicit none
      integer(4),intent(in) :: a(1-ng:pv_iml+ng,1:1,1-ng:pv_kml+ng)
      integer(4),intent(in) :: ng
      character(len=*),intent(in) :: fnamebase
      integer(4),intent(in) :: timestep
      real(8),intent(in) :: timeval
      integer(4),intent(in) :: gpid
      character(len=*),intent(in) :: description

      integer(4) :: i,j,k,nn
      character(len=300) :: fnameg, fnamel, subdir
      character(len=8) :: ts2s, pid2s
      integer(4) :: pidx, pidy, pidz
      integer(4) :: il, iu, jl, ju, kl, ku
      integer(4) :: itmp, jtmp, ktmp
      character(len=300) :: str_extent_g, str_extent_l
      integer(4) :: ierr

      call get_pidxyz(pidx, pidy, pidz, gpid)
      call get_ijkrange(il, iu, jl, ju, kl, ku, gpid)

      subdir = trim(fnamebase)//'/'

      ! Make a subdir
      if(gpid.eq.0)then
      call execute_command_line(
     & 'mkdir -p '//trim(CHECKPOINTPATH))
      call execute_command_line(
     & 'mkdir -p '//trim(CHECKPOINTPATH)//trim(subdir))
      endif
      call mpi_barrier(MPI_COMM_WORLD, ierr)

      write(ts2s,"(i8.8)") int(timestep)
      write(pid2s,"(i4.4)") int(gpid)
      write(str_extent_g,"(6(i7,1x))")
     &                            0,pv_img-1,0,pv_jmg-1,0,pv_kmg-1
      fnameg = trim(CHECKPOINTPATH)//trim(subdir)//
     & trim(fnamebase)//'_'//trim(ts2s)//'.pvtr'

      ! Output of pvtr file (meta file to combine each local data(vtr))
      if(gpid.eq.0)then
        call output_pvd(CHECKPOINTPATH,fnamebase,subdir,timestep,
     &                                              timeval,description)
        open(1001,FILE=trim(fnameg),status='unknown',action='write')
        write(1001,'(A)')'<?xml version="1.0"?>'
        write(1001,'(A)')'<VTKFile type="PRectilinearGrid" '//
     &   'version="0.1" byte_order="LittleEndian">'
        write(1001,'(A)')'<PRectilinearGrid WholeExtent="'//
     &   trim(adjustl(str_extent_g))//'">'

        write(1001,'(A)')'<PCoordinates>'
        write(1001,'(A)')
     &   '<DataArray Name="x" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="y" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="z" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PCoordinates>'

        write(1001,'(A)')'<PPointData Scalars="'//
     &   trim(fnamebase)//'">'
        write(1001,'(A)')'<DataArray type="Int32" Name="'//
     &   trim(fnamebase)//'" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PPointData>'

        do nn=0,nprocs_o-1
          write(pid2s,"(i4.4)") int(nn)
          fnamel=trim(fnamebase)//'_'//trim(ts2s)//'_'//
     &     trim(pid2s)//'.vtr'
          call get_ijkrange(il, iu, jl, ju, kl, ku, nn)
          call get_pidxyz(pidx, pidy, pidz, nn)

          itmp=il-1
          jtmp=jl-1
          ktmp=kl-1
          if(pidx.ge.1)itmp=il-2
          if(pidy.ge.1)jtmp=jl-2
          if(pidz.ge.1)ktmp=kl-2
          write(str_extent_l,"(6(i7,1x))")
     &      itmp,iu-1,jtmp,ju-1,ktmp,ku-1

          write(1001,'(A)')
     &     '<Piece Extent="'//trim(adjustl(str_extent_l))//
     &                '" Source="'//trim(fnamel)//'">'
          write(1001,'(A)')'</Piece>'
        enddo !do nn=0,nprocs_o-1

        write(1001,'(A)')'</PRectilinearGrid>'
        write(1001,'(A)')'</VTKFile>'
        close(1001)      
      endif !if(gpid.eq.0)then


      ! Output of vtr file (actual local data)
      write(pid2s,"(i4.4)") int(gpid)
      call get_ijkrange(il, iu, jl, ju, kl, ku, gpid)
      call get_pidxyz(pidx, pidy, pidz, gpid)

      itmp=il-1
      jtmp=jl-1
      ktmp=kl-1
      if(pidx.ge.1)itmp=il-2
      if(pidy.ge.1)jtmp=jl-2
      if(pidz.ge.1)ktmp=kl-2
      write(str_extent_l,"(6(i7,1x))") itmp,iu-1,jtmp,ju-1,ktmp,ku-1
      fnamel=trim(CHECKPOINTPATH)//trim(subdir)//
     & trim(fnamebase)//'_'//trim(ts2s)//'_'//trim(pid2s)//'.vtr'
      open(1002,FILE=trim(fnamel),status='unknown',action='write')
      write(1002,'(A)')'<?xml version="1.0"?>'
      write(1002,'(A)')'<VTKFile type="RectilinearGrid" '//
     & 'version="0.1" byte_order="LittleEndian">'
      write(1002,'(A)')'<RectilinearGrid WholeExtent="'//
     & trim(adjustl(str_extent_g))//'">'
      write(1002,'(A)')'<Piece Extent="'//
     & trim(adjustl(str_extent_l))//'">'
      write(1002,'(A)')'<Coordinates>'
      write(1002,'(A)')
     & '<DataArray Name="x" type="Float32" format="ascii">'
      do i=itmp,iu-1
            write(1002,fmt_r) x0_o + i*dx_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="y" type="Float32" format="ascii">'
      do j=jtmp,ju-1
            write(1002,fmt_r) y0_o + j*dy_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="z" type="Float32" format="ascii">'
      do k=ktmp,ku-1
            write(1002,fmt_r) z0_o + k*dz_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</Coordinates>'

      write(1002,'(A)')'<PointData Scalars="'//trim(fnamebase)//'">'
      write(1002,'(A)')'<DataArray type="Int32" Name="'//
     & trim(fnamebase)//'" format="ascii">'
      itmp=1-1
      jtmp=1-1
      ktmp=1-1
      if(pidx.ge.1)itmp=1-2
      if(pidy.ge.1)jtmp=1-2
      if(pidz.ge.1)ktmp=1-2
      do k=ktmp,pv_kml-1
      do j=1,1
      do i=itmp,pv_iml-1
            write(1002,'(I0)')a(i+1,j,k+1)
      enddo
      enddo
      enddo

      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</PointData>'
      write(1002,'(A)')'</Piece>'
      write(1002,'(A)')'</RectilinearGrid>'
      write(1002,'(A)')'</VTKFile>'
      close(1002)

      end subroutine checkpoint_i4_pvtr_2d

      !-----------------------------------------------------------------

      subroutine checkpoint_r8_pvtr(a,ng,fnamebase,timestep,timeval,
     & gpid, description)
      ! Output of pvtr requires at least one ghost layer (ng>=1)

      implicit none
      real(8),intent(in) :: 
     &                   a(1-ng:pv_iml+ng,1-ng:pv_jml+ng,1-ng:pv_kml+ng)
      integer(4),intent(in) :: ng
      character(len=*),intent(in) :: fnamebase
      integer(4),intent(in) :: timestep
      real(8),intent(in) :: timeval
      integer(4),intent(in) :: gpid
      character(len=*),intent(in) :: description

      integer(4) :: i,j,k,nn
      character(len=300) :: fnameg, fnamel, subdir
      character(len=8) :: ts2s, pid2s
      integer(4) :: pidx, pidy, pidz
      integer(4) :: il, iu, jl, ju, kl, ku
      integer(4) :: itmp, jtmp, ktmp
      character(len=300) :: str_extent_g, str_extent_l
      integer(4) :: ierr

      call get_pidxyz(pidx, pidy, pidz, gpid)
      call get_ijkrange(il, iu, jl, ju, kl, ku, gpid)

      subdir = trim(fnamebase)//'/'

      ! Make a subdir
      if(gpid.eq.0)then
      call execute_command_line(
     & 'mkdir -p '//trim(CHECKPOINTPATH))
      call execute_command_line(
     & 'mkdir -p '//trim(CHECKPOINTPATH)//trim(subdir))
      endif
      call mpi_barrier(MPI_COMM_WORLD, ierr)

      write(ts2s,"(i8.8)") int(timestep)
      write(pid2s,"(i4.4)") int(gpid)
      write(str_extent_g,"(6(i7,1x))")
     &                            0,pv_img-1,0,pv_jmg-1,0,pv_kmg-1
      fnameg = trim(CHECKPOINTPATH)//trim(subdir)//
     & trim(fnamebase)//'_'//trim(ts2s)//'.pvtr'

      ! Output of pvtr file (meta file to combine each local data(vtr))
      if(gpid.eq.0)then
        call output_pvd(CHECKPOINTPATH,fnamebase,subdir,timestep,
     &                                              timeval,description)
        open(1001,FILE=trim(fnameg),status='unknown',action='write')
        write(1001,'(A)')'<?xml version="1.0"?>'
        write(1001,'(A)')'<VTKFile type="PRectilinearGrid" '//
     &   'version="0.1" byte_order="LittleEndian">'
        write(1001,'(A)')'<PRectilinearGrid WholeExtent="'//
     &   trim(adjustl(str_extent_g))//'">'

        write(1001,'(A)')'<PCoordinates>'
        write(1001,'(A)')
     &   '<DataArray Name="x" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="y" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="z" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PCoordinates>'

        write(1001,'(A)')'<PPointData Scalars="'//
     &   trim(fnamebase)//'">'
        write(1001,'(A)')'<DataArray type="Float32" Name="'//
     &   trim(fnamebase)//'" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PPointData>'

        do nn=0,nprocs_o-1
          write(pid2s,"(i4.4)") int(nn)
          fnamel=trim(fnamebase)//'_'//trim(ts2s)//'_'//
     &     trim(pid2s)//'.vtr'
          call get_ijkrange(il, iu, jl, ju, kl, ku, nn)
          call get_pidxyz(pidx, pidy, pidz, nn)

          itmp=il-1
          jtmp=jl-1
          ktmp=kl-1
          if(pidx.ge.1)itmp=il-2
          if(pidy.ge.1)jtmp=jl-2
          if(pidz.ge.1)ktmp=kl-2
          write(str_extent_l,"(6(i7,1x))")
     &      itmp,iu-1,jtmp,ju-1,ktmp,ku-1

          write(1001,'(A)')
     &     '<Piece Extent="'//trim(adjustl(str_extent_l))//
     &                '" Source="'//trim(fnamel)//'">'
          write(1001,'(A)')'</Piece>'
        enddo !do nn=0,nprocs_o-1

        write(1001,'(A)')'</PRectilinearGrid>'
        write(1001,'(A)')'</VTKFile>'
        close(1001)      
      endif !if(gpid.eq.0)then


      ! Output of vtr file (actual local data)
      write(pid2s,"(i4.4)") int(gpid)
      call get_ijkrange(il, iu, jl, ju, kl, ku, gpid)
      call get_pidxyz(pidx, pidy, pidz, gpid)

      itmp=il-1
      jtmp=jl-1
      ktmp=kl-1
      if(pidx.ge.1)itmp=il-2
      if(pidy.ge.1)jtmp=jl-2
      if(pidz.ge.1)ktmp=kl-2
      write(str_extent_l,"(6(i7,1x))") itmp,iu-1,jtmp,ju-1,ktmp,ku-1
      fnamel=trim(CHECKPOINTPATH)//trim(subdir)//
     & trim(fnamebase)//'_'//trim(ts2s)//'_'//trim(pid2s)//'.vtr'
      open(1002,FILE=trim(fnamel),status='unknown',action='write')
      write(1002,'(A)')'<?xml version="1.0"?>'
      write(1002,'(A)')'<VTKFile type="RectilinearGrid" '//
     & 'version="0.1" byte_order="LittleEndian">'
      write(1002,'(A)')'<RectilinearGrid WholeExtent="'//
     & trim(adjustl(str_extent_g))//'">'
      write(1002,'(A)')'<Piece Extent="'//
     & trim(adjustl(str_extent_l))//'">'
      write(1002,'(A)')'<Coordinates>'
      write(1002,'(A)')
     & '<DataArray Name="x" type="Float32" format="ascii">'
      do i=itmp,iu-1
            write(1002,fmt_r) x0_o + i*dx_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="y" type="Float32" format="ascii">'
      do j=jtmp,ju-1
            write(1002,fmt_r) y0_o + j*dy_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="z" type="Float32" format="ascii">'
      do k=ktmp,ku-1
            write(1002,fmt_r) z0_o + k*dz_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</Coordinates>'

      write(1002,'(A)')'<PointData Scalars="'//trim(fnamebase)//'">'
      write(1002,'(A)')'<DataArray type="Float32" Name="'//
     & trim(fnamebase)//'" format="ascii">'
      itmp=1-1
      jtmp=1-1
      ktmp=1-1
      if(pidx.ge.1)itmp=1-2
      if(pidy.ge.1)jtmp=1-2
      if(pidz.ge.1)ktmp=1-2
      do k=ktmp,pv_kml-1
      do j=jtmp,pv_jml-1
      do i=itmp,pv_iml-1
            write(1002,fmt_r16)a(i+1,j+1,k+1)
      enddo
      enddo
      enddo

      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</PointData>'
      write(1002,'(A)')'</Piece>'
      write(1002,'(A)')'</RectilinearGrid>'
      write(1002,'(A)')'</VTKFile>'
      close(1002)

      end subroutine checkpoint_r8_pvtr

      !-----------------------------------------------------------------

      subroutine checkpoint_r8_pvtr_2d(a,ng,fnamebase,timestep,timeval, 
     & gpid,description)
      ! Output of pvtr requires at least one ghost layer (ng>=1)
      use boundary_conditions
      implicit none
      real(8),intent(in) :: a(1-ng:pv_iml+ng,1:1,1-ng:pv_kml+ng)
      integer(4),intent(in) :: ng
      character(len=*),intent(in) :: fnamebase
      integer(4),intent(in) :: timestep
      real(8),intent(in) :: timeval
      integer(4),intent(in) :: gpid
      character(len=*),intent(in) :: description

      integer(4) :: i,j,k,nn
      character(len=300) :: fnameg, fnamel, subdir
      character(len=8) :: ts2s, pid2s
      integer(4) :: pidx, pidy, pidz
      integer(4) :: il, iu, jl, ju, kl, ku
      integer(4) :: itmp, jtmp, ktmp
      character(len=300) :: str_extent_g, str_extent_l
      integer(4) :: ierr

      call get_pidxyz(pidx, pidy, pidz, gpid)
      call get_ijkrange(il, iu, jl, ju, kl, ku, gpid)

      subdir = trim(fnamebase)//'/'

      ! Make a subdir
      if(gpid.eq.0)then
      call execute_command_line(
     & 'mkdir -p '//trim(PVOPATH)//'checkpoint/')
      call execute_command_line(
     & 'mkdir -p '//trim(PVOPATH)//'checkpoint/'//trim(subdir))
      endif
      call mpi_barrier(MPI_COMM_WORLD, ierr)

      write(ts2s,"(i8.8)") int(timestep)
      write(pid2s,"(i4.4)") int(gpid)
      write(str_extent_g,"(6(i7,1x))")
     &                            0,pv_img-1,0,pv_jmg-1,0,pv_kmg-1
      fnameg = trim(PVOPATH)//'checkpoint/'//trim(subdir)//
     & trim(fnamebase)//'_'//trim(ts2s)//'.pvtr'

      ! Output of pvtr file (meta file to combine each local data(vtr))
      if(gpid.eq.0)then
        call output_pvd(CHECKPOINTPATH,fnamebase,subdir,timestep,
     &                                              timeval,description)
        open(1001,FILE=trim(fnameg),status='unknown',action='write')
        write(1001,'(A)')'<?xml version="1.0"?>'
        write(1001,'(A)')'<VTKFile type="PRectilinearGrid" '//
     &   'version="0.1" byte_order="LittleEndian">'
        write(1001,'(A)')'<PRectilinearGrid WholeExtent="'//
     &   trim(adjustl(str_extent_g))//'">'

        write(1001,'(A)')'<PCoordinates>'
        write(1001,'(A)')
     &   '<DataArray Name="x" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="y" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="z" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PCoordinates>'

        write(1001,'(A)')'<PPointData Scalars="'//
     &   trim(fnamebase)//'">'
        write(1001,'(A)')'<DataArray type="Float32" Name="'//
     &   trim(fnamebase)//'" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PPointData>'

        do nn=0,nprocs_o-1
          write(pid2s,"(i4.4)") int(nn)
          fnamel=trim(fnamebase)//'_'//trim(ts2s)//'_'//
     &     trim(pid2s)//'.vtr'
          call get_ijkrange(il, iu, jl, ju, kl, ku, nn)
          call get_pidxyz(pidx, pidy, pidz, nn)

          itmp=il-1
          jtmp=jl-1
          ktmp=kl-1
          if(pidx.ge.1)itmp=il-2
          if(pidy.ge.1)jtmp=jl-2
          if(pidz.ge.1)ktmp=kl-2
          write(str_extent_l,"(6(i7,1x))")
     &      itmp,iu-1,jtmp,ju-1,ktmp,ku-1

          write(1001,'(A)')
     &     '<Piece Extent="'//trim(adjustl(str_extent_l))//
     &                '" Source="'//trim(fnamel)//'">'
          write(1001,'(A)')'</Piece>'
        enddo !do nn=0,nprocs_o-1

        write(1001,'(A)')'</PRectilinearGrid>'
        write(1001,'(A)')'</VTKFile>'
        close(1001)      
      endif !if(gpid.eq.0)then


      ! Output of vtr file (actual local data)
      write(pid2s,"(i4.4)") int(gpid)
      call get_ijkrange(il, iu, jl, ju, kl, ku, gpid)
      call get_pidxyz(pidx, pidy, pidz, gpid)

      itmp=il-1
      jtmp=jl-1
      ktmp=kl-1
      if(pidx.ge.1)itmp=il-2
      if(pidy.ge.1)jtmp=jl-2
      if(pidz.ge.1)ktmp=kl-2
      write(str_extent_l,"(6(i7,1x))") itmp,iu-1,jtmp,ju-1,ktmp,ku-1
      fnamel=trim(PVOPATH)//'checkpoint/'//trim(subdir)//
     & trim(fnamebase)//'_'//trim(ts2s)//'_'//trim(pid2s)//'.vtr'
      open(1002,FILE=trim(fnamel),status='unknown',action='write')
      write(1002,'(A)')'<?xml version="1.0"?>'
      write(1002,'(A)')'<VTKFile type="RectilinearGrid" '//
     & 'version="0.1" byte_order="LittleEndian">'
      write(1002,'(A)')'<RectilinearGrid WholeExtent="'//
     & trim(adjustl(str_extent_g))//'">'
      write(1002,'(A)')'<Piece Extent="'//
     & trim(adjustl(str_extent_l))//'">'
      write(1002,'(A)')'<Coordinates>'
      write(1002,'(A)')
     & '<DataArray Name="x" type="Float32" format="ascii">'
      do i=itmp,iu-1
            write(1002,fmt_r) x0_o + i*dx_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="y" type="Float32" format="ascii">'
      do j=jtmp,ju-1
            write(1002,fmt_r) y0_o + j*dy_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="z" type="Float32" format="ascii">'
      do k=ktmp,ku-1
            write(1002,fmt_r) z0_o + k*dz_o
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</Coordinates>'

      write(1002,'(A)')'<PointData Scalars="'//trim(fnamebase)//'">'
      write(1002,'(A)')'<DataArray type="Float32" Name="'//
     & trim(fnamebase)//'" format="ascii">'
      itmp=1-1
      jtmp=1-1
      ktmp=1-1
      if(pidx.ge.1)itmp=1-2
      if(pidy.ge.1)jtmp=1-2
      if(pidz.ge.1)ktmp=1-2
      do k=ktmp,pv_kml-1
      do j=1,1
      do i=itmp,pv_iml-1
            write(1002,fmt_r16)a(i+1,j,k+1)
      enddo
      enddo
      enddo

      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</PointData>'
      write(1002,'(A)')'</Piece>'
      write(1002,'(A)')'</RectilinearGrid>'
      write(1002,'(A)')'</VTKFile>'
      close(1002)

      end subroutine checkpoint_r8_pvtr_2d

      ! ----------------------------------------------------------------

      subroutine checkpoint_i4_pvtr_2d_with_grid(gridobj,
     &     iml,jml,kml,ng,a,fnamebase,timestep,timeval,gpid,description)
      ! Output of pvtr requires at least one ghost layer (ng>=1)

      implicit none
      type(grid),intent(in) :: gridobj
      integer(4),intent(in) :: iml,jml,kml
      integer(4),intent(in) :: ng
      integer(4),intent(in) :: a(1-ng:iml+ng,1:1,1-ng:kml+ng)
      character(len=*),intent(in) :: fnamebase
      integer(4),intent(in) :: timestep
      real(8),intent(in) :: timeval
      integer(4),intent(in) :: gpid
      character(len=*),intent(in) :: description

      integer(4) :: i,j,k,nn
      character(len=300) :: fnameg, fnamel, subdir
      character(len=8) :: ts2s, pid2s
      integer(4) :: pidx, pidy, pidz
      integer(4) :: il, iu, jl, ju, kl, ku
      integer(4) :: itmp, jtmp, ktmp
      character(len=300) :: str_extent_g, str_extent_l
      integer(4) :: ierr
      integer(4) :: img, jmg, kmg
      real(8) :: x0, y0, z0
      real(8) :: dxl, dyl, dzl 

      if(gridobj%is_initialized .eqv. .false.)then
        write(*,*)"Error, M801 pvoutputs, output_i4_pvtr_with_grid()"
        write(*,*)"  Uninitialized grid object delivered."
        write(*,*)"  Aborting..."
        call abort()
      endif

      x0 = gridobj%x0 
      y0 = gridobj%y0 
      z0 = gridobj%z0 
      dxl = gridobj%dxl 
      dyl = gridobj%dyl 
      dzl = gridobj%dzl 

      img = gridobj%img 
      jmg = gridobj%jmg 
      kmg = gridobj%kmg 

      pidx = gridobj%pidx
      pidy = gridobj%pidy
      pidz = gridobj%pidz
      il = gridobj%il 
      iu = gridobj%iu 
      jl = gridobj%jl 
      ju = gridobj%ju 
      kl = gridobj%kl 
      ku = gridobj%ku 

      subdir = trim(fnamebase)//'/'

      ! Make a subdir
      if(gpid.eq.0)then
        call execute_command_line(
     & 'mkdir -p '//trim(CHECKPOINTPATH))
        call execute_command_line(
     & 'mkdir -p '//trim(CHECKPOINTPATH)//trim(subdir))
      endif
      call mpi_barrier(MPI_COMM_WORLD, ierr)

      write(ts2s,"(i8.8)") int(timestep)
      write(pid2s,"(i4.4)") int(gpid)
      write(str_extent_g,"(6(i7,1x))")
     &                            0,img-1,0,jmg-1,0,kmg-1
      fnameg = trim(CHECKPOINTPATH)//trim(subdir)//
     & trim(fnamebase)//'_'//trim(ts2s)//'.pvtr'

      ! Output of pvtr file (meta file to combine each local data(vtr))
      if(gpid.eq.0)then
        call output_pvd(CHECKPOINTPATH,fnamebase,subdir,timestep,
     &                                              timeval,description)
        open(1001,FILE=trim(fnameg),status='unknown',action='write')
        write(1001,'(A)')'<?xml version="1.0"?>'
        write(1001,'(A)')'<VTKFile type="PRectilinearGrid" '//
     &   'version="0.1" byte_order="LittleEndian">'
        write(1001,'(A)')'<PRectilinearGrid WholeExtent="'//
     &   trim(adjustl(str_extent_g))//'">'

        write(1001,'(A)')'<PCoordinates>'
        write(1001,'(A)')
     &   '<DataArray Name="x" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="y" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="z" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PCoordinates>'

        write(1001,'(A)')'<PPointData Scalars="'//
     &   trim(fnamebase)//'">'
        write(1001,'(A)')'<DataArray type="Int32" Name="'//
     &   trim(fnamebase)//'" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PPointData>'

        do nn=0,nprocs_o-1
          write(pid2s,"(i4.4)") int(nn)
          fnamel=trim(fnamebase)//'_'//trim(ts2s)//'_'//
     &     trim(pid2s)//'.vtr'
          call grid_get_ijkrange(il, iu, jl, ju, kl, ku, gridobj, nn)
          call grid_get_pidxyz(pidx, pidy, pidz, nn, gridobj)

          itmp=il-1
          jtmp=jl-1
          ktmp=kl-1
          if(pidx.ge.1)itmp=il-2
          if(pidy.ge.1)jtmp=jl-2
          if(pidz.ge.1)ktmp=kl-2
          write(str_extent_l,"(6(i7,1x))")
     &      itmp,iu-1,jtmp,ju-1,ktmp,ku-1

          write(1001,'(A)')
     &     '<Piece Extent="'//trim(adjustl(str_extent_l))//
     &                '" Source="'//trim(fnamel)//'">'
          write(1001,'(A)')'</Piece>'
        enddo !do nn=0,nprocs_o-1

        write(1001,'(A)')'</PRectilinearGrid>'
        write(1001,'(A)')'</VTKFile>'
        close(1001)      
      endif !if(gpid.eq.0)then


      ! Output of vtr file (actual local data)
      write(pid2s,"(i4.4)") int(gpid)
      call grid_get_ijkrange(il, iu, jl, ju, kl, ku, gridobj, gpid)
      call grid_get_pidxyz(pidx, pidy, pidz, gpid, gridobj)

      itmp=il-1
      jtmp=jl-1
      ktmp=kl-1
      if(pidx.ge.1)itmp=il-2
      if(pidy.ge.1)jtmp=jl-2
      if(pidz.ge.1)ktmp=kl-2
      write(str_extent_l,"(6(i7,1x))") itmp,iu-1,jtmp,ju-1,ktmp,ku-1
      fnamel=trim(CHECKPOINTPATH)//trim(subdir)//
     & trim(fnamebase)//'_'//trim(ts2s)//'_'//trim(pid2s)//'.vtr'
      open(1002,FILE=trim(fnamel),status='unknown',action='write')
      write(1002,'(A)')'<?xml version="1.0"?>'
      write(1002,'(A)')'<VTKFile type="RectilinearGrid" '//
     & 'version="0.1" byte_order="LittleEndian">'
      write(1002,'(A)')'<RectilinearGrid WholeExtent="'//
     & trim(adjustl(str_extent_g))//'">'
      write(1002,'(A)')'<Piece Extent="'//
     & trim(adjustl(str_extent_l))//'">'
      write(1002,'(A)')'<Coordinates>'
      write(1002,'(A)')
     & '<DataArray Name="x" type="Float32" format="ascii">'
      do i=itmp,iu-1
            write(1002,fmt_r) x0 + i*dxl
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="y" type="Float32" format="ascii">'
      do j=jtmp,ju-1
            write(1002,fmt_r) y0 + j*dyl
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="z" type="Float32" format="ascii">'
      do k=ktmp,ku-1
            write(1002,fmt_r) z0 + k*dzl
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</Coordinates>'

      write(1002,'(A)')'<PointData Scalars="'//trim(fnamebase)//'">'
      write(1002,'(A)')'<DataArray type="Int32" Name="'//
     & trim(fnamebase)//'" format="ascii">'
      itmp=1-1
      jtmp=1-1
      ktmp=1-1
      if(pidx.ge.1)itmp=1-2
      if(pidy.ge.1)jtmp=1-2
      if(pidz.ge.1)ktmp=1-2
      do k=ktmp,kml-1
      do j=1,1
      do i=itmp,iml-1
            write(1002,'(I0)')a(i+1,j,k+1)
      enddo
      enddo
      enddo

      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</PointData>'
      write(1002,'(A)')'</Piece>'
      write(1002,'(A)')'</RectilinearGrid>'
      write(1002,'(A)')'</VTKFile>'
      close(1002)

      end subroutine checkpoint_i4_pvtr_2d_with_grid

      ! ----------------------------------------------------------------

      subroutine checkpoint_r8_pvtr_2d_with_grid(gridobj,
     &     iml,jml,kml,ng,a,fnamebase,timestep,timeval,gpid,description)
      ! Output of pvtr requires at least one ghost layer (ng>=1)

      implicit none
      type(grid),intent(in) :: gridobj
      integer(4),intent(in) :: iml,jml,kml
      integer(4),intent(in) :: ng
      real(8),intent(in) :: a(1-ng:iml+ng,1:1,1-ng:kml+ng)
      character(len=*),intent(in) :: fnamebase
      integer(4),intent(in) :: timestep
      real(8),intent(in) :: timeval
      integer(4),intent(in) :: gpid
      character(len=*),intent(in) :: description

      integer(4) :: i,j,k,nn
      character(len=300) :: fnameg, fnamel, subdir
      character(len=8) :: ts2s, pid2s
      integer(4) :: pidx, pidy, pidz
      integer(4) :: il, iu, jl, ju, kl, ku
      integer(4) :: itmp, jtmp, ktmp
      character(len=300) :: str_extent_g, str_extent_l
      integer(4) :: ierr
      integer(4) :: img, jmg, kmg
      real(8) :: x0, y0, z0
      real(8) :: dxl, dyl, dzl 

      if(gridobj%is_initialized .eqv. .false.)then
        write(*,*)"Error, M801 pvoutputs, output_i4_pvtr_with_grid()"
        write(*,*)"  Uninitialized grid object delivered."
        write(*,*)"  Aborting..."
        call abort()
      endif

      x0 = gridobj%x0 
      y0 = gridobj%y0 
      z0 = gridobj%z0 
      dxl = gridobj%dxl 
      dyl = gridobj%dyl 
      dzl = gridobj%dzl 

      img = gridobj%img 
      jmg = gridobj%jmg 
      kmg = gridobj%kmg 

      pidx = gridobj%pidx
      pidy = gridobj%pidy
      pidz = gridobj%pidz
      il = gridobj%il 
      iu = gridobj%iu 
      jl = gridobj%jl 
      ju = gridobj%ju 
      kl = gridobj%kl 
      ku = gridobj%ku 

      subdir = trim(fnamebase)//'/'

      ! Make a subdir
      if(gpid.eq.0)then
        call execute_command_line(
     & 'mkdir -p '//trim(CHECKPOINTPATH))
        call execute_command_line(
     & 'mkdir -p '//trim(CHECKPOINTPATH)//trim(subdir))
      endif
      call mpi_barrier(MPI_COMM_WORLD, ierr)

      write(ts2s,"(i8.8)") int(timestep)
      write(pid2s,"(i4.4)") int(gpid)
      write(str_extent_g,"(6(i7,1x))")
     &                            0,img-1,0,jmg-1,0,kmg-1
      fnameg = trim(CHECKPOINTPATH)//trim(subdir)//
     & trim(fnamebase)//'_'//trim(ts2s)//'.pvtr'


      ! Output of pvtr file (meta file to combine each local data(vtr))
      if(gpid.eq.0)then
        call output_pvd(CHECKPOINTPATH,fnamebase,subdir,timestep,
     &                                              timeval,description)
        open(1001,FILE=trim(fnameg),status='unknown',action='write')
        write(1001,'(A)')'<?xml version="1.0"?>'
        write(1001,'(A)')'<VTKFile type="PRectilinearGrid" '//
     &   'version="0.1" byte_order="LittleEndian">'
        write(1001,'(A)')'<PRectilinearGrid WholeExtent="'//
     &   trim(adjustl(str_extent_g))//'">'

        write(1001,'(A)')'<PCoordinates>'
        write(1001,'(A)')
     &   '<DataArray Name="x" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="y" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')
     &   '<DataArray Name="z" type="Float32" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PCoordinates>'

        write(1001,'(A)')'<PPointData Scalars="'//
     &   trim(fnamebase)//'">'
        write(1001,'(A)')'<DataArray type="Float32" Name="'//
     &   trim(fnamebase)//'" format="ascii">'
        write(1001,'(A)')'</DataArray>'
        write(1001,'(A)')'</PPointData>'

        do nn=0,nprocs_o-1
          write(pid2s,"(i4.4)") int(nn)
          fnamel=trim(fnamebase)//'_'//trim(ts2s)//'_'//
     &     trim(pid2s)//'.vtr'
          call grid_get_ijkrange(il, iu, jl, ju, kl, ku, gridobj, nn)
          call grid_get_pidxyz(pidx, pidy, pidz, nn, gridobj)

          itmp=il-1
          jtmp=jl-1
          ktmp=kl-1
          if(pidx.ge.1)itmp=il-2
          if(pidy.ge.1)jtmp=jl-2
          if(pidz.ge.1)ktmp=kl-2
          write(str_extent_l,"(6(i7,1x))")
     &      itmp,iu-1,jtmp,ju-1,ktmp,ku-1

          write(1001,'(A)')
     &     '<Piece Extent="'//trim(adjustl(str_extent_l))//
     &                '" Source="'//trim(fnamel)//'">'
          write(1001,'(A)')'</Piece>'
        enddo !do nn=0,nprocs_o-1

        write(1001,'(A)')'</PRectilinearGrid>'
        write(1001,'(A)')'</VTKFile>'
        close(1001)      
      endif !if(gpid.eq.0)then


      ! Output of vtr file (actual local data)
      write(pid2s,"(i4.4)") int(gpid)
      call grid_get_ijkrange(il, iu, jl, ju, kl, ku, gridobj, gpid)
      call grid_get_pidxyz(pidx, pidy, pidz, gpid, gridobj)

      itmp=il-1
      jtmp=jl-1
      ktmp=kl-1
      if(pidx.ge.1)itmp=il-2
      if(pidy.ge.1)jtmp=jl-2
      if(pidz.ge.1)ktmp=kl-2
      write(str_extent_l,"(6(i7,1x))") itmp,iu-1,jtmp,ju-1,ktmp,ku-1
      fnamel=trim(CHECKPOINTPATH)//trim(subdir)//
     & trim(fnamebase)//'_'//trim(ts2s)//'_'//trim(pid2s)//'.vtr'
      open(1002,FILE=trim(fnamel),status='unknown',action='write')
      write(1002,'(A)')'<?xml version="1.0"?>'
      write(1002,'(A)')'<VTKFile type="RectilinearGrid" '//
     & 'version="0.1" byte_order="LittleEndian">'
      write(1002,'(A)')'<RectilinearGrid WholeExtent="'//
     & trim(adjustl(str_extent_g))//'">'
      write(1002,'(A)')'<Piece Extent="'//
     & trim(adjustl(str_extent_l))//'">'
      write(1002,'(A)')'<Coordinates>'
      write(1002,'(A)')
     & '<DataArray Name="x" type="Float32" format="ascii">'
      do i=itmp,iu-1
            write(1002,fmt_r) x0 + i*dxl
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="y" type="Float32" format="ascii">'
      do j=jtmp,ju-1
            write(1002,fmt_r) y0 + j*dyl
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')
     & '<DataArray Name="z" type="Float32" format="ascii">'
      do k=ktmp,ku-1
            write(1002,fmt_r) z0 + k*dzl
      enddo
      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</Coordinates>'

      write(1002,'(A)')'<PointData Scalars="'//trim(fnamebase)//'">'
      write(1002,'(A)')'<DataArray type="Float32" Name="'//
     & trim(fnamebase)//'" format="ascii">'
      itmp=1-1
      jtmp=1-1
      ktmp=1-1
      if(pidx.ge.1)itmp=1-2
      if(pidy.ge.1)jtmp=1-2
      if(pidz.ge.1)ktmp=1-2
      do k=ktmp,gridobj%kml-1
      do j=1,1
      do i=itmp,gridobj%iml-1
            write(1002,fmt_r)a(i+1,j,k+1)
      enddo
      enddo
      enddo

      write(1002,'(A)')'</DataArray>'
      write(1002,'(A)')'</PointData>'
      write(1002,'(A)')'</Piece>'
      write(1002,'(A)')'</RectilinearGrid>'
      write(1002,'(A)')'</VTKFile>'
      close(1002)

      end subroutine checkpoint_r8_pvtr_2d_with_grid


      ! ----------------------------------------------------------------

      subroutine read_iternum_from_pvd(timestep,path,fnamebase,time)
      implicit none
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

      subroutine pvd_data_parser(pvdobj, inputstr)
      implicit none 
      type(pvddata),intent(inout) :: pvdobj 
      character(len=*),intent(in) :: inputstr 

      character(len=100) :: timestepstr
      character(len=100) :: groupstr
      character(len=100) :: partstr
      character(len=100) :: filestr
      character(len=100) :: namestr 
      logical :: found 

      ! type :: pvddata
      !   real(8) :: timestep
      !   integer(4) :: group 
      !   integer(4) :: part 
      !   character(len=300) :: file 
      !   character(len=30) :: name 
      ! end type

      ! Parsing input string and get the pvd data 
      call get_option_value(inputstr,"timestep",timestepstr,found)
      call get_option_value(inputstr,"group",groupstr,found)
      call get_option_value(inputstr,"part",partstr,found)
      call get_option_value(inputstr,"file",filestr,found)
      call get_option_value(inputstr,"name",namestr,found)

      read(timestepstr,*) pvdobj%timestep
      read(groupstr,*) pvdobj%group
      read(partstr,*) pvdobj%part 
      pvdobj%file=trim(filestr)
      pvdobj%name=trim(namestr)

      return 
      end subroutine pvd_data_parser

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
        print *, "Error reading integers from the string:", str_to_read
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
      integer(4) :: io

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
        open(newunit=io,FILE=trim(fnamevtr),status='old',action='read')
        read(io,*) ! Skipping <?xml version="1.0"?>
        read(io,*) ! Skipping <VTKFile type="RectilinearGrid" version="0.1" byte_order="LittleEndian">
        read(io,*) ! Skipping Whole Extent line which we already know
        read(io,*) ! Skipping Piece Extent line which we already know 
        read(io,*) ! Skipping <Coordinate> line
        read(io,*) ! Skipping <DataArray Name="x" type="Float32" format="ascii">
        do l=1,numdat_x 
          read(io,*) ! Now x coordinates comes (all skip)
        enddo
        read(io,*) ! Skipping </DataArray>
        read(io,*) ! Skipping <DataArray Name="y" type="Float32" format="ascii">
        do l=1,numdat_y 
          read(io,*) ! Now y coordinates comes (all skip)
        enddo 
        read(io,*) ! Skipping </DataArray>
        read(io,*) ! Skipping <DataArray Name="z" type="Float32" format="ascii">
        do l=1,numdat_z 
          read(io,*) ! Now z coordinates comes (all skip)
        enddo 
        read(io,*) ! Skipping </DataArray>
        read(io,*) ! Skipping </Coordinates>
        read(io,*) ! Skipping <PointData Scalars="SI">
        read(io,*) ! Skipping <DataArray type="Float32" Name="SI" format="ascii">
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
          read(io,*) a(i,j,k) 
        enddo
        enddo
        enddo

        read(io,*) ! Skipping </DataArray>
        read(io,*) ! Skipping </PointData>
        read(io,*) ! Skipping </Piece>
        read(io,*) ! Skipping </RectilinearGrid>
        read(io,*) ! Skipping </VTKFile>

        close(io)
      endif !if (fexist .eqv. .false.) then ... else
        
      return 
      end subroutine read_vtr_r8_2d

      ! ----------------------------------------------------------------


      subroutine read_vtr_i4_2d(gridobj, ng, a, path, fnamebase, iter)
      type(grid),intent(in) :: gridobj 
      integer(4),intent(in) :: ng
      integer(4),intent(inout) :: 
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
      integer(4) :: io 

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
        open(newunit=io,FILE=trim(fnamevtr),status='old',action='read')
        read(io,*) ! Skipping <?xml version="1.0"?>
        read(io,*) ! Skipping <VTKFile type="RectilinearGrid" version="0.1" byte_order="LittleEndian">
        read(io,*) ! Skipping Whole Extent line which we already know
        read(io,*) ! Skipping Piece Extent line which we already know 
        read(io,*) ! Skipping <Coordinate> line
        read(io,*) ! Skipping <DataArray Name="x" type="Float32" format="ascii">
        do l=1,numdat_x 
          read(io,*) ! Now x coordinates comes (all skip)
        enddo
        read(io,*) ! Skipping </DataArray>
        read(io,*) ! Skipping <DataArray Name="y" type="Float32" format="ascii">
        do l=1,numdat_y 
          read(io,*) ! Now y coordinates comes (all skip)
        enddo 
        read(io,*) ! Skipping </DataArray>
        read(io,*) ! Skipping <DataArray Name="z" type="Float32" format="ascii">
        do l=1,numdat_z 
          read(io,*) ! Now z coordinates comes (all skip)
        enddo 
        read(io,*) ! Skipping </DataArray>
        read(io,*) ! Skipping </Coordinates>
        read(io,*) ! Skipping <PointData Scalars="SI">
        read(io,*) ! Skipping <DataArray type="Float32" Name="SI" format="ascii">
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
          read(io,*) a(i,j,k) 
        enddo
        enddo
        enddo

        read(io,*) ! Skipping </DataArray>
        read(io,*) ! Skipping </PointData>
        read(io,*) ! Skipping </Piece>
        read(io,*) ! Skipping </RectilinearGrid>
        read(io,*) ! Skipping </VTKFile>

        close(io)
      endif !if (fexist .eqv. .false.) then ... else
        
      return 
      end subroutine read_vtr_i4_2d

      ! ----------------------------------------------------------------

      end module pvoutputs  