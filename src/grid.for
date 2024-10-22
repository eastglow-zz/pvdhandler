      module grid_module
      use mpi 

      implicit none 

      type :: grid
        character(len=30) :: name
        real(8) :: x0  ! in meter unit
        real(8) :: y0  ! in meter unit
        real(8) :: z0  ! in meter unit
        real(8) :: dxl ! in meter unit
        real(8) :: dyl ! in meter unit
        real(8) :: dzl ! in meter unit
        integer(4) :: img
        integer(4) :: jmg
        integer(4) :: kmg
        integer(4) :: nprocs
        integer(4) :: npx
        integer(4) :: npy
        integer(4) :: npz
        integer(4) :: pid
        integer(4) :: pidx
        integer(4) :: pidy 
        integer(4) :: pidz
        integer(4) :: iml
        integer(4) :: jml
        integer(4) :: kml
        integer(4) :: il
        integer(4) :: iu
        integer(4) :: jl
        integer(4) :: ju
        integer(4) :: kl
        integer(4) :: ku
        logical :: is_initialized = .false.
      end type  !type :: grid

      contains

      ! ----------------------------------------------------------------

      subroutine grid_init(gridobj,gridname,x0,y0,z0,dxl,dyl,dzl,
     &                               img,jmg,kmg,nprocs,npx,npy,npz,pid)
      implicit none

      type(grid),intent(inout) :: gridobj
      character(len=*),intent(in) :: gridname
      real(8),intent(in) :: x0
      real(8),intent(in) :: y0
      real(8),intent(in) :: z0
      real(8),intent(in) :: dxl
      real(8),intent(in) :: dyl
      real(8),intent(in) :: dzl
      integer(4),intent(in) :: img
      integer(4),intent(in) :: jmg
      integer(4),intent(in) :: kmg
      integer(4),intent(in) :: nprocs
      integer(4),intent(in) :: npx
      integer(4),intent(in) :: npy
      integer(4),intent(in) :: npz
      integer(4),intent(in) :: pid

      integer(4) :: pidx
      integer(4) :: pidy
      integer(4) :: pidz
      integer(4) :: il,iu,jl,ju,kl,ku
      
      ! Items in gridobj, gridobj%*
      !character(len=30) :: name
      !real(8) :: x0
      !real(8) :: y0
      !real(8) :: z0
      !real(8) :: dxl
      !real(8) :: dyl
      !real(8) :: dzl
      !integer(4) :: img  ! Directly given
      !integer(4) :: jmg  ! Directly given
      !integer(4) :: kmg  ! Directly given
      !integer(4) :: nprocs  ! Directly given
      !integer(4) :: npx  ! Directly given
      !integer(4) :: npy  ! Directly given
      !integer(4) :: npz  ! Directly given
      !integer(4) :: pid  ! Directly given
      !integer(4) :: pidx ! Computed here
      !integer(4) :: pidy ! Computed here
      !integer(4) :: pidz ! Computed here
      !integer(4) :: iml  ! Computed here
      !integer(4) :: jml  ! Computed here
      !integer(4) :: kml  ! Computed here
      !integer(4) :: il  ! Computed here
      !integer(4) :: iu  ! Computed here
      !integer(4) :: jl  ! Computed here
      !integer(4) :: ju  ! Computed here
      !integer(4) :: kl  ! Computed here
      !integer(4) :: ku  ! Computed here
      !logical :: is_initialized = .false. ! Computed here

      if(gridobj%is_initialized .eqv. .true.) then
        if(pid .eq. 0) then
          write(*,*)"Error, grid module, grid_init()"
          write(*,*)"  Cannot init. the grid already initialized"
          write(*,*)"  grid obj. name: ",trim(gridobj%name)
          write(*,*)"  grid name tried to be overwriteen: ", gridname 
          write(*,*)"  Please use grid_reinit() to overwrite the grid."
          write(*,*)"  Aborting..."
          call abort()
        endif
      endif

      gridobj%name = trim(gridname)

      gridobj%x0 = x0 
      gridobj%y0 = y0 
      gridobj%z0 = z0
      gridobj%dxl = dxl 
      gridobj%dyl = dyl 
      gridobj%dzl = dzl 

      gridobj%img = img
      gridobj%jmg = jmg
      gridobj%kmg = kmg 
      gridobj%nprocs = nprocs 
      gridobj%npx = npx 
      gridobj%npy = npy 
      gridobj%npz = npz 
      gridobj%pid = pid
      
      ! Calc. of pidx,pidy,pidz,il,iu,jl,ju,kl,ku,iml,jml,kml
      ! calc. gridobj%pidx, ...pidy, ...pidz
      call grid_get_pidxyz(pidx,pidy,pidz, pid,gridobj)
      gridobj%pidx = pidx
      gridobj%pidy = pidy
      gridobj%pidz = pidz

      call grid_get_ijkrange(il,iu,jl,ju,kl,ku,gridobj,pid)
      gridobj%il = il 
      gridobj%iu = iu 
      gridobj%jl = jl 
      gridobj%ju = ju 
      gridobj%kl = kl 
      gridobj%ku = ku 

      gridobj%iml=iu-il+1
      gridobj%jml=ju-jl+1
      gridobj%kml=ku-kl+1

      gridobj%is_initialized = .true.
      
      return 
      end subroutine grid_init

      ! ----------------------------------------------------------------

      subroutine grid_reinit(gridobj,gridname,img,jmg,kmg,nprocs,
     &                                                  npx,npy,npz,pid)
      implicit none

      type(grid),intent(inout) :: gridobj
      character(len=*),intent(in) :: gridname
      integer(4),intent(in) :: img
      integer(4),intent(in) :: jmg
      integer(4),intent(in) :: kmg
      integer(4),intent(in) :: nprocs
      integer(4),intent(in) :: npx
      integer(4),intent(in) :: npy
      integer(4),intent(in) :: npz
      integer(4),intent(in) :: pid

      integer(4) :: pidx
      integer(4) :: pidy
      integer(4) :: pidz
      integer(4) :: il,iu,jl,ju,kl,ku
      
      ! Items in gridobj, gridobj%*
      !character(len=30) :: name
      !integer(4) :: img  ! Directly given
      !integer(4) :: jmg  ! Directly given
      !integer(4) :: kmg  ! Directly given
      !integer(4) :: nprocs  ! Directly given
      !integer(4) :: npx  ! Directly given
      !integer(4) :: npy  ! Directly given
      !integer(4) :: npz  ! Directly given
      !integer(4) :: pid  ! Directly given
      !integer(4) :: pidx ! Computed here
      !integer(4) :: pidy ! Computed here
      !integer(4) :: pidz ! Computed here
      !integer(4) :: iml  ! Computed here
      !integer(4) :: jml  ! Computed here
      !integer(4) :: kml  ! Computed here
      !integer(4) :: il  ! Computed here
      !integer(4) :: iu  ! Computed here
      !integer(4) :: jl  ! Computed here
      !integer(4) :: ju  ! Computed here
      !integer(4) :: kl  ! Computed here
      !integer(4) :: ku  ! Computed here
      !logical :: is_initialized = .false. ! Computed here

      if(gridobj%is_initialized .eqv. .false.) then
        if(pid .eq. 0) then
          write(*,*)"Error, grid module, grid_reinit()"
          write(*,*)"  Cannot reinit. the grid uninitialized"
          write(*,*)"  grid name tried to be overwriteen: ", gridname 
          write(*,*)"  Please use grid_init() to init. the grid."
          write(*,*)"  Aborting..."
          call abort()
        endif
      endif

      gridobj%name = trim(gridname)

      gridobj%img = img
      gridobj%jmg = jmg
      gridobj%kmg = kmg 
      gridobj%nprocs = nprocs 
      gridobj%npx = npx 
      gridobj%npy = npy 
      gridobj%npz = npz 
      gridobj%pid = pid
      
      call grid_get_pidxyz(pidx,pidy,pidz, pid,gridobj)
      gridobj%pidx = pidx
      gridobj%pidy = pidy
      gridobj%pidz = pidz

      call grid_get_ijkrange(il,iu,jl,ju,kl,ku,gridobj,pid)
      gridobj%il = il 
      gridobj%iu = iu 
      gridobj%jl = jl 
      gridobj%ju = ju 
      gridobj%kl = kl 
      gridobj%ku = ku 

      gridobj%iml=iu-il+1
      gridobj%jml=ju-jl+1
      gridobj%kml=ku-kl+1

      gridobj%is_initialized = .true.
      
      return 
      end subroutine grid_reinit

      ! ----------------------------------------------------------------

      subroutine grid_get_ijkrange(il,iu,jl,ju,kl,ku,gridobj,pid)
      implicit none
      integer(4),intent(inout) :: il,iu,jl,ju,kl,ku
      type(grid),intent(in) :: gridobj
      integer(4),intent(in) :: pid

      integer(4) :: pidx, pidy, pidz

      call grid_get_pidxyz(pidx,pidy,pidz, pid,gridobj)
      
      ! x bounds
      call grid_para_range(1,gridobj%img,gridobj%npx,pidx,il,iu)
      ! y bounds
      call grid_para_range(1,gridobj%jmg,gridobj%npy,pidy,jl,ju)
      ! z bounds
      call grid_para_range(1,gridobj%kmg,gridobj%npz,pidz,kl,ku)


      return
      end subroutine grid_get_ijkrange

      ! ----------------------------------------------------------------

      subroutine grid_get_pidxyz(pidx,pidy,pidz,pid,gridobj)
      implicit none
      integer(4),intent(inout):: pidx, pidy, pidz
      integer(4),intent(in) :: pid
      type(grid),intent(in) :: gridobj

      integer(4) :: gnpx, gnpy, gnpz
      

      gnpx = gridobj%npx 
      gnpy = gridobj%npy 
      gnpz = gridobj%npz 

      pidz = int(pid/gnpx/gnpy)
      pidy = int((pid - pidz*gnpx*gnpy)/gnpx)
      pidx = pid - pidz*gnpx*gnpy - pidy*gnpx

      return
      end subroutine grid_get_pidxyz

      ! ----------------------------------------------------------------

      subroutine grid_get_global_ijk(ig,jg,kg, 
     &                                ilocal,jlocal,klocal, gridobj,pid)
      implicit none
      integer(4),intent(inout) :: ig, jg, kg !global ijk index
      integer(4),intent(in) :: ilocal, jlocal, klocal !local ijk index
      type(grid),intent(in) :: gridobj
      integer(4),intent(in) :: pid

      integer(4) :: il,iu, jl,ju, kl,ku

      call grid_get_ijkrange(il,iu,jl,ju,kl,ku,gridobj,pid)
      ig = il + (ilocal-1)  !for spatial index starts from 1, not 0
      jg = jl + (jlocal-1)  !for spatial index starts from 1, not 0
      kg = kl + (klocal-1)  !for spatial index starts from 1, not 0

      return
      end subroutine grid_get_global_ijk

      ! ----------------------------------------------------------------

      subroutine grid_para_range(n1,n2,npart,irank,ista,iend)
      ! Deviding system size bc_kmori into KM = bc_kmori/(porcess number)
      ! npart: total number of porcesses (npart=48 in my computer)
      ! irank: current process name (irank=1~47)
      ! ista and iend: first and last grids (J-direction) on each process in the total grid system
      ! use: CALL para_range(1,bc_kmori,nprocs,myrank,ista,iend)
      implicit none
      integer(4),intent(in) :: n1,n2,npart
      integer(4),intent(inout) :: irank,ista,iend

      integer(4) :: iwork1, iwork2

      iwork1=(n2-n1+1)/npart
      iwork2=MOD(n2-n1+1,npart)
      ista=irank*iwork1+n1+MIN(irank,iwork2)
      iend=ista+iwork1-1
      if(iwork2.gt.irank) iend=iend+1

      return
      end subroutine grid_para_range

      ! ----------------------------------------------------------------

      integer(4) function grid_get_pid_from_globalijk(gridobj,ig,jg,kg)
      implicit none
      type(grid),intent(in) :: gridobj
      integer(4),intent(in) :: ig,jg,kg

      integer(4) :: igf, jgf, kgf
      integer(4) :: pidx, pidy, pidz

      if(gridobj%is_initialized .eqv. .false.)then
        write(*,*)
     &    "Error, grid module, integer(4) grid_get_pid_from_globalijk()"
        write(*,*)"  Uninitialized grid obj. given."
        write(*,*)"  Aborting..."
        call abort()
      endif

      if(ig.gt.gridobj%img.or.ig.lt.1)then
        write(*,*)
     &    "Error, grid module, integer(4) grid_get_pid_from_globalijk()"
        write(*,*)"  Out of range in x grid index. ig: ",ig
        write(*,*)"  grid name: ",trim(gridobj%name)
        write(*,*)"  y grid range limit: ",gridobj%img
        write(*,*)"  Aborting..."
        call abort()
      endif

      if(jg.gt.gridobj%jmg.or.jg.lt.1)then
        write(*,*)
     &    "Error, grid module, integer(4) grid_get_pid_from_globalijk()"
        write(*,*)"  Out of range in y grid index. jg: ",jg
        write(*,*)"  grid name: ",trim(gridobj%name)
        write(*,*)"  z grid range limit: ",gridobj%jmg
        write(*,*)"  Aborting..."
        call abort()
      endif

      if(kg.gt.gridobj%kmg.or.kg.lt.1)then
        write(*,*)
     &    "Error, grid module, integer(4) grid_get_pid_from_globalijk()"
        write(*,*)"  Out of range in x grid index. ig: ",kg
        write(*,*)"  grid name: ",trim(gridobj%name)
        write(*,*)"  x grid range limit: ",gridobj%kmg
        write(*,*)"  Aborting..."
        call abort()
      endif

      pidx = int((ig-1)/gridobj%iml)
      pidy = int((jg-1)/gridobj%jml)
      pidz = int((kg-1)/gridobj%kml)

      grid_get_pid_from_globalijk = 
     &                            grid_get_pid(gridobj,pidx, pidy, pidz)
      return
      end function grid_get_pid_from_globalijk


      ! ----------------------------------------------------------------

      integer(4) function grid_get_pid(gridobj,pidx, pidy, pidz)
      implicit none
      type(grid),intent(in) :: gridobj
      integer(4),intent(in) :: pidx, pidy, pidz

      integer(4) :: pidxw, pidyw, pidzw

      if(gridobj%is_initialized .eqv. .false.)then
        write(*,*)
     &    "Error, grid module, integer(4) grid_get_pid()"
        write(*,*)"  Uninitialized grid obj. given."
        write(*,*)"  Aborting..."
        call abort()
      endif

      pidxw = pidx
      if(pidxw.gt.gridobj%npx-1)pidxw=0
      if(pidxw.lt.0)pidxw=gridobj%npx-1

      pidyw = pidy
      if(pidyw.gt. gridobj%npy-1)pidyw=0
      if(pidyw.lt.0)pidyw= gridobj%npy-1

      pidzw = pidz
      if(pidzw.gt.gridobj%npz-1)pidzw=0
      if(pidzw.lt.0)pidzw=gridobj%npz-1

      grid_get_pid = pidxw + pidyw*gridobj%npx + 
     &                                     pidzw*gridobj%npx*gridobj%npy
      return
      end function grid_get_pid

      ! ----------------------------------------------------------------

      subroutine grid_get_localijk_from_globalijk(iloc,jloc,kloc,
     &                                                 ig,jg,kg,gridobj)
      implicit none
      integer(4),intent(inout) :: iloc,jloc,kloc
      integer(4),intent(in) :: ig,jg,kg
      type(grid),intent(in) :: gridobj

      integer(4) :: pid, pidx, pidy, pidz

      if(gridobj%is_initialized .eqv. .false.)then
        write(*,*)
     &    "Error, grid module, grid_get_localijk_from_globalijk()"
        write(*,*)"  Uninitialized grid obj. given."
        write(*,*)"  Aborting..."
        call abort()
      endif

      if(ig.gt.gridobj%img.or.ig.lt.1)then
        write(*,*)
     &    "Error, grid module, grid_get_localijk_from_globalijk()"
        write(*,*)"  Out of range in x grid index. ig: ",ig
        write(*,*)"  grid name: ",trim(gridobj%name)
        write(*,*)"  x grid range limit: ",gridobj%img
        write(*,*)"  Aborting..."
        call abort()
      endif

      if(jg.gt.gridobj%jmg.or.jg.lt.1)then
        write(*,*)
     &    "Error, grid module, grid_get_localijk_from_globalijk()"
        write(*,*)"  Out of range in y grid index. jg: ",jg
        write(*,*)"  grid name: ",trim(gridobj%name)
        write(*,*)"  y grid range limit: ",gridobj%jmg
        write(*,*)"  Aborting..."
        call abort()
      endif

      if(kg.gt.gridobj%kmg.or.kg.lt.1)then
        write(*,*)
     &    "Error, grid module, grid_get_localijk_from_globalijk()"
        write(*,*)"  Out of range in y grid index. jg: ",kg
        write(*,*)"  grid name: ",trim(gridobj%name)
        write(*,*)"  z grid range limit: ",gridobj%kmg
        write(*,*)"  Aborting..."
        call abort()
      endif

      pid = grid_get_pid_from_globalijk(gridobj,ig,jg,kg)
      call grid_get_pidxyz(pidx, pidy, pidz, pid, gridobj)

      iloc=ig-pidx*gridobj%iml
      jloc=jg-pidy*gridobj%jml
      kloc=kg-pidz*gridobj%kml

      return
      end subroutine grid_get_localijk_from_globalijk

      ! ----------------------------------------------------------------

      subroutine grid_write_in_file(gridobj, path, fname, pid)
      implicit none

      type(grid),intent(in) :: gridobj
      character(len=*),intent(in) :: path 
      character(len=*),intent(in) :: fname
      integer(4),intent(in) :: pid

      character(len=4) ::tpid
      integer(4) :: ierr 

      ! Items in gridobj, gridobj%*
      !character(len=30) :: name
      !real(8) :: x0
      !real(8) :: y0
      !real(8) :: z0
      !real(8) :: dxl
      !real(8) :: dyl
      !real(8) :: dzl
      !integer(4) :: img  ! Directly given
      !integer(4) :: jmg  ! Directly given
      !integer(4) :: kmg  ! Directly given
      !integer(4) :: nprocs  ! Directly given
      !integer(4) :: npx  ! Directly given
      !integer(4) :: npy  ! Directly given
      !integer(4) :: npz  ! Directly given
      !integer(4) :: pid  ! Directly given
      !integer(4) :: pidx ! Computed here
      !integer(4) :: pidy ! Computed here
      !integer(4) :: pidz ! Computed here
      !integer(4) :: iml  ! Computed here
      !integer(4) :: jml  ! Computed here
      !integer(4) :: kml  ! Computed here
      !integer(4) :: il  ! Computed here
      !integer(4) :: iu  ! Computed here
      !integer(4) :: jl  ! Computed here
      !integer(4) :: ju  ! Computed here
      !integer(4) :: kl  ! Computed here
      !integer(4) :: ku  ! Computed here
      !logical :: is_initialized = .false. ! Computed here

      if(gridobj%is_initialized .eqv. .false.)then
        write(*,*)"Error, grid module, grid_write_in_file()"
        write(*,*)"  Uninitialized grid object delivered."
        write(*,*)"  Aborting..."
        call abort()
      endif

      ! Make a subdir
      if(pid.eq.0)then
        call execute_command_line(
     &                 'mkdir -p '//trim(path)//'/'//trim(gridobj%name))
      endif
      call mpi_barrier(MPI_COMM_WORLD, ierr)

      write(tpid,'(i4.4)')int(pid)

      ! JSON compatible format
      open(100,file=trim(path)//'/'//trim(gridobj%name)//'/'
     &               //fname//'_'//tpid,status='unknown',action='write')
      write(100,*)'{'
      write(100,*)'"name":',trim(gridobj%name),','
      write(100,*)'"x0":',gridobj%x0,','
      write(100,*)'"y0":',gridobj%y0,','
      write(100,*)'"z0":',gridobj%z0,','
      write(100,*)'"dxl":',gridobj%dxl,','
      write(100,*)'"dyl":',gridobj%dyl,','
      write(100,*)'"dzl":',gridobj%dzl,','
      write(100,*)'"img":',gridobj%img,','
      write(100,*)'"jmg":',gridobj%jmg,','
      write(100,*)'"kmg":',gridobj%kmg,','
      write(100,*)'"nprocs":',gridobj%nprocs,','
      write(100,*)'"npx":',gridobj%npx,','
      write(100,*)'"npy":',gridobj%npy,','
      write(100,*)'"npz":',gridobj%npz,','
      write(100,*)'"pid":',gridobj%pid,','
      write(100,*)'"pidx":',gridobj%pidx,','
      write(100,*)'"pidy":',gridobj%pidy,','
      write(100,*)'"pidz":',gridobj%pidz,','
      write(100,*)'"iml":',gridobj%iml,','
      write(100,*)'"jml":',gridobj%jml,','
      write(100,*)'"kml":',gridobj%kml,','
      write(100,*)'"il":',gridobj%il,','
      write(100,*)'"iu":',gridobj%iu,','
      write(100,*)'"jl":',gridobj%jl,','
      write(100,*)'"ju":',gridobj%ju,','
      write(100,*)'"kl":',gridobj%kl,','
      write(100,*)'"ku":',gridobj%ku,','
      write(100,*)'"vectorization_ordering":',"fortran"
      write(100,*)'}'
      close(100)

      return
      end subroutine grid_write_in_file

      ! ----------------------------------------------------------------

      subroutine grid_verbose(gridobj)
      implicit none

      type(grid),intent(in) :: gridobj

      ! Items in gridobj, gridobj%*
      !character(len=30) :: name
      !real(8) :: x0
      !real(8) :: y0
      !real(8) :: z0
      !real(8) :: dxl
      !real(8) :: dyl
      !real(8) :: dzl
      !integer(4) :: img  ! Directly given
      !integer(4) :: jmg  ! Directly given
      !integer(4) :: kmg  ! Directly given
      !integer(4) :: nprocs  ! Directly given
      !integer(4) :: npx  ! Directly given
      !integer(4) :: npy  ! Directly given
      !integer(4) :: npz  ! Directly given
      !integer(4) :: pid  ! Directly given
      !integer(4) :: pidx ! Computed here
      !integer(4) :: pidy ! Computed here
      !integer(4) :: pidz ! Computed here
      !integer(4) :: iml  ! Computed here
      !integer(4) :: jml  ! Computed here
      !integer(4) :: kml  ! Computed here
      !integer(4) :: il  ! Computed here
      !integer(4) :: iu  ! Computed here
      !integer(4) :: jl  ! Computed here
      !integer(4) :: ju  ! Computed here
      !integer(4) :: kl  ! Computed here
      !integer(4) :: ku  ! Computed here
      !logical :: is_initialized = .false. ! Computed here

      if(gridobj%is_initialized .eqv. .false.)then
        write(*,*)"Error, grid module, grid_verbose()"
        write(*,*)"  Uninitialized grid object delivered."
        write(*,*)"  Aborting..."
        call abort()
      endif


      ! JSON compatible format
      write(*,*)'{'
      write(*,*)'"name":',trim(gridobj%name),','
      write(*,*)'"x0":',gridobj%x0,','
      write(*,*)'"y0":',gridobj%y0,','
      write(*,*)'"z0":',gridobj%z0,','
      write(*,*)'"dxl":',gridobj%dxl,','
      write(*,*)'"dyl":',gridobj%dyl,','
      write(*,*)'"dzl":',gridobj%dzl,','
      write(*,*)'"img":',gridobj%img,','
      write(*,*)'"jmg":',gridobj%jmg,','
      write(*,*)'"kmg":',gridobj%kmg,','
      write(*,*)'"nprocs":',gridobj%nprocs,','
      write(*,*)'"npx":',gridobj%npx,','
      write(*,*)'"npy":',gridobj%npy,','
      write(*,*)'"npz":',gridobj%npz,','
      write(*,*)'"pid":',gridobj%pid,','
      write(*,*)'"pidx":',gridobj%pidx,','
      write(*,*)'"pidy":',gridobj%pidy,','
      write(*,*)'"pidz":',gridobj%pidz,','
      write(*,*)'"iml":',gridobj%iml,','
      write(*,*)'"jml":',gridobj%jml,','
      write(*,*)'"kml":',gridobj%kml,','
      write(*,*)'"il":',gridobj%il,','
      write(*,*)'"iu":',gridobj%iu,','
      write(*,*)'"jl":',gridobj%jl,','
      write(*,*)'"ju":',gridobj%ju,','
      write(*,*)'"kl":',gridobj%kl,','
      write(*,*)'"ku":',gridobj%ku,','
      write(*,*)'"vectorization_ordering":',"fortran"
      write(*,*)'}'

      return
      end subroutine grid_verbose
      
      ! ----------------------------------------------------------------

      subroutine grid_update_origin(gridobj, x0_in, y0_in, z0_in)
      implicit none
      type(grid),intent(inout) :: gridobj 
      real(8),intent(in) :: x0_in 
      real(8),intent(in) :: y0_in 
      real(8),intent(in) :: z0_in 

      if(gridobj%is_initialized .eqv. .false.)then
        write(*,*)"Error, grid module, grid_update_origin()"
        write(*,*)"  Uninitialized grid object delivered."
        write(*,*)"  Aborting..."
        call abort()
      endif

      gridobj%x0 = x0_in 
      gridobj%y0 = y0_in 
      gridobj%z0 = z0_in 

      return
      end subroutine grid_update_origin

      ! ----------------------------------------------------------------

      subroutine grid_update_grid_sizes(gridobj, dxl_in, dyl_in, dzl_in)
      implicit none
      type(grid),intent(inout) :: gridobj 
      real(8),intent(in) :: dxl_in
      real(8),intent(in) :: dyl_in
      real(8),intent(in) :: dzl_in

      if(gridobj%is_initialized .eqv. .false.)then
        write(*,*)"Error, grid module, grid_update_grid_sizes()"
        write(*,*)"  Uninitialized grid object delivered."
        write(*,*)"  Aborting..."
        call abort()
      endif

      gridobj%dxl = dxl_in 
      gridobj%dyl = dyl_in 
      gridobj%dzl = dzl_in 

      return
      end subroutine grid_update_grid_sizes

      ! ----------------------------------------------------------------

      subroutine grid_copy(gridobj, gridname_target, grid_source)
      implicit none 
      type(grid),intent(inout) :: gridobj 
      character(len=*),intent(in) :: gridname_target
      type(grid),intent(in) :: grid_source

      if(grid_source%is_initialized .eqv. .false.)then
        write(*,*)"Error, grid module, grid_copy()"
        write(*,*)"  Uninitialized source grid object delivered."
        write(*,*)"  Aborting..."
        call abort()
      endif

      gridobj = grid_source 
      gridobj%name = trim(gridname_target)

      return
      end subroutine grid_copy

      ! ----------------------------------------------------------------

      subroutine initialize_grid_from_file(grid_var, filename)
      implicit none
    
      type(grid), intent(out) :: grid_var
      character(len=*), intent(in) :: filename
    
      integer :: iunit, ios
      character(len=256) :: line, key, value
    
      ! Open the file for reading
      open(newunit=iunit, file=filename, status='old', action='read', 
     &                                                       iostat=ios)
      if (ios /= 0) then
        print *, "Error: Could not open file ", filename
        return
      end if
    
      ! Read the file line by line
      do
        read(iunit, '(A)', iostat=ios) line
        if (ios /= 0) exit  ! End of file or error
    
        ! Skip lines that don't contain key-value pairs
        if (index(line, ':') == 0) cycle
    
        ! Extract key and value from the line
        call parse_line_json(line, key, value)
        !write(*,*)'key=',trim(key),', value=',trim(value)
    
        ! Assign the value to the corresponding grid field
        select case(trim(adjustl(key)))
          case('name')
            grid_var%name = trim(adjustl(value))
          case('x0')
            read(value, *) grid_var%x0
          case('y0')
            read(value, *) grid_var%y0
          case('z0')
            read(value, *) grid_var%z0
          case('dxl')
            read(value, *) grid_var%dxl
          case('dyl')
            read(value, *) grid_var%dyl
          case('dzl')
            read(value, *) grid_var%dzl
          case('img')
            read(value, *) grid_var%img
          case('jmg')
            read(value, *) grid_var%jmg
          case('kmg')
            read(value, *) grid_var%kmg
          case('nprocs')
            read(value, *) grid_var%nprocs
          case('npx')
            read(value, *) grid_var%npx
          case('npy')
            read(value, *) grid_var%npy
          case('npz')
            read(value, *) grid_var%npz
          case('pid')
            read(value, *) grid_var%pid
          case('pidx')
            read(value, *) grid_var%pidx
          case('pidy')
            read(value, *) grid_var%pidy
          case('pidz')
            read(value, *) grid_var%pidz
          case('iml')
            read(value, *) grid_var%iml
          case('jml')
            read(value, *) grid_var%jml
          case('kml')
            read(value, *) grid_var%kml
          case('il')
            read(value, *) grid_var%il
          case('iu')
            read(value, *) grid_var%iu
          case('jl')
            read(value, *) grid_var%jl
          case('ju')
            read(value, *) grid_var%ju
          case('kl')
            read(value, *) grid_var%kl
          case('ku')
            read(value, *) grid_var%ku
          case default
            ! Do nothing for unrecognized keys
        end select
      end do
    
      ! Mark grid as initialized
      grid_var%is_initialized = .true.
    
      ! Close the file
      close(iunit)
      
      return 
      end subroutine initialize_grid_from_file

      ! ----------------------------------------------------------------
      
      ! Helper subroutine to parse a line into a key and value
      subroutine parse_line_json(line, key, value)
      implicit none
      character(len=*), intent(in) :: line
      character(len=256), intent(out) :: key, value
      integer :: sep_pos
    
      ! Find the position of the colon separator
      sep_pos = index(line, ':')
    
      ! Extract the key and value
      key = line(1:sep_pos-1)
      value = line(sep_pos+1:)
    
      ! Trim quotes and spaces from key and value
      key = adjustl(trim(key))
      value = adjustl(trim(value))
    
      ! Remove any surrounding double quotes from the key
      if (key(1:1) == '"') then
        key = key(2:)  ! Remove the starting double quote
      end if
      if (key(len_trim(key):len_trim(key)) == '"') then
        key = key(:len_trim(key)-1)  ! Remove the ending double quote
      end if
    
      ! Remove any surrounding double quotes from the value (if it's a string)
      if (value(1:1) == '"') then
        value = value(2:)  ! Remove the starting double quote
      end if
      if (value(len_trim(value):len_trim(value)) == '"') then
        value = value(:len_trim(value)-1)  ! Remove the ending double quote
      end if
    
      ! Remove any trailing commas from the value
      if (len_trim(value) > 0 .and. 
     &               value(len_trim(value):len_trim(value)) == ',') then
        value = value(1:len_trim(value)-1)
      end if

      return 
      end subroutine parse_line_json 
      
      
      ! ----------------------------------------------------------------

      end module grid_module