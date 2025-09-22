      module boundary_conditions

      use mpi

      implicit none

      integer(4) :: bc_call_count
      integer(4) :: bc_max_mpi_tag = 0
      integer(4),parameter :: bc_mpitag_jumpstep = 13

      integer(4) :: bc_npx=1, bc_npy=1, bc_npz=1 !Number of partitioning along x,y, and z axis
      integer(4) :: bc_im,bc_jm,bc_km
      integer(4) :: bc_imori, bc_jmori, bc_kmori
      integer(4) :: bc_pid

      logical :: bc_initialized = .false.

      integer(4),parameter :: namelen = 30
      type :: bctype
        character(len=100) :: handlename
        character(len=namelen) :: xp
        character(len=namelen) :: xn
        character(len=namelen) :: yp
        character(len=namelen) :: yn
        character(len=namelen) :: zp
        character(len=namelen) :: zn
        character(len=namelen) :: legacy
        logical :: is_initialized

        ! Available keys for xp, xn, yp, yn, zp, zn:
        ! PERIODIC
        ! ADIABATIC
        ! NEUMANN
        ! NEUMANN_CONST
        ! DIRICHLET
        ! DIRICHLET_CONST

        ! Available keys for legacy:
        ! PERIODIC
        ! ADIABATIC
        ! ADIABATIC_X
        ! ADIABATIC_Y
        ! ADIABATIC_Z
        ! ADIABATIC_XY
        ! ADIABATIC_YZ
        ! ADIABATIC_ZX
      end type

      contains

      ! ----------------------------------------------------------------

      subroutine init_bc_grid(inpx,inpy,inpz,i_all,j_all,k_all,irank)
      implicit none
      integer(4),intent(in) :: inpx,inpy,inpz !saved in bc_npx, bc_npy, bc_npz
      integer(4),intent(in) :: i_all, j_all, k_all !saved in bc_imori,bc_jmori,bc_kmori
      integer(4),intent(in) :: irank !saved in bc_pid

      integer(4) :: ilb,iub,jlb,jub,klb,kub

      bc_npx = inpx
      bc_npy = inpy
      bc_npz = inpz

      bc_imori = i_all
      bc_jmori = j_all
      bc_kmori = k_all
      
      bc_pid = irank

      call get_ijkrange(ilb,iub,jlb,jub,klb,kub,bc_pid)
      bc_im=iub-ilb+1
      bc_jm=jub-jlb+1
      bc_km=kub-klb+1

      bc_initialized = .true.

      end subroutine init_bc_grid

      ! ----------------------------------------------------------------

      subroutine init_bc_type(bctype_handle, bc_name, 
     &                                             bctype_xn, bctype_xp, 
     &                                             bctype_yn, bctype_yp, 
     &                                             bctype_zn, bctype_zp)
      implicit none
      type(bctype),intent(inout) :: bctype_handle
      character(len=*),intent(in) :: bc_name
      character(len=*),intent(in) :: bctype_xn
      character(len=*),intent(in) :: bctype_xp
      character(len=*),intent(in) :: bctype_yn
      character(len=*),intent(in) :: bctype_yp
      character(len=*),intent(in) :: bctype_zn
      character(len=*),intent(in) :: bctype_zp
      
      logical :: is_xn_periodic
      logical :: is_xp_periodic
      logical :: is_yn_periodic
      logical :: is_yp_periodic
      logical :: is_zn_periodic
      logical :: is_zp_periodic

      bctype_handle%handlename = trim(bc_name)
      bctype_handle%is_initialized = .false. ! Not initialized yet

      ! Parity check. PERIODIC BC should be paired for both boundaries
      is_xn_periodic = trim(bctype_xn).eq.'PERIODIC'
      is_xp_periodic = trim(bctype_xp).eq.'PERIODIC'
      is_yn_periodic = trim(bctype_yn).eq.'PERIODIC'
      is_yp_periodic = trim(bctype_yp).eq.'PERIODIC'
      is_zn_periodic = trim(bctype_zn).eq.'PERIODIC'
      is_zp_periodic = trim(bctype_zp).eq.'PERIODIC'

      if(is_xn_periodic .or. is_xp_periodic) then
        if(is_xn_periodic .neqv. is_xp_periodic) then
          write(*,*)'Error, init_bc_type(), parity error'
          write(*,*)'  PERIODIC BC should be paired for both boundaries'
          write(*,*)'  BC name:', trim(bctype_handle%handlename)
          write(*,*)'  bctype_xn:',bctype_xn
          write(*,*)'  bctype_xp:',bctype_xp
          write(*,*)'  Aborting...'
          call abort()
        endif
      endif

      if(is_yn_periodic .or. is_yp_periodic) then
        if(is_yn_periodic .neqv. is_yp_periodic) then
          write(*,*)'Error, init_bc_type(), parity error'
          write(*,*)'  PERIODIC BC should be paired for both boundaries'
          write(*,*)'  BC name:', trim(bctype_handle%handlename)
          write(*,*)'  bctype_xn:',bctype_yn
          write(*,*)'  bctype_xp:',bctype_yp
          write(*,*)'  Aborting...'
          call abort()
        endif
      endif

      if(is_zn_periodic .or. is_zp_periodic) then
        if(is_zn_periodic .neqv. is_zp_periodic) then
          write(*,*)'Error, init_bc_type(), parity error'
          write(*,*)'  PERIODIC BC should be paired for both boundaries'
          write(*,*)'  BC name:', trim(bctype_handle%handlename)
          write(*,*)'  bctype_xn:',bctype_zn
          write(*,*)'  bctype_xp:',bctype_zp
          write(*,*)'  Aborting...'
          call abort()
        endif
      endif

      bctype_handle%xn = trim(bctype_xn)
      bctype_handle%xp = trim(bctype_xp)
      bctype_handle%yn = trim(bctype_yn)
      bctype_handle%yp = trim(bctype_yp)
      bctype_handle%zn = trim(bctype_zn)
      bctype_handle%zp = trim(bctype_zp)
      
      bctype_handle%is_initialized = .true.

      return
      end subroutine init_bc_type

      ! ----------------------------------------------------------------

      subroutine init_bc_type_legacy(bctype_handle, bc_name,
     &                                               bctype_legacy,ndim)
      implicit none
      type(bctype),intent(inout) :: bctype_handle
      character(len=*),intent(in) :: bc_name
      character(len=*),intent(in) :: bctype_legacy
      integer(4),intent(in) :: ndim

      bctype_handle%handlename = trim(bc_name)
      bctype_handle%is_initialized = .false.

      if (ndim .eq. 2) then
        select case(trim(bctype_legacy))
        case('PERIODIC')
          bctype_handle%xn = 'PERIODIC'
          bctype_handle%xp = 'PERIODIC'
          bctype_handle%yn = 'NA'
          bctype_handle%yp = 'NA'
          bctype_handle%zn = 'PERIODIC'
          bctype_handle%zp = 'PERIODIC'
        case('ADIABATIC')
          bctype_handle%xn = 'ADIABATIC'
          bctype_handle%xp = 'ADIABATIC'
          bctype_handle%yn = 'NA'
          bctype_handle%yp = 'NA'
          bctype_handle%zn = 'ADIABATIC'
          bctype_handle%zp = 'ADIABATIC'
        case('ADIABATIC_X')
          bctype_handle%xn = 'ADIABATIC'
          bctype_handle%xp = 'ADIABATIC'
          bctype_handle%yn = 'NA'
          bctype_handle%yp = 'NA'
          bctype_handle%zn = 'PERIODIC'
          bctype_handle%zp = 'PERIODIC'
        case('ADIABATIC_Z')
          bctype_handle%xn = 'PERIODIC'
          bctype_handle%xp = 'PERIODIC'
          bctype_handle%yn = 'NA'
          bctype_handle%yp = 'NA'
          bctype_handle%zn = 'ADIABATIC'
          bctype_handle%zp = 'ADIABATIC'
        case default
          write(*,*)'Error, init_bc_type_legacy(), '
          write(*,*)'  Unknown legacy BC type: ',trim(bctype_legacy)
          write(*,*)'  bc name: ', trim(bctype_handle%handlename)
          write(*,*)'  Aborting...'
          call abort()
        end select

      elseif (ndim .eq. 3) then
        select case(trim(bctype_legacy))
        case('PERIODIC')
          bctype_handle%xn = 'PERIODIC'
          bctype_handle%xp = 'PERIODIC'
          bctype_handle%yn = 'PERIODIC'
          bctype_handle%yp = 'PERIODIC'
          bctype_handle%zn = 'PERIODIC'
          bctype_handle%zp = 'PERIODIC'
        case('ADIABATIC')
          bctype_handle%xn = 'ADIABATIC'
          bctype_handle%xp = 'ADIABATIC'
          bctype_handle%yn = 'ADIABATIC'
          bctype_handle%yp = 'ADIABATIC'
          bctype_handle%zn = 'ADIABATIC'
          bctype_handle%zp = 'ADIABATIC'
        case('ADIABATIC_X')
          bctype_handle%xn = 'ADIABATIC'
          bctype_handle%xp = 'ADIABATIC'
          bctype_handle%yn = 'PERIODIC'
          bctype_handle%yp = 'PERIODIC'
          bctype_handle%zn = 'PERIODIC'
          bctype_handle%zp = 'PERIODIC'
        case('ADIABATIC_Y')
          bctype_handle%xn = 'PERIODIC'
          bctype_handle%xp = 'PERIODIC'
          bctype_handle%yn = 'ADIABATIC'
          bctype_handle%yp = 'ADIABATIC'
          bctype_handle%zn = 'PERIODIC'
          bctype_handle%zp = 'PERIODIC'
        case('ADIABATIC_Z')
          bctype_handle%xn = 'PERIODIC'
          bctype_handle%xp = 'PERIODIC'
          bctype_handle%yn = 'PERIODIC'
          bctype_handle%yp = 'PERIODIC'
          bctype_handle%zn = 'ADIABATIC'
          bctype_handle%zp = 'ADIABATIC'
        case('ADIABATIC_XY')
          bctype_handle%xn = 'ADIABATIC'
          bctype_handle%xp = 'ADIABATIC'
          bctype_handle%yn = 'ADIABATIC'
          bctype_handle%yp = 'ADIABATIC'
          bctype_handle%zn = 'PERIODIC'
          bctype_handle%zp = 'PERIODIC'
        case('ADIABATIC_YZ')
          bctype_handle%xn = 'PERIODIC'
          bctype_handle%xp = 'PERIODIC'
          bctype_handle%yn = 'ADIABATIC'
          bctype_handle%yp = 'ADIABATIC'
          bctype_handle%zn = 'ADIABATIC'
          bctype_handle%zp = 'ADIABATIC'
        case('ADIABATIC_ZX')
          bctype_handle%xn = 'ADIABATIC'
          bctype_handle%xp = 'ADIABATIC'
          bctype_handle%yn = 'PERIODIC'
          bctype_handle%yp = 'PERIODIC'
          bctype_handle%zn = 'ADIABATIC'
          bctype_handle%zp = 'ADIABATIC'
        case default
          write(*,*)'Error, init_bc_type_legacy(), '
          write(*,*)'  Unknown legacy BC type: ',trim(bctype_legacy)
          write(*,*)'  bc name: ', trim(bctype_handle%handlename)
          write(*,*)'  Aborting...'
          call abort()
        end select

      endif !if (ndim .eq. 2) then ... elseif (ndim .eq. 3) then

      bctype_handle%legacy = trim(bctype_legacy)

      bctype_handle%is_initialized = .true.

      return
      end subroutine init_bc_type_legacy

      ! ----------------------------------------------------------------

      subroutine init_bc_call_count(val)
      implicit none
      integer(4),intent(in) :: val

       bc_call_count = val

      end subroutine init_bc_call_count

      ! ----------------------------------------------------------------

      subroutine verbose_bc_call_count()
      implicit none
      
        write(*,*) 'bc_call_count:', bc_call_count

      end subroutine verbose_bc_call_count

      ! ----------------------------------------------------------------

      subroutine verbose_bc_max_mpi_tag()
      implicit none

      write(*,*) 'bc_max_mpi_tag:', bc_max_mpi_tag

      end subroutine verbose_bc_max_mpi_tag

      ! ----------------------------------------------------------------

      subroutine monitor_bc_max_mpi_tag()
      implicit none
      if(bc_max_mpi_tag.ge.32767)then
            write(*,*)'Reminder: bc_max_mpi_tag exceeded 32767.',
     &           bc_max_mpi_tag
            if(bc_max_mpi_tag.ge.MPI_TAG_UB)then
            write(*,*)'ERROR: bc_max_mpi_tag exceeded MPI_TAG_UB.',
     &        bc_max_mpi_tag
            write(*,*)'MPI_TAG_UB for this system:',MPI_TAG_UB
            write(*,*)'Aborting due to high bc_max_mpi_tag.'
            call abort
            endif
      endif
      end subroutine monitor_bc_max_mpi_tag

      ! ----------------------------------------------------------------

      subroutine bound_check(flagid, i,j,k, ng)
      implicit none
      character(len=*),intent(in) :: flagid
      integer(4),intent(in) :: i,j,k, ng

      logical :: go_nogo
      character(len=100) :: err_sign

      if ( bc_initialized .eqv. .false.) then
            write(*,*)'Error, bound_check()'
            write(*,*)'  boundary condition grid uninitialized.'
            write(*,*)'  Aborting...'
            call abort
      endif


      err_sign=" "
      go_nogo = .true.

      if(i.lt.1-ng) then
            go_nogo = .false.
            err_sign=trim(err_sign)//"(iL)"
      endif
      if(i.gt.bc_im+ng) then
            go_nogo = .false.
            err_sign=trim(err_sign)//"(iU)"
      endif
      if(j.lt.1-ng) then
            go_nogo = .false.
            err_sign=trim(err_sign)//"(jL)"
      endif
      if(j.gt.bc_jm+ng) then
            go_nogo = .false.
            err_sign=trim(err_sign)//"(jU)"
      endif
      if(k.lt.1-ng) then
            go_nogo = .false.
            err_sign=trim(err_sign)//"(kL)"
      endif
      if(k.gt.bc_km+ng) then
            go_nogo = .false.
            err_sign=trim(err_sign)//"(kU)"
      endif

      if(go_nogo.eqv. .false.) then
            write(*,*)"Tried to access the out of array bound:"
            write(*,*)"flagid=",trim(flagid),", i,j,k = ",i,j,k
            write(*,*)"err_sign: ",trim(err_sign)," ,rank=",bc_pid
            call abort
      endif

      end subroutine bound_check

      ! ----------------------------------------------------------------

      integer(4) function get_pid(pidx, pidy, pidz)
      implicit none
      integer(4),intent(in) :: pidx, pidy, pidz

      integer(4) :: pidxw, pidyw, pidzw

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, get_pid()'
        write(*,*)'  boundary condition grid uninitialized.'
        write(*,*)'  Aborting...'
        call abort
      endif


      pidxw = pidx
      if(pidxw.gt.bc_npx-1)pidxw=0
      if(pidxw.lt.0)pidxw=bc_npx-1

      pidyw = pidy
      if(pidyw.gt.bc_npy-1)pidyw=0
      if(pidyw.lt.0)pidyw=bc_npy-1

      pidzw = pidz
      if(pidzw.gt.bc_npz-1)pidzw=0
      if(pidzw.lt.0)pidzw=bc_npz-1

      get_pid = pidxw + pidyw*bc_npx + pidzw*bc_npx*bc_npy
      return
      end function get_pid

      ! ----------------------------------------------------------------

      subroutine get_pidxyz(pidx, pidy, pidz, pid)
      implicit none
      integer(4),intent(inout) :: pidx, pidy, pidz
      integer(4),intent(in) :: pid

      pidz = int(pid/bc_npx/bc_npy)
      pidy = int((pid - pidz*bc_npx*bc_npy)/bc_npx)
      pidx = pid - pidz*bc_npx*bc_npy - pidy*bc_npx

      end subroutine get_pidxyz

      ! ----------------------------------------------------------------

      subroutine get_ijkrange(il,iu,jl,ju,kl,ku,pid)
      implicit none
      integer(4),intent(inout) :: il,iu,jl,ju,kl,ku
      integer(4),intent(in) :: pid

      integer(4) :: pidx, pidy, pidz


      call get_pidxyz(pidx, pidy, pidz, pid)

      ! x bounds
      call para_range(1,bc_imori,bc_npx,pidx,il,iu)
      ! y bounds
      call para_range(1,bc_jmori,bc_npy,pidy,jl,ju)
      ! z bounds
      call para_range(1,bc_kmori,bc_npz,pidz,kl,ku)

      end subroutine get_ijkrange

      ! ----------------------------------------------------------------

      subroutine get_global_ijk(ig,jg,kg, ilocal,jlocal,klocal, pid)
      implicit none
      integer(4),intent(inout) :: ig, jg, kg !global ijk index
      integer(4),intent(in) :: ilocal, jlocal, klocal !local ijk index
      integer(4),intent(in) :: pid

      integer(4) :: il,iu, jl,ju, kl,ku

      call get_ijkrange(il,iu,jl,ju,kl,ku,pid)
      ig = il + (ilocal-1)  !for spatial index starts from 1, not 0
      jg = jl + (jlocal-1)  !for spatial index starts from 1, not 0
      kg = kl + (klocal-1)  !for spatial index starts from 1, not 0

      end subroutine get_global_ijk

      ! ----------------------------------------------------------------

      SUBROUTINE para_range(n1,n2,npart,irank,ista,iend)
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

      END SUBROUTINE para_range

      ! ----------------------------------------------------------------

      integer(4) function get_pid_from_globalijk(ig,jg,kg,bc)
      implicit none
      integer(4),intent(in) :: ig,jg,kg
      !character(len=*),intent(in) :: bctype
      type(bctype),intent(in) :: bc
      
      integer(4) :: igf, jgf, kgf
      integer(4) :: pidx, pidy, pidz

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, get_pid_from_globalijk()'
        write(*,*)'  boundary condition grid uninitialized.'
        write(*,*)'  Aborting...'
        call abort
      endif

      if ( bc%is_initialized .eqv. .false.) then
        write(*,*)'Error, get_pid_from_globalijk()'
        write(*,*)'  bc handle uninitialized.'
        write(*,*)'  bc name:',trim(bc%handlename)
        write(*,*)'  Aborting...'
        call abort
      endif


      call clamp_globalijk(igf,jgf,kgf,ig,jg,kg,bc)

      pidx = int((igf-1)/bc_im)
      pidy = int((jgf-1)/bc_jm)
      pidz = int((kgf-1)/bc_km)

      get_pid_from_globalijk = get_pid(pidx, pidy, pidz)
      return
      end function get_pid_from_globalijk

      ! ----------------------------------------------------------------

      subroutine get_localijk_from_globalijk(il,jl,kl,ig,jg,kg,bc)
      implicit none
      integer(4),intent(inout) :: il,jl,kl
      integer(4),intent(in) :: ig,jg,kg
      !character(len=*) :: bctype
      type(bctype),intent(in) :: bc

      integer(4) :: igf, jgf, kgf
      integer(4) :: pid, pidx, pidy, pidz

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, get_localijk_from_globalijk()'
        write(*,*)'  boundary condition grid uninitialized.'
        write(*,*)'  Aborting...'
        call abort
      endif

      if ( bc%is_initialized .eqv. .false.) then
        write(*,*)'Error, get_localijk_from_globalijk()'
        write(*,*)'  bc handle uninitialized.'
        write(*,*)'  bc name:',trim(bc%handlename)
        write(*,*)'  Aborting...'
        call abort
      endif


      pid = get_pid_from_globalijk(ig,jg,kg,bc)
      call get_pidxyz(pidx, pidy, pidz, pid)

      call clamp_globalijk(igf,jgf,kgf,ig,jg,kg,bc)

      il=igf-pidx*bc_im
      jl=jgf-pidy*bc_jm
      kl=kgf-pidz*bc_km

      end subroutine get_localijk_from_globalijk

      ! ----------------------------------------------------------------

      subroutine clamp_globalijk(igo,jgo,kgo,ig,jg,kg,bc)
      integer(4),intent(inout) :: igo, jgo, kgo
      integer(4),intent(in) :: ig,jg,kg
      !character(len=*) :: bctype
      type(bctype),intent(in) :: bc

      if ( bc_initialized .eqv. .false.) then
            write(*,*)'Error, clamp_globalijk()'
            write(*,*)'  boundary condition grid uninitialized.'
            write(*,*)'  Aborting...'
            call abort
      endif

      if ( bc%is_initialized .eqv. .false.) then
        write(*,*)'Error, clamp_globalijk()'
        write(*,*)'  bc handle uninitialized.'
        write(*,*)'  bc name:',trim(bc%handlename)
        write(*,*)'  Aborting...'
        call abort
      endif


      igo=ig
      jgo=jg
      kgo=kg

      if (igo.lt.1 .or. igo.gt.bc_imori ) then
        if(trim(bc%xn).eq.'PERIODIC' .or. 
     &                               trim(bc%xp).eq.'PERIODIC') then
          igo=mod(ig+min(ig-1,0)*(-bc_imori)-1,bc_imori)+1
        else
          igo=min(max(ig,1),bc_imori)
        endif
      endif !if (igo.lt.1 .or. igo.gt.bc_imori ) then

      if (jgo.lt.1 .or. jgo.gt.bc_jmori ) then
        if(trim(bc%yn).eq.'PERIODIC' .or. 
     &                               trim(bc%yp).eq.'PERIODIC') then
          jgo=mod(jg+min(jg-1,0)*(-bc_jmori)-1,bc_jmori)+1
        else
          jgo=min(max(jg,1),bc_jmori)
        endif
      endif !if (jgo.lt.1 .or. jgo.gt.bc_jmori ) then

      if (kgo.lt.1 .or. kgo.gt.bc_kmori ) then
        if(trim(bc%zn).eq.'PERIODIC' .or. 
     &                               trim(bc%zp).eq.'PERIODIC') then
          kgo=mod(kg+min(kg-1,0)*(-bc_kmori)-1,bc_kmori)+1
        else
          kgo=min(max(kg,1),bc_kmori)
        endif
      endif !if (kgo.lt.1 .or. kgo.gt.bc_kmori ) then

      end subroutine clamp_globalijk

      ! ----------------------------------------------------------------

      subroutine boundary_condition_2d_i4(bc,a,ng)
      implicit none

      !Subroutine arguments
      !character(len=*),intent(in) :: bctype
      type(bctype),intent(in) :: bc
      integer(4),intent(inout) :: a(1-ng:bc_im+ng,1,1-ng:bc_km+ng)
      integer(4),intent(in) :: ng

      !Local variables
      integer(4) :: i,j,k,n1,n2
      integer(4) :: il,iu,jl,ju,kl,ku

      ! MPI parallelization variables
      integer(4) :: pidx, pidy, pidz
      integer(4) :: pidto01, pidto02, pidto03, pidto04
      integer(4) :: pidto05, pidto06, pidto07, pidto08
      integer(4) :: pidto09, pidto10, pidto11, pidto12
      integer(4) :: pidfrom01, pidfrom02, pidfrom03, pidfrom04
      integer(4) :: pidfrom05, pidfrom06, pidfrom07, pidfrom08
      integer(4) :: pidfrom09, pidfrom10, pidfrom11, pidfrom12

      integer(4) :: isend01, isend02, isend03, isend04
      integer(4) :: isend05, isend06, isend07, isend08
      integer(4) :: isend09, isend10, isend11, isend12
      integer(4) :: irecv01, irecv02, irecv03, irecv04
      integer(4) :: irecv05, irecv06, irecv07, irecv08
      integer(4) :: irecv09, irecv10, irecv11, irecv12


      !MPI_HANDLES for MPI types for non-contiguous array data communication
      ! MPI type is a kind of MPI_HANDLE in C, but it is equivalent to integer in fortran
      integer(4) :: xslab, yslab, zslab

      integer(4) :: arrsize(3)
      integer(4) :: subsize(3)
      integer(4) :: starts(3)

      integer(4) :: istatus(MPI_STATUS_SIZE), ierr
      integer(4) :: basetag

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, boundary_condition_2d_i4(),'
        write(*,*)'  boundary condition grid uninitialized. Aborting.'
        call abort
      endif

      if (bc%is_initialized .eqv. .false.) then
        write(*,*)'Error, boundary_condition_2d_i4(),'
        write(*,*)'  bc handle uninitialized.'
        write(*,*)'  bc name:',trim(bc%handlename)
        write(*,*)'  Aborting...'
        call abort
      endif


      ! Resister the call count to determine the base MPI_TAG value
      bc_call_count = 1 + bc_call_count
      basetag = bc_call_count*bc_mpitag_jumpstep
      bc_max_mpi_tag = max(bc_max_mpi_tag,basetag+bc_mpitag_jumpstep)


      arrsize(1) = 2*ng + bc_im
      arrsize(2) = 1
      arrsize(3) = 2*ng + bc_km

      il = 1-ng
      iu = bc_im+ng
      jl = 1
      ju = 1
      kl = 1-ng
      ku = bc_km+ng

      ! MPI derived data type definition for non-cnotiguous transfer
      ! xslab
      starts(1) = 0
      starts(2) = 0
      starts(3) = 0

      subsize(1) = ng
      subsize(2) = arrsize(2)
      subsize(3) = arrsize(3)
      call mpi_type_create_subarray(3, arrsize, subsize, starts,
     &  MPI_ORDER_FORTRAN, MPI_INTEGER, xslab, ierr)
      call mpi_type_commit(xslab, ierr)

      ! zslab
      subsize(1) = arrsize(1)
      subsize(2) = arrsize(2)
      subsize(3) = ng
      call mpi_type_create_subarray(3, arrsize, subsize, starts,
     &  MPI_ORDER_FORTRAN, MPI_INTEGER, zslab, ierr)
      call mpi_type_commit(zslab, ierr)

      select case(trim(bc%legacy))
      case('PERIODIC')
        !Process coordinate
        call get_pidxyz(pidx, pidy, pidz, bc_pid)

        ! Layer send/recv
        ! i planes
        if(bc_npx.gt.1)then
          pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
          call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)

          pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
          call MPI_IRECV(a(bc_im+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)

          pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
          call MPI_ISEND(a(bc_im-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      

          pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
          call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
          
          call MPI_WAIT(isend01,istatus,ierr)
          call MPI_WAIT(isend02,istatus,ierr)
          call MPI_WAIT(irecv01,istatus,ierr)
          call MPI_WAIT(irecv02,istatus,ierr)
        endif
        ! Done only for x direction - not complete yet

        ! k planes
        if(bc_npz.gt.1)then
          pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
          call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)

          pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
          call MPI_IRECV(a(il,jl,bc_km+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
        
          pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
          call MPI_ISEND(a(il,jl,bc_km-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)

          pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
          call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)

          call MPI_WAIT(isend05,istatus,ierr)
          call MPI_WAIT(isend06,istatus,ierr)
          call MPI_WAIT(irecv05,istatus,ierr)
          call MPI_WAIT(irecv06,istatus,ierr)
        endif

        call apply_periodic_bc_2d_i4(a, ng)
        ! Done for all directions - complete

      case('ADIABATIC')
        ! Necessary data transfer among processes
        !Process coordinate
        call get_pidxyz(pidx, pidy, pidz, bc_pid)

        ! Layer send/recv
        ! i planes
        if(pidx.ge.0+1)then
          pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
          call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
        endif

        if(pidx.le.bc_npx-1-1)then
          pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
          call MPI_IRECV(a(bc_im+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
        endif

        if(pidx.le.bc_npx-1-1)then
          pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
          call MPI_ISEND(a(bc_im-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
        endif

        if(pidx.ge.0+1)then
          pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
          call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
        endif
        
        if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
        if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
        if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
        if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
        ! Done only for x direction - not complete yet

        ! k planes
        if(pidz.ge.0+1)then
          pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
          call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
        endif

        if(pidz.le.bc_npz-1-1)then
          pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
          call MPI_IRECV(a(il,jl,bc_km+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
        endif
        
        if(pidz.le.bc_npz-1-1)then
          pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
          call MPI_ISEND(a(il,jl,bc_km-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
        endif

        if(pidz.ge.0+1)then
          pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
          call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
        endif

        if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
        if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
        if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
        if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
        ! Done for x and y direction - complete

        ! Apply ADIABATIC BC
        call apply_nonperiodic_bc_2d_i4(bc, a, ng)

      case('ADIABATIC_X')
        ! Necessary data transfer among processes
        !Process coordinate
        call get_pidxyz(pidx, pidy, pidz, bc_pid)

        ! Layer send/recv
        ! i planes
        if(pidx.ge.0+1)then
          pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
          call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
        endif

        if(pidx.le.bc_npx-1-1)then
          pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
          call MPI_IRECV(a(bc_im+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
        endif

        if(pidx.le.bc_npx-1-1)then
          pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
          call MPI_ISEND(a(bc_im-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
        endif

        if(pidx.ge.0+1)then
          pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
          call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
        endif
        
        if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
        if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
        if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
        if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
        ! Done only for x direction - not complete yet

        ! k planes
        if(.true.)then
          pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
          call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
        endif

        if(.true.)then
          pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
          call MPI_IRECV(a(il,jl,bc_km+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
        endif
        
        if(.true.)then
          pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
          call MPI_ISEND(a(il,jl,bc_km-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
        endif

        if(.true.)then
          pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
          call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
        endif

        if(.true.) call MPI_WAIT(isend05,istatus,ierr)
        if(.true.) call MPI_WAIT(isend06,istatus,ierr)
        if(.true.) call MPI_WAIT(irecv05,istatus,ierr)
        if(.true.) call MPI_WAIT(irecv06,istatus,ierr)
        ! Done for x and y direction - complete

        ! Apply ADIABATIC_X BC
        call apply_nonperiodic_bc_2d_i4(bc, a, ng)

      case('ADIABATIC_Z')
        ! Necessary data transfer among processes
        !Process coordinate
        call get_pidxyz(pidx, pidy, pidz, bc_pid)

        ! Layer send/recv
        ! i planes
        if(.true.)then
          pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
          call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
        endif

        if(.true.)then
          pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
          call MPI_IRECV(a(bc_im+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
        endif

        if(.true.)then
          pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
          call MPI_ISEND(a(bc_im-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
        endif

        if(.true.)then
          pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
          call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
        endif
        
        if(.true.) call MPI_WAIT(isend01,istatus,ierr)
        if(.true.) call MPI_WAIT(isend02,istatus,ierr)
        if(.true.) call MPI_WAIT(irecv01,istatus,ierr)
        if(.true.) call MPI_WAIT(irecv02,istatus,ierr)
        ! Done only for x direction - not complete yet

        ! k planes
        if(pidz.ge.0+1)then
          pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
          call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
        endif

        if(pidz.le.bc_npz-1-1)then
          pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
          call MPI_IRECV(a(il,jl,bc_km+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
        endif
        
        if(pidz.le.bc_npz-1-1)then
          pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
          call MPI_ISEND(a(il,jl,bc_km-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
        endif

        if(pidz.ge.0+1)then
          pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
          call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
        endif

        if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
        if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
        if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
        if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
        ! Done for x and y direction - complete

        ! Apply ADIABATIC_Z BC
        call apply_nonperiodic_bc_2d_i4(bc, a, ng)

      end select !select case(trim(bc))

      call mpi_barrier(MPI_COMM_WORLD, ierr)
      call mpi_type_free(xslab, ierr)
      call mpi_type_free(zslab, ierr)

      end subroutine boundary_condition_2d_i4

      ! ----------------------------------------------------------------

      subroutine apply_periodic_bc_2d_i4(a,ng)
      implicit none
      integer(4),intent(inout) :: a(1-ng:bc_im+ng,1,1-ng:bc_km+ng)
      integer(4),intent(in) :: ng

      integer(4) :: i,j,k,lg

      if(bc_npx.gt.1.and.bc_npz.gt.1)then
            return
      endif

      if(bc_npx.eq.1)then
            do k=1-ng,bc_km+ng
                  do lg=1,ng
                        a(-ng+lg,1,k) = a(bc_im-ng+lg,1,k)
                        a(bc_im+lg,1,k) = a(lg,1,k)
                  enddo
            enddo
      endif

      if(bc_npz.eq.1)then
            do i=1-ng,bc_im+ng
                  do lg=1,ng
                        a(i,1,-ng+lg) = a(i,1,bc_km-ng+lg)
                        a(i,1,bc_km+lg) = a(i,1,lg)
                  enddo
            enddo
      endif

      end subroutine apply_periodic_bc_2d_i4

      ! ----------------------------------------------------------------

      subroutine apply_nonperiodic_bc_2d_i4(bc, a, ng)
      implicit none
      !Subroutine arguments
      !character(len=*),intent(in) :: bctype
      type(bctype),intent(in) :: bc
      integer(4),intent(inout) :: a(1-ng:bc_im+ng,1,1-ng:bc_km+ng)
      integer(4),intent(in) :: ng

      !Local variables
      integer(4) :: i,j,k,lg
      integer(4) :: l1, l2, lh
      logical :: gconxn, gconzn
      logical :: gconxp, gconzp
      logical :: lconil, lconiu
      logical :: lconkl, lconku
      integer(4) :: pidx,pidy,pidz

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, apply_nonperiodic_bc_2d_i4(),'
        write(*,*)'  boundary condition grid uninitialized. Aborting.'
        call abort
      endif

      if (bc%is_initialized .eqv. .false.) then
        write(*,*)'Error, apply_nonperiodic_bc_2d_i4(),'
        write(*,*)'  bc handle uninitialized.'
        write(*,*)'  bc name:',trim(bc%handlename)
        write(*,*)'  Aborting...'
        call abort
      endif


      gconxn=.false.; gconzn=.false.
      gconxp=.false.; gconzp=.false.
      lconil=.false.; lconkl=.false.
      lconiu=.false.; lconku=.false.

      if(trim(bc%xn) .ne. 'PERIODIC') gconxn = .true.
      if(trim(bc%xp) .ne. 'PERIODIC') gconxp = .true.
      if(trim(bc%zn) .ne. 'PERIODIC') gconzn = .true.
      if(trim(bc%zp) .ne. 'PERIODIC') gconzp = .true.

      call get_pidxyz(pidx, pidy, pidz, bc_pid)
      if(pidx.eq.0)     lconil = .true.
      if(pidx.eq.bc_npx-1) lconiu = .true.
      if(pidz.eq.0)     lconkl = .true.
      if(pidz.eq.bc_npz-1) lconku = .true.

      ! Combining global condition and local condition
      lconil=lconil.and.gconxn
      lconiu=lconiu.and.gconxp
      lconkl=lconkl.and.gconzn
      lconku=lconku.and.gconzp

      ! Terminate if this local process is nothing to do with BC
      if(lconil.or.lconiu.or.lconkl.or.lconku
     &    .eqv. .false.) then
            return
      endif
 
      lh = max(bc_im,1,bc_km)  ! The highest # of grids

      do l1=1-ng,lh+ng

      ! x-layers: l1=j=1, l2=k
      if(lconil.or.lconiu .eqv. .true.)then
      if(l1.le.bc_km+ng)then  
        do lg=1,ng
          ! il
          if(lconil.eqv..true.)a(1-lg,1,l1)=a(lg,1,l1)
          ! iu
          if(lconiu.eqv..true.)a(bc_im+lg,1,l1)=a(bc_im+1-lg,1,l1)
        enddo
      endif
      endif !if(lconil.or.lconiu .eqv. .true.)then

      ! z-layers: l1=i, l2=j=1
      if(lconkl.or.lconku .eqv. .true.)then
      if(l1.le.bc_im+ng)then  
        do lg=1,ng
          ! kl
          if(lconkl.eqv..true.)a(l1,1,1-lg)=a(l1,1,lg)
          ! ku
          if(lconku.eqv..true.)a(l1,1,bc_km+lg)=a(l1,1,bc_km+1-lg)
        enddo
      endif
      endif !if(lconkl.or.lconku .eqv. .true.)then
      

      ! Edges along y-axis: l2=j=1, l1=k, lg=i
      ! klil
      if(lconkl.and.lconil .eqv. .true.) then
            if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("klil",1-ng+l1,1,1-ng+lg,ng)
                  a(1-ng+l1,1,1-ng+lg)=a(ng-lg,1,ng-l1)
            enddo
            endif
      endif
      ! kuil
      if(lconku.and.lconil .eqv. .true.) then
            if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("kuil",bc_im+1+l1,1,bc_km+1+lg,ng)
                  a(bc_im+1+l1,1,bc_km+1+lg)=a(bc_im-lg,1,bc_km-l1)
            enddo
            endif
      endif
      ! kliu
      if(lconkl.and.lconiu .eqv. .true.) then
            if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("kliu",1-ng+l1,1,bc_km+1+lg,ng)
                  a(1-ng+l1,1,bc_km+1+lg)=a(1+lg,1,bc_km+1-ng+l1)
            enddo
            endif
      endif
      ! kuiu
      if(lconku.and.lconiu .eqv. .true.) then
            if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("kuiu",bc_im+1+l1,1,bc_km+1+lg,ng)
                  a(bc_im+1+l1,1,bc_km+1+lg)=a(bc_im-lg,1,bc_km-l1)
            enddo
            endif
      endif
      
      enddo !do l1=1-ng,lh+ng

      end subroutine apply_nonperiodic_bc_2d_i4

      ! ----------------------------------------------------------------

      subroutine boundary_condition_i4(bc,a,ng)
      implicit none

      !Subroutine arguments
      !character(len=*),intent(in) :: bctype
      type(bctype),intent(in) :: bc
      integer(4),intent(inout) :: 
     &                      a(1-ng:bc_im+ng,1-ng:bc_jm+ng,1-ng:bc_km+ng)
      integer(4),intent(in) :: ng

      !Local variables
      integer(4) :: i,j,k,n1,n2
      integer(4) :: il,iu,jl,ju,kl,ku

      ! MPI parallelization variables
      integer(4) :: pidx, pidy, pidz
      integer(4) :: pidto01, pidto02, pidto03, pidto04
      integer(4) :: pidto05, pidto06, pidto07, pidto08
      integer(4) :: pidto09, pidto10, pidto11, pidto12
      integer(4) :: pidfrom01, pidfrom02, pidfrom03, pidfrom04
      integer(4) :: pidfrom05, pidfrom06, pidfrom07, pidfrom08
      integer(4) :: pidfrom09, pidfrom10, pidfrom11, pidfrom12

      integer(4) :: isend01, isend02, isend03, isend04
      integer(4) :: isend05, isend06, isend07, isend08
      integer(4) :: isend09, isend10, isend11, isend12
      integer(4) :: irecv01, irecv02, irecv03, irecv04
      integer(4) :: irecv05, irecv06, irecv07, irecv08
      integer(4) :: irecv09, irecv10, irecv11, irecv12


      !MPI_HANDLES for MPI types for non-contiguous array data communication
      ! MPI type is a kind of MPI_HANDLE in C, but it is equivalent to integer in fortran
      integer(4) :: xslab, yslab, zslab

      integer(4) :: arrsize(3)
      integer(4) :: subsize(3)
      integer(4) :: starts(3)

      integer(4) :: istatus(MPI_STATUS_SIZE), ierr
      integer(4) :: basetag

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, boundary_condition_i4(),'
        write(*,*)'  boundary condition grid uninitialized. Aborting.'
        call abort
      endif

      if (bc%is_initialized .eqv. .false.) then
        write(*,*)'Error, boundary_condition_i4(),'
        write(*,*)'  bc handle uninitialized.'
        write(*,*)'  bc name:',trim(bc%handlename)
        write(*,*)'  Aborting...'
        call abort
      endif


      ! Resister the call count to determine the base MPI_TAG value
      bc_call_count = 1 + bc_call_count
      basetag = bc_call_count*bc_mpitag_jumpstep
      bc_max_mpi_tag = max(bc_max_mpi_tag,basetag+bc_mpitag_jumpstep)


      arrsize(1) = 2*ng + bc_im
      arrsize(2) = 2*ng + bc_jm
      arrsize(3) = 2*ng + bc_km

      il = 1-ng
      iu = bc_im+ng
      jl = 1-ng
      ju = bc_jm+ng
      kl = 1-ng
      ku = bc_km+ng

      ! MPI derived data type definition for non-cnotiguous transfer
      ! xslab
      starts(1) = 0
      starts(2) = 0
      starts(3) = 0

      subsize(1) = ng
      subsize(2) = arrsize(2)
      subsize(3) = arrsize(3)
      call mpi_type_create_subarray(3, arrsize, subsize, starts,
     &  MPI_ORDER_FORTRAN, MPI_INTEGER, xslab, ierr)
      call mpi_type_commit(xslab, ierr)

      ! yslab
      subsize(1) = arrsize(1)
      subsize(2) = ng
      subsize(3) = arrsize(3)
      call mpi_type_create_subarray(3, arrsize, subsize, starts,
     &  MPI_ORDER_FORTRAN, MPI_INTEGER, yslab, ierr)
      call mpi_type_commit(yslab, ierr)

      ! zslab
      subsize(1) = arrsize(1)
      subsize(2) = arrsize(2)
      subsize(3) = ng
      call mpi_type_create_subarray(3, arrsize, subsize, starts,
     &  MPI_ORDER_FORTRAN, MPI_INTEGER, zslab, ierr)
      call mpi_type_commit(zslab, ierr)



      select case(trim(bc%legacy))
      case('PERIODIC')
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            
            ! Layer send/recv
            ! i planes
            if(bc_npx.gt.1)then
            pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
            call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)

            pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
            call MPI_IRECV(a(bc_im+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)

            pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
            call MPI_ISEND(a(bc_im-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      

            pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
            call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            
            call MPI_WAIT(isend01,istatus,ierr)
            call MPI_WAIT(isend02,istatus,ierr)
            call MPI_WAIT(irecv01,istatus,ierr)
            call MPI_WAIT(irecv02,istatus,ierr)
            endif
            ! Done only for x direction - not complete yet

            ! j planes
            if(bc_npy.gt.1)then
            pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
            call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)

            pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
            call MPI_IRECV(a(il,bc_jm+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)

            pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
            call MPI_ISEND(a(il,bc_jm-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)

            pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
            call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)

            call MPI_WAIT(isend03,istatus,ierr)
            call MPI_WAIT(isend04,istatus,ierr)
            call MPI_WAIT(irecv03,istatus,ierr)
            call MPI_WAIT(irecv04,istatus,ierr)
            endif
            ! Done for x and y direction - not completed yet

            ! k planes
            if(bc_npz.gt.1)then
            pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
            call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)

            pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
            call MPI_IRECV(a(il,jl,bc_km+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            
            pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
            call MPI_ISEND(a(il,jl,bc_km-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)

            pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
            call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)

            call MPI_WAIT(isend05,istatus,ierr)
            call MPI_WAIT(isend06,istatus,ierr)
            call MPI_WAIT(irecv05,istatus,ierr)
            call MPI_WAIT(irecv06,istatus,ierr)
            endif

            call apply_periodic_bc_i4(a, ng)
            ! Done for all directions - complete

      case('ADIABATIC')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(pidx.ge.0+1)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(bc_im+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(bc_im-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(pidx.ge.0+1)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
            if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(pidy.ge.0+1)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,bc_jm+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,bc_jm-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(pidy.ge.0+1)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(pidy.ge.0+1)     call MPI_WAIT(isend03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(isend04,istatus,ierr)
            if(pidy.ge.0+1)     call MPI_WAIT(irecv03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(pidz.ge.0+1)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(pidz.le.bc_npz-1-1)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,bc_km+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(pidz.le.bc_npz-1-1)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,bc_km-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(pidz.ge.0+1)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
            if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC BC
            call apply_nonperiodic_bc_i4(bc, a, ng)

      case('ADIABATIC_XY')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(pidx.ge.0+1)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(bc_im+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(bc_im-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(pidx.ge.0+1)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
            if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(pidy.ge.0+1)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,bc_jm+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,bc_jm-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(pidy.ge.0+1)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(pidy.ge.0+1)     call MPI_WAIT(isend03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(isend04,istatus,ierr)
            if(pidy.ge.0+1)     call MPI_WAIT(irecv03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(.true.)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(.true.)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,bc_km+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(.true.)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,bc_km-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(.true.)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(.true.) call MPI_WAIT(isend05,istatus,ierr)
            if(.true.) call MPI_WAIT(isend06,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv05,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_XY BC
            call apply_nonperiodic_bc_i4(bc, a, ng)

      case('ADIABATIC_YZ')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(.true.)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(.true.)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(bc_im+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(.true.)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(bc_im-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(.true.)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(.true.) call MPI_WAIT(isend01,istatus,ierr)
            if(.true.) call MPI_WAIT(isend02,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv01,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(pidy.ge.0+1)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,bc_jm+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,bc_jm-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(pidy.ge.0+1)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(pidy.ge.0+1)     call MPI_WAIT(isend03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(isend04,istatus,ierr)
            if(pidy.ge.0+1)     call MPI_WAIT(irecv03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(pidz.ge.0+1)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(pidz.le.bc_npz-1-1)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,bc_km+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(pidz.le.bc_npz-1-1)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,bc_km-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(pidz.ge.0+1)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
            if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_YZ BC
            call apply_nonperiodic_bc_i4(bc, a, ng)

      case('ADIABATIC_ZX')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(pidx.ge.0+1)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(bc_im+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(bc_im-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(pidx.ge.0+1)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
            if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(.true.)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(.true.)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,bc_jm+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(.true.)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,bc_jm-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(.true.)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(.true.) call MPI_WAIT(isend03,istatus,ierr)
            if(.true.) call MPI_WAIT(isend04,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv03,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(pidz.ge.0+1)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(pidz.le.bc_npz-1-1)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,bc_km+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(pidz.le.bc_npz-1-1)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,bc_km-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(pidz.ge.0+1)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
            if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_ZX BC
            call apply_nonperiodic_bc_i4(bc, a, ng)

      case('ADIABATIC_X')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(pidx.ge.0+1)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(bc_im+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(bc_im-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(pidx.ge.0+1)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
            if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(.true.)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(.true.)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,bc_jm+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(.true.)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,bc_jm-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(.true.)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(.true.) call MPI_WAIT(isend03,istatus,ierr)
            if(.true.) call MPI_WAIT(isend04,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv03,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(.true.)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(.true.)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,bc_km+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(.true.)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,bc_km-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(.true.)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(.true.) call MPI_WAIT(isend05,istatus,ierr)
            if(.true.) call MPI_WAIT(isend06,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv05,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_X BC
            call apply_nonperiodic_bc_i4(bc, a, ng)

      case('ADIABATIC_Y')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(.true.)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(.true.)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(bc_im+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(.true.)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(bc_im-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(.true.)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(.true.) call MPI_WAIT(isend01,istatus,ierr)
            if(.true.) call MPI_WAIT(isend02,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv01,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(pidy.ge.0+1)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,bc_jm+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,bc_jm-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(pidy.ge.0+1)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(pidy.ge.0+1)     call MPI_WAIT(isend03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(isend04,istatus,ierr)
            if(pidy.ge.0+1)     call MPI_WAIT(irecv03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(.true.)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(.true.)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,bc_km+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(.true.)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,bc_km-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(.true.)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(.true.) call MPI_WAIT(isend05,istatus,ierr)
            if(.true.) call MPI_WAIT(isend06,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv05,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_Y BC
            call apply_nonperiodic_bc_i4(bc, a, ng)

      case('ADIABATIC_Z')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(.true.)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(.true.)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(bc_im+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(.true.)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(bc_im-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(.true.)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(.true.) call MPI_WAIT(isend01,istatus,ierr)
            if(.true.) call MPI_WAIT(isend02,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv01,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(.true.)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(.true.)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,bc_jm+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(.true.)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,bc_jm-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(.true.)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(.true.) call MPI_WAIT(isend03,istatus,ierr)
            if(.true.) call MPI_WAIT(isend04,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv03,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(pidz.ge.0+1)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(pidz.le.bc_npz-1-1)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,bc_km+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(pidz.le.bc_npz-1-1)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,bc_km-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(pidz.ge.0+1)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
            if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_Z BC
            call apply_nonperiodic_bc_i4(bc, a, ng)

      end select !select case(trim(bc))

      call mpi_barrier(MPI_COMM_WORLD, ierr)
      call mpi_type_free(xslab, ierr)
      call mpi_type_free(yslab, ierr)
      call mpi_type_free(zslab, ierr)

      end subroutine boundary_condition_i4

      ! ----------------------------------------------------------------

      subroutine apply_periodic_bc_i4(a,ng)
      implicit none
      integer(4),intent(inout) :: 
     &                      a(1-ng:bc_im+ng,1-ng:bc_jm+ng,1-ng:bc_km+ng)
      integer(4),intent(in) :: ng

      integer(4) :: i,j,k,lg

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, apply_periodic_bc_i4(),'
        write(*,*)'  boundary condition grid uninitialized. Aborting.'
        call abort
      endif


      if(bc_npx.gt.1.and.bc_npy.gt.1.and.bc_npz.gt.1)then
            return
      endif

      if(bc_npx.eq.1)then
            do k=1-ng,bc_km+ng
            do j=1-ng,bc_jm+ng
                  do lg=1,ng
                        a(-ng+lg,j,k) = a(bc_im-ng+lg,j,k)
                        a(bc_im+lg,j,k) = a(lg,j,k)
                  enddo
            enddo
            enddo
      endif
      if(bc_npy.eq.1)then
            do i=1-ng,bc_im+ng
            do k=1-ng,bc_km+ng
                  do lg=1,ng
                        a(i,-ng+lg,k) = a(i,bc_jm-ng+lg,k)
                        a(i,bc_jm+lg,k) = a(i,lg,k)
                  enddo
            enddo
            enddo
      endif
      if(bc_npz.eq.1)then
            do j=1-ng,bc_jm+ng
            do i=1-ng,bc_im+ng
                  do lg=1,ng
                        a(i,j,-ng+lg) = a(i,j,bc_km-ng+lg)
                        a(i,j,bc_km+lg) = a(i,j,lg)
                  enddo
            enddo
            enddo
      endif

      end subroutine apply_periodic_bc_i4

      ! ----------------------------------------------------------------
      
      subroutine apply_nonperiodic_bc_i4(bc, a, ng)
      implicit none
      !Subroutine arguments
      !character(len=*),intent(in) :: bctype
      type(bctype),intent(in) :: bc
      integer(4),intent(inout) :: 
     &                      a(1-ng:bc_im+ng,1-ng:bc_jm+ng,1-ng:bc_km+ng)
      integer(4),intent(in) :: ng

      !Local variables
      integer(4) :: i,j,k,lg
      integer(4) :: l1, l2, lh
      logical :: gconxn, gconyn, gconzn
      logical :: gconxp, gconyp, gconzp
      logical :: lconil, lconiu
      logical :: lconjl, lconju
      logical :: lconkl, lconku
      integer(4) :: pidx,pidy,pidz

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, apply_nonperiodic_bc_i4(),'
        write(*,*)'  boundary condition grid uninitialized. Aborting.'
        call abort
      endif

      if (bc%is_initialized .eqv. .false.) then
        write(*,*)'Error, apply_nonperiodic_bc_i4(),'
        write(*,*)'  bc handle uninitialized.'
        write(*,*)'  bc name:',trim(bc%handlename)
        write(*,*)'  Aborting...'
        call abort
      endif


      gconxn=.false.; gconyn=.false.; gconzn=.false.
      gconxp=.false.; gconyp=.false.; gconzp=.false.
      lconil=.false.; lconjl=.false.; lconkl=.false.
      lconiu=.false.; lconju=.false.; lconku=.false.

      if(trim(bc%xn) .ne. 'PERIODIC') gconxn = .true.
      if(trim(bc%xp) .ne. 'PERIODIC') gconxp = .true.
      if(trim(bc%yn) .ne. 'PERIODIC') gconyn = .true.
      if(trim(bc%yp) .ne. 'PERIODIC') gconyp = .true.
      if(trim(bc%zn) .ne. 'PERIODIC') gconzn = .true.
      if(trim(bc%zp) .ne. 'PERIODIC') gconzp = .true.

      call get_pidxyz(pidx, pidy, pidz, bc_pid)
      if(pidx.eq.0)     lconil = .true.
      if(pidx.eq.bc_npx-1) lconiu = .true.
      if(pidy.eq.0)     lconjl = .true.
      if(pidy.eq.bc_npy-1) lconju = .true.
      if(pidz.eq.0)     lconkl = .true.
      if(pidz.eq.bc_npz-1) lconku = .true.

      ! Combining global condition and local condition
      lconil=lconil.and.gconxn
      lconiu=lconiu.and.gconxp
      lconjl=lconjl.and.gconyn
      lconju=lconju.and.gconyp
      lconkl=lconkl.and.gconzn
      lconku=lconku.and.gconzp

      ! Terminate if this local process is nothing to do with BC
      if(lconil.or.lconiu.or.lconjl.or.lconju.or.lconkl.or.lconku
     &    .eqv. .false.) then
            return
      endif
 
      lh = max(bc_im,bc_jm,bc_km)  ! The highest # of grids

      do l2=1-ng,lh+ng
      do l1=1-ng,lh+ng

      ! x-layers: l1=j, l2=k
      if(lconil.or.lconiu .eqv. .true.)then
        if(l1.le.bc_jm+ng.and.l2.le.bc_km+ng)then  
          do lg=1,ng
            ! il
            if(lconil.eqv..true.)a(1-lg,l1,l2)=a(lg,l1,l2)
            ! iu
            if(lconiu.eqv..true.)a(bc_im+lg,l1,l2)=a(bc_im+1-lg,l1,l2)
          enddo
        endif
      endif !if(lconil.or.lconiu .eqv. .true.)then

      ! y-layers: l1=k, l2=i
      if(lconjl.or.lconju .eqv. .true.)then
        if(l1.le.bc_km+ng.and.l2.le.bc_im+ng)then  
          do lg=1,ng
            ! jl
            if(lconjl.eqv..true.)a(l2,1-lg,l1)=a(l2,lg,l1)
            ! ju
            if(lconju.eqv..true.)a(l2,bc_jm+lg,l1)=a(l2,bc_jm+1-lg,l1)
          enddo
        endif
      endif !if(lconjl.or.lconju .eqv. .true.)then

      ! z-layers: l1=i, l2=j
      if(lconkl.or.lconku .eqv. .true.)then
        if(l1.le.bc_im+ng.and.l2.le.bc_jm+ng)then  
          do lg=1,ng
            ! kl
            if(lconkl.eqv..true.)a(l1,l2,1-lg)=a(l1,l2,lg)
            ! ku
            if(lconku.eqv..true.)a(l1,l2,bc_km+lg)=a(l1,l2,bc_km+1-lg)
          enddo
        endif
      endif !if(lconkl.or.lconku .eqv. .true.)then
      

      ! Edges along x-axis: l2=i, l1=j, lg=k
      ! jlkl
      if(lconjl.and.lconkl .eqv. .true.) then
        if(l2.ge.1.and.l2.le.bc_im) then
          if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("jlkl",l2,1-ng+lg,1-ng+l1,ng)
                  a(l2,1-ng+lg,1-ng+l1)=a(l2,ng-l1,ng-lg)
            enddo
          endif
        endif
      endif
      ! jukl
      if(lconju.and.lconkl .eqv. .true.) then
        if(l2.ge.1.and.l2.le.bc_im) then
          if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("jukl",l2,bc_jm+1+lg,1-ng+l1,ng)
                  a(l2,bc_jm+1+lg,1-ng+l1)=a(l2,bc_jm+1-ng+l1,1+lg)
            enddo
          endif
        endif
      endif
      ! jlku
      if(lconjl.and.lconku .eqv. .true.) then
        if(l2.ge.1.and.l2.le.bc_im) then
          if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("jlku",l2,1-ng+lg,bc_km+1+l1,ng)
                  a(l2,1-ng+lg,bc_km+1+l1)=a(l2,1+l1,bc_km+1-ng+lg)
            enddo
          endif
        endif
      endif
      ! juku
      if(lconju.and.lconku .eqv. .true.) then
        if(l2.ge.1.and.l2.le.bc_im) then
          if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("juku",l2,bc_jm+1+lg,bc_km+1+l1,ng)
                  a(l2,bc_jm+1+lg,bc_km+1+l1)=a(l2,bc_jm-l1,bc_km-lg)
            enddo
          endif
        endif
      endif

      
      ! Edges along y-axis: l2=j, l1=k, lg=i
      ! klil
      if(lconkl.and.lconil .eqv. .true.) then
        if(l2.ge.1.and.l2.le.bc_jm) then
          if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("klil",1-ng+l1,l2,1-ng+lg,ng)
                  a(1-ng+l1,l2,1-ng+lg)=a(ng-lg,l2,ng-l1)
            enddo
          endif
        endif
      endif
      ! kuil
      if(lconku.and.lconil .eqv. .true.) then
        if(l2.ge.1.and.l2.le.bc_jm) then
          if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("kuil",bc_im+1+l1,l2,bc_km+1+lg,ng)
                  a(bc_im+1+l1,l2,bc_km+1+lg)=a(bc_im-lg,l2,bc_km-l1)
            enddo
          endif
        endif
      endif
      ! kliu
      if(lconkl.and.lconiu .eqv. .true.) then
        if(l2.ge.1.and.l2.le.bc_jm) then
          if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("kliu",1-ng+l1,l2,bc_km+1+lg,ng)
                  a(1-ng+l1,l2,bc_km+1+lg)=a(1+lg,l2,bc_km+1-ng+l1)
            enddo
          endif
        endif
      endif
      ! kuiu
      if(lconku.and.lconiu .eqv. .true.) then
        if(l2.ge.1.and.l2.le.bc_jm) then
          if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("kuiu",bc_im+1+l1,l2,bc_km+1+lg,ng)
                  a(bc_im+1+l1,l2,bc_km+1+lg)=a(bc_im-lg,l2,bc_km-l1)
            enddo
          endif
        endif
      endif
      
      ! Edges along z-axis: l2=k, l1=i, lg=j
      ! iljl
      if(lconil.and.lconjl .eqv. .true.) then
        if(l2.ge.1.and.l2.le.bc_km) then
          if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("iljl",1-ng+lg,1-ng+l1,l2,ng)
                  a(1-ng+lg,1-ng+l1,l2)=a(ng-l1,ng-lg,l2)
            enddo
          endif
        endif
      endif
      ! iujl
      if(lconiu.and.lconjl .eqv. .true.) then
        if(l2.ge.1.and.l2.le.bc_km) then
          if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("iujl",bc_im+1+lg,1-ng+l1,l2,ng)
                  a(bc_im+1+lg,1-ng+l1,l2)=a(bc_im+1-ng+l1,1+lg,l2)
            enddo
          endif
        endif
      endif
      ! ilju
      if(lconil.and.lconju .eqv. .true.) then
        if(l2.ge.1.and.l2.le.bc_km) then
          if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("ilju",1-ng+lg,bc_jm+1+l1,l2,ng)
                  a(1-ng+lg,bc_jm+1+l1,l2)=a(1+l1,bc_jm+1-ng+lg,l2)
            enddo
          endif
        endif
      endif
      ! iuju
      if(lconiu.and.lconju .eqv. .true.) then
        if(l2.ge.1.and.l2.le.bc_km) then
          if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("iuju",bc_im+1+lg,bc_jm+1+l1,l2,ng)
                  a(bc_im+1+lg,bc_jm+1+l1,l2)=a(bc_im-l1,bc_jm-lg,l2)
            enddo
          endif
        endif
      endif
      
      ! Corners
      ! lll ! l1=i, l2=j, lg=k
      if(lconil.and.lconjl.and.lconkl .eqv. .true.) then
        if(l2.ge.0.and.l2.le.ng-1) then
          if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
              a(1-ng+l1,1-ng+l2,1-ng+lg)=a(ng-lg,ng-l2,ng-l1)
            enddo
          endif
        endif
      endif

      ! ull ! l1=i, l2=j, lg=k
      if(lconiu.and.lconjl.and.lconkl .eqv. .true.) then
        if(l2.ge.0.and.l2.le.ng-1) then
          if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
              a(bc_im+1+l1,1-ng+l2,1-ng+lg)=a(bc_im-ng+1-lg,ng-l2,ng-l1)
            enddo
          endif
        endif
      endif

      ! lul ! l1=i, l2=j, lg=k
      if(lconil.and.lconju.and.lconkl .eqv. .true.) then
        if(l2.ge.0.and.l2.le.ng-1) then
          if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
              a(1-ng+l1,bc_jm+1+l2,1-ng+lg)=a(ng-lg,bc_jm-ng+1-l2,ng-l1)
            enddo
          endif
        endif
      endif

      ! uul ! l1=i, l2=j, lg=k
      if(lconiu.and.lconju.and.lconkl .eqv. .true.) then
        if(l2.ge.0.and.l2.le.ng-1) then
          if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
              a(bc_im+1+l1,bc_jm+1+l2,1-ng+lg)=
     &                              a(bc_im-ng+1-lg,bc_jm-ng+1-l2,ng-l1)
            enddo
          endif
        endif
      endif

      ! llu ! l1=i, l2=j, lg=k
      if(lconil.and.lconjl.and.lconku .eqv. .true.) then
        if(l2.ge.0.and.l2.le.ng-1) then
          if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
              a(1-ng+l1,1-ng+l2,bc_km+1+lg)=a(ng-lg,ng-l2,bc_km-ng+1-l1)
            enddo
          endif
        endif
      endif

      ! ulu ! l1=i, l2=j, lg=k
      if(lconiu.and.lconjl.and.lconku .eqv. .true.) then
        if(l2.ge.0.and.l2.le.ng-1) then
          if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
              a(bc_im+1+l1,1-ng+l2,bc_km+1+lg)=
     &                              a(bc_im-ng+1-lg,ng-l2,bc_km-ng+1-l1)
            enddo
          endif
        endif
      endif

      ! luu ! l1=i, l2=j, lg=k
      if(lconil.and.lconju.and.lconku .eqv. .true.) then
        if(l2.ge.0.and.l2.le.ng-1) then
          if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
              a(1-ng+l1,bc_jm+1+l2,bc_km+1+lg)=
     &                              a(ng-lg,bc_jm-ng+1-l2,bc_km-ng+1-l1)
            enddo
          endif
        endif
      endif

      ! uuu ! l1=i, l2=j, lg=k
      if(lconiu.and.lconju.and.lconku .eqv. .true.) then
        if(l2.ge.0.and.l2.le.ng-1) then
          if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
              a(bc_im+1+l1,bc_jm+1+l2,bc_km+1+lg)=
     &                      a(bc_im-ng+1-lg,bc_jm-ng+1-l2,bc_km-ng+1-l1)
            enddo
          endif
        endif
      endif

      enddo !do l1=1-ng,lh+ng
      enddo !do l2=1-ng,lh+ng

      end subroutine apply_nonperiodic_bc_i4

      ! ----------------------------------------------------------------

      subroutine boundary_condition_2d_r8(bc,a,ng)
      implicit none

      !Subroutine arguments
      !character(len=*),intent(in) :: bctype
      type(bctype),intent(in) :: bc
      real(8),intent(inout) :: a(1-ng:bc_im+ng,1,1-ng:bc_km+ng)
      integer(4),intent(in) :: ng

      !Local variables
      integer(4) :: i,j,k,n1,n2
      integer(4) :: il,iu,jl,ju,kl,ku

      ! MPI parallelization variables
      integer(4) :: pidx, pidy, pidz
      integer(4) :: pidto01, pidto02, pidto03, pidto04
      integer(4) :: pidto05, pidto06, pidto07, pidto08
      integer(4) :: pidto09, pidto10, pidto11, pidto12
      integer(4) :: pidfrom01, pidfrom02, pidfrom03, pidfrom04
      integer(4) :: pidfrom05, pidfrom06, pidfrom07, pidfrom08
      integer(4) :: pidfrom09, pidfrom10, pidfrom11, pidfrom12

      integer(4) :: isend01, isend02, isend03, isend04
      integer(4) :: isend05, isend06, isend07, isend08
      integer(4) :: isend09, isend10, isend11, isend12
      integer(4) :: irecv01, irecv02, irecv03, irecv04
      integer(4) :: irecv05, irecv06, irecv07, irecv08
      integer(4) :: irecv09, irecv10, irecv11, irecv12


      !MPI_HANDLES for MPI types for non-contiguous array data communication
      ! MPI type is a kind of MPI_HANDLE in C, but it is equivalent to integer in fortran
      integer(4) :: xslab, yslab, zslab

      integer(4) :: arrsize(3)
      integer(4) :: subsize(3)
      integer(4) :: starts(3)

      integer(4) :: istatus(MPI_STATUS_SIZE), ierr
      integer(4) :: basetag

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, boundary_condition_2d_r8(),'
        write(*,*)'  boundary condition grid uninitialized. Aborting.'
        call abort
      endif

      if (bc%is_initialized .eqv. .false.) then
        write(*,*)'Error, boundary_condition_2d_r8(),'
        write(*,*)'  bc handle uninitialized.'
        write(*,*)'  bc name:',trim(bc%handlename)
        write(*,*)'  Aborting...'
        call abort
      endif

      ! Resister the call count to determine the base MPI_TAG value
      bc_call_count = 1 + bc_call_count
      basetag = bc_call_count*bc_mpitag_jumpstep
      bc_max_mpi_tag = max(bc_max_mpi_tag,basetag+bc_mpitag_jumpstep)

      arrsize(1) = 2*ng + bc_im
      arrsize(2) = 1
      arrsize(3) = 2*ng + bc_km

      il = 1-ng
      iu = bc_im+ng
      jl = 1
      ju = 1
      kl = 1-ng
      ku = bc_km+ng

      ! MPI derived data type definition for non-cnotiguous transfer
      ! xslab
      starts(1) = 0
      starts(2) = 0
      starts(3) = 0

      subsize(1) = ng
      subsize(2) = arrsize(2)
      subsize(3) = arrsize(3)
      call mpi_type_create_subarray(3, arrsize, subsize, starts,
     &  MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, xslab, ierr)
      call mpi_type_commit(xslab, ierr)

      ! zslab
      subsize(1) = arrsize(1)
      subsize(2) = arrsize(2)
      subsize(3) = ng
      call mpi_type_create_subarray(3, arrsize, subsize, starts,
     &  MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, zslab, ierr)
      call mpi_type_commit(zslab, ierr)

      select case(trim(bc%legacy))
      case('PERIODIC')
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(bc_npx.gt.1)then
            pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
            call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)

            pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
            call MPI_IRECV(a(bc_im+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)

            pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
            call MPI_ISEND(a(bc_im-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      

            pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
            call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            
            call MPI_WAIT(isend01,istatus,ierr)
            call MPI_WAIT(isend02,istatus,ierr)
            call MPI_WAIT(irecv01,istatus,ierr)
            call MPI_WAIT(irecv02,istatus,ierr)
            endif
            ! Done only for x direction - not complete yet

            ! k planes
            if(bc_npz.gt.1)then
            pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
            call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)

            pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
            call MPI_IRECV(a(il,jl,bc_km+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            
            pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
            call MPI_ISEND(a(il,jl,bc_km-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)

            pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
            call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)

            call MPI_WAIT(isend05,istatus,ierr)
            call MPI_WAIT(isend06,istatus,ierr)
            call MPI_WAIT(irecv05,istatus,ierr)
            call MPI_WAIT(irecv06,istatus,ierr)
            endif

            call apply_periodic_bc_2d_r8(a, ng)
            ! Done for all directions - complete

      case('ADIABATIC')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(pidx.ge.0+1)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(bc_im+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(bc_im-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(pidx.ge.0+1)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
            if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! k planes
            if(pidz.ge.0+1)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(pidz.le.bc_npz-1-1)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,bc_km+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(pidz.le.bc_npz-1-1)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,bc_km-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(pidz.ge.0+1)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
            if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC BC
            call apply_nonperiodic_bc_2d_r8(bc, a, ng)

      case('ADIABATIC_X')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(pidx.ge.0+1)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(bc_im+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(bc_im-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(pidx.ge.0+1)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
            if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! k planes
            if(.true.)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(.true.)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,bc_km+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(.true.)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,bc_km-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(.true.)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(.true.) call MPI_WAIT(isend05,istatus,ierr)
            if(.true.) call MPI_WAIT(isend06,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv05,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_X BC
            call apply_nonperiodic_bc_2d_r8(bc, a, ng)

      case('ADIABATIC_Z')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(.true.)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(.true.)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(bc_im+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(.true.)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(bc_im-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(.true.)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(.true.) call MPI_WAIT(isend01,istatus,ierr)
            if(.true.) call MPI_WAIT(isend02,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv01,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! k planes
            if(pidz.ge.0+1)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(pidz.le.bc_npz-1-1)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,bc_km+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(pidz.le.bc_npz-1-1)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,bc_km-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(pidz.ge.0+1)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
            if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_Z BC
            call apply_nonperiodic_bc_2d_r8(bc, a, ng)

      end select !select case(trim(bc))

      call mpi_barrier(MPI_COMM_WORLD, ierr)
      call mpi_type_free(xslab, ierr)
      call mpi_type_free(zslab, ierr)

      end subroutine boundary_condition_2d_r8

      ! ----------------------------------------------------------------

      subroutine apply_periodic_bc_2d_r8(a,ng)
      implicit none
      real(8),intent(inout) :: a(1-ng:bc_im+ng,1,1-ng:bc_km+ng)
      integer(4),intent(in) :: ng

      integer(4) :: i,j,k,lg

      if(bc_npx.gt.1.and.bc_npz.gt.1)then
            return
      endif

      if(bc_npx.eq.1)then
            do k=1-ng,bc_km+ng
                  do lg=1,ng
                        a(-ng+lg,1,k) = a(bc_im-ng+lg,1,k)
                        a(bc_im+lg,1,k) = a(lg,1,k)
                  enddo
            enddo
      endif

      if(bc_npz.eq.1)then
            do i=1-ng,bc_im+ng
                  do lg=1,ng
                        a(i,1,-ng+lg) = a(i,1,bc_km-ng+lg)
                        a(i,1,bc_km+lg) = a(i,1,lg)
                  enddo
            enddo
      endif

      end subroutine apply_periodic_bc_2d_r8

      ! ----------------------------------------------------------------

      subroutine apply_nonperiodic_bc_2d_r8(bc, a, ng)
      implicit none
      !Subroutine arguments
      !character(len=*),intent(in) :: bctype
      type(bctype),intent(in) :: bc
      real(8),intent(inout) :: a(1-ng:bc_im+ng,1,1-ng:bc_km+ng)
      integer(4),intent(in) :: ng

      !Local variables
      integer(4) :: i,j,k,lg
      integer(4) :: l1, l2, lh
      logical :: gconxn, gconzn
      logical :: gconxp, gconzp
      logical :: lconil, lconiu
      logical :: lconkl, lconku
      integer(4) :: pidx,pidy,pidz


      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, apply_nonperiodic_bc_2d_r8(),'
        write(*,*)'  boundary condition grid uninitialized. Aborting.'
        call abort
      endif

      if (bc%is_initialized .eqv. .false.) then
        write(*,*)'Error, apply_nonperiodic_bc_2d_r8(),'
        write(*,*)'  bc handle uninitialized.'
        write(*,*)'  bc name:',trim(bc%handlename)
        write(*,*)'  Aborting...'
        call abort
      endif


      gconxn=.false.; gconzn=.false.
      gconxp=.false.; gconzp=.false.
      lconil=.false.; lconkl=.false.
      lconiu=.false.; lconku=.false.

      if(trim(bc%xn) .ne. 'PERIODIC') gconxn = .true.
      if(trim(bc%xp) .ne. 'PERIODIC') gconxp = .true.
      if(trim(bc%zn) .ne. 'PERIODIC') gconzn = .true.
      if(trim(bc%zp) .ne. 'PERIODIC') gconzp = .true.

      call get_pidxyz(pidx, pidy, pidz, bc_pid)
      if(pidx.eq.0)     lconil = .true.
      if(pidx.eq.bc_npx-1) lconiu = .true.
      if(pidz.eq.0)     lconkl = .true.
      if(pidz.eq.bc_npz-1) lconku = .true.

      ! Combining global condition and local condition
      lconil=lconil.and.gconxn
      lconiu=lconiu.and.gconxp
      lconkl=lconkl.and.gconzn
      lconku=lconku.and.gconzp

      ! Terminate if this local process is nothing to do with BC
      if(lconil.or.lconiu.or.lconkl.or.lconku
     &    .eqv. .false.) then
            return
      endif
 
      lh = max(bc_im,1,bc_km)  ! The highest # of grids

      do l1=1-ng,lh+ng

      ! x-layers: l1=j=1, l2=k
      if(lconil.or.lconiu .eqv. .true.)then
        if(l1.le.bc_km+ng)then  
          do lg=1,ng
            ! il
            if(lconil.eqv..true.)a(1-lg,1,l1)=a(lg,1,l1)
            ! iu
            if(lconiu.eqv..true.)a(bc_im+lg,1,l1)=a(bc_im+1-lg,1,l1)
          enddo
        endif
      endif !if(lconil.or.lconiu .eqv. .true.)then

      ! z-layers: l1=i, l2=j=1
      if(lconkl.or.lconku .eqv. .true.)then
        if(l1.le.bc_im+ng)then  
          do lg=1,ng
            ! kl
            if(lconkl.eqv..true.)a(l1,1,1-lg)=a(l1,1,lg)
            ! ku
            if(lconku.eqv..true.)a(l1,1,bc_km+lg)=a(l1,1,bc_km+1-lg)
          enddo
        endif
      endif !if(lconkl.or.lconku .eqv. .true.)then
      

      ! Edges along y-axis: l2=j=1, l1=k, lg=i
      ! klil
      if(lconkl.and.lconil .eqv. .true.) then
        if(l1.ge.0.and.l1.le.ng-1) then
          do lg=0,ng-1
            call bound_check("klil",1-ng+l1,1,1-ng+lg,ng)
            a(1-ng+l1,1,1-ng+lg)=a(ng-lg,1,ng-l1)
          enddo
        endif
      endif
      ! kuil
      if(lconku.and.lconil .eqv. .true.) then
        if(l1.ge.0.and.l1.le.ng-1) then
          do lg=0,ng-1
            call bound_check("kuil",bc_im+1+l1,1,bc_km+1+lg,ng)
            a(bc_im+1+l1,1,bc_km+1+lg)=a(bc_im-lg,1,bc_km-l1)
          enddo
        endif
      endif
      ! kliu
      if(lconkl.and.lconiu .eqv. .true.) then
        if(l1.ge.0.and.l1.le.ng-1) then
          do lg=0,ng-1
            call bound_check("kliu",1-ng+l1,1,bc_km+1+lg,ng)
            a(1-ng+l1,1,bc_km+1+lg)=a(1+lg,1,bc_km+1-ng+l1)
          enddo
        endif
      endif
      ! kuiu
      if(lconku.and.lconiu .eqv. .true.) then
        if(l1.ge.0.and.l1.le.ng-1) then
          do lg=0,ng-1
            call bound_check("kuiu",bc_im+1+l1,1,bc_km+1+lg,ng)
            a(bc_im+1+l1,1,bc_km+1+lg)=a(bc_im-lg,1,bc_km-l1)
          enddo
        endif
      endif
      
      
      enddo !do l1=1-ng,lh+ng

      end subroutine apply_nonperiodic_bc_2d_r8

      ! ----------------------------------------------------------------

      subroutine boundary_condition_r8(bc,a,ng)
      implicit none

      !Subroutine arguments
      !character(len=*),intent(in) :: bctype
      type(bctype),intent(in) :: bc
      real(8),intent(inout) :: 
     &                      a(1-ng:bc_im+ng,1-ng:bc_jm+ng,1-ng:bc_km+ng)
      integer(4),intent(in) :: ng

      !Local variables
      integer(4) :: i,j,k,n1,n2
      integer(4) :: il,iu,jl,ju,kl,ku

      ! MPI parallelization variables
      integer(4) :: pidx, pidy, pidz
      integer(4) :: pidto01, pidto02, pidto03, pidto04
      integer(4) :: pidto05, pidto06, pidto07, pidto08
      integer(4) :: pidto09, pidto10, pidto11, pidto12
      integer(4) :: pidfrom01, pidfrom02, pidfrom03, pidfrom04
      integer(4) :: pidfrom05, pidfrom06, pidfrom07, pidfrom08
      integer(4) :: pidfrom09, pidfrom10, pidfrom11, pidfrom12

      integer(4) :: isend01, isend02, isend03, isend04
      integer(4) :: isend05, isend06, isend07, isend08
      integer(4) :: isend09, isend10, isend11, isend12
      integer(4) :: irecv01, irecv02, irecv03, irecv04
      integer(4) :: irecv05, irecv06, irecv07, irecv08
      integer(4) :: irecv09, irecv10, irecv11, irecv12


      !MPI_HANDLES for MPI types for non-contiguous array data communication
      ! MPI type is a kind of MPI_HANDLE in C, but it is equivalent to integer in fortran
      integer(4) :: xslab, yslab, zslab

      integer(4) :: arrsize(3)
      integer(4) :: subsize(3)
      integer(4) :: starts(3)

      integer(4) :: istatus(MPI_STATUS_SIZE), ierr
      integer(4) :: basetag

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, boundary_condition_r8(),'
        write(*,*)'  boundary condition grid uninitialized. Aborting.'
        call abort
      endif

      if (bc%is_initialized .eqv. .false.) then
        write(*,*)'Error, boundary_condition_r8(),'
        write(*,*)'  bc handle uninitialized.'
        write(*,*)'  bc name:',trim(bc%handlename)
        write(*,*)'  Aborting...'
        call abort
      endif


      ! Resister the call count to determine the base MPI_TAG value
      bc_call_count = 1 + bc_call_count
      basetag = bc_call_count*bc_mpitag_jumpstep
      bc_max_mpi_tag = max(bc_max_mpi_tag,basetag+bc_mpitag_jumpstep)


      arrsize(1) = 2*ng + bc_im
      arrsize(2) = 2*ng + bc_jm
      arrsize(3) = 2*ng + bc_km

      il = 1-ng
      iu = bc_im+ng
      jl = 1-ng
      ju = bc_jm+ng
      kl = 1-ng
      ku = bc_km+ng

      ! MPI derived data type definition for non-cnotiguous transfer
      ! xslab
      starts(1) = 0
      starts(2) = 0
      starts(3) = 0

      subsize(1) = ng
      subsize(2) = arrsize(2)
      subsize(3) = arrsize(3)
      call mpi_type_create_subarray(3, arrsize, subsize, starts,
     &  MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, xslab, ierr)
      call mpi_type_commit(xslab, ierr)

      ! yslab
      subsize(1) = arrsize(1)
      subsize(2) = ng
      subsize(3) = arrsize(3)
      call mpi_type_create_subarray(3, arrsize, subsize, starts,
     &  MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, yslab, ierr)
      call mpi_type_commit(yslab, ierr)

      ! zslab
      subsize(1) = arrsize(1)
      subsize(2) = arrsize(2)
      subsize(3) = ng
      call mpi_type_create_subarray(3, arrsize, subsize, starts,
     &  MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, zslab, ierr)
      call mpi_type_commit(zslab, ierr)



      select case(trim(bc%legacy))
      case('PERIODIC')
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            
            ! Layer send/recv
            ! i planes
            if(bc_npx.gt.1)then
            pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
            call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)

            pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
            call MPI_IRECV(a(bc_im+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)

            pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
            call MPI_ISEND(a(bc_im-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      

            pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
            call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            
            call MPI_WAIT(isend01,istatus,ierr)
            call MPI_WAIT(isend02,istatus,ierr)
            call MPI_WAIT(irecv01,istatus,ierr)
            call MPI_WAIT(irecv02,istatus,ierr)
            endif
            ! Done only for x direction - not complete yet

            ! j planes
            if(bc_npy.gt.1)then
            pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
            call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)

            pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
            call MPI_IRECV(a(il,bc_jm+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)

            pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
            call MPI_ISEND(a(il,bc_jm-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)

            pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
            call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)

            call MPI_WAIT(isend03,istatus,ierr)
            call MPI_WAIT(isend04,istatus,ierr)
            call MPI_WAIT(irecv03,istatus,ierr)
            call MPI_WAIT(irecv04,istatus,ierr)
            endif
            ! Done for x and y direction - not completed yet

            ! k planes
            if(bc_npz.gt.1)then
            pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
            call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)

            pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
            call MPI_IRECV(a(il,jl,bc_km+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            
            pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
            call MPI_ISEND(a(il,jl,bc_km-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)

            pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
            call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)

            call MPI_WAIT(isend05,istatus,ierr)
            call MPI_WAIT(isend06,istatus,ierr)
            call MPI_WAIT(irecv05,istatus,ierr)
            call MPI_WAIT(irecv06,istatus,ierr)
            endif

            call apply_periodic_bc_r8(a, ng)
            ! Done for all directions - complete

      case('ADIABATIC')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(pidx.ge.0+1)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(bc_im+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(bc_im-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(pidx.ge.0+1)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
            if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(pidy.ge.0+1)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,bc_jm+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,bc_jm-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(pidy.ge.0+1)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(pidy.ge.0+1)     call MPI_WAIT(isend03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(isend04,istatus,ierr)
            if(pidy.ge.0+1)     call MPI_WAIT(irecv03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(pidz.ge.0+1)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(pidz.le.bc_npz-1-1)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,bc_km+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(pidz.le.bc_npz-1-1)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,bc_km-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(pidz.ge.0+1)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
            if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC BC
            call apply_nonperiodic_bc_r8(bc, a, ng)

      case('ADIABATIC_XY')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(pidx.ge.0+1)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(bc_im+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(bc_im-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(pidx.ge.0+1)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
            if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(pidy.ge.0+1)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,bc_jm+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,bc_jm-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(pidy.ge.0+1)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(pidy.ge.0+1)     call MPI_WAIT(isend03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(isend04,istatus,ierr)
            if(pidy.ge.0+1)     call MPI_WAIT(irecv03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(.true.)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(.true.)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,bc_km+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(.true.)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,bc_km-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(.true.)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(.true.) call MPI_WAIT(isend05,istatus,ierr)
            if(.true.) call MPI_WAIT(isend06,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv05,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_XY BC
            call apply_nonperiodic_bc_r8(bc, a, ng)

      case('ADIABATIC_YZ')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(.true.)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(.true.)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(bc_im+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(.true.)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(bc_im-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(.true.)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(.true.) call MPI_WAIT(isend01,istatus,ierr)
            if(.true.) call MPI_WAIT(isend02,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv01,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(pidy.ge.0+1)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,bc_jm+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,bc_jm-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(pidy.ge.0+1)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(pidy.ge.0+1)     call MPI_WAIT(isend03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(isend04,istatus,ierr)
            if(pidy.ge.0+1)     call MPI_WAIT(irecv03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(pidz.ge.0+1)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(pidz.le.bc_npz-1-1)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,bc_km+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(pidz.le.bc_npz-1-1)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,bc_km-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(pidz.ge.0+1)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
            if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_YZ BC
            call apply_nonperiodic_bc_r8(bc, a, ng)

      case('ADIABATIC_ZX')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(pidx.ge.0+1)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(bc_im+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(bc_im-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(pidx.ge.0+1)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
            if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(.true.)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(.true.)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,bc_jm+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(.true.)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,bc_jm-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(.true.)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(.true.) call MPI_WAIT(isend03,istatus,ierr)
            if(.true.) call MPI_WAIT(isend04,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv03,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(pidz.ge.0+1)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(pidz.le.bc_npz-1-1)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,bc_km+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(pidz.le.bc_npz-1-1)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,bc_km-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(pidz.ge.0+1)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
            if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_ZX BC
            call apply_nonperiodic_bc_r8(bc, a, ng)

      case('ADIABATIC_X')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(pidx.ge.0+1)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(bc_im+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(pidx.le.bc_npx-1-1)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(bc_im-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(pidx.ge.0+1)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(pidx.ge.0+1)     call MPI_WAIT(isend01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(isend02,istatus,ierr)
            if(pidx.ge.0+1)     call MPI_WAIT(irecv01,istatus,ierr)
            if(pidx.le.bc_npx-1-1) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(.true.)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(.true.)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,bc_jm+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(.true.)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,bc_jm-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(.true.)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(.true.) call MPI_WAIT(isend03,istatus,ierr)
            if(.true.) call MPI_WAIT(isend04,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv03,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(.true.)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(.true.)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,bc_km+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(.true.)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,bc_km-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(.true.)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(.true.) call MPI_WAIT(isend05,istatus,ierr)
            if(.true.) call MPI_WAIT(isend06,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv05,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_X BC
            call apply_nonperiodic_bc_r8(bc, a, ng)

      case('ADIABATIC_Y')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(.true.)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(.true.)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(bc_im+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(.true.)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(bc_im-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(.true.)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(.true.) call MPI_WAIT(isend01,istatus,ierr)
            if(.true.) call MPI_WAIT(isend02,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv01,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(pidy.ge.0+1)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,bc_jm+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(pidy.le.bc_npy-1-1)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,bc_jm-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(pidy.ge.0+1)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(pidy.ge.0+1)     call MPI_WAIT(isend03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(isend04,istatus,ierr)
            if(pidy.ge.0+1)     call MPI_WAIT(irecv03,istatus,ierr)
            if(pidy.le.bc_npy-1-1) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(.true.)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(.true.)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,bc_km+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(.true.)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,bc_km-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(.true.)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(.true.) call MPI_WAIT(isend05,istatus,ierr)
            if(.true.) call MPI_WAIT(isend06,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv05,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_Y BC
            call apply_nonperiodic_bc_r8(bc, a, ng)

      case('ADIABATIC_Z')
            ! Necessary data transfer among processes
            !Process coordinate
            call get_pidxyz(pidx, pidy, pidz, bc_pid)

            ! Layer send/recv
            ! i planes
            if(.true.)then
                  pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
                  call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)
            endif

            if(.true.)then
                  pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
                  call MPI_IRECV(a(bc_im+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)
            endif

            if(.true.)then
                  pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
                  call MPI_ISEND(a(bc_im-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      
            endif

            if(.true.)then
                  pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
            endif
            
            if(.true.) call MPI_WAIT(isend01,istatus,ierr)
            if(.true.) call MPI_WAIT(isend02,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv01,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv02,istatus,ierr)
            ! Done only for x direction - not complete yet

            ! j planes
            if(.true.)then
                  pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
                  call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)
            endif

            if(.true.)then
                  pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
                  call MPI_IRECV(a(il,bc_jm+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)
            endif

            if(.true.)then
                  pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
                  call MPI_ISEND(a(il,bc_jm-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)
            endif

            if(.true.)then
                  pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)
            endif

            if(.true.) call MPI_WAIT(isend03,istatus,ierr)
            if(.true.) call MPI_WAIT(isend04,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv03,istatus,ierr)
            if(.true.) call MPI_WAIT(irecv04,istatus,ierr)
            ! Done for x and y direction - complete

            ! k planes
            if(pidz.ge.0+1)then
                  pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
                  call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)
            endif

            if(pidz.le.bc_npz-1-1)then
                  pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
                  call MPI_IRECV(a(il,jl,bc_km+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
            endif
            
            if(pidz.le.bc_npz-1-1)then
                  pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
                  call MPI_ISEND(a(il,jl,bc_km-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)
            endif

            if(pidz.ge.0+1)then
                  pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
                  call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)
            endif

            if(pidz.ge.0+1)     call MPI_WAIT(isend05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(isend06,istatus,ierr)
            if(pidz.ge.0+1)     call MPI_WAIT(irecv05,istatus,ierr)
            if(pidz.le.bc_npz-1-1) call MPI_WAIT(irecv06,istatus,ierr)
            ! Done for x and y direction - complete

            ! Apply ADIABATIC_Z BC
            call apply_nonperiodic_bc_r8(bc, a, ng)

      end select !select case(trim(bc))

      call mpi_barrier(MPI_COMM_WORLD, ierr)
      call mpi_type_free(xslab, ierr)
      call mpi_type_free(yslab, ierr)
      call mpi_type_free(zslab, ierr)

      end subroutine boundary_condition_r8

      ! ----------------------------------------------------------------

      subroutine apply_periodic_bc_r8(a,ng)
      implicit none
      real(8),intent(inout) :: 
     &                      a(1-ng:bc_im+ng,1-ng:bc_jm+ng,1-ng:bc_km+ng)
      integer(4),intent(in) :: ng

      integer(4) :: i,j,k,lg

      if(bc_npx.gt.1.and.bc_npy.gt.1.and.bc_npz.gt.1)then
            return
      endif

      if(bc_npx.eq.1)then
            do k=1-ng,bc_km+ng
            do j=1-ng,bc_jm+ng
                  do lg=1,ng
                        a(-ng+lg,j,k) = a(bc_im-ng+lg,j,k)
                        a(bc_im+lg,j,k) = a(lg,j,k)
                  enddo
            enddo
            enddo
      endif
      if(bc_npy.eq.1)then
            do i=1-ng,bc_im+ng
            do k=1-ng,bc_km+ng
                  do lg=1,ng
                        a(i,-ng+lg,k) = a(i,bc_jm-ng+lg,k)
                        a(i,bc_jm+lg,k) = a(i,lg,k)
                  enddo
            enddo
            enddo
      endif
      if(bc_npz.eq.1)then
            do j=1-ng,bc_jm+ng
            do i=1-ng,bc_im+ng
                  do lg=1,ng
                        a(i,j,-ng+lg) = a(i,j,bc_km-ng+lg)
                        a(i,j,bc_km+lg) = a(i,j,lg)
                  enddo
            enddo
            enddo
      endif

      end subroutine apply_periodic_bc_r8

      ! ----------------------------------------------------------------
      
      subroutine apply_nonperiodic_bc_r8(bc, a, ng)
      implicit none
      !Subroutine arguments
      !character(len=*),intent(in) :: bctype
      type(bctype),intent(in) :: bc
      real(8),intent(inout) :: 
     &                      a(1-ng:bc_im+ng,1-ng:bc_jm+ng,1-ng:bc_km+ng)
      integer(4),intent(in) :: ng

      !Local variables
      integer(4) :: i,j,k,lg
      integer(4) :: l1, l2, lh
      logical :: gconxn, gconyn, gconzn
      logical :: gconxp, gconyp, gconzp
      logical :: lconil, lconiu
      logical :: lconjl, lconju
      logical :: lconkl, lconku
      integer(4) :: pidx,pidy,pidz

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, apply_nonperiodic_bc_r8(),'
        write(*,*)'  boundary condition grid uninitialized. Aborting.'
        call abort
      endif

      if (bc%is_initialized .eqv. .false.) then
        write(*,*)'Error, apply_nonperiodic_bc_r8(),'
        write(*,*)'  bc handle uninitialized.'
        write(*,*)'  bc name:',trim(bc%handlename)
        write(*,*)'  Aborting...'
        call abort
      endif


      gconxn=.false.; gconyn=.false.; gconzn=.false.
      gconxp=.false.; gconyp=.false.; gconzp=.false.
      lconil=.false.; lconjl=.false.; lconkl=.false.
      lconiu=.false.; lconju=.false.; lconku=.false.

      if(trim(bc%xn) .ne. 'PERIODIC') gconxn = .true.
      if(trim(bc%xp) .ne. 'PERIODIC') gconxp = .true.
      if(trim(bc%yn) .ne. 'PERIODIC') gconyn = .true.
      if(trim(bc%yp) .ne. 'PERIODIC') gconyp = .true.
      if(trim(bc%zn) .ne. 'PERIODIC') gconzn = .true.
      if(trim(bc%zp) .ne. 'PERIODIC') gconzp = .true.

      call get_pidxyz(pidx, pidy, pidz, bc_pid)
      if(pidx.eq.0)     lconil = .true.
      if(pidx.eq.bc_npx-1) lconiu = .true.
      if(pidy.eq.0)     lconjl = .true.
      if(pidy.eq.bc_npy-1) lconju = .true.
      if(pidz.eq.0)     lconkl = .true.
      if(pidz.eq.bc_npz-1) lconku = .true.

      ! Combining global condition and local condition
      lconil=lconil.and.gconxn
      lconiu=lconiu.and.gconxp
      lconjl=lconjl.and.gconyn
      lconju=lconju.and.gconyp
      lconkl=lconkl.and.gconzn
      lconku=lconku.and.gconzp

      ! Terminate if this local process is nothing to do with BC
      if(lconil.or.lconiu.or.lconjl.or.lconju.or.lconkl.or.lconku
     &    .eqv. .false.) then
            return
      endif
 
      lh = max(bc_im,bc_jm,bc_km)  ! The highest # of grids

      do l2=1-ng,lh+ng
      do l1=1-ng,lh+ng

        ! x-layers: l1=j, l2=k
        if(lconil.or.lconiu .eqv. .true.)then
          if(l1.le.bc_jm+ng.and.l2.le.bc_km+ng)then  
            do lg=1,ng
              ! il
              if(lconil.eqv..true.)a(1-lg,l1,l2)=a(lg,l1,l2)
              ! iu
              if(lconiu.eqv..true.)a(bc_im+lg,l1,l2)=a(bc_im+1-lg,l1,l2)
            enddo
          endif
        endif !if(lconil.or.lconiu .eqv. .true.)then

        ! y-layers: l1=k, l2=i
        if(lconjl.or.lconju .eqv. .true.)then
          if(l1.le.bc_km+ng.and.l2.le.bc_im+ng)then  
            do lg=1,ng
              ! jl
              if(lconjl.eqv..true.)a(l2,1-lg,l1)=a(l2,lg,l1)
              ! ju
              if(lconju.eqv..true.)a(l2,bc_jm+lg,l1)=a(l2,bc_jm+1-lg,l1)
            enddo
          endif
        endif !if(lconjl.or.lconju .eqv. .true.)then

        ! z-layers: l1=i, l2=j
        if(lconkl.or.lconku .eqv. .true.)then
          if(l1.le.bc_im+ng.and.l2.le.bc_jm+ng)then  
            do lg=1,ng
              ! kl
              if(lconkl.eqv..true.)a(l1,l2,1-lg)=a(l1,l2,lg)
              ! ku
              if(lconku.eqv..true.)a(l1,l2,bc_jm+lg)=a(l1,l2,bc_jm+1-lg)
            enddo
          endif
        endif !if(lconkl.or.lconku .eqv. .true.)then
      

        ! Edges along x-axis: l2=i, l1=j, lg=k
        ! jlkl
        if(lconjl.and.lconkl .eqv. .true.) then
          if(l2.ge.1.and.l2.le.bc_im) then
            if(l1.ge.0.and.l1.le.ng-1) then
              do lg=0,ng-1
                    call bound_check("jlkl",l2,1-ng+lg,1-ng+l1,ng)
                    a(l2,1-ng+lg,1-ng+l1)=a(l2,ng-l1,ng-lg)
              enddo
            endif
          endif
        endif
        ! jukl
        if(lconju.and.lconkl .eqv. .true.) then
          if(l2.ge.1.and.l2.le.bc_im) then
            if(l1.ge.0.and.l1.le.ng-1) then
              do lg=0,ng-1
                    call bound_check("jukl",l2,bc_jm+1+lg,1-ng+l1,ng)
                    a(l2,bc_jm+1+lg,1-ng+l1)=a(l2,bc_jm+1-ng+l1,1+lg)
              enddo
            endif
          endif
        endif
        ! jlku
        if(lconjl.and.lconku .eqv. .true.) then
          if(l2.ge.1.and.l2.le.bc_im) then
            if(l1.ge.0.and.l1.le.ng-1) then
              do lg=0,ng-1
                    call bound_check("jlku",l2,1-ng+lg,bc_km+1+l1,ng)
                    a(l2,1-ng+lg,bc_km+1+l1)=a(l2,1+l1,bc_km+1-ng+lg)
              enddo
            endif
          endif
        endif
        ! juku
        if(lconju.and.lconku .eqv. .true.) then
          if(l2.ge.1.and.l2.le.bc_im) then
            if(l1.ge.0.and.l1.le.ng-1) then
              do lg=0,ng-1
                    call bound_check("juku",l2,bc_jm+1+lg,bc_km+1+l1,ng)
                    a(l2,bc_jm+1+lg,bc_km+1+l1)=a(l2,bc_jm-l1,bc_km-lg)
              enddo
            endif
          endif
        endif

      
        ! Edges along y-axis: l2=j, l1=k, lg=i
        ! klil
        if(lconkl.and.lconil .eqv. .true.) then
          if(l2.ge.1.and.l2.le.bc_jm) then
            if(l1.ge.0.and.l1.le.ng-1) then
              do lg=0,ng-1
                    call bound_check("klil",1-ng+l1,l2,1-ng+lg,ng)
                    a(1-ng+l1,l2,1-ng+lg)=a(ng-lg,l2,ng-l1)
              enddo
            endif
          endif
        endif
        ! kuil
        if(lconku.and.lconil .eqv. .true.) then
          if(l2.ge.1.and.l2.le.bc_jm) then
            if(l1.ge.0.and.l1.le.ng-1) then
              do lg=0,ng-1
                    call bound_check("kuil",bc_im+1+l1,l2,bc_km+1+lg,ng)
                    a(bc_im+1+l1,l2,bc_km+1+lg)=a(bc_im-lg,l2,bc_km-l1)
              enddo
            endif
          endif
        endif
        ! kliu
        if(lconkl.and.lconiu .eqv. .true.) then
          if(l2.ge.1.and.l2.le.bc_jm) then
            if(l1.ge.0.and.l1.le.ng-1) then
              do lg=0,ng-1
                    call bound_check("kliu",1-ng+l1,l2,bc_km+1+lg,ng)
                    a(1-ng+l1,l2,bc_km+1+lg)=a(1+lg,l2,bc_km+1-ng+l1)
              enddo
            endif
          endif
        endif
        ! kuiu
        if(lconku.and.lconiu .eqv. .true.) then
          if(l2.ge.1.and.l2.le.bc_jm) then
            if(l1.ge.0.and.l1.le.ng-1) then
              do lg=0,ng-1
                    call bound_check("kuiu",bc_im+1+l1,l2,bc_km+1+lg,ng)
                    a(bc_im+1+l1,l2,bc_km+1+lg)=a(bc_im-lg,l2,bc_km-l1)
              enddo
            endif
          endif
        endif
      
        ! Edges along z-axis: l2=k, l1=i, lg=j
        ! iljl
        if(lconil.and.lconjl .eqv. .true.) then
          if(l2.ge.1.and.l2.le.bc_km) then
            if(l1.ge.0.and.l1.le.ng-1) then
              do lg=0,ng-1
                    call bound_check("iljl",1-ng+lg,1-ng+l1,l2,ng)
                    a(1-ng+lg,1-ng+l1,l2)=a(ng-l1,ng-lg,l2)
              enddo
            endif
          endif
        endif
        ! iujl
        if(lconiu.and.lconjl .eqv. .true.) then
          if(l2.ge.1.and.l2.le.bc_km) then
            if(l1.ge.0.and.l1.le.ng-1) then
              do lg=0,ng-1
                    call bound_check("iujl",bc_im+1+lg,1-ng+l1,l2,ng)
                    a(bc_im+1+lg,1-ng+l1,l2)=a(bc_im+1-ng+l1,1+lg,l2)
              enddo
            endif
          endif
        endif
        ! ilju
        if(lconil.and.lconju .eqv. .true.) then
          if(l2.ge.1.and.l2.le.bc_km) then
            if(l1.ge.0.and.l1.le.ng-1) then
              do lg=0,ng-1
                    call bound_check("ilju",1-ng+lg,bc_jm+1+l1,l2,ng)
                    a(1-ng+lg,bc_jm+1+l1,l2)=a(1+l1,bc_jm+1-ng+lg,l2)
              enddo
            endif
          endif
        endif
        ! iuju
        if(lconiu.and.lconju .eqv. .true.) then
          if(l2.ge.1.and.l2.le.bc_km) then
            if(l1.ge.0.and.l1.le.ng-1) then
              do lg=0,ng-1
                    call bound_check("iuju",bc_im+1+lg,bc_jm+1+l1,l2,ng)
                    a(bc_im+1+lg,bc_jm+1+l1,l2)=a(bc_im-l1,bc_jm-lg,l2)
              enddo
            endif
          endif
        endif
      
        ! Corners
        ! lll ! l1=i, l2=j, lg=k
        if(lconil.and.lconjl.and.lconkl .eqv. .true.) then
              if(l2.ge.0.and.l2.le.ng-1) then
              if(l1.ge.0.and.l1.le.ng-1)then
              do lg=0,ng-1
              a(1-ng+l1,1-ng+l2,1-ng+lg)=a(ng-lg,ng-l2,ng-l1)
              enddo
              endif
              endif
        endif

        ! ull ! l1=i, l2=j, lg=k
        if(lconiu.and.lconjl.and.lconkl .eqv. .true.) then
          if(l2.ge.0.and.l2.le.ng-1) then
          if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
              a(bc_im+1+l1,1-ng+l2,1-ng+lg)=a(bc_im-ng+1-lg,ng-l2,ng-l1)
            enddo
          endif
          endif
        endif

        ! lul ! l1=i, l2=j, lg=k
        if(lconil.and.lconju.and.lconkl .eqv. .true.) then
          if(l2.ge.0.and.l2.le.ng-1) then
          if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
              a(1-ng+l1,bc_jm+1+l2,1-ng+lg)=a(ng-lg,bc_jm-ng+1-l2,ng-l1)
            enddo
          endif
          endif
        endif

        ! uul ! l1=i, l2=j, lg=k
        if(lconiu.and.lconju.and.lconkl .eqv. .true.) then
          if(l2.ge.0.and.l2.le.ng-1) then
          if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
              a(bc_im+1+l1,bc_jm+1+l2,1-ng+lg)=
     &                              a(bc_im-ng+1-lg,bc_jm-ng+1-l2,ng-l1)
            enddo
          endif
          endif
        endif

        ! llu ! l1=i, l2=j, lg=k
        if(lconil.and.lconjl.and.lconku .eqv. .true.) then
          if(l2.ge.0.and.l2.le.ng-1) then
          if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
              a(1-ng+l1,1-ng+l2,bc_km+1+lg)=a(ng-lg,ng-l2,bc_km-ng+1-l1)
            enddo
          endif
          endif
        endif

        ! ulu ! l1=i, l2=j, lg=k
        if(lconiu.and.lconjl.and.lconku .eqv. .true.) then
          if(l2.ge.0.and.l2.le.ng-1) then
          if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
              a(bc_im+1+l1,1-ng+l2,bc_km+1+lg)=
     &                              a(bc_im-ng+1-lg,ng-l2,bc_km-ng+1-l1)
            enddo
          endif
          endif
        endif

        ! luu ! l1=i, l2=j, lg=k
        if(lconil.and.lconju.and.lconku .eqv. .true.) then
          if(l2.ge.0.and.l2.le.ng-1) then
          if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
              a(1-ng+l1,bc_jm+1+l2,bc_km+1+lg)=
     &                              a(ng-lg,bc_jm-ng+1-l2,bc_km-ng+1-l1)
            enddo
          endif
          endif
        endif

        ! uuu ! l1=i, l2=j, lg=k
        if(lconiu.and.lconju.and.lconku .eqv. .true.) then
          if(l2.ge.0.and.l2.le.ng-1) then
          if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
              a(bc_im+1+l1,bc_jm+1+l2,bc_km+1+lg)=
     &                      a(bc_im-ng+1-lg,bc_jm-ng+1-l2,bc_km-ng+1-l1)
            enddo
          endif
          endif
        endif

      enddo !do l1=1-ng,lh+ng
      enddo !do l2=1-ng,lh+ng

      end subroutine apply_nonperiodic_bc_r8

      ! ----------------------------------------------------------------


      subroutine boundary_condition_i4_m1(bc,a,ng)
      implicit none

      !Subroutine arguments
      !character(len=*),intent(in) :: bctype
      type(bctype),intent(in) :: bc
      integer(4),intent(inout) :: 
     &                      a(1-ng:bc_im+ng,1-ng:bc_jm+ng,1-ng:bc_km+ng)
      integer(4),intent(in) :: ng

      !Local variables
      integer(4) :: i,j,k,n1,n2
      integer(4) :: il,iu,jl,ju,kl,ku

      ! MPI parallelization variables
      integer(4) :: pidx, pidy, pidz
      integer(4) :: pidto01, pidto02, pidto03, pidto04
      integer(4) :: pidto05, pidto06, pidto07, pidto08
      integer(4) :: pidto09, pidto10, pidto11, pidto12
      integer(4) :: pidfrom01, pidfrom02, pidfrom03, pidfrom04
      integer(4) :: pidfrom05, pidfrom06, pidfrom07, pidfrom08
      integer(4) :: pidfrom09, pidfrom10, pidfrom11, pidfrom12

      integer(4) :: isend01, isend02, isend03, isend04
      integer(4) :: isend05, isend06, isend07, isend08
      integer(4) :: isend09, isend10, isend11, isend12
      integer(4) :: irecv01, irecv02, irecv03, irecv04
      integer(4) :: irecv05, irecv06, irecv07, irecv08
      integer(4) :: irecv09, irecv10, irecv11, irecv12


      !MPI_HANDLES for MPI types for non-contiguous array data communication
      ! MPI type is a kind of MPI_HANDLE in C, but it is equivalent to integer in fortran
      integer(4) :: xslab, yslab, zslab

      integer(4) :: arrsize(3)
      integer(4) :: subsize(3)
      integer(4) :: starts(3)

      integer(4) :: istatus(MPI_STATUS_SIZE), ierr
      integer(4) :: basetag

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, boundary_condition_i4_new(),'
        write(*,*)'  boundary condition grid uninitialized. Aborting.'
        call abort
      endif

      if (bc%is_initialized .eqv. .false.) then
        write(*,*)'Error, boundary_condition_i4_new(),'
        write(*,*)'  bc handle uninitialized.'
        write(*,*)'  bc name:',trim(bc%handlename)
        write(*,*)'  Aborting...'
        call abort
      endif


      ! Resister the call count to determine the base MPI_TAG value
      bc_call_count = 1 + bc_call_count
      basetag = bc_call_count*bc_mpitag_jumpstep
      bc_max_mpi_tag = max(bc_max_mpi_tag,basetag+bc_mpitag_jumpstep)


      arrsize(1) = 2*ng + bc_im
      arrsize(2) = 2*ng + bc_jm
      arrsize(3) = 2*ng + bc_km

      il = 1-ng
      iu = bc_im+ng
      jl = 1-ng
      ju = bc_jm+ng
      kl = 1-ng
      ku = bc_km+ng

      ! MPI derived data type definition for non-cnotiguous transfer
      ! xslab
      starts(1) = 0
      starts(2) = 0
      starts(3) = 0

      subsize(1) = ng
      subsize(2) = arrsize(2)
      subsize(3) = arrsize(3)
      call mpi_type_create_subarray(3, arrsize, subsize, starts,
     &  MPI_ORDER_FORTRAN, MPI_INTEGER, xslab, ierr)
      call mpi_type_commit(xslab, ierr)

      ! yslab
      subsize(1) = arrsize(1)
      subsize(2) = ng
      subsize(3) = arrsize(3)
      call mpi_type_create_subarray(3, arrsize, subsize, starts,
     &  MPI_ORDER_FORTRAN, MPI_INTEGER, yslab, ierr)
      call mpi_type_commit(yslab, ierr)

      ! zslab
      subsize(1) = arrsize(1)
      subsize(2) = arrsize(2)
      subsize(3) = ng
      call mpi_type_create_subarray(3, arrsize, subsize, starts,
     &  MPI_ORDER_FORTRAN, MPI_INTEGER, zslab, ierr)
      call mpi_type_commit(zslab, ierr)



      ! Necessary transfer among all processes - periodic BC
      !Process coordinate
      call get_pidxyz(pidx, pidy, pidz, bc_pid)

      ! Layer send/recv
      ! i planes
      if(bc_npx.gt.1)then
      pidto01=get_pid(pidx-1,pidy,pidz) ! il_s
      call MPI_ISEND(a(1,jl,kl),
     &      1, xslab, 
     &      pidto01,basetag+1,MPI_COMM_WORLD,isend01, ierr)

      pidfrom02=get_pid(pidx+1,pidy,pidz) ! iu_r
      call MPI_IRECV(a(bc_im+1,jl,kl), 
     &      1, xslab, 
     &      pidfrom02,basetag+1,MPI_COMM_WORLD,irecv02,ierr)

      pidto02=get_pid(pidx+1,pidy,pidz) ! iu_s
      call MPI_ISEND(a(bc_im-ng+1,jl,kl), 
     &      1, xslab,  
     &      pidto02,basetag+2,MPI_COMM_WORLD,isend02,ierr)      

      pidfrom01=get_pid(pidx-1,pidy,pidz) ! il_r
      call MPI_IRECV(a(il,jl,kl), 
     &      1, xslab, 
     &      pidfrom01,basetag+2,MPI_COMM_WORLD,irecv01,ierr)
      
      call MPI_WAIT(isend01,istatus,ierr)
      call MPI_WAIT(isend02,istatus,ierr)
      call MPI_WAIT(irecv01,istatus,ierr)
      call MPI_WAIT(irecv02,istatus,ierr)
      endif
      ! Done only for x direction - not complete yet

      ! j planes
      if(bc_npy.gt.1)then
      pidto03=get_pid(pidx,pidy-1,pidz) ! jl_s
      call MPI_ISEND(a(il,1,kl),
     &      1, yslab, 
     &      pidto03,basetag+3,MPI_COMM_WORLD,isend03,ierr)

      pidfrom04=get_pid(pidx,pidy+1,pidz) ! ju_r
      call MPI_IRECV(a(il,bc_jm+1,kl),
     &      1, yslab, 
     &      pidfrom04,basetag+3,MPI_COMM_WORLD,irecv04,ierr)

      pidto04=get_pid(pidx,pidy+1,pidz) ! ju_s
      call MPI_ISEND(a(il,bc_jm-ng+1,kl), 
     &      1, yslab, 
     &      pidto04,basetag+4,MPI_COMM_WORLD,isend04,ierr)

      pidfrom03=get_pid(pidx,pidy-1,pidz) ! jl_r      
      call MPI_IRECV(a(il,jl,kl), 
     &      1, yslab, 
     &      pidfrom03,basetag+4,MPI_COMM_WORLD,irecv03,ierr)

      call MPI_WAIT(isend03,istatus,ierr)
      call MPI_WAIT(isend04,istatus,ierr)
      call MPI_WAIT(irecv03,istatus,ierr)
      call MPI_WAIT(irecv04,istatus,ierr)
      endif
      ! Done for x and y direction - not completed yet

      ! k planes
      if(bc_npz.gt.1)then
      pidto05=get_pid(pidx,pidy,pidz-1) ! kl_s
      call MPI_ISEND(a(il,jl,1),
     &      1, zslab, 
     &      pidto05,basetag+5,MPI_COMM_WORLD,isend05,ierr)

      pidfrom06=get_pid(pidx,pidy,pidz+1) ! ku_r
      call MPI_IRECV(a(il,jl,bc_km+1),
     &      1, zslab, 
     &      pidfrom06,basetag+5,MPI_COMM_WORLD,irecv06,ierr)
      
      pidto06=get_pid(pidx,pidy,pidz+1) ! ku_s
      call MPI_ISEND(a(il,jl,bc_km-ng+1), 
     &      1, zslab, 
     &      pidto06,basetag+6,MPI_COMM_WORLD,isend06,ierr)

      pidfrom05=get_pid(pidx,pidy,pidz-1) ! kl_r
      call MPI_IRECV(a(il,jl,kl), 
     &      1, zslab, 
     &      pidfrom05,basetag+6,MPI_COMM_WORLD,irecv05,ierr)

      call MPI_WAIT(isend05,istatus,ierr)
      call MPI_WAIT(isend06,istatus,ierr)
      call MPI_WAIT(irecv05,istatus,ierr)
      call MPI_WAIT(irecv06,istatus,ierr)
      endif

      call apply_periodic_bc_i4(a, ng)
      ! Done for all directions - complete

       
      ! Overwrite the boundary data as per given BC
      call apply_adiabatic_bc_i4(bc, a, ng)
      

      call mpi_barrier(MPI_COMM_WORLD, ierr)
      call mpi_type_free(xslab, ierr)
      call mpi_type_free(yslab, ierr)
      call mpi_type_free(zslab, ierr)

      end subroutine boundary_condition_i4_m1

      ! ----------------------------------------------------------------


      subroutine apply_adiabatic_bc_i4(bc, a, ng)
      implicit none
      !Subroutine arguments
      !character(len=*),intent(in) :: bctype
      type(bctype),intent(in) :: bc
      integer(4),intent(inout) :: 
     &                      a(1-ng:bc_im+ng,1-ng:bc_jm+ng,1-ng:bc_km+ng)
      integer(4),intent(in) :: ng

      !Local variables
      integer(4) :: i,j,k,lg
      integer(4) :: l1, l2, lh
      logical :: gconxn, gconyn, gconzn
      logical :: gconxp, gconyp, gconzp
      logical :: lconil, lconiu
      logical :: lconjl, lconju
      logical :: lconkl, lconku
      integer(4) :: pidx,pidy,pidz

      if ( bc_initialized .eqv. .false.) then
        write(*,*)'Error, apply_adiabatic_bc_i4(),'
        write(*,*)'  boundary condition grid uninitialized. Aborting.'
        call abort
      endif

      if (bc%is_initialized .eqv. .false.) then
        write(*,*)'Error, apply_adiabatic_bc_i4(),'
        write(*,*)'  bc handle uninitialized.'
        write(*,*)'  bc name:',trim(bc%handlename)
        write(*,*)'  Aborting...'
        call abort
      endif


      gconxn=.false.; gconyn=.false.; gconzn=.false.
      gconxp=.false.; gconyp=.false.; gconzp=.false.
      lconil=.false.; lconjl=.false.; lconkl=.false.
      lconiu=.false.; lconju=.false.; lconku=.false.

      if(trim(bc%xn) .ne. 'PERIODIC') gconxn = .true.
      if(trim(bc%xp) .ne. 'PERIODIC') gconxp = .true.
      if(trim(bc%yn) .ne. 'PERIODIC') gconyn = .true.
      if(trim(bc%yp) .ne. 'PERIODIC') gconyp = .true.
      if(trim(bc%zn) .ne. 'PERIODIC') gconzn = .true.
      if(trim(bc%zp) .ne. 'PERIODIC') gconzp = .true.

      call get_pidxyz(pidx, pidy, pidz, bc_pid)
      if(pidx.eq.0)     lconil = .true.
      if(pidx.eq.bc_npx-1) lconiu = .true.
      if(pidy.eq.0)     lconjl = .true.
      if(pidy.eq.bc_npy-1) lconju = .true.
      if(pidz.eq.0)     lconkl = .true.
      if(pidz.eq.bc_npz-1) lconku = .true.

      ! Combining global condition and local condition
      lconil=lconil.and.gconxn
      lconiu=lconiu.and.gconxp
      lconjl=lconjl.and.gconyn
      lconju=lconju.and.gconyp
      lconkl=lconkl.and.gconzn
      lconku=lconku.and.gconzp

      ! Terminate if this local process is nothing to do with BC
      if(lconil.or.lconiu.or.lconjl.or.lconju.or.lconkl.or.lconku
     &    .eqv. .false.) then
            return
      endif
 
      lh = max(bc_im,bc_jm,bc_km)  ! The highest # of grids

      do l2=1-ng,lh+ng
      do l1=1-ng,lh+ng

      ! x-layers: l1=j, l2=k
      if(lconil.or.lconiu .eqv. .true.)then
        if(l1.le.bc_jm+ng.and.l2.le.bc_km+ng)then  
          do lg=1,ng
            ! il
            if(lconil.eqv..true.)a(1-lg,l1,l2)=a(lg,l1,l2)
            ! iu
            if(lconiu.eqv..true.)a(bc_im+lg,l1,l2)=a(bc_im+1-lg,l1,l2)
          enddo
        endif
      endif !if(lconil.or.lconiu .eqv. .true.)then

      ! y-layers: l1=k, l2=i
      if(lconjl.or.lconju .eqv. .true.)then
        if(l1.le.bc_km+ng.and.l2.le.bc_im+ng)then  
          do lg=1,ng
            ! jl
            if(lconjl.eqv..true.)a(l2,1-lg,l1)=a(l2,lg,l1)
            ! ju
            if(lconju.eqv..true.)a(l2,bc_jm+lg,l1)=a(l2,bc_jm+1-lg,l1)
          enddo
        endif
      endif !if(lconjl.or.lconju .eqv. .true.)then

      ! z-layers: l1=i, l2=j
      if(lconkl.or.lconku .eqv. .true.)then
        if(l1.le.bc_im+ng.and.l2.le.bc_jm+ng)then  
          do lg=1,ng
            ! kl
            if(lconkl.eqv..true.)a(l1,l2,1-lg)=a(l1,l2,lg)
            ! ku
            if(lconku.eqv..true.)a(l1,l2,bc_km+lg)=a(l1,l2,bc_km+1-lg)
          enddo
        endif
      endif !if(lconkl.or.lconku .eqv. .true.)then
      

      ! Edges along x-axis: l2=i, l1=j, lg=k
      ! jlkl
      if(lconjl.and.lconkl .eqv. .true.) then
        if(l2.ge.1.and.l2.le.bc_im) then
          if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("jlkl",l2,1-ng+lg,1-ng+l1,ng)
                  a(l2,1-ng+lg,1-ng+l1)=a(l2,ng-l1,ng-lg)
            enddo
          endif
        endif
      endif
      ! jukl
      if(lconju.and.lconkl .eqv. .true.) then
        if(l2.ge.1.and.l2.le.bc_im) then
          if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("jukl",l2,bc_jm+1+lg,1-ng+l1,ng)
                  a(l2,bc_jm+1+lg,1-ng+l1)=a(l2,bc_jm+1-ng+l1,1+lg)
            enddo
          endif
        endif
      endif
      ! jlku
      if(lconjl.and.lconku .eqv. .true.) then
        if(l2.ge.1.and.l2.le.bc_im) then
          if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("jlku",l2,1-ng+lg,bc_km+1+l1,ng)
                  a(l2,1-ng+lg,bc_km+1+l1)=a(l2,1+l1,bc_km+1-ng+lg)
            enddo
          endif
        endif
      endif
      ! juku
      if(lconju.and.lconku .eqv. .true.) then
        if(l2.ge.1.and.l2.le.bc_im) then
          if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("juku",l2,bc_jm+1+lg,bc_km+1+l1,ng)
                  a(l2,bc_jm+1+lg,bc_km+1+l1)=a(l2,bc_jm-l1,bc_km-lg)
            enddo
          endif
        endif
      endif

      
      ! Edges along y-axis: l2=j, l1=k, lg=i
      ! klil
      if(lconkl.and.lconil .eqv. .true.) then
        if(l2.ge.1.and.l2.le.bc_jm) then
          if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("klil",1-ng+l1,l2,1-ng+lg,ng)
                  a(1-ng+l1,l2,1-ng+lg)=a(ng-lg,l2,ng-l1)
            enddo
          endif
        endif
      endif
      ! kuil
      if(lconku.and.lconil .eqv. .true.) then
        if(l2.ge.1.and.l2.le.bc_jm) then
          if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("kuil",bc_im+1+l1,l2,bc_km+1+lg,ng)
                  a(bc_im+1+l1,l2,bc_km+1+lg)=a(bc_im-lg,l2,bc_km-l1)
            enddo
          endif
        endif
      endif
      ! kliu
      if(lconkl.and.lconiu .eqv. .true.) then
        if(l2.ge.1.and.l2.le.bc_jm) then
          if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("kliu",1-ng+l1,l2,bc_km+1+lg,ng)
                  a(1-ng+l1,l2,bc_km+1+lg)=a(1+lg,l2,bc_km+1-ng+l1)
            enddo
          endif
        endif
      endif
      ! kuiu
      if(lconku.and.lconiu .eqv. .true.) then
        if(l2.ge.1.and.l2.le.bc_jm) then
          if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("kuiu",bc_im+1+l1,l2,bc_km+1+lg,ng)
                  a(bc_im+1+l1,l2,bc_km+1+lg)=a(bc_im-lg,l2,bc_km-l1)
            enddo
          endif
        endif
      endif
      
      ! Edges along z-axis: l2=k, l1=i, lg=j
      ! iljl
      if(lconil.and.lconjl .eqv. .true.) then
        if(l2.ge.1.and.l2.le.bc_km) then
          if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("iljl",1-ng+lg,1-ng+l1,l2,ng)
                  a(1-ng+lg,1-ng+l1,l2)=a(ng-l1,ng-lg,l2)
            enddo
          endif
        endif
      endif
      ! iujl
      if(lconiu.and.lconjl .eqv. .true.) then
        if(l2.ge.1.and.l2.le.bc_km) then
          if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("iujl",bc_im+1+lg,1-ng+l1,l2,ng)
                  a(bc_im+1+lg,1-ng+l1,l2)=a(bc_im+1-ng+l1,1+lg,l2)
            enddo
          endif
        endif
      endif
      ! ilju
      if(lconil.and.lconju .eqv. .true.) then
        if(l2.ge.1.and.l2.le.bc_km) then
          if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("ilju",1-ng+lg,bc_jm+1+l1,l2,ng)
                  a(1-ng+lg,bc_jm+1+l1,l2)=a(1+l1,bc_jm+1-ng+lg,l2)
            enddo
          endif
        endif
      endif
      ! iuju
      if(lconiu.and.lconju .eqv. .true.) then
        if(l2.ge.1.and.l2.le.bc_km) then
          if(l1.ge.0.and.l1.le.ng-1) then
            do lg=0,ng-1
                  call bound_check("iuju",bc_im+1+lg,bc_jm+1+l1,l2,ng)
                  a(bc_im+1+lg,bc_jm+1+l1,l2)=a(bc_im-l1,bc_jm-lg,l2)
            enddo
          endif
        endif
      endif
      
      ! Corners
      ! lll ! l1=i, l2=j, lg=k
      if(lconil.and.lconjl.and.lconkl .eqv. .true.) then
        if(l2.ge.0.and.l2.le.ng-1) then
          if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
              a(1-ng+l1,1-ng+l2,1-ng+lg)=a(ng-lg,ng-l2,ng-l1)
            enddo
          endif
        endif
      endif

      ! ull ! l1=i, l2=j, lg=k
      if(lconiu.and.lconjl.and.lconkl .eqv. .true.) then
        if(l2.ge.0.and.l2.le.ng-1) then
          if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
              a(bc_im+1+l1,1-ng+l2,1-ng+lg)=a(bc_im-ng+1-lg,ng-l2,ng-l1)
            enddo
          endif
        endif
      endif

      ! lul ! l1=i, l2=j, lg=k
      if(lconil.and.lconju.and.lconkl .eqv. .true.) then
        if(l2.ge.0.and.l2.le.ng-1) then
          if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
              a(1-ng+l1,bc_jm+1+l2,1-ng+lg)=a(ng-lg,bc_jm-ng+1-l2,ng-l1)
            enddo
          endif
        endif
      endif

      ! uul ! l1=i, l2=j, lg=k
      if(lconiu.and.lconju.and.lconkl .eqv. .true.) then
        if(l2.ge.0.and.l2.le.ng-1) then
          if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
              a(bc_im+1+l1,bc_jm+1+l2,1-ng+lg)=
     &                              a(bc_im-ng+1-lg,bc_jm-ng+1-l2,ng-l1)
            enddo
          endif
        endif
      endif

      ! llu ! l1=i, l2=j, lg=k
      if(lconil.and.lconjl.and.lconku .eqv. .true.) then
        if(l2.ge.0.and.l2.le.ng-1) then
          if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
              a(1-ng+l1,1-ng+l2,bc_km+1+lg)=a(ng-lg,ng-l2,bc_km-ng+1-l1)
            enddo
          endif
        endif
      endif

      ! ulu ! l1=i, l2=j, lg=k
      if(lconiu.and.lconjl.and.lconku .eqv. .true.) then
        if(l2.ge.0.and.l2.le.ng-1) then
          if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
              a(bc_im+1+l1,1-ng+l2,bc_km+1+lg)=
     &                              a(bc_im-ng+1-lg,ng-l2,bc_km-ng+1-l1)
            enddo
          endif
        endif
      endif

      ! luu ! l1=i, l2=j, lg=k
      if(lconil.and.lconju.and.lconku .eqv. .true.) then
        if(l2.ge.0.and.l2.le.ng-1) then
          if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
              a(1-ng+l1,bc_jm+1+l2,bc_km+1+lg)=
     &                              a(ng-lg,bc_jm-ng+1-l2,bc_km-ng+1-l1)
            enddo
          endif
        endif
      endif

      ! uuu ! l1=i, l2=j, lg=k
      if(lconiu.and.lconju.and.lconku .eqv. .true.) then
        if(l2.ge.0.and.l2.le.ng-1) then
          if(l1.ge.0.and.l1.le.ng-1)then
            do lg=0,ng-1
              a(bc_im+1+l1,bc_jm+1+l2,bc_km+1+lg)=
     &                      a(bc_im-ng+1-lg,bc_jm-ng+1-l2,bc_km-ng+1-l1)
            enddo
          endif
        endif
      endif

      enddo !do l1=1-ng,lh+ng
      enddo !do l2=1-ng,lh+ng

      end subroutine apply_adiabatic_bc_i4


      end module boundary_conditions