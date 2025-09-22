      module interpolation
      implicit none

      !-----------------------------------------------------------------
      contains
      real(8) function interpolation_L1_1D(x, x0,f0, x1,f1, dbgtag)
      implicit none
      real(8),intent(in) :: x
      real(8),intent(in) :: x0, f0
      real(8),intent(in) :: x1, f1
      character(len=*),intent(in) :: dbgtag

      real(8) :: xi
      real(8) :: sf0, sf1

      xi = 2.0*(x - x0)/(x1 - x0) - 1.0

      ! Sanitation check
      if(xi.lt.-1.0) then
        write(*,*)'dbgtag: ',dbgtag
        write(*,*)'  FATAL ERROR: interpolation_L1_1D: invalid xi.'
        write(*,*)'               Double-check the x coordinate inputs.'
        write(*,*)'               x: ',x
        write(*,*)'               xi: ',xi
        write(*,*)'               Aborting'
        write(*,*)'            x0 = ',x0
        write(*,*)'            x1 = ',x1
        call abort
      endif
      if(xi.gt.1.0) then
        write(*,*)'dbgtag: ',dbgtag
        write(*,*)'  FATAL ERROR: interpolation_L1_1D: invalid xi.'
        write(*,*)'               Double-check the x coordinate inputs.'
        write(*,*)'               x: ',x
        write(*,*)'               xi: ',xi
        write(*,*)'               Aborting'
        call abort
      endif

      sf0 = 0.5 * (1.0 - xi) * f0
      sf1 = 0.5 * (1.0 + xi) * f1

      interpolation_L1_1D = sf0 + sf1
      return 
      end function interpolation_L1_1D
      
      !-----------------------------------------------------------------

      real(8) function interpolation_L1_2D(x, y,
     &                                    x0,y0,f0, 
     &                                    x1,y1,f1, 
     &                                    x2,y2,f2, 
     &                                    x3,y3,f3 , dbgtag)
      implicit none
      real(8),intent(in) :: x, y
      real(8),intent(in) :: x0, y0, f0
      real(8),intent(in) :: x1, y1, f1
      real(8),intent(in) :: x2, y2, f2
      real(8),intent(in) :: x3, y3, f3 
      character(len=*),intent(in) :: dbgtag
      !Assumed 4 coordinates are ordered properly.
      !Assumed the mesh geometry is rectilinear, 
      !                             expecting x2=x0, x3=x1, y2=y0, y3=y1
      ! @ (x0, y0), (xi,aeta) = (-1,-1)
      ! @ (x1, y1), (xi,aeta) = (+1,-1)
      ! @ (x2, y2), (xi,aeta) = (-1,+1)
      ! @ (x3, y3), (xi,aeta) = (+1,+1)

      real(8) :: xi, aeta !parametric coordinates for interpolation
      real(8) :: sf0, sf1, sf2, sf3 !shape function values

      xi = 2.0*(x - x0)/(x3 - x0) - 1.0
      aeta = 2.0*(y - y0)/(y3 - y0) - 1.0

      ! Sanitation check
      if(xi.lt.-1.0) then
        write(*,*)'dbgtag: ',dbgtag
        write(*,*)'  FATAL ERROR: interpolation_L1_2D: invalid xi.'
        write(*,*)'               Double-check the x coordinate inputs.'
        write(*,*)'               x: ',x
        write(*,*)'               xi: ',xi
        write(*,*)'  Data grid coordinates: '
        write(*,*)'            x0,y0 = ',x0,y0
        write(*,*)'            x1,y1 = ',x1,y1
        write(*,*)'            x2,y2 = ',x2,y2
        write(*,*)'            x3,y3 = ',x3,y3
        write(*,*)'               Aborting'
        call abort
      endif
      if(xi.gt.1.0) then
        write(*,*)'dbgtag: ',dbgtag
        write(*,*)'  FATAL ERROR: interpolation_L1_2D: invalid xi.'
        write(*,*)'               Double-check the x coordinate inputs.'
        write(*,*)'               x: ',x
        write(*,*)'               xi: ',xi
        write(*,*)'            x0,y0 = ',x0,y0
        write(*,*)'            x1,y1 = ',x1,y1
        write(*,*)'            x2,y2 = ',x2,y2
        write(*,*)'            x3,y3 = ',x3,y3
        write(*,*)'               Aborting'
        call abort
      endif
      if(aeta.lt.-1.0) then
        write(*,*)'dbgtag: ',dbgtag
        write(*,*)'  FATAL ERROR: interpolation_L1_2D: invalid aeta.'
        write(*,*)'               Double-check the y coordinate inputs.'
        write(*,*)'               y: ',y
        write(*,*)'               aeta: ',aeta
        write(*,*)'            x0,y0 = ',x0,y0
        write(*,*)'            x1,y1 = ',x1,y1
        write(*,*)'            x2,y2 = ',x2,y2
        write(*,*)'            x3,y3 = ',x3,y3
        write(*,*)'               Aborting'
        call abort
      endif
      if(aeta.gt.1.0) then
        write(*,*)'dbgtag: ',dbgtag
        write(*,*)'  FATAL ERROR: interpolation_L1_2D: invalid aeta.'
        write(*,*)'               Double-check the y coordinate inputs.'
        write(*,*)'               x: ',y
        write(*,*)'               xi: ',aeta
        write(*,*)'            x0,y0 = ',x0,y0
        write(*,*)'            x1,y1 = ',x1,y1
        write(*,*)'            x2,y2 = ',x2,y2
        write(*,*)'            x3,y3 = ',x3,y3
        write(*,*)'               Aborting'
        call abort
      endif

      sf0 = 0.25 * (1.0 - xi)*(1.0 - aeta) * f0
      sf1 = 0.25 * (1.0 + xi)*(1.0 - aeta) * f1
      sf2 = 0.25 * (1.0 - xi)*(1.0 + aeta) * f2
      sf3 = 0.25 * (1.0 + xi)*(1.0 + aeta) * f3
      
      interpolation_L1_2D = sf0 + sf1 + sf2 + sf3
      return
      end function interpolation_L1_2D

      !-----------------------------------------------------------------

      real(8) function interpolation_L1_3D(x, y, z,
     &                                    x0,y0,z0,f0, 
     &                                    x1,y1,z1,f1,
     &                                    x2,y2,z2,f2, 
     &                                    x3,y3,z3,f3,
     &                                    x4,y4,z4,f4,
     &                                    x5,y5,z5,f5,
     &                                    x6,y6,z6,f6,
     &                                    x7,y7,z7,f7, dbgtag )
      implicit none
      real(8),intent(in) :: x, y, z
      real(8),intent(in) :: x0,y0,z0,f0
      real(8),intent(in) :: x1,y1,z1,f1
      real(8),intent(in) :: x2,y2,z2,f2
      real(8),intent(in) :: x3,y3,z3,f3 
      real(8),intent(in) :: x4,y4,z4,f4
      real(8),intent(in) :: x5,y5,z5,f5
      real(8),intent(in) :: x6,y6,z6,f6 
      real(8),intent(in) :: x7,y7,z7,f7 
      character(len=*),intent(in) :: dbgtag
      !Assumed 8 coordinates are ordered properly.
      !Assumed the mesh geometry is rectilinear, 
      ! @ (x0, y0, z0), (xi,aeta,zeta) = (-1,-1,-1)
      ! @ (x1, y1, z1), (xi,aeta,zeta) = (+1,-1,-1)
      ! @ (x2, y2, z2), (xi,aeta,zeta) = (-1,+1,-1)
      ! @ (x3, y3, z3), (xi,aeta,zeta) = (+1,+1,-1)
      ! @ (x4, y4, z4), (xi,aeta,zeta) = (-1,-1,+1)
      ! @ (x5, y5, z5), (xi,aeta,zeta) = (+1,-1,+1)
      ! @ (x6, y6, z6), (xi,aeta,zeta) = (-1,+1,+1)
      ! @ (x7, y7, z7), (xi,aeta,zeta) = (+1,+1,+1)

      real(8) :: xi, aeta, zeta !parametric coordinates for interpolation
      real(8) :: sf0,sf1,sf2,sf3,sf4,sf5,sf6,sf7 !shape function values

      xi = 2.0*(x - x0)/(x7 - x0) - 1.0
      aeta = 2.0*(y - y0)/(y7 - y0) - 1.0
      zeta = 2.0*(z - z0)/(z7 - z0) - 1.0

      ! Sanitation check
      if(xi.lt.-1.0) then
        write(*,*)'dbgtag: ',dbgtag
        write(*,*)'  FATAL ERROR: interpolation_L1_3D: invalid xi.'
        write(*,*)'               Double-check the x coordinate inputs.'
        write(*,*)'               x: ',x
        write(*,*)'               xi: ',xi
        write(*,*)'               Aborting'
        call abort
      endif
      if(xi.gt.1.0) then
        write(*,*)'dbgtag: ',dbgtag
        write(*,*)'  FATAL ERROR: interpolation_L1_3D: invalid xi.'
        write(*,*)'               Double-check the x coordinate inputs.'
        write(*,*)'               x: ',x
        write(*,*)'               xi: ',xi
        write(*,*)'               Aborting'
        call abort
      endif
      if(aeta.lt.-1.0) then
        write(*,*)'dbgtag: ',dbgtag
        write(*,*)'  FATAL ERROR: interpolation_L1_3D: invalid aeta.'
        write(*,*)'               Double-check the y coordinate inputs.'
        write(*,*)'               y: ',y
        write(*,*)'               aeta: ',aeta
        write(*,*)'               Aborting'
        call abort
      endif
      if(aeta.gt.1.0) then
        write(*,*)'dbgtag: ',dbgtag
        write(*,*)'  FATAL ERROR: interpolation_L1_3D: invalid aeta.'
        write(*,*)'               Double-check the y coordinate inputs.'
        write(*,*)'               y: ',y
        write(*,*)'               aeta: ',aeta
        write(*,*)'               Aborting'
        call abort
      endif
      if(zeta.lt.-1.0) then
        write(*,*)'dbgtag: ',dbgtag
        write(*,*)'  FATAL ERROR: interpolation_L1_3D: invalid aeta.'
        write(*,*)'               Double-check the z coordinate inputs.'
        write(*,*)'               z: ',z
        write(*,*)'               zeta: ',zeta
        write(*,*)'               Aborting'
        call abort
      endif
      if(zeta.gt.1.0) then
        write(*,*)'dbgtag: ',dbgtag
        write(*,*)'  FATAL ERROR: interpolation_L1_3D: invalid aeta.'
        write(*,*)'               Double-check the z coordinate inputs.'
        write(*,*)'               z: ',z
        write(*,*)'               zeta: ',zeta
        write(*,*)'               Aborting'
        call abort
      endif

      sf0 = 0.125 * (1.0 - xi)*(1.0 - aeta)*(1.0 - zeta) * f0
      sf1 = 0.125 * (1.0 + xi)*(1.0 - aeta)*(1.0 - zeta) * f1
      sf2 = 0.125 * (1.0 - xi)*(1.0 + aeta)*(1.0 - zeta) * f2
      sf3 = 0.125 * (1.0 + xi)*(1.0 + aeta)*(1.0 - zeta) * f3
      sf4 = 0.125 * (1.0 - xi)*(1.0 - aeta)*(1.0 + zeta) * f4
      sf5 = 0.125 * (1.0 + xi)*(1.0 - aeta)*(1.0 + zeta) * f5
      sf6 = 0.125 * (1.0 - xi)*(1.0 + aeta)*(1.0 + zeta) * f6
      sf7 = 0.125 * (1.0 + xi)*(1.0 + aeta)*(1.0 + zeta) * f7
      
      interpolation_L1_3D = sf0 + sf1 + sf2 + sf3 +
     &                      sf4 + sf5 + sf6 + sf7
      return
      end function interpolation_L1_3D

      !-----------------------------------------------------------------

      subroutine test_interpolation()
      implicit none
      real(8) :: x,y,z,f
      real(8) :: x0, x1, y0, y1, z0, z1
      real(8) :: f0, f1, f2, f3, f4, f5, f6, f7

      x = 1.0
      y = 1.0
      x0 = 0.0
      x1 = 1.0
      y0 = 0.0
      y1 = 1.0
      f0 = -2.0
      f1 = 1.0
      f2 = 5.0
      f3 = 0.0

      f = interpolation_L1_1D(x,
     &                        x0,f0,  
     &                        x1,f1, 'test1d')
      write(*,*)'interpolation 1D test value:',f

      f = interpolation_L1_2D(x,y,
     &                        x0,y0,f0,  
     &                        x1,y0,f1,
     &                        x0,y1,f2,
     &                        x1,y1,f3, 'test2d')

      write(*,*)'interpolation 2D test value:',f

      z = 1.0
      z0 = 0.0
      z1 = 1.0
      f4 = 3.0
      f5 = 2.0
      f6 = -10.0
      f7 = -5.0

      f = interpolation_L1_3D(x,y,z,
     &                        x0,y0,z0,f0,  
     &                        x1,y0,z0,f1,
     &                        x0,y1,z0,f2,
     &                        x1,y1,z0,f3,
     &                        x0,y0,z1,f4,
     &                        x1,y0,z1,f5,
     &                        x0,y1,z1,f6,
     &                        x1,y1,z1,f7, 'test1d' )

      write(*,*)'interpolation 3D test value:',f

      return
      end subroutine test_interpolation

      !-----------------------------------------------------------------

      end module interpolation