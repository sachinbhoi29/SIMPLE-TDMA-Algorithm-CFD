!     Two dimensional pressure correction assignment 
program SIMPLE_TDMA_2D
! variables
! dimensions:
      integer :: iNxCells, iNyCells, iNxNodes, iNyNodes, iNxUNodes, iNyUnodes, & 
                 iNxVNodes,iNyVNodes,iNxPNodes,iNyPNodes, iNxMax, iNyMax
!iterations
      integer :: iMaxIter, i,j,n
!convergence limit
      real :: rConvergence
!residuals
      real :: rURes, rVRes, rPRes, rMaxRes, rCRes
!underrelaxation
      real :: rAlphaP, rAlphaU, rAlphaV, rAlphaC
!Dimensionless numbers
      real :: rReynolds, rPeclet, rInjStart,rInjEnd
!massflow
      real :: rTotFlow, rMRatio, rCInFlow, rCOutFlow
!grid
      real :: rHeight, rWidth, rDX, rDY
! current auxiliary variables
      real :: rPCur, rUCur,rVCur,rCCur
! storage for variables
      real, allocatable, dimension (:,:) :: rU,rV,rP,rC,rDU,rDV,rDP


! read parameters from file
      
      open (unit=1, file="parameters.par",status="old")
      read (1,*) rHeight
      read (1,*) rWidth
      read (1,*) rReynolds
      read (1,*) rPeclet
      read (1,*) rInjStart
      read (1,*) rInjEnd
      read (1,*) iNxCells
      read (1,*) iNyCells
      read (1,*) iMaxIter
      read (1,*) rConvergence
      read (1,*) rAlphaP
      read (1,*) rAlphaU
      read (1,*) rAlphaV
      read (1,*) rAlphaC
      close(1)
! grid cells
      rDx=rWidth/iNxCells
      rDy=rHeight/iNyCells
! grid nodes
      iNxNodes=iNxCells+1
      iNyNodes=iNyCells+1
! nodes for P - control volumes + 2 for BCs
      iNxPNodes=iNxCells+2
      iNyPNodes=iNyCells+2
! nodes for U - control volumes + 1 in X and 2 in Y
      iNxUNodes=iNxCells+1
      iNyUNodes=iNyCells+2
! nodes for V - control volumes + 2 in X and 1 in Y
      iNxVNodes=iNxCells+2
      iNyVNodes=iNyCells+1
!
      iNxMax=iNxCells+2
      iNyMax=iNyCells+2
! allocate memory
! scalar - same control volume as P
      allocate(rU(iNxMax,iNyMax),rV(iNxMax,iNyMax),rP(iNxMax,iNyMax),rC(iNxMax,iNyMax))
      allocate(rDU(iNxMax,iNyMax),rDV(iNxMax,iNyMax),rDP(iNxMax,iNyMax))
! Initialise
      rP=0.
      rU=0.
      rV=0.
      rC=0.
      rMaxRes=10.e+3
      n=1
! start time loop   ! It is the stopping criteria; imaxiter is maximum steps and reaxres is the residual
      do while(n.le.iMaxIter.and.rMaxRes.gt.rConvergence)
! compute u-momentum equation
       call UMomentum(iNxUNodes, iNyUNodes,iNxMax,iNyMax,rDx,rDy,rU,rV,rP,rDU,rAlphaU,rURes,rReynolds)
! compute v-momentum equation
       call VMomentum(iNxVNodes, iNyVNodes,iNxMax,iNyMax,rDx,rDy,rU,rV,rP,rDV,rAlphaV,rVRes,rReynolds)
! correct mass flux
       rTotFlow=0.
       do j = 2 , iNyUNodes-1
          rTotFlow=rTotFlow+rDY*rU(iNxUNodes,j)  
       end do
       if(rTotFlow.gt.0.0) then
            rMRatio=rHeight/rTotFlow
         else
            rMRatio=1.0
         end if
       do j = 2 , iNyUNodes-1
          rU(iNxUNodes,j)=rU(iNxUNodes,j)*rMRatio
       end do
! check scalar flux
       rCOutFlow=0.
       do j = 2 , iNyUNodes-1
          rCOutFlow=rCOutFlow+rDY*rC(iNxUNodes,j)  
       end do
! check scalar flux
       rCInFlow=0.
       do j = 2 , iNyUNodes-1
          rCInFlow=rCInFlow+rDY*rC(2,j)  
       end do
!       print *, "Scalar: ", rCInFlow,rCOutFlow
! compute pressure correction
       call PCorrection(iNxPNodes, iNyPNodes,iNxMax,iNyMax,rDx,rDy,rU,rV,rDU,rDV,rDP,rPRes)
! correct velocities 
       do i=2,iNxUNodes-1
        do j=2,iNyUNodes-1
          rU(i,j)=rU(i,j)-rDU(i,j)*(rDP(i+1,j)-rDP(i,j))
        end do
       end do
       do i=2,iNxVNodes-1
        do j=2,iNyVNodes-1
          rV(i,j)=rV(i,j)-rDV(i,j)*(rDP(i,j+1)-rDP(i,j))
        end do
       end do
! correct pressure
       rP=rP+rAlphaP*rDP
! scale pressure
       rP=rP-rP(2,2)
! Solve the scalar equation for corrected flow field:
       call Scalar(iNxPNodes, iNyPNodes,iNxMax,iNyMax,rDx,rDy,rU,rV,rC,rAlphaC,rCRes,rPeclet,rInjStart,rInjEnd)
! maximum residual
       rMaxRes=max(rURes,rVRes,rPRes,rCRes)
       write(*,"(2X,I4,2X,5(2X,E10.3))") n,rURes,rVRes,rPRes,rCRes,rP(2,2)
! advance iterations
       n=n+1
      end do
! check convergence3
      if(rMaxRes.gt.rConvergence) then 
       write(*,*) "Stat: Convergence not reached. Exiting..."
      else
       write(*,*) "Stat: Converged."
      end if
!output tecplot file:
      open (unit=1,file="resultn.dat",status="replace",action="write")
      write(1,*) 'Variables="X","Y","P","U","V", "Scalar"'
      write(1,'("ZONE DATAPACKING=POINT, I=",I4,", J=",I4)') iNxNodes,iNyNodes
      do j=1,iNyNodes
       do i=1,iNxNodes
         rPCur=0.25*(rP(i,j)+rP(i+1,j)+rP(i,j+1)+rP(i+1,j+1))
         rCCur=0.25*(rC(i,j)+rC(i+1,j)+rC(i,j+1)+rC(i+1,j+1))
         rUCur=0.5 *(rU(i,j)+rU(i,j+1))
         rVCur=0.5 *(rV(i,j)+rV(i+1,j))
         write(1,"(8E15.7)") (i-1)*rDX,(j-1)*rDY,rPCur,rUCur,rVCur,rCCur
       end do
      end do
      close(1)
end program SIMPLE_TDMA_2D


subroutine UMomentum(iNxUNodes, iNyUNodes,iNxMax,iNyMax,rDx,rDy,rU,rV,rP,rDU,rAlphaU,rURes,rReynolds)
!arguments
    integer, intent(in) :: iNxUNodes, iNyUNodes,iNxMax,iNyMax
    real, intent(in)    :: rDx,rDy, rAlphaU,rReynolds
    real, intent(inout) :: rURes
    real, dimension (iNxMax, iNyMax), intent (in) :: rV,rP
    real, dimension (iNxMax, iNyMax), intent (inout) :: rU,rDU
! local variables
    real ::  de, dw, dn, ds, ge, gw, gn, gs,rTemp
    real, dimension(iNxMax,iNyMax) ::aw,ae,as,an,ap,su
    integer :: i,j

! set boundary conditions for coefficients in the equations
! south
    aw(:,1)=0.
    ae(:,1)=0.
    as(:,1)=0.
    an(:,1)=0.
    ap(:,1)=1.
    su(:,1)=0.
! west
    aw(1,:)=0.
    ae(1,:)=0.
    as(1,:)=0.
    an(1,:)=0.
    ap(1,:)=1.
    su(1,:)=0.
! north
    aw(:,iNyUNodes)=0.
    ae(:,iNyUNodes)=0.
    as(:,iNyUNodes)=0.
    an(:,iNyUNodes)=0.
    ap(:,iNyUNodes)=1.
    su(:,iNyUNodes)=1.
! east
    aw(iNxUNodes,:)=0.
    ae(iNxUNodes,:)=0.
    as(iNxUNodes,:)=0.
    an(iNxUNodes,:)=0.
    ap(iNxUNodes,:)=1.
    su(iNxUNodes,:)=0.
! set the coefficeints in the interior staggered cells

      do i = 2 , iNxUNodes-1
         do J = 2 , iNyUNodes-1
            dw = rDY/(rDX*rReynolds)
            gw = rDY*(rU(i-1, j) + rU(i, j)) / 2.
            de = rDY/(rDX*rReynolds)
            ge = rDY*(rU(i, j) + rU(i+1, j)) / 2.
            if ( j .eq. 2 ) then
               ds = 2.*rDX/(rDY*rReynolds)
               gs = 0.0
            else
               ds = rDX/(rDY*rReynolds)
               gs = rDX*(rV(i, j-1) + rV(i+1, j-1)) / 2.
            end if 
            if ( j .eq. (iNyUNodes - 1)) then   
               dn = 2.*rDX/(rDY*rReynolds)
               gn = 0.0
            else   
               dn = rDX/(rDY*rReynolds)
               gn = rDX*(rV(i, j) + rV(i+1, j)) / 2.
            end if
            aw(i, j) = dw + max( 0.0 ,  gw)
            ae(i, j) = de + max( 0.0 , -ge)
            as(i, j) = ds + max( 0.0 ,  gs)
            an(i, j) = dn + max( 0.0 , -gn)
            ap(i, j) = dw+de+ds+dn+max( 0.0 ,  -gw) + max( 0.0 , ge) +max( 0.0 ,  -gs)+max( 0.0 ,  gn)
            su(i, j) = -(rP(i+1, j) - rP(i, j))*rDY
         end do
      end do


! evaluate the residual in the interior
      rURes=0.
      do i = 2 , iNxUNodes - 1
         do j = 2 , iNyUNodes - 1 
           rTemp = abs( ap(i, j)*rU(i, j) & 
                      - ae(i, j)*rU(i+1, j) - aw(i, j)*rU(i-1, j) &
                      - an(i, j)*rU(i, j+1) - as(i, j)*rU(i, j-1) &
                      - su(i, j))
           if(rTemp.ge.rURes) rURes=rTemp
         end do
      end do
! apply underrelaxation to interior cells
      ap(2:iNxUNodes-1,2:iNyUNodes-1)=ap(2:iNxUNodes-1,2:iNyUNodes-1)/rAlphaU
      su(2:iNxUNodes-1,2:iNyUNodes-1)=su(2:iNxUNodes-1,2:iNyUNodes-1)+ &
         (1.-rAlphaU)*ap(2:iNxUNodes-1,2:iNyUNodes-1)*rU(2:iNxUNodes-1,2:iNyUNodes-1)
! store dU for pressure correction:
      rDU(2:iNxUNodes-1,2:iNyUNodes-1)=rDY/ap(2:iNxUNodes-1,2:iNyUNodes-1)
! solve the tridiagonal systems
!  for each i-direction line
      call tridag_i(as,aw,ap,ae,an,su,rU,iNxUNodes,iNyUNodes,iNxMax,iNyMax)
!  for each j-direction line
      call tridag_j(aw,as,ap,an,ae,su,rU,iNyUNodes,iNxUNodes,iNxMax,iNyMax)

      return 
end subroutine UMomentum



subroutine VMomentum(iNxVNodes, iNyVNodes,iNxMax,iNyMax,rDx,rDy,rU,rV,rP,rDV,rAlphaV,rVRes,rReynolds)
!arguments
    integer, intent(in) :: iNxVNodes, iNyVNodes,iNxMax,iNyMax
    real, intent(in)    :: rDx,rDy, rAlphaV,rReynolds
    real, intent(inout) :: rVRes
    real, dimension (iNxMax, iNyMax), intent (in) :: rU,rP
    real, dimension (iNxMax, iNyMax), intent (inout) :: rV,rDV
! local variables
    real ::  de, dw, dn, ds, ge, gw, gn, gs, rTemp
    real, dimension(iNxMax,iNyMax) ::aw,ae,as,an,ap,su
    integer :: i,j

! set boundary conditions for coefficients in the equations
! south
    aw(:,1)=0.
    ae(:,1)=0.
    as(:,1)=0.
    an(:,1)=0.
    ap(:,1)=1.
    su(:,1)=0.
! west
    aw(1,:)=0.
    ae(1,:)=0.
    as(1,:)=0.
    an(1,:)=0.
    ap(1,:)=1.
    su(1,:)=0.
! north
    aw(:,iNyVNodes)=0.
    ae(:,iNyVNodes)=0.
    as(:,iNyVNodes)=0.
    an(:,iNyVNodes)=0.
    ap(:,iNyVNodes)=1.
    su(:,iNyVNodes)=0.
! east
    aw(iNxVNodes,:)=0.
    ae(iNxVNodes,:)=0.
    as(iNxVNodes,:)=0.
    an(iNxVNodes,:)=0.
    ap(iNxVNodes,:)=1.
    su(iNxVNodes,:)=0.
! set the coefficeints in the interior staggered cells

      do i = 2 , iNxVNodes-1
         do j = 2 , iNyVNodes-1
            if ( i .eq. 2 ) then
             dw = 2.*rDY/(rDX*rReynolds)
             gw = rDY*(rU(i-1, j) + rU(i-1, j+1)) / 2.
            else
             dw = rDY/(rDX*rReynolds)
             gw = rDY*(rU(i-1, j) + rU(i-1, j+1)) / 2.
            end if
            if ( i .eq. (iNxVNodes - 1)) then   
             de = 2.*rDY/(rDX*rReynolds)
             ge = rDY*(rU(i, j) + rU(i, j+1)) / 2.
	    else
             de = rDY/(rDX*rReynolds)
             ge = rDY*(rU(i, j) + rU(i, j+1)) / 2.
	    end if
            ds = rDX/(rDY*rReynolds)
            gs = rDX*(rV(i, j-1) + rV(i, j)) / 2.
            dn = rDX/(rDY*rReynolds)
            gn = rDX*(rV(i, j) + rV(i, j+1)) / 2.

            aw(i, j) = dw + max( 0.0 ,  gw)
            ae(i, j) = de + max( 0.0 , -ge)
            as(i, j) = ds + max( 0.0 ,  gs)
            an(i, j) = dn + max( 0.0 , -gn)
            ap(i, j) = dw+de+ds+dn+max( 0.0 ,  -gw) + max( 0.0 , ge) +max( 0.0 ,  -gs)+max( 0.0 ,  gn)
            su(i, j) = -(rP(i, j+1) - rP(i, J))*rDX
         end do
      end do
! evaluate the residual in the interior
      rVRes=0.
      do i = 2 , iNxVNodes - 1
         do j = 2 , iNyVNodes - 1 
           rTemp=abs( ap(i, j)*rV(i, j) & 
                    - ae(i, j)*rV(i+1, j) - aw(i, j)*rV(i-1, j) &
                    - an(i, j)*rV(i, j+1) - as(i, j)*rV(i, j-1) &
                    - su(i, j))
           if(rTemp.ge.rVRes) rVRes=rTemp
         end do
      end do

!output tecplot file:
! apply underrelaxation to interior cells
      ap(2:iNxVNodes-1,2:iNyVNodes-1)=ap(2:iNxVNodes-1,2:iNyVNodes-1)/rAlphaV
      su(2:iNxVNodes-1,2:iNyVNodes-1)=su(2:iNxVNodes-1,2:iNyVNodes-1)+ &
        (1.-rAlphaV)*ap(2:iNxVNodes-1,2:iNyVNodes-1)*rV(2:iNxVNodes-1,2:iNyVNodes-1)
! store dV for pressure correction:
      rDV(2:iNxVNodes-1,2:iNyVNodes-1)=rDX/ap(2:iNxVNodes-1,2:iNyVNodes-1)
! solve the tridiagonal systems
!  for each i-direction line
      call tridag_i(as,aw,ap,ae,an,su,rV,iNxVNodes,iNyVNodes,iNxMax,iNyMax)
!  for each j-direction line
      call tridag_j(aw,as,ap,an,ae,su,rV,iNyVNodes,iNxVNodes,iNxMax,iNyMax)
      return 
end subroutine VMomentum

subroutine PCorrection(iNxPNodes,iNyPNodes,iNxMax,iNyMax,rDx,rDy,rU,rV,rDU,rDV,rDP,rPRes)
!arguments
    integer, intent(in) :: iNxPNodes, iNyPNodes,iNxMax,iNyMax
    real, intent(in)    :: rDx,rDy
    real, intent(inout) :: rPRes
    real, dimension (iNxMax, iNyMax), intent (in) :: rU,rV,rDU,rDV
    real, dimension (iNxMax, iNyMax), intent (inout) :: rDP
! local variables
    real, dimension(iNxMax,iNyMax) ::aw,ae,as,an,ap,su
    integer :: i,j
    real :: rTemp

! initalise pressure correction
    rDP=0.
! set boundary conditions for coefficients in the equations
    aw=0.
    ae=0.
    as=0.
    an=0.
    ap=0.
    su=0.
! south
    aw(2:iNxPNodes-1,1)=0.
    ae(2:iNxPNodes-1,1)=0.
    as(2:iNxPNodes-1,1)=0.
    an(2:iNxPNodes-1,1)=0.
    ap(2:iNxPNodes-1,1)=1.
    su(2:iNxPNodes-1,1)=0.
! west
    aw(1,2:iNyPNodes-1)=0.
    ae(1,2:iNyPNodes-1)=0.
    as(1,2:iNyPNodes-1)=0.
    an(1,2:iNyPNodes-1)=0.
    ap(1,2:iNyPNodes-1)=1.
    su(1,2:iNyPNodes-1)=0.
! north
    aw(2:iNxPNodes-1,iNyPNodes)=0.
    ae(2:iNxPNodes-1,iNyPNodes)=0.
    as(2:iNxPNodes-1,iNyPNodes)=0.
    an(2:iNxPNodes-1,iNyPNodes)=0.
    ap(2:iNxPNodes-1,iNyPNodes)=1.
    su(2:iNxPNodes-1,iNyPNodes)=0.
! east
    aw(iNxPNodes,2:iNyPNodes-1)=0.
    ae(iNxPNodes,2:iNyPNodes-1)=0.
    as(iNxPNodes,2:iNyPNodes-1)=0.
    an(iNxPNodes,2:iNyPNodes-1)=0.
    ap(iNxPNodes,2:iNyPNodes-1)=1.
    su(iNxPNodes,2:iNyPNodes-1)=0.
! set the coefficeints in the interior staggered cells

      do j = 2 , iNyPNodes-1
       do i = 2 , iNxPNodes-1
           aw(i, j) = rDY*rDU(i-1, j)
           ae(i, j) = rDY*rDU(i, j)
           as(i, j) = rDX*rDV(i, j-1)
           an(i, j) = rDX*rDV(i, j)            
           ap(i,j)  = aw(i, j)+ae(i, j)+as(i, j)+an(i, j)
           su(i, j) = rDY*(rU(i-1, j)-rU(i, j))+ &
                      rDX*(rV(i, j-1)-rV(i, j))
         end do
      end do
! sweep through the domain repeatedly and solve the algebraic equations in each line
      do i=1,10
        call tridag_i(as,aw,ap,ae,an,su,rDP,iNxPNodes,iNyPNodes,iNxMax,iNyMax)
        call tridag_j(aw,as,ap,an,ae,su,rDP,iNyPNodes,iNxPNodes,iNxMax,iNyMax)
      end do
! compute the pressure residual 
      rPRes=0.
      do j=2,iNyPNodes-1
       do i=2,iNxPNodes-1
        rTemp=abs(su(i,j))
        if(rTemp.ge.rPRes) rPRes=rTemp

       end do
      end do
      return 
end subroutine PCorrection


subroutine Scalar(iNxPNodes, iNyPNodes,iNxMax,iNyMax,rDx,rDy,rU,rV,rC,rAlphaC,rCRes,rPeclet,rInjStart,rInjEnd)
!arguments
    integer, intent(in) :: iNxPNodes, iNyPNodes,iNxMax,iNyMax
    real, intent(in)    :: rDx,rDy, rAlphaC,rPeclet,rInjStart,rInjEnd
    real, intent(inout) :: rCRes
    real, dimension (iNxMax, iNyMax), intent (in) :: rU,rV
    real, dimension (iNxMax, iNyMax), intent (inout) :: rC
! local variables
    real ::  de, dw, dn, ds, ge, gw, gn, gs, rTemp
    real, dimension(iNxMax,iNyMax) ::aw,ae,as,an,ap,su
    integer :: i,j

! set boundary conditions for coefficients in the equations
! south
    aw(:,1)=0.
    ae(:,1)=0.
    as(:,1)=0.
    an(:,1)=1.
    ap(:,1)=1.
    su(:,1)=0.
! west
    aw(1,:)=0.
    ae(1,:)=0.
    as(1,:)=0.
    an(1,:)=0.
    do j=2,iNyPNodes-1
     rTemp=(j-2)*rDY+0.5*rDY
     if (rTemp.ge.rInjStart.and.rTemp.le.rInjEnd) then
      ap(1,j)=1.
      su(1,j)=1.
     else
      ap(1,j)=1.
      su(1,j)=0.
     end if
    end do
! north
    aw(:,iNyPNodes)=0.
    ae(:,iNyPNodes)=0.
    as(:,iNyPNodes)=1.
    an(:,iNyPNodes)=0.
    ap(:,iNyPNodes)=1.
    su(:,iNyPNodes)=0.
! east
    aw(iNxPNodes,:)=0.
    ae(iNxPNodes,:)=0.
    as(iNxPNodes,:)=0.
    an(iNxPNodes,:)=0.
    ap(iNxPNodes,:)=1.
    su(iNxPNodes,:)=0.
! set the coefficeints in the interior staggered cells

      do i = 2 , iNxPNodes-1
         do j = 2 , iNyPNodes-1
          dw = rDY/(rDX*rPeclet)
          gw = rDY*rU(i-1, j)
          de = rDY/(rDX*rPeclet)
          ge = rDY*rU(i, j)
          if (j.eq.2) then
            ds = 0.0
            gs = 0.0
          else
            ds = rDX/(rDY*rPeclet)
            gs = rDX*rV(i, j-1)
          end if 
          if (j.eq.(iNyPNodes-1)) then   
           dn = 0.0
           gn = 0.0
          else   
           dn = rDX/(rDY*rPeclet)
           gn = rDX*rV(i, j)
          end if
           aw(i, j) = dw + gw/2
	   ae(i, j) = de - ge/2
	   as(i, j) = ds + gs/2
	   an(i, j) = dn - gn/2
	   ap(i, j) =  dw+de+ds+dn+-gw/2 + ge/2 + -gs/2 + gn/2
	   su(i, j) = 0.0
         end do
      end do
! evaluate the residual in the interior
      rCRes=0.
      do i = 2 , iNxPNodes - 1
         do j = 2 , iNyPNodes - 1 
           rTemp=abs( ap(i, j)*rC(i, j) & 
                    - ae(i, j)*rC(i+1, j) - aw(i, j)*rC(i-1, j) &
                    - an(i, j)*rC(i, j+1) - as(i, j)*rC(i, j-1) &
                    - su(i, j))
           if(rTemp.ge.rCRes) rCRes=rTemp
         end do
      end do

!output tecplot file:
! apply underrelaxation to interior cells
      ap(2:iNxPNodes-1,2:iNyPNodes-1)=ap(2:iNxPNodes-1,2:iNyPNodes-1)/rAlphaC
      su(2:iNxPNodes-1,2:iNyPNodes-1)=su(2:iNxPNodes-1,2:iNyPNodes-1)+ &
        (1.-rAlphaC)*ap(2:iNxPNodes-1,2:iNyPNodes-1)*rC(2:iNxPNodes-1,2:iNyPNodes-1)
! solve the tridiagonal systems
!  for each i-direction line
      call tridag_i(as,aw,ap,ae,an,su,rC,iNxPNodes,iNyPNodes,iNxMax,iNyMax)
!  for each j-direction linE
      call tridag_j(aw,as,ap,an,ae,su,rC,iNyPNodes,iNxPNodes,iNxMax,iNyMax)
      return 
end subroutine Scalar

!-----------------------------------------------------------------------
!     Tri-diagonal matrix solver (based on Numerical Recipes)  
!
!     Solving i-lines, looping over domain in j-direction
!
!     Note signs of Coefficients 
! note - r - righthand side 
!-----------------------------------------------------------------------
      subroutine tridag_i(d, a, b, c, e, r, u, n, l,iNxMax,iNyMax)
!-----------------------------------------------------------------------
!     Argument type declarations
!-----------------------------------------------------------------------
      integer, intent(in) ::  n, l, iNxMax, iNyMax
      real, dimension(iNxMax,iNyMax) ::     a, b, c,d,e,r,u
!-----------------------------------------------------------------------
!     Local variable type declarations
!-----------------------------------------------------------------------
      integer j, m
      real    bet,gam(n),rhs
!-----------------------------------------------------------------------
!     Loop in j-direction
!-----------------------------------------------------------------------
      do j = 2 , l-1
!-----------------------------------------------------------------------
!     Forward recurrence
!-----------------------------------------------------------------------
         bet  = b(1, j)
         u(1, j) = r(1, j)/bet
         do m = 2 , n
            gam(m) = -c(m-1, j) / bet
            bet    = b(m, j) + a(m, j)*gam(m)
            rhs    = r(m,j) + d(m,j)*u(m,j-1) + e(m,j)*u(m,j+1)
            u(m, j)= (rhs + a(m, j)*u(m-1, j))/bet
         end do
!-----------------------------------------------------------------------
!     Backward recurrence
!-----------------------------------------------------------------------
         do m = n-1, 1, -1
            u(m, j) = u(m, j) - gam(m+1)*u(m+1, j)
         end do
!-----------------------------------------------------------------------
!     Next j line
!-----------------------------------------------------------------------
      end do
!-----------------------------------------------------------------------
!     finished so return
!-----------------------------------------------------------------------
      return
      end



!-----------------------------------------------------------------------
!     Tri-diagonal matrix solver (based on Numerical Recipes)  
!
!     Solving j-lines, looping over domain in i-direction
!
!     Note signs of Coefficients 
! note - r - righthand side 
!-----------------------------------------------------------------------
      subroutine tridag_j(d, a, b, c, e, r, u, n, l,iNxMax,iNyMax)
!-----------------------------------------------------------------------
!     Argument type declarations
!-----------------------------------------------------------------------
      integer, intent(in) ::  l, n,iNxMax,iNyMax
      real, dimension (iNxMax,iNyMax) ::  a,b,c,d,e,r,u

!-----------------------------------------------------------------------
!     Local variable type declarations
!-----------------------------------------------------------------------
      integer i, m
      real    bet,gam(n),rhs
!-----------------------------------------------------------------------
!     Loop in i-diection
!-----------------------------------------------------------------------
      do i = 2, l-1
!
!-----------------------------------------------------------------------
!     Forward recurrence
!-----------------------------------------------------------------------
         bet  = b(i, 1)
         u(i, 1) = r(i, 1)/bet
         do m = 2, n
            gam(m) = -c(i, m-1) / bet
            bet    = b(i, m) + a(i, m)*gam(m)
            rhs    = r(i,m) + d(i,m)*u(i-1,m) + e(i,m)*u(i+1,m)
            u(i, m)= (rhs + a(i, m)*u(i, m-1))/bet
         end do
!-----------------------------------------------------------------------
!     Backward recurrence
!-----------------------------------------------------------------------
         do m = n-1, 1, -1
            u(i, m) = u(i, m) - gam(m+1)*u(i, m+1)
         end do
!-----------------------------------------------------------------------
!     Next i line
!-----------------------------------------------------------------------
      end do
!-----------------------------------------------------------------------
!     finished so return
!-----------------------------------------------------------------------
      return
      end
