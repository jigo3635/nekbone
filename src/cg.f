c-----------------------------------------------------------------------
      subroutine cg(x,f,g,c,r,w,p,z,n,niter,flop_cg)


c     Solve Ax=f where A is SPD and is invoked by ax()
c
c     Output:  x - vector of length n
c
c     Input:   f - vector of length n
c     Input:   g - geometric factors for SEM operator
c     Input:   c - inverse of the counting matrix
c
c     Work arrays:   r,w,p,z  - vectors of length n
c
c     User-provided ax(w,z,n) returns  w := Az,  
c
c     User-provided solveM(z,r,n) ) returns  z := M^-1 r,  
c

#ifdef XSMM
      USE :: LIBXSMM
      USE :: STREAM_UPDATE_KERNELS
#endif

      include 'SIZE'
      include 'DXYZ'
C      include 'INPUT'
      common /mymask/cmask(-1:lx1*ly1*lz1*lelt)
      parameter (lt=lx1*ly1*lz1*lelt)


      real x(n),f(n),r(n),z(n),c(n)

      character*1 ans

#ifdef XSMM
      real w(lx1,ly1,lz1,lelt)
      real p(lx1,ly1,lz1,lelt)
      real wk(lx1,ly1,lz1,lelt)
      real g(6,lx1,ly1,lz1,lelt)

      INTEGER, PARAMETER :: T = KIND(0D0)
      REAL, PARAMETER :: alpha = 1, beta0 = 0, beta1 = 1
      
      REAL, allocatable, dimension(:,:,:), target :: ur, us, ut
      REAL, allocatable, target :: dx(:,:), dxt(:,:)
      REAL, ALLOCATABLE,TARGET,SAVE :: tm1(:,:,:), tm2(:,:,:),
     $     tm3(:,:,:)

      TYPE(LIBXSMM_DMMFUNCTION) :: xmm1, xmm2, xmm3, xmm4, xmm5

      DOUBLE PRECISION :: max_diff
      INTEGER :: argc, m, n, k, routine, check
      INTEGER(8) :: i, j, ix, iy, iz,  s, size0, size1, size, 
     $     start, it
      CHARACTER(32) :: argv
      s = lelt
      size = s
      ALLOCATE(ur(lx1,ly1,lz1))
      ALLOCATE(us(lx1,ly1,lz1))
      ALLOCATE(ut(lx1,ly1,lz1))
      ALLOCATE(dx(lx1,lx1), dxt(ly1,ly1))


      lxyz = lx1*ly1*lz1
! Initialize LIBXSMM
      CALL libxsmm_init()

!  Initialize 
      do j = 1,ly1
         do i = 1, lx1
            dx(i,j)  = dxm1(i,j)
            dxt(i,j) = dxtm1(i,j)
         enddo
      enddo

!      WRITE(*, "(A)") " Streamed... (specialized)"
      CALL libxsmm_dispatch(xmm1,lx1,ly1*lz1,lx1,alpha=alpha,beta=beta0)
      CALL libxsmm_dispatch(xmm2,lx1,ly1,ly1,alpha=alpha,beta=beta0)
      CALL libxsmm_dispatch(xmm3,lx1*ly1,lz1,lz1,alpha=alpha,beta=beta0)
      CALL libxsmm_dispatch(xmm4,lx1,ly1,ly1,alpha=alpha,beta=beta1)
      CALL libxsmm_dispatch(xmm5,lx1*ly1,lz1,lz1,alpha=alpha,beta=beta1)

      IF (libxsmm_available(xmm1).AND.libxsmm_available(xmm2) 
     $     .AND. libxsmm_available(xmm3) .AND. libxsmm_available(xmm4) 
     $     .AND. libxsmm_available(xmm5) ) THEN

         ALLOCATE(tm1(lx1,ly1,lz1), tm2(lx1,ly1,lz1), tm3(lx1,ly1,lz1))
         tm1 = 0; tm2 = 0; tm3 = 0
         start = libxsmm_timer_tick()
      ELSE
         WRITE(*,*) "Could not build specialized function(s)!"
      END IF

#else
      real w(n),p(n)
      real wk(lt)
      real ur(lt),us(lt),ut(lt)
      real g
#endif

      pap = 0.0

c     set machine tolerances
      one = 1.
      eps = 1.e-20
      if (one+eps .eq. one) eps = 1.e-14
      if (one+eps .eq. one) eps = 1.e-7

      rtz1=1.0

      call rzero(x,n)
      call copy (r,f,n)
      call maskit (r,cmask,nx1,ny1,nz1) ! Zero out Dirichlet conditions

      rnorm = sqrt(glsc3(r,c,r,n))
      iter = 0
      if (nid.eq.0)  write(6,6) iter,rnorm

      miter = niter
c     call tester(z,r,n)  
      do iter=1,miter
         call solveM(z,r,n)    ! preconditioner here

         rtz2=rtz1                                                       ! OPS
         rtz1=glsc3(r,c,z,n)   ! parallel weighted inner product r^T C z ! 3n

         beta_cg = rtz1/rtz2
         if (iter.eq.1) beta_cg=0.0
         call add2s1(p,z,beta_cg,n)                                         ! 2n

#ifdef XSMM
C         call ax_xsmm(w,p,g,ur,us,ut,wk,n,
C     $                dx1,dxt,xmm1,xmm2,xmm3,tm1,tm2,tm3 ) ! flopa


C Local_grad3
      DO e = 1, nelt
         CALL libxsmm_call(xmm1, C_LOC(dx), C_LOC(p(1,1,1,e)),
     $        C_LOC(ur(1,1,1)))
         DO j = 1, ly1
            CALL libxsmm_call(xmm2, C_LOC(p(1,1,j,e)), 
     $           C_LOC(dxt), C_LOC(us(1,1,j)))
         END DO
         CALL libxsmm_call(xmm3, C_LOC(p(1,1,1,e)), C_LOC(dxt),
     $        C_LOC(ut(1,1,1)))

C Geometric multiplication

         DO k = 1,lz1
         DO j = 1,ly1
         DO i = 1,lx1
            wr = g(1,i,j,k,e)*ur(i,j,k) + 
     $           g(2,i,j,k,e)*us(i,j,k) +
     $           g(3,i,j,k,e)*ut(i,j,k)
            ws = g(2,i,j,k,e)*ur(i,j,k) + 
     $           g(4,i,j,k,e)*us(i,j,k) +
     $           g(5,i,j,k,e)*ut(i,j,k)
            wt = g(3,i,j,k,e)*ur(i,j,k) +
     $           g(5,i,j,k,e)*us(i,j,k) +
     $           g(6,i,j,k,e)*ut(i,j,k)
            ur(i,j,k) = wr
            us(i,j,k) = ws
            ut(i,j,k) = wt
         ENDDO
         ENDDO
         ENDDO            
                  
C local_grad3_t

         CALL libxsmm_call(xmm1,  C_LOC(dxt), C_LOC(ur(1,1,1)),
     $        C_LOC(tm1(1,1,1)))

         DO j = 1, ly1
            CALL libxsmm_call(xmm4, C_LOC(us(1,1,j)), 
     $           C_LOC(dx), C_LOC(tm1(1,1,j)))
         END DO

         CALL libxsmm_call(xmm5, C_LOC(ut(1,1,1)), C_LOC(dx),
     $        C_LOC(tm1(1,1,1)))

         CALL stream_vector_copy(tm1(1,1,1),w(1,1,1,e),lxyz)

      END DO

      call dssum(w)         ! Gather-scatter operation  ! w   = QQ  w
                                                           !            L
      call add2s2(w,p,.1,n)   !2n
      call maskit(w,cmask,nx1,ny1,nz1)  ! Zero out Dirichlet conditions

#else
         call ax(w,p,g,ur,us,ut,wk,n)                                    ! flopa
#endif
         pap=glsc3(w,c,p,n)                                              ! 3n

         alpha_cg=rtz1/pap
         alphm=-alpha_cg
         call add2s2(x,p,alpha_cg,n)                                        ! 2n
         call add2s2(r,w,alphm,n)                                        ! 2n

         rtr = glsc3(r,c,r,n)                                            ! 3n
         if (iter.eq.1) rlim2 = rtr*eps**2
         if (iter.eq.1) rtr0  = rtr
         rnorm = sqrt(rtr)
c        if (nid.eq.0.and.mod(iter,100).eq.0) 
c    $      write(6,6) iter,rnorm,alpha_cg,beta_cg,pap
    6    format('cg:',i4,1p4e12.4)
c        if (rtr.le.rlim2) goto 1001

      enddo

 1001 continue

#ifdef XSMM
      duration = libxsmm_timer_duration(start, libxsmm_timer_tick())
      
      DEALLOCATE(tm1, tm2, tm3)
      DEALLOCATE(ur, us, ut)

! finalize LIBXSMM
      CALL libxsmm_finalize()
#endif

      if (nid.eq.0) write(6,6) iter,rnorm,alpha_cg,beta_cg,pap
      flop_cg = flop_cg + niter*15.*n

      return
      end
c-----------------------------------------------------------------------
      subroutine solveM(z,r,n)
      include 'INPUT'
      real z(n),r(n)

      nn = n
      call h1mg_solve(z,r,nn)

      return
      end
c-----------------------------------------------------------------------
      subroutine ax(w,u,gxyz,ur,us,ut,wk,n) ! Matrix-vector product: w=A*u

      include 'SIZE'
      include 'TOTAL'

      real w(nx1*ny1*nz1,nelt),u(nx1*ny1*nz1,nelt)
      real gxyz(2*ldim,nx1*ny1*nz1,nelt)

      parameter (lt=lx1*ly1*lz1*lelt)
      real ur(lt),us(lt),ut(lt),wk(lt)
      common /mymask/cmask(-1:lx1*ly1*lz1*lelt)

      integer e


      do e=1,nelt                                ! ~
         call ax_e( w(1,e),u(1,e),gxyz(1,1,e)    ! w   = A  u
     $                             ,ur,us,ut,wk) !  L     L  L
      enddo                                      ! 

      call dssum(w)         ! Gather-scatter operation  ! w   = QQ  w
                                                           !            L
      call add2s2(w,u,.1,n)   !2n
      call maskit(w,cmask,nx1,ny1,nz1)  ! Zero out Dirichlet conditions

      nxyz=nx1*ny1*nz1
      flop_a = flop_a + (19*nxyz+12*nx1*nxyz)*nelt

      return
      end
c-------------------------------------------------------------------------
      subroutine ax1(w,u,n)
      include 'SIZE'
      real w(n),u(n)
      real h2i
  
      h2i = (n+1)*(n+1)  
      do i = 2,n-1
         w(i)=h2i*(2*u(i)-u(i-1)-u(i+1))
      enddo
      w(1)  = h2i*(2*u(1)-u(2  ))
      w(n)  = h2i*(2*u(n)-u(n-1))

      return
      end
c-------------------------------------------------------------------------
      subroutine ax_e(w,u,g,ur,us,ut,wk) ! Local matrix-vector product
      include 'SIZE'
      include 'TOTAL'

      parameter (lxyz=lx1*ly1*lz1)
      real ur(lxyz),us(lxyz),ut(lxyz),wk(lxyz)
      real w(nx1*ny1*nz1),u(nx1*ny1*nz1),g(2*ldim,nx1*ny1*nz1)


      nxyz = nx1*ny1*nz1
      n    = nx1-1

      call local_grad3(ur,us,ut,u,n,dxm1,dxtm1)

      do i=1,nxyz
         wr = g(1,i)*ur(i) + g(2,i)*us(i) + g(3,i)*ut(i)
         ws = g(2,i)*ur(i) + g(4,i)*us(i) + g(5,i)*ut(i)
         wt = g(3,i)*ur(i) + g(5,i)*us(i) + g(6,i)*ut(i)
         ur(i) = wr
         us(i) = ws
         ut(i) = wt
      enddo

      call local_grad3_t(w,ur,us,ut,n,dxm1,dxtm1,wk)

      return
      end
c-------------------------------------------------------------------------
      subroutine local_grad3(ur,us,ut,u,n,D,Dt)
c     Output: ur,us,ut         Input:u,n,D,Dt
      real ur(0:n,0:n,0:n),us(0:n,0:n,0:n),ut(0:n,0:n,0:n)
      real u (0:n,0:n,0:n)
      real D (0:n,0:n),Dt(0:n,0:n)
      integer e

      m1 = n+1
      m2 = m1*m1

      call mxm(D ,m1,u,m1,ur,m2)
      do k=0,n
         call mxm(u(0,0,k),m1,Dt,m1,us(0,0,k),m1)
      enddo
      call mxm(u,m2,Dt,m1,ut,m1)

      return
      end
c-----------------------------------------------------------------------
      subroutine local_grad3_t(u,ur,us,ut,N,D,Dt,w)
c     Output: ur,us,ut         Input:u,N,D,Dt
      real u (0:N,0:N,0:N)
      real ur(0:N,0:N,0:N),us(0:N,0:N,0:N),ut(0:N,0:N,0:N)
      real D (0:N,0:N),Dt(0:N,0:N)
      real w (0:N,0:N,0:N)
      integer e

      m1 = N+1
      m2 = m1*m1
      m3 = m1*m1*m1

      call mxm(Dt,m1,ur,m1,u,m2)

      do k=0,N
         call mxm(us(0,0,k),m1,D ,m1,w(0,0,k),m1)
      enddo
      call add2(u,w,m3)

      call mxm(ut,m2,D ,m1,w,m1)
      call add2(u,w,m3)

      return
      end
c-----------------------------------------------------------------------
      subroutine maskit(w,pmask,nx,ny,nz)   ! Zero out Dirichlet conditions
      include 'SIZE'
      include 'PARALLEL'

      real pmask(-1:lx1*ly1*lz1*lelt)
      real w(1)
      integer e

      nxyz = nx*ny*nz
      nxy  = nx*ny
      if(pmask(-1).lt.0) then
        j=pmask(0)
        do i = 1,j
           k = pmask(i)
           w(k)=0.0
        enddo
      else
c         Zero out Dirichlet boundaries.
c
c                      +------+     ^ Y
c                     /   3  /|     |
c               4--> /      / |     |
c                   +------+ 2 +    +----> X
c                   |   5  |  /    /
c                   |      | /    /
c                   +------+     Z   
c

        nn = 0
        do e  = 1,nelt
          call get_face(w,nx,e)
          do i = 1,nxyz
             if(w(i).eq.0) then
               nn=nn+1
               pmask(nn)=i
             endif
          enddo
        enddo     
        pmask(-1) = -1.
        pmask(0) = nn
      endif


      return
      end
c-----------------------------------------------------------------------
      subroutine masko(w)   ! Old 'mask'
      include 'SIZE'
      real w(1)

      if (nid.eq.0) w(1) = 0.  ! suitable for solvability

      return
      end
c-----------------------------------------------------------------------
      subroutine masking(w,nx,e,x0,x1,y0,y1,z0,z1)
c     Zeros out boundary
      include 'SIZE'
      integer e,x0,x1,y0,y1,z0,z1
      real w(nx,nx,nx,nelt)
      
c       write(6,*) x0,x1,y0,y1,z0,z1
      do k=z0,z1
      do j=y0,y1
      do i=x0,x1
          w(i,j,k,e)=0.0
      enddo
      enddo
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine tester(z,r,n)
c     Used to test if solution to precond. is SPD
      real r(n),z(n)

      do j=1,n
         call rzero(r,n)
         r(j) = 1.0
         call solveM(z,r,n)
         do i=1,n
            write(79,*) z(i)
         enddo
      enddo
      call exitt0
      return
      end
c-----------------------------------------------------------------------
      subroutine get_face(w,nx,ie)
c     zero out all boundaries as Dirichlet
c     to change, change this routine to only zero out 
c     the nodes desired to be Dirichlet, and leave others alone.
      include 'SIZE'
      include 'PARALLEL'
      real w(1)
      integer nx,ie,nelx,nely,nelz
      integer x0,x1,y0,y1,z0,z1
      
      x0=1
      y0=1
      z0=1
      x1=nx
      y1=nx
      z1=nx
      
      nelxy=nelx*nely
      ngl = lglel(ie)        !global element number

      ir = 1+(ngl-1)/nelxy   !global z-count
      iq = mod1(ngl,nelxy)   !global y-count
      iq = 1+(iq-1)/nelx     
      ip = mod1(ngl,nelx)    !global x-count

c     write(6,1) ip,iq,ir,nelx,nely,nelz, nelt,' test it'
c  1  format(7i7,a8)

      if(mod(ip,nelx).eq.1.or.nelx.eq.1)   then  ! Face4
         x0=1
         x1=1
         call masking(w,nx,ie,x0,x1,y0,y1,z0,z1)
      endif
      if(mod(ip,nelx).eq.0)               then   ! Face2
         x0=nx
         x1=nx
         call masking(w,nx,ie,x0,x1,y0,y1,z0,z1)
      endif

      x0=1
      x1=nx
      if(mod(iq,nely).eq.1.or.nely.eq.1) then    ! Face1
         y0=1
         y1=1
         call masking(w,nx,ie,x0,x1,y0,y1,z0,z1)
      endif
      if(mod(iq,nely).eq.0)              then    ! Face3
         y0=nx
         y1=nx
         call masking(w,nx,ie,x0,x1,y0,y1,z0,z1)
      endif

      y0=1
      y1=nx
      if(mod(ir,nelz).eq.1.or.nelz.eq.1) then    ! Face5
         z0=1
         z1=1
         call masking(w,nx,ie,x0,x1,y0,y1,z0,z1)
      endif
      if(mod(ir,nelz).eq.0)              then    ! Face6
         z1=nx
         z0=nx
         call masking(w,nx,ie,x0,x1,y0,y1,z0,z1)
      endif

      return
      end
c-----------------------------------------------------------------------


#if 0
#ifdef XSMM
c-----------------------------------------------------------------------
      subroutine ax_xsmm(w,u,gxyz,ur,us,ut,wk,n,
     $                   dx,dxt,xmm1,xmm2,xmm3,tm1,tm2,tm3) ! Matrix-vector product: w=A*u

      include 'SIZE'
      include 'TOTAL'

!      real w(nx1*ny1*nz1,nelt),u(nx1*ny1*nz1,nelt)
      REAL, allocatable, target :: dx(:,:), dxt(:,:)

c      TYPE(LIBXSMM_DMMFUNCTION) :: xmm1, xmm2, xmm3

      REAL, allocatable, dimension(:,:,:,:), target :: ur, us, ut
      REAL, allocatable, dimension(:,:,:), target :: tm1, tm2, tm3
      
      real gxyz(2*ldim,lx1,ly1,lz1,lelt)

      parameter (lt=lx1*ly1*lz1*lelt)
      real wk(lx1,ly1,lz1,lelt)
      common /mymask/cmask(-1:lx1*ly1*lz1*lelt)

      integer i,j,k, e, lxyz
      lxyz = lx1*ly1*lz1

C Local_grad3
      DO e = 1, nelt
         CALL libxsmm_call(xmm1,  C_LOC(dx), C_LOC(u(1,1,1,e)),
     $        C_LOC(tm1(1,1,1)))
         CALL stream_vector_copy(tm1(1,1,1), ur(1,1,1,e),lxyz)
         DO j = 1, ly1
            CALL libxsmm_call(xmm2, C_LOC(u(1,1,j,e)), 
     $           C_LOC(dxt), C_LOC(tm2(1,1,j)))
         END DO
         CALL stream_vector_copy(tm2(1,1,1),us(1,1,1,e),lxyz)
         CALL libxsmm_call(xmm3, C_LOC(u(1,1,1,e)), C_LOC(dxt),
     $        C_LOC(tm3(1,1,1)))
         CALL stream_vector_copy(tm3(1,1,1),ut(1,1,1,e),lxyz)
      ENDDO

C Geometric multiplication
      DO e = 1,nelt
         DO k = 1,lz1
         DO j = 1,ly1
         DO i = 1,lx1
            wr = gxyz(1,i,j,k,e)*ur(i,j,k,e) + 
     $           gxyz(2,i,j,k,e)*us(i,j,k,e) +
     $           gxyz(3,i,j,k,e)*ut(i,j,k,e)
            ws = gxyz(2,i,j,k,e)*ur(i,j,k,e) + 
     $           gxyz(4,i,j,k,e)*us(i,j,k,e) +
     $           gxyz(5,i,j,k,e)*ut(i,j,k,e)
            wt = gxyz(3,i,j,k,e)*ur(i,j,k,e) +
     $           gxyz(5,i,j,k,e)*us(i,j,k,e) +
     $           gxyz(6,i,j,k,e)*ut(i,j,k,e)
            ur(i,j,k,e) = wr
            us(i,j,k,e) = ws
            ut(i,j,k,e) = wt
         ENDDO
         ENDDO
         ENDDO            
      ENDDO
                  
C local_grad3_t
      DO e = 1, nelt
         CALL libxsmm_call(xmm1,  C_LOC(dxt), C_LOC(ur(1,1,1,e)),
     $        C_LOC(tm1(1,1,1)))

         CALL stream_vector_copy(tm1(1,1,1),w(1,1,1,e),lxyz)

         DO j = 1, ly1
            CALL libxsmm_call(xmm2, C_LOC(us(1,1,j,e)), 
     $           C_LOC(dx), C_LOC(tm2(1,1,j)))
         END DO
         CALL stream_vector_copy(tm2(1,1,1),wk(1,1,1,e),lxyz)

         CALL add2(w(1,1,1,e), wk(1,1,1,e), lxyz)
         CALL libxsmm_call(xmm3, C_LOC(ut(1,1,1,e)), C_LOC(dx),
     $        C_LOC(tm3(1,1,1)))
         CALL stream_vector_copy(tm3(1,1,1),wk(1,1,1,e),lxyz)
         CALL add2(w(1,1,1,e), wk(1,1,1,e), lxyz)
      END DO

      return
      end


c-----------------------------------------------------------------------
      SUBROUTINE performance(duration, iters, m, n, k, size)
      DOUBLE PRECISION, INTENT(IN) :: duration
      INTEGER, INTENT(IN)    :: m, n, k
C      INTEGER(8), INTENT(IN) :: size
      integer size
      real T
      T = 8.0
      IF (0.LT.duration) THEN
         WRITE(*, 2) CHAR(9), "performance:", 
     $         (1D-9 * iters * size * m * n * k * (4*(m+n+k) - 4) / 
     $        duration),     " GFLOPS/s"
         WRITE(*, 2) CHAR(9), "bandwidth:  ", 
     $        (size*m*n*k*(2)*T*iters / (duration * ISHFT(1_8, 30)))
     $        , " GB/s"
      END IF
      WRITE(*, 2) CHAR(9), "duration:   ", (1D3 * duration), " ms"

 2    format(1A,A,F10.5,A)

      return 
      END

#endif

#endif
