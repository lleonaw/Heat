c-----------------------------------------------------------------------
c-----| Discontinuous Galerkin for fluid |-----
c-----------------------------------------------------------------------
c     - Combine with dgf3_stpf, dgf3_3dlm.  
c       . Done            -  Thu Oct  8 10:06:19 CDT 2015
c     - Add in Zhang, Shu limiter (2010). !!
c       . Coded 2D (?3D)  -  Wed Oct 14 18:37:47 CDT 2015
c       . Isen ok.        -  Wed Oct 14 19:53:01 CDT 2015
c     - Add in Mach number and CFL computation. !!
c       . Done, isen ok   -  Mon Oct 12 18:12:12 CDT 2015
c     - Do de-aliasing for surface integrals. !!
c       . Coded. Testing. -  Fri Oct 16 13:42:13 CDT 2015
c     - Diffusion operator, 
c       . Bassi, Rebay 1997 Coded 
c       .        Testing  -  Mon Oct 26 23:22:24 CDT 2015
c         Seems shear ok. -  Wed Oct 28 14:00:08 CDT 2015
c         eddy_uv does not converge for velocity,
c           density seems fine. 1e-9 
c         Do more tests 
c       . H&W page 288, Central, LDG, IP 
c       . LDG to be coded -  Tue Oct 27 20:14:10 CDT 2015
c         Seems to be ok? -  Wed Nov  4 18:37:15 CST 2015
c       . IP  to be coded -  Wed Nov  4 18:37:44 CST 2015
c           IP coded in shear case, did not do well      
c       . SIPG to do      -  Thu Mar 31 12:25:57 CDT 2016
c           For shear case I guess 
c     - Heat problem, 2D  -  Wed Feb 17 09:55:52 CST 2016
c       . Coded now! Testing 
c     - Heat, primal, 2D  -  Tue Mar 22 22:09:56 CDT 2016
c       . Done, convergence testing Done
c     - Redo bc, v and f  -  Tue Feb 16 13:06:35 CST 2016
c       . Make v the previous f for flow dirichlet 
c       . Use f as the flux condition for temp dirichlet 
c       . Starting from this .f file
c     . Filter. !
c     . Mask of boundary. !
c     . Roe flux. !
c-----------------------------------------------------------------------
      subroutine dg_flow
c      /----------------------------------\
c-----|  dg method for fluid flow problem  |-----
c      \----------------------------------/
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'

c     miu = 0.01 
c     miu = 1.0e-20  ! 0, inviscid, Euler equations
      miu = param(2) 
      prt = 2.0e20  ! large prandtl number => neglib. thermal 
                    ! or just comment out lines in vis_flx 
      gama = 1.4

      call dg_advect_setup

      call dg_advect_execute

      return
      end
c-----------------------------------------------------------------------
c-----| Setting up communication |-----
c-----------------------------------------------------------------------
      subroutine izero8(a,n)
      integer*8 a(1)
      do i=1,n
         a(i)=0
      enddo
      return
      end
c-----------------------------------------------------------------------
      subroutine iface_vert_int8(fa,va,jz0,jz1,nel)
      include 'SIZE'
      integer*8 fa(nx1*nz1,2*ndim,nel),va(0:nx1+1,0:ny1+1,jz0:jz1,nel)
      integer e,f

      n = nx1*nz1*2*ndim*nel
      call izero8(fa,n)

      mx1 = nx1+2
      my1 = ny1+2
      mz1 = nz1+2
      if (ndim.eq.2) mz1=1

      nface = 2*ndim
      do e=1,nel
      do f=1,nface
         call facind (kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,f)

         if     (f.eq.1) then ! EB notation
            ky1=ky1-1
            ky2=ky1
         elseif (f.eq.2) then
            kx1=kx1+1
            kx2=kx1
         elseif (f.eq.3) then
            ky1=ky1+1
            ky2=ky1
         elseif (f.eq.4) then
            kx1=kx1-1
            kx2=kx1
         elseif (f.eq.5) then
            kz1=kz1-1
            kz2=kz1
         elseif (f.eq.6) then
            kz1=kz1+1
            kz2=kz1
         endif

         i = 0
         do iz=kz1,kz2
         do iy=ky1,ky2
         do ix=kx1,kx2
            i = i+1
            fa(i,f,e)=va(ix,iy,iz,e)
c           write(6,*) 'fa:',fa(i,f,e),i,f,e
         enddo
         enddo
         enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setup_dg_gs(dg_hndl,nx,ny,nz,nel,melg,vertex)

c     Global-to-local mapping for gs

      include 'SIZE'
      include 'TOTAL'

      integer   dg_hndl
      integer   vertex(1)

      parameter(lf=lx1*lz1*2*ldim*lelt)
      common /c_is1/ glo_num_face(lf)
     $             , glo_num_vol((lx1+2)*(ly1+2)*(lz1+2)*lelt)
      integer*8 glo_num_face,glo_num_vol,ngv,nf

      common /nekmpi/ mid,mp,nekcomm,nekgroup,nekreal

      mx = nx+2
      call set_vert(glo_num_vol,ngv,mx,nel,vertex,.false.)

      mz0 = 1
      mz1 = 1
      if (if3d) mz0 = 0
      if (if3d) mz1 = nz1+1
      call iface_vert_int8 (glo_num_face,glo_num_vol,mz0,mz1,nelt) 

      nf = nx1*nz1*2*ndim*nelt !total number of points on faces
      call gs_setup(dg_hndl,glo_num_face,nf,nekcomm,np)

      return
      end
c-----------------------------------------------------------------------
      subroutine dg_advect_setup
      include 'SIZE'
      include 'TOTAL'
      include 'DG'

      common /ivrtx/ vertex ((2**ldim)*lelt)
      integer vertex

      call setup_dg_gs(dg_hndl,nx1,ny1,nz1,nelt,nelgt,vertex)

      call dg_set_fc_ptr

      return
      end
c-----------------------------------------------------------------------
      subroutine dg_set_fc_ptr
c
c     Set up pointer to restrict u to faces ! NOTE: compact
c
      include 'SIZE'
      include 'TOTAL'
      include 'DG'

      integer e,f,ef

      call dsset(nx1,ny1,nz1) ! set skpdat

      nxyz  = nx1*ny1*nz1
      nxz   = nx1*nz1
      nface = 2*ndim
      nxzf  = nx1*nz1*nface ! red'd mod to area, unx, etc.

      k = 0

      do e=1,nelv
      do f=1,nface

         ef     = eface(f)
         js1    = skpdat(1,f)
         jf1    = skpdat(2,f)
         jskip1 = skpdat(3,f)
         js2    = skpdat(4,f)
         jf2    = skpdat(5,f)
         jskip2 = skpdat(6,f)

         i = 0
         do j2=js2,jf2,jskip2
         do j1=js1,jf1,jskip1

            i = i+1
            k = i+nxz*(ef-1)+nxzf*(e-1)           ! face   numbering
            dg_face(k) = j1+nx1*(j2-1)+nxyz*(e-1) ! global numbering

         enddo
         enddo

      enddo
      enddo
      ndg_face = nxzf*nelv

      return
      end
c-----------------------------------------------------------------------
c-----| Time stepping driver |-----
c-----------------------------------------------------------------------
      subroutine dg_advect_execute   ! do time stepping here? 2+3
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      real cfl 
      logical ifout,ifexc, ifoup 
      real  rhex(lx1,ly1,lz1,lelt)
     $   , rvxex(lx1,ly1,lz1,lelt), rvyex(lx1,ly1,lz1,lelt)
     $   , rvzex(lx1,ly1,lz1,lelt)
     $   ,  Enex(lx1,ly1,lz1,lelt)
      real  mach(lx1,ly1,lz1,lelt), mxma, mxmn
      integer n 

      n=nx1*ny1*nz1*nelt
      ifout = .false.
      ifout = .true.
      ifoup = .true.   ! if output prim 
      ifoup = .false.  ! if output prim 
      ifexc = .true. ! if there is exact solution in usr file 
      ifexc = .false. ! if there is exact solution in usr file 

c    Init for dg
      call setlogi
      call setallic ! 
      if(nid.eq.0) write(6,*) ' miu',miu,', gama',gama,', Pr',prt

      kstep = 0
      call get_duxyz (dux, duy, duz, rh, rvx, rvy, rvz) 
      call printmax(rh,rvx,rvy,rvz,En,kstep)
      call dg_cfl  (cfl,kstep)  ! compute cfl 
      if(nid.eq.0) write(6,*) 'CFL = ', cfl 
      call cmp_mach(mxma,mxmn,mach)  
      if(nid.eq.0) write(6,*) 'Mach number extrema'
      call printmax_1(mach)  

      ifxyo=.true.
      if(ifout) then 
        call outpost(rvx,rvy,rvz,rh,En,'   ')
        if(nid.eq.0) write(6,*) 'write complete', kstep
      endif

c     write(6,*) 'gamma= ', gama, ', miu=', miu
      do kstep=1,nsteps

         call dg_advect(kstep) 

         time = time + dt

         call get_duxyz (dux, duy, duz, rh, rvx, rvy, rvz) 
         call dg_cfl    (cfl,kstep)  ! compute cfl 
         call cmp_mach  (mxma,mxmn,mach)  ! compute mach number  

c        if (mod(kstep,1).eq.0) then
         if (mod(kstep,iostep).eq.0) then
             call printmax(rh,rvx,rvy,rvz,En,kstep)

             if(nid.eq.0) write(6,*) 'Conv. CFL = ', cfl 
             if(nid.eq.0) 
     $         write(6,*) 'Spd of sound CFL = ',2.*cfl/(mxma+mxmn)
             if(nid.eq.0) write(6,*) 'Mach number extrema'
             call printmax_1(mach)  ! print max for Mach 
c   error values during the process is not meaningful
c   vortex already out of domain. Just look at the last one.
             if(ifout) then
                 ifxyo=.true.
                 if (kstep.gt.iostep) ifxyo=.false.
                 call outpost(rvx,rvy,rvz,rh,En,'   ')
                 if(ifoup) then
                   if(if3d) then
                    call get_pres3(ppr,rh,rvx,rvy,rvz,En,gama,n)  ! 3
                   else 
                    call get_pres(ppr,rh,rvx,rvy,En,gama,n)  ! 2
                   endif
                   call outpost(dux,duy,duz,rh,ppr,'prm')
                 endif
                 if(nid.eq.0) write(6,*) 'write complete', kstep
             endif
         endif
      enddo


      return
      end
c-----------------------------------------------------------------------
c----- Cases ----- 
c-----------------------------------------------------------------------
c----- Cylinder ----- 
c-----------------------------------------------------------------------
      subroutine estimate_strouhal

      include 'SIZE'
      include 'TOTAL'

      real tlast,vlast,tcurr,vcurr,t0,t1
      save tlast,vlast,tcurr,vcurr,t0,t1
      data tlast,vlast,tcurr,vcurr,t0,t1 / 6*0 /

      integer e,eg,eg0,e0

      eg0 = 622          ! Identify element/processor in wake
      mid = gllnid(eg0)
      e0  = gllel (eg0)

      st  = 0

      if (nid.eq.mid) then

         tlast = tcurr
         vlast = vcurr

         tcurr = time
         vcurr = vy (1,ny1,1,e0)

         xcurr = xm1(1,ny1,1,e0)
         ycurr = ym1(1,ny1,1,e0)

         write(6,2) istep,time,vcurr,xcurr,ycurr
    2    format(i9,1p4e13.5,' vcurr')

         if (vlast.gt.0.and.vcurr.le.0) then ! zero crossing w/ negative slope
            t0  = t1
            t1  = tlast + (tcurr-tlast)*(vlast-0)/(vlast-vcurr)
            per = t1-t0
            if (per.gt.0) st = 1./per
         endif
      endif

      st = glmax(st,1)

      n  = nx1*ny1*nz1*nelv
      ux = glamax(vx,n)
      uy = glamax(vy,n)

      if (nid.eq.0.and.st.gt.0) write(6,1) istep,time,st,ux,uy
    1 format(i5,1p4e12.4,' Strouhal')

      return
      end
c-----------------------------------------------------------------------
      subroutine set_obj  ! define objects for surface integrals
c
      include 'SIZE'
      include 'TOTAL'

      integer e,f,eg

      nobj = 1
      iobj = 0
      do ii=nhis+1,nhis+nobj
         iobj = iobj+1
         hcode(10,ii) = 'I'
         hcode( 1,ii) = 'F'
         hcode( 2,ii) = 'F'
         hcode( 3,ii) = 'F'
         lochis(1,ii) = iobj
      enddo
      nhis = nhis + nobj

      if (maxobj.lt.nobj) call exitti('increase maxobj in SIZE$',nobj)

      nxyz  = nx1*ny1*nz1
      nface = 2*ndim

      do e=1,nelv
      do f=1,nface
         if (cbc(f,e,1).eq.'W  ') then
            iobj  = 1
            if (iobj.gt.0) then
               nmember(iobj) = nmember(iobj) + 1
               mem = nmember(iobj)
               eg  = lglel(e)
               object(iobj,mem,1) = eg
               object(iobj,mem,2) = f
c              write(6,1) iobj,mem,f,eg,e,nid,' OBJ'
c   1          format(6i9,a4)

            endif
         endif
      enddo
      enddo

c     write(6,*) 'number',(nmember(k),k=1,4)
c
      return
      end
c-----------------------------------------------------------------------
c----- For eddy_uv case -----
c-----------------------------------------------------------------------
      subroutine err_sol_vel(evx,evy,ivx,ivy)
      include 'SIZE'
      include 'TOTAL'
c 
c     Evaluate error between num. & exact solutions vel. fields
c     Input:
c     Outpu:
c        .  Inf norm of errors 
c        .  Second norm of errors 
c    
      real    ivx(1), ivy(1)
      real    evx(1), evy(1)

      real    er1(lx1*ly1*lz1*lelt), er2(lx1*ly1*lz1*lelt)
      real    er0(lx1*ly1*lz1*lelt) 
      real    erni, ern2 

      n    = nx1*ny1*nz1*nelt

      call sub3(er1, evx, ivx, n)
      call sub3(er2, evy, ivy, n)
      call rzero(er0, n)

      call outpost(er1,er2,er0,er0,er0,'err')

c  -- Inf norm -- 
      erni = glamax  (er1,n)
      if(nid.eq.0) write(6,*) 'Inf norm error, ux, ', erni
      erni = glamax  (er2,n)
      if(nid.eq.0) write(6,*) 'Inf norm error, uy, ', erni

c  -- 2nd norm -- 
      ern2 = gl2norm (er1,n)  ! math.f, ln 1146
      if(nid.eq.0) write(6,*) '2nd norm error, ux, ', ern2
      ern2 = gl2norm (er2,n)  ! 
      if(nid.eq.0) write(6,*) '2nd norm error, uy, ', ern2

      if(nid.eq.0) write(6,*) 'done :: Error compute'

      return
      end
c-----------------------------------------------------------------------
c-----| Generic codes for dg flow |----- 
c-----------------------------------------------------------------------
c-----| Init |----- 
c-----------------------------------------------------------------------
      subroutine setlogi ! 2+3
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'

      tmstpp = 1 
      tmstpp = 2      ! second order, check with MATLAB
      tmstpp = 3      ! third order for comparing with Matlab 
      flxtyp = 2      ! 1 - LF , 2 - Roe( not there yet )
      flxtyp = 1      ! 1 - LF , 2 - Roe( not there yet )
c     
      ifstr = .true.  ! strong form 
      ifstr = .false. ! weak  form 

c     de-aliase, currently only the weak form has de-aliasing option
      ifdeal = .true. 
      ifdeal = .false. 

c     diffusion operator 
      ifdifu = .false. ! 
      ifdifu = .true.  ! 

c     debugging, print out function call history 
      ifdbg = .true.   ! for debugging, write out routine stack
      ifdbg = .false.  ! for debugging, write out routine stack

c     limiter 
      iflim = .true.  
      iflim = .false. 

c     for outpost 
      ifvo = .true.
      ifpo = .true.
      ifto = .true.

      if(nid.eq.0) then 
      write(6,*) 'o------------------------------------------------o'
      write(6,*) '| Time-stepper order is           ', tmstpp
      write(6,*) 'o------------------------------------------------o'
      write(6,*) '| Flux type is ( 1 - LF, 2 - Roe) ', flxtyp
      write(6,*) 'o------------------------------------------------o'
      if(ifstr) then 
      write(6,*) '| Strong form                     '
      write(6,*) 'o------------------------------------------------o'
      else
      write(6,*) '| Weak form                       '
      write(6,*) 'o------------------------------------------------o'
      if(ifdeal) then 
        write(6,*) '| De-aliase volume & surface term '
        write(6,*) 'o------------------------------------------------o'
      else
        write(6,*) '| Not de-aliased                  '
        write(6,*) 'o------------------------------------------------o'
      endif
      endif
      if(iflim) then 
        write(6,*) '| Limiter on                      '
        write(6,*) 'o------------------------------------------------o'
      else
        write(6,*) '| Limiter off                     '
        write(6,*) 'o------------------------------------------------o'
      endif
      if(ifdifu) then 
        write(6,*) '| Viscous on (Weak form)          '
        write(6,*) 'o------------------------------------------------o'
      else
        write(6,*) '| Viscous off                     '
        write(6,*) 'o------------------------------------------------o'
      endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine setallic ! 2+3
      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'
      include 'DGUSE'

      n = nx1*ny1*nz1*nelt
c     write(6,*) ifgetx, ifgetu, ifgetp, ifgett 
c     this if is working. for ifdg flag turned on 
      if(       ifgetx             ! if true means mesh restart
     $    .and. ifgetu             ! if true means velo restart
     $    .and. ifgetp             ! if true means dens restart
     $    .and. ifgett      ) then ! if true means ener restart
          ! copy into my data fields, time is set in restarting
          if(if3d) then 
              if(nid.eq.0) write(6,*) 'Done :: DG Restarting 3d...'
              call copy_all(rvx,rvy,rvz, rh,En
     $                    ,  vx, vy, vz, pr, t,n)  ! pres -> density
          else
              call copy_all2(rvx,rvy, rh,En
     $                      , vx, vy, pr, t,n)  ! pres -> density
              if(nid.eq.0) write(6,*) 'in restart '
              call printmax(rh,rvx,rvy,rvz,En,0)
c             if(iflim) then ! Slope limiter on initial field
c                 call slopelim(rh,rvx,rvy,rvz,En)
c             endif
c             stop
              if(nid.eq.0) write(6,*) 'Done :: DG Restarting 2d...'
          endif
      else ! not restarting 
          call setic
          time = 0.
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine setic ! 2+3
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      integer n, nel, nxyz, ie, i, j, k

      n = nx1*ny1*nz1*nelt
      nel = nelfld(ifield)
      nxyz=nx1*ny1*nz1
      do ie=1,nel
          ieg  = lglel(ie)
          do k=1,nz1
          do j=1,ny1
          do i=1,nx1
              call nekasgn(i,j,k,ie)
              call useric_fl(i,j,k,ie) ! get rho, enr, rux,.. ruz
              rh(i,j,k,ie) = rho 
              En(i,j,k,ie) = enr 
              rvx(i,j,k,ie)= rux
              rvy(i,j,k,ie)= ruy
              dux(i,j,k,ie)= rvx(i,j,k,ie)/rh(i,j,k,ie)
              duy(i,j,k,ie)= rvy(i,j,k,ie)/rh(i,j,k,ie)
              if(if3d) then
                  rvz(i,j,k,ie)= ruz
                  duz(i,j,k,ie)= rvy(i,j,k,ie)/rh(i,j,k,ie)
              endif
          enddo 
          enddo 
          enddo 
      enddo 
      if(nid.eq.0) write(6,*) 'done DG :: Set init '

      return
      end
c-----------------------------------------------------------------------
c----- Compute Mach number ----- 
c-----------------------------------------------------------------------
      subroutine cmp_mach(mxa,mxn,mch)  ! compute mach number  
c     Assuming that dux, duy, duz have velo
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      parameter(le=lx1*ly1*lz1) 
      integer kstep, e, i
      real mch(lx1,ly1,lz1,1), mxa, mxn 
      real pres, spd 
      real tmp 

c     Get pressure, p = (gama - 1)*(En - rho ( u^2 + v^2 )/2) 
      mxa = -1.
      mxn = 1.e20 
      do e=1,nelt
      do i=1,le
c       get pres 
         if(if3d) then 
         pres = (gama-1.)*(En(i,1,1,e) - 
     $                  0.5*rh(i,1,1,e)*
     $       (dux(i,1,1,e)**2+duy(i,1,1,e)**2+duz(i,1,1,e)**2) ) 
         else 
         pres = (gama-1.)*(En(i,1,1,e) - 
     $                  0.5*rh(i,1,1,e)*
     $       (dux(i,1,1,e)**2+duy(i,1,1,e)**2) ) 
         endif 
c       Idea gas c = sqrt( gamma * p / rho) 
         if(rh(i,1,1,e) .gt. 1e-10 ) then 
             spd = (gama*pres)/(rh(i,1,1,e))
             if( spd .gt. 1e-13) then
                 spd = sqrt(spd) 
             else
                 write(6,*) 'Near zero or neg. pressure, a^2', spd
                 write(6,*) 'MPI rank = ',nid,' , Elem #. = ', e 
                 write(6,*) 'E',En(i,1,1,e),',rh',rh(i,1,1,e) 
                 write(6,*) 'dux',dux(i,1,1,e),',duy',duy(i,1,1,e) 
                 write(6,*) 'duz',duz(i,1,1,e)
                 call exitt
                 stop 
             endif
         else 
             write(6,*) 'Near zero density, rh', rh(i,1,1,e) 
             write(6,*) 'MPI rank = ',nid,' , Elem #. = ', e 
             write(6,*) 'i = ',i 
             write(6,*) 'E',En(i,1,1,e),',rh',rh(i,1,1,e) 
             write(6,*) 'dux',dux(i,1,1,e),',duy',duy(i,1,1,e) 
             write(6,*) 'duz',duz(i,1,1,e)
             call exitt
             stop 
         endif
         mch(i,1,1,e) = 
     $       sqrt( dux(i,1,1,e)*dux(i,1,1,e)
     $           + duy(i,1,1,e)*duy(i,1,1,e)
     $           + duz(i,1,1,e)*duz(i,1,1,e) ) / spd
         if(mch(i,1,1,e).ge.mxa) mxa = mch(i,1,1,e) 
         if(mch(i,1,1,e).le.mxn) mxn = mch(i,1,1,e) 
      enddo 
      enddo

      return
      end
c-----------------------------------------------------------------------
c----- Compute CFL  ----- 
c-----------------------------------------------------------------------
      subroutine  dg_cfl(cfl,kstep)  ! compute cfl 
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      integer kstep
      real cfl_rt

      call compute_cfl(cfl_rt,dux,duy,duz,dt)  ! cfl from routine 

c     Important: how to compute the CFL based on speed of sound ! 

      cfl = cfl_rt 

      return
      end
c-----------------------------------------------------------------------
c----- Velo. from conservative variables ----- 
c-----------------------------------------------------------------------
      subroutine get_duxyz(dux, duy, duz, rh, rvx, rvy, rvz) 
c     Get velocities data from rh, and momentum 
      include 'SIZE'
      include 'TOTAL'
      parameter(lu = lx1*ly1*lz1*lelt) 
      parameter(le = lx1*ly1*lz1)
      real     dux(1), duy(1), duz(1) 
      real     rvx(1), rvy(1), rvz(1), rh(1) 

      n  = nx1*ny1*nz1*nelt
c     get velocities 
      call invcol3 ( dux, rvx, rh, n)
      call invcol3 ( duy, rvy, rh, n)
      if(if3d) call invcol3 ( duz, rvz, rh, n)

      return
      end
c-----------------------------------------------------------------------
c-----| Time stepper |-----
c-----------------------------------------------------------------------
      subroutine dg_advect(kstep) ! Advance u^{n-1} --> u^n
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'

      integer kstep

      if(tmstpp.eq.1) then      ! euler 

         call dg_advect_euler1(kstep)

      else if(tmstpp.eq.2) then ! rk2

         call dg_advect_ssprk2(kstep)

      else if(tmstpp.eq.3) then ! rk3

         call dg_advect_ssprk3(kstep)

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine dg_advect_euler1(kstep) ! Advance u^{n-1} --> u^n
c     First order time stepper 
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'

      real    tm
      real    rhs1(lt), rhs2(lt), rhs3(lt), rhs4(lt), rhs5(lt)

      n = nx1*ny1*nz1*nelt
      tm = time + dt

      call copyall(tmrh, tmrx, tmry, tmrz, tmen   ! copy into scratch
     $            , rh, rvx, rvy, rvz, En, n)

      call pre_comp      ! use scrdg to find stress, pressure

      call comp_dg_rhs_all(rhs1,rhs2,rhs3,rhs4,rhs5,tm) ! compute rhs 

      do i=1,n       ! update solution 
          rh (i,1,1,1) = rh (i,1,1,1) + rhs1(i)*dt
          rvx(i,1,1,1) = rvx(i,1,1,1) + rhs2(i)*dt
          rvy(i,1,1,1) = rvy(i,1,1,1) + rhs3(i)*dt
          rvz(i,1,1,1) = rvz(i,1,1,1) + rhs4(i)*dt
          En (i,1,1,1) = En (i,1,1,1) + rhs5(i)*dt
      enddo

      if(iflim) then ! Slope limiter on solution field 
          call slopelim(rh,rvx,rvy,rvz,En)
      endif

      if(ifdbg) then
          if(nid.eq.0) 
     $       write(6,*) 'done :: advance Euler time step', kstep
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine dg_advect_ssprk2(kstep) ! Advance u^{n-1} --> u^n
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'

      real    t1, t2
      real    rh1(lt), rvx1(lt), rvy1(lt), rvz1(lt), En1(lt)
     $      , rhs1(lt), rhs2(lt), rhs3(lt), rhs4(lt), rhs5(lt)

      n = nx1*ny1*nz1*nelt
c     Stage 1 
      t1 = time
      call copyall(tmrh, tmrx, tmry, tmrz, tmen   ! copy into scratch
     $            , rh, rvx, rvy, rvz, En, n)
      call pre_comp  ! find pressure, put into scr_dg
      call comp_dg_rhs_all(rhs1,rhs2,rhs3,rhs4,rhs5,t1) ! compute rhs 
      do i=1,n       ! update solution 
          rh1  (i) = rh (i,1,1,1) + rhs1(i)*dt
          rvx1 (i) = rvx(i,1,1,1) + rhs2(i)*dt
          rvy1 (i) = rvy(i,1,1,1) + rhs3(i)*dt
          if(if3d) then 
              rvz1 (i) = rvz(i,1,1,1) + rhs4(i)*dt
          else
              rvz1 (i) = 0. 
          endif
          En1  (i) = En (i,1,1,1) + rhs5(i)*dt
      enddo
      if(iflim) then ! Slope limiter on solution field 
          call slopelim(rh1,rvx1,rvy1,rvz1,En1)
      endif
      if(ifdbg) then
          if(nid.eq.0) 
     $       write(6,*) 'done :: stage 1 after lim', kstep
      endif

c     Stage 2 
      t2 = time + dt
      call copyall(tmrh, tmrx, tmry, tmrz, tmen   ! copy into scratch
     $           , rh1,rvx1,rvy1,rvz1, En1, n)
      call pre_comp               ! find stress, pressure, put into scr_dg
      call comp_dg_rhs_all(rhs1,rhs2,rhs3,rhs4,rhs5,t2) ! compute rhs 
      do i=1,n       ! update solution 
          rh  (i,1,1,1) = (rh (i,1,1,1) + rh1 (i) + rhs1(i)*dt)/2.
          rvx (i,1,1,1) = (rvx(i,1,1,1) + rvx1(i) + rhs2(i)*dt)/2.
          rvy (i,1,1,1) = (rvy(i,1,1,1) + rvy1(i) + rhs3(i)*dt)/2.
          if(if3d) rvz (i,1,1,1) = (rvz(i,1,1,1) + rvz1(i) + rhs4(i)*dt)/2.
          En  (i,1,1,1) = (En (i,1,1,1) + En1 (i) + rhs5(i)*dt)/2.
      enddo
      if(iflim) then ! Slope limiter on solution field 
          call slopelim(rh,rvx,rvy,rvz,En)
      endif

      if(ifdbg) then
          if(nid.eq.0) 
     $       write(6,*) 'done :: stage 2 after lim', kstep
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine dg_advect_ssprk3(kstep) ! Advance u^{n-1} --> u^n
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
c
      real    t1, t2, t3
      real    rh1(lt), rvx1(lt), rvy1(lt), rvz1(lt), En1(lt)
     $     ,  rh2(lt), rvx2(lt), rvy2(lt), rvz2(lt), En2(lt)
     $     , rhs1(lt), rhs2(lt), rhs3(lt), rhs4(lt), rhs5(lt)
      n = nx1*ny1*nz1*nelt
c
c     Stage 1 
c     write(6,*) 'Stage 1 '
      t1 = time
      call copyall(tmrh, tmrx, tmry, tmrz, tmen   ! copy into scratch
     $            , rh, rvx, rvy, rvz, En, n)
      call pre_comp               ! find stress, pressure, put into scr_dg
      call comp_dg_rhs_all(rhs1,rhs2,rhs3,rhs4,rhs5,t1) ! compute rhs 

      do i=1,n       ! update solution 
          rh1  (i) = rh (i,1,1,1) + rhs1(i)*dt
          rvx1 (i) = rvx(i,1,1,1) + rhs2(i)*dt
          rvy1 (i) = rvy(i,1,1,1) + rhs3(i)*dt
          if(if3d) then
              rvz1 (i) = rvz(i,1,1,1) + rhs4(i)*dt
          else
              rvz1 (i) = 0. 
          endif
          En1  (i) = En (i,1,1,1) + rhs5(i)*dt
      enddo
      if(iflim) then ! Slope limiter on solution field 
          call slopelim(rh1,rvx1,rvy1,rvz1,En1)
      endif
c
c     Stage 2 
c     write(6,*) 'Stage 2 ' ! stage 2 breaks down 
      t2 = time + dt
      call copyall(tmrh, tmrx, tmry, tmrz, tmen   ! copy into scratch
     $           , rh1,rvx1,rvy1,rvz1, En1, n)
      call pre_comp               ! find stress, pressure, put into scr_dg
      call comp_dg_rhs_all(rhs1,rhs2,rhs3,rhs4,rhs5,t2) ! compute rhs 
      do i=1,n       ! update solution 
          rh2  (i) = (rh (i,1,1,1)*3. + rh1 (i) + rhs1(i)*dt)/4.
          rvx2 (i) = (rvx(i,1,1,1)*3. + rvx1(i) + rhs2(i)*dt)/4.
          rvy2 (i) = (rvy(i,1,1,1)*3. + rvy1(i) + rhs3(i)*dt)/4.
          if(if3d) then 
              rvz2(i) = (rvz(i,1,1,1)*3. + rvz1(i) + rhs4(i)*dt)/4.
          else 
              rvz2(i) = 0. 
          endif
          En2  (i) = (En (i,1,1,1)*3. + En1 (i) + rhs5(i)*dt)/4.
      enddo
c     OK, problematic rhs for rvy rvz at least. 
c     Thu Aug  6 09:49:13 CDT 2015
      if(iflim) then ! Slope limiter on solution field 
          call slopelim(rh2,rvx2,rvy2,rvz2,En2)
      endif
c
c     Stage 3 
c     write(6,*) 'Stage 3 '
      t3 = time + dt/2.
      call copyall(tmrh, tmrx, tmry, tmrz, tmen   ! copy into scratch
     $           , rh2,rvx2,rvy2,rvz2, En2, n)
      call pre_comp               ! find stress, pressure, put into scr_dg
      call comp_dg_rhs_all(rhs1,rhs2,rhs3,rhs4,rhs5,t3) ! compute rhs 
      do i=1,n       ! update solution 
          rh  (i,1,1,1)=(rh (i,1,1,1)+rh2 (i)*2.+rhs1(i)*dt*2.)/3.
          rvx (i,1,1,1)=(rvx(i,1,1,1)+rvx2(i)*2.+rhs2(i)*dt*2.)/3.
          rvy (i,1,1,1)=(rvy(i,1,1,1)+rvy2(i)*2.+rhs3(i)*dt*2.)/3.
          if(if3d) then
          rvz (i,1,1,1)=(rvz(i,1,1,1)+rvz2(i)*2.+rhs4(i)*dt*2.)/3.
          else
          rvz (i,1,1,1)=0. 
          endif
          En  (i,1,1,1)=(En (i,1,1,1)+En2 (i)*2.+rhs5(i)*dt*2.)/3.
      enddo
      if(iflim) then ! Slope limiter on solution field 
          call slopelim(rh,rvx,rvy,rvz,En)
      endif
c
      return
      end
c-----------------------------------------------------------------------
c-----| End time stepper |-----
c-----------------------------------------------------------------------
c-----| Service routines 1 |-----
c-----------------------------------------------------------------------
      subroutine copyall(ro1,ro2,ro3,ro4,ro5
     $                 , ri1,ri2,ri3,ri4,ri5, n)
c     Copy 5 matrices from riN to roN , size n
      real     ro1(1), ro2(1), ro3(1), ro4(1), ro5(1)
     $       , ri1(1), ri2(1), ri3(1), ri4(1), ri5(1)
      integer  i, n

      call copy(ro1,ri1,n)
      call copy(ro2,ri2,n)
      call copy(ro3,ri3,n)
      call copy(ro4,ri4,n)
      call copy(ro5,ri5,n)

c     do i=1,n
c        ro1(i) = ri1(i)
c        ro2(i) = ri2(i)
c        ro3(i) = ri3(i)
c        ro4(i) = ri4(i)
c        ro5(i) = ri5(i)
c     enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine ctr2full(vol_ary, ctr_ary)

      include 'SIZE'
      include 'TOTAL'
      include 'DG'

      real     vol_ary(lx1,ly1,lz1,lelt)
      real     ctr_ary(lelt)
      integer  i

      ne = nx1*ny1*nz1
      do ie=1,nelt ! fill to volumetric array 
         call cfill(vol_ary(1,1,1,ie),ctr_ary(ie),ne) 
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine ctr2face(faceary, ctr_ary)
c     Input : ctr_ary
c     Output: faceary
c     Temp  : vol_ary
      include 'SIZE'
      include 'TOTAL'

      real     faceary(lx1*lz1,2*ldim,lelt)
      real     vol_ary(lx1,ly1,lz1,lelt)
      real     ctr_ary(lelt)

      call ctr2full (vol_ary, ctr_ary)
      call full2face(faceary, vol_ary)

      return
      end
c-----------------------------------------------------------------------
      subroutine full2face(faceary, vol_ary)

      include 'SIZE'
      include 'TOTAL'
      include 'DG'

      real     faceary(lx1*lz1,2*ldim,lelt)
      real     vol_ary(lx1,ly1,lz1,lelt)
      integer  i,j

      do j=1,ndg_face
         i=dg_face(j)
         faceary(j,1,1) = vol_ary(i,1,1,1)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine face2full(vol_ary, faceary)

      include 'SIZE'
      include 'TOTAL'
      include 'DG'

      real     faceary(lx1*lz1,2*ldim,lelt)
      real     vol_ary(lx1,ly1,lz1,lelt)
      integer  i,j

      n=nx1*ny1*nz1*nelfld(ifield)
      call rzero(vol_ary,n)

      do j=1,ndg_face
         i=dg_face(j)
         vol_ary(i,1,1,1) = vol_ary(i,1,1,1)+faceary(j,1,1)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine add_face2full(vol_ary, faceary)

      include 'SIZE'
      include 'TOTAL'
      include 'DG'      ! dg_face is stored

      real     faceary(lx1*lz1,2*ldim,lelt)
      real     vol_ary(lx1,ly1,lz1,lelt)
      integer  i,j

      do j=1,ndg_face
         i=dg_face(j)
         vol_ary(i,1,1,1) = vol_ary(i,1,1,1)+faceary(j,1,1)
      enddo

      return
      end
c-----------------------------------------------------------------------
c----- Evaluate right hand side  -----
c-----------------------------------------------------------------------
      subroutine comp_dg_rhs_all(rhs1,rhs2,rhs3,rhs4,rhs5,tm)
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      real    rhs1(1), rhs2(1), rhs3(1), rhs4(1), rhs5(1)
      real    tm

      if(ifdifu) then
          call pre_difu_dg        ! solve aux var. 
      endif

      call comp_dg_rhs(rhs1,1,tm) ! compute rhs 
      call comp_dg_rhs(rhs2,2,tm)
      call comp_dg_rhs(rhs3,3,tm)
      if(if3d) call comp_dg_rhs(rhs4,4,tm)
      call comp_dg_rhs(rhs5,5,tm)

      if(ifdbg) then
          if(nid.eq.0) 
     $    write(6,*) 'done :: compute all rhs'
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine comp_dg_rhs(rhs,flg,tm)
c
c    rhs =  volume + surface
c
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      real     rhs(1) , tm
      real     rsf(lt), res(lt)
      real     rhd(lt)
      common /srcns/ uf(lf)
      integer  flg

      n  = nx1*ny1*nz1*nelt
      nf = nx1*nz1*2*ndim*nelt
c 
c output 
      call rzero(rhs,n)
c Temp 
      call rzero(res,n)
      call rzero(rsf,n) 
c 
c     Volume
      call vol_res (rhs,flg)
c 
c     Surface
      call surf_res     ( uf,flg)
      call add_face2full(res,uf)  ! to volumetric array
      call invbf        (rsf,res) ! 3D 
c 
c     Adv done 
      call add2 (rhs,rsf,n) 
c 
c     Diffusion operator 
      if(ifdifu) then
          call difu_dg(rhd, flg)  ! right hand side of diffusion 
          call add2   (rhs,rhd,n) 
      endif
c 
      return
      end
c-----------------------------------------------------------------------
c----- Diffusion operator, BR2 (hopefully) 
c-----------------------------------------------------------------------
      subroutine pre_difu_dg        ! solve aux var. 
c
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      integer i, e 
      n =nx1*ny1*nz1*nelt
      ne=nx1*ny1*nz1
c 
c   Solve for auxi 
      call aux_slv (shr,1) ! 
      call aux_slv (shu,2) ! 
      call aux_slv (shv,3) ! 
      if(if3d) call aux_slv (shw,4) !
      call aux_slv (she,5) ! 
c 
      if(if3d) then 
        write(6,*) 'need equations ' 
      else 
c Find derivatives of velocity 
        do i=1,n
          udx(i,1)= (shu(i,1)- shr(i,1)*tmux(i,1,1,1))/tmrh(i,1,1,1) 
          udy(i,1)= (shu(i,2)- shr(i,2)*tmux(i,1,1,1))/tmrh(i,1,1,1) 
          vdx(i,1)= (shv(i,1)- shr(i,1)*tmuy(i,1,1,1))/tmrh(i,1,1,1) 
          vdy(i,1)= (shv(i,2)- shr(i,2)*tmuy(i,1,1,1))/tmrh(i,1,1,1) 
          edx(i,1)= (she(i,1)- shr(i,1)*tme (i,1,1,1))/tmrh(i,1,1,1) 
          edy(i,1)= (she(i,2)- shr(i,2)*tme (i,1,1,1))/tmrh(i,1,1,1) 
        enddo 
      endif
c 
      return
      end
c-----------------------------------------------------------------------
      subroutine difu_dg(rhs,flg)
c
c    rhs =  volume + surface
c
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      real     rhs(1) 
      integer  flg

      n  = nx1*ny1*nz1*nelt
      nf = nx1*nz1*2*ndim*nelt
c
c   output 
      call rzero(rhs,n)
c
c   Eval F_v, viscous flux term  
      call vis_flx (flg) ! fdf, gdf, hdf  !! 
c
c   Vol integral 
      call vis_vol (rhs, flg) ! use fd, gd, hd  !
c
c   Surf integral 
      call vis_srf (rhs, flg) ! use fdf, gdf, hdf  
c
      return
      end
c-----------------------------------------------------------------------
      subroutine vis_srf(rhs,flg) ! 3+2
c     Central flux, only weak form 
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      integer  flg, i, f, e, k
      real     flx(1) ! defined on surface 
      real     ucf(lf,1) 
      real     rhf(lf), tm1(lt), tm2(lt) 
      real     taud(lf) 
      logical  ifctr  ! if central flux 
      logical  ifldg  ! if central flux 
      real     btax(lf), btay(lf), btaz(lf) 
      ifctr = .true.  ! Central needs penalty on u
      ifctr = .false. ! Central needs penalty on u
c   ifctr, ifldg = false -> default: BR 1, real central 
c     ifldg = .true.  ! LDG basedo on central, need ifctr for ifldg
      ifldg = .false. ! LDG basedo on central, need ifctr for ifldg

      nf=nx1*nz1*2*ndim*nelt
      n=nx1*ny1*nz1*nelt
      nxz=nx1*nz1
      nfaces=2*ndim

      if(flg.ge.2) then 
        call vis_ctr_f(flg) ! fcd, gcd, hcd 
        call rzero(dud,nf)  ! dud 
        if(ifctr) call vis_dif_u(taud,flg) ! taud, dud 
        if(ifldg) call cmp_beta (btax,btay,btaz) !
        if(ifldg) call vis_dif_f(flg) ! fdd, gdd, hdd

        k = 0
        do e=1,nelt
        do f=1,nfaces
        do i=1,nxz
           k=k+1
           if(if3d) then
              rhf(k) = ( fcd(k)*unx(i,1,f,e) 
     $                 + gcd(k)*uny(i,1,f,e) 
     $                 + hcd(k)*unz(i,1,f,e) ) *area(i,1,f,e) 
              if(ifctr) then 
                rhf(k) = rhf(k) - taud(k)* ! change to - 
     $                 ( dud(k)*unx(i,1,f,e) 
     $                 + dud(k)*uny(i,1,f,e) 
     $                 + dud(k)*unz(i,1,f,e) ) *area(i,1,f,e) 
              endif
              if(ifldg) then 
                rhf(k) = rhf(k) - 
     $                 ( fdd(k)*btax(k) 
     $                 + gdd(k)*btay(k) 
     $                 + hdd(k)*btaz(k) ) *area(i,1,f,e) 
              endif
           else ! 2D 
              rhf(k) = ( fcd(k)*unx(i,1,f,e) 
     $                 + gcd(k)*uny(i,1,f,e) ) *area(i,1,f,e) 
              if(ifctr) then 
                rhf(k) = rhf(k) - taud(k)* ! change to - 
     $                 ( dud(k)*unx(i,1,f,e) 
     $                 + dud(k)*uny(i,1,f,e) ) *area(i,1,f,e) 
c     Is this different from ? 
c               rhf(k) = rhf(k) - taud(k)*dud(k)*area(i,1,f,e) 
              endif
              if(ifldg) then 
                rhf(k) = rhf(k) - 
     $                 ( fdd(k)*btax(k) 
     $                 + gdd(k)*btay(k) ) *area(i,1,f,e) 
              endif
           endif
        enddo
        enddo
        enddo

        call face2full (tm1,rhf) 
        call invbf     (tm2,tm1) ! all elements 

        call add2      (rhs,tm2,n) 
      endif ! if flg ge 2 

      return
      end
c-----------------------------------------------------------------------
      subroutine vis_dif_f(flg)
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      integer  flg
      real     ua (lt)
      real     bf1(lf),bg1(lf),bh1(lf) ! fluxes on surfaces 
      real     tmf(lf) 

      n=nx1*ny1*nz1*nelt
      nf=nx1*nz1*2*ndim*nelt !total number of points on faces
c 
c    Use fd, gd, hd

      if(flg.ge.2) then 
        call full2face(fdf, fd) 
        call full2face(gdf, gd) 
        if(if3d) call full2face(hdf, hd) 

c    boundary conditions on viscous flux term 
        call bc_flx_vis(bf1,bg1,bh1,flg) ! Page 249 

        call copy  (tmf,fdf,nf)         ! a 
        call cmult (tmf,2.0,nf)         ! times 2.
        call add2  (fdf,bf1,nf)         ! first add bc , a- + a+(bc)
        call gs_op (dg_hndl,fdf,1,1,0)  ! 1 ==> +      , a- + a+(itr)
        call chsign(fdf,nf)             ! - ( a + b ) 
        call add2  (fdf,tmf,nf)         !  

        call copy  (tmf,gdf,nf)         ! a 
        call cmult (tmf,2.0,nf)         ! times 2.
        call add2  (gdf,bf1,nf)         ! first add bc , a- + a+(bc)
        call gs_op (dg_hndl,gdf,1,1,0)  ! 1 ==> +      , a- + a+(itr)
        call chsign(gdf,nf)             ! - ( a + b ) 
        call add2  (gdf,tmf,nf)         ! 

        if(if3d) then
          call copy  (tmf,hdf,nf)         ! a 
          call cmult (tmf,2.0,nf)         ! times 2.
          call add2  (hdf,bf1,nf)         ! first add bc , a- + a+(bc)
          call gs_op (dg_hndl,hdf,1,1,0)  ! 1 ==> +      , a- + a+(itr)
          call chsign(hdf,nf)             ! - ( a + b ) 
          call add2  (hdf,tmf,nf)         ! 
        endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine vis_dif_u(tau,flg)
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      integer  flg
      real     ua (lt)
      real     buf(lf), ufsv(lf) ! conservative variable on surfaces 
      real     tau(lf) 

      n=nx1*ny1*nz1*nelt
      nf=nx1*nz1*2*ndim*nelt !total number of points on faces
c 
c    Eval ua based on flg 
      call aux_u     ( ua, flg)  !
      call full2face (dud, ua)  !
c    
      call cmp_tau (tau, flg)    !
c
      call rzero(buf, nf)
      call bc_aux_eval(buf,flg) ! bc for u
c 
      call copy  (ufsv,dud,nf)  ! a 
      call cmult (ufsv,2.,nf)   ! 2 a  
c
      call add2  (dud,buf,nf)   ! first add bc 
      call gs_op (dg_hndl,dud,1,1,0)  ! 1 ==> + , a + b 
      call chsign(dud,nf)       ! - ( a + b ) 
      call add2  (dud,ufsv,nf)  ! 2 a - ( a + b ) = a - b

      return
      end
c-----------------------------------------------------------------------
      subroutine cmp_beta(btx,bty,btz) ! 3+2
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      parameter(le = lx1*ly1*lz1)
      integer i, f, e, k
      real    btx(1), bty(1), btz(1) 

      nxz = nx1*nz1
      ne = nx1*ny1*nz1
      nfaces=2*ndim
      nf = nx1*nz1*2*ndim*nelt 
      n  = nx1*ny1*nz1*nelt

c     Get beta =  unx, uny, unz ?  
      if(if3d) then 
        call copy(btx, unx, nf) ! for 2d not correct 
        call copy(bty, uny, nf)
        call copy(btz, unz, nf)
      else
c  -- 1.
c  -- 2.
        k=0
        do e=1,nelt
        do f=1,nfaces
          do i=1,nxz
            k = k + 1 
c           btx(k) = unx(i,1,f,e) 
c           bty(k) = uny(i,1,f,e) 
            btx(k) = -1.*unx(i,1,f,e) 
            bty(k) = -1.*uny(i,1,f,e) 
          enddo 
        enddo 
        enddo 
      endif
c 
      return
      end
c-----------------------------------------------------------------------
      subroutine cmp_tau(tau,flg) ! 3+2
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      parameter(le = lx1*ly1*lz1)
      real    tav(lt)
      integer e 

      ne = nx1*ny1*nz1
      n  = nx1*ny1*nz1*nelt

c     Get tau =  ? 
c  -- 1.
c     call rzero(tav,n) ! tau = 0, reduces to central 
c  -- 2. 
      call rone (tav,n) ! tau = 1, control jump 
      call cmult(tav,0.5,n) ! tau = .5, control jump 
c 
      call full2face(tau,tav)         ! lm always > 0 
c     call gs_op (dg_hndl,lm,1,4,0)  ! 4 ==> max, same on interfaces
c     all zeros or ones, no need to do max 

      return
      end
c-----------------------------------------------------------------------
      subroutine vis_ctr_f(flg) 
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      integer  flg 
      real     bf1(lf),bg1(lf),bh1(lf) ! fluxes on surfaces 

      n=nx1*ny1*nz1*nelt
      nf=nx1*nz1*2*ndim*nelt !total number of points on faces

      if(flg.ge.2) then 
        call full2face(fdf, fd) 
        call full2face(gdf, gd) 
        if(if3d) call full2face(hdf, hd) 

c    boundary conditions on viscous flux term 
        call bc_flx_vis(bf1,bg1,bh1,flg) ! Page 249 

        call copy  (fcd,fdf,nf)         ! a 
        call add2  (fcd,bf1,nf)         ! first add bc , a- + a+(bc)
        call gs_op (dg_hndl,fcd,1,1,0)  ! 1 ==> +      , a- + a+(itr)
        call cmult (fcd,0.5,nf)         ! times 1/2, (a- + a+)/2 

        call copy  (gcd,gdf,nf)         ! a 
        call add2  (gcd,bg1,nf)         ! first add bc , a- + a+(bc)
        call gs_op (dg_hndl,gcd,1,1,0)  ! 1 ==> +      , a- + a+(itr)
        call cmult (gcd,0.5,nf)         ! times 1/2, (a- + a+)/2 

        if(if3d) then
          call copy  (hcd,hdf,nf)         ! a 
          call add2  (hcd,bh1,nf)         ! first add bc , a- + a+(bc)
          call gs_op (dg_hndl,hcd,1,1,0)  ! 1 ==> +      , a- + a+(itr)
          call cmult (hcd,0.5,nf)         ! times 1/2, (a- + a+)/2 
        endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine bc_flx_vis(bf1,bg1,bh1,flg) ! 3 + 2
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      common /nekcb/ cb
      character cb*3
      integer  flg, i, f, e, k
      real     bf1(lf), bg1(1), bh1(1) 
      n  = nx1*ny1*nz1*nelt
      nf = nx1*nz1*2*ndim*nelt 
      nxz= nx1*nz1
      nfaces=2*ndim
      call rzero(bf1,nf)
      call rzero(bg1,nf)
      call rzero(bh1,nf)
      k=0
      do e=1,nelt
      do f=1,nfaces
         ieg=lglel(e)
         cb =cbc(f,e,ifield)
         if(cb.eq.'f  ' .or. cb.eq.'F  ' .or. 
     $      cb.eq.'o  ' .or. cb.eq.'SYM' .or.  
     $      cb.eq.'O  ' .or. 
     $      cb.eq.'w  ' .or. cb.eq.'W  ' ) then
            ia=0
            call facind(kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,f)
            do iz=kz1,kz2
            do iy=ky1,ky2
            do ix=kx1,kx2
               ia = ia + 1
               k  = ia + (f-1)*nxz + (e-1)*nfaces*nxz

               if(cb.eq.'f  ') then  ! qh+ = qh- 

c                call nekasgn(ix,iy,iz,e)
c                call userbc_f(ix,iy,iz,f,e,ieg) 

                 bf1(k) = fdf(k) 
                 bg1(k) = gdf(k)  ! BC for F_v
                 if(if3d) bh1(k) = hdf(k) 
c              else if(cb.eq.'o  ') then
c                  call nekasgn(ix,iy,iz,e)
c                  call userbc_o(ix,iy,iz,f,e,ieg)
               else if(cb.eq.'SYM') then ! slip wall 
                   call nekasgn(ix,iy,iz,e)
                   call userbc_sym(ix,iy,iz,f,e,ieg)
               else if(cb.eq.'W  '.or.cb.eq.'w  ') then ! slip wall 
c if no boundary conditions on \nabla u \cdot n 
c set bf1 = fdf 
                 bf1(k) = fdf(k) 
                 bg1(k) = gdf(k)  ! BC for F_v
                 if(if3d) bh1(k) = hdf(k) 
c                  call nekasgn(ix,iy,iz,e)
c                  call userbc_w_vis(ix,iy,iz,f,e,ieg)
               endif
c              if(if3d) then 
c                  call surf_pnt(bf1(k),bg1(k),bh1(k)
c    $             , rho, rux, ruy, ruz, enr, pp, gama, flg)
c              else  
c                  call surf_pnt2(bf1(k),bg1(k)
c    $             , rho, rux, ruy, enr, pp, gama, flg)
c              endif
            enddo 
            enddo 
            enddo 
         endif
      enddo 
      enddo 

      return
      end
c-----------------------------------------------------------------------
      subroutine vis_vol(rhs,flg) ! ! 3+2? 
c
c  weak form: 1. times mass matrix 
c             2. transpose of derivatives 
c             3. times mass matrix inverse 
c
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      parameter(le=lx1*ly1*lz1)
      integer  flg, e
      real     rhs(le,1), tm2(le,lelt)

      n=nx1*ny1*nz1*nelt
      ne=nx1*ny1*nz1
      
      if(if3d) then
c       write(6,*) '3d not ther e' 
        call rzero(tm2,n)
        call grad_bm1_t(tm2,fd,gd,hd)
      else  ! first use code that's already there, will it work? 
        call rzero(tm2,n)
        call rzero(hd,n)
        if(flg.ge.2) then 
          call grad_bm1_t(tm2,fd,gd,hd) ! no de-aliasing now 
        endif
      endif

      if(flg.ge.2) then 
        call invbf (rhs,tm2) ! all elements 
        call chsign(rhs,n) ! all elements 
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine vis_flx(flg)  ! 3 + 2 
c     Assume : shr, shu, shv, shw, she already up to date
c     Assume : tmux, tmuy, tmuz, tme already up to date
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      integer flg 
      ne = nx1*ny1*nz1
      n  = nx1*ny1*nz1*nelt

      if(if3d) then 
        write(6,*) ' not there yet' 
        call exitt
        stop 
      else  ! 2 
        if (flg .eq. 1) then 
          call rzero ( fd,n)
          call rzero ( gd,n)
        else if(flg .eq. 2) then 
          do i=1,n 
            fd(i) = 2.*udx(i,1) - 2.*(udx(i,1) + vdy(i,1))/3. 
            gd(i) =                   udy(i,1) + vdx(i,1)
          enddo 
        else if(flg .eq. 3) then 
          do i=1,n 
            fd(i) =                   udy(i,1) + vdx(i,1)
            gd(i) = 2.*vdy(i,1) - 2.*(udx(i,1) + vdy(i,1))/3. 
          enddo 
        else if(flg .eq. 5) then  ! 
          do i=1,n 
            fd(i) = tmux(i,1,1,1)*
     $             (2.*udx(i,1) - 2.*(udx(i,1) + vdy(i,1))/3.)
     $            + tmuy(i,1,1,1)*( udy(i,1) + vdx(i,1) )
c    $            + gama*edx(i,1) / prt  
            gd(i) = tmux(i,1,1,1)*( udy(i,1) + vdx(i,1) )
     $            + tmuy(i,1,1,1)*
     $             (2.*vdy(i,1) - 2.*(udx(i,1) + vdy(i,1))/3.)
c    $            + gama*edy(i,1) / prt  
          enddo 
        endif
        call cmult(fd,miu,n) 
        call cmult(gd,miu,n) 
      endif

      return
      end
c-----------------------------------------------------------------------
c----- Auxillary variable for diffusion 
c-----------------------------------------------------------------------
      subroutine aux_slv(sh,flg) ! 
c    Only do weak form ? 
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      real     sh (lt,1) 
      integer  flg
      real     ua (lt), vua(lt,ldim) 
      real     uaf(lf,ldim) 

      n  = nx1*ny1*nz1*nelt
      nf = nx1*nz1*2*ndim*nelt

c    Eval ua based on flg 
      call aux_u   ( ua, flg)  !

c    Vol integral 
      call aux_vol ( vua, ua)  ! 

c    Surf integral 
      call aux_srf ( uaf, ua,flg)  ! 

c    Get sh 
      call aux_s   ( sh, vua, uaf)  ! 

      return
      end
c-----------------------------------------------------------------------
      subroutine aux_s(sh,vua,uaf) ! 3+2
c     
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      integer  flg, i, f, e, k
      real     ftmp(lt) ! surf flux, def on volume 
      real     sh(lt,1), vua(lt,1)
      real     uaf(lf,1)
c
      n=nx1*ny1*nz1*nelt
      nxz=nx1*nz1
      nfaces=2*ndim
c
      do id=1,ndim  ! 3, 2 
        call face2full (ftmp,uaf(1,id)) 
        call sub2      (ftmp,vua(1,id),n)  ! surf - vol 
c    Multi. mass matrix inverse 
        call invbf     (sh(1,id), ftmp)  ! inv bm1 
      enddo 
c
      return
      end
c-----------------------------------------------------------------------
      subroutine aux_srf(ucf,ua,flg) ! 3+2
c     Central flux, only weak form 
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      integer  flg, i, f, e, k
      real     flx(1) ! defined on surface 
      real     ucf(lf,1) 
      real     ua (lt)
c 
      logical  ifldg  ! if local dg 
      real     btax(lf), btay(lf), btaz(lf) 

      ifldg = .true. 
      ifldg = .false. 
      nf=nx1*nz1*2*ndim*nelt !total number of points on faces
      n=nx1*ny1*nz1*nelt
      nxz=nx1*nz1
      nfaces=2*ndim

      call rzero (btax, nf) 
      call rzero (btay, nf) 
      if(if3d) call rzero(btaz,nf)  ! beta 
      if(ifldg) call cmp_beta (btax,btay,btaz)  ! beta, page 267, H&W

      call aux_ctr_u(ua,flg) !!  duc , flg for bc 
      call aux_dif_u(flg)    !!  dud 

      k = 0
      do e=1,nelt
      do f=1,nfaces
      do i=1,nxz
         k=k+1
         if(if3d) then
            ucf(k,1) = duc(k)*unx(i,1,f,e)*area(i,1,f,e)
            ucf(k,2) = duc(k)*uny(i,1,f,e)*area(i,1,f,e)
            ucf(k,3) = duc(k)*unz(i,1,f,e)*area(i,1,f,e)
            if(ifldg) then
              ucf(k,1) = ucf(k,1) 
     $                 + dud(k)*btax(k)*area(i,1,f,e)
              ucf(k,2) = ucf(k,2) 
     $                 + dud(k)*btay(k)*area(i,1,f,e)
              ucf(k,3) = ucf(k,3) 
     $                 + dud(k)*btaz(k)*area(i,1,f,e)
            endif
         else ! 2D 
            ucf(k,1) = duc(k)*unx(i,1,f,e)*area(i,1,f,e)
            ucf(k,2) = duc(k)*uny(i,1,f,e)*area(i,1,f,e)
            if(ifldg) then
              ucf(k,1) = ucf(k,1) 
     $                 + dud(k)*btax(k)*area(i,1,f,e)
              ucf(k,2) = ucf(k,2) 
     $                 + dud(k)*btay(k)*area(i,1,f,e)
            endif
         endif
      enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine aux_dif_u(flg)
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      integer  flg
      real     ua (lt)
      real     buf(lf), ufsv(lf) ! conservative variable on surfaces 
      real     tau(lf) 

      n=nx1*ny1*nz1*nelt
      nf=nx1*nz1*2*ndim*nelt !total number of points on faces
c 
c    Eval ua based on flg 
      call aux_u     ( ua, flg)  !
      call full2face (dud, ua)  !
c
      call rzero(buf, nf)
      call bc_aux_eval(buf,flg) ! bc for u
c 
      call copy  (ufsv,dud,nf)  ! a 
      call cmult (ufsv,2.,nf)   ! 2 a  
c
      call add2  (dud,buf,nf)   ! first add bc 
      call gs_op (dg_hndl,dud,1,1,0)  ! 1 ==> + , a + b 
      call chsign(dud,nf)       ! - ( a + b ) 
      call add2  (dud,ufsv,nf)  ! 2 a - ( a + b ) = a - b
c 
      return
      end
c-----------------------------------------------------------------------
      subroutine aux_ctr_u(ua,flg)
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      integer  flg 
      real     buf(lf), uf (lf) 
      real     ua(lt) 

      n=nx1*ny1*nz1*nelt
      nf=nx1*nz1*2*ndim*nelt !total number of points on faces

      call full2face   (duc, ua) 
      call bc_aux_eval (buf,flg) ! Page 249, H&W for dirichlet 

      call add2  (duc,buf,nf)         ! first add bc , a- + a+(bc) 
      call gs_op (dg_hndl,duc,1,1,0)  ! 1 ==> +      , a- + a+(itr)
      call cmult (duc,0.5,nf)         ! times 1/2, (a- + a+)/2 

      return
      end
c-----------------------------------------------------------------------
      subroutine bc_aux_eval(buf,flg) ! 3+2
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      common /nekcb/ cb
      character cb*3
      integer  flg, i, f, e, k
      real     buf(lf)
      n   =nx1*ny1*nz1*nelt
      nf  =nx1*nz1*2*ndim*nelt 
      nxz =nx1*nz1
      nfaces=2*ndim
      k=0
      do e=1,nelt
      do f=1,nfaces
         ieg=lglel(e)
         cb =cbc(f,e,ifield)
         if(cb.eq.'f  ' .or. cb.eq.'F  ' .or. 
     $      cb.eq.'w  ' .or. cb.eq.'W  ' .or. 
     $      cb.eq.'o  ' .or. cb.eq.'O  ') then
            ia=0
            call facind(kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,f)
            do iz=kz1,kz2
            do iy=ky1,ky2
            do ix=kx1,kx2
               ia = ia + 1
               k  = ia + (f-1)*nxz + (e-1)*nfaces*nxz
               if(cb.eq.'f  '.or. cb .eq.'F  ') then 
                   call nekasgn(ix,iy,iz,e)
                   call userbc_f(ix,iy,iz,f,e,ieg) ! rho, rux, ruy..
                   if(if3d) then
                       call surf_pnt_conv(buf(k)
     $             , rho, rux, ruy, ruz, enr, pp, gama, flg)
                   else 
                       call surf_pnt_conv2(buf(k)
     $             , rho, rux, ruy, enr, pp, gama, flg)
                   endif
                   if(flg.eq.1) then
                     buf(k) = 2.*buf(k) - tmrh(ix,iy,iz,e) 
                   else if(flg.eq.2) then 
                     buf(k) = 2.*buf(k) - tmrx(ix,iy,iz,e) 
                   else if(flg.eq.3) then 
                     buf(k) = 2.*buf(k) - tmry(ix,iy,iz,e) 
                   else if(flg.eq.4) then 
                     buf(k) = 2.*buf(k) - tmrz(ix,iy,iz,e) 
                   else if(flg.eq.5) then 
                     buf(k) = 2.*buf(k) - tmen(ix,iy,iz,e) 
                   endif
               else if(cb.eq.'W  ') then ! ! solid wall 
                   call nekasgn(ix,iy,iz,e)
                   call userbc_w(ix,iy,iz,f,e,ieg) ! rho, rux, ruy..
                   if(if3d) then
                       call surf_pnt_conv(buf(k)
     $             , rho, rux, ruy, ruz, enr, pp, gama, flg)
                   else 
                       call surf_pnt_conv2(buf(k)
     $             , rho, rux, ruy, enr, pp, gama, flg)
                   endif
c                  if(flg.eq.1) then
c                    buf(k) = 2.*buf(k) - tmrh(ix,iy,iz,e) 
c                  else if(flg.eq.2) then 
c                    buf(k) = 2.*buf(k) - tmrx(ix,iy,iz,e) 
c                  else if(flg.eq.3) then 
c                    buf(k) = 2.*buf(k) - tmry(ix,iy,iz,e) 
c                  else if(flg.eq.4) then 
c                    buf(k) = 2.*buf(k) - tmrz(ix,iy,iz,e) 
c                  else if(flg.eq.5) then 
c                    buf(k) = 2.*buf(k) - tmen(ix,iy,iz,e) 
c                  endif
               else if(cb.eq.'o  ' .or. cb.eq.'O  ') then ! 
                   write(6,*) 'What happens at outflow? '
                   call exitt 
               endif
            enddo 
            enddo 
            enddo 
         endif
      enddo 
      enddo 

      return
      end
c-----------------------------------------------------------------------
      subroutine surf_pnt_aux(cu 
     $           , rho, rux, ruy, ruz, en, pp, gama, flg)
c     
c     -  Input  :  rho, rux, ruy, ruz, en, flg
c     -  Output :  cu, pp
      integer  flg
      real     cu
      real     rho, rux, ruy, ruz, en, pp, gama 

      if(flg.eq.1) then
          cu = rho   ! rho 
      else if(flg.eq.2) then
          cu = rux   ! rh u
      else if(flg.eq.3) then
          cu = ruy   ! rh v
      else if(flg.eq.4) then
          cu = ruz   ! rh w
      else if(flg.eq.5) then
          cu = en    ! en  
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine surf_pnt_aux2(cu 
     $           , rho, rux, ruy, en, pp, gama, flg)
c     Evaluate conservative variables to cu 
c     No common blocks 
c     -  Input  :  rho, rux, ruy, ruz, en, flg
c     -  Output :  cu, pp
      integer  flg
      real     cu
      real     rho, rux, ruy, ruz, en, pp, gama 

      pp = (gama-1.)*(en - 0.5*(rux*rux+ruy*ruy+ruz*ruz)/rho) 
      if(flg.eq.1) then
          cu = rho   ! rho 
      else if(flg.eq.2) then
          cu = rux   ! rh u
      else if(flg.eq.3) then
          cu = ruy   ! rh v
c     else if(flg.eq.4) then
c         cu = ruz   ! rh w
      else if(flg.eq.5) then
          cu = en    ! en  
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine aux_vol(vu,ua) ! ! 3+2? 
c
c  weak form: 1. times mass matrix 
c             2. transpose of derivatives 
c
c     uses w3m1, the diagnonal 'vanilla' mass matrix
c
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      parameter(le=lx1*ly1*lz1)
      integer  flg, e
      real     vu(lt,1)
      real     ua(le,1), mua(le,lelt) 

      n=nx1*ny1*nz1*nelt
      ne=nx1*ny1*nz1
      
      do e=1,nelt 
        call col3 (mua(1,e),ua(1,e),w3m1,le) ! multi mass matrix 
      enddo 

      if(if3d) then 
        call grad_t_u (vu(1,1),vu(1,2),vu(1,3),mua) ! transpose of grad
      else 
        call grad_t_u2(vu(1,1),vu(1,2),mua)
      endif

      return
      end
c-----------------------------------------------------------------------
c----- Service routine for diffusion operator 
c-----------------------------------------------------------------------
      subroutine grad_t_u(ux,uy,uz,mua) !  3
c
c     Input : M*u 
c     Output: ux, uy, uz 
c   DxT (M fx), (M fy) Dy , DzT (M fz) 
      include 'SIZE'
      include 'TOTAL'
c
      parameter (lxyz=lx1*ly1*lz1)
      real mua(lxyz,1)
      real ux (lxyz,1),uy(lx1,ly1,lz1,1),uz(lxyz,1)
      real fr (lxyz), fs(lx1,ly1,lz1), ft(lxyz) 
      integer e, i, k
      integer m1, m2 
      m1 = lx1 
      m2 = m1*m1 
c
      do e = 1,nelt
c   
         do i=1,lx1*ly1*lz1
             fr(i) = ( rxm1(i,1,1,e)
     $             +   rym1(i,1,1,e)
     $             +   rzm1(i,1,1,e) )*mua(i,e) ! jacm1 
             fs(i,1,1) = ( sxm1(i,1,1,e)
     $             +   sym1(i,1,1,e)
     $             +   szm1(i,1,1,e) )*mua(i,e) ! jacm1 
             ft(i) = ( txm1(i,1,1,e)
     $             +   tym1(i,1,1,e)
     $             +   tzm1(i,1,1,e) )*mua(i,e) ! jacm1 
         enddo 
c   
c      dxtm1 fr, (do loop dym1) fs, ft dzm1
c      Similar to local_grad3_t(v,fr,fs,ft,lx1-1,1,dxm1,dxtm1,w2)
c 
         call mxm(dxtm1,m1,fr,m1,ux(1,e),m2)  ! x 
         do k=1,m1
           call mxm(fs(1,1,k),m1,dym1,m1,uy(1,1,k,e),m1)  ! y 
         enddo 
         call mxm(ft,m2,dzm1,m1,uz(1,e),m1)  ! z 
c
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine grad_t_u2(ux,uy,mua) !  2
c
c     Input : M*u 
c     Output: ux, uy
c   DxT (M fx), (M fy) Dy
      include 'SIZE'
      include 'TOTAL'
c
      parameter (lxy=lx1*ly1)
      real mua(lxy,1)
      real ux (lxy,1), uy(lxy,1)
      real fr (lxy)  , fs(lxy)
      integer e, i, m1
      m1 = lx1 
c
      do e = 1,nelt
c   
         do i=1,lx1*ly1
             fr(i) = ( rxm1(i,1,1,e)
     $             +   rym1(i,1,1,e) )*mua(i,e) ! jacm1 
             fs(i) = ( sxm1(i,1,1,e)
     $             +   sym1(i,1,1,e) )*mua(i,e) ! jacm1 
         enddo 
c   
         call mxm(dxtm1,m1,fr,m1,ux(1,e),m1)  ! x 
         call mxm(fs,m1,dym1,m1,uy(1,e),m1)   ! y 
c
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine aux_u(ua, flg)
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      integer  flg
      real     ua(1) 

      n=nx1*ny1*nz1*nelt
      if(flg.eq.1) then
          call copy    (ua,tmrh,n)    ! rho
      else if(flg.eq.2) then
          call copy    (ua,tmrx,n)    ! rho u
      else if(flg.eq.3) then
          call copy    (ua,tmry,n)    ! rho v
      else if(flg.eq.4) then
          call copy    (ua,tmrz,n)    ! rho w
      else if(flg.eq.5) then
          call copy    (ua,tmen,n)    ! energy
      endif

      return
      end
c-----------------------------------------------------------------------
c----- End diffusion 
c-----------------------------------------------------------------------
c----- Service routine for main loop 
c-----------------------------------------------------------------------
      subroutine pre_comp
c     Assume : tmrh, tmrx, tmry, tmrz, tmen already up to date
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      parameter(lu = lx1*ly1*lz1*lelt) 
      parameter(le = lx1*ly1*lz1)
      real     ur(le), us(le), ut(le), vr(le), vs(le), vt(le)
     $       , wr(le), ws(le), wt(le)
      integer  e 
      ne = nx1*ny1*nz1
      n  = nx1*ny1*nz1*nelt

c     energy 
      call invcol3 ( tme , tmen, tmrh, n)

c     velocities 
      call invcol3 ( tmux, tmrx, tmrh, n)
      call invcol3 ( tmuy, tmry, tmrh, n)
      if(if3d) call invcol3 ( tmuz, tmrz, tmrh, n)
      if(if3d) then 
c     Get pressure, p = (gama - 1)*(En - rho ( u^2 + v^2 )/2) 
      do e=1,nelt
      do i=1,le
         tmpr(i,1,1,e) = (gama-1.)*(tmen(i,1,1,e) - 
     $                  0.5*tmrh(i,1,1,e)*
     $       (tmux(i,1,1,e)**2+tmuy(i,1,1,e)**2+tmuz(i,1,1,e)**2) ) 
      enddo 
      enddo
      else if(ifaxis) then ! for 2d axial 
          write(6,*) 'Not axial yet, quit in pre_comp'
          call exitt
      else                 ! for 2d, no axial though
      do e=1,nelt
      do i=1,le
         tmpr(i,1,1,e) = (gama-1.)*(tmen(i,1,1,e) - 
     $       0.5*tmrh(i,1,1,e)*(tmux(i,1,1,e)**2+tmuy(i,1,1,e)**2)) 
      enddo 
      enddo
      endif
      if(ifdbg) then
          if(nid.eq.0) write(6,*) 'done :: pre compute' 
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine get_fpr(fpr) ! 3
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      parameter(ld=lxd*lyd*lzd)
      parameter(ltd=lxd*lyd*lzd*lelt)
      real     fpr(ld,1)

      mx = lx1
      md = lxd
! is this the correct way? feel like this is already aliased
      do ie=1,nelt 
         call intp_rstd(fpr(1,ie),tmpr(1,1,1,ie),mx,md,if3d,0) 
c forward, GLL to GL
      enddo 

      return
      end
c-----------------------------------------------------------------------
      subroutine vol_eval_d(fx1,fx2,gx1,gx2,hx1,hx2
     $                     , flg) ! 3
c 
c         fx1    fx2    gx1    gx2     hx1    hx2   
c  r      r u     1     r v     1      r w     1    
c  ru      u     r u     v     r u      w     r u   
c  rv      u     r v     v     r v      w     r v   
c  rw      u     r w     v     r w      w     r w   
c  en      u     E+p     v     E+p      w     E+p   
c 
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      integer  flg
      real     f1(lt), g1(lt), h1(lt)
      real     fx1(lt), fx2(lt), gx1(lt), gx2(lt), hx1(lt), hx2(lt) 
      real     uc(lt), vc(lt), wc(lt), cpr(lt)

      n=nx1*ny1*nz1*nelt
      if(flg.eq.1) then
          call copy    (fx1,tmrx,n)    ! rh u
          call rone    (fx2, n) 
          call copy    (gx1,tmry,n)    ! rh v
          call rone    (gx2, n) 
          call copy    (hx1,tmrz,n)    ! rh w
          call rone    (hx2, n) 
      else if(flg.eq.2) then
          call copy    (fx1, tmux, n)  ! u
          call copy    (fx2, tmrx, n)  ! rh u
          call copy    (gx1, tmuy, n)  ! v
          call copy    (gx2, tmrx, n)  ! rh u
          call copy    (hx1, tmuz, n)  ! w
          call copy    (hx2, tmrx, n)  ! rh u
      else if(flg.eq.3) then
          call copy    (fx1, tmux, n)  ! u
          call copy    (fx2, tmry, n)  ! rh v
          call copy    (gx1, tmuy, n)  ! v
          call copy    (gx2, tmry, n)  ! rh v
          call copy    (hx1, tmuz, n)  ! w
          call copy    (hx2, tmry, n)  ! rh v
      else if(flg.eq.4) then
          call copy    (fx1, tmux, n)  ! u
          call copy    (fx2, tmrz, n)  ! rh w
          call copy    (gx1, tmuy, n)  ! v
          call copy    (gx2, tmrz, n)  ! rh w
          call copy    (hx1, tmuz, n)  ! w
          call copy    (hx2, tmrz, n)  ! rh w
      else if(flg.eq.5) then
          call copy    (fx1, tmux, n)  ! u
          call copy    (fx2, tmen, n)  ! tmen
          call add2    (fx2, tmpr, n)  ! E + p
          call copy    (gx1, tmuy, n)  ! v
          call copy    (gx2, fx2, n)   ! rh w
          call copy    (hx1, tmuz, n)  ! w
          call copy    (hx2, fx2, n)   ! rh w
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine vol_eval2(f1,g1,f2,g2,flg) ! 2
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      integer  flg
      real     f1(lt), g1(lt), f2(lt), g2(lt)
      real     tm1(lt), tm2(lt), tm3(lt)

      n=nx1*ny1*nz1*nelt
      call rzero(f2,n)
      call rzero(g2,n)
      if(flg.eq.1) then
          call copy    (f1,tmrx,n)    ! rh u
          call copy    (g1,tmry,n)    ! rh v
      else if(flg.eq.2) then
          call col3    (f1,tmrx,tmrx,n)
          call invcol2 (f1,tmrh,n)
          call add2    (f1,tmpr,n)    ! p + rh u^2 
          call col3    (g1,tmrx,tmry,n)
          call invcol2 (g1,tmrh,n)    ! rh u v
      else if(flg.eq.3) then
          call col3    (f1,tmry,tmrx,n)
          call invcol2 (f1,tmrh,n)    ! rh v u
          call col3    (g1,tmry,tmry,n)
          call invcol2 (g1,tmrh,n) 
          call add2    (g1,tmpr,n)    ! p + rh v^2 
c     else if(flg.eq.4) then
c         call col3    (f1,tmrz,tmrx,n)
c         call invcol2 (f1,tmrh,n)    ! rh w u
c         call col3    (g1,tmrz,tmry,n)
c         call invcol2 (g1,tmrh,n)    ! rh w v
c         call col3    (h1,tmrz,tmrz,n) 
c         call invcol2 (h1,tmrh,n)    
c         call add2    (h1,tmpr,n)    ! p + rh w^2 
      else if(flg.eq.5) then
          call add3    (f1,tmen,tmpr,n)
          call copy    (g1,f1,n)
          call col2    (f1,tmrx,n)
          call invcol2 (f1,tmrh,n)    ! u (E + p)
          call col2    (g1,tmry,n)
          call invcol2 (g1,tmrh,n)    ! v (E + p)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine vol_eval(f1,g1,h1,f2,g2,h2,flg) ! 3
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      integer  flg
      real     rhs(1)
      real     f1(lt), g1(lt), h1(lt), f2(lt), g2(lt), h2(lt)
      real     tm1(lt), tm2(lt), tm3(lt)

      n=nx1*ny1*nz1*nelt
      call rzero(f2,n)
      call rzero(g2,n)
      call rzero(h2,n)
      if(flg.eq.1) then
          call copy    (f1,tmrx,n)    ! rh u
          call copy    (g1,tmry,n)    ! rh v
          call copy    (h1,tmrz,n)    ! rh w
      else if(flg.eq.2) then
          call col3    (f1,tmrx,tmrx,n)
          call invcol2 (f1,tmrh,n)
          call add2    (f1,tmpr,n)    ! p + rh u^2 
          call col3    (g1,tmrx,tmry,n)
          call invcol2 (g1,tmrh,n)    ! rh u v
          call col3    (h1,tmrx,tmrz,n)
          call invcol2 (h1,tmrh,n)    ! rh u w
      else if(flg.eq.3) then
          call col3    (f1,tmry,tmrx,n)
          call invcol2 (f1,tmrh,n)    ! rh v u
          call col3    (g1,tmry,tmry,n)
          call invcol2 (g1,tmrh,n) 
          call add2    (g1,tmpr,n)    ! p + rh v^2 
          call col3    (h1,tmry,tmrz,n) 
          call invcol2 (h1,tmrh,n)    ! rh v w
      else if(flg.eq.4) then
          call col3    (f1,tmrz,tmrx,n)
          call invcol2 (f1,tmrh,n)    ! rh w u
          call col3    (g1,tmrz,tmry,n)
          call invcol2 (g1,tmrh,n)    ! rh w v
          call col3    (h1,tmrz,tmrz,n) 
          call invcol2 (h1,tmrh,n)    
          call add2    (h1,tmpr,n)    ! p + rh w^2 
      else if(flg.eq.5) then
          call add3    (f1,tmen,tmpr,n)
          call copy    (g1,f1,n)
          call copy    (h1,f1,n)
          call col2    (f1,tmrx,n)
          call invcol2 (f1,tmrh,n)    ! u (E + p)
          call col2    (g1,tmry,n)
          call invcol2 (g1,tmrh,n)    ! v (E + p)
          call col2    (h1,tmrz,n)
          call invcol2 (h1,tmrh,n)    ! w (E + p)
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine vol_res(rhs,flg)
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE' ! need logical switch ifstr
      integer  flg
      real     rhs(1)

      if(ifstr) then 
          call vol_strg(rhs, flg) ! ! strong form 
      else           
          if(ifdeal) then
              call vol_weak_d(rhs, flg) ! de-aliased 
          else
              call vol_weak(rhs, flg) ! ! weak form 
          endif
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine vol_strg(rhs,flg) ! 
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      integer  flg
      real     rhs(1)
      real     f1(lt), g1(lt), h1(lt)
     $       , f2(lt), g2(lt), h2(lt)
      real     tm1(lt), tm2(lt), tm3(lt)
     $       , f1x(lt), g1y(lt), h1z(lt)

      n=nx1*ny1*nz1*nelt
      if(if3d) then
        call vol_eval(f1,g1,h1,f2,g2,h2,flg) ! get f1, g1 ... 
          
        call gradm1(f1x,tm2,tm3,f1) ! 
        call gradm1(tm1,g1y,tm3,g1) ! 
        call gradm1(tm1,tm2,h1z,h1) ! 

        call copy(rhs,f1x,n)
        call add2(rhs,g1y,n)
        call add2(rhs,h1z,n)

        call chsign(rhs, n)  ! minus sign here 

        if(flg.gt.1) then
          call gradm1(f1x,tm2,tm3,f2)
          call gradm1(tm1,g1y,tm3,g2)
          call gradm1(tm1,tm2,h1z,h2)
          call copy(tm1,f1x,n)
          call add2(tm1,g1y,n) ! positive in strong form 
          call add2(tm1,h1z,n)

          call add2(rhs,tm1,n)
        endif
      else  ! 2D, h1 h2 should be 0 
        call vol_eval(f1,g1,h1,f2,g2,h2,flg) ! get f1, g1 ... 
          
        call gradm1(f1x,tm2,tm3,f1) ! 
        call gradm1(tm1,g1y,tm3,g1) ! 

        call copy(rhs,f1x,n)
        call add2(rhs,g1y,n)
        call chsign(rhs, n)  ! minus sign here 

        if(flg.gt.1) then
          call gradm1(f1x,tm2,tm3,f2) ! df2/dx
          call gradm1(tm1,g1y,tm3,g2) ! dg2/dy
          call copy(tm1,f1x,n)
          call add2(tm1,g1y,n) ! positive in strong form 

          call add2(rhs,tm1,n)
        endif
      endif

      return
      end
c-----------------------------------------------------------------------
c----- Evaluate volume routines 
c-----------------------------------------------------------------------
      subroutine vol_weak_d(rhs,flg) ! ! 3+2? deal 
c
c  Dealiased weak form: 1. times mass matrix 
c                       2. transpose of derivatives 
c                       3. times mass matrix inverse 
c
c     Multiplicative de-aliasing for vol routine 
c
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      parameter(ldt=lxd*lyd*lzd*lelt)
      parameter(le=lx1*ly1*lz1)
      integer  flg, e
      real     rhs(le,1), tm1(le,lelt), tm2(le,lelt)
      real     f1(le,lelt), g1(le,lelt), h1(le,lelt)
     $       , f2(le,lelt), g2(le,lelt), h2(le,lelt)
      real fx1(lt), fx2(lt), gx1(lt), gx2(lt), hx1(lt), hx2(lt)
     $   , uc (lt),  vc(lt),  wc(lt), cpr(lt), fpr(ldt) 

      n=nx1*ny1*nz1*nelt
      ne=nx1*ny1*nz1
      
c     components of fluxes , only for Euler 
c     coarse and fine pressure 

      if(if3d) then
          call vol_eval(f1,g1,h1,f2,g2,h2,flg) ! get f1, g1 ... 
          call get_fpr(fpr) ! get fpr
          call vol_eval_d(fx1,fx2,gx1,gx2,hx1,hx2
     $                   , flg) ! ! ! get f1, g1 ... 
            
          call rzero(tm2,n)

          if(flg.eq.1) then      ! no multiply needed 
              call grad_dbm1_t         (tm2,f1,g1,h1,ifdeal)
          else if(flg.eq.2) then ! 2,3,4
              call grad_dbm1_t_col2add1(tm2,fx1,fx2,fpr,gx1,gx2,hx1,hx2)
          else if(flg.eq.3) then ! 2,3,4
              call grad_dbm1_t_col2add2(tm2,fx1,fx2,gx1,gx2,fpr,hx1,hx2)
          else if(flg.eq.4) then ! 2,3,4
              call grad_dbm1_t_col2add3(tm2,fx1,fx2,gx1,gx2,hx1,hx2,fpr)
          else if(flg.eq.5) then ! multiply needed 
              call grad_dbm1_t_col2    (tm2,fx1,fx2,gx1,gx2,hx1,hx2)
          endif

      else if(ifaxis) then
          write(6,*) 'Not axial yet, quit in vol_res'
          call exitt
      else  ! first use code that's already there, will it work? 

          call vol_eval(f1,g1,h1,f2,g2,h2,flg)
          call get_fpr(fpr) ! get fpr
          call vol_eval_d(fx1,fx2,gx1,gx2,hx1,hx2
     $                   , flg) ! ! ! get f1, g1 ... for 3D though...

          call rzero(h1,n)
          call rzero(h2,n)
          call rzero(hx1,n)
          call rzero(hx2,n)

          call rzero(tm2,n)

          if(flg.eq.1) then      ! no multiply needed 
              call grad_dbm1_t         (tm2,f1,g1,h1,ifdeal)
          else if(flg.eq.2) then ! 2,3,4
              call grad_dbm1_t_col2add1(tm2,fx1,fx2,fpr,gx1,gx2,hx1,hx2)
          else if(flg.eq.3) then ! 2,3,4
              call grad_dbm1_t_col2add2(tm2,fx1,fx2,gx1,gx2,fpr,hx1,hx2)
          else if(flg.eq.5) then ! multiply needed 
              call grad_dbm1_t_col2    (tm2,fx1,fx2,gx1,gx2,hx1,hx2)
        endif

      endif

      call invbf(rhs,tm2) ! all elements 

      if(ifdbg) then
          write(6,*) 'done :: compute weak vol (deal) flg = ', flg
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine vol_weak(rhs,flg) ! ! 3+2? 
c
c  weak form: 1. times mass matrix 
c             2. transpose of derivatives 
c             3. times mass matrix inverse 
c
c     refer to routine gradm11ts in navier2.f 
c     elem based, uses w3m1, the diagnonal 'vanilla' mass matrix
c         sum to the first array ! not evaluate, but sum 
c
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      parameter(le=lx1*ly1*lz1)
      integer  flg, e
      real     rhs(le,1), tm1(le,lelt), tm2(le,lelt)
      real     f1(le,lelt), g1(le,lelt), h1(le,lelt)
     $       , f2(le,lelt), g2(le,lelt), h2(le,lelt)

      n=nx1*ny1*nz1*nelt
      ne=nx1*ny1*nz1
      
      if(if3d) then
        call vol_eval(f1,g1,h1,f2,g2,h2,flg) ! get f1, g1 ... 
        call rzero(tm2,n)
        call grad_m1_t(tm2,f1,g1,h1,ifdeal)

      else if(ifaxis) then
          write(6,*) 'Not axial yet, quit in vol_res'
          call exitt
      else  ! first use code that's already there, will it work? 
c       call vol_eval2(f1,g1,f2,g2,flg)
        call vol_eval(f1,g1,h1,f2,g2,h2,flg)
        call rzero(h1,n)
        call rzero(h2,n)
        call rzero(tm2,n)

c       call grad_m1_t2(tm2,f1,g1,ifdeal)
        call grad_m1_t(tm2,f1,g1,h1,ifdeal)
      endif

      call invbf(rhs,tm2) ! all elements 

      if(ifdbg) then
          write(6,*) 'done :: compute weak vol flg = ', flg
      endif
      return
      end
c-----------------------------------------------------------------------
c-----| Surface |-----
c-----------------------------------------------------------------------
      subroutine surf_res(flx,flg)
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      integer  flg
      real     flx(1) ! defined on surface 

      if(flxtyp.eq.1) then ! Lax Friedrichs
          if(ifdeal) then
              call flx_lf_d(flx,flg)
          else
              call flx_lf  (flx,flg)
          endif
      else if(flxtyp.eq.2) then ! Roe fluxes
          call flx_roe_3d(flx,flg) ! only 3D
      endif

      return
      end
c-----------------------------------------------------------------------
c----- Surface de-aliasing 
c-----------------------------------------------------------------------
      subroutine flx_lf_d(flx,flg) ! 3+2
c     Lax-Friedrichs flux, distinguishes between strong and weak form 
c     Add in de-aliasing for surface fluxes 
c     1. Even if weak form works strong form may still be wrong
c     b/c I never interpolated and integrated f1f terms. 
c     2. lm is also not computed using interpolation 

      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      integer  flg, i, f, e, k
      real     flx(1) ! defined on surface 

      n=nx1*ny1*nz1*nelt
      nxz=nx1*nz1
      nfaces=2*ndim

      call surf_eval  (flg) ! get f1f, g1f, h1f, cuf, no need change
c     write(6,*) 'f1f max' 
c     call printmax_fac1(f1f)  ! 2+3
c     write(6,*) 'g1f max' 
c     call printmax_fac1(g1f)  ! 2+3

      call ctr_inter_d(flg) ! eval fcf, gcf, hcf
c     write(6,*) 'fcf max' 
c     call printmax_fac1(fcf)  ! 2+3
c     write(6,*) 'gcf max' 
c     call printmax_fac1(gcf)  ! 2+3

      call dif_inter_d(flg) ! eval udf, lm 
c     write(6,*) 'udf max' 
c     call printmax_fac1(udf)  ! 2+3
c     write(6,*) 'lmd max' 
c     call printmax_fac1(lm)  ! 2+3

      k = 0
      do e=1,nelt
      do f=1,nfaces
      do i=1,nxz
         k=k+1
         if(if3d) then
           if(ifstr) then ! strong form, f- - f*
            flx(k) =  ( f1f(k) - fcf(k) )*unx(i,1,f,e)
     $              + ( g1f(k) - gcf(k) )*uny(i,1,f,e)
     $              + ( h1f(k) - hcf(k) )*unz(i,1,f,e)
     $              - lm(k)*udf(k)/2. ! lambda / 2 
c           flx(k) = flx(k)*area(i,1,f,e)
           else           ! weak form, no changes 
            flx(k) = -( fcf(k)*unx(i,1,f,e)
     $              + gcf(k)*uny(i,1,f,e)
     $              + hcf(k)*unz(i,1,f,e)
     $              + lm(k)*udf(k)/2. )
c           flx(k) = flx(k)*area(i,1,f,e)
           endif 
         else ! 2D 
           if(ifstr) then ! strong form, f- - f*
            flx(k) =  ( f1f(k) - fcf(k) )*unx(i,1,f,e)
     $              + ( g1f(k) - gcf(k) )*uny(i,1,f,e)
     $              - lm(k)*udf(k)/2. ! lambda / 2 
c           flx(k) = flx(k)*area(i,1,f,e)
           else           ! weak form, no changes 
            flx(k) = -( fcf(k)*unx(i,1,f,e)
     $                + gcf(k)*uny(i,1,f,e)
     $                + lm(k)*udf(k)/2. )
c           flx(k) = flx(k)*area(i,1,f,e)
           endif
         endif
      enddo
      enddo
      enddo

c     write(6,*) 'flx max' 
c     call printmax_fac1(flx)  ! 2+3
c     stop 

      return
      end
c-----------------------------------------------------------------------
      subroutine ctr_inter_d(flg) 
c mult. with area on fine mesh  !! 
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      parameter(lfd=lxd*lzd)              ! fine surface 
      integer  flg 
      integer  e, f, kf, fdim
      real     bf1(lf),bg1(lf),bh1(lf) ! fluxes on surfaces 

      real     jfcf(lfd,6,lelt)
     $       , jgcf(lfd,6,lelt)
     $       , jhcf(lfd,6,lelt) 
      real     zptf(lxd), wgtf(lxd), wghtc(lx1*lz1), wghtf(lxd*lzd)
      real     jaco_c(lx1*lz1),jaco_f(lxd*lzd) 

      n=nx1*ny1*nz1*nelt
      nf=nx1*nz1*2*ndim*nelt !total number of points on faces
      nfaces = 2*ndim 

      call bc_flx_eval(bf1,bg1,bh1,flg)

      call copy  (fcf,f1f,nf)         ! a 
      call add2  (fcf,bf1,nf)         ! first add bc , a- + a+(bc)
      call gs_op (dg_hndl,fcf,1,1,0)  ! 1 ==> +      , a- + a+(itr)
      call cmult (fcf,0.5,nf)         ! times 1/2, (a- + a+)/2 

      call copy  (gcf,g1f,nf)
      call add2  (gcf,bg1,nf)         ! first add bc 
      call gs_op (dg_hndl,gcf,1,1,0)  ! 1 ==> + 
      call cmult (gcf,0.5,nf)         ! times 1/2

      if(if3d) then
          call copy  (hcf,h1f,nf)
          call add2  (hcf,bh1,nf)         ! first add bc 
          call gs_op (dg_hndl,hcf,1,1,0)  ! 1 ==> + 
          call cmult (hcf,0.5,nf)         ! times 1/2, not correct bc
      endif

c     Set fine-scale surface weights
      call zwgl(zptf,wgtf,nxd)
      if (if3d) then
         k=0
         do j=1,ny1
         do i=1,nx1
            k=k+1
            wghtc(k)=wxm1(i)*wzm1(j)
         enddo
         enddo
         k=0
         do j=1,nyd
         do i=1,nxd
            k=k+1
            wghtf(k)=wgtf(i)*wgtf(j)
         enddo
         enddo
      else
         call copy(wghtc,wxm1,nx1)
         call copy(wghtf,wgtf,nxd)
      endif

c     Interpolate and sum 
      nxz  = nx1*nz1  
      nxzd = nxd*nzd  
      kf = 1 
      fdim = ndim-1 ! Dimension of face 
      do e=1,nelt 
      do f=1,nfaces

        do i=1,nxz 
          jaco_c(i) = area(i,1,f,e)/wghtc(i)
        enddo 
        call map_faced(jaco_f,jaco_c,nx1,nxd,fdim,0) ! 

        call map_faced(jfcf(1,f,e),fcf(kf),nx1,nxd,fdim,0) ! Dealiased quadrature,
        call map_faced(jgcf(1,f,e),gcf(kf),nx1,nxd,fdim,0) ! Dealiased quadrature,
        if(if3d) call map_faced(jhcf(1,f,e),hcf(kf),nx1,nxd,fdim,0) 

        do i=1,nxzd 
          jfcf(i,f,e) = wghtf(i)*jaco_f(i)*jfcf(i,f,e)
          jgcf(i,f,e) = wghtf(i)*jaco_f(i)*jgcf(i,f,e)
          if(if3d) jhcf(i,f,e) = wghtf(i)*jaco_f(i)*jhcf(i,f,e)
        enddo 

        call map_faced(fcf(kf),jfcf(1,f,e),nx1,nxd,fdim,1)  ! 1 --> uf = J^T ufine
        call map_faced(gcf(kf),jgcf(1,f,e),nx1,nxd,fdim,1)  ! 1 --> uf = J^T ufine
        if(if3d) call map_faced(hcf(kf),jhcf(1,f,e),nx1,nxd,fdim,1)  !

        kf = kf + nxz 
      enddo 
      kf = (e)*(nxz*2*ldim) + 1 
      enddo 

      return
      end
c-----------------------------------------------------------------------
      subroutine dif_inter_d(flg) ! 3+2
c mult. with area on fine mesh  !! 
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      parameter(lfd=lxd*lzd)              ! fine surface 
      integer  flg
      real     buf(lf), ufsv(lf) ! conservative variable on surfaces 

      real     judf(lfd,2*ldim,lelt) 
      integer  e, f, kf 
      real     zptf(lxd), wgtf(lxd), wghtc(lx1*lz1), wghtf(lxd*lzd)
      real     jaco_c(lx1*lz1),jaco_f(lxd*lzd) 

      n=nx1*ny1*nz1*nelt
      nf=nx1*nz1*2*ndim*nelt !total number of points on faces
      nfaces=2*ndim

c    compute lambda only for flg = 1, it's the same for later equations 
      if(flg.eq.1) then 
         call cmp_lm  ! actually can be put into pre_comp routine 
      endif 

      call rzero(buf, nf)
      call bc_cons_eval(buf,flg)

      call copy  (ufsv,cuf,nf)  ! a 
      call cmult (ufsv,2.,nf)   ! 2 a  
      call copy  (udf,cuf,nf)   ! a 
      call add2  (udf,buf,nf)   ! first add bc 
      call gs_op (dg_hndl,udf,1,1,0)  ! 1 ==> + , a + b 
      call chsign(udf,nf)       ! - ( a + b ) 
      call add2  (udf,ufsv,nf)  ! 2 a - ( a + b ) = a - b

c     Set fine-scale surface weights
      call zwgl(zptf,wgtf,nxd)
      if (if3d) then
         k=0
         do j=1,ny1
         do i=1,nx1
            k=k+1
            wghtc(k)=wxm1(i)*wzm1(j)
         enddo
         enddo
         k=0
         do j=1,nyd
         do i=1,nxd
            k=k+1
            wghtf(k)=wgtf(i)*wgtf(j)
         enddo
         enddo
      else
         call copy(wghtc,wxm1,nx1)
         call copy(wghtf,wgtf,nxd)
      endif

c     Interpolate and sum 
      nxz  = nx1*nz1  
      nxzd = nxd*nzd  
      kf = 1 
      fdim = ndim-1 ! Dimension of face 

      do e=1,nelt 
      do f=1,nfaces

        do i=1,nxz 
          jaco_c(i) = area(i,1,f,e)/wghtc(i)
        enddo 
        call map_faced(jaco_f,jaco_c,nx1,nxd,fdim,0) ! 
        call map_faced(judf(1,f,e),udf(kf),nx1,nxd,fdim,0) !
        do i=1,nxzd 
          judf(i,f,e) = wghtf(i)*jaco_f(i)*judf(i,f,e)
        enddo 
        call map_faced(udf(kf),judf(1,f,e),nx1,nxd,fdim,1)  ! 
        kf = kf + nxz 
      enddo 
      kf = (e)*(nxz*2*ldim) + 1 
      enddo 

      return
      end
c-----------------------------------------------------------------------
c---- From convect.f      - Fri Oct 16 12:36:59 CDT 2015
c-----------------------------------------------------------------------
      subroutine map_faced(ju,u,mx,md,fdim,idir) ! GLL->GL interpolation

c     GLL interpolation from mx to md for a face array of size (nx,nz)

c     If idir ^= 0, then apply transpose operator  (md to mx)

      include 'SIZE'

      real    ju(1),u(1)
      integer fdim

      parameter (ldg=lxd**3,lwkd=4*lxd*lxd)
      common /dgrad/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real jgl,jgt

      parameter (ld=2*lxd)
      common /ctmp0/ w(ld**ldim,2)

      call lim_chk(md,ld,'md   ','ld   ','map_faced ')
      call lim_chk(mx,ld,'mx   ','ld   ','map_faced ')

c     call copy(ju,u,mx)
c     return

      call get_int_ptr (i,mx,md)

      if (idir.eq.0) then
         if (fdim.eq.2) then
            call mxm(jgl(i),md,u,mx,wkd,mx)
            call mxm(wkd,md,jgt(i),mx,ju,md)
         else
            call mxm(jgl(i),md,u,mx,ju,1)
         endif
      else
         if (fdim.eq.2) then
            call mxm(jgt(i),mx,u,md,wkd,md)
            call mxm(wkd,md,jgl(i),md,ju,mx)
         else
            call mxm(jgt(i),mx,u,md,ju,1)
         endif
      endif

      return
      end
c-----------------------------------------------------------------------
c----- Surface flx routines 
c-----------------------------------------------------------------------
      subroutine flx_lf(flx,flg) ! 3+2
c     Lax-Friedrichs flux, distinguishes between strong and weak form 
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      integer  flg, i, f, e, k
      real     flx(1) ! defined on surface 

      n=nx1*ny1*nz1*nelt
      nxz=nx1*nz1
      nfaces=2*ndim

      call surf_eval(flg) ! get f1f, g1f, h1f, cuf 
c     write(6,*) 'f1f max' 
c     call printmax_fac1(f1f)  ! 2+3
c     write(6,*) 'g1f max' 
c     call printmax_fac1(g1f)  ! 2+3

      call ctr_inter(flg) ! eval fcf, gcf, hcf
c     write(6,*) 'fcf max' 
c     call printmax_fac1(fcf)  ! 2+3
c     write(6,*) 'gcf max' 
c     call printmax_fac1(gcf)  ! 2+3

      call dif_inter(flg) ! eval udf, lm 
c     write(6,*) 'udf max' 
c     call printmax_fac1(udf)  ! 2+3
c     write(6,*) 'lmd max' 
c     call printmax_fac1(lm)  ! 2+3

      k = 0
      do e=1,nelt
      do f=1,nfaces
      do i=1,nxz
         k=k+1
         if(if3d) then
           if(ifstr) then ! strong form, f- - f*
            flx(k) =  ( f1f(k) - fcf(k) )*unx(i,1,f,e)
     $              + ( g1f(k) - gcf(k) )*uny(i,1,f,e)
     $              + ( h1f(k) - hcf(k) )*unz(i,1,f,e)
     $              - lm(k)*udf(k)/2. ! lambda / 2 
            flx(k) = flx(k)*area(i,1,f,e)
           else           ! weak form, no changes 
            flx(k) = -( fcf(k)*unx(i,1,f,e)
     $              + gcf(k)*uny(i,1,f,e)
     $              + hcf(k)*unz(i,1,f,e)
     $              + lm(k)*udf(k)/2. )
            flx(k) = flx(k)*area(i,1,f,e)
           endif 
         else ! 2D 
           if(ifstr) then ! strong form, f- - f*
            flx(k) =  ( f1f(k) - fcf(k) )*unx(i,1,f,e)
     $              + ( g1f(k) - gcf(k) )*uny(i,1,f,e)
     $              - lm(k)*udf(k)/2. ! lambda / 2 
            flx(k) = flx(k)*area(i,1,f,e)
           else           ! weak form, no changes 
            flx(k) = -( fcf(k)*unx(i,1,f,e)
     $                + gcf(k)*uny(i,1,f,e)
     $                + lm(k)*udf(k)/2. )
            flx(k) = flx(k)*area(i,1,f,e)
           endif
         endif
      enddo
      enddo
      enddo

c     write(6,*) 'flx max' 
c     call printmax_fac1(flx)  ! 2+3
c     stop 
      if(ifdbg) then
          if(nid.eq.0) write(6,*) 'done :: L-F flux for flg = ', flg
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine ctr_inter(flg)
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      integer  flg 
      real     bf1(lf),bg1(lf),bh1(lf) ! fluxes on surfaces 

      n=nx1*ny1*nz1*nelt
      nf=nx1*nz1*2*ndim*nelt !total number of points on faces

      call bc_flx_eval(bf1,bg1,bh1,flg)

      call copy  (fcf,f1f,nf)         ! a 
      call add2  (fcf,bf1,nf)         ! first add bc , a- + a+(bc)
      call gs_op (dg_hndl,fcf,1,1,0)  ! 1 ==> +      , a- + a+(itr)
      call cmult (fcf,0.5,nf)         ! times 1/2, (a- + a+)/2 

      call copy  (gcf,g1f,nf)
      call add2  (gcf,bg1,nf)         ! first add bc 
      call gs_op (dg_hndl,gcf,1,1,0)  ! 1 ==> + 
      call cmult (gcf,0.5,nf)         ! times 1/2

      if(if3d) then
          call copy  (hcf,h1f,nf)
          call add2  (hcf,bh1,nf)         ! first add bc 
          call gs_op (dg_hndl,hcf,1,1,0)  ! 1 ==> + 
          call cmult (hcf,0.5,nf)         ! times 1/2, not correct bc
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine dif_inter(flg) ! 3+2
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      integer  flg
      real     buf(lf), ufsv(lf) ! conservative variable on surfaces 

      n=nx1*ny1*nz1*nelt
      nf=nx1*nz1*2*ndim*nelt !total number of points on faces

c    compute lambda only for flg = 1, it's the same for later equations 
      if(flg.eq.1) then 
         call cmp_lm  ! actually can be put into pre_comp routine 
      endif 

      call rzero(buf, nf)
      call bc_cons_eval(buf,flg)

      call copy  (ufsv,cuf,nf)  ! a 
      call cmult (ufsv,2.,nf)   ! 2 a  
      call copy  (udf,cuf,nf)   ! a 
      call add2  (udf,buf,nf)   ! first add bc 
      call gs_op (dg_hndl,udf,1,1,0)  ! 1 ==> + , a + b 
      call chsign(udf,nf)       ! - ( a + b ) 
      call add2  (udf,ufsv,nf)  ! 2 a - ( a + b ) = a - b

      return
      end
c-----------------------------------------------------------------------
      subroutine bc_cons_eval(buf,flg) ! 3+2
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      common /nekcb/ cb
      character cb*3
      integer  flg, i, f, e, k
      real     buf(lf)
      n   =nx1*ny1*nz1*nelt
      nf  =nx1*nz1*2*ndim*nelt 
      nxz =nx1*nz1
      nfaces=2*ndim
      k=0
      do e=1,nelt
      do f=1,nfaces
         ieg=lglel(e)
         cb =cbc(f,e,ifield)
         if(cb.eq.'f  ' .or. cb.eq.'F  ' .or. 
     $      cb.eq.'o  ' .or. cb.eq.'SYM' .or.  
     $                       cb.eq.'W  ' ) then
            ia=0
            call facind(kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,f)
            do iz=kz1,kz2
            do iy=ky1,ky2
            do ix=kx1,kx2
               ia = ia + 1
               k  = ia + (f-1)*nxz + (e-1)*nfaces*nxz
               if(cb.eq.'f  ') then
                   call nekasgn(ix,iy,iz,e)
                   call userbc_f(ix,iy,iz,f,e,ieg) ! rho, rux, ruy..
               else if(cb.eq.'o  ') then
                   call nekasgn(ix,iy,iz,e)
                   call userbc_o(ix,iy,iz,f,e,ieg)
               else if(cb.eq.'W  '.or.cb.eq.'SYM') then ! reflective wall
                   call nekasgn(ix,iy,iz,e)
                   call userbc_sym(ix,iy,iz,f,e,ieg)
c              else if(cb.eq.'W  ') then ! no-slip wall 
c                  call nekasgn(ix,iy,iz,e)
c                  call userbc_w(ix,iy,iz,f,e,ieg)
               endif
               if(if3d) then
                   call surf_pnt_conv(buf(k)
     $             , rho, rux, ruy, ruz, enr, pp, gama, flg)
               else 
                   call surf_pnt_conv2(buf(k)
     $             , rho, rux, ruy, enr, pp, gama, flg)
               endif
            enddo 
            enddo 
            enddo 
         endif
      enddo 
      enddo 

      return
      end
c-----------------------------------------------------------------------
      subroutine bc_flx_eval(bf1,bg1,bh1,flg) ! 3 + 2
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      common /nekcb/ cb
      character cb*3
      integer  flg, i, f, e, k
      real     bf1(lf), bg1(1), bh1(1) 
      n  = nx1*ny1*nz1*nelt
      nf = nx1*nz1*2*ndim*nelt 
      nxz= nx1*nz1
      nfaces=2*ndim
      call rzero(bf1,nf)
      call rzero(bg1,nf)
      call rzero(bh1,nf)
      k=0
      do e=1,nelt
      do f=1,nfaces
         ieg=lglel(e)
         cb =cbc(f,e,ifield)
         if(cb.eq.'f  ' .or. cb.eq.'F  ' .or. 
     $      cb.eq.'o  ' .or. cb.eq.'SYM' .or.  
     $                       cb.eq.'W  ' ) then
            ia=0
            call facind(kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,f)
            do iz=kz1,kz2
            do iy=ky1,ky2
            do ix=kx1,kx2
               ia = ia + 1
               k  = ia + (f-1)*nxz + (e-1)*nfaces*nxz
               if(cb.eq.'f  ') then  ! rho, rux, ruy..
                   call nekasgn(ix,iy,iz,e)
                   call userbc_f(ix,iy,iz,f,e,ieg) 
               else if(cb.eq.'o  ') then
                   call nekasgn(ix,iy,iz,e)
                   call userbc_o(ix,iy,iz,f,e,ieg)
               else if(cb.eq.'W  '.or.cb.eq.'SYM') then ! slip wall 
                   call nekasgn(ix,iy,iz,e)
                   call userbc_sym(ix,iy,iz,f,e,ieg)
c              else if(cb.eq.'W  ') then ! no-slip wall 
c                  call nekasgn(ix,iy,iz,e)
c                  call userbc_w(ix,iy,iz,f,e,ieg)
               endif
               if(if3d) then 
                   call surf_pnt(bf1(k),bg1(k),bh1(k)
     $             , rho, rux, ruy, ruz, enr, pp, gama, flg)
               else  
                   call surf_pnt2(bf1(k),bg1(k)
     $             , rho, rux, ruy, enr, pp, gama, flg)
               endif
            enddo 
            enddo 
            enddo 
         endif
      enddo 
      enddo 

      return
      end
c-----------------------------------------------------------------------
c----|  Roe flux  |-----
c-----------------------------------------------------------------------
      subroutine flx_roe_3d(flx,flg)
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      integer  flg
      real     flx(1) ! defined on surface 

c     For here i am only completing the 3d Roe-Pike scheme
c     Refer to Chapter 11 of Toro (Riemann solvers)

      
      nf=nx1*nz1*2*ndim*nelt 
      if(flg.eq.1) then
          call cmp_roe_3d

          call copy(flx,rflx(1,1),nf)
      else if(flg.eq.2) then
          call copy(flx,rflx(1,2),nf)
      else if(flg.eq.3) then
          call copy(flx,rflx(1,3),nf)
      else if(flg.eq.4) then
          call copy(flx,rflx(1,4),nf)
      else if(flg.eq.5) then
          call copy(flx,rflx(1,5),nf)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine cmp_roe_3d ! for 3D 
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
c
      call eval_fn_3d  ! evaluation onto surfaces, boundary values

      call roe_star_3d ! star variables, avg

      call cmp_dw_3d   ! get dw arrays 

      call cmp_rflx_3d ! get final flux term, in rflx(lf,5)
c --
      return
      end
c-----------------------------------------------------------------------
      subroutine eval_fn_3d !  !
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      real     rhf(lf), ruf(lf), rvf(lf), rwf(lf), enf(lf) 

      nf=nx1*nz1*2*ndim*nelt 

c     from solution field to face, bd + inter
      call full2face_all( rhf,  ruf,  rvf,  rwf,  enf
     $                 , tmrh, tmrx, tmry, tmrz, tmen) 

      call copy_all    ( brh, brx, bry, brz, ben
     $                ,  rhf, ruf, rvf, rwf, enf, nf) 

      call bc_cons_eval(brh, 1)
      call bc_cons_eval(brx, 2) ! modify boundary values, 
      call bc_cons_eval(bry, 3) ! nonzero interior values, 
      call bc_cons_eval(brz, 4) ! the same as rhf 
      call bc_cons_eval(ben, 5)

      call rot2nt_3d(ruf,rvf,rwf,ruf,rvf,rwf)  ! rot 2 n-t 
      call rot2nt_3d(brx,bry,brz,brx,bry,brz)  ! rot 2 n-t 

c     Get surface terms and prim. variables 
      call euler_flx_3d(fn, prv, prs
     $            , rhf, ruf, rvf, rwf, enf, gama)
c     
      call euler_flx_3d(bfn, bprv, bprs ! 
     $            , brh, brx, bry, brz, ben, gama)

      call cmp_enth(  eth,  prv,  prs)  ! 
      call cmp_enth( beth, bprv, bprs)  ! enthalpy

      return
      end
c-----------------------------------------------------------------------
      subroutine prim2cons_2( orh, orx, ory, oen
     $                     ,  irh, iux, iuy, ien, gama) !  3 
      include 'SIZE'
      include 'TOTAL'
      parameter(lf=lx1*lz1*2*ldim*lelt)
c       Get prim variables out of conservative input
c       output: prim varb term, pressure
c        input: conservative variables(surface arrays though) 
c
      real     orh(1), orx(1), ory(1), oen(1)
     $       , irh(1), iux(1), iuy(1), ien(1)
      real     gama
      real     tmp, tm1, tm2, tm3, tm5 

      nf=nx1*nz1*2*ndim*nelt 
      do i=1,nf
          orh(i)   = irh(i)  ! could have directly evaluated 
          orx(i)   = iux(i)*orh(i) 
          ory(i)   = iuy(i)*orh(i) 
          oen(i)   = ien(i) 
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine cons2prim_2( orh, oux, ouy, oen, opr
     $                     ,  irh, irx, iry, ien, gama) !  3 
      include 'SIZE'
      include 'TOTAL'
      parameter(lf=lx1*lz1*2*ldim*lelt)
c       Get prim variables out of conservative input
c       output: prim varb term, pressure
c        input: conservative variables(surface arrays though) 
c
      real     orh(1), oux(1), ouy(1), oen(1), opr(1)
     $       , irh(1), irx(1), iry(1), ien(1)
      real     gama
      real     tmp, tm1, tm2, tm3, tm5 

      nf=nx1*nz1*2*ndim*nelt 
      do i=1,nf

          call surf_pnt_conv(tm1 !  = rho 
     $       , irh(i), irx(i), iry(i), 0., ien(i), tmp, gama, 1)
c         call surf_pnt_conv(tm2 !  = rho u 
c    $       , irh(i), irx(i), iry(i), irz(i), ien(i), tmp, gama, 2)
c         call surf_pnt_conv(tm3 !  = rho v 
c    $       , irh(i), irx(i), iry(i), irz(i), ien(i), tmp, gama, 3)
c         call surf_pnt_conv(tm5 !  = en
c    $       , irh(i), irx(i), iry(i), irz(i), ien(i), tmp, gama, 5)

          opr(i)   = tmp
          orh(i)   = irh(i)  ! could have directly evaluated 
          oux(i)   = irx(i)/orh(i) 
          ouy(i)   = iry(i)/orh(i) 
          oen(i)   = ien(i) 
      enddo

      if(nid.eq.0) write(6,*) 'done :: conservative to primitive ' 
      return
      end
c-----------------------------------------------------------------------
      subroutine cons2prim( orh, oux, ouy, ouz, oen, opr
     $                   ,  irh, irx, iry, irz, ien, gama) !  3 
      include 'SIZE'
      include 'TOTAL'
      parameter(lf=lx1*lz1*2*ldim*lelt)
c       Get prim variables out of conservative input
c       output: prim varb term, pressure
c        input: conservative variables(surface arrays though) 
c
      real     orh(1), oux(1), ouy(1), ouz(1), oen(1), opr(1)
     $       , irh(1), irx(1), iry(1), irz(1), ien(1)
      real     gama
      real     tmp, tm1, tm2, tm3, tm4, tm5 

      nf=nx1*nz1*2*ndim*nelt 
      do i=1,nf

          call surf_pnt_conv(tm1 !  = rho 
     $       , irh(i), irx(i), iry(i), irz(i), ien(i), tmp, gama, 1)
          call surf_pnt_conv(tm2 !  = rho u 
     $       , irh(i), irx(i), iry(i), irz(i), ien(i), tmp, gama, 2)
          call surf_pnt_conv(tm3 !  = rho v 
     $       , irh(i), irx(i), iry(i), irz(i), ien(i), tmp, gama, 3)
          call surf_pnt_conv(tm4 !  = rho w, = 0? 
     $       , irh(i), irx(i), iry(i), irz(i), ien(i), tmp, gama, 4)
          call surf_pnt_conv(tm5 !  = en
     $       , irh(i), irx(i), iry(i), irz(i), ien(i), tmp, gama, 5)

          opr(i)   = tmp
          orh(i)   = irh(i)  ! could have directly evaluated 
          oux(i)   = tm2/orh(i) 
          ouy(i)   = tm3/orh(i) 
          ouz(i)   = tm4/orh(i) 
          oen(i)   = ien(i) 

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine euler_flx_3d( ofx, opv, opr
     $         ,  irh, irx, iry, irz, ien, gama) !  !
      include 'SIZE'
      include 'TOTAL'
      parameter(lf=lx1*lz1*2*ldim*lelt)
c       Get prim variables out of conservative input
c       output: flux term, prim varb term, pressure
c        input: conservative variables(surface arrays though) 
c
      real     ofx(lf,1), opv(lf,1), opr(1) 
     $       , irh(1), irx(1), iry(1), irz(1), ien(1)
      real     gama
      real     tmp, tm1, tm2, tm3, tm4, tm5 
      real     ofytmp, ofztmp  ! not gonna be returned 

      nf=nx1*nz1*2*ndim*nelt 

c     2 --  Does the same as 1 part 
c         call cons2prim(opv(1,1),opv(1,2),opv(1,3),opv(1,4),opv(1,5)
c    $         opr, irh, irx, iry, irz, ien, gama) 

      do i=1,nf

c     1 -- 
          call surf_pnt_conv(tm1 !  = rho 
     $       , irh(i), irx(i), iry(i), irz(i), ien(i), tmp, gama, 1)
          call surf_pnt_conv(tm2 !  = rho u 
     $       , irh(i), irx(i), iry(i), irz(i), ien(i), tmp, gama, 2)
          call surf_pnt_conv(tm3 !  = rho v 
     $       , irh(i), irx(i), iry(i), irz(i), ien(i), tmp, gama, 3)
          call surf_pnt_conv(tm4 !  = rho w, = 0? 
     $       , irh(i), irx(i), iry(i), irz(i), ien(i), tmp, gama, 4)
          call surf_pnt_conv(tm5 !  = en
     $       , irh(i), irx(i), iry(i), irz(i), ien(i), tmp, gama, 5)

          opr(i)   = tmp
          opv(i,1) = tm1  ! could have directly evaluated 
          opv(i,2) = tm2/opv(i,1) 
          opv(i,3) = tm3/opv(i,1) 
          opv(i,4) = tm4/opv(i,1) 
          opv(i,5) = tm5 
c  -- 1 -- 

          call surf_pnt     (ofx(i,1), ofytmp, ofztmp 
     $       , irh(i), irx(i), iry(i), irz(i), ien(i), tmp, gama, 1)
          call surf_pnt     (ofx(i,2), ofytmp, ofztmp 
     $       , irh(i), irx(i), iry(i), irz(i), ien(i), tmp, gama, 2)
          call surf_pnt     (ofx(i,3), ofytmp, ofztmp 
     $       , irh(i), irx(i), iry(i), irz(i), ien(i), tmp, gama, 3)
          call surf_pnt     (ofx(i,4), ofytmp, ofztmp ! ofx (, 4) ?= 0 
     $       , irh(i), irx(i), iry(i), irz(i), ien(i), tmp, gama, 4)
          call surf_pnt     (ofx(i,5), ofytmp, ofztmp 
     $       , irh(i), irx(i), iry(i), irz(i), ien(i), tmp, gama, 5)

      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine roe_star_3d !  !
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      real    smsrh(lf), mlsrh(lf)
     $      , usrh(lf) , vsrh(lf) , wsrh(lf), hsrh(lf)
      integer i,nf

      nf=nx1*nz1*2*ndim*nelt 

      call gs_roe_star_3d( mlsrh, smsrh
     $                 , usrh, vsrh, wsrh, hsrh) 

      call copy    ( rhsr, mlsrh, nf)        ! rho star 
      call invcol3 ( uxsr, usrh, smsrh, nf)  !  u  star 
      call invcol3 ( uysr, vsrh, smsrh, nf)  !  v  star 
      call invcol3 ( uzsr, wsrh, smsrh, nf)  !  v  star 
      call invcol3 ( ehsr, hsrh, smsrh, nf)  !  h  star 

      do i=1,nf ! speed of sound 
          cpd2(i) = (gama-1.)*(ehsr(i) 
     $             - 0.5*(uxsr(i)**2 + uysr(i)**2 + uzsr(i)**2))
          cpd (i) = sqrt(cpd2(i)) 
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine gs_roe_star_3d( mlsrh, smsrh
     $         , usrh, vsrh, wsrh, hsrh)  ! 
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      real sqrrh(lf), bsqrrh(lf) ! scratch boundary arrays 
     $   , smsrh(lf), bsmsrh(lf) 
     $   , mlsrh(lf), bmlsrh(lf) ! out mlsrh, smsrh, usrh, vsrh, hsrh
     $   ,  usrh(lf), busrh(lf) 
     $   ,  vsrh(lf), bvsrh(lf) 
     $   ,  wsrh(lf), bwsrh(lf) 
     $   ,  hsrh(lf), bhsrh(lf) 

      nf=nx1*nz1*2*ndim*nelt 

      call copy (  sqrrh,  prv(1,1), nf) 
      call copy ( bsqrrh, bprv(1,1), nf) 
      call vsqrt(  sqrrh, nf)  ! vector sqrt root of rho 
      call vsqrt( bsqrrh, nf)  ! vector sqrt root of rho 

      call copy (  smsrh,  sqrrh, nf)  ! for sum 
      call copy ( bsmsrh, bsqrrh, nf)  ! for sum 
      call copy (  mlsrh,  sqrrh, nf)  ! for mult
      call copy ( bmlsrh, bsqrrh, nf)  ! for mult

      call fil_inter( bsmsrh, 0.) ! fill 0 for inter
      call gs_add_bc(smsrh, bsmsrh) 

      call fil_inter( bmlsrh, 1.) ! fill 1 for inter
      call gs_mul_bc(mlsrh, bmlsrh) 

      call copy (  usrh,  sqrrh, nf)     ! u sqrt(rho) 
      call col2 (  usrh,  prv(1,2), nf)  ! prmt varb 2, u
      call copy ( busrh, bsqrrh, nf)
      call col2 ( busrh, bprv(1,2), nf) 
      call fil_inter( busrh, 0.)
      call gs_add_bc( usrh, busrh) 

      call copy ( vsrh, sqrrh, nf)      ! v sqrt(rho) 
      call col2 ( vsrh, prv(1,3), nf)   ! prmt varb 3, v
      call copy ( bvsrh, bsqrrh, nf)
      call col2 ( bvsrh, bprv(1,3), nf)  
      call fil_inter( bvsrh, 0.)        ! 
      call gs_add_bc( vsrh, bvsrh) 

      call copy ( wsrh, sqrrh, nf)      ! w sqrt(rho) 
      call col2 ( wsrh, prv(1,4), nf)   ! prmt varb 4, w
      call copy ( bwsrh, bsqrrh, nf)
      call col2 ( bwsrh, bprv(1,4), nf)  
      call fil_inter( bwsrh, 0.)        ! 
      call gs_add_bc( wsrh, bwsrh) 

      call copy ( hsrh, sqrrh, nf)      ! eh sqrt(rho) 
      call col2 ( hsrh, eth, nf) 
      call copy ( bhsrh, bsqrrh, nf)
      call col2 ( bhsrh, beth, nf) 
      call fil_inter( bhsrh, 0.)        ! 
      call gs_add_bc( hsrh, bhsrh) 

      return
      end
c-----------------------------------------------------------------------
      subroutine cmp_dw_3d !  !
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'

c     Refer to Ch 11, page 364 of Toro 
c dimension: [dw1] = [rho] * [u]
c            [dw2] = [rho] * [u]
c            [dw3] = [rho] * [u] * [u] , 1 higher 
c            [dw5] = [rho] * [u]
      nf=nx1*nz1*2*ndim*nelt 

      call gs_roe_dif_3d  ! get difference in rho, u, v, w, p

      do i=1,nf
          dw1(i) = 0.5*rhsr(i)*difu(i)/cpd(i) - 0.5*difp(i)/cpd2(i)
          dw1(i) = dw1(i) * abs(uxsr(i) - cpd(i))

          dw2(i) = -difr(i) + difp(i)/cpd2(i)
          dw2(i) = dw2(i) * abs(uxsr(i))

          dw3(i) = -rhsr(i)*difv(i)
          dw3(i) = dw3(i) * abs(uxsr(i))

          dw4(i) = -rhsr(i)*difw(i)
          dw4(i) = dw4(i) * abs(uxsr(i))

          dw5(i) = -0.5*rhsr(i)*difu(i)/cpd(i) - 0.5*difp(i)/cpd2(i) 
          dw5(i) = dw5(i) * abs(uxsr(i) + cpd(i))
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine gs_roe_dif_3d !  !
c     Get u- - u+
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      real     rsv(lf), usv(lf), vsv(lf), wsv(lf), psv(lf)
     $      , bdfr(lf), bdfu(lf), bdfv(lf), bdfw(lf), bdfp(lf) 

      nf=nx1*nz1*2*ndim*nelt 

      call copy(rsv, prv(1,1),nf)
      call copy(usv, prv(1,2),nf)
      call copy(vsv, prv(1,3),nf)
      call copy(wsv, prv(1,4),nf)
      call copy(psv, prs,nf)

      call copy(difr, prv(1,1),nf)   ! face 
      call copy(difu, prv(1,2),nf) 
      call copy(difv, prv(1,3),nf) 
      call copy(difw, prv(1,4),nf) 
      call copy(difp, prs,nf) 

      call copy(bdfr, bprv(1,1),nf)  ! bdry 
      call copy(bdfu, bprv(1,2),nf) 
      call copy(bdfv, bprv(1,3),nf) 
      call copy(bdfw, bprv(1,4),nf) 
      call copy(bdfp, bprs,nf) 

      call fil_inter(bdfr, 0.)
      call gs_add_bc(difr, bdfr)  ! a + b 

      call fil_inter(bdfu, 0.)
      call gs_add_bc(difu, bdfu) 

      call fil_inter(bdfv, 0.)
      call gs_add_bc(difv, bdfv) 

      call fil_inter(bdfw, 0.)
      call gs_add_bc(difw, bdfw) 

      call fil_inter(bdfp, 0.)
      call gs_add_bc(difp, bdfp) 

      call chsign_5(difr, difu, difv, difw, difp, nf) 

      call cmult   (rsv, 2., nf) 
      call cmult   (usv, 2., nf) 
      call cmult   (vsv, 2., nf) 
      call cmult   (wsv, 2., nf) 
      call cmult   (psv, 2., nf) 

      call add2    (difr, rsv, nf)  ! 2 a - ( a + b )
      call add2    (difu, usv, nf)  ! = a - b
      call add2    (difv, vsv, nf) 
      call add2    (difw, wsv, nf) 
      call add2    (difp, psv, nf) 

      return
      end
c-----------------------------------------------------------------------
      subroutine cmp_rflx_3d ! !
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      real     fnsv(lf,5) 

      nf=nx1*nz1*2*ndim*nelt 

      if(ifstr) then
          call copy_all (
     $   fnsv(1,1), fnsv(1,2), fnsv(1,3), fnsv(1,4), fnsv(1,5)
     $  ,  fn(1,1),   fn(1,2),   fn(1,3),   fn(1,4),   fn(1,5), nf)
      endif

c sum and divide by 2, central flux 

      call gs_add_bc_all(fn, bfn) 
      call cmult_all    (fn, 0.5, nf)  ! (f- + f+) / 2 

c use dw_ values 
      do i=1,nf
          rflx(i,1) = fn(i,1) 
     $      - (dw1(i)*1. + dw2(i)*1. + dw3(i)*0. + dw4(i)*0.
     $                   + dw5(i)*1. )/2.
          rflx(i,2) = fn(i,2) 
     $      - ( dw1(i)*(uxsr(i) - cpd(i)) 
     $        + dw2(i)*uxsr(i) + dw3(i)*0. + dw4(i)*0. 
     $        + dw5(i)*(uxsr(i) + cpd(i)) )/2.
          rflx(i,3) = fn(i,3) 
     $      - ( dw1(i)*uysr(i) 
     $        + dw2(i)*uysr(i) + dw3(i)*1. + dw4(i)*0.
     $        + dw5(i)*uysr(i) )/2.
          rflx(i,4) = fn(i,4) 
     $      - ( dw1(i)*uzsr(i) + dw2(i)*uzsr(i)
     $        + dw3(i)*0.  + dw4(i)*1.
     $        + dw5(i)*uzsr(i) )/2.
          rflx(i,5) = fn(i,5) 
     $      - ( dw1(i)*(ehsr(i) - uxsr(i)*cpd(i))
     $        + dw2(i)*(uxsr(i)**2 + uysr(i)**2 + uzsr(i)**2)/2. 
     $        + dw3(i)*uysr(i) + dw4(i)*uzsr(i)
     $        + dw5(i)*(ehsr(i) + uxsr(i)*cpd(i)) )/2.
      enddo

      call chsign_all(rflx, nf)  ! complete flux here

      if(ifstr) then
          call add2_all  (rflx, fnsv, nf)  ! f- - f* 
      endif

c from n-t to x-y 
      call rot2xy_3d(rflx(1,2),rflx(1,3),rflx(1,4)
     $              ,rflx(1,2),rflx(1,3),rflx(1,4))
c     
      call col_area_all(rflx) ! times area 

      return
      end
c-----------------------------------------------------------------------
c----- Service routines for Roe flux 
c-----------------------------------------------------------------------
      subroutine cmp_enth( enth, prmv, prsu)  ! 
c     Enthalpy = ( Energy + Pressure ) / Rho
      include 'SIZE'
      include 'TOTAL'
      parameter(lf=lx1*lz1*2*ldim*lelt)
      real     enth(1), prmv(lf,1), prsu(1) 
      integer  i, nf

      nf=nx1*nz1*2*ndim*nelt 
      do i=1,nf
          enth(i) = ( prmv(i,5) + prsu(i)) / prmv(i,1) 
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine col_area_all(flx)
c     
      include 'SIZE'
      include 'TOTAL'
      parameter(lf=lx1*lz1*2*ldim*lelt)
      integer  i, f, e, k
      real     flx(lf,1) ! defined on surface 

      nxz=nx1*nz1
      nfaces=2*ndim

      k = 0
      do e=1,nelt
      do f=1,nfaces
      do i=1,nxz
         k=k+1
         flx(k,1) = flx(k,1)*area(i,1,f,e)
         flx(k,2) = flx(k,2)*area(i,1,f,e)
         flx(k,3) = flx(k,3)*area(i,1,f,e)
         flx(k,4) = flx(k,4)*area(i,1,f,e)
         flx(k,5) = flx(k,5)*area(i,1,f,e)
      enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine rot2nt_3d(rn,rt,rt2,arx,ary,arz) ! from x-y to n-t 
      include 'SIZE'
      include 'TOTAL'
      real    tm1, tm2 , tm3
      integer i,f,e,k 
      real    rn(1), rt(1), rt2(1), arx(1), ary(1), arz(1)

      nfaces=2*ndim
      nxz = nx1*nz1
      k=0
      do e=1,nelt
      do f=1,nfaces
      do i=1,nxz
          k = k + 1
          tm1 =   arx(k)*unx(i,1,f,e) + ary(k)*uny(i,1,f,e)
     $          + arz(k)*unz(i,1,f,e)
          tm2 =   arx(k)*t1x(i,1,f,e) + ary(k)*t1y(i,1,f,e)
     $          + arz(k)*t1z(i,1,f,e)
          tm3 =   arx(k)*t2x(i,1,f,e) + ary(k)*t2y(i,1,f,e)
     $          + arz(k)*t2z(i,1,f,e)
          rn (k) = tm1 
          rt (k) = tm2  ! should be able to deal with arx = rn
          rt2(k) = tm3 
      enddo 
      enddo 
      enddo 

      return
      end
c-----------------------------------------------------------------------
      subroutine rot2xy_3d(arx,ary,arz,rn,rt,rt2) ! from n-t to x-y 
      include 'SIZE'
      include 'TOTAL'
c     Orthogonal matrix, inverse equals to transpose
      real    tm1, tm2 , tm3
      integer i,f,e,k 
      real    rn(1), rt(1), rt2(1), arx(1), ary(1), arz(1)

      nfaces=2*ndim
      nxz = nx1*nz1
      k=0
      do e=1,nelt
      do f=1,nfaces
      do i=1,nxz
          k = k + 1
          tm1 =   rn (k)*unx(i,1,f,e) + rt (k)*t1x(i,1,f,e)
     $          + rt2(k)*t2x(i,1,f,e)
          tm2 =   rn (k)*uny(i,1,f,e) + rt (k)*t1y(i,1,f,e)
     $          + rt2(k)*t2y(i,1,f,e)
          tm3 =   rn (k)*unz(i,1,f,e) + rt (k)*t1z(i,1,f,e)
     $          + rt2(k)*t2z(i,1,f,e)
          arx(k) = tm1 
          ary(k) = tm2  ! should be able to deal with arx = rn
          arz(k) = tm3 
      enddo 
      enddo 
      enddo 

      return
      end
c-----------------------------------------------------------------------
c----- Service routines 2 
c-----------------------------------------------------------------------
      subroutine full2face_all2( rhf,  ruf,  rvf,  enf
     $                        , tmrh, tmrx, tmry, tmen) 
      include 'SIZE'
      include 'TOTAL' 
      real rhf(1),  ruf(1),  rvf(1),  enf(1)
     $  , tmrh(1), tmrx(1), tmry(1), tmen(1)

      call full2face(rhf, tmrh) 
      call full2face(ruf, tmrx) 
      call full2face(rvf, tmry) 
      call full2face(enf, tmen) 

      return
      end
c-----------------------------------------------------------------------
      subroutine full2face_all( rhf,  ruf,  rvf,  rwf,  enf
     $                       , tmrh, tmrx, tmry, tmrz, tmen) 
      include 'SIZE'
      include 'TOTAL' 
      real rhf(1),  ruf(1),  rvf(1),  rwf(1),  enf(1)
     $  , tmrh(1), tmrx(1), tmry(1), tmrz(1), tmen(1)

      call full2face(rhf, tmrh) 
      call full2face(ruf, tmrx) 
      call full2face(rvf, tmry) 
      call full2face(rwf, tmrz) 
      call full2face(enf, tmen) 

      return
      end
c-----------------------------------------------------------------------
      subroutine copy_all2(fo1,fo2,fo3,fo4
     $                   , fl1,fl2,fl3,fl4,n) 
      include 'SIZE'
      include 'TOTAL' 
      real    fo1(1),fo2(1),fo3(1),fo4(1)
     $     ,  fl1(1),fl2(1),fl3(1),fl4(1)
      integer n 

      call copy (fo1, fl1, n) 
      call copy (fo2, fl2, n) 
      call copy (fo3, fl3, n) 
      call copy (fo4, fl4, n) 

      return
      end
c-----------------------------------------------------------------------
      subroutine copy_all(fo1,fo2,fo3,fo4,fo5
     $                  , fl1,fl2,fl3,fl4,fl5,n) 
      include 'SIZE'
      include 'TOTAL' 
      real    fo1(1),fo2(1),fo3(1),fo4(1),fo5(1)
     $     ,  fl1(1),fl2(1),fl3(1),fl4(1),fl5(1) 
      integer n 

      call copy (fo1, fl1, n) 
      call copy (fo2, fl2, n) 
      call copy (fo3, fl3, n) 
      call copy (fo4, fl4, n) 
      call copy (fo5, fl5, n) 

      return
      end
c-----------------------------------------------------------------------
      subroutine chsign_5(fl1,fl2,fl3,fl4,fl5,n) 
c     include 'SIZE'
c     include 'TOTAL' 
      real    fl1(1), fl2(1), fl3(1), fl4(1), fl5(1)
      integer n  ! should be nf 

      call chsign (fl1, n) 
      call chsign (fl2, n) 
      call chsign (fl3, n) 
      call chsign (fl4, n) 
      call chsign (fl5, n) 

      return
      end
c-----------------------------------------------------------------------
      subroutine chsign_all(fl, n) 
      include 'SIZE'
      include 'TOTAL' 
      parameter(lf=lx1*lz1*2*ldim*lelt)
      real    fl(lf,1) 
      integer n  ! should be nf 

      call chsign (fl(1,1), n) 
      call chsign (fl(1,2), n) 
      call chsign (fl(1,3), n) 
      call chsign (fl(1,4), n) 
      call chsign (fl(1,5), n) 

      return
      end
c-----------------------------------------------------------------------
      subroutine add2_all(fl, bfl, n) 
      include 'SIZE'
      include 'TOTAL' 
      parameter(lf=lx1*lz1*2*ldim*lelt)
      real    fl(lf,1) , bfl(lf,1) 
      integer n  ! should be nf 

      call add2 (fl(1,1), bfl(1,1), n) 
      call add2 (fl(1,2), bfl(1,2), n) 
      call add2 (fl(1,3), bfl(1,3), n) 
      call add2 (fl(1,4), bfl(1,4), n) 
      call add2 (fl(1,5), bfl(1,5), n) 

      return
      end
c-----------------------------------------------------------------------
      subroutine cmult_all2(fl1,fl2,fl3,fl4, c, n) 
      include 'SIZE'
      include 'TOTAL' 
      parameter(lf=lx1*lz1*2*ldim*lelt)
      real    fl1(1), fl2(1),fl3(1),fl4(1)
      real    c
      integer n 

      call cmult (fl1, c, n) 
      call cmult (fl2, c, n) 
      call cmult (fl3, c, n) 
      call cmult (fl4, c, n) 

      return
      end
c-----------------------------------------------------------------------
c----- Gather scatter service routines 
c-----------------------------------------------------------------------
      subroutine gs_add_bc_all2(flx1, flx2, flx3, flx4
     $                        , bfx1, bfx2, bfx3, bfx4) 
      include 'SIZE'
      include 'TOTAL' 
      parameter(lf=lx1*lz1*2*ldim*lelt)
      real    flx1(1), bfx1(1)
     $      , flx2(1), bfx2(1)
     $      , flx3(1), bfx3(1)
     $      , flx4(1), bfx4(1)
      real    bftm1(lf)

      nf=nx1*nz1*2*ndim*nelt !total number of points on faces

      call copy     ( bftm1, bfx1, nf)  ! copy into tmp 
      call fil_inter( bftm1, 0.)        ! zero out inter
      call gs_add_bc( flx1, bftm1) 

      call copy     ( bftm1, bfx2, nf)  ! copy into tmp 
      call fil_inter( bftm1, 0.)        ! zero out inter
      call gs_add_bc( flx2, bftm1) 

      call copy     ( bftm1, bfx3, nf)  ! copy into tmp 
      call fil_inter( bftm1, 0.)        ! zero out inter
      call gs_add_bc( flx3, bftm1) 

      call copy     ( bftm1, bfx4, nf)  ! copy into tmp 
      call fil_inter( bftm1, 0.)        ! zero out inter
      call gs_add_bc( flx4, bftm1) 

      return
      end
c-----------------------------------------------------------------------
      subroutine cmult_all(fl, c, n)  ! does not belong here 
      include 'SIZE'
      include 'TOTAL' 
      parameter(lf=lx1*lz1*2*ldim*lelt)
      real    fl(lf,1)
      real    c
      integer n 

      call cmult (fl(1,1), c, n) 
      call cmult (fl(1,2), c, n) 
      call cmult (fl(1,3), c, n) 
      call cmult (fl(1,4), c, n) 
      call cmult (fl(1,5), c, n) 

      return
      end
c-----------------------------------------------------------------------
      subroutine gs_add_bc_all(flxa, bfxa) 
      include 'SIZE'
      include 'TOTAL' 
      parameter(lf=lx1*lz1*2*ldim*lelt)
      real    flxa(lf,1), bfxa(lf,1)
      real    bftm1(lf), bftm2(lf), bftm3(lf)
     $      , bftm4(lf), bftm5(lf) 

      nf=nx1*nz1*2*ndim*nelt !total number of points on faces

      call copy     ( bftm1, bfxa(1,1), nf)  ! copy into tmp 
      call fil_inter( bftm1, 0.)             ! zero out inter
      call gs_add_bc(flxa(1,1), bftm1) 

      call copy     ( bftm2, bfxa(1,2), nf)  ! copy into tmp 
      call fil_inter( bftm2, 0.)             ! zero out inter
      call gs_add_bc(flxa(1,2), bftm2) 

      call copy     ( bftm3, bfxa(1,3), nf)  ! copy into tmp 
      call fil_inter( bftm3, 0.)             ! zero out inter
      call gs_add_bc(flxa(1,3), bftm3) 

      call copy     ( bftm4, bfxa(1,4), nf)  ! copy into tmp 
      call fil_inter( bftm4, 0.)             ! zero out inter
      call gs_add_bc(flxa(1,4), bftm4) 

      call copy     ( bftm5, bfxa(1,5), nf)  ! copy into tmp 
      call fil_inter( bftm5, 0.)             ! zero out inter
      call gs_add_bc(flxa(1,5), bftm5) 

      return
      end
c-----------------------------------------------------------------------
      subroutine gs_add_bc(flx, bflx) 
      include 'SIZE'
      include 'TOTAL' 
      include 'DG'
      real    flx(1), bflx(1)

      nf=nx1*nz1*2*ndim*nelt !total number of points on faces

      call gs_op(dg_hndl,flx,1,1,0) ! 1 ==> + , a + b , inter 
      call add2 ( flx, bflx, nf)     ! add boundary values 

      return
      end
c-----------------------------------------------------------------------
      subroutine gs_mul_bc(flx, bflx) 
      include 'SIZE'
      include 'TOTAL' 
      include 'DG'
      real    flx(1), bflx(1)

      nf=nx1*nz1*2*ndim*nelt !total number of points on faces

      call gs_op(dg_hndl,flx,1,2,0) ! 2 ==> * , a * b , inter 
      call col2 ( flx, bflx, nf)     ! mul boundary values 

      return
      end
c-----------------------------------------------------------------------
c----- Service 
c-----------------------------------------------------------------------
      subroutine fil_inter( bf, cons) ! fill 0 for inter
      include 'SIZE'
      include 'TOTAL'
      parameter(lf=lx1*lz1*2*ldim*lelt)
      common /nekcb/ cb
      character cb*3
      integer  flg, i, f, e, k
      real     bf(lf)
      n   =nx1*ny1*nz1*nelt
      nf  =nx1*nz1*2*ndim*nelt 
      nxz =nx1*nz1
      nfaces=2*ndim
      k=0
      do e=1,nelt
      do f=1,nfaces
         ieg=lglel(e)
         cb =cbc(f,e,ifield)
         if(cb.eq.'e  ' .or. cb.eq.'E  ') then ! inter 
            ia=0
            call facind(kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,f)
            do iz=kz1,kz2
            do iy=ky1,ky2
            do ix=kx1,kx2
               ia = ia + 1
               k  = ia + (f-1)*nxz + (e-1)*nfaces*nxz
               bf(k) = cons
            enddo 
            enddo 
            enddo 
         endif
      enddo 
      enddo 

      return
      end
c-----------------------------------------------------------------------
c-----| Boundary conditions 
c-----------------------------------------------------------------------
      subroutine userbc_o(ix,iy,iz,iside,e,eg)
      include 'SIZE'
      include 'TSTEP'
      include 'NEKUSE'
      include 'DGUSE'
      integer e,eg

      rho = tmrh(ix,iy,iz,e)  ! copy interior values 
      rux = tmrx(ix,iy,iz,e)
      ruy = tmry(ix,iy,iz,e)
      ruz = tmrz(ix,iy,iz,e)
      enr = tmen(ix,iy,iz,e)

      return
      end
c-----------------------------------------------------------------------
      subroutine userbc_sym(ix,iy,iz,iside,e,eg)
      include 'SIZE'
      include 'TSTEP'
      include 'GEOM'
      include 'NEKUSE'
      include 'DGUSE'
      include 'INPUT' ! if3d is in INPUT 
      integer e, eg, i 
      real    tm1, tm2, tm3 

c page 234, reflective boundary , 2D 

      rho = tmrh(ix,iy,iz,e)  ! was wrong because e was using f's value
      enr = tmen(ix,iy,iz,e)

      tm1 = tmrx(ix,iy,iz,e) 
      tm2 = tmry(ix,iy,iz,e) 
      tm3 = tmrz(ix,iy,iz,e) 

      if(iside.le.4) then ! 3+2 ? 
          i = mod(iside,2)*ix + (1-mod(iside,2))*iy + (iz-1)*nx1 
      else ! iside = 5, 6
          i = ix + (iy - 1)* nx1
      endif

      if(if3d) then 
          rux = tm1 - 2.*( tm1*unx(i,1,iside,e) + tm2*uny(i,1,iside,e) 
     $                   + tm3*unz(i,1,iside,e) ) * unx(i,1,iside,e) 
          ruy = tm2 - 2.*( tm1*unx(i,1,iside,e) + tm2*uny(i,1,iside,e) 
     $                   + tm3*unz(i,1,iside,e) ) * uny(i,1,iside,e) 
          ruz = tm3 - 2.*( tm1*unx(i,1,iside,e) + tm2*uny(i,1,iside,e) 
     $                   + tm3*unz(i,1,iside,e) ) * unz(i,1,iside,e) 
      else
          rux = tm1 - 2.*( tm1*unx(i,1,iside,e) + tm2*uny(i,1,iside,e) 
     $                    ) * unx(i,1,iside,e) 
          ruy = tm2 - 2.*( tm1*unx(i,1,iside,e) + tm2*uny(i,1,iside,e) 
     $                    ) * uny(i,1,iside,e) 
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine userbc_w_vis(ix,iy,iz,iside,e,eg)
      include 'SIZE'
      include 'TSTEP'
      include 'NEKUSE'
      include 'DGUSE'
      integer e,eg
c     Trying to do : Bassi Rebay 1997 
c  This routine is for imposing bc on the viscous flux 

      write(6,*) 'quit in routine bc_w_vis'
      call exitt 
      stop 
c     rho = tmrh(ix,iy,iz,e)  !  rh+ = rh- 
c     rux = - tmrx(ix,iy,iz,e)
c     ruy = - tmry(ix,iy,iz,e)
c     ruz = - tmrz(ix,iy,iz,e)
c     enr = tmen(ix,iy,iz,e)

      return
      end
c-----------------------------------------------------------------------
      subroutine userbc_w(ix,iy,iz,iside,e,eg)
      include 'SIZE'
      include 'TSTEP'
      include 'NEKUSE'
      include 'DGUSE'
      integer e,eg
c     Trying to do : Bassi Rebay 1997 
c   Problem: isentropic process, does E (total energy) change? 
c   Known: stagnation temperature should not change, given no 
c     heat coming in and no work. 

      rho = tmrh(ix,iy,iz,e)  !  rh+ = rh- 
      rux = - tmrx(ix,iy,iz,e)
      ruy = - tmry(ix,iy,iz,e)
      ruz = - tmrz(ix,iy,iz,e)
      enr = tmen(ix,iy,iz,e)

      return
      end
c-----------------------------------------------------------------------
      subroutine surf_pnt_conv2(cu 
     $           , rho, rux, ruy, en, pp, gama, flg)
c     Evaluate conservative variables to cu 
c     No common blocks 
c     -  Input  :  rho, rux, ruy, ruz, en, flg
c     -  Output :  cu, pp
      integer  flg
      real     cu
      real     rho, rux, ruy, ruz, en, pp, gama 

      pp = (gama-1.)*(en - 0.5*(rux*rux+ruy*ruy+ruz*ruz)/rho) 
      if(flg.eq.1) then
          cu = rho   ! rho 
      else if(flg.eq.2) then
          cu = rux   ! rh u
      else if(flg.eq.3) then
          cu = ruy   ! rh v
c     else if(flg.eq.4) then
c         cu = ruz   ! rh w
      else if(flg.eq.5) then
          cu = en    ! en  
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine surf_pnt2(f1,g1
     $           , rho, rux, ruy, en, pp, gama, flg)
c     Evaluate surface flux at point 
c     No common blocks 
c     -  Input  :  rho, rux, ruy, ruz, en, flg
c     -  Output :  f1, g1, h1, pp
      integer  flg
      real     f1, g1
      real     rho, rux, ruy, en, pp, gama 

      pp = (gama-1.)*(en - 0.5*(rux*rux+ruy*ruy)/rho) 
      if(flg.eq.1) then
          f1 = rux   ! rh u
          g1 = ruy   ! rh v
      else if(flg.eq.2) then
          f1 = rux*rux/rho + pp   ! rh u ^2 + p
          g1 = ruy*rux/rho        ! rh v u
      else if(flg.eq.3) then
          f1 = rux*ruy/rho        ! rh u v
          g1 = ruy*ruy/rho + pp   ! rh v ^2 + p
c     else if(flg.eq.4) then
c         f1 = rux*ruz/rho        ! rh u w
c         g1 = ruy*ruz/rho        ! rh v w
c         h1 = ruz*ruz/rho + pp   ! rh w ^2 + p
      else if(flg.eq.5) then
          f1 = rux*(en + pp)/rho  ! u ( E + p) 
          g1 = ruy*(en + pp)/rho  ! v ( E + p) 
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine surf_pnt_conv(cu 
     $           , rho, rux, ruy, ruz, en, pp, gama, flg)
c     Evaluate conservative variables to cu 
c     No common blocks 
c     -  Input  :  rho, rux, ruy, ruz, en, flg
c     -  Output :  cu, pp
      integer  flg
      real     cu
      real     rho, rux, ruy, ruz, en, pp, gama 

      pp = (gama-1.)*(en - 0.5*(rux*rux+ruy*ruy+ruz*ruz)/rho) 
      if(flg.eq.1) then
          cu = rho   ! rho 
      else if(flg.eq.2) then
          cu = rux   ! rh u
      else if(flg.eq.3) then
          cu = ruy   ! rh v
      else if(flg.eq.4) then
          cu = ruz   ! rh w
      else if(flg.eq.5) then
          cu = en    ! en  
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine surf_pnt(f1,g1,h1
     $           , rho, rux, ruy, ruz, en, pp, gama, flg)
c     Evaluate surface flux at point 
c     No common blocks 
c     -  Input  :  rho, rux, ruy, ruz, en, flg
c     -  Output :  f1, g1, h1, pp
      integer  flg
      real     f1, g1, h1
      real     rho, rux, ruy, ruz, en, pp, gama 

      pp = (gama-1.)*(en - 0.5*(rux*rux+ruy*ruy+ruz*ruz)/rho) 
      if(flg.eq.1) then
          f1 = rux   ! rh u
          g1 = ruy   ! rh v
          h1 = ruz   ! rh w
      else if(flg.eq.2) then
          f1 = rux*rux/rho + pp   ! rh u ^2 + p
          g1 = ruy*rux/rho        ! rh v u
          h1 = ruz*rux/rho        ! rh w u
      else if(flg.eq.3) then
          f1 = rux*ruy/rho        ! rh u v
          g1 = ruy*ruy/rho + pp   ! rh v ^2 + p
          h1 = ruz*ruy/rho        ! rh w v
      else if(flg.eq.4) then
          f1 = rux*ruz/rho        ! rh u w
          g1 = ruy*ruz/rho        ! rh v w
          h1 = ruz*ruz/rho + pp   ! rh w ^2 + p
      else if(flg.eq.5) then
          f1 = rux*(en + pp)/rho  ! u ( E + p) 
          g1 = ruy*(en + pp)/rho  ! v ( E + p) 
          h1 = ruz*(en + pp)/rho  ! w ( E + p) 
      endif
      return
      end
c-----------------------------------------------------------------------
      subroutine cmp_lm  ! 3+2
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      parameter(lu = lx1*ly1*lz1*lelt) 
      parameter(le = lx1*ly1*lz1)
      real    lmv(lu)
      integer e 

      ne = nx1*ny1*nz1
      n  = nx1*ny1*nz1*nelt

c     Get lambda = abs(rhou/ rho) + sqrt( gamma* pres./rho) 
      do e=1,nelt
      do i=1,le
       if(if3d) then
         lmv(i+(e-1)*le) = sqrt(gama*tmpr(i,1,1,e)/tmrh(i,1,1,e)) 
     $  + sqrt(tmux(i,1,1,e)**2 + tmuy(i,1,1,e)**2 + tmuz(i,1,1,e)**2) 
       else 
         lmv(i+(e-1)*le) = sqrt(gama*tmpr(i,1,1,e)/tmrh(i,1,1,e)) 
     $  + sqrt(tmux(i,1,1,e)**2 + tmuy(i,1,1,e)**2) 
       endif
      enddo 
      enddo
      call full2face(lm,lmv)         ! lm always > 0 
      call gs_op (dg_hndl,lm,1,4,0)  ! 4 ==> max, same on interfaces
c     ---
c     this can find the maximum at element interface, but totally 
c     ignores boundary conditions, is it problematic ?
c     I would assume it matters but not a lot; still try to fix it with
c     correct procedure ( get correct boundary conditions ) 
c     Mon Aug 10 20:52:45 CDT 2015
c     ---

      return
      end
c-----------------------------------------------------------------------
      subroutine surf_eval2(flg)
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      integer  flg
      real     f1(lt), g1(lt), u1(lt)

      n=nx1*ny1*nz1*nelt
      if(flg.eq.1) then
          call copy    (u1,tmrh,n)    ! rho
          call copy    (f1,tmrx,n)    ! rh u
          call copy    (g1,tmry,n)    ! rh v
      else if(flg.eq.2) then
          call copy    (u1,tmrx,n)    ! rho u
          call col3    (f1,tmrx,tmrx,n)
          call invcol2 (f1,tmrh,n)
          call add2    (f1,tmpr,n)    ! p + rh u^2 
          call col3    (g1,tmrx,tmry,n)
          call invcol2 (g1,tmrh,n)    ! rh u v
      else if(flg.eq.3) then
          call copy    (u1,tmry,n)    ! rho v
          call col3    (f1,tmry,tmrx,n)
          call invcol2 (f1,tmrh,n)    ! rh v u
          call col3    (g1,tmry,tmry,n)
          call invcol2 (g1,tmrh,n) 
          call add2    (g1,tmpr,n)    ! p + rh v^2 
c     else if(flg.eq.4) then
c         call copy    (u1,tmrz,n)    ! rho w
c         call col3    (f1,tmrz,tmrx,n)
c         call invcol2 (f1,tmrh,n)    ! rh w u
c         call col3    (g1,tmrz,tmry,n)
c         call invcol2 (g1,tmrh,n)    ! rh w v
c         call col3    (h1,tmrz,tmrz,n) 
c         call invcol2 (h1,tmrh,n)    
c         call add2    (h1,tmpr,n)    ! p + rh w^2 
      else if(flg.eq.5) then
          call copy    (u1,tmen,n)    ! energy
          call add3    (f1,tmen,tmpr,n)
          call copy    (g1,f1,n)
          call copy    (h1,f1,n)      ! E + p
          call col2    (f1,tmrx,n)
          call invcol2 (f1,tmrh,n)    ! u (E + p)
          call col2    (g1,tmry,n)
          call invcol2 (g1,tmrh,n)    ! v (E + p)
      endif
      call full2face(f1f,f1)
      call full2face(g1f,g1)
      call full2face(cuf,u1)
      return
      end
c-----------------------------------------------------------------------
      subroutine surf_eval(flg)
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      integer  flg
      real     f1(lt), g1(lt), h1(lt), u1(lt)

      n=nx1*ny1*nz1*nelt
      if(flg.eq.1) then
          call copy    (u1,tmrh,n)    ! rho
          call copy    (f1,tmrx,n)    ! rh u
          call copy    (g1,tmry,n)    ! rh v
          call copy    (h1,tmrz,n)    ! rh w
      else if(flg.eq.2) then
          call copy    (u1,tmrx,n)    ! rho u
          call col3    (f1,tmrx,tmrx,n)
          call invcol2 (f1,tmrh,n)
          call add2    (f1,tmpr,n)    ! p + rh u^2 
          call col3    (g1,tmrx,tmry,n)
          call invcol2 (g1,tmrh,n)    ! rh u v
          call col3    (h1,tmrx,tmrz,n)
          call invcol2 (h1,tmrh,n)    ! rh u w
      else if(flg.eq.3) then
          call copy    (u1,tmry,n)    ! rho v
          call col3    (f1,tmry,tmrx,n)
          call invcol2 (f1,tmrh,n)    ! rh v u
          call col3    (g1,tmry,tmry,n)
          call invcol2 (g1,tmrh,n) 
          call add2    (g1,tmpr,n)    ! p + rh v^2 
          call col3    (h1,tmry,tmrz,n) 
          call invcol2 (h1,tmrh,n)    ! rh v w
      else if(flg.eq.4) then
          call copy    (u1,tmrz,n)    ! rho w
          call col3    (f1,tmrz,tmrx,n)
          call invcol2 (f1,tmrh,n)    ! rh w u
          call col3    (g1,tmrz,tmry,n)
          call invcol2 (g1,tmrh,n)    ! rh w v
          call col3    (h1,tmrz,tmrz,n) 
          call invcol2 (h1,tmrh,n)    
          call add2    (h1,tmpr,n)    ! p + rh w^2 
      else if(flg.eq.5) then
          call copy    (u1,tmen,n)    ! energy
          call add3    (f1,tmen,tmpr,n)
          call copy    (g1,f1,n)
          call copy    (h1,f1,n)      ! E + p
          call col2    (f1,tmrx,n)
          call invcol2 (f1,tmrh,n)    ! u (E + p)
          call col2    (g1,tmry,n)
          call invcol2 (g1,tmrh,n)    ! v (E + p)
          call col2    (h1,tmrz,n)
          call invcol2 (h1,tmrh,n)    ! w (E + p)
      endif
      call full2face(f1f,f1)
      call full2face(g1f,g1)
      call full2face(h1f,h1)
      call full2face(cuf,u1)
      return
      end
c-----------------------------------------------------------------------
c-----| Service routines 3, printing |-----
c-----------------------------------------------------------------------
      subroutine outfld3d_deal_e1(x,txt10)  ! print out 3d de-al. array 
      include 'SIZE'
      real x(lxd,lyd,lzd,1)
      character*10 txt10
      integer      iz 

      call outfldzd_e1(x,txt10,1) ! de-aliased print
      if(lzd/2 .gt. 1) call outfldzd_e1(x,txt10,lzd/2)
      return
      end
c-----------------------------------------------------------------------
      subroutine outfld3d_deal(x,txt10)  ! print out 3d de-al. array 
      include 'SIZE'
      real x(lxd,lyd,lzd,lelt)
      character*10 txt10
      integer      iz 

      call outfldzd(x,txt10,1) ! de-aliased print
      if(lzd/2 .gt. 1) call outfldzd(x,txt10,lzd/2)
      return
      end
c-----------------------------------------------------------------------
      subroutine outfldzd_e1(x,txt10,iz)  ! print out one elem 2d array on a level
      include 'SIZE'
      real x(lxd,lyd,lzd,1)
      character*10 txt10
      character*18 s(20,20)
      integer      ie, iz 
      call blank(s,18*20*20)
      if(iz.eq.1) write(6,*) 'Dealiased array'
      write(6,108) txt10,iz,lzd
  108   FORMAT(  /,5X,'    ^          ',/,
     $             5X,'  Y |          ',A10,/,
     $             5X,'    |          ',/, ! 'elem.=',I2,'/',I2,/,
     $             5X,'    +---->     ','#z   =',I2,'/',I2,/,
     $             5X,'      X        ')
      if(lxd .gt. 6) then
          write(6,*) 'poly order too high for printing out'
          return
      endif
      do ie=1,1 !nelv
          if(ie.eq.1) then
              istart  = 1
              jstart  = 1
              istride = 1
          elseif(ie.eq.2) then
              istart  = 1
              jstart  = 1 + lxd
              istride = 1
          elseif(ie.eq.3) then
              istart  = 1 + lxd
              jstart  = 1
              istride = 1
          elseif(ie.eq.4) then
              istart  = 1 + lxd
              jstart  = 1 + lxd
              istride = 1
          endif
          i=istart
          do iy=nyd,1,-1
              j = jstart
              do ix=1,nxd
                  write(s(i,j),16) x(ix,iy,iz,ie)
                  j = j+1
              enddo
              i = i+istride
           enddo
   16      format(f14.10)
      enddo 
      do i=1,6
          write(6,17) (s(i,l),l=1,2*lxd)
      enddo
   17 format(12a18)

      write(6,*)

      return
      end
c-----------------------------------------------------------------------
      subroutine outfldzd(x,txt10,iz)  ! print out 2d array on a level
      include 'SIZE'
      real x(lxd,lyd,lzd,lelt)
      character*10 txt10
      character*18 s(20,20)
      integer      ie, iz 
      call blank(s,18*20*20)
      if(iz.eq.1) write(6,*) 'Dealiased array'
      write(6,107) txt10,nelv,nelv,iz,lzd
  107   FORMAT(  /,5X,'    ^          ',/,
     $             5X,'  Y |          ',A10,/,
     $             5X,'    |          ','elem.=',I2,'/',I2,/,
     $             5X,'    +---->     ','#z   =',I2,'/',I2,/,
     $             5X,'      X        ')
      if(nelv .gt. 6) then
          write(6,*) 'too many elements for printing out'
          return
      endif
      if(lxd .gt. 6) then
          write(6,*) 'poly order too high for printing out'
          return
      endif
      do ie=1,nelv
          if(ie.eq.1) then
              istart  = 1
              jstart  = 1
              istride = 1
          elseif(ie.eq.2) then
              istart  = 1
              jstart  = 1 + lxd
              istride = 1
          elseif(ie.eq.3) then
              istart  = 1 + lxd
              jstart  = 1
              istride = 1
          elseif(ie.eq.4) then
              istart  = 1 + lxd
              jstart  = 1 + lxd
              istride = 1
          endif
          i=istart
          do iy=nyd,1,-1
              j = jstart
              do ix=1,nxd
                  write(s(i,j),16) x(ix,iy,iz,ie)
                  j = j+1
              enddo
              i = i+istride
           enddo
   16      format(f14.10)
      enddo 
      do i=1,6
          write(6,17) (s(i,l),l=1,2*lxd)
      enddo
   17 format(12a18)

      write(6,*)

      return
      end
c-----------------------------------------------------------------------
      subroutine outface3d_all(x,txt10)  ! print out 3d array level by levl
      include 'SIZE'
      real x(1), volx(lx1,ly1,lz1,lelt)
      character*10 txt10
      integer      iz 

      call face2full(volx,x)

      call outfld3d_all(volx,txt10,iz)

      return
      end
c-----------------------------------------------------------------------
      subroutine outfacez1(x,txt10)  ! print out 3d array level by levl
      include 'SIZE'
      real x(1), volx(lx1,ly1,lz1,lelt)
      character*10 txt10

      call face2full(volx,x)

      call outfldz1(volx,txt10,1)

      return
      end
c-----------------------------------------------------------------------
      subroutine outface3d(x,txt10)  ! print out 3d array level by levl
      include 'SIZE'
      real x(1), volx(lx1,ly1,lz1,lelt)
      character*10 txt10

      call face2full(volx,x)

      call outfld3d(volx,txt10)

      return
      end
c-----------------------------------------------------------------------
      subroutine outfld3d_all(x,txt10)  ! print out 3d array all level
      include 'SIZE'
      real x(lx1,ly1,lz1,lelt)
      character*10 txt10
      integer      iz 

      do iz=1,lz1
          call outfldz1(x,txt10,iz)
      enddo 

      return
      end
c-----------------------------------------------------------------------
      subroutine outfld3d(x,txt10)  ! print out 3d array 2 slices
      include 'SIZE'
      real x(lx1,ly1,lz1,lelt)
      character*10 txt10
      integer      iz 

c     do iz=1,lz1
c         call outfldz1(x,txt10,iz)
c     enddo 

      call outfldz1(x,txt10,1)
      if(lz1/2 .gt. 1) call outfldz1(x,txt10,lz1/2)
      return
      end
c-----------------------------------------------------------------------
      subroutine outfldz1_bg(x,txt10,iz)  
      include 'SIZE'
      real x(lx1,ly1,lz1,lelt)
      character*10 txt10
      character*18 s(20,20)
      integer      ie, iz 
      call blank(s,18*20*20)
      write(6,106) txt10,nelv,nelv,iz,lz1
  106   FORMAT(  /,5X,'    ^          ',/,
     $             5X,'  Y |          ',A10,/,
     $             5X,'    |          ','elem.=',I2,'/',I2,/,
     $             5X,'    +---->     ','#z   =',I2,'/',I2,/,
     $             5X,'      X        ')
      if(nelv .gt. 9) then
          write(6,*) 'too many elements for printing out'
          return
      endif
      if(lx1 .gt. 6) then
          write(6,*) 'poly order too high for printing out'
          return
      endif
      do ie=1,nelv
          if(ie.eq.1) then
              istart  = 1
              jstart  = 1
              istride = 1
          elseif(ie.eq.2) then
              istart  = 1
              jstart  = 1 + lx1
              istride = 1
          elseif(ie.eq.3) then
              istart  = 1
              jstart  = 1 + 2*lx1
              istride = 1
          elseif(ie.eq.4) then
              istart  = 1 + lx1
              jstart  = 1
              istride = 1
          elseif(ie.eq.5) then
              istart  = 1 + lx1
              jstart  = 1 + lx1
              istride = 1
          elseif(ie.eq.6) then
              istart  = 1 + lx1
              jstart  = 1 + 2*lx1
              istride = 1
          elseif(ie.eq.7) then
              istart  = 1 + 2*lx1
              jstart  = 1
              istride = 1
          elseif(ie.eq.8) then
              istart  = 1 + 2*lx1
              jstart  = 1 + lx1
              istride = 1
          elseif(ie.eq.9) then
              istart  = 1 + 2*lx1
              jstart  = 1 + 2*lx1
              istride = 1
          endif
          i=istart
          do iy=ny1,1,-1
              j = jstart
              do ix=1,nx1
                  write(s(i,j),6) x(ix,iy,iz,ie)
                  j = j+1
              enddo
              i = i+istride
           enddo
    6      format(f14.10)
      enddo 
      do i=1,3*lx1 !6
          write(6,7) (s(i,l),l=1,3*lx1)
      enddo
    7 format(12a18)

      write(6,*)

      return
      end
c-----------------------------------------------------------------------
      subroutine outfldz1_e1(x,txt10)  ! print out 2d array on a level
      include 'SIZE'
      real x(lx1,ly1,lz1,lelt)
      character*10 txt10
      character*18 s(20,20)
      integer      ie, iz 
      call blank(s,18*20*20)
      iz  = 1 
      write(6,106) txt10,nelv,nelv,iz,lz1
  106   FORMAT(  /,5X,'    ^          ',/,
     $             5X,'  Y |          ',A10,/,
     $             5X,'    |          ','elem.=',I2,'/',I2,/,
     $             5X,'    +---->     ','#z   =',I2,'/',I2,/,
     $             5X,'      X        ')
      if(nelv .gt. 6) then
          write(6,*) 'too many elements for printing out'
          return
      endif
      if(lx1 .gt. 6) then
          write(6,*) 'poly order too high for printing out'
          return
      endif
      ie=1
      istart  = 1
      jstart  = 1
      istride = 1

      i=istart
      do iy=ny1,1,-1
         j = jstart
         do ix=1,nx1
            write(s(i,j),6) x(ix,iy,iz,ie)
            j = j+1
         enddo
         i = i+istride
      enddo
    6 format(f14.10)

c     do i=1,6
      do i=1,lx1
          write(6,7) (s(i,l),l=1,ly1)
      enddo
    7 format(12a18)

      write(6,*)

      return
      end
c-----------------------------------------------------------------------
      subroutine outfldz1(x,txt10,iz)  ! print out 2d array on a level
      include 'SIZE'
      real x(lx1,ly1,lz1,lelt)
      character*10 txt10
      character*18 s(20,20)
      integer      ie, iz 
      call blank(s,18*20*20)
      write(6,106) txt10,nelv,nelv,iz,lz1
  106   FORMAT(  /,5X,'    ^          ',/,
     $             5X,'  Y |          ',A10,/,
     $             5X,'    |          ','elem.=',I2,'/',I2,/,
     $             5X,'    +---->     ','#z   =',I2,'/',I2,/,
     $             5X,'      X        ')
      if(nelv .gt. 6) then
          write(6,*) 'too many elements for printing out'
          return
      endif
      if(lx1 .gt. 6) then
          write(6,*) 'poly order too high for printing out'
          return
      endif
      do ie=1,nelv
          if(ie.eq.1) then
              istart  = 1 + lx1
              jstart  = 1
              istride = 1
          elseif(ie.eq.2) then
              istart  = 1 + lx1
              jstart  = 1 + lx1
              istride = 1
          elseif(ie.eq.3) then
              istart  = 1
              jstart  = 1
              istride = 1
          elseif(ie.eq.4) then
              istart  = 1
              jstart  = 1 + lx1
              istride = 1
          endif
          i=istart
          do iy=ny1,1,-1
              j = jstart
              do ix=1,nx1
                  write(s(i,j),6) x(ix,iy,iz,ie)
                  j = j+1
              enddo
              i = i+istride
           enddo
    6      format(f14.10)
      enddo 
      if(nx1.lt.4) then
        do i=1,6
            write(6,7) (s(i,l),l=1,2*lx1)
        enddo
      else
        do i=1,12
            write(6,7) (s(i,l),l=1,2*lx1)
        enddo
      endif
    7 format(12a18)

      write(6,*)

      return
      end
c-----------------------------------------------------------------------
      subroutine printmax_ctr(rhc, rvxc, rvyc, rvzc, Enc, k)  ! 2+3
      include 'SIZE'
      include 'TOTAL'

      n=nelt 
      rhmax = glmax( rhc,n)
      umax  = glmax(rvxc,n)
      vmax  = glmax(rvyc,n)
      wmax  = glmax(rvzc,n)
      enmax = glmax( Enc,n)
      if(nid .eq. 0) write(6,*) 'Center, Step = ', k
      if(nid .eq. 0) write(6,*) 'rh  max = ', rhmax
      if(nid .eq. 0) write(6,*) 'rvx max = ', umax
      if(nid .eq. 0) write(6,*) 'rvy max = ', vmax
      if(nid .eq. 0) write(6,*) 'rvz max = ', wmax
      if(nid .eq. 0) write(6,*) ' En max = ', enmax
      rhmax = glmin( rhc,n)
      umax  = glmin(rvxc,n)
      vmax  = glmin(rvyc,n)
      wmax  = glmin(rvzc,n)
      enmax = glmin( Enc,n)
      if(nid .eq. 0) write(6,*) 'rh  min = ', rhmax
      if(nid .eq. 0) write(6,*) 'rvx min = ', umax
      if(nid .eq. 0) write(6,*) 'rvy min = ', vmax
      if(nid .eq. 0) write(6,*) 'rvz min = ', wmax
      if(nid .eq. 0) write(6,*) ' En min = ', enmax

      return
      end
c-----------------------------------------------------------------------
      subroutine printmax_fac1(r1)  ! 2+3
      include 'SIZE'
      include 'TOTAL'
      real rmax 

      n=nx1*nz1*2*ndim*nelt 
      rmax = glmax(r1 ,n)
      if(nid .eq. 0) write(6,*) 'max = ', rmax
      rmax = glmin(r1 ,n)
      if(nid .eq. 0) write(6,*) 'min = ', rmax

      return
      end
c-----------------------------------------------------------------------
      subroutine printmax_fac(rhf, rvxf, rvyf, rvzf, Enf, k)  ! 2+3
      include 'SIZE'
      include 'TOTAL'

      n=nx1*nz1*2*ndim*nelt 
      rhmax = glmax(rhf ,n)
      umax  = glmax(rvxf,n)
      vmax  = glmax(rvyf,n)
      wmax  = glmax(rvzf,n)
      enmax = glmax( Enf,n)
      if(nid .eq. 0) write(6,*) 'Surface, Step = ', k
      if(nid .eq. 0) write(6,*) 'rh  max = ', rhmax
      if(nid .eq. 0) write(6,*) 'rvx max = ', umax
      if(nid .eq. 0) write(6,*) 'rvy max = ', vmax
      if(nid .eq. 0) write(6,*) 'rvz max = ', wmax
      if(nid .eq. 0) write(6,*) ' En max = ', enmax
      rhmax = glmin(rhf ,n)
      umax  = glmin(rvxf,n)
      vmax  = glmin(rvyf,n)
      wmax  = glmin(rvzf,n)
      enmax = glmin( Enf,n)
      if(nid .eq. 0) write(6,*) 'rh  min = ', rhmax
      if(nid .eq. 0) write(6,*) 'rvx min = ', umax
      if(nid .eq. 0) write(6,*) 'rvy min = ', vmax
      if(nid .eq. 0) write(6,*) 'rvz min = ', wmax
      if(nid .eq. 0) write(6,*) ' En min = ', enmax

      return
      end
c-----------------------------------------------------------------------
      subroutine print_terr(terr,kstep,tmmm,dlt,np1)  ! 2+3
      include 'SIZE'
      include 'TOTAL'
      integer  kstep, np1
      real     terr(1), tmax, tmin, tmmm, dlt

      n = nx1*ny1*nz1*nelt
      tmax = glamax(terr,n)
      tmin = glamin(terr,n)
      if(nid .eq. 0) then 
         write(6,*) 'Step, Time, max error, dt, nx1' ! , Step ', kstep 
         write(6,*) kstep, tmmm, tmax, dlt, np1
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine printamax_1(r1)  ! 2+3
      include 'SIZE'
      include 'TOTAL'

      n = nx1*ny1*nz1*nelt
      rhmax = glamax(r1,n)
      if(nid .eq. 0) write(6,*) 'max = ', rhmax
      rhmax = glamin(r1,n)
      if(nid .eq. 0) write(6,*) 'min = ', rhmax

      return
      end
c-----------------------------------------------------------------------
      subroutine printmax_1(r1)  ! 2+3
      include 'SIZE'
      include 'TOTAL'

      n = nx1*ny1*nz1*nelt
      rhmax = glmax(r1,n)
      if(nid .eq. 0) write(6,*) 'max = ', rhmax
      rhmax = glmin(r1,n)
      if(nid .eq. 0) write(6,*) 'min = ', rhmax

      return
      end
c-----------------------------------------------------------------------
      subroutine printmax(rh, rvx, rvy, rvz, En, k)  ! 2+3
      include 'SIZE'
      include 'TOTAL'

      n = nx1*ny1*nz1*nelt
      rhmax = glmax(rh ,n)
      umax  = glmax(rvx,n)
      vmax  = glmax(rvy,n)
      wmax  = glmax(rvz,n)
      enmax = glmax( En,n)
      if(nid .eq. 0) write(6,*) ' Step = ', k
      if(nid .eq. 0) write(6,*) 'rh  max = ', rhmax
      if(nid .eq. 0) write(6,*) 'rvx max = ', umax
      if(nid .eq. 0) write(6,*) 'rvy max = ', vmax
      if(nid .eq. 0) write(6,*) 'rvz max = ', wmax
      if(nid .eq. 0) write(6,*) ' En max = ', enmax
      rhmax = glmin(rh ,n)
      umax  = glmin(rvx,n)
      vmax  = glmin(rvy,n)
      wmax  = glmin(rvz,n)
      enmax = glmin( En,n)
      if(nid .eq. 0) write(6,*) 'rh  min = ', rhmax
      if(nid .eq. 0) write(6,*) 'rvx min = ', umax
      if(nid .eq. 0) write(6,*) 'rvy min = ', vmax
      if(nid .eq. 0) write(6,*) 'rvz min = ', wmax
      if(nid .eq. 0) write(6,*) ' En min = ', enmax

      return
      end
c-----------------------------------------------------------------------
c----- Stress 
c-----------------------------------------------------------------------
      subroutine comp_sij_e(sij,nij,u,v,w,ur,us,ut,vr,vs,vt,wr,ws,wt,e)
c                                       du_i       du_j
c     Compute the stress tensor S_ij := ----   +   ----
c                                       dx_j       dx_i
c
      include 'SIZE'
      include 'TOTAL'
c
      integer e
c
      real sij(lx1*ly1*lz1,nij)
      real u  (lx1*ly1*lz1,1)
      real v  (lx1*ly1*lz1,1)
      real w  (lx1*ly1*lz1,1)
      real ur (1) , us (1) , ut (1)
     $   , vr (1) , vs (1) , vt (1)
     $   , wr (1) , ws (1) , wt (1)

      real j ! Inverse Jacobian

      n    = nx1-1      ! Polynomial degree
      nxyz = nx1*ny1*nz1

      if (if3d) then     ! 3D CASE
        call local_grad3(ur,us,ut,u,N,e,dxm1,dxtm1)
        call local_grad3(vr,vs,vt,v,N,e,dxm1,dxtm1)
        call local_grad3(wr,ws,wt,w,N,e,dxm1,dxtm1)

        do i=1,nxyz

         j = jacmi(i,e)

         sij(i,1) = j*  ! du/dx + du/dx
     $   2*(ur(i)*rxm1(i,1,1,e)+us(i)*sxm1(i,1,1,e)+ut(i)*txm1(i,1,1,e))

         sij(i,2) = j*  ! dv/dy + dv/dy
     $   2*(vr(i)*rym1(i,1,1,e)+vs(i)*sym1(i,1,1,e)+vt(i)*tym1(i,1,1,e))

         sij(i,3) = j*  ! dw/dz + dw/dz
     $   2*(wr(i)*rzm1(i,1,1,e)+ws(i)*szm1(i,1,1,e)+wt(i)*tzm1(i,1,1,e))

         sij(i,4) = j*  ! du/dy + dv/dx
     $   (ur(i)*rym1(i,1,1,e)+us(i)*sym1(i,1,1,e)+ut(i)*tym1(i,1,1,e) +
     $    vr(i)*rxm1(i,1,1,e)+vs(i)*sxm1(i,1,1,e)+vt(i)*txm1(i,1,1,e) )

         sij(i,5) = j*  ! dv/dz + dw/dy
     $   (wr(i)*rym1(i,1,1,e)+ws(i)*sym1(i,1,1,e)+wt(i)*tym1(i,1,1,e) +
     $    vr(i)*rzm1(i,1,1,e)+vs(i)*szm1(i,1,1,e)+vt(i)*tzm1(i,1,1,e) )

         sij(i,6) = j*  ! du/dz + dw/dx
     $   (ur(i)*rzm1(i,1,1,e)+us(i)*szm1(i,1,1,e)+ut(i)*tzm1(i,1,1,e) +
     $    wr(i)*rxm1(i,1,1,e)+ws(i)*sxm1(i,1,1,e)+wt(i)*txm1(i,1,1,e) )

        enddo

      elseif (ifaxis) then  ! AXISYMMETRIC CASE  

c
c        Notation:                       ( 2  x  Acheson, p. 353)
c                     Cylindrical
c            Nek5k    Coordinates
c
c              x          z
c              y          r
c              z          theta
c

            call setaxdy ( ifrzer(e) )  ! change dytm1 if on-axis
            call local_grad2(ur,us,u,N,1,dxm1,dytm1)
            call local_grad2(vr,vs,v,N,1,dxm1,dytm1)
            call local_grad2(wr,ws,w,N,1,dxm1,dytm1)

            do i=1,nxyz
               j = jacmi(i,e)
               r = ym1(i,1,1,e)                              ! Cyl. Coord:

               sij(i,1) = j*  ! du/dx + du/dx              ! e_zz
     $           2*(ur(i)*rxm1(i,1,1,e)+us(i)*sxm1(i,1,1,e))

               sij(i,2) = j*  ! dv/dy + dv/dy              ! e_rr
     $           2*(vr(i)*rym1(i,1,1,e)+vs(i)*sym1(i,1,1,e))

               if (r.gt.0) then                              ! e_@@
                  sij(i,3) = v(i,e)/r  ! v / r 
               else
                  sij(i,3) = j*  ! L'Hopital's rule: e_@@ = dv/dr
     $            2*(vr(i)*rym1(i,1,1,e)+vs(i)*sym1(i,1,1,e))
               endif

               sij(i,4) = j*  ! du/dy + dv/dx             ! e_zr
     $            ( ur(i)*rym1(i,1,1,e)+us(i)*sym1(i,1,1,e) +
     $              vr(i)*rxm1(i,1,1,e)+vs(i)*sxm1(i,1,1,e) )

               if (yyyr.gt.0) then                             ! e_r@
                  sij(i,5) = j*  ! dw/dy 
     $              ( wr(i)*rym1(i,1,1,e)+ws(i)*sym1(i,1,1,e) )
     $              - w(i,e) / r
               else
                  sij(i,5) = 0
               endif

               sij(i,6) = j*  ! dw/dx                     ! e_@z
     $            ( wr(i)*rxm1(i,1,1,e)+ws(i)*sxm1(i,1,1,e) )

            enddo

      else              ! 2D CASE

            call local_grad2(ur,us,u,N,1,dxm1,dxtm1)
            call local_grad2(vr,vs,v,N,1,dxm1,dxtm1)

            do i=1,nxyz
               j = jacmi(i,e)

               sij(i,1) = j*  ! du/dx + du/dx
     $           2*(ur(i)*rxm1(i,1,1,e)+us(i)*sxm1(i,1,1,e))

               sij(i,2) = j*  ! dv/dy + dv/dy
     $           2*(vr(i)*rym1(i,1,1,e)+vs(i)*sym1(i,1,1,e))

               sij(i,3) = j*  ! du/dy + dv/dx
     $           (ur(i)*rym1(i,1,1,e)+us(i)*sym1(i,1,1,e) +
     $            vr(i)*rxm1(i,1,1,e)+vs(i)*sxm1(i,1,1,e) )

            enddo
      endif
      return
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c-----| Limiters |-----
c------ Thu Oct  8 10:13:52 CDT 2015
c------ Code in Zhang Shu limiter both for 2D and 3D 
c-----------------------------------------------------------------------
      subroutine slopelim(rh,rvx,rvy,rvz,En)
      include 'SIZE'
      include 'TOTAL'
      real    rh(1), rvx(1), rvy(1), rvz(1), en(1)
      logical ifsod

      ifsod=.true.  ! unless for sod problem 
      ifsod=.false.  ! unless for sod problem 

      if(ifsod) then
          call slopelim1d(rh,rvx,rvy,rvz,En)
      else if(if3d) then
          write(6,*) '3D limiters...'
c         call slopelim3d_zs(rh,rvx,rvy,rvz,En) !! zhang shu 
c         call slopelim3d(rh,rvx,rvy,rvz,En) !! Luo, baum
          call exitt
      else ! 2d limiter
c     Wed Aug 12 16:56:13 CDT 2015
c     For forward facing step problem it seems I have to 
c     get 2D limiters to work. Otherwise it just won't compute 
c         write(6,*) '2D limiters...'
c         call slopelim2d(rh,rvx,rvy,rvz,En) !! no confidence
          call slopelim2d_zs(rh,rvx,rvy,rvz,En) !! zhang shu 
      endif

      return
      end
c-----------------------------------------------------------------------
c------ 3D limiter, Zhang, Shu, 2010 
c------ Sun Oct 11 14:20:46 CDT 2015
c-----------------------------------------------------------------------
      subroutine slopelim3d_zs(rha,rvxa,rvya,rvza,enra) ! 2 
      include 'SIZE'
      include 'TOTAL'
      include 'DG'
      include 'DGUSE'
      parameter(le=lx1*ly1*lz1)
      parameter(l3f= lx1*lz1*6*lelt)
      real      rha(le,lelt),  rvxa(le,lelt),  rvya(le,lelt)
     $      ,  rvza(le,lelt),  enra(le,lelt)

      return
      end
c-----------------------------------------------------------------------
c------ 3D limiter, Luo, Baum, Lohner, 2008
c------ Tue Aug 18 21:29:17 CDT 2015
c-----------------------------------------------------------------------
      subroutine slopelim3d(rha,rvxa,rvya,rvza,enra) ! 3
c this should be able to easily transform into 2d case 
      include 'SIZE'
      include 'TOTAL'
      include 'DG'
      include 'DGUSE'
      parameter(le=lx1*ly1*lz1)
      parameter(l3f= lx1*lz1*6*lelt)
      real      rha(1),  rvxa(1),  rvya(1),  rvza(1),  enra(1)
      common /tdlimr/   lrh(lt), lrvx(lt), lrvy(lt), lrvz(lt)
     $                , lenr(lt), lmp(lt)
      real       lrh, lrvx, lrvy, lrvz, lenr, lmp

      common /tdlimi/  ngh(6,lelt)
      integer*8  ngh

      real rhc(lelt), rhuc(lelt), rhvc(lelt), rhwc(lelt)
     $   , enc(lelt), pprc(lelt) ! elem avg 
      real x_c(lelt), y_c(lelt), z_c(lelt)
      real rhcv(le,lelt),rhuv(le,lelt),rhvv(le,lelt)
     $   , rhwv(le,lelt),encv(le,lelt)
      real rhf(lf), rxf(lf), ryf(lf), rzf(lf), enf(lf) 
      real mxrh(lf), mxru(lf), mxrv(lf),mxrw(lf), mxen(lf) 
      real mnrh(lf), mnru(lf), mnrv(lf),mnrw(lf), mnen(lf) 
      real mxc1(lelt),mxc2(lelt),mxc3(lelt),mxc4(lelt),mxc5(lelt)
     $   , mnc1(lelt),mnc2(lelt),mnc3(lelt),mnc4(lelt),mnc5(lelt)
      real  al1(le,lelt),al2(le,lelt),al3(le,lelt)
     $    , al4(le,lelt),al5(le,lelt)
      real  alc1(lelt),alc2(lelt),alc3(lelt)
     $    , alc4(lelt),alc5(lelt)

!if in range, default .true.
      logical ifin1(lelt),ifin2(lelt),ifin3(lelt)
     $       ,ifin4(lelt),ifin5(lelt) 

      nf=nx1*nz1*2*ndim*nelt !total number of points on faces
      n=nx1*ny1*nz1*nelt

c    1. Find element avg. coord., rh, ru, rv, en 
      call elm_avg3(x_c, y_c, z_c, rhc, rhuc, rhvc, rhwc, enc
     $                  , rha,  rvxa,  rvya, rvza, enra)

c    2. Find avg(U)max and min, boundary ? 
      call fill2full_all(rhcv,rhuv,rhvv,rhwv,encv
     $                   ,rhc,rhuc,rhvc,rhwc,enc)
c      Facial central value
      call full2face_all(rhf, rxf, ryf, rzf, enf
     $                 ,rhcv,rhuv,rhvv,rhwv,encv) 

c    3. Find elements that need limiting /Shock detector ? 
      call copy_all(mxrh,mxru,mxrv,mxrw,mxen
     $             , rhf, rxf, ryf, rzf, enf, nf)
      call copy_all(mnrh,mnru,mnrv,mnrw,mnen
     $             , rhf, rxf, ryf, rzf, enf, nf)

      call gs_op_mx( mxrh, mxru, mxrv, mxrw, mxen) ! no boundary ? 
      call gs_op_mn( mnrh, mnru, mnrv, mnrw, mnen) ! no boundary ? 

      call exm_frm_fac( mxc1, mxc2, mxc3, mxc4, mxc5
     $                , mnc1, mnc2, mnc3, mnc4, mnc5
     $                , mxrh, mxru, mxrv, mxrw, mxen
     $                , mnrh, mnru, mnrv, mnrw, mnen ) ! 3 

      call init_ifrng(ifin1, ifin2, ifin3, ifin4, ifin5)
      call set_ifrng (ifin1, ifin2, ifin3, ifin4, ifin5
     $              , mxc1, mxc2, mxc3, mxc4, mxc5
     $              , mnc1, mnc2, mnc3, mnc4, mnc5
     $              , rha,  rvxa,  rvya, rvza, enra)

c    4. Find factor alpha 
c   Get alpha equation in here , only for elements with ifin? false? 
c    Tue Aug 18 22:54:21 CDT 2015
c    alpha = 0 means evaluating with the mean of the element
c    maybe not a good idea for flagging do not limit 
c    alpha = 1 means linearizing 
c    also dose not seem like a good idea  ??? what to use? 
      call init_alc(alc1, alc2, alc3, alc4, alc5)
c     write(6,*) 'Centered alpha values init '
c     call printmax_ctr(alc1, alc2, alc3, alc4, alc5,0)
c     do ie =1, nelt
c         write(6,*) ie, ' = ie', alc1(ie), ' = alc1(ie) ',ifin1(ie)
c         write(6,*) mnc1(ie), 'min center '
c     enddo 
      call eval_alp(al1,al2,al3,al4,al5
     $            , alc1, alc2, alc3, alc4, alc5
     $            ,ifin1, ifin2, ifin3, ifin4, ifin5
     $            , rhc,rhuc,rhvc,rhwc,enc
     $            , rha, rvxa, rvya, rvza, enra  
     $            , mxc1, mxc2, mxc3, mxc4, mxc5
     $            , mnc1, mnc2, mnc3, mnc4, mnc5)
c     write(6,*) 'Centered alpha values '
c     call printmax_ctr(alc1, alc2, alc3, alc4, alc5,0)
c     do ie =1, nelt
c         write(6,*) ie, ' = ie', alc1(ie), ' = alc1(ie) ' 
c     enddo 
c     do i  =1, le
c         write(6,*) i , ' = i ', al1(i,16), ' = al1(i,16) ' 
c     enddo 
c    Seems problematic to have 0 values here 
c    Wed Aug 26 23:28:42 CDT 2015

      call copy_all(lrh, lrvx, lrvy, lrvz, lenr
     $            , rha, rvxa, rvya, rvza, enra,n) 
c     write(6,*) 'Before lim grad  '
c     call printmax_ctr(rhc, rhuc, rhvc, rhwc, enc,0)
c    5. Find gradient, then limit 
      call lim_soln(lrh, lrvx, lrvy, lrvz, lenr
     $            ,ifin1, ifin2, ifin3, ifin4, ifin5
     $            , alc1, alc2, alc3, alc4, alc5
     $            , rhc, rhuc, rhvc, rhwc,  enc
     $            , x_c,  y_c,  z_c
     $            , rha, rvxa, rvya, rvza, enra) 
c     write(6,*) 'After lim grad  '
c     call printmax_ctr(lrh, lrvx, lrvy, lrvz, lenr,0)
c     stop
c      copy back into solution fields
      call copy_all(rha,rvxa,rvya,rvza,enra
     $            , lrh,lrvx,lrvy,lrvz,lenr,n)
c     
      return
      end
c-----------------------------------------------------------------------
      subroutine fill2full_all(fcv1,fcv2,fcv3,fcv4,fcv5
     $                       , fc1,fc2,fc3,fc4,fc5)
      include 'SIZE'
      include 'TOTAL' 
      parameter(le=lx1*ly1*lz1)
      real  fcv1(le,1),fcv2(le,1),fcv3(le,1)
     $    , fcv4(le,1),fcv5(le,1)
      real  fc1(1),fc2(1),fc3(1),fc4(1),fc5(1)
      integer ie, ne
      ne=nx1*ny1*nz1

      do ie=1,nelt ! e is negative ??? 
c         write(6,*) 'ie = ',ie 
          call cfill(fcv1(1,ie),fc1(ie),ne)
          call cfill(fcv2(1,ie),fc2(ie),ne)
          call cfill(fcv3(1,ie),fc3(ie),ne)
          call cfill(fcv4(1,ie),fc4(ie),ne)
          call cfill(fcv5(1,ie),fc5(ie),ne)
      enddo 

      return
      end
c-----------------------------------------------------------------------
      subroutine lim_soln(lrh, lrvx, lrvy, lrvz, lenr
     $                 ,ifin1, ifin2, ifin3, ifin4, ifin5
     $                 , alc1, alc2, alc3, alc4, alc5
     $                 , rhc, rhuc, rhvc, rhwc,  enc
     $                 , x_c,  y_c,  z_c
     $                 , rha, rvxa, rvya, rvza, enra) 
c  Need the gradient of each solution, seems expensive 
c  also if the solution has very huge derivative ( in the case 
c  of discontinuities), will this limiting process take care of 
c  it? 
      include 'SIZE'
      include 'TOTAL'
      include 'DG'
      include 'DGUSE'
      parameter(le=lx1*ly1*lz1)
      logical   ifin1(1), ifin2(1), ifin3(1)
     $        , ifin4(1), ifin5(1) 
      real  rhc(1), rhuc(1), rhvc(1), rhwc(1), enc(1) ! elem avg 
     $    , x_c(1), y_c(1), z_c(1)
      real  rha(le,1),rvxa(le,1),rvya(le,1)
     $    , rvza(le,1),enra(le,1)
      real  alc1(1),alc2(1),alc3(1)
     $    , alc4(1),alc5(1)
      real  lrh (le,1), lrvx(le,1), lrvy(le,1), lrvz(le,1)
     $    , lenr(le,1) 
      real  tmx(le,lelt), tmy(le,lelt), tmz(le,lelt) 

      integer e, i, j, k 

      ne = nx1*ny1*nz1
      nf = nx1*nz1*2*ndim*nelt

      call gradm1(tmx,tmy,tmz,rha)
      do e=1,nelt
          if(.not.ifin1(e)) then ! rho 
              do i = 1, ne
                  lrh(i,e) = rhc(e) 
     $        + alc1(e)*( tmx(i,e)*(xm1(i,1,1,e)-x_c(e))
     $                  + tmy(i,e)*(ym1(i,1,1,e)-y_c(e))
     $                  + tmz(i,e)*(zm1(i,1,1,e)-z_c(e)))
                  if(lrh(i,e).le.1e-8) then
                      write(6,*) 'in limiting routine' ! check this  
                      write(6,*) 'i =',i,'e =',e,'lrh(i,e) =',lrh(i,e)
                  endif
              enddo 
          endif
      enddo  

      call gradm1(tmx,tmy,tmz,rvxa)
      do e=1,nelt
          if(.not.ifin2(e)) then ! rho u 
              do i = 1, ne
                  lrvx(i,e) = rhuc(e) 
     $        + alc2(e)*( tmx(i,e)*(xm1(i,1,1,e)-x_c(e))
     $                  + tmy(i,e)*(ym1(i,1,1,e)-y_c(e))
     $                  + tmz(i,e)*(zm1(i,1,1,e)-z_c(e)))
              enddo 
          endif
      enddo  

      call gradm1(tmx,tmy,tmz,rvya)
      do e=1,nelt
          if(.not.ifin3(e)) then ! rho v 
              do i = 1, ne
                  lrvy(i,e) = rhvc(e) 
     $        + alc3(e)*( tmx(i,e)*(xm1(i,1,1,e)-x_c(e))
     $                  + tmy(i,e)*(ym1(i,1,1,e)-y_c(e))
     $                  + tmz(i,e)*(zm1(i,1,1,e)-z_c(e)))
              enddo 
          endif
      enddo  

      call gradm1(tmx,tmy,tmz,rvza)
      do e=1,nelt
          if(.not.ifin4(e)) then ! rho w 
              do i = 1, ne
                  lrvz(i,e) = rhwc(e) 
     $        + alc4(e)*( tmx(i,e)*(xm1(i,1,1,e)-x_c(e))
     $                  + tmy(i,e)*(ym1(i,1,1,e)-y_c(e))
     $                  + tmz(i,e)*(zm1(i,1,1,e)-z_c(e)))
              enddo 
          endif
      enddo  

      call gradm1(tmx,tmy,tmz,enra)
      do e=1,nelt
          if(.not.ifin5(e)) then ! Ener
              do i = 1, ne
                  lenr(i,e) = enc(e) 
     $        + alc5(e)*( tmx(i,e)*(xm1(i,1,1,e)-x_c(e))
     $                  + tmy(i,e)*(ym1(i,1,1,e)-y_c(e))
     $                  + tmz(i,e)*(zm1(i,1,1,e)-z_c(e)))
              enddo 
          endif
      enddo  

      return
      end
c-----------------------------------------------------------------------
      subroutine eval_alp(al1,al2,al3,al4,al5
     $                 , alc1, alc2, alc3, alc4, alc5
     $                 ,ifin1, ifin2, ifin3, ifin4, ifin5
     $                 ,  rhc,rhuc,rhvc,rhwc,enc
     $                 , rha, rvxa, rvya, rvza, enra  
     $                 , mxc1, mxc2, mxc3, mxc4, mxc5
     $                 , mnc1, mnc2, mnc3, mnc4, mnc5)
      include 'SIZE'
      include 'TOTAL'
      include 'DG'
      include 'DGUSE'
      parameter(le=lx1*ly1*lz1)
      logical   ifin1(1), ifin2(1), ifin3(1)
     $        , ifin4(1), ifin5(1) 
      real  rhc(1), rhuc(1), rhvc(1), rhwc(1), enc(1) ! elem avg 
      real  rha(le,1), rvxa(le,1),rvya(le,1)
     $    , rvza(le,1),enra(le,1)
      real  mxc1(1),mxc2(1),mxc3(1),mxc4(1),mxc5(1)
     $   ,  mnc1(1),mnc2(1),mnc3(1),mnc4(1),mnc5(1)
      real  al1(le,1),al2(le,1),al3(le,1)
     $    , al4(le,1),al5(le,1)
      real  alc1(1),  alc2(1),  alc3(1)
     $    , alc4(1),  alc5(1)

      integer e, i, j, k 

      ne = nx1*ny1*nz1

      do e=1,nelt
          if(.not.ifin1(e)) then ! rho 
              do i = 1, ne
                  if(rha(i,e) - rhc(e).ge. 1.e-8 ) then
                      al1(i,e)=
     $              min(1.,(mxc1(e)-rhc(e))/(rha(i,e)-rhc(e)))
                  else if(rha(i,e) - rhc(e).le. -1.e-8) then
                      al1(i,e)=
     $              min(1.,(mnc1(e)-rhc(e))/(rha(i,e)-rhc(e)))
                  else 
                      al1(i,e)= 1.
                  endif 
c this is the problem , value very close to 0 
c  Thu Aug 27 10:52:05 CDT 2015
c                 if(i.eq.6.and.e.eq.16) then
c                     write(6,*) 'i,',i,'e,',e
c                     write(6,*) 'mxc1(e),',mxc1(e),'mnc1(e),',mnc1(e)
c                     write(6,*) 'rhc (e),',rhc(e),'rha(i,e),',rha(i,e)
c                     write(6,*) 'al1(i,e),',
c    $              (mxc1(e)-rhc(e))/(rha(i,e)-rhc(e))
c                 endif

              enddo 
          endif
          if(.not.ifin2(e)) then ! rvx 
              do i = 1, ne
                  if(rvxa(i,e) .gt. rhuc(e)) then
                      al2(i,e)=
     $              min(1.,(mxc2(e)-rhuc(e))/(rvxa(i,e)-rhuc(e)))
                  else if(rvxa(i,e) .lt. rhuc(e)) then
                      al2(i,e)=
     $              min(1.,(mnc2(e)-rhuc(e))/(rvxa(i,e)-rhuc(e)))
                  else 
                      al2(i,e)= 1.
                  endif 
              enddo 
          endif
          if(.not.ifin3(e)) then ! rvy 
              do i = 1, ne
                  if(rvya(i,e) .gt. rhvc(e)) then
                      al3(i,e)=
     $            min(1.,(mxc3(e)-rhvc(e))/(rvya(i,e)-rhvc(e)))
                  else if(rvya(i,e) .lt. rhvc(e)) then
                      al3(i,e)=
     $            min(1.,(mnc3(e)-rhvc(e))/(rvya(i,e)-rhvc(e)))
                  else 
                      al3(i,e)= 1.
                  endif 
              enddo 
          endif
          if(.not.ifin4(e)) then ! rvz 
              do i = 1, ne
                  if(rvza(i,e) .gt. rhwc(e)) then
                      al4(i,e)=
     $            min(1.,(mxc4(e)-rhwc(e))/(rvza(i,e)-rhwc(e)))
                  else if(rvza(i,e) .lt. rhwc(e)) then
                      al4(i,e)=
     $            min(1.,(mnc4(e)-rhwc(e))/(rvza(i,e)-rhwc(e)))
                  else 
                      al4(i,e)= 1.
                  endif 
              enddo 
          endif
          if(.not.ifin5(e)) then ! enr 
              do i = 1, ne
                  if(enra(i,e) .gt. enc(e)) then
                      al5(i,e)=
     $            min(1.,(mxc5(e)-enc(e))/(enra(i,e)-enc(e)))
                  else if(enra(i,e) .lt. enc(e)) then
                      al5(i,e)=
     $            min(1.,(mnc5(e)-enc(e))/(enra(i,e)-enc(e)))
                  else 
                      al5(i,e)= 1.
                  endif 
              enddo 
          endif
      enddo  ! this is for each point 
c now each element only has one alpha that get used, which will be 
c  the minimum out of all the alpha's within one element 
      do e=1,nelt
          if(.not.ifin1(e)) then 
              alc1(e) = vlmin(al1(1,e),ne)
          endif
          if(.not.ifin2(e)) then 
              alc2(e) = vlmin(al2(1,e),ne)
          endif
          if(.not.ifin3(e)) then 
              alc3(e) = vlmin(al3(1,e),ne)
          endif
          if(.not.ifin4(e)) then 
              alc4(e) = vlmin(al4(1,e),ne)
          endif
          if(.not.ifin5(e)) then 
              alc5(e) = vlmin(al5(1,e),ne)
          endif
      enddo  ! minimum within one element 

      return
      end
c-----------------------------------------------------------------------
      subroutine init_alc(alc1, alc2, alc3, alc4, alc5)
      include 'SIZE'
      include 'TOTAL'
      include 'DG'
      include 'DGUSE'
      parameter(le=lx1*ly1*lz1)
      real    alc1(1), alc2(1), alc3(1), alc4(1), alc5(1) 
      integer e

      do e=1,nelt
          alc1(e) = 1.  
          alc2(e) = 1.  
          alc3(e) = 1.  
          alc4(e) = 1.  
          alc5(e) = 1.
      enddo 

      return
      end
c-----------------------------------------------------------------------
      subroutine init_ifrng(ifin1, ifin2, ifin3, ifin4, ifin5)
      include 'SIZE'
      include 'TOTAL'
      include 'DG'
      include 'DGUSE'
      parameter(le=lx1*ly1*lz1)
      logical   ifin1(1), ifin2(1), ifin3(1)
     $        , ifin4(1), ifin5(1) 
      integer e

      do e=1,nelt
          ifin1(e) = .true.
          ifin2(e) = .true.
          ifin3(e) = .true.
          ifin4(e) = .true.
          ifin5(e) = .true.
      enddo 

      return
      end
c-----------------------------------------------------------------------
      subroutine set_ifrng(ifin1, ifin2, ifin3, ifin4, ifin5
     $                    , mxrh, mxru, mxrv, mxrw, mxen
     $                    , mnrh, mnru, mnrv, mnrw, mnen
     $                    , rha,  rvxa,  rvya, rvza, enra) ! 3 
      include 'SIZE'
      include 'TOTAL'
      include 'DG'
      include 'DGUSE'
      parameter(le=lx1*ly1*lz1)
      logical   ifin1(1), ifin2(1), ifin3(1)
     $        , ifin4(1), ifin5(1) 
      real      rha(le,1), rvxa(le,1), rvya(le,1), rvza(le,1)
     $        , enra(le,1)
      real      mxrh(lelt), mxru(lelt), mxrv(lelt), mxrw(lelt)
     $        , mxen(lelt), mnrh(lelt), mnru(lelt), mnrv(lelt)
     $        , mnrw(lelt), mnen(lelt) 

      integer e, i, j, k 

      ne = nx1*ny1*nz1
      nxz    = nx1*nz1
      nfaces = 2*ndim
      nver = 2**(ndim)
      nf = nx1*nz1*2*ndim*nelt

      do e=1,nelt
      do i = 1, ne
          if( rha(i,e) .lt. mnrh(e) .or.  rha(i,e) .gt. mxrh(e)) then
              ifin1(e) = .false. 
          endif
          if(rvxa(i,e) .lt. mnru(e) .or. rvxa(i,e) .gt. mxru(e)) then
              ifin2(e) = .false. 
          endif
          if(rvya(i,e) .lt. mnrv(e) .or. rvya(i,e) .gt. mxrv(e)) then
              ifin3(e) = .false. 
          endif
          if(rvza(i,e) .lt. mnrw(e) .or. rvza(i,e) .gt. mxrw(e)) then
              ifin4(e) = .false. 
          endif
          if(enra(i,e) .lt. mnen(e) .or. enra(i,e) .gt. mxen(e)) then
              ifin5(e) = .false. 
          endif
      enddo 
      enddo 

      return
      end
c-----------------------------------------------------------------------
      subroutine exm_frm_fac( max1, max2, max3, max4, max5
     $                      , min1, min2, min3, min4, min5
     $                      , mx1, mx2, mx3, mx4, mx5
     $                      , mn1, mn2, mn3, mn4, mn5) ! 3 
      include 'SIZE'
      include 'TOTAL'
      include 'DG'
      include 'DGUSE'
      parameter(le=lx1*ly1*lz1)
      real      mx1(lf), mx2(lf), mx3(lf), mx4(lf), mx5(lf)
     $        , mn1(lf), mn2(lf), mn3(lf), mn4(lf), mn5(lf)
      real      max1(lelt), max2(lelt), max3(lelt), max4(lelt)
     $        , max5(lelt), min1(lelt), min2(lelt), min3(lelt)
     $        , min4(lelt), min5(lelt)

      integer e, i, j, k 
      nxz    = nx1*nz1
      nfaces = 2*ndim
      nf = nx1*nz1*2*ndim*nelt
      nfe= nx1*nz1*2*ndim

      do e=1,nelt
          max1(e)=max(mx1(1+      nfe*(e-1)),mx1(1+  nxz+nfe*(e-1))
     $             ,  mx1(1+2*nxz+nfe*(e-1)),mx1(1+3*nxz+nfe*(e-1))
     $             ,  mx1(1+4*nxz+nfe*(e-1)),mx1(1+5*nxz+nfe*(e-1)))
          min1(e)=min(mn1(1+      nfe*(e-1)),mn1(1+  nxz+nfe*(e-1))
     $             ,  mn1(1+2*nxz+nfe*(e-1)),mn1(1+3*nxz+nfe*(e-1))
     $             ,  mn1(1+4*nxz+nfe*(e-1)),mn1(1+5*nxz+nfe*(e-1)))
          max2(e)=max(mx2(1+      nfe*(e-1)),mx2(1+  nxz+nfe*(e-1))
     $             ,  mx2(1+2*nxz+nfe*(e-1)),mx2(1+3*nxz+nfe*(e-1))
     $             ,  mx2(1+4*nxz+nfe*(e-1)),mx2(1+5*nxz+nfe*(e-1)))
          min2(e)=min(mn2(1+      nfe*(e-1)),mn2(1+  nxz+nfe*(e-1))
     $             ,  mn2(1+2*nxz+nfe*(e-1)),mn2(1+3*nxz+nfe*(e-1))
     $             ,  mn2(1+4*nxz+nfe*(e-1)),mn2(1+5*nxz+nfe*(e-1)))
          max3(e)=max(mx3(1+      nfe*(e-1)),mx3(1+  nxz+nfe*(e-1))
     $             ,  mx3(1+2*nxz+nfe*(e-1)),mx3(1+3*nxz+nfe*(e-1))
     $             ,  mx3(1+4*nxz+nfe*(e-1)),mx3(1+5*nxz+nfe*(e-1)))
          min3(e)=min(mn3(1+      nfe*(e-1)),mn3(1+  nxz+nfe*(e-1))
     $             ,  mn3(1+2*nxz+nfe*(e-1)),mn3(1+3*nxz+nfe*(e-1))
     $             ,  mn3(1+4*nxz+nfe*(e-1)),mn3(1+5*nxz+nfe*(e-1)))
          max4(e)=max(mx4(1+      nfe*(e-1)),mx4(1+  nxz+nfe*(e-1))
     $             ,  mx4(1+2*nxz+nfe*(e-1)),mx4(1+3*nxz+nfe*(e-1))
     $             ,  mx4(1+4*nxz+nfe*(e-1)),mx4(1+5*nxz+nfe*(e-1)))
          min4(e)=min(mn4(1+      nfe*(e-1)),mn4(1+  nxz+nfe*(e-1))
     $             ,  mn4(1+2*nxz+nfe*(e-1)),mn4(1+3*nxz+nfe*(e-1))
     $             ,  mn4(1+4*nxz+nfe*(e-1)),mn4(1+5*nxz+nfe*(e-1)))
          max5(e)=max(mx5(1+      nfe*(e-1)),mx5(1+  nxz+nfe*(e-1))
     $             ,  mx5(1+2*nxz+nfe*(e-1)),mx5(1+3*nxz+nfe*(e-1))
     $             ,  mx5(1+4*nxz+nfe*(e-1)),mx5(1+5*nxz+nfe*(e-1)))
          min5(e)=min(mn5(1+      nfe*(e-1)),mn5(1+  nxz+nfe*(e-1))
     $             ,  mn5(1+2*nxz+nfe*(e-1)),mn5(1+3*nxz+nfe*(e-1))
     $             ,  mn5(1+4*nxz+nfe*(e-1)),mn5(1+5*nxz+nfe*(e-1)))
      enddo 

      return
      end
c-----------------------------------------------------------------------
c----- Service routines 4, finding max min 
c-----------------------------------------------------------------------
      real function min_6( m1, m2, m3, m4, m5, m6)  ! does not work 
      real m1, m2, m3, m4, m5, m6

      min_6 = min(m1,m2,m3,m4,m5,m6) ! could it be ? 

      return
      end
c-----------------------------------------------------------------------
      real function max_6( m1, m2, m3, m4, m5, m6)  ! does not work 
      real m1, m2, m3, m4, m5, m6

      max_6 = max(m1,m2,m3,m4,m5,m6) ! could it be ? 

      return
      end
c-----------------------------------------------------------------------
      subroutine gs_op_mx(mx1,mx2,mx3,mx4,mx5) ! 3
c     on boundaries the values shouldn't change 
      include 'SIZE'
      include 'TOTAL'
      include 'DG'
      real mx1(1),  mx2(1),  mx3(1),  mx4(1),  mx5(1)

      call gs_op (dg_hndl,mx1,1,4,0)  ! 4 ==> max
      call gs_op (dg_hndl,mx2,1,4,0)  ! 4 ==> max
      call gs_op (dg_hndl,mx3,1,4,0)  ! 4 ==> max
      call gs_op (dg_hndl,mx4,1,4,0)  ! 4 ==> max
      call gs_op (dg_hndl,mx5,1,4,0)  ! 4 ==> max

      return
      end
c-----------------------------------------------------------------------
      subroutine gs_op_mn(mn1,mn2,mn3,mn4,mn5) ! 3
c     on boundaries the values shouldn't change 
      include 'SIZE'
      include 'TOTAL'
      include 'DG'
      real      mn1(1),  mn2(1),  mn3(1),  mn4(1),  mn5(1)

      call gs_op (dg_hndl,mn1,1,3,0)  ! 3 ==> min
      call gs_op (dg_hndl,mn2,1,3,0)  ! 3 ==> min
      call gs_op (dg_hndl,mn3,1,3,0)  ! 3 ==> min
      call gs_op (dg_hndl,mn4,1,3,0)  ! 3 ==> min
      call gs_op (dg_hndl,mn5,1,3,0)  ! 3 ==> min

      return
      end
c-----------------------------------------------------------------------
      subroutine elm_avg3(x_c, y_c, z_c, rhc, rhuc, rhvc, rhwc, enc
     $                  , rh,  rvx,  rvy, rvz, enr) ! no rvz 
      include 'SIZE'
      include 'TOTAL'
c
c   Get element averaged values
c     Input:  rh(le,nelt), rvx(le,nelt), rvy(le,nelt), enr(le,nelt)
c     Output: x_c(nelt),y_c(nelt), z_c(lelt) 
c            ,rhc(nelt), rhvc(nelt), rhvc(nelt), enc(nelt)
c     
      parameter(le = lx1*ly1*lz1)
      parameter(lf = lx1*lz1*2*ldim*lelt)
      real      x_c(1), y_c(1), z_c(1), rhc(1), rhuc(1), rhvc(1)
     $        , rhwc(1), enc(1)
     $        , rh(le,1), rvx(le,1), rvy(le,1), rvz(le,1), enr(le,1)
      real      ones(le), tmp(le) 
      real      sum1,sum2,sum3,sum4,sum5,sum6,sum7,sum8,sum9

      nxz    = nx1*nz1
      ne     = nx1*ny1*nz1
      nfaces = 2*ndim 
      nf = nx1*nz1*2*ndim*nelt

      do ie = 1, nelt
c - - 
          call cfill(ones,1.,ne)
          call col2(ones,bm1(1,1,1,ie),ne)
          sum1 = vlsum(ones,ne) ! volume for elem #ie
c - - 
          call col3   (tmp,bm1(1,1,1,ie),xm1(1,1,1,ie),ne)
          sum2 = vlsum(tmp ,ne) ! vol. sum for x
c - - 
          call col3   (tmp,bm1(1,1,1,ie),ym1(1,1,1,ie),ne)
          sum3 = vlsum(tmp ,ne) ! vol. sum for y
c - - 
          call col3   (tmp,bm1(1,1,1,ie),zm1(1,1,1,ie),ne)
          sum9 = vlsum(tmp ,ne) ! vol. sum for z
c - - 
          call col3   (tmp,bm1(1,1,1,ie),rh (1,ie),ne)
          sum4 = vlsum(tmp ,ne) ! vol. sum for rh
c - - 
          call col3   (tmp,bm1(1,1,1,ie),rvx(1,ie),ne)
          sum5 = vlsum(tmp ,ne) ! vol. sum for rvx
c - - 
          call col3   (tmp,bm1(1,1,1,ie),rvy(1,ie),ne)
          sum6 = vlsum(tmp ,ne) ! vol. sum for rvy
c - - 
          call col3   (tmp,bm1(1,1,1,ie),rvz(1,ie),ne)
          sum7 = vlsum(tmp ,ne) ! vol. sum for rvz
c - - 
          call col3   (tmp,bm1(1,1,1,ie),enr(1,ie),ne)
          sum8 = vlsum(tmp ,ne) ! vol. sum for rvy
c 
          x_c (ie) = sum2/sum1
          y_c (ie) = sum3/sum1
          z_c (ie) = sum9/sum1
          rhc (ie) = sum4/sum1
          rhuc(ie) = sum5/sum1
          rhvc(ie) = sum6/sum1
          rhwc(ie) = sum7/sum1
          enc (ie) = sum8/sum1

      enddo 

      return
      end
c-----------------------------------------------------------------------
c------ 2D limiter
c-----------------------------------------------------------------------
c------ Zhang, Shu 2009 
c------ Sun Oct 11 14:21:40 CDT 2015
c-----------------------------------------------------------------------
      subroutine slopelim2d_zs(rha,rvxa,rvya,rvza,enra) ! 2 
      include 'SIZE'
      include 'TOTAL'
      include 'DG'
      include 'DGUSE'
      parameter(le=lx1*ly1*lz1)
      parameter(l3f= lx1*lz1*6*lelt)

      real      rha(le,lelt),  rvxa(le,lelt),  rvya(le,lelt)
     $      ,  rvza(le,lelt),  enra(le,lelt)
      real      lrh(le,lelt),  lru(le,lelt),  lrv(le,lelt)
     $      ,   lrw(le,lelt),  lenr(le,lelt)

      real rhc(lelt), rhuc(lelt), rhvc(lelt), enc(lelt), pprc(lelt)
     $   , uxc(lelt), uyc(lelt) 
      real p_rh(le,lelt)

      real eps, theta1, theta2 
      real minp, minr 
      integer e, nall 

c     Find elem avg
      if(if3d) then 
        call elm_avg_3(rhc, rhuc, rhvc, rhwc, enc
     $        , uxc, uyc,uzc, rha, rvxa, rvya, rvza, enra) ! 
        call get_pres3(pprc,rhc,rhuc,rhvc,rhwc,enc,gama, nelt) ! 
      else
        call elm_avg(x_c, y_c, rhc, rhuc, rhvc, enc, uxc, uyc
     $                     ,  rha,  rvxa,  rvya, enra) ! ck
        call get_pres(pprc,rhc,rhuc,rhvc,enc,gama, nelt) ! 
      endif

      do e = 1,nelt
c     Define eps 
        eps = 1.e-13
        eps = min (eps,rhc(e))
        eps = min (eps,pprc(e))

c     Limiting density 
        call pos_rh(p_rh,rha,rhc,eps,e) ! 
c       write(6,*) 'rh' ,e

c     Find t for limiting press
        call pos_pr(theta2,p_rh,rvxa,rvya,rvza,enra
     $                    , rhc,rhuc,rhvc,rhwc,enc,gama,eps,e) !  ! 
c       write(6,*) 'theta 2 ' ,theta2,'e',e 

c     Limiting everything based on theta2 
        call lm_all( lrh,  lru,  lrv,  lrw, lenr
     $            , p_rh, rvxa, rvya, rvza, enra
     $            , rhc, rhuc, rhvc, rhwc, enc,theta2,e) ! 
c       write(6,*) 'round ' , e 

      enddo 

      nall = nx1*ny1*nz1*nelt
      if(if3d) then
        call copy_all ( rha, rvxa, rvya, rvza, enra
     $                , lrh,  lru,  lrv,  lrw, lenr,nall) 
      else
        call copy_all2( rha, rvxa, rvya,  enra
     $                , lrh,  lru,  lrv,  lenr,nall) 
      endif 

      return
      end
c-----------------------------------------------------------------------
      subroutine pos_pr(theta2,p_rh,rvxa,rvya,rvza,enra
     $                        , rhc,rhuc,rhvc,rhwc,enc,gama,eps,e) !  ! 
      include 'SIZE'
      include 'TOTAL'
      parameter(lg=lx1)
      parameter(lge=lg*ly1*lz1)
      parameter(le=lx1*ly1*lz1)
      parameter(lt=lx1*ly1*lz1*lelt)
      parameter(l3f= lx1*lz1*6*lelt)

      real p_rh(le,1),rvxa(le,1),rvya(le,1),rvza(le,1)  
     $   , enra(le,1) 
      real rhc(1),rhuc(1),rhvc(1),rhwc(1),enc(1) 
      real rhgx(lge), rhgy(lge), rhgz(lge)
     $   , rugx(lge), rugy(lge), rugz(lge)
     $   , rvgx(lge), rvgy(lge), rvgz(lge)
     $   , rwgx(lge), rwgy(lge), rwgz(lge)
     $   , engx(lge), engy(lge), engz(lge)

      real ptmp(lge)

      real tx(lge), ty(lge), tz(lge), mnx,mny,mnz 

      real gama, eps, theta2
      integer i, e

c     Expand in three directions 
      call intp_expd(rhgx,rhgy,rhgz,lg,lx1,p_rh,if3d,0,e)  ! 
      call intp_expd(rugx,rugy,rugz,lg,lx1,rvxa,if3d,0,e)  ! 
      call intp_expd(rvgx,rvgy,rvgz,lg,lx1,rvya,if3d,0,e)  ! 
      if(if3d) then
        call intp_expd(rwgx,rwgy,rwgz,lg,lx1,rvza,if3d,0,e)  ! 
      endif
      call intp_expd(engx,engy,engz,lg,lx1,enra,if3d,0,e)  ! 

c     Get pressure, Gauss in x 
c     write(6,*) 'lge',lge
      if(if3d) then
        call get_pres3(ptmp,rhgx,rugx,rvgx,rwgx,engx,gama,lge) ! 
      else 
        call get_pres (ptmp,rhgx,rugx,rvgx,engx,gama,lge) ! 
c       write(6,*) 'back from get pres 2 ' 
      endif

      do i=1,lge
          if(ptmp(i) .gt. eps) then
            tx(i) = 1.
          else 
           if(if3d) then
            call get_t_quad3(tx(i)
     $                   ,rhgx(i),rugx(i),rvgx(i),rwgx(i),engx(i)
     $                   ,rhc(e),rhuc(e),rhvc(e),rhwc(e),enc(e)  
     $                   , gama, eps) 
           else
            call get_t_quad2(tx(i)
     $                   ,rhgx(i),rugx(i),rvgx(i),engx(i)
     $                   ,rhc(e),rhuc(e),rhvc(e),enc(e)  
     $                   , gama, eps) 
           endif
          endif
      enddo  ! tx done 

c    Gauss in y 
      if(if3d) then
        call get_pres3(ptmp,rhgx,rugx,rvgx,rwgx,engx,gama,lge) ! 
      else 
        call get_pres(ptmp,rhgx,rugx,rvgx,engx,gama,lge) ! 
      endif
      do i=1,lge
          if(ptmp(i) .gt. eps) then
            ty(i) = 1.
          else 
           if(if3d) then
            call get_t_quad3(ty(i)
     $                   ,rhgy(i),rugy(i),rvgy(i),rwgy(i),engy(i)
     $                   ,rhc(e),rhuc(e),rhvc(e),rhwc(e),enc(e)  
     $                   , gama, eps) 
           else
            call get_t_quad2(ty(i)
     $                   ,rhgy(i),rugy(i),rvgy(i),engy(i)
     $                   ,rhc(e),rhuc(e),rhvc(e),enc(e)  
     $                   , gama, eps) 
           endif
          endif
      enddo  ! ty done 

c    Gauss in z 
      if(if3d) then 
        call get_pres3(ptmp,rhgz,rugz,rvgz,rwgz,engz,gama,lge) ! 
        do i=1,lge
            if(ptmp(i) .gt. eps) then
              tz(i) = 1.
            else 
              if(if3d) then
               call get_t_quad3(tz(i)
     $                      ,rhgz(i),rugz(i),rvgz(i),rwgz(i),engz(i)
     $                      ,rhc(e),rhuc(e),rhvc(e),rhwc(e),enc(e)  
     $                      , gama, eps) 
              else
               call get_t_quad2(tz(i)
     $                      ,rhgz(i),rugz(i),rvgz(i),engz(i)
     $                      ,rhc(e),rhuc(e),rhvc(e),enc(e)  
     $                      , gama, eps) 
              endif

            endif
        enddo  ! tz done 
      endif

      mnx = vlmin(tx,lge) 
      mny = vlmin(ty,lge) 
      mnz = vlmin(tz,lge) 

      theta2 = min(mnx,mny) 
      if(if3d) theta2 = min(theta2,mnz)  ! DONE 

      return
      end
c-----------------------------------------------------------------------
      subroutine get_t_quad2(tx
     $                   ,rhgx,rugx,rvgx,engx
     $                   ,rhc,rhuc,rhvc,enc  
     $                   , gama, eps) 
c     Solve quadratic equation 

      real t, eps, gama, g1 
     $   ,    rh, ru, rv, rw, en 
     $   ,   rhc, ruc, rvc, rwc, enc 
      real a, b, c, dl, t1, t2 

      g1 = gama - 1.
      a =  enc*rhc + rh*en - rhc*en - rc*enc
     $   - .5*(ru*ru + rv*rv + ruc*ruc + rvc*rvc)

      b =  -2.*enc*rhc + rhc*en + rc*enc
     $   + ruc*ruc + rvc*rvc + eps*rhc/g1 - eps*rh/g1 

      c =  -.5*(ruc*ruc + rvc*rvc) - eps*rhc/g1

      dl = b*b - 4.*a*c 
      if(dl .gt. 1e-15) then
          dl = sqrt(dl) 
      else
          write(6,*) 'In solving quad eq., detm. < 0',dl  
          call exitt
      endif

c     Formula 1 
      t1 = (-1.*b + dl)/(2.*a) 
      t2 = (-1.*b - dl)/(2.*a) 

c     Formula 2 
c     t1 = 2.*c/(-1.*b - dl)
c     t2 = 2.*c/(-1.*b + dl)

      if(t1.gt.1e-16) then  ! pick the positive one 
          t = t1
      else if(t2.gt.1e-16) then
          t = t2
      else
          write(6,*) '2 roots are negative, quit',t1,t2 
          call exitt
      endif 

      return
      end
c-----------------------------------------------------------------------
      subroutine get_t_quad3(tx
     $                   ,rhgx,rugx,rvgx,rwgx,engx
     $                   ,rhc,rhuc,rhvc,rhwc,enc  
     $                   , gama, eps) 
c     Solve quadratic equation 

      real t, eps, gama, g1 
     $   ,    rh, ru, rv, rw, en 
     $   ,   rhc, ruc, rvc, rwc, enc 
      real a, b, c, dl, t1, t2 

      g1 = gama - 1.
      a =  enc*rhc + rh*en - rhc*en - rc*enc
     $   - .5*(ru*ru + rv*rv + rw*rw + ruc*ruc + rvc*rvc + rwc*rwc)

      b =  -2.*enc*rhc + rhc*en + rc*enc
     $   + ruc*ruc + rvc*rvc + rwc*rwc + eps*rhc/g1 - eps*rh/g1 

      c =  -.5*(ruc*ruc + rvc*rvc + rwc*rwc) - eps*rhc/g1

      dl = b*b - 4.*a*c 
      if(dl .gt. 1e-15) then
          dl = sqrt(dl) 
      else
          write(6,*) 'In solving quad eq., detm. < 0',dl  
          call exitt
      endif

c     Formula 1 
      t1 = (-1.*b + dl)/(2.*a) 
      t2 = (-1.*b - dl)/(2.*a) 

c     Formula 2 
c     t1 = 2.*c/(-1.*b - dl)
c     t2 = 2.*c/(-1.*b + dl)

      if(t1.gt.1e-16) then  ! pick the positive one 
          t = t1
      else if(t2.gt.1e-16) then
          t = t2
      else
          write(6,*) '2 roots are negative, quit',t1,t2 
          call exitt
      endif 

      return
      end
c-----------------------------------------------------------------------
      subroutine lm_all( lrh,  lru,  lrv,  lrw, lenr
     $                 , prh, rvxa, rvya, rvza, enra
     $                 , rhc, rhuc, rhvc, rhwc, enc,theta2,e) ! 
      include 'SIZE'
      include 'TOTAL'
      parameter(le=lx1*ly1*lz1)

      real lrh(le,1), lru(le,1), lrv(le,1), lrw(le,1) 
     $   , lenr(le,1)
      real prh(le,1), rvxa(le,1), rvya(le,1), rvza(le,1) 
     $   , enra(le,1) 
      real rhc(1), rhuc(1), rhvc(1), rhwc(1), enc(1) 

      real theta2
      integer i, e

      do i=1,le
         lrh(i,e)  = theta2*( prh(i,e) -  rhc(e)) +  rhc(e) 
         lru(i,e)  = theta2*(rvxa(i,e) - rhuc(e)) + rhuc(e) 
         lrv(i,e)  = theta2*(rvya(i,e) - rhvc(e)) + rhvc(e) 
         if(if3d) then
           lrw(i,e)  = theta2*(rvza(i,e) - rhwc(e)) + rhwc(e) 
         endif
         lenr(i,e) = theta2*(enra(i,e) -  enc(e)) +  enc(e) 
      enddo 

      return
      end
c-----------------------------------------------------------------------
      subroutine elm_avg_3(rhc, rhuc, rhvc, rhwc, enc
     $           , uxc, uyc,uzc, rha, rvxa, rvya, rvza, enra) ! 
c     There is another routine earlier called elm_avg3 
      include 'SIZE'
      include 'TOTAL'
c   only local, works for MPI 
c
c   Get element averaged values
c     Input : rh(le,nelt), rvx(le,nelt), rvy(le,nelt), enr(le,nelt)
c     Output: rhc(nelt), rhvc(nelt), rhvc(nelt), enc(nelt)
c     
      parameter(le = lx1*ly1*lz1)
      parameter(lf = lx1*lz1*2*ldim*lelt)
      real      rhc(1), rhuc(1), rhvc(1), rhwc(1), enc(1)
     $        , uxc(1), uyc(1), uzc(1) 
     $        , rh(le,1),rvx(le,1),rvy(le,1),rvz(le,1),enr(le,1)
      real      ones(le), tmp(le) 

      nxz    = nx1*nz1
      ne     = nx1*ny1*nz1
      nfaces = 2*ndim 
      nf = nx1*nz1*2*ndim*nelt

      do ie = 1, nelt
          call cfill   (ones,1.,ne)
          call col2    (ones,bm1(1,1,1,ie),ne)
          sum1 = vlsum (ones,ne) ! volume for elem #ie
          call col3    (tmp,bm1(1,1,1,ie),rh(1,ie),ne)
          sum4 = vlsum (tmp ,ne) ! vol. sum for rh
          call col3    (tmp,bm1(1,1,1,ie),rvx(1,ie),ne)
          sum5 = vlsum (tmp ,ne) ! vol. sum for rvx
          call col3    (tmp,bm1(1,1,1,ie),rvy(1,ie),ne)
          sum6 = vlsum (tmp ,ne) ! vol. sum for rvy
          call col3    (tmp,bm1(1,1,1,ie),rvz(1,ie),ne)
          sum2 = vlsum (tmp ,ne) ! vol. sum for rvz
          call col3    (tmp,bm1(1,1,1,ie),enr(1,ie),ne)
          sum7 = vlsum (tmp ,ne) ! vol. sum for enr
         
          rhc (ie) = sum4/sum1
          rhuc(ie) = sum5/sum1
          rhvc(ie) = sum6/sum1
          rhwc(ie) = sum2/sum1
          enc (ie) = sum7/sum1

          uxc (ie) = rhuc(ie)/rhc(ie)
          uyc (ie) = rhvc(ie)/rhc(ie)
          uzc (ie) = rhwc(ie)/rhc(ie)
      enddo 

c     if(nid.eq.0) write(6,*) 'done :: finding elem avg3' 

      return
      end
c-----------------------------------------------------------------------
      subroutine pos_rh(prh,rho,rhc,eps,e) ! 
      include 'SIZE'
      include 'TOTAL'
      parameter(le=lx1*ly1*lz1)
      parameter(lt=lx1*ly1*lz1*lelt)
      parameter(l3f= lx1*lz1*6*lelt)

      real rho(le,1), prh(le,1)  
      real rhc(1)

      real eps, theta1, minrh 
      integer i, e

c     Find min on expanded set of points, (GLL x G)U(G x GLL)
c     write(6,*) 'e ',  e 
c     write(6,*) 'rhc(e) ',  rhc(e) 
      minrh = 1.e-13

      call mn_sij(minrh,rho,e) !

c     Find theta
c     write(6,*) 'eps ', eps
      theta1 = min(1.,(rhc(e) - eps)/(rhc(e) - minrh) )  
      
      do i = 1,nx1*ny1*nz1
         prh(i,e) = theta1*(rho(i,e) - rhc(e)) + rhc(e) 
      enddo 
c     write(6,*) 'rho limit finished' 
c     write(6,*) 'Theta 1 found ', theta1 

      return
      end
c-----------------------------------------------------------------------
      subroutine mn_sij(mnu,u,e) ! for elem # e 
c     Find min on expanded set of points, (GLL x G)U(G x GLL)
      include 'SIZE'
      include 'TOTAL'
      parameter(lg=lx1) ! gauss order + 1 
      parameter(le=lx1*ly1*lz1)
      parameter(lge=lx1*lg*lz1) ! in case lg ~= lx1 

      real u(le,lelt) ! 2 or 3 D 
      real jgz(lx1,ly1,lg) ! Gauss z 
     $   , jgy(lx1,lg,ly1) ! Gauss y 
     $   , jgx(lg,lx1,ly1) ! Gauss x 

      real mnu, mux, muy, muz 
      integer i, e
      integer ng, nl 

      mnu = 1.e-13
      ng = lg ! Gauss order + 1 
      nl = nx1 

      mux = vlmin(u,le) 
c     problem here 
      call intp_expd(jgx,jgy,jgz,ng,nl,u,if3d,0,e) ! 0 = forward

      mux = vlmin(jgx,lge) 
      muy = vlmin(jgy,lge) 
      muz = vlmin(jgz,lge) 

      mnu = min(mux,mnu) 
      mnu = min(muy,mnu) 
      if(if3d) then 
          mnu = min(muz,mnu) 
      endif
c     write(6,*) mnu,mux,muy,muz 

      if(mnu .lt. 1.e-13) then 
          write(6,*) 'Min rho less than 1e-13',mnu
c         call exitt 
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine intp_expd(jgx,jgy,jgz,ng,nl,u,if3d,idir,e) 
c     
c     Interpolate from (nl x nl x nl) to (ng x nl x nl)
c                                        (nl x ng x nl)
c                                        (nl x nl x ng)
c 
      include 'SIZE'

      parameter (ldg=lxd**3,lwkd=4*lxd*lxd)
      common /dgrad/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real jgl,jgt

      real    jgx(1),jgy(1),jgz(1) 
      real    u(1) 
      integer e, idir 
      integer ng,nl
      logical if3d

      call get_int_ptr (i,nl,ng) ! from nl to ng 
c
      if (idir.eq.0) then
         call local_intp(jgx,jgy,jgz,ng,u,nl,jgl(i),jgt(i)
     $                  , e ,if3d)
      else
         call local_intp(jgx,jgy,jgz,nl,u,ng,jgl(i),jgt(i)
     $                  , e ,if3d)
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine local_intp(jgx,jgy,jgz,nl,u,ng,jg,jgt
     $                  , e ,if3d)
c     
      include 'SIZE'

      real    jgx(1),jgy(1),jgz(1) 
      real    jg(1),jgt(1) 
      real    u(1) 
      integer e
      integer ng,nl
      logical if3d

      if (if3d) then
         call local_intp3(jgx,jgy,jgz,ng,u,nl,jg,jgt
     $                  , e )
      else
         call local_intp2(jgx,jgy,    ng,u,nl,jg,jgt
     $                  , e )
c        write(6,*) 'Local interp 2' 
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine local_intp3(jgx,jgy,jgz,ng,u,nl,jg,jgt
     $                 ,e)
c     Output: jgx,jgy,jgz      Input:ng,u,nl,jg,jgt
c     only u is global, meaning have more than 1 elem 
      real jgx(ng,nl,nl),jgy(nl,ng,nl),jgz(nl,nl,ng)
      real u (nl,nl,nl,1)
      real jg(ng,nl),jgt(nl,ng)
      integer e
c
      n1 = nl 
      n2 = n1*n1
      m1 = ng 
      m2 = m1*m1
c
      call mxm(jg,m1,u(1,1,1,e),n1,jgx,n2)
      do k=1,n1
          call mxm(u(1,1,k,e),n1,jgt,n1,jgy,m1)
      enddo
      call mxm(u(1,1,1,e),n2,jgt,n1,jgz,m1)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine local_intp2(jgx,jgy,ng,u,nl,jg,jgt
     $                 ,e)
c     Output: jgx,jgy,jgz      Input:ng,u,nl,jg,jgt
c     only u is global, meaning have more than 1 elem 
      real jgx(ng,nl),jgy(nl,ng)
      real u (nl,nl,1,1)
      real jg(ng,nl),jgt(nl,ng)
      integer e
      logical if3d
c
      n1 = nl 
      m1 = ng 
c
      call mxm(jg,m1,u(1,1,1,e),n1,jgx,n1)
      call mxm(u(1,1,1,e),n1,jgt,n1,jgy,m1)
c
      return
      end
c-----------------------------------------------------------------------
c------ Mon Aug 10 21:57:05 CDT 2015
c------ Change get_ngh, use current handle, just fill surface arrays
c------ with the wanted value 
c------ Thu Sep 10 09:52:40 CDT 2015
c-----------------------------------------------------------------------
      subroutine slopelim2d(rha,rvxa,rvya,rvza,enra) ! 2
      include 'SIZE'
      include 'TOTAL'
      include 'DG'
      include 'DGUSE'
      parameter(le=lx1*ly1*lz1)
      parameter(l3f= lx1*lz1*6*lelt)
      real      rha(le,lelt),  rvxa(le,lelt),  rvya(le,lelt)
     $      ,  rvza(le,lelt),  enra(le,lelt)
      common /tdlimr/   lrh(lt), lrvx(lt), lrvy(lt), lrvz(lt)
     $                , lenr(lt), lmp(lt)
      real       lrh, lrvx, lrvy, lrvz, lenr, lmp
c     common /tdlimi/  ngh(6,lelt)
c     integer*8  ngh
      common /tdlimi/  ngb(4,lelt)
      logical ngb

      real rhc(lelt), rhuc(lelt), rhvc(lelt), enc(lelt), pprc(lelt)
     $   , uxc(lelt), uyc(lelt) 
      real x_c(lelt), y_c(lelt)
      real vcod(2,4,lelt) ! (x,y),(1,2,3,4),(nelt), (3,8,nelt) for 3d
      real ccod(2,4,lelt) 
      real arf(4,lelt)    ! (1,2,3,4),(nelt)
      real rhf(lf), rxf(lf), ryf(lf), enf(lf) 
      real brhc(4,lelt), buc(4,lelt), bvc(4,lelt)
     $   , bprc(4,lelt) 
      real prh(lf), pux(lf), puy(lf), pen(lf), prsu(lf) 
      real brh3(l3f),brx3(l3f),bry3(l3f), bpr3(l3f) ! match area def.

      real vrh(4,2,lelt), vux(4,2,lelt)
     $   , vuy(4,2,lelt), vpr(4,2,lelt) 
      real prpx(4,lelt), prpy(4,lelt), pppx(4,lelt), pppy(4,lelt)
     $   , pupx(4,lelt), pupy(4,lelt), pvpx(4,lelt), pvpy(4,lelt)
      real prpxc(lelt), prpyc(lelt), pupxc(lelt), pupyc(lelt)
     $   , pvpxc(lelt), pvpyc(lelt), pppxc(lelt), pppyc(lelt)
      real lrxc(lelt), lryc(lelt), luxc(lelt), luyc(lelt) 
     $   , lvxc(lelt), lvyc(lelt), lpxc(lelt), lpyc(lelt) 
      real drh(le,lelt),du(le,lelt),dv(le,lelt),dp(le,lelt)

      nf=nx1*nz1*2*ndim*nelt !total number of points on faces
      n=nx1*ny1*nz1*nelt

      call lg_false(ngb,4*nelt)  ! init to .false. 
c    1. Find element avg. coord., rh, ru, rv, en 
      call elm_avg(x_c, y_c, rhc, rhuc, rhvc, enc, uxc, uyc
     $                     ,  rha,  rvxa,  rvya, enra) ! ck
      call get_pres(pprc,rhc,rhuc,rhvc,enc,gama, nelt) ! ck, lelt vector
c                                                      , on each proc 
c    2. Primitive variables at vertices, deal with boundary
      call vert_coor(vcod)     !  ck, Get vertices coord
      call cntr_coor(ccod,x_c,y_c,vcod,ngb) ! ck, Get neighbor center coor
      call area_face(arf,x_c,y_c,ccod,vcod) ! ck, area, for each face 

      call full2face_all2(rhf, rxf, ryf, enf
     $                  , rha,rvxa,rvya,enra) ! ck 

      call copy_all2(brh,brx,bry,ben
     $             , rhf,rxf,ryf,enf,nf)
      call bc_cons_eval(brh, 1) ! ck, only on face though
      call bc_cons_eval(brx, 2) ! inter the same as rhf
      call bc_cons_eval(bry, 3) ! 
      call bc_cons_eval(ben, 5)

      call get_pres(bpr ,brh,brx ,bry ,ben,gama, nf) ! ck, get pressure 

      call get_ctr_b(brhc,buc, bvc,bprc   ! get center values, 
     $             , rhc ,uxc, uyc,pprc   ! with ghost element 
     $             , brh, brx, bry, bpr,ngb) ! ck, prim var 

      call gs_add_bc_all2( rhf, rxf, ryf, enf
     $                   , brh, brx, bry, ben)  ! ck, a + b
      call cmult_all2 (rhf,rxf,ryf,enf,0.5, nf)   ! ck,(a + b)/ 2

      call cons2prim_2( prh, pux, puy, pen, prsu
     $               ,  rhf, rxf, ryf, enf,gama) ! ck, prim  

      call vt_prim(vrh, vux, vuy, vpr
     $           , prh, pux, puy, prsu)  ! ck 

c    3. Face gradient 
      call fac_grad(prpx,prpy,pupx,pupy,pvpx,pvpy,pppx,pppy
     $            ,  vrh, vux, vuy, vpr, rhc,uxc, uyc, pprc
     $            , brhc,buc,bvc,bprc,arf,x_c,y_c,ccod,vcod ) ! ! ck,

c    4. Center gradient 
      call ctr_grad(prpxc, prpyc, pupxc, pupyc, pvpxc, pvpyc
     $            , pppxc, pppyc, prpx, prpy, pupx, pupy
     $            , pvpx, pvpy, pppx, pppy, arf ) ! ! ck,

c    5. Limited gradient 
      eps = 1e-8
      call lim_grad(lrxc, lryc, luxc, luyc, lvxc, lvyc
     $            , lpxc, lpyc, prpxc, prpyc, pupxc, pupyc
     $            , pvpxc, pvpyc, pppxc, pppyc,ngh,eps) ! ! ck,

c    6. Construct difference term
      call diff_sl(drh, du, dv, dp
     $          , lrxc, lryc, luxc, luyc
     $          , lvxc, lvyc, lpxc, lpyc,x_c,y_c ) ! ! ck, 2
c    7. Conservative variabls & 
c    8. Check for negative densities and pressure 
      call limt_cons(lrh, lrvx, lrvy, lenr, lmp
     $             , drh,   du,   dv,   dp
     $             , rhc, rhuc, rhvc, enc  
     $             , uxc, uyc, pprc, gama) ! ! 
c    
      call copy_all2(rha,rvxa,rvya,enra
     $             , lrh,lrvx,lrvy,lenr,n)
      return
      end
c-----------------------------------------------------------------------
      subroutine get_pres(pr, rh, ru, rv, en, gama, n)  ! 2
c         Input : rh, ru, rv, en, gama, n
c        Output : pr 
      real pr(1), rh(1), ru(1), rv(1), en(1)
      real gama
      integer n, i
      do i=1,n
          pr(i) = (gama-1)*(en(i) 
     $            - .5*(ru(i)*ru(i) + rv(i)*rv(i))/rh(i))
      enddo 

      if(nid.eq.0) write(6,*) 'done :: return presu ' 

      return
      end
c-----------------------------------------------------------------------
      subroutine get_pres3(pr, rh, ru, rv, rw, en, gama, n)  ! 2
c         Input : rh, ru, rv, en, gama, n
c        Output : pr 
      real pr(1), rh(1), ru(1), rv(1), en(1)
      real gama
      integer n, i
      do i=1,n
          pr(i) = (gama-1)*(en(i) 
     $    - .5*(ru(i)*ru(i) + rv(i)*rv(i) + rw(i)*rw(i))/rh(i))
      enddo 

      if(nid.eq.0) write(6,*) 'done :: return presu 3 ' 
      return
      end
c-----------------------------------------------------------------------
      subroutine limt_cons(lrh, lrvx, lrvy, lenr, lmp
     $                   , drh,   du,   dv,   dp
     $                   , rhc, rhuc, rhvc, enc 
     $                   , uxc, uyc, pprc, gama)  ! !
      include 'SIZE'
      include 'TOTAL'
      include 'DG'
c
c   Get limited conservative variables
c     AND check negative densities and pressures 
c     Input : drh, du, dv, dp, rhc, rhuc, rhvc, enc 
c    Output : lrh, lrvx, lrvy, lenr, lmp
c                                                                   
      parameter(le = lx1*ly1*lz1)
      parameter(lf = lx1*lz1*2*ldim*lelt)
      real drh(le,1), du(le,1), dv(le,1), dp(le,1)
      real de (le,lelt)  ! temp, used only here
      real rhc(lelt), uxc(lelt), uyc(lelt), pprc(lelt)
     $   , rhuc(lelt), rhvc(lelt), enc(lelt)
      real lrh(le,lelt), lrvx(le,lelt), lrvy(le,lelt)
     $   , lenr(le,lelt), lmp(le,lelt) 
      real gama
      real tol

      nxz    = nx1*nz1
      nfaces = 2*ndim
      nver = 2**(ndim)
      nf = nx1*nz1*2*ndim*nelt

      tol = 1e-9
      do ie = 1, nelt
      do i = 1,le
          lrh(i,ie) = rhc(ie) + drh(i,ie)
          do while(lrh(i,ie) .lt. tol) ! negative density 
              write(6,*) 'Correcting neg. density'
              write(6,*) 'i = ', i,'ie = ', ie 
              drh(i,ie) = drh(i,ie)*0.5
              lrh(i,ie) = rhc(ie) + drh(i,ie)
          enddo 
      enddo 
      enddo 
 
      do ie = 1, nelt
      do i = 1,le
          lrvx(i,ie) = rhuc(ie) + drh(i,ie)*uxc(ie)
     $                          + du (i,ie)*rhc(ie)
          lrvy(i,ie) = rhvc(ie) + drh(i,ie)*uyc(ie)
     $                          + dv (i,ie)*rhc(ie)
          de  (i,ie) = dp(i,ie)/(gama - 1.) 
     $         + .5*drh(i,ie)*(uxc(ie)**2 + uyc(ie)**2)
     $         + rhc(ie)*(uxc(ie)*du(i,ie) + uyc(ie)*dv(i,ie) )
          lenr(i,ie) = enc(ie) + de(i,ie)
          lmp (i,ie) = (gama - 1.)
     $         * ( -0.5*(lrvx(i,ie)**2 + lrvy(i,ie)**2)/lrh(i,ie)
     $            + lenr(i,ie) )
!!!      
          if(lmp(i,ie) .lt. tol) then ! negative pressure
              write(6,*) 'Correcting neg. pressure'
              write(6,*) 'i = ', i,'ie = ', ie
              lmp(i,ie) = enc(ie)
          endif ! Problem being after correction the energy term is 
c         still gonna produce nagative pressure b/c energy is not 
c         changed
      enddo 
      enddo 

      if(nid.eq.0) write(6,*) 'done :: limited solution ' 
      return
      end
c-----------------------------------------------------------------------
      subroutine diff_sl(drh, du, dv, dp
     $                , lrxc, lryc, luxc, luyc
     $                , lvxc, lvyc, lpxc, lpyc,x_c,y_c ) ! ! 2
      include 'SIZE'
      include 'TOTAL'
      include 'DG'
c
c   Get difference terms ( used as blocks for cons variables)
c     Input : lrxc, lryc, luxc, luyc, lvxc, lvyc, lpxc, lpyc 
c    Output : drh, du, dv, dp
c                                                                   
      parameter(le = lx1*ly1*lz1)
      parameter(lf = lx1*lz1*2*ldim*lelt)
      real lrxc(lelt), lryc(lelt), luxc(lelt), luyc(lelt) 
     $   , lvxc(lelt), lvyc(lelt), lpxc(lelt), lpyc(lelt) 
      real x_c(1), y_c(1)
      real drh(le,1), du(le,1), dv(le,1), dp(le,1)

      real eps

      nxz    = nx1*nz1
      nfaces = 2*ndim
      nver = 2**(ndim)
      nf = nx1*nz1*2*ndim*nelt

      do ie = 1, nelt
      do i =1, le
          x = xm1(i,1,1,ie)
          y = ym1(i,1,1,ie)
          drh(i,ie) = lrxc(ie)*(x - x_c(ie))
     $              + lryc(ie)*(y - y_c(ie))
          du (i,ie) = luxc(ie)*(x - x_c(ie))
     $              + luyc(ie)*(y - y_c(ie))
          dv (i,ie) = lvxc(ie)*(x - x_c(ie))
     $              + lvyc(ie)*(y - y_c(ie))
          dp (i,ie) = lpxc(ie)*(x - x_c(ie))
     $              + lpyc(ie)*(y - y_c(ie))
      enddo 
      enddo 
 
      if(nid.eq.0) write(6,*) 'done :: difference terms ' 
      return
      end
c-----------------------------------------------------------------------
      subroutine lim_grad(lrxc, lryc, luxc, luyc, lvxc, lvyc
     $            , lpxc, lpyc, prpxc, prpyc, pupxc, pupyc
     $            , pvpxc, pvpyc, pppxc, pppyc, ngb,eps) ! ! 2
      include 'SIZE'
      include 'TOTAL'
      include 'DG'
c
c   Get limited gradient 
c     Input : prpxc, prpyc, pupxc, pupyc, pvpxc, pvpyc, pppxc, pppyc
c    Output : lrxc, lryc, luxc, luyc, lvxc, lvyc, lpxc, lpyc 
c                                                                   
      parameter(le = lx1*ly1*lz1)
      parameter(lf = lx1*lz1*2*ldim*lelt)
      real lrxc(lelt), lryc(lelt), luxc(lelt), luyc(lelt) 
     $   , lvxc(lelt), lvyc(lelt), lpxc(lelt), lpyc(lelt) 

      logical ngb(4,lelt)
      real prpxc(lelt), prpyc(lelt), pupxc(lelt), pupyc(lelt)
     $   , pvpxc(lelt), pvpyc(lelt), pppxc(lelt), pppyc(lelt)
      real eps
      real g1, g2, g3, g4, gsm 
c obtain g1,g2,g3,g4 from gs_op 
c  Thu Sep 10 12:52:58 CDT 2015
      real   fdxrh(lx1*lz1,2*ldim,lelt), fdxrhs(lx1*lz1,2*ldim,lelt)
     $     , fdxpr(lx1*lz1,2*ldim,lelt), fdxprs(lx1*lz1,2*ldim,lelt)
     $     , fdxux(lx1*lz1,2*ldim,lelt), fdxuxs(lx1*lz1,2*ldim,lelt)
     $     , fdxuy(lx1*lz1,2*ldim,lelt), fdxuys(lx1*lz1,2*ldim,lelt)
     $     , fdyrh(lx1*lz1,2*ldim,lelt), fdyrhs(lx1*lz1,2*ldim,lelt)
     $     , fdypr(lx1*lz1,2*ldim,lelt), fdyprs(lx1*lz1,2*ldim,lelt)
     $     , fdyux(lx1*lz1,2*ldim,lelt), fdyuxs(lx1*lz1,2*ldim,lelt)
     $     , fdyuy(lx1*lz1,2*ldim,lelt), fdyuys(lx1*lz1,2*ldim,lelt)

      nxz    = nx1*nz1
      nfaces = 2*ndim
      nver = 2**(ndim)
      nf = nx1*nz1*2*ndim*nelt
      ns = (nx1*nz1)/2 ! pointer on faces 

      call ctr2face(fdxrh,prpxc) ! from center to face 
      call copy    (fdxrhs,fdxrh,nf)
      call gs_op   (dg_hndl,fdxrh,1,1,0)  ! 1 ==> +
      call sub2    (fdxrh,fdxrhs,nf) ! a + b - a 
      call ctr2face(fdyrh,prpyc) ! from center to face 
      call copy    (fdyrhs,fdyrh,nf)
      call gs_op   (dg_hndl,fdyrh,1,1,0)  ! 1 ==> +
      call sub2    (fdyrh,fdyrhs,nf) ! a + b - a 

      call ctr2face(fdxux,pupxc) ! from center to face 
      call copy    (fdxuxs,fdxux,nf)
      call gs_op   (dg_hndl,fdxux,1,1,0)  ! 1 ==> +
      call sub2    (fdxux,fdxuxs,nf) ! a + b - a 
      call ctr2face(fdyux,pupyc) ! from center to face 
      call copy    (fdyuxs,fdyux,nf)
      call gs_op   (dg_hndl,fdyux,1,1,0)  ! 1 ==> +
      call sub2    (fdyux,fdyuxs,nf) ! a + b - a 

      call ctr2face(fdxuy,pvpxc) ! from center to face 
      call copy    (fdxuys,fdxuy,nf)
      call gs_op   (dg_hndl,fdxuy,1,1,0)  ! 1 ==> +
      call sub2    (fdxuy,fdxuys,nf) ! a + b - a 
      call ctr2face(fdyuy,pvpyc) ! from center to face 
      call copy    (fdyuys,fdyuy,nf)
      call gs_op   (dg_hndl,fdyuy,1,1,0)  ! 1 ==> +
      call sub2    (fdyuy,fdyuys,nf) ! a + b - a 

      call ctr2face(fdxpr,pppxc) ! from center to face 
      call copy    (fdxprs,fdxpr,nf)
      call gs_op   (dg_hndl,fdxpr,1,1,0)  ! 1 ==> +
      call sub2    (fdxpr,fdxprs,nf) ! a + b - a 
      call ctr2face(fdypr,pppyc) ! from center to face 
      call copy    (fdyprs,fdypr,nf)
      call gs_op   (dg_hndl,fdypr,1,1,0)  ! 1 ==> +
      call sub2    (fdypr,fdyprs,nf) ! a + b - a 

      do ie = 1, nelt
          if( ngb(ie,1) .or. ngb(ie,2)
     $   .or. ngb(ie,3) .or. ngb(ie,4) ) then ! ? 
c           write(6,*) 'What happens at boundary' 
c     now do the same thing? The value would be just 0. 
c         else
          endif 
c   -- rho --
            g1 = fdxrh(ns,1,ie)**2 + fdyrh(ns,1,ie)**2
            g2 = fdxrh(ns,2,ie)**2 + fdyrh(ns,2,ie)**2
            g3 = fdxrh(ns,3,ie)**2 + fdyrh(ns,3,ie)**2
            g4 = fdxrh(ns,4,ie)**2 + fdyrh(ns,4,ie)**2
            gsm = g1**2 + g2**2 + g3**2 + g4**2  !! L2 norm sum
            gsm = (gsm)**(1.5)  ! ! for dimension purpose
            w1 = (g2*g3*g4 + eps) / ( gsm + 4.*eps)
            w2 = (g1*g3*g4 + eps) / ( gsm + 4.*eps)
            w3 = (g1*g2*g4 + eps) / ( gsm + 4.*eps)
            w4 = (g1*g2*g3 + eps) / ( gsm + 4.*eps)
            lrxc(ie) = w1*fdxrh(ns,1,ie) + w2*fdxrh(ns,2,ie)
     $               + w3*fdxrh(ns,3,ie) + w4*fdxrh(ns,4,ie)
            lryc(ie) = w1*fdyrh(ns,1,ie) + w2*fdyrh(ns,2,ie)
     $               + w3*fdyrh(ns,3,ie) + w4*fdyrh(ns,4,ie)
c   -- ux --
            g1 = fdxux(ns,1,ie)**2 + fdyux(ns,1,ie)**2
            g2 = fdxux(ns,2,ie)**2 + fdyux(ns,2,ie)**2
            g3 = fdxux(ns,3,ie)**2 + fdyux(ns,3,ie)**2
            g4 = fdxux(ns,4,ie)**2 + fdyux(ns,4,ie)**2
            gsm = g1**2 + g2**2 + g3**2 + g4**2  !! L2 norm sum
            gsm = (gsm)**(1.5)  ! ! for dimension purpose raise to 3
            w1 = (g2*g3*g4 + eps) / ( gsm + 4.*eps)
            w2 = (g1*g3*g4 + eps) / ( gsm + 4.*eps)
            w3 = (g1*g2*g4 + eps) / ( gsm + 4.*eps)
            w4 = (g1*g2*g3 + eps) / ( gsm + 4.*eps)
            luxc(ie) = w1*fdxux(ns,1,ie) + w2*fdxux(ns,2,ie)
     $               + w3*fdxux(ns,3,ie) + w4*fdxux(ns,4,ie)
            luyc(ie) = w1*fdyux(ns,1,ie) + w2*fdyux(ns,2,ie)
     $               + w3*fdyux(ns,3,ie) + w4*fdyux(ns,4,ie)
c   -- uy --
            g1 = fdxuy(ns,1,ie)**2 + fdyuy(ns,1,ie)**2
            g2 = fdxuy(ns,2,ie)**2 + fdyuy(ns,2,ie)**2
            g3 = fdxuy(ns,3,ie)**2 + fdyuy(ns,3,ie)**2
            g4 = fdxuy(ns,4,ie)**2 + fdyuy(ns,4,ie)**2
            gsm = g1**2 + g2**2 + g3**2 + g4**2  !! L2 norm sum
            gsm = (gsm)**(1.5)  ! ! for dimension purpose raise to 3
            w1 = (g2*g3*g4 + eps) / ( gsm + 4.*eps)
            w2 = (g1*g3*g4 + eps) / ( gsm + 4.*eps)
            w3 = (g1*g2*g4 + eps) / ( gsm + 4.*eps)
            w4 = (g1*g2*g3 + eps) / ( gsm + 4.*eps)
            lvxc(ie) = w1*fdxuy(ns,1,ie) + w2*fdxuy(ns,2,ie)
     $               + w3*fdxuy(ns,3,ie) + w4*fdxuy(ns,4,ie)
            lvyc(ie) = w1*fdyuy(ns,1,ie) + w2*fdyuy(ns,2,ie)
     $               + w3*fdyuy(ns,3,ie) + w4*fdyuy(ns,4,ie)
c   -- p  --
            g1 = fdxpr(ns,1,ie)**2 + fdypr(ns,1,ie)**2
            g2 = fdxpr(ns,2,ie)**2 + fdypr(ns,2,ie)**2
            g3 = fdxpr(ns,3,ie)**2 + fdypr(ns,3,ie)**2
            g4 = fdxpr(ns,4,ie)**2 + fdypr(ns,4,ie)**2
            gsm = g1**2 + g2**2 + g3**2 + g4**2  !! L2 norm sum
            gsm = (gsm)**(1.5)  ! ! for dimension purpose
            w1 = (g2*g3*g4 + eps) / ( gsm + 4.*eps)
            w2 = (g1*g3*g4 + eps) / ( gsm + 4.*eps)
            w3 = (g1*g2*g4 + eps) / ( gsm + 4.*eps)
            w4 = (g1*g2*g3 + eps) / ( gsm + 4.*eps)
            lpxc(ie) = w1*fdxpr(ns,1,ie) + w2*fdxpr(ns,2,ie)
     $               + w3*fdxpr(ns,3,ie) + w4*fdxpr(ns,4,ie)
            lpyc(ie) = w1*fdypr(ns,1,ie) + w2*fdypr(ns,2,ie)
     $               + w3*fdypr(ns,3,ie) + w4*fdypr(ns,4,ie)
c         endif
      enddo 
 
      if(nid.eq.0) write(6,*) 'done :: limited gradient ' 
      return
      end
c-----------------------------------------------------------------------
      subroutine ctr_grad(prpxc, prpyc, pupxc, pupyc, pvpxc, pvpyc
     $                  , pppxc, pppyc, prpx, prpy, pupx, pupy
     $                  , pvpx, pvpy, pppx, pppy, arf ) ! !
      include 'SIZE'
      include 'TOTAL'
      include 'DG'
c
c   Get center gradients 
c     Input : prpx, prpy, pupx, pupy, pvpx, pvpy, pppx, pppy, arf
c    Output : prpxc, prpyc, pupxc, pupyc, pvpxc, pvpyc, pppxc, pppyc
c                                                                   
c     vertice should aligh with the order of increasing 
c     vertex number
c
      parameter(le = lx1*ly1*lz1)
      parameter(lf = lx1*lz1*2*ldim*lelt)
      real prpx(4,lelt), prpy(4,lelt), pupx(4,lelt), pupy(4,lelt)
     $   , pppx(4,lelt), pppy(4,lelt), pvpx(4,lelt), pvpy(4,lelt)
      real  arf(4,lelt) 

      real prpxc(lelt), prpyc(lelt), pupxc(lelt), pupyc(lelt)
     $   , pvpxc(lelt), pvpyc(lelt), pppxc(lelt), pppyc(lelt)

      nxz    = nx1*nz1
      nfaces = 2*ndim
      nver = 2**(ndim)
      nf = nx1*nz1*2*ndim*nelt

      k = 1
      do ie = 1, nelt
      arsum = arf(1,ie) + arf(2,ie) + arf(3,ie) + arf(4,ie)
      prpxc(ie) = (arf(1,ie)*prpx(1,ie) + arf(2,ie)*prpx(2,ie) 
     $           + arf(3,ie)*prpx(3,ie) + arf(4,ie)*prpx(4,ie))
     $           / arsum 
      prpyc(ie) = (arf(1,ie)*prpy(1,ie) + arf(2,ie)*prpy(2,ie) 
     $           + arf(3,ie)*prpy(3,ie) + arf(4,ie)*prpy(4,ie))
     $           / arsum 
      pupxc(ie) = (arf(1,ie)*pupx(1,ie) + arf(2,ie)*pupx(2,ie) 
     $           + arf(3,ie)*pupx(3,ie) + arf(4,ie)*pupx(4,ie))
     $           / arsum 
      pupyc(ie) = (arf(1,ie)*pupy(1,ie) + arf(2,ie)*pupy(2,ie) 
     $           + arf(3,ie)*pupy(3,ie) + arf(4,ie)*pupy(4,ie))
     $           / arsum 
      pvpxc(ie) = (arf(1,ie)*pvpx(1,ie) + arf(2,ie)*pvpx(2,ie) 
     $           + arf(3,ie)*pvpx(3,ie) + arf(4,ie)*pvpx(4,ie))
     $           / arsum 
      pvpyc(ie) = (arf(1,ie)*pvpy(1,ie) + arf(2,ie)*pvpy(2,ie) 
     $           + arf(3,ie)*pvpy(3,ie) + arf(4,ie)*pvpy(4,ie))
     $           / arsum 
      pppxc(ie) = (arf(1,ie)*pppx(1,ie) + arf(2,ie)*pppx(2,ie) 
     $           + arf(3,ie)*pppx(3,ie) + arf(4,ie)*pppx(4,ie))
     $           / arsum 
      pppyc(ie) = (arf(1,ie)*pppy(1,ie) + arf(2,ie)*pppy(2,ie) 
     $           + arf(3,ie)*pppy(3,ie) + arf(4,ie)*pppy(4,ie))
     $           / arsum 
c    
      enddo 
 
      if(nid.eq.0) write(6,*) 'done :: cneter gradient ' 
      return
      end
c-----------------------------------------------------------------------
      subroutine fac_grad(prpx,prpy,pupx,pupy,pvpx,pvpy,pppx,pppy
     $            ,  vrh, vux, vuy, vpr, rhc,uxc, uyc, pprc
     $            , brhc,buc,bvc,bprc,arf,x_c,y_c,ccod,vcod ) ! ! 2 
      include 'SIZE'
      include 'TOTAL'
      include 'DG'
c
c   Get face gradients 
c     Input : vrh, vux, vuy, vpr, rhc,uxc, uyc, pprc 
c           , brhc,buc,bvc,bprc
c    Output : prpx, prpy, pupx, pupy, pvpx, pvpy, pppx, pppy
c                                                                   
c           ifc =(1,2,3,4)                                          
c  \1/  -------------------------   \2/                              
c                                                                   
c                                                                   
c     vertice should aligh with the order of increasing 
c     vertex number
      parameter(le = lx1*ly1*lz1)
      parameter(lf = lx1*lz1*2*ldim*lelt)
      real prpx(4,lelt), prpy(4,lelt), pupx(4,lelt), pupy(4,lelt)
     $   , pppx(4,lelt), pppy(4,lelt), pvpx(4,lelt), pvpy(4,lelt)

      real vrh(4,2,lelt), vux(4,2,lelt)
     $   , vuy(4,2,lelt), vpr(4,2,lelt) 
      real rhc(lelt), uxc(lelt), uyc(lelt), pprc(lelt)
      real brhc(4,lelt), buc(4,lelt), bvc(4,lelt), bprc(4,lelt)
      real ccod(2,4,lelt), vcod(2,4,lelt), arf(4,lelt)
      real x_c(1), y_c(1)
      real tmp

      nxz    = nx1*nz1
      nfaces = 2*ndim
      nver = 2**(ndim)
      nf = nx1*nz1*2*ndim*nelt

      k = 1
      do ie = 1, nelt
      xc0 = x_c(ie)
      yc0 = y_c(ie)
      rc0 = rhc(ie)
      uc0 = uxc(ie)
      vc0 = uyc(ie)
      pc0 =pprc(ie)
      do iface = 1, nfaces
          iv1 = iface
          iv2 = iface + 1
          if(iface.eq.4) iv2 = 1
          xc1 = ccod(1,iface,ie)
          yc1 = ccod(2,iface,ie)
          xv1 = vcod(1,iv1,ie)
          yv1 = vcod(2,iv1,ie)
          xv2 = vcod(1,iv2,ie)
          yv2 = vcod(2,iv2,ie)

          rc1 = brhc(iface,ie)
          uc1 = buc (iface,ie)
          vc1 = bvc (iface,ie)
          pc1 = bprc(iface,ie)
          rv1 = vrh (iface,1,ie)
          uv1 = vux (iface,1,ie)
          vv1 = vuy (iface,1,ie)
          pv1 = vpr (iface,1,ie)
          rv2 = vrh (iface,2,ie)
          uv2 = vux (iface,2,ie)
          vv2 = vuy (iface,2,ie)
          pv2 = vpr (iface,2,ie)
c -- r
          tmp =   ((rc1 - rc0)*(yv2-yv1) - (rv2 - rv1)*(yc1 - yc0))
          prpx(iface,ie) = tmp/(2.*arf(iface,ie) )
          tmp = - ((rc1 - rc0)*(xv2-xv1) - (rv2 - rv1)*(xc1 - xc0))
          prpy(iface,ie) = tmp/(2.*arf(iface,ie) )
c -- u          
          tmp =   ((uc1 - uc0)*(yv2-yv1) - (uv2 - uv1)*(yc1 - yc0))
          pupx(iface,ie) = tmp/(2.*arf(iface,ie) )
          tmp = - ((uc1 - uc0)*(xv2-xv1) - (uv2 - uv1)*(xc1 - xc0))
          pupy(iface,ie) = tmp/(2.*arf(iface,ie) )
c -- v          
          tmp =   ((vc1 - vc0)*(yv2-yv1) - (vv2 - vv1)*(yc1 - yc0))
          pvpx(iface,ie) = tmp/(2.*arf(iface,ie) )
          tmp = - ((vc1 - vc0)*(xv2-xv1) - (vv2 - vv1)*(xc1 - xc0))
          pvpy(iface,ie) = tmp/(2.*arf(iface,ie) )
c -- p          
          tmp =   ((pc1 - pc0)*(yv2-yv1) - (pv2 - pv1)*(yc1 - yc0))
          pppx(iface,ie) = tmp/(2.*arf(iface,ie) )
          tmp = - ((pc1 - pc0)*(xv2-xv1) - (pv2 - pv1)*(xc1 - xc0))
          pppy(iface,ie) = tmp/(2.*arf(iface,ie) )
          
      enddo 
      enddo 
 
      if(nid.eq.0) write(6,*) 'done :: face gradient ' 
      return
      end
c-----------------------------------------------------------------------
      subroutine vt_prim(vrh, vux, vuy, vpr ! 2 
     $                 , prh, pux, puy, prsu) 
      include 'SIZE'
      include 'TOTAL'
      include 'DG'
c
c   Get verticies values with ghost element, prim var 
c   - Assumeing boundary values are already added and /2 
c     in the computation of prh, pux, puy, prsu  -
c     Input : prh, pux, puy, prsu ( surface values )
c    Output : vrh, vux, vuy, vpr
c                                                                   
c           ifc =(1,2,3,4)                                          
c  \1/  -------------------------   \2/                              
c                                                                   
c                                                                   
c     vertice always start with lower index to higher index 
      parameter(le = lx1*ly1*lz1)
      parameter(lf = lx1*lz1*2*ldim*lelt)
      real vrh(4,2,lelt), vux(4,2,lelt)
     $   , vuy(4,2,lelt), vpr(4,2,lelt) 
      real      prh(lf), pux(lf), puy(lf), prsu(lf)

      real      ones(lf), tmp, tm1
      nxz    = nx1*nz1
      nfaces = 2*ndim
      nver = 2**(ndim)
      nf = nx1*nz1*2*ndim*nelt

      k = 1
      do ie = 1, nelt
      do iface = 1, nfaces
        if(iface.le.2) then
          vrh(iface,1,ie) = prh(k)
          vux(iface,1,ie) = pux(k)
          vuy(iface,1,ie) = puy(k)
          vpr(iface,1,ie) = prsu(k)
          k = k + nxz - 1 ! nx1 
          vrh(iface,2,ie) = prh(k)
          vux(iface,2,ie) = pux(k)
          vuy(iface,2,ie) = puy(k)
          vpr(iface,2,ie) = prsu(k)
          k = k + 1
        else 
          vrh(iface,2,ie) = prh(k)
          vux(iface,2,ie) = pux(k)
          vuy(iface,2,ie) = puy(k)
          vpr(iface,2,ie) = prsu(k)
          k = k + nxz - 1 ! nx1 
          vrh(iface,1,ie) = prh(k)
          vux(iface,1,ie) = pux(k)
          vuy(iface,1,ie) = puy(k)
          vpr(iface,1,ie) = prsu(k)
          k = k + 1
        endif
      enddo 
      enddo 
 
      if(nid.eq.0) write(6,*) 'done :: vertices prim. variables ' 
      return
      end
c-----------------------------------------------------------------------
      subroutine lift_3d_surf( brh3, brx3, bry3, bpr3
     $                        , brh , brx, bry, bpr) ! from lf to l3f 
      include 'SIZE'
      include 'TOTAL'
      include 'DG'
c
c   Get surface values from lf to l3f 
c     Input : brh , brx, bry, bpr
c    Output : brh3, brx3, bry3, bpr3
c     
      parameter(le = lx1*ly1*lz1)
      parameter(lf = lx1*lz1*2*ldim*lelt)  ! lf == 2/3 l3f in 2D 
      parameter(l3f= lx1*lz1*6*lelt)       !,lf ==     l3f in 3D
      real      brh (lx1,lz1,4,lelt),  brx(lx1,lz1,4,lelt)
     $        , bry (lx1,lz1,4,lelt),  bpr(lx1,lz1,4,lelt)
      real      brh3(lx1,lz1,6,lelt), brx3(lx1,lz1,6,lelt)
     $        , bry3(lx1,lz1,6,lelt), bpr3(lx1,lz1,6,lelt)

      nxz    = nx1*nz1
      nfaces = 2*ndim 
      nf = nx1*nz1*2*ndim*nelt
      n3f = nx1*nz1*6*nelt

c     !-!-! !-!-! !-!-! !-!-! !-!-! !-!-! !-!-! !-!-! 
c     It seems that to use facint_a routine, the array 
c     has to be defined as brh(lx1,lz1,6,lelt) 
c     regardless of 2D or 3D 
c     but brh is in DGUSE...
c     Let's just don't use brh, define something new 
c     Wed Aug 12 16:44:23 CDT 2015
c     !-!-! !-!-! !-!-! !-!-! !-!-! !-!-! !-!-! !-!-! 
      if(if3d) then  
          write(6,*) '3D, no need here.'
      else
        do ie = 1, nelt
        do iface=1,nfaces ! 1 to 4 
          call copy(brh3(1,1,iface,ie), brh(1,1,iface,ie),nxz)
          call copy(brx3(1,1,iface,ie), brx(1,1,iface,ie),nxz)
          call copy(bry3(1,1,iface,ie), bry(1,1,iface,ie),nxz)
          call copy(bpr3(1,1,iface,ie), bpr(1,1,iface,ie),nxz)
        enddo 
        call rzero(brh3(1,1,5,ie),nxz)
        call rzero(brh3(1,1,6,ie),nxz)
        call rzero(brx3(1,1,5,ie),nxz)
        call rzero(brx3(1,1,6,ie),nxz)
        call rzero(bry3(1,1,5,ie),nxz)
        call rzero(bry3(1,1,6,ie),nxz) ! shouldn't matter though
        call rzero(bpr3(1,1,5,ie),nxz)
        call rzero(bpr3(1,1,6,ie),nxz)
        enddo 
      endif
 
      return
      end
c-----------------------------------------------------------------------
      subroutine get_ctr_b(brhc,buc, bvc,bprc   ! ! 2 
     $                   , rhc ,uxc, uyc,pprc  
     $                   , brh , brx, bry, bpr,ngb)
      include 'SIZE'
      include 'TOTAL'
      include 'DG'
c
c   Get vol. center values with ghost element, prim var 
c     Input : rhc, uxc, uyc, pprc (center values) 
c    Output : brhc, buc, bvc, bprc
c     
      parameter(le = lx1*ly1*lz1)
      parameter(lf = lx1*lz1*2*ldim*lelt)
      parameter(l3f= lx1*lz1*6*lelt)
      logical   ngb(4,lelt)
      real      brhc(4,lelt), buc(4,lelt), bvc(4,lelt)
     $        , bprc(4,lelt)
      real      rhc (1),  uxc(1), uyc(1), pprc(1)
      real      brh(lf), brx(lf), bry(lf), bpr(lf)
      real      brh3(l3f), brx3(l3f), bry3(l3f), bpr3(l3f)
      real      ones(lx1*lz1*6*lelt), tmp, tm1
      real      frhc(lx1*lz1,2*ldim,lelt), frhcs(lx1*lz1,2*ldim,lelt)
     $        , fprc(lx1*lz1,2*ldim,lelt), fprcs(lx1*lz1,2*ldim,lelt)
     $        , fuxc(lx1*lz1,2*ldim,lelt), fuxcs(lx1*lz1,2*ldim,lelt)
     $        , fuyc(lx1*lz1,2*ldim,lelt), fuycs(lx1*lz1,2*ldim,lelt)

c     !-!-! 
c     It seems that to use facint_a routine, the array 
c     has to be defined as brh(lx1,lz1,6,lelt) 
c     but brh is in DGUSE...
c     !-!-! 
c     filling 0's at positions where faces 5 and 6 are
      call lift_3d_surf( brh3, brx3, bry3, bpr3
     $                 , brh , brx, bry, bpr) ! ck, from lf to l3f 

      nxz    = nx1*nz1
      nfaces = 2*ndim 
      nf = nx1*nz1*2*ndim*nelt
      n3f = nx1*nz1*6*nelt
      ns = (nx1*nz1)/2 ! pointer on faces 

      call ctr2face(frhc,rhc)
      call copy    (frhcs,frhc,nf)
      call gs_op   (dg_hndl,frhc,1,1,0)  ! 1 ==> +
      call sub2    (frhc,frhcs,nf) ! a + b - a 
      call ctr2face(fprc,pprc)
      call copy    (fprcs,fprc,nf)
      call gs_op   (dg_hndl,fprc,1,1,0)  ! 1 ==> +
      call sub2    (fprc,fprcs,nf) ! a + b - a 
      call ctr2face(fuxc,uxc)
      call copy    (fuxcs,fuxc,nf)
      call gs_op   (dg_hndl,fuxc,1,1,0)  ! 1 ==> +
      call sub2    (fuxc,fuxcs,nf) ! a + b - a 
      call ctr2face(fuyc,uyc)
      call copy    (fuycs,fuyc,nf)
      call gs_op   (dg_hndl,fuyc,1,1,0)  ! 1 ==> +
      call sub2    (fuyc,fuycs,nf) ! a + b - a 

      call cfill(ones,1.,n3f)
      do ie = 1, nelt
      do iface=1,nfaces
         if(ngb(iface,ie)) then ! on boundary, use element number 
            tm1 = facint_a(ones,area,iface,ie)  !?? doubt right 
            tmp = facint_a(brh3,area,iface,ie)  !?? doubt right 
            brhc(iface,ie) =  tmp/tm1                ! rho c
            tmp = facint_a(brx3,area,iface,ie)  !?? doubt right 
            buc (iface,ie) =  tmp/(tm1*brhc(iface,ie)) ! uc 
            tmp = facint_a(bry3,area,iface,ie)  !?? doubt right 
            bvc (iface,ie) =  tmp/(tm1*brhc(iface,ie)) ! vc 
            tmp = facint_a(bpr3,area,iface,ie)  !?? doubt right 
            bprc(iface,ie) =  tmp/tm1                ! pres c
         else 
            brhc(iface,ie) =  frhc(ns,iface,ie)
            bprc(iface,ie) =  fprc(ns,iface,ie)
            buc (iface,ie) =  fuxc(ns,iface,ie)
            bvc (iface,ie) =  fuyc(ns,iface,ie)
         endif
      enddo 
      enddo 
 
      if(nid.eq.0) write(6,*) 'done :: center val. from nghb' 
      return
      end
c-----------------------------------------------------------------------
      subroutine get_ngh(ngh) ! 2+3? not called anywhere
      include 'SIZE'
      include 'TOTAL'
      include 'DG'
c
c   Get neighboring element number 
c     
      parameter(le = lx1*ly1*lz1)
      parameter(lf = lx1*lz1*2*ldim*lelt)
c   real version 
      integer*8 ngh(6,lelt) ! same size as real 
      real      env(le,lelt), renf(lf), enfsv(lf)

      nxz    = nx1*nz1
      nfaces = 2*ndim 
      nf = nx1*nz1*2*ndim*nelt

c   Wed Sep  9 12:55:43 CDT 2015
c this approach will return the local elem. number on another proc.
c if the elem. of interest has a neighbor on another proc. 
      call rzero    (ngh,6*nelt)
      do ie = 1, nelt
c -- Approach 1. -- 
c        call cfill (env(1,ie),real(ie),le) ! element number 
c -- Approach 2. -- 
c   Find everywhere ngh is used, use gllel to get local elem. #
         ig = lglel(ie)  ! global elem number, can fix it? 
         call cfill (env(1,ie),real(ig),le) ! gl element number 
      enddo 
      call full2face(renf,env)
      call copy     (enfsv,renf,nf)      ! a 
      call chsign   (enfsv,nf)           ! - a 
      call gs_op    (dg_hndl,renf,1,1,0) ! 1 ==> + , a + b 
      call add2     (renf,enfsv,nf)      ! (a + b) - a = b 
c     0 on boundary

      k     = 1 
      do ie = 1, nelt
         do iface=1,nfaces
            ngh(iface,ie) = int(renf(k))
            if(abs(renf(k) + 0.) .le. 1e-10 ) then 
c     == 0, on boundary, use element number 
c -- Approach 1. -- 
c              ngh(iface,ie) = ie
c -- Approach 2. -- 
               ig = lglel(ie)  ! global elem number
               ngh(iface,ie) = ig
            endif
            k = k + nxz ! only use the first value on each face
         enddo 
      enddo 
 
      return
      end
c-----------------------------------------------------------------------
      subroutine area_face(arf,x_c,y_c,ccod,vcod) ! 
      include 'SIZE'
      include 'TOTAL'
c   Get area, for each face 
c     Input:  x-c, y-c, ccod(xy,1-4,nelt), vcod(xy,1-4,nelt)
c     Output: arf(1-4,nelt) 
c     Area : 
c 
c       (xv2,yv2) 
c   +-------+-------+                           
c   | ec0   ^  ec1  |                                    
c   |      /|\      |  A = 1/2|(xc1-xc0,yc1-yc0)x(xv2-xv1,yv2-yv1)| 
c   |     / | \     |                           
c   | xc0/__|__\xc1 |                           
c   | yc0\  |  /yc1 |                           
c   |     \ | /     |                           
c   |      \|/      |                           
c   |       v       |                           
c   +-------+-------+                            
c       (xv1,yv1)                           
c
      parameter(le = lx1*ly1*lz1)
c     in 3d ccod(3,6,nelt), vcod(3,8,nelt) 
      real      ccod(2,4,lelt), vcod(2,4,lelt)
      real      x_c(1), y_c(1) 
      real      arf(4,lelt) 
      real      xc0, yc0, xc1, yc1 
     $        , xv2, yv2, xv1, yv1 , tmp

      ne     = nx1*ny1*nz1
      nfaces = 2*ndim 

      do ie = 1, nelt
      xc0 = x_c(ie)
      yc0 = y_c(ie)
      do iface = 1,nfaces
          iv1 = iface
          iv2 = iface + 1
          if(iface.eq.4) iv2 = 1
          xc1 = ccod(1,iface,ie)
          yc1 = ccod(2,iface,ie)
          xv1 = vcod(1,iv1,ie)
          yv1 = vcod(2,iv1,ie)
          xv2 = vcod(1,iv2,ie)
          yv2 = vcod(2,iv2,ie)
          tmp = ((xc1 - xc0)*(yv2-yv1) - (xv2 - xv1)*(yc1 - yc0))
          tmp = abs(tmp)/2.
          arf(iface,ie) = tmp 
      enddo 
      enddo 

      if(nid.eq.0) write(6,*) 'done :: compute area ' 
      return
      end
c-----------------------------------------------------------------------
      subroutine lg_false(larry,n)  ! init to .false. 
      include 'SIZE'
      include 'TOTAL'
      logical larry(1)
      integer i,n 
      do i=1,n 
          larry(i) = .false. 
      enddo 
      if(nid.eq.0) write(6,*) 'done :: lg false return ' 
      return
      end
c-----------------------------------------------------------------------
      subroutine cntr_coor(ccod,x_c,y_c,vcod,ngb) ! 2 
      include 'SIZE'
      include 'TOTAL'
      include 'DG'
c     parameter(lf=lx1*lz1*2*ldim*lelt)
c   Get neighbor center coor
c     Input : x-c, y-c, vcod(xy,1-4,nelt)
c     Output: ccod(xy,1-4,nelt)
c     
c   \4/ 1        2  \3/             .            3, n     
c     2 +---------+ 2  Vertex list  .       +----------+    Face list 
c       |         |                 .       |          |            
c       |         |                 .  4, w |          |             
c       |         |                 .       |          |  2, e      
c       |         |                 .       |          |             
c     1 +---------+ 1               .       +----------+            
c   \1/ 1        2  \2/             .            1, s                 
c
      parameter(le = lx1*ly1*lz1)
c     in 3d ccod(3,6,nelt), vcod(3,8,nelt) 
      real   ccod(2,4,lelt), vcod(2,4,lelt)
      real     fx_c(lx1*lz1,2*ldim,lelt),   fy_c(lx1*lz1,2*ldim,lelt)
     $     , fx_csv(lx1*lz1,2*ldim,lelt), fy_csv(lx1*lz1,2*ldim,lelt) 
      real   x_c(lelt), y_c(lelt)
      logical ngb(4,lelt)

      ne = nx1*ny1*nz1
      nf = nx1*nz1*2*ndim*nelt
      nfaces = 2*ndim 
      ns = (nx1*nz1)/2 ! pointer on faces 

      call ctr2face(fx_c,x_c) ! center point to surface
      call ctr2face(fy_c,y_c)
      call copy(fx_csv,fx_c,nf)
      call copy(fy_csv,fy_c,nf)
      call gs_op (dg_hndl,fx_c,1,1,0)  ! 1 ==> +
      call gs_op (dg_hndl,fy_c,1,1,0)  ! 1 ==> + , a + b 
      call sub2(fx_c,fx_csv,nf) ! a + b - a 
      call sub2(fy_c,fy_csv,nf) ! 

      do ie = 1,nelt
      do iface = 1,nfaces
          if(abs(fx_c(ns,iface,ie)).lt.1e-8 
     $ .and. abs(fy_c(ns,iface,ie)).lt.1e-8) then 
c ! on boundary 
              iv1 = iface       !   1, 2 - 1; 2, 3 - 2
              iv2 = iface + 1   !   3, 4 - 3; 4, 1 - 4
              if(iface.eq.4) iv2 = 1
              ccod(1,iface,ie) = vcod(1,iv1,ie) 
     $              + vcod(1,iv2,ie) - x_c(ie) 
              ccod(2,iface,ie) = vcod(2,iv1,ie) 
     $              + vcod(2,iv2,ie) - y_c(ie) 
              ngb(iface,ie) = .true. ! .true. for on bound
          else                  ! has neighbor, just grab
              ccod(1,iface,ie) = fx_c(ns,iface,ie)
              ccod(2,iface,ie) = fy_c(ns,iface,ie)
          endif
      enddo 
      enddo 

      if(nid.eq.0) write(6,*) 'done :: center coord ' 
      return
      end
c-----------------------------------------------------------------------
      subroutine vert_coor(vcoor) ! ! 2
      include 'SIZE'
      include 'TOTAL'
c
c   Get element vertices coord
c     Input:  xm1,ym1 in TOTAL
c     Output: vcoor(xy,1-4,nelt)
c     
c   \4/ 1        2  \3/             .            3, n     
c     2 +---------+ 2  Vertex list  .       +----------+    Face list 
c       |         |                 .       |          |            
c       |         |                 .  4, w |          |             
c       |         |                 .       |          |  2, e      
c       |         |                 .       |          |             
c     1 +---------+ 1               .       +----------+            
c   \1/ 1        2  \2/             .            1, s                 
c
      parameter(le = lx1*ly1*lz1)
      real      vcoor(2,4,lelt) ! (x:y,1:4,1:nelt)

      ne     = nx1*ny1*nz1
      nfaces = 2*ndim 

      do ie = 1, nelt
          vcoor(1,1,ie) = xm1(  1,  1,1,ie) 
          vcoor(2,1,ie) = ym1(  1,  1,1,ie) 
          vcoor(1,2,ie) = xm1(lx1,  1,1,ie) 
          vcoor(2,2,ie) = ym1(lx1,  1,1,ie) 
          vcoor(1,3,ie) = xm1(lx1,ly1,1,ie) 
          vcoor(2,3,ie) = ym1(lx1,ly1,1,ie) 
          vcoor(1,4,ie) = xm1(  1,ly1,1,ie) 
          vcoor(2,4,ie) = ym1(  1,ly1,1,ie) 
      enddo 

      if(nid.eq.0) write(6,*) 'done :: vertices coord ' 
      return
      end
c-----------------------------------------------------------------------
      subroutine elm_avg(x_c, y_c, rhc, rhuc, rhvc, enc, uxc, uyc
     $                     ,  rh,  rvx,  rvy, enr) ! no rvz 
c     There are already two routines elm_avg3, elm_avg_3 
      include 'SIZE'
      include 'TOTAL'
c   only local, works for MPI 
c
c   Get element averaged values
c     Input : rh(le,nelt), rvx(le,nelt), rvy(le,nelt), enr(le,nelt)
c     Output: x_c(nelt),  y_c(nelt)
c           , rhc(nelt), rhvc(nelt), rhvc(nelt), enc(nelt)
c     
      parameter(le = lx1*ly1*lz1)
      parameter(lf = lx1*lz1*2*ldim*lelt)
      real      x_c(1), y_c(1), rhc(1), rhuc(1), rhvc(1), enc(1)
     $        , uxc(1), uyc(1) 
     $        , rh(le,1), rvx(le,1), rvy(le,1), enr(le,1)
      real      ones(le), tmp(le) 

      nxz    = nx1*nz1
      ne     = nx1*ny1*nz1
      nfaces = 2*ndim 
      nf = nx1*nz1*2*ndim*nelt

      do ie = 1, nelt
          call cfill   (ones,1.,ne)
          call col2    (ones,bm1(1,1,1,ie),ne)
          sum1 = vlsum (ones,ne) ! volume for elem #ie
          call col3    (tmp,bm1(1,1,1,ie),xm1(1,1,1,ie),ne)
          sum2 = vlsum (tmp ,ne) ! vol. sum for x
          call col3    (tmp,bm1(1,1,1,ie),ym1(1,1,1,ie),ne)
          sum3 = vlsum (tmp ,ne) ! vol. sum for y
          call col3    (tmp,bm1(1,1,1,ie),rh(1,ie),ne)
          sum4 = vlsum (tmp ,ne) ! vol. sum for rh
          call col3    (tmp,bm1(1,1,1,ie),rvx(1,ie),ne)
          sum5 = vlsum (tmp ,ne) ! vol. sum for rvx
          call col3    (tmp,bm1(1,1,1,ie),rvy(1,ie),ne)
          sum6 = vlsum (tmp ,ne) ! vol. sum for rvy
          call col3    (tmp,bm1(1,1,1,ie),enr(1,ie),ne)
          sum7 = vlsum (tmp ,ne) ! vol. sum for enr
         
          x_c (ie) = sum2/sum1
          y_c (ie) = sum3/sum1
          rhc (ie) = sum4/sum1
          rhuc(ie) = sum5/sum1
          rhvc(ie) = sum6/sum1
          enc (ie) = sum7/sum1

          uxc (ie) = rhuc(ie)/rhc(ie)
          uyc (ie) = rhvc(ie)/rhc(ie)
      enddo 

c     if(nid.eq.0) write(6,*) 'done :: finding elem avg ' 
      return
      end
c-----------------------------------------------------------------------
c----- 1D limiter, works for elem. aligned on 1 line ( x axis)
c-----------------------------------------------------------------------
      subroutine slopelim1d(rh,rvx,rvy,rvz,enr)
c     Only for 1d special cases, all elements aligned on x axis 
      include 'SIZE'
      include 'TOTAL'
      parameter (lt=lx1*ly1*lz1*lelt)
      real         rh(1),  rvx(1),  rvy(1),  rvz(1),  enr(1)
      real       lrh(lt), lrvx(lt), lrvy(lt), lrvz(lt), lenr(lt)

      n=nx1*ny1*nz1*nelt
      call slopelimN(lrh , rh)
      call slopelimN(lrvx,rvx)
      call slopelimN(lrvy,rvy)
      call slopelimN(lrvz,rvz)
      call slopelimN(lenr,enr)

      call copy(rh , lrh ,n)
      call copy(rvx, lrvx,n)
      call copy(rvy, lrvy,n)
      call copy(rvz, lrvz,n)
      call copy(enr, lenr,n)

      return
      end
c-----------------------------------------------------------------------
      subroutine slopelimN(ulm,u)  
      include 'SIZE'
      include 'TOTAL'
c    Apply slope limiter to u, choose elements that are in need 
c
c       .  u      =   input conservative variable 
c
c       .  ulm    =   limited conservative variable 
c 
c     Simple slope limiter, destorys solution in smooth regions
c 
      parameter(le = lx1*ly1*lz1)
      parameter(lu = lx1*ly1*lz1*lelt)

      real    u (lx1,ly1,lz1,lelt), ulm(lx1,ly1,lz1,lelt)
      real    du(lx1,ly1,lz1,lelt), ul (lx1,ly1,lz1,lelt)
      integer i, ie, iem, iep, n
c scratch 
      real    uave(lelt), xave(lelt), yave(lelt) ,zave(lelt)
      real    us(lelt), ue(lelt), un(lelt)
     $      , uw(lelt), ub(lelt), ut(lelt)

      real    vkm(lelt), vkp(lelt)  ! uave moved left or right by 1 
      real    vel(lelv), ver(lelv)  ! fluxes 

      real    varray(3), difv
      real    er1, er2, eps0
      real    vanf(lx1,lx1), vinv(lx1,lx1)  ! full vandermonde matrix, inv
      real    uhl(lx1) , uhl2(lx1) 
      integer soddir 

      n = le*nelt    
      eps0 = 10.**(-8)       ! epsilon 

      call copy  (ulm, u, n)   ! initialize 

c    1. Find element avg, element end values  
      call cmp_avg_vs(uave, xave, yave, zave, us, ue, un, uw
     $                      , ub,   ut, u)

c    2. Get Vandermonde matrix
      call cmp_van_1d (vanf,vinv)  ! 1D Vandermonde 

c    3. Find interface fluxes 
      do ie = 1, nelv
         iep = ie + 1
         iem = ie - 1
         if (ie .eq. 1) then 
            iem = ie            ! ghost elem left
         else if (ie .eq. nelt) then 
            iep = ie            ! ghost elem right
         endif
         vkm(ie) = uave(iem)    ! get cell avg shifted left and right
         vkp(ie) = uave(iep) 

         varray(1) =  uave(ie) - uw (ie)  ! left flux 
         varray(2) =  uave(ie) - vkm(ie)  ! 
         varray(3) =  vkp(ie) - uave(ie)  ! 
         call minmod ( difv,  varray, 3)  ! get difference in v 
         vel(ie)   =  uave(ie) - difv     ! 

         varray(1) =  ue (ie) - uave(ie)  ! right flux 
         varray(2) =  uave(ie) - vkm(ie)  ! 
         varray(3) =  vkp(ie) - uave(ie)  ! 
         call minmod ( difv,  varray, 3)  ! get difference in v 
         ver(ie)   =  uave(ie) + difv     ! 

c    4. Determine which elements need limiting
c    5. Do limiting to those elements

         er1 = abs(vel(ie) - uw (ie)) 
         er2 = abs(ver(ie) - ue (ie)) 

         if (er1 .gt. eps0 .or. er2 .gt. eps0 ) then

c      5.1 Get linear polynomial -> solve a linear least square problem 
c     ---- vandermonde inverse method ----
c      Evaluate p_{j-1} ( z_i), full Vandermonde matrix, in speclib.f 

            do ju = 1, ly1*lz1 ! for 3d !!! 2d is --> ju = 1,ly1
               call mxm     (vinv,lx1,u(1,ju,1,ie),lx1,uhl,1)   ! vinv * u = uhl
               call rzero   (uhl2,lx1)
               uhl2(1) =  uhl(1) 
               uhl2(2) =  uhl(2)   ! linear 
               call mxm     (vanf,lx1,uhl2,lx1,ul(1,ju,1,ie),1) ! ul = v * uhl
            enddo 

c      5.2  Do limiting 
            v0  = uave(ie) 
            iep = ie + 1
            iem = ie - 1
            if (ie .eq. 1) then 
               iem = ie         ! ghost elem left
            else if (ie .eq. nelv) then 
               iep = ie         ! ghost elem right
            endif
            vm1 = uave(iem) 
            vp1 = uave(iep) 

            call slopelimLin( ulm(1,1,1,ie), ul(1,1,1,ie)
     $              , xm1(1,1,1,ie), xave(ie), vm1, v0, vp1, ie) 

         endif 
      enddo 
      return
      end
c-----------------------------------------------------------------------
      subroutine slopelimLin(ulm,ul,xl,x0,vm1,v0,vp1,ie) 
      include 'SIZE'
      include 'TOTAL'
c 
c     Apply slope limiter to linear function u on coordinate xl 
c vm1 (v minus 1) v0 vp1 (v plus 1) are element averages left, center, right
c
      parameter(le = lx1*ly1*lz1)

      real    u(lx1,ly1,lz1), ulm(lx1,ly1,lz1)
      real    xl(lx1,ly1,lz1), x0  ! x0 ave value 
      real    varray(3), vm1, vp1, v0 
      integer i, ie 

      real du 
      real dudx(lx1,ly1,lz1)

      call copy( ulm, ul, le) 

      xmin = vlmin(xl,le) 
      xmax = vlmax(xl,le) 
      h    = xmax - xmin

      dvm  = (v0 - vm1) / h
      dvp  = (vp1 - v0) / h

      varray(2) = dvm
      varray(3) = dvp

c     Get du dx 
      call mxm  (dxm1, lx1, ul, lx1, dudx, ly1*lz1)  ! dxm1: lx1 x lx1 
      call col2 (dudx, rxm1(1,1,1,ie), le) ! RXM1 = dr/dx * jacm1 
      call col2 (dudx, jacmi(1,ie), le)    ! RXM1 = dr/dx * jacm1 
      do i=1,le

         varray(1) = dudx(i,1,1) 
         call minmod(du,  varray, 3)  ! get new du dx 

         ulm(i,1,1) = v0 + (xl(i,1,1) - x0) * du
      enddo 

      return
      end
c-----------------------------------------------------------------------
      subroutine cmp_avg_vs(uave, xave, yave, zave, us, ue, un, uw, 
     $                        ub,   ut, u)
      include 'SIZE'
      include 'TOTAL'
c  Compute u average, x y avg, ( volume averages) 
c     east, west, south, north, bottom, top face avg  ( surface averages) 
c
      parameter(le = lx1*ly1*lz1)
      real    u(lx1,ly1,lz1,lelt), ons(lx1,ly1,lz1,lelt) 
      real    us(1), ue(1), un(1), uw(1), ub(1), ut(1)
      real    uave(1), xave(1), yave(1), zave(1)
      integer i, ie, iface, nfaces
c     scratch 
      real    sum1, sum2, sum3, sum4 , sum5

      call rone(ons,nx1*ny1*nz1*nelt)
      nfaces = 2*ndim 
      nxz    = nx1*nz1
      do ie = 1, nelv
         sum1 = 0.
         sum2 = 0.
         sum3 = 0.
         sum4 = 0.
         sum5 = 0.
         do i  = 1, le
            sum1 = sum1 + bm1(i,1,1,ie)*u(i,1,1,ie) 
            sum2 = sum2 + bm1(i,1,1,ie)
            sum3 = sum3 + bm1(i,1,1,ie)*xm1(i,1,1,ie) 
            sum4 = sum4 + bm1(i,1,1,ie)*ym1(i,1,1,ie) 
            sum5 = sum5 + bm1(i,1,1,ie)*zm1(i,1,1,ie) 
         enddo 
         uave(ie) = sum1/sum2      ! element avg u 
         xave(ie) = sum3/sum2      ! element avg x 
         yave(ie) = sum4/sum2      ! element avg y 
         zave(ie) = sum5/sum2      ! element avg z 

         do iface = 1, nfaces 
            sum2     = facint_v(ons, area, iface, ie) ! all ones 
            sum1     = facint_v(  u, area, iface, ie)
            if(iface.eq.1) then        ! south 
                 us(ie) = sum1/sum2   
            else if(iface.eq.2) then   ! east  
                 ue(ie) = sum1/sum2   
            else if(iface.eq.3) then   ! north 
                 un(ie) = sum1/sum2   
            else if(iface.eq.4) then   ! west  
                 uw(ie) = sum1/sum2   
            else if(iface.eq.5) then   ! west  
                 ub(ie) = sum1/sum2   
            else if(iface.eq.6) then   ! west  
                 ut(ie) = sum1/sum2   
            endif
         enddo 
      enddo 

      return
      end
c-----------------------------------------------------------------------
      subroutine cmp_van_1d(vanf,vinv) !  1D Vandermonde 
      include 'SIZE'
      include 'TOTAL'
c
      parameter(le = lx1*ly1*lz1)
      real    vanf(lx1,lx1), vinv(lx1,lx1)  ! full vandermonde matrix, inv
      real    z1(lx1), w1(lx1)         ! working variables 
c    Scratch variables for gauss jordan 
      integer indr(lx1), indc(lx1), ipiv(lx1)
      real    rmult(lx1) 
c     real    dv, pm1, pdm1, pm2, pdm2 ! working variables 

      integer iv, jv

      call zwgll (z1,w1,lx1)  ! np = lx1, poly degree = np - 1
      do jv=1,nx1
      do iv=1,nx1  ! index for x-coord 

c         call jacobf(vanf(iv,jv),dv
c    $         ,pm1,pdm1,pm2,pdm2,jv-1,0.,0.,z1(iv)) 

          vanf(iv,jv) = pnleg(z1(iv), jv-1) 

      enddo 
      enddo 

c       Next, find v inv 

      call copy    (vinv,vanf,nx1*nx1)   ! copy into vinv 
      call gaujordf(vinv,nx1,nx1,indr,indc,ipiv,ierr,rmult) 
c      Gauss Jordan w/ full pivoting, full Vandermonde inverse 

      return
      end
c-----------------------------------------------------------------------
      subroutine minmod(sv,v,nv) 
c
c    minmod function: 
c     s =  sum_{i=1,..n} sign(v_i) / n , integer divide
c 
c                       | s *min| v_i | , if |s| = 1
c     m(v1, ... , vn) = 
c                       | 0             , otherwise 

      real    v(1), sv
      real    mnv, abv
      integer nv, i
      integer sum   ! note integer 

      sum = 0 
      mnv = abs(v(1))

      do i=1,nv
         sum = sum + sign(1.,v(i))
         abv = abs(v(i))
         mnv = min(mnv, abv)  ! find the absolute minimum 
      enddo 
      sum = sum / nv ! 1, -1 means the same signs for all v_i  

      sv = 0. 
      if (sum .eq. 1 .or. sum .eq. -1) then
         sv = sum * mnv
      endif 

      return
      end
c-----------------------------------------------------------------------
c-----| End limiters |-----
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c-----| Mass matrix |-----
c     TODO: Full mass matrix oneday maybe? 
c-----------------------------------------------------------------------
      subroutine massf(ub,u)  
      include 'SIZE'
      include 'TOTAL'
      real     ub(1), u(1)
      integer  btyp

      btyp = 3
      if(btyp.eq.1) then      ! btype = 1, rectangular
          call massf_rect(ub,u)
      else if(btyp.eq.2) then ! btype = 2, curved side ! ! 
          call massf_curv(ub,u)
      else if(btyp.eq.3) then ! btype = 3, diagonal 
          call massf_diag(ub,u)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine massf_diag(ub,u)  
c     Diagonal mass matrix 
      include 'SIZE'
      include 'TOTAL'
      real     ub(1), u(1)

      n=nx1*ny1*nz1*nelt
      do i=1,n ! is this right? 
         ub(i) = u(i) * bm1(i,1,1,1)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine massf_rect(ub,u)   ! 
c     Rectangle mass matrix with full mass matrix
      include 'SIZE'
      include 'TOTAL'
      parameter(le=lx1*ly1*lz1)
      parameter(ln=lx1)
      parameter(ln2=ln*ln)
      parameter(lm=lx1+1)
      real     ub(le,1), u(le,1)
      real     b1f(ln*ln), b1i(ln*ln), b1it(ln*ln) ! 1d mass
     $       , b1ft(ln*ln)
      real     jgl(lm,ln), jgt(ln,lm), w(ln*ln)
      real     cons, xl, yl, xmin, xmax

      call full_mass_1d(b1f,b1i,lm,ln,jgl,jgt,w) 

      call transpose (b1ft,ln,b1f,ln) 

      do ie=1,nelt
         xmin = vlmin(xm1(i,1,1,ie),le)
         xmax = vlmax(xm1(i,1,1,ie),le)
         xl   = xmax - xmin
         xmin = vlmin(ym1(i,1,1,ie),le)
         xmax = vlmax(ym1(i,1,1,ie),le)
         yl   = xmax - xmin
         cons = (xl*yl)/4.
         call tensr2(ub(1,ie),ln,u(1,ie),ln,b1f,b1ft,w)
         call cmult (ub(1,ie),cons,le)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine massf_curv(ub,u)  ! ! 
c     Curve sided mass matrix  , not coded up yet 
      include 'SIZE'
      include 'TOTAL'
      real     ub(1), u(1)

      n=nx1*ny1*nz1*nelt
      do i=1,n ! is this right? 
         ub(i) = u(i) * bm1(i,1,1,1)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine invbf(uivb,u)   ! 3 + 2
      include 'SIZE'
      include 'TOTAL'
      real     uivb(1), u(1)
      integer  invbtyp

      invbtyp = 3 ! currently only diagonal works
      if(invbtyp.eq.1) then       ! invbtype = 1, rectangular
          call invbf_rect(uivb,u) ! not working yet 
      else if(invbtyp.eq.2) then  ! invbtype = 2, curved side
          call invbf_curv(uivb,u) ! also not working
      else if(invbtyp.eq.3) then  ! invbtype = 3, diagonal 
          call invbf_diag(uivb,u) ! 
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine invbf_diag(uivb,u)  
c     Inverse diagonal mass matrix 
      include 'SIZE'
      include 'TOTAL'
      real     uivb(1), u(1)

      n=nx1*ny1*nz1*nelt
      do i=1,n ! is this right? 
         uivb(i) = u(i) / bm1(i,1,1,1)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine invbf_rect(uivb,u)  
c     Inverse full mass matrix for rectangular elements
      include 'SIZE'
      include 'TOTAL'
      parameter(le=lx1*ly1*lz1)
      parameter(ln=lx1)
      parameter(ln2=ln*ln)
      parameter(lm=lx1+1)
      real     uivb(le,1), u(le,1)
      real     b1f(ln*ln), b1i(ln*ln), b1it(ln*ln) ! 1d mass
      real     jgl(lm,ln), jgt(ln,lm), w(ln*ln)
      real     cons, xl, yl, xmin, xmax

      call full_mass_1d(b1f,b1i,lm,ln,jgl,jgt,w) 

      call transpose (b1it,ln,b1i,ln)

      do ie=1,nelt
         xmin = vlmin(xm1(i,1,1,ie),le)
         xmax = vlmax(xm1(i,1,1,ie),le)
         xl   = xmax - xmin
         xmin = vlmin(ym1(i,1,1,ie),le)
         xmax = vlmax(ym1(i,1,1,ie),le)
         yl   = xmax - xmin
         cons = 4./(xl*yl)
         call tensr2(uivb(1,ie),ln,u(1,ie),ln,b1i,b1it,w)
         call cmult(uivb(1,ie),cons,le)
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine tensr2(v,nv,u,nu,A,Bt,w)
C
C     Originiated from tensr3 routine in fasts.f
C 
C     -  Tensor product application of v = (B x A) u
C        NOTE -- the transpose of B must be input, rather than B.
C
C     -  scratch arrays: w(nu*nv)
C
C
      include 'SIZE'
      include 'INPUT'
      real v(nv,nv),u(nu,nu)
      real A(1),Bt(1)
      real w(1)

      if (nu.gt.nv) then
         write(6,*) nid,nu,nv,' ERROR in tensr3. Contact P.Fischer.'
         write(6,*) nid,nu,nv,' Memory problem.'
         call exitt
      endif

c     v = A u Bt 
      call mxm(A,nv,u,nu,w,nu)
      call mxm(w,nv,Bt,nu,v,nv)

      return
      end

c-----------------------------------------------------------------------
      subroutine invbf_curv(uivb,u)  
c     Inverse full mass matrix for curved elements 
c        not tested yet 
      include 'SIZE'
      include 'TOTAL'
      parameter(le=lx1*ly1*lz1)
      parameter(ln=lx1)
      parameter(lm=lx1+1)
      real     uivb(le,1), u(le,1)
      ne = nx1*ny1*nz1
      n  = nx1*ny1*nz1*nelt

c     write(6,*) 'not there yet '
      return
      end

c-----------------------------------------------------------------------
      subroutine full_mass_1d(bf,bi,mp,np,jgl,jgt,w)
c
c     Generate full mass matrix
c
c      jgl  = interpolation matrix, mapping from velocity nodes to pressure
c      jgt  = transpose of interpolation matrix
c      w    = work array of size (np+mp)
c
c      np   = number of points on GLL grid
c      mp   = number of points on GL  grid
c
c
      real bf(np*np),bi(np*np)

      real    jgl(mp,np),jgt(np*mp),w(1)

      integer idr(np), idc(np), ipiv(np)

      iz = 1
      id = iz + np

      call zwgll (w(iz),bi,np)
      call zwgl  (w(id),bi,mp) ! bi contains mass matrix on Gauss points

      n  = np-1
      do i=1,mp
         call fd_weights_full(w(id+i-1),w(iz),n,0,jgt)
         do j=1,np
            jgl(i,j) = jgt(j)                  !  Interpolation matrix
         enddo
      enddo

      call transpose(jgt,np,jgl,mp)

      do j=1,np
      do i=1,mp
         jgl(i,j) =  bi(i)*jgl(i,j)   ! W*Jt
      enddo
      enddo

      call mxm      (jgt,np,jgl,mp,bf,np)  ! Bf = J*W*Jt
      call copy     (bi,bf,np*np)

c     call gaujordf (bi,np,np,jgt,jgl,w,ierr,jgt(1+np))

      call gaujordf (bi,np,np,idr,idc,ipiv,ierr,w)

c    both call work
c    working arrays don't really matter, just to be clear 

      return
      end
c-----------------------------------------------------------------------
c=======================================================================
c----- Source: ~fischer/nek5_svn/trunk/nek_fullb/hmholtz.f
c----- Deleted for now, if need it go back and grab it 
c=======================================================================
c-----------------------------------------------------------------------
c----- Volume de-aliasing 
c-----------------------------------------------------------------------
      subroutine grad_m1_t2(u,fx,fy,ifdal)
c     New task: get de-aliasing for volume term, July 29th
c                 T
c     Compute grad  of (fx,fy) and set to u.
c
c   - Input : fx, fy, e
c   - Output: u
c   DxT (M fx)  + (M fy) Dy  
c
      real      u(1),fx(1),fy(1)
      logical   ifdal
c
      if (ifdal) then    ! de-aliased
          call grad_dbm1_t2(u,fx,fy) ! de-aliased
      else               ! not de-aliased
          call grad_bm1_t2 (u,fx,fy) ! not de-aliased
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine grad_dbm1_t2(u,fx,fy) ! ! 
c     De-aliased form not correct, y -term is problematic 
c     Fri Jul 31 11:27:12 CDT 2015
c     
c     Get de-aliasing for volume term, July 29th
c                 T
c     Compute grad  of (fx,fy,fz) and set to u.
c
c     Input : fx, fy, fz  |   Output: u
      include 'SIZE'
      include 'GEOM'
      include 'INPUT' ! need if3d from here 
c
      parameter (lxyz=lx1*ly1*lz1)
      parameter (ltd =lxd*lyd*lzd)
      parameter (ldg =lxd**3,lwkd=4*lxd*lxd) ! ln 541, from convect.f
      common /ctmp1/ v(lxyz), w2(ltd)
      common /ctmp2/ jfx(ltd), jfy(ltd), jfz(ltd)
     $             , tm1(ltd), tm2(ltd), tm3(ltd)
     $             , jv (ltd), w(lxd), z(lxd)
     $             , dxmd(lxd,lyd), dxtmd(lxd,lyd)
      real     jfx, jfy, jfz, jv
c     dt really bothers me here, causes me not being able to 
      common /dgrad/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real     jgl, jgt
c
      real     u(lxyz,1),fx(lxyz,1),fy(lxyz,1),fz(lxyz,1)
      real     fr(ltd), fs(ltd), ft(ltd)
      integer  e, nxyz1, nxyzd, mx, md 
c
      nxyz1 = nx1*ny1*nz1
      nxyzd = nxd*nyd*nzd
      nxyd  = nxd*nyd
      mx = nx1 
      md = lxd 
c
c     in induct.f 
      call set_dealias_rx ! get rx ( in GEOM) , zwgl is used
c
c     alpha = 0, beta = 0 could be either GL or GLL ??
c     Fri Jul 31 17:02:46 CDT 2015
      call zwgl ( z, w, md)  ! ! GL points , alpha = 0, beta = 0 ?
c
c  Works, but a bit messy because this routine uses common block
c     dgrad ( dg and dgt) 
      call get_dgl_ptr(ip,md,md) 
      call copy ( dxmd, dg(ip),nxyd)
      call copy (dxtmd,dgt(ip),nxyd) ! !

      if (if3d) then
c   
c     -  in navier5.f, convect_cons is for: v div (C u)
c        I need (grad v) . F though, namely J^T D^T Bf JF w/ jaco
c
         do e=1,nelt
c        
         call rzero(u(1,e),lxyz) ! set to 0 
         call intp_rstd(jfx,fx(1,e),mx,md,if3d,0) ! forward, GLL to GL
         call intp_rstd(jfy,fy(1,e),mx,md,if3d,0) ! forward 
         call intp_rstd(jfz,fz(1,e),mx,md,if3d,0) ! forward 
c   
         do i=1,nxyzd
             fr(i) = rx(i,1,e)*jfx(i) ! rxm1
     $             + rx(i,2,e)*jfy(i) ! rym1, has mass embeded in
     $             + rx(i,3,e)*jfz(i) ! rzm1, rx matrix 
             fs(i) = rx(i,4,e)*jfx(i) ! sxm1, refer to induct.f
     $             + rx(i,5,e)*jfy(i) ! sym1
     $             + rx(i,6,e)*jfz(i) ! szm1
             ft(i) = rx(i,7,e)*jfx(i) ! txm1
     $             + rx(i,8,e)*jfy(i) ! tym1
     $             + rx(i,9,e)*jfz(i) ! tzm1
         enddo 
c
c   Input: fr, fs, ft   Output: v
c        dxtm1 fr + (do loop dym1) fs + ft dzm1
         call local_grad3_t(jv,fr,fs,ft,lxd-1,1,dxmd,dxtmd,w2) !!??
c
         call intp_rstd(v, jv, mx, md, if3d, 1) ! backward ( transpose) 
c
         call add2 ( u(1,e),v,lxyz)
c
         enddo 
      else 
         call exitt
         write(6,*) 'exit :: 2d in grad_dbm1_t'
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine grad_bm1_t2(u,fx,fy) ! 
c                 T
c     Compute grad  of (fx,fy,fz) and set to u.
c
c     Single element case, not de-aliased 
c
c     Input : fx, fy, e
c     Output: u
c   DxT (M fx)  + (M fy) Dy  
      include 'SIZE'
      include 'TOTAL'
c
      parameter (lxyz=lx1*ly1*lz1)
      common /ctmp1/ v(lxyz),w(lxyz),w2(lxyz),w3(lxyz)
      real u(lxyz,1),fx(lxyz,1),fy(lxyz,1),fz(lxyz,1)
      real fr(lxyz), fs(lxyz), ft(lxyz)
      integer e
c
      if (if3d) then

         do e = 1,nelt

         call rzero(u(1,e),lxyz) ! set to 0 
         call copy (w, fx(1,e), lxyz) 
         call col2 (w, w3m1,lxyz) ! 
         call copy (w2, fy(1,e), lxyz) 
         call col2 (w2, w3m1,lxyz) ! 
         call copy (w3, fz(1,e), lxyz) 
         call col2 (w3, w3m1,lxyz) ! 

c   what should fr, fs, ft be for fx, fy, fz ? 
         do i=1,lx1*ly1*lz1
             fr(i) = rxm1(i,1,1,e)*w (i) 
     $             + rym1(i,1,1,e)*w2(i)
     $             + rzm1(i,1,1,e)*w3(i) ! jacm1 
             fs(i) = sxm1(i,1,1,e)*w (i) 
     $             + sym1(i,1,1,e)*w2(i)
     $             + szm1(i,1,1,e)*w3(i)
             ft(i) = txm1(i,1,1,e)*w (i) 
     $             + tym1(i,1,1,e)*w2(i)
     $             + tzm1(i,1,1,e)*w3(i)
         enddo 
c   Input: fr, fs, ft   Output: w
c        dxtm1 fr + (do loop dym1) fs + ft dzm1
         call local_grad3_t(v,fr,fs,ft,lx1-1,1,dxm1,dxtm1,w2)
c
         call add2 ( u(1,e),v,lxyz)
c
         enddo
      else 
         call exitt
         write(6,*) 'exit :: 2d in grad_bm1_t'
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine grad_m1_t(u,fx,fy,fz,ifdal) ! 2+3???
c     New task: get de-aliasing for volume term, July 29th
c                 T
c     Compute grad  of (fx,fy,fz) and set to u.
c
c   - Input : fx, fy, fz, e
c   - Output: u
c   DxT (M fx)  + (M fy) Dy  + DzT (M fz) 
c
      real      u(1),fx(1),fy(1),fz(1)
      logical   ifdal
c
      if (ifdal) then    ! de-aliased
          call grad_dbm1_t(u,fx,fy,fz) ! de-aliased
      else               ! not de-aliased
          call grad_bm1_t (u,fx,fy,fz) ! not de-aliased
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine grad_dbm1_t_col2add3(u
     $                             , fx1,fx2,fy1,fy2,fz1,fz2,fpr) ! ! 
c     
c     Multiple input arrays so fx1*fx2, fy1*fy2, fz1*fz2 is not aliased
c                 T
c     Compute grad  of (fx,fy,fz) and set to u.
c
c     Input : fx1, fy1, fz1  |   Output: u
c           , fx2, fy2, fz2  |   
c                     , fpr  |   
c     fx = fx1*fx2 , fy = fy1*fy2, fz = fz1*fz2 + fpr  |   

      include 'SIZE'
      include 'GEOM'
      include 'INPUT' ! need if3d from here 
c
      parameter (lxyz=lx1*ly1*lz1)
      parameter (ltd =lxd*lyd*lzd)
      parameter (ldg =lxd**3,lwkd=4*lxd*lxd) ! ln 541, from convect.f
      common /ctmp1/ v(lxyz), w2(ltd)
      common /ctmp2/ jfx1(ltd), jfy1(ltd), jfz1(ltd)
     $             , jfx2(ltd), jfy2(ltd), jfz2(ltd)
     $             , jfpr(ltd)
     $             , tm1(ltd), tm2(ltd), tm3(ltd)
     $             , jv (ltd), w(lxd), z(lxd)
     $             , dxmd(lxd,lyd), dxtmd(lxd,lyd)
      real     jfx1, jfy1, jfz1
     $       , jfx2, jfy2, jfz2, jfpr, jv
c     
      common /dgrad/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real     jgl, jgt
c
      real     u  (lxyz,1)
     $       , fx1(lxyz,1),fy1(lxyz,1),fz1(lxyz,1)
     $       , fx2(lxyz,1),fy2(lxyz,1),fz2(lxyz,1)
      real     fpr(ltd,1) ! already on fine mesh 
      real     fr(ltd), fs(ltd), ft(ltd)
      integer  e, nxyz1, nxyzd, mx, md 
c
      nxyz1 = nx1*ny1*nz1
      nxyzd = nxd*nyd*nzd
      nxyd  = nxd*nyd ! for derivative matricies
      mx = nx1 
      md = nxd 
c
c     in induct.f 
      call set_dealias_rx ! get rx ( in GEOM) , zwgl is used
c
      call zwgl ( z, w, md)  ! ! GL points , alpha = 0, beta = 0 ?
c     jgl, jgt are used in intp_rstd routine 
      call get_dgl_ptr(ip,md,md) 
      call copy ( dxmd, dg(ip),nxyd)
      call copy (dxtmd,dgt(ip),nxyd) ! !
c
      if (if3d) then
c
         do e=1,nelt
c        
         call rzero(u(1,e),lxyz) ! set to 0 
         call intp_rstd(jfx1,fx1(1,e),mx,md,if3d,0) ! forward, GLL to GL
         call intp_rstd(jfx2,fx2(1,e),mx,md,if3d,0) ! forward, GLL to GL
         call intp_rstd(jfy1,fy1(1,e),mx,md,if3d,0) ! forward 
         call intp_rstd(jfy2,fy2(1,e),mx,md,if3d,0) ! forward 
         call intp_rstd(jfz1,fz1(1,e),mx,md,if3d,0) ! forward 
         call intp_rstd(jfz2,fz2(1,e),mx,md,if3d,0) ! forward 
c
         call col2     (jfx1, jfx2, ltd)
         call col2     (jfy1, jfy2, ltd) ! 2 parts mult. together
         call col2     (jfz1, jfz2, ltd)
         call add2     (jfz1,  fpr(1,e), ltd) ! fz1*fz2 + p
c   
         do i=1,nxyzd
             fr(i) = rx(i,1,e)*jfx1(i) ! rxm1
     $             + rx(i,2,e)*jfy1(i) ! rym1, has mass embeded in
     $             + rx(i,3,e)*jfz1(i) ! rzm1, rx matrix 
             fs(i) = rx(i,4,e)*jfx1(i) ! sxm1, refer to induct.f
     $             + rx(i,5,e)*jfy1(i) ! sym1
     $             + rx(i,6,e)*jfz1(i) ! szm1
             ft(i) = rx(i,7,e)*jfx1(i) ! txm1
     $             + rx(i,8,e)*jfy1(i) ! tym1
     $             + rx(i,9,e)*jfz1(i) ! tzm1
         enddo 
c
         call local_grad3_t(jv,fr,fs,ft,nxd-1,1,dxmd,dxtmd,w2) !!??
c
         call intp_rstd(v, jv, mx, md, if3d, 1) ! backward ( transpose) 
c
         call add2 ( u(1,e),v,lxyz)
c
         enddo 
      else  
c     this routine shouldn't have 2d because it deals with rho w eq. 
         call exitt
         write(6,*) 'exit :: 2d in grad_dbm1_t_col2add3'
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine grad_dbm1_t_col2add2(u
     $                             , fx1,fx2,fy1,fy2,fpr,fz1,fz2) ! ! 
c     
c     Multiple input arrays so fx1*fx2, fy1*fy2, fz1*fz2 is not aliased
c                 T
c     Compute grad  of (fx,fy,fz) and set to u.
c
c     Input : fx1, fy1, fz1  |   Output: u
c           , fx2, fy2, fz2  |   
c                , fpr       |   
c     fx = fx1*fx2 , fy = fy1*fy2 + fpr, fz = fz1*fz2  |   

      include 'SIZE'
      include 'GEOM'
      include 'INPUT' ! need if3d from here 
c
      parameter (lxyz=lx1*ly1*lz1)
      parameter (ltd =lxd*lyd*lzd)
      parameter (ldg =lxd**3,lwkd=4*lxd*lxd) ! ln 541, from convect.f
      common /ctmp1/ v(lxyz), w2(ltd)
      common /ctmp2/ jfx1(ltd), jfy1(ltd), jfz1(ltd)
     $             , jfx2(ltd), jfy2(ltd), jfz2(ltd)
     $             , jfpr(ltd)
     $             , tm1(ltd), tm2(ltd), tm3(ltd)
     $             , jv (ltd), w(lxd), z(lxd)
     $             , dxmd(lxd,lyd), dxtmd(lxd,lyd)
      real     jfx1, jfy1, jfz1
     $       , jfx2, jfy2, jfz2, jfpr, jv
c     
      common /dgrad/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real     jgl, jgt
c
      real     u  (lxyz,1)
     $       , fx1(lxyz,1),fy1(lxyz,1),fz1(lxyz,1)
     $       , fx2(lxyz,1),fy2(lxyz,1),fz2(lxyz,1)
      real     fpr(ltd,1) ! already on fine mesh 
      real     fr(ltd), fs(ltd), ft(ltd)
      integer  e, nxyz1, nxyzd, mx, md 
c
      nxyz1 = nx1*ny1*nz1
      nxyzd = nxd*nyd*nzd
      nxyd  = nxd*nyd ! for derivative matricies
      mx = nx1 
      md = nxd 
c
c     in induct.f 
      call set_dealias_rx ! get rx ( in GEOM) , zwgl is used
c
      call zwgl ( z, w, md)  ! ! GL points , alpha = 0, beta = 0 ?
c     jgl, jgt are used in intp_rstd routine 
      call get_dgl_ptr(ip,md,md) 
      call copy ( dxmd, dg(ip),nxyd)
      call copy (dxtmd,dgt(ip),nxyd) ! !
c
      if (if3d) then
c
         do e=1,nelt
c        
         call rzero(u(1,e),lxyz) ! set to 0 
         call intp_rstd(jfx1,fx1(1,e),mx,md,if3d,0) ! forward, GLL to GL
         call intp_rstd(jfx2,fx2(1,e),mx,md,if3d,0) ! forward, GLL to GL
         call intp_rstd(jfy1,fy1(1,e),mx,md,if3d,0) ! forward 
         call intp_rstd(jfy2,fy2(1,e),mx,md,if3d,0) ! forward 
         call intp_rstd(jfz1,fz1(1,e),mx,md,if3d,0) ! forward 
         call intp_rstd(jfz2,fz2(1,e),mx,md,if3d,0) ! forward 
c
         call col2     (jfx1, jfx2, ltd)
         call col2     (jfy1, jfy2, ltd) ! 2 parts mult. together
         call add2     (jfy1,  fpr(1,e), ltd) ! fy1*fy2 + p
         call col2     (jfz1, jfz2, ltd)
c   
         do i=1,nxyzd
             fr(i) = rx(i,1,e)*jfx1(i) ! rxm1
     $             + rx(i,2,e)*jfy1(i) ! rym1, has mass embeded in
     $             + rx(i,3,e)*jfz1(i) ! rzm1, rx matrix 
             fs(i) = rx(i,4,e)*jfx1(i) ! sxm1, refer to induct.f
     $             + rx(i,5,e)*jfy1(i) ! sym1
     $             + rx(i,6,e)*jfz1(i) ! szm1
             ft(i) = rx(i,7,e)*jfx1(i) ! txm1
     $             + rx(i,8,e)*jfy1(i) ! tym1
     $             + rx(i,9,e)*jfz1(i) ! tzm1
         enddo 
c
c   Input: fr, fs, ft   Output: v
c        dxtm1 fr + (do loop dym1) fs + ft dzm1
         call local_grad3_t(jv,fr,fs,ft,nxd-1,1,dxmd,dxtmd,w2) !!??
c
         call intp_rstd(v, jv, mx, md, if3d, 1) ! backward ( transpose) 
c
         call add2 ( u(1,e),v,lxyz)
c
         enddo 
      else  ! 2D 
         do e=1,nelt
         call rzero(u(1,e),lxyz) ! set to 0 
         call intp_rstd(jfx1,fx1(1,e),mx,md,if3d,0) ! forward, GLL to GL
         call intp_rstd(jfx2,fx2(1,e),mx,md,if3d,0) ! forward, GLL to GL
         call intp_rstd(jfy1,fy1(1,e),mx,md,if3d,0) ! forward 
         call intp_rstd(jfy2,fy2(1,e),mx,md,if3d,0) ! forward 
c 
         call col2     (jfx1, jfx2, ltd)
         call col2     (jfy1, jfy2, ltd) ! 2 parts mult. together
         call add2     (jfy1,  fpr(1,e), ltd) ! fx1*fx2 + p
c 
         do i=1,nxyzd
             fr(i) = rx(i,1,e)*jfx1(i) ! rxm1, has mass embeded in
     $             + rx(i,2,e)*jfy1(i) ! rym1, rx matrix 
             fs(i) = rx(i,3,e)*jfx1(i) ! sxm1, refer to induct.f
     $             + rx(i,4,e)*jfy1(i) ! sym1
         enddo 
c
         call local_grad2_t(jv,fr,fs,lxd-1,1,dxmd,dxtmd,w2) !!??
c
         call intp_rstd(v, jv, mx, md, if3d, 1) ! backward ( transpose) 
c
         call add2 ( u(1,e),v,lxyz)

         enddo

      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine grad_dbm1_t_col2add1(u
     $                             , fx1,fx2,fpr,fy1,fy2,fz1,fz2) ! ! 
c     
c     Multiple input arrays so fx1*fx2, fy1*fy2, fz1*fz2 is not aliased
c                 T
c     Compute grad  of (fx,fy,fz) and set to u.
c
c     Input : fx1, fy1, fz1  |   Output: u
c           , fx2, fy2, fz2  |   
c           , fpr            |   
c     fx = fx1*fx2 + fpr, fy = fy1*fy2, fz = fz1*fz2  |   

      include 'SIZE'
      include 'GEOM'
      include 'INPUT' ! need if3d from here 
c
      parameter (lxyz=lx1*ly1*lz1)
      parameter (ltd =lxd*lyd*lzd)
      parameter (ldg =lxd**3,lwkd=4*lxd*lxd) ! ln 541, from convect.f
      common /ctmp1/ v(lxyz), w2(ltd)
      common /ctmp2/ jfx1(ltd), jfy1(ltd), jfz1(ltd)
     $             , jfx2(ltd), jfy2(ltd), jfz2(ltd)
     $             , jfpr(ltd)
     $             , tm1(ltd), tm2(ltd), tm3(ltd)
     $             , jv (ltd), w(lxd), z(lxd)
     $             , dxmd(lxd,lyd), dxtmd(lxd,lyd)
      real     jfx1, jfy1, jfz1
     $       , jfx2, jfy2, jfz2, jfpr, jv
c     
      common /dgrad/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real     jgl, jgt
c
      real     u  (lxyz,1)
     $       , fx1(lxyz,1),fy1(lxyz,1),fz1(lxyz,1)
     $       , fx2(lxyz,1),fy2(lxyz,1),fz2(lxyz,1)
      real     fpr(ltd,1) ! already on fine mesh 
      real     fr(ltd), fs(ltd), ft(ltd)
      integer  e, nxyz1, nxyzd, mx, md 
c
      nxyz1 = nx1*ny1*nz1
      nxyzd = nxd*nyd*nzd
      nxyd  = nxd*nyd ! for derivative matricies
      mx = nx1 
      md = nxd 
c
c     in induct.f 
      call set_dealias_rx ! get rx ( in GEOM) , zwgl is used
c
      call zwgl ( z, w, md)  ! ! GL points , alpha = 0, beta = 0 ?
c     jgl, jgt are used in intp_rstd routine 
      call get_dgl_ptr(ip,md,md) 
      call copy ( dxmd, dg(ip),nxyd)
      call copy (dxtmd,dgt(ip),nxyd) ! !
c
      if (if3d) then
c
         do e=1,nelt
c        
         call rzero(u(1,e),lxyz) ! set to 0 
         call intp_rstd(jfx1,fx1(1,e),mx,md,if3d,0) ! forward, GLL to GL
         call intp_rstd(jfx2,fx2(1,e),mx,md,if3d,0) ! forward, GLL to GL
         call intp_rstd(jfy1,fy1(1,e),mx,md,if3d,0) ! forward 
         call intp_rstd(jfy2,fy2(1,e),mx,md,if3d,0) ! forward 
         call intp_rstd(jfz1,fz1(1,e),mx,md,if3d,0) ! forward 
         call intp_rstd(jfz2,fz2(1,e),mx,md,if3d,0) ! forward 
c
         call col2     (jfx1, jfx2, ltd)
         call add2     (jfx1, fpr(1,e), ltd) ! fx1*fx2 + p
         call col2     (jfy1, jfy2, ltd) ! 2 parts mult. together
         call col2     (jfz1, jfz2, ltd)
c   
         do i=1,nxyzd
             fr(i) = rx(i,1,e)*jfx1(i) ! rxm1
     $             + rx(i,2,e)*jfy1(i) ! rym1, has mass embeded in
     $             + rx(i,3,e)*jfz1(i) ! rzm1, rx matrix 
             fs(i) = rx(i,4,e)*jfx1(i) ! sxm1, refer to induct.f
     $             + rx(i,5,e)*jfy1(i) ! sym1
     $             + rx(i,6,e)*jfz1(i) ! szm1
             ft(i) = rx(i,7,e)*jfx1(i) ! txm1
     $             + rx(i,8,e)*jfy1(i) ! tym1
     $             + rx(i,9,e)*jfz1(i) ! tzm1
         enddo 
c
c   Input: fr, fs, ft   Output: v
c        dxtm1 fr + (do loop dym1) fs + ft dzm1
         call local_grad3_t(jv,fr,fs,ft,nxd-1,1,dxmd,dxtmd,w2) !!??
c
         call intp_rstd(v, jv, mx, md, if3d, 1) ! backward ( transpose) 
c
         call add2 ( u(1,e),v,lxyz)
c
         enddo 
      else 
         do e=1,nelt
         call rzero(u(1,e),lxyz) ! set to 0 
         call intp_rstd(jfx1,fx1(1,e),mx,md,if3d,0) ! forward, GLL to GL
         call intp_rstd(jfx2,fx2(1,e),mx,md,if3d,0) ! forward, GLL to GL
         call intp_rstd(jfy1,fy1(1,e),mx,md,if3d,0) ! forward 
         call intp_rstd(jfy2,fy2(1,e),mx,md,if3d,0) ! forward 
c 
         call col2     (jfx1, jfx2, ltd)
         call add2     (jfx1,  fpr(1,e), ltd) ! fx1*fx2 + p
         call col2     (jfy1, jfy2, ltd) ! 2 parts mult. together
c 
         do i=1,nxyzd
             fr(i) = rx(i,1,e)*jfx1(i) ! rxm1, has mass embeded in
     $             + rx(i,2,e)*jfy1(i) ! rym1, rx matrix 
             fs(i) = rx(i,3,e)*jfx1(i) ! sxm1, refer to induct.f
     $             + rx(i,4,e)*jfy1(i) ! sym1
         enddo 
c
         call local_grad2_t(jv,fr,fs,lxd-1,1,dxmd,dxtmd,w2) !!??
c
         call intp_rstd(v, jv, mx, md, if3d, 1) ! backward ( transpose) 
c
         call add2 ( u(1,e),v,lxyz)

         enddo
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine grad_dbm1_t_col2(u,fx1,fx2,fy1,fy2,fz1,fz2) ! ! 
c     Won't compute x^3 exactly. Don't know why yet.
c     Fri Aug  7 12:32:06 CDT 2015
c     
c     2 input arrays so fx1*fx2, fy1*fy2, fz1*fz2 is not aliased
c                 T
c     Compute grad  of (fx,fy,fz) and set to u.
c
c     Input : fx1, fy1, fz1  |   Output: u
c           , fx2, fy2, fz2  |   
c     fx = fx1*fx2, fy = fy1*fy2, fz = fz1*fz2  |   
      include 'SIZE'
      include 'GEOM'
      include 'INPUT' ! need if3d from here 
c
      parameter (lxyz=lx1*ly1*lz1)
      parameter (ltd =lxd*lyd*lzd)
      parameter (ldg =lxd**3,lwkd=4*lxd*lxd) ! ln 541, from convect.f
      common /ctmp1/ v(lxyz), w2(ltd)
      common /ctmp2/ jfx1(ltd), jfy1(ltd), jfz1(ltd)
     $             , jfx2(ltd), jfy2(ltd), jfz2(ltd)
     $             , jfpr(ltd)
     $             , tm1(ltd), tm2(ltd), tm3(ltd)
     $             , jv (ltd), w(lxd), z(lxd)
     $             , dxmd(lxd,lyd), dxtmd(lxd,lyd)
      real     jfx1, jfy1, jfz1
     $       , jfx2, jfy2, jfz2, jfpr, jv
c     
      common /dgrad/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real     jgl, jgt
c
      real     u  (lxyz,1)
     $       , fx1(lxyz,1),fy1(lxyz,1),fz1(lxyz,1)
     $       , fx2(lxyz,1),fy2(lxyz,1),fz2(lxyz,1)
      real     fr(ltd), fs(ltd), ft(ltd)
      integer  e, nxyz1, nxyzd, mx, md 
c
      nxyz1 = nx1*ny1*nz1
      nxyzd = nxd*nyd*nzd
      nxyd  = nxd*nyd ! for derivative matricies
      mx = nx1 
      md = nxd 
c
c     in induct.f 
      call set_dealias_rx ! get rx ( in GEOM) , zwgl is used
c
      call zwgl ( z, w, md)  ! ! GL points , alpha = 0, beta = 0 ?
c     jgl, jgt are used in intp_rstd routine 
      call get_dgl_ptr(ip,md,md) 
      call copy ( dxmd, dg(ip),nxyd)
      call copy (dxtmd,dgt(ip),nxyd) ! !
c
      if (if3d) then
c
         do e=1,nelt
c        
         call rzero(u(1,e),lxyz) ! set to 0 
         call intp_rstd(jfx1,fx1(1,e),mx,md,if3d,0) ! forward, GLL to GL
         call intp_rstd(jfx2,fx2(1,e),mx,md,if3d,0) ! forward, GLL to GL
         call intp_rstd(jfy1,fy1(1,e),mx,md,if3d,0) ! forward 
         call intp_rstd(jfy2,fy2(1,e),mx,md,if3d,0) ! forward 
         call intp_rstd(jfz1,fz1(1,e),mx,md,if3d,0) ! forward 
         call intp_rstd(jfz2,fz2(1,e),mx,md,if3d,0) ! forward 
c
         call col2     (jfx1, jfx2, ltd)
         call col2     (jfy1, jfy2, ltd) ! 2 parts mult. together
         call col2     (jfz1, jfz2, ltd)
c   
         do i=1,nxyzd
             fr(i) = rx(i,1,e)*jfx1(i) ! rxm1
     $             + rx(i,2,e)*jfy1(i) ! rym1, has mass embeded in
     $             + rx(i,3,e)*jfz1(i) ! rzm1, rx matrix 
             fs(i) = rx(i,4,e)*jfx1(i) ! sxm1, refer to induct.f
     $             + rx(i,5,e)*jfy1(i) ! sym1
     $             + rx(i,6,e)*jfz1(i) ! szm1
             ft(i) = rx(i,7,e)*jfx1(i) ! txm1
     $             + rx(i,8,e)*jfy1(i) ! tym1
     $             + rx(i,9,e)*jfz1(i) ! tzm1
         enddo 
c
c   Input: fr, fs, ft   Output: v
c        dxtm1 fr + (do loop dym1) fs + ft dzm1
         call local_grad3_t(jv,fr,fs,ft,nxd-1,1,dxmd,dxtmd,w2) !!??
c
         call intp_rstd(v, jv, mx, md, if3d, 1) ! backward ( transpose) 
c
         call add2 ( u(1,e),v,lxyz)
c
         enddo 
      else 
         do e=1,nelt
         call rzero(u(1,e),lxyz) ! set to 0 
         call intp_rstd(jfx1,fx1(1,e),mx,md,if3d,0) ! forward, GLL to GL
         call intp_rstd(jfx2,fx2(1,e),mx,md,if3d,0) ! forward, GLL to GL
         call intp_rstd(jfy1,fy1(1,e),mx,md,if3d,0) ! forward 
         call intp_rstd(jfy2,fy2(1,e),mx,md,if3d,0) ! forward 
c 
         call col2     (jfx1, jfx2, ltd)
         call col2     (jfy1, jfy2, ltd) ! 2 parts mult. together
c 
         do i=1,nxyzd
             fr(i) = rx(i,1,e)*jfx1(i) ! rxm1, has mass embeded in
     $             + rx(i,2,e)*jfy1(i) ! rym1, rx matrix 
             fs(i) = rx(i,3,e)*jfx1(i) ! sxm1, refer to induct.f
     $             + rx(i,4,e)*jfy1(i) ! sym1
         enddo 
c
         call local_grad2_t(jv,fr,fs,lxd-1,1,dxmd,dxtmd,w2) !!??
c
         call intp_rstd(v, jv, mx, md, if3d, 1) ! backward ( transpose) 
c
         call add2 ( u(1,e),v,lxyz)

         enddo
c        call exitt
c        call grad_dbm1_t2_col2(tm2,f1,g1,h1,ifdeal)
c        write(6,*) 'exit :: 2d in grad_dbm1_t_2'
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine grad_dbm1_t(u,fx,fy,fz) ! ! 
c     De-aliased form not correct, y -term is problematic 
c     Fri Jul 31 11:27:12 CDT 2015
c     
c     Get de-aliasing for volume term, July 29th
c                 T
c     Compute grad  of (fx,fy,fz) and set to u.
c
c     Input : fx, fy, fz  |   Output: u
      include 'SIZE'
      include 'GEOM'
      include 'INPUT' ! need if3d from here 
c
      parameter (lxyz=lx1*ly1*lz1)
      parameter (ltd =lxd*lyd*lzd)
      parameter (ldg =lxd**3,lwkd=4*lxd*lxd) ! ln 541, from convect.f
      common /ctmp1/ v(lxyz), w2(ltd)
      common /ctmp2/ jfx (ltd), jfy (ltd), jfz (ltd)
     $             , jfpr(ltd)
     $             , tm1(ltd), tm2(ltd), tm3(ltd)
     $             , jv (ltd), w(lxd), z(lxd)
     $             , dxmd(lxd,lyd), dxtmd(lxd,lyd)
      real     jfx , jfy , jfz 
     $       , jfpr, jv
c     dt really bothers me here, causes me not being able to 
      common /dgrad/ d(ldg),dt(ldg),dg(ldg),dgt(ldg),jgl(ldg),jgt(ldg)
     $             , wkd(lwkd)
      real     jgl, jgt
c
      real     u(lxyz,1),fx(lxyz,1),fy(lxyz,1),fz(lxyz,1)
      real     fr(ltd), fs(ltd), ft(ltd)
      integer  e, nxyz1, nxyzd, mx, md 
      real ons(lxyz)
c
      nxyz1 = nx1*ny1*nz1
      nxyzd = nxd*nyd*nzd
      nxyd  = nxd*nyd
      mx = nx1 
      md = lxd 
c
c     in induct.f 
      call set_dealias_rx ! get rx ( in GEOM) , zwgl is used
c
c     alpha = 0, beta = 0 could be either GL or GLL ??
c     Fri Jul 31 17:02:46 CDT 2015
      call zwgl ( z, w, md)  ! ! GL points , alpha = 0, beta = 0 ?
c
c  Works, but a bit messy because this routine uses common block
c     dgrad ( dg and dgt) 
      call get_dgl_ptr(ip,md,md) 
      call copy ( dxmd, dg(ip),nxyd)
      call copy (dxtmd,dgt(ip),nxyd) ! !

      if (if3d) then
c   
c     -  in navier5.f, convect_cons is for: v div (C u)
c        I need (grad v) . F though, namely J^T D^T Bf JF w/ jaco
c
         do e=1,nelt
c        
         call rzero(u(1,e),lxyz) ! set to 0 
         call intp_rstd(jfx,fx(1,e),mx,md,if3d,0) ! forward, GLL to GL
         call intp_rstd(jfy,fy(1,e),mx,md,if3d,0) ! forward 
         call intp_rstd(jfz,fz(1,e),mx,md,if3d,0) ! forward 
c   
         do i=1,nxyzd
             fr(i) = rx(i,1,e)*jfx(i) ! rxm1
     $             + rx(i,2,e)*jfy(i) ! rym1, has mass embeded in
     $             + rx(i,3,e)*jfz(i) ! rzm1, rx matrix 
             fs(i) = rx(i,4,e)*jfx(i) ! sxm1, refer to induct.f
     $             + rx(i,5,e)*jfy(i) ! sym1
     $             + rx(i,6,e)*jfz(i) ! szm1
             ft(i) = rx(i,7,e)*jfx(i) ! txm1
     $             + rx(i,8,e)*jfy(i) ! tym1
     $             + rx(i,9,e)*jfz(i) ! tzm1
         enddo 
c
c   Input: fr, fs, ft   Output: v
c        dxtm1 fr + (do loop dym1) fs + ft dzm1
         call local_grad3_t(jv,fr,fs,ft,lxd-1,1,dxmd,dxtmd,w2) !!??
c
         call intp_rstd(v, jv, mx, md, if3d, 1) ! backward ( transpose) 
c
         call add2 ( u(1,e),v,lxyz)
c
         enddo 
      else 
c 
         do e=1,nelt
         call rzero(u(1,e),lxyz) ! set to 0 
         call intp_rstd(jfx,fx(1,e),mx,md,if3d,0) ! forward, GLL to GL
         call intp_rstd(jfy,fy(1,e),mx,md,if3d,0) ! forward 
c 
         do i=1,nxyzd
             fr(i) = rx(i,1,e)*jfx(i) ! rxm1, has mass embeded in
     $             + rx(i,2,e)*jfy(i) ! rym1, rx matrix 
             fs(i) = rx(i,3,e)*jfx(i) ! sxm1, refer to induct.f
     $             + rx(i,4,e)*jfy(i) ! sym1
         enddo 
c
         call local_grad2_t(jv,fr,fs,lxd-1,1,dxmd,dxtmd,w2) !!??
c
         call intp_rstd(v, jv, mx, md, if3d, 1) ! backward ( transpose) 
c
         call add2 ( u(1,e),v,lxyz)

         enddo
      endif
c
      return
      end
c-----------------------------------------------------------------------
      subroutine grad_bm1_t(u,fx,fy,fz) !  3+2
c     Alright surface term is not figured out, July 28th
c     New task: get de-aliasing for volume term, July 29th
c                 T
c     Compute grad  of (fx,fy,fz) and set to u.
c
c     Single element case, not de-aliased 
c
c     Input : fx, fy, fz
c     Output: u
c   DxT (M fx)  + (M fy) Dy  + DzT (M fz) 
      include 'SIZE'
      include 'TOTAL'
c
      parameter (lxyz=lx1*ly1*lz1)
      common /ctmp1/ v(lxyz),w(lxyz),w2(lxyz),w3(lxyz)
      real u(lxyz,1),fx(lxyz,1),fy(lxyz,1),fz(lxyz,1)
      real fr(lxyz), fs(lxyz), ft(lxyz)
      integer e
c
      if(if3d) then

         do e = 1,nelt

         call rzero(u(1,e),lxyz) ! set to 0 
         call copy (w, fx(1,e), lxyz) 
         call col2 (w, w3m1,lxyz) ! 
         call copy (w2, fy(1,e), lxyz) 
         call col2 (w2, w3m1,lxyz) ! 
         call copy (w3, fz(1,e), lxyz) 
         call col2 (w3, w3m1,lxyz) ! 

c   what should fr, fs, ft be for fx, fy, fz ? 
         do i=1,lx1*ly1*lz1
             fr(i) = rxm1(i,1,1,e)*w (i) 
     $             + rym1(i,1,1,e)*w2(i)
     $             + rzm1(i,1,1,e)*w3(i) ! jacm1 
             fs(i) = sxm1(i,1,1,e)*w (i) 
     $             + sym1(i,1,1,e)*w2(i)
     $             + szm1(i,1,1,e)*w3(i)
             ft(i) = txm1(i,1,1,e)*w (i) 
     $             + tym1(i,1,1,e)*w2(i)
     $             + tzm1(i,1,1,e)*w3(i)
         enddo 
c   Input: fr, fs, ft   Output: w
c        dxtm1 fr + (do loop dym1) fs + ft dzm1
         call local_grad3_t(v,fr,fs,ft,lx1-1,1,dxm1,dxtm1,w2)
c
         call add2 ( u(1,e),v,lxyz)
c
         enddo
      else 
         do e = 1,nelt

         call rzero(u(1,e),lxyz)       ! set to 0 
         call copy (w,  fx(1,e), lxyz) 
         call col2 (w,  w3m1,lxyz)     ! 
         call copy (w2, fy(1,e), lxyz) 
         call col2 (w2, w3m1,lxyz)     ! 
c 
         do i=1,lx1*ly1*lz1            ! lz1 = 1
             fr(i) = rxm1(i,1,1,e)*w (i) 
     $             + rym1(i,1,1,e)*w2(i)
             fs(i) = sxm1(i,1,1,e)*w (i) 
     $             + sym1(i,1,1,e)*w2(i)
         enddo 
c   Input: fr, fs   Output: w
c        dxtm1 fr + fs dym1
         call local_grad2_t(v,fr,fs,lx1-1,1,dxm1,dxtm1,w3)
c
         call add2 ( u(1,e),v,lxyz)
c
         enddo
      endif
c
      return
      end
c-----------------------------------------------------------------------
c---- End de-aliasing  
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      subroutine dg_heat
c      /----------------------------\
c-----|  dg method for heat problem  |-----
c      \----------------------------/
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'

      rcp = 1.000  ! air: specific heat 1.005 at 25 C, density 1.1~1.2 
      kcd = 0.026  ! air: thermal conductivity 
      kcd = 1.000  ! model problem: take 1 
      kcd = 0.100  ! model problem: take .1 
      prd = 0.711  ! air: prandtl number? can't remember now  

      nu = kcd / rcp       ! k/(rho cp) 

      call dg_advect_setup

      call dg_advance_execute_ht 

      return
      end
c-----------------------------------------------------------------------
      subroutine dg_advance_execute_ht   ! 
c     heat problem 
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'

      real      cfl ! "Diffusion" CFL 
      logical   ifout, ifexc, ifoup 
      real      thex(lx1,ly1,lz1,lelt), terr(lx1,ly1,lz1,lelt) 
      real      temx
      integer   n 

      n=nx1*ny1*nz1*nelt
      ifout = .true.
      ifout = .false.
      ifexc = .false.      ! if there is exact solution in usr file 
      ifexc = .true.       ! if there is exact solution in usr file 

c    Init for dg
      call setlogi_ht      ! 
      call setallic_ht     ! 
      if(nid.eq.0) write(6,*) ' rho C_p',rcp,', k',kcd,', Pr',prd
      if(nid.eq.0) write(6,*) ' nu = k/(rho C_p) ', nu 

      kstep = 0
      if(nid.eq.0) write(6,*) 'Temp field, Step ',kstep
      call printmax_1  (tep) 
      call ht_cfl      (cfl,nu,dt)  ! ! compute cfl 
      if(nid.eq.0) write(6,*) 'Difu. CFL = ', cfl,' , dt', dt   !  

      ifxyo=.true.
      if(ifout) then 
        call outpost(vx,vy,vz,pr,tep,'   ') ! temperature in temperature
        if(nid.eq.0) write(6,*) 'write complete', kstep
        if(ifexc) then ! if there is exact solution in usr file 
          call exc_fld(thex) 
          call sub3(terr,tep,thex,n) 
          ifxyo=.true.
          if (kstep.gt.iostep) ifxyo=.false.
          call outpost(vx,vy,vz,pr,terr,'err') ! 
          if(nid.eq.0) write(6,*) 'err fld write complete'
        endif
      endif

      do kstep=1,nsteps

         call dg_advance_ht(kstep)  ! 
         time = time + dt
         call ht_cfl       (cfl,nu,dt)  ! ! compute cfl 

c        if(mod(kstep,1).eq.0) then
         if(mod(kstep,iostep).eq.0) then

           if(nid.eq.0) write(6,*) 'Temp field, Step ',kstep
     $                            ,', T = ',time
           call printmax_1(tep) 

           if(nid.eq.0) write(6,*) 'Difu. CFL = ', cfl  ! ! 
           if(ifexc) then  ! if there is exact solution in usr file 
             call exc_fld(thex) 
             if(nid.eq.0) write(6,*) 'Exact Temp field, Step ',kstep
     $                              ,', T = ',time
             call printmax_1(thex) 
             call sub3(terr,tep,thex,n)  ! abs error - relative?
             temx = glamax(thex,n)
             temx = 1./temx
c            write(6,*) 'temx', temx,' ,kstep',kstep
             call cmult(terr,temx,n)      ! relative error 
             call print_terr  (terr,kstep,time,dt,nx1) 

             if(ifout) then
               ifxyo=.true.
               if (kstep.gt.iostep) ifxyo=.false.
               call outpost(vx,vy,vz,pr,tep,'   ') ! 
               if(nid.eq.0) write(6,*) 'write complete', kstep
               ifxyo=.true.
               if (kstep.gt.iostep) ifxyo=.false.
               call outpost(vx,vy,vz,pr,terr,'err') ! 
               if(nid.eq.0) write(6,*) 'err fld write complete',kstep
             endif
           endif
         endif
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine setlogi_ht ! 2+3
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'

      tmstpp = 1      ! first order, euler forward 
      tmstpp = 3      ! third order for comparing with Matlab 
      tmstpp = 2      ! second order, check with MATLAB
      htflx  = 1      ! 1 - Central, 2 - ? 
      htflx  = 2      ! 2 - Sym. Interior Penalty, namely primal form 
c     
c     for outpost 
      ifvo = .true.
      ifpo = .true.
      ifto = .true.
c     

      if(nid.eq.0) then 
        write(6,*) 'o------------------------------------------------o'
        write(6,*) '| Heat problem                    '
        write(6,*) 'o------------------------------------------------o'
        write(6,*) '| Time-stepper order is           ', tmstpp
        write(6,*) 'o------------------------------------------------o'
        write(6,*) '| Flux type is ( 1-Ctr, 2-SIPG)   ', htflx
        write(6,*) 'o------------------------------------------------o'
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine setallic_ht ! 2+3
      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'
      include 'DGUSE'

      n = nx1*ny1*nz1*nelt
c     write(6,*) ifgetx, ifgetu, ifgetp, ifgett 
c     this if is working. for ifdg flag turned on 
      if(       ifgetx             ! if true means mesh restart
     $    .or.  ifgett      ) then ! if true means ener restart
          ! copy into my data fields, time is set in restarting
          if(if3d) then 
              if(nid.eq.0) write(6,*) 'Done :: DG Restarting 3d...'
              call copy(tep,t,n) 
          else
              call copy(tep,t,n) 
              if(nid.eq.0) write(6,*) 'in restart '
              call printmax_1(tep) 
              if(nid.eq.0) write(6,*) 'Done :: DG Restarting 2d...'
          endif
      else ! not restarting 
          call setic_ht ! ! 
          time = 0.
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine setic_ht ! 2+3
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      integer n, nel, nxyz, ie, i, j, k

      n = nx1*ny1*nz1*nelt
      nel = nelfld(ifield)
      nxyz=nx1*ny1*nz1
      do ie=1,nel
          ieg  = lglel(ie)
          do k=1,nz1
          do j=1,ny1
          do i=1,nx1
              call nekasgn(i,j,k,ie)
              call useric_ht(i,j,k,ie) ! get rho, enr, rux,.. ruz
              tep(i,j,k,ie)= tem 
          enddo 
          enddo 
          enddo 
      enddo 
      if(nid.eq.0) write(6,*) 'done DG :: Set init for ht'

      return
      end
c-----------------------------------------------------------------------
c----- Copied from induct.f, subroutine compute_cfl 
c--------------------------------------------------------------------
      subroutine ht_cfl(cfl,nu,dt)  ! ! compute diffusion cfl!NOT RIGHT 
c
c     Given kinematic viscosity nu, compute cfl = nu/dx^2 
c
      include 'SIZE'
      include 'GEOM'
      include 'INPUT'
      include 'WZ'
      include 'SOLN'
c
      real nu, dt  
      real cfl, cflr, cfls, cflt, cflm 
c
c     Store the inverse jacobian to speed up this operation
c
      common /cfldx/ dri(lx1),dsi(ly1),dti(lz1)
c
      integer e
c
      integer icalld
      save    icalld
      data    icalld /0/
c
      if (icalld.eq.0) then
         icalld=1
         call getdr(dri,zgm1(1,1),nx1) ! inverse dr 
         call getdr(dsi,zgm1(1,2),ny1)
         if (if3d) call getdr(dti,zgm1(1,3),nz1)
      endif

      cfl = 0.
      l   = 0

      if (if3d) then
         nxyz = nx1*ny1*nz1
         do e=1,nelv
            do k=1,nz1
            do j=1,ny1
            do i=1,nx1
               l = l+1

               ur = ( 1.*rxm1(i,j,k,e)
     $            +   1.*rym1(i,j,k,e)
     $            +   1.*rzm1(i,j,k,e) ) * jacmi(l,1)
               us = ( 1.*sxm1(i,j,k,e)
     $            +   1.*sym1(i,j,k,e)
     $            +   1.*szm1(i,j,k,e) ) * jacmi(l,1)
               ut = ( 1.*txm1(i,j,k,e)
     $            +   1.*tym1(i,j,k,e)
     $            +   1.*tzm1(i,j,k,e) ) * jacmi(l,1)

               cflr = abs(dt*nu*(ur*dri(i))**2)  !! 
               cfls = abs(dt*nu*(us*dsi(j))**2)  !! 
               cflt = abs(dt*nu*(ut*dti(k))**2)  !! 

               cflm = cflr + cfls + cflt 
               cfl  = max(cfl,cflm)

               cflf(i,j,k,e) = cflm ! in SOLN 
 
            enddo
            enddo
            enddo
         enddo
      else
         nxyz = nx1*ny1
         do e=1,nelv
            do j=1,ny1
            do i=1,nx1
               l = l+1

               ur = ( 1.*rxm1(i,j,1,e)
     $            +   1.*rym1(i,j,1,e) ) * jacmi(l,1)
               us = ( 1.*sxm1(i,j,1,e)
     $            +   1.*sym1(i,j,1,e) ) * jacmi(l,1)

               cflr = abs(dt*nu*(ur*dri(i))**2)  !! 
               cfls = abs(dt*nu*(us*dsi(j))**2)  !! 

               cflm = cflr + cfls 
               cfl  = max(cfl,cflm)

               cflf(i,j,1,e) = cflm   ! in SOLN 

            enddo
            enddo
         enddo
      endif
c
      cfl = glmax(cfl,1)
c
      return
      end
c-----------------------------------------------------------------------
      subroutine dg_advance_ht(kstep) ! Advance u^{n-1} --> u^n
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'

      integer kstep

      if(tmstpp.eq.1) then      ! euler 

         call dg_advance_ht_euler1(kstep)

      else if(tmstpp.eq.2) then ! rk2

         call dg_advance_ht_ssprk2(kstep)

      else if(tmstpp.eq.3) then ! rk3

         call dg_advance_ht_ssprk3(kstep)

      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine dg_advance_ht_euler1(kstep) ! Advance u^{n-1} --> u^n
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
c
      real    t1
      real    rhst(lt)
c
      n  = nx1*ny1*nz1*nelt
      t1 = time
c 
      call copy           (tmtp, tep, n)   ! copy into scratch
      call comp_dg_rhs_ht (rhst, t1)       ! compute rhs 
c
      do i=1,n                             ! update solution 
          tep(i,1,1,1) = tep(i,1,1,1) + rhst(i)*dt 
      enddo
c     write(6,*) '7' 
c     call outfldz1(tep,'tep fin   ',1)  
c     stop 
c
      return
      end
c-----------------------------------------------------------------------
      subroutine dg_advance_ht_ssprk2(kstep) ! Advance u^{n-1} --> u^n
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
c
      real    t1, t2
      real    tep1(lt), rhst(lt)
      n = nx1*ny1*nz1*nelt
c
c     Stage 1 
c     write(6,*) 'Stage 1 '
      t1 = time
      call copy(tmtp, tep, n)                     ! copy into scratch
      call comp_dg_rhs_ht (rhst,t1) ! compute rhs 
      do i=1,n                                    ! update solution 
          tep1 (i) = tep(i,1,1,1) + rhst(i)*dt
      enddo
c
c     Stage 2 
c     write(6,*) 'Stage 2 ' ! 
      t2 = time + dt
      call copy(tmtp, tep1, n)                    ! copy into scratch
      call comp_dg_rhs_ht (rhst,t2) ! compute rhs 
      do i=1,n                                    ! update solution 
          tep(i,1,1,1) = (tep(i,1,1,1) + tep1(i) + rhst(i)*dt)/2.
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine dg_advance_ht_ssprk3(kstep) ! Advance u^{n-1} --> u^n
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
c
      real    t1, t2, t3
      real    tep1(lt), tep2(lt), rhst(lt)
      n = nx1*ny1*nz1*nelt
c
c     Stage 1 
c     write(6,*) 'Stage 1 '
      t1 = time
      call copy(tmtp, tep, n)                     ! copy into scratch
      call comp_dg_rhs_ht (rhst,t1) ! compute rhs 
      write(6,*) 'before stage 1 in rk3'

      do i=1,n       ! update solution 
          tep1 (i) = tep(i,1,1,1) + rhst(i)*dt
      enddo
      write(6,*) 'stage 1 in rk3'
c
c     Stage 2 
c     write(6,*) 'Stage 2 ' ! stage 2 breaks down 
      t2 = time + dt
      call copy(tmtp, tep1, n)                     ! copy into scratch
      call comp_dg_rhs_ht (rhst,t2) ! compute rhs 
      do i=1,n       ! update solution 
          tep2 (i) = (tep(i,1,1,1)*3. + tep1(i) + rhst(i)*dt)/4.
      enddo
      write(6,*) 'stage 2 in rk3'
c
c     Stage 3 
c     write(6,*) 'Stage 3 '
      t3 = time + dt/2.
      call copy(tmtp, tep2, n)                     ! copy into scratch
      call comp_dg_rhs_ht (rhst,t3) ! compute rhs 
      do i=1,n       ! update solution 
          tep (i,1,1,1)=(tep(i,1,1,1)+tep2(i)*2.+rhst(i)*dt*2.)/3.
      enddo
      write(6,*) 'stage 3 in rk3'
      stop 
c
      return
      end
c-----------------------------------------------------------------------
      subroutine comp_dg_rhs_ht(rhs,tm) ! find right hand side 
c     Input: 
c         . tmtp 
c     Output: 
c         . rhs 
c
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'   ! Read htflx from here
      include 'DG'      ! dg_face is stored
c
      real rhs(1), tm 
      n  = nx1*ny1*nz1*nelt
c 
      call rzero(rhs,n)
c 
      if(htflx.eq.1) then         ! Central 
          call comp_dg_rhs_ht_ctr(rhs,tm) 
      elseif(htflx.eq.2) then     ! Primal 
c         write(6,*) 'before getting into'
          call comp_dg_rhs_ht_prm(rhs,tm) 
c         write(6,*) 'at least u should be here'
      endif
c 
      return
      end
c-----------------------------------------------------------------------
c----| SIP/ primal form |-----------------------------------------------
c-----------------------------------------------------------------------
      subroutine comp_dg_rhs_ht_prm(rhs,tm) ! rhs for htflx 2 
c
c    rhs =  I1 - I2 - I3 + I4 
c     Input: 
c         . tmtp 
c     Output: 
c         . rhs - first (neg) volume res, then vol + surf 
c     Interm: 
c         . rsf - surf 
c
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      real     rhs(lt), tm
      real     ku(lt) , hu(lt)
     $       , gtu(lt), gu(lt)
c
      n  = nx1*ny1*nz1*nelt
      nf = nx1*nz1*2*ndim*nelt
c   
      call iku (ku)      ! + grad v \cdot grad u    ! 
c     write(6,*) ' 1 .'
c   
      call igtu(gtu)     ! - [u] \cdot grad v       ! 
c     write(6,*) ' 2 .'
c     
      call igu (gu)      ! - [v] \cdot grad u       ! 
c     write(6,*) ' 3 .'  
c     
c         gu is missing factor of 4 
c        gtu prbbly missing that as well
c        reason: dxm1 does not have geom/jacob
c        in axhelm, it is dealt with by 
c        g1m1, etc, line 153 hmholtz.f, axhelm 
c 
      call ihu (hu)      ! + eta/ h ([v] \cdot [u]) ! 
c     write(6,*) ' 4 .'
c     
c     write(6,*) 'before summing them'
      call add3  (rhs,ku,hu,n) !
      call chsign(gtu,n)       !
      call add2  (rhs,gtu,n)   !
      call chsign(gu,n)        !
      call add2  (rhs,gu,n)    !
c                              !
      call chsign(rhs,n)       ! - sign
      call cmult (rhs,nu,n)    ! everything times nu 
c
      call invbf (rhs,rhs)     !
c     write(6,*) 'check everythin'
c     stop 
c 
      return
      end
c-----------------------------------------------------------------------
      subroutine iku(ku) ! 
c     Vol integral 
c     Input: 
c         . tmtp 
c     Output: 
c         . ku
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      parameter(le=lx1*ly1*lz1)
      real     grdu(le,lelt,ldim)  ! second index for dimension 
      real     ku  (le,lelt)  ! second index for dimension 
      real     h1(le,lelt), h2(le,lelt) 
      integer  e, n, npl, nf 
      integer  imesh, isd 

      npl = nx1 - 1
      n  = nx1*ny1*nz1*nelt
      nf = nx1*nz1*2*ndim*nelt
      nfaces=2*ndim
      nxz=nx1*nz1

c     grad u 
c     Src: . navier5.f : local_grad3; gradm1 
c          . hmholtz.f : axhelm ? 
c   . Approach 1 for step 1 
c      call gradm1(grdu(1,1,1),grdu(1,1,2),grdu(1,1,3),tmtp) 
c     ! this converts to xyz coord., not directly useful for local_grad3_t 
c   . Approach 2 for step 1 
c     do e=1,nelt 
c       if(if3d) then
c         call local_grad3(grdu(1,e,1),grdu(1,e,2),grdu(1,e,3) 
c    $                    ,tmtp,np,e,dxm1,dxtm1) 
c       else 
c         call local_grad2(grdu(1,e,1),grdu(1,e,2),
c    $                    ,tmtp,np,e,dxm1,dxtm1) 
c       endif 
c     enddo
c    B grad u 
c     do e=1,nelt 
c       call ! multi mass matrix 
c     enddo 
c    div B grad u 

c   . Approach 3 
c     call axhelm, which should yield:  h1 A u + h2 B u
      imesh = 1 
      call rone (h1,n) 
      call rzero(h2,n) 

      isd   = 1   ! What is isd? 
      call axhelm(ku,tmtp,h1,h2,imesh,isd) 

      return
      end
c-----------------------------------------------------------------------
      subroutine igtu(gtu) ! 
c     Vol integral 
c     Input: 
c         . tmtp 
c     Output: 
c         . gtu
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      parameter(le=lx1*ly1*lz1)
      real     gtu  (le,lelt)         
      real     uf   (lf), ufsv (lf)
      real     ufxyz(lf,ldim), uxyz(le,lelt,ldim) 
      real     tmp  (le)
      integer  e, n, npl, nf, f, i, k
c
      n  = nx1*ny1*nz1*nelt
      npl = nx1 - 1
      nf = nx1*nz1*2*ndim*nelt
      nfaces=2*ndim
      nxz=nx1*nz1
c     
c     Src: . navier1.f : local_grad3_t 
c   . R
      call full2face(uf,tmtp)
c   . diff 
      call copy(ufsv,uf,nf) 
c     
c     P(eriodic), do nothing 
c     call add2  (dvgf,bf1,nf)         ! P, implied in next step, a- + a+(bc)
      call gs_op (dg_hndl,uf,1,1,0)    ! 1 ==> +      , a- + a+(itr)
      call cmult (uf,0.5,nf)           ! times 2.
      call chsign(uf,nf)               ! - ( a + b ) 
      call add2  (ufsv,uf,nf)          !  
c   . nhat, area 
      do idm=1,ldim
         call rzero(ufxyz(1,idm),nf) 
      enddo
      k = 0
      do e=1,nelt
      do f=1,nfaces
      do i=1,nxz
         k=k+1
         if(if3d) then
            ufxyz(k,1) = ufsv(k)*unx(i,1,f,e)*area(i,1,f,e)
            ufxyz(k,2) = ufsv(k)*uny(i,1,f,e)*area(i,1,f,e)
            ufxyz(k,3) = ufsv(k)*unz(i,1,f,e)*area(i,1,f,e)
         else ! 2D, 
            ufxyz(k,1) = ufsv(k)*unx(i,1,f,e)*area(i,1,f,e)
            ufxyz(k,2) = ufsv(k)*uny(i,1,f,e)*area(i,1,f,e)
         endif
      enddo
      enddo
      enddo
c   . back to volume 
      do idm=1,ndim
        call face2full(uxyz(1,1,idm),ufxyz(1,idm))        !  
      enddo
c   . grad T (ux, uy, uz) 
c    . Approach 1 
c     do e=1,nelt 
c       if(if3d) then
c         call local_grad3_t(gtu,uxyz(1,e,1),uxyz(1,e,2),uxyz(1,e,3)
c    $                    ,npl,e,dxm1,dxtm1,tmp) 
c       else 
c         call local_grad2_t(gtu,uxyz(1,e,1),uxyz(1,e,2)
c    $                    ,npl,e,dxm1,dxtm1,tmp) 
c       endif 
c     enddo
c    . Not done yet! Multiple with rxm1, divide by Jacm1 
c     - so we are back to gradm1 - too bad there is no gram1_T 
c    . Approach 2 
      call gradm1_t(gtu,uxyz(1,1,1),uxyz(1,1,2),uxyz(1,1,3)) 

c     write(6,*) 'done in gtu' 

      return
      end
c-----------------------------------------------------------------------
      subroutine gradm1_t(u,ux,uy,uz)
c     source: . gradm1, from navier5.f 
c             . gradm1, from navier5.f 
c
c     Compute divergence of ux,uy,uz -- mesh 1 to mesh 1 (vel. to vel.)
c
      include 'SIZE'
      include 'DXYZ'
      include 'GEOM'
      include 'INPUT'
      include 'TSTEP'
c
      parameter (lxyz=lx1*ly1*lz1)
      real ux(lxyz,1),uy(lxyz,1),uz(lxyz,1),u(lxyz,1)

c     common /ctmp1/ ur(lxyz),us(lxyz),ut(lxyz),tmp(lxyz)
      real ur(lxyz),us(lxyz),ut(lxyz),tmp(lxyz)

      integer e

      nxyz = nx1*ny1*nz1
      ntot = nxyz*nelt

      N = nx1-1
      do e=1,nelt
         if (if3d) then

            do i=1,lxyz
               ur(i) = jacmi(i,e)*(ux(i,e)*rxm1(i,1,1,e)
     $                           + uy(i,e)*rym1(i,1,1,e)
     $                           + uz(i,e)*rzm1(i,1,1,e) )
               us(i) = jacmi(i,e)*(ux(i,e)*sxm1(i,1,1,e)
     $                           + uy(i,e)*sym1(i,1,1,e)
     $                           + uz(i,e)*szm1(i,1,1,e) )
               ut(i) = jacmi(i,e)*(ux(i,e)*txm1(i,1,1,e)
     $                           + uy(i,e)*tym1(i,1,1,e)
     $                           + uz(i,e)*tzm1(i,1,1,e) )
            enddo
            call local_grad3_t(u,ur,us,ut,N,e,dxm1,dxtm1,tmp)
         else
            do i=1,lxyz
               ur(i) =jacmi(i,e)*(ux(i,e)*rxm1(i,1,1,e)
     $                          + uy(i,e)*rym1(i,1,1,e) )
               us(i) =jacmi(i,e)*(ux(i,e)*sxm1(i,1,1,e)
     $                          + uy(i,e)*sym1(i,1,1,e) )
            enddo
            call local_grad2_t(u,ur,us,N,e,dxm1,dytm1,tmp)
         endif
      enddo
c
      return
      end
c-----------------------------------------------------------------------
      subroutine igu(gu) ! 
c     Vol integral 
c     Input: 
c         . tmtp 
c     Output: 
c         . gu
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      parameter(le=lx1*ly1*lz1)
      real     grdu(le,lelt,3), grduf(lf,3) 
      real     dvgf(lf), dvsv(lf)     ! div \cdot grad u, on surface
      real     gu  (le,lelt)         ! second index for dimension 
      real     h1(le,lelt), h2(le,lelt) 
      integer  e, n, npl, nf, f, i, k, idm
c
      n  = nx1*ny1*nz1*nelt
      nf = nx1*nz1*2*ndim*nelt
      nfaces=2*ndim
      nxz=nx1*nz1
      npl = nx1 - 1
c     
c     Src: . navier5.f : local_grad3 
c   . grad u 
c    . Approach 1 - this works ? NO! y derivative not right 
      call gradm1(grdu(1,1,1),grdu(1,1,2),grdu(1,1,3),tmtp) 
c    . Approach 2 
c     do e=1,nelt 
c       if(if3d) then
c         call local_grad3(grdu(1,e,1),grdu(1,e,2),grdu(1,e,3) 
c    $                    ,tmtp,npl,e,dxm1,dxtm1) 
c       else 
c         call local_grad2(grdu(1,e,1),grdu(1,e,2)
c    $                    ,tmtp,npl,e,dxm1,dxtm1) 
c       endif 
c     enddo
c     not done yet! Need to multi. rxm1 then divide by jacm1
c 
c   . R 
      do idm=1,ndim
        call full2face(grduf(1,idm),grdu(1,1,idm)) 
      enddo
c   . dot nhat 
      call rzero(dvgf,nf) 
      k = 0
      do e=1,nelt
      do f=1,nfaces
      do i=1,nxz
         k=k+1
         if(if3d) then
            dvgf(k) = grduf(k,1)*unx(i,1,f,e)*area(i,1,f,e)
     $              + grduf(k,2)*uny(i,1,f,e)*area(i,1,f,e)
     $              + grduf(k,3)*unz(i,1,f,e)*area(i,1,f,e)
         else ! 2D, 
            dvgf(k) = grduf(k,1)*unx(i,1,f,e)*area(i,1,f,e)
     $              + grduf(k,2)*uny(i,1,f,e)*area(i,1,f,e)
         endif
      enddo
      enddo
      enddo
c   . diff 
      call copy  (dvsv,dvgf,nf) 
c     
c     the boundary thing is still really ugly 
c     ok for periodic bc we should be able to get away with 
c     do nothing 
c     call add2  (dvgf,bf1,nf)         ! P, implied in next step, a- + a+(bc)
      call gs_op (dg_hndl,dvgf,1,1,0)  ! 1 ==> +      , a- + a+(itr)
      call cmult (dvgf,0.5,nf)         ! times 2.
      call chsign(dvgf,nf)             ! - ( a + b ) 
      call add2  (dvsv,dvgf,nf)        !  
c   . back to volume 
      call face2full(gu,dvsv)        !  
c     write(6,*) 'done in gu' 

      return
      end
c-----------------------------------------------------------------------
      subroutine set_eta_he(eta,he) ! 
c     Vol integral 
c     Input: 
c         . geometry and weights 
c     Output: 
c         . eta, he 
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      parameter(le=lx1*ly1*lz1)
      real     eta (lf), he(lf), tmpw 
      integer  e, n, npl, nf, f, i, k
c
      n  = nx1*ny1*nz1*nelt
      nfaces=2*ndim
      nxz=nx1*nz1
      nf = nx1*nz1*2*ndim*nelt
c     
c   . Set eta and h ! Yep this awkward, but seems to work 
      tmpw = 1./sqrt(bm1(1,1,1,1)) 
c     write(6,*) 'bm1 (1) ', bm1(1,1,1,1), ', tmpw', tmpw
c     stop 
      call rone (eta,nf) 
      call cmult(eta,tmpw,nf) 

      k = 0
      do e=1,nelt
      do f=1,nfaces

      sumle = 0.
      do i1=1,nxz
         sumle = sumle + area(i1,1,f,e) 
      enddo

      do i=1,nxz
         k=k+1
         he(k) = sumle 
      enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine ihu(hu) ! 
c     Vol integral 
c     Input: 
c         . tmtp 
c     Output: 
c         . hu
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      parameter(le=lx1*ly1*lz1)
      real     hu  (le,lelt), huf(lf)
      real     uf  (lf), ufsv(lf)
      real     eta (lf), he(lf)
      integer  e, n, npl, nf, f, i, k
c
      n  = nx1*ny1*nz1*nelt
      nfaces=2*ndim
      nxz=nx1*nz1
      npl = nx1 - 1
      nf = nx1*nz1*2*ndim*nelt
c     
c   . Set eta and h 
      call set_eta_he(eta,he) ! ! constants 
c   . R
      call full2face(uf,tmtp)
c   . diff 
      call copy  (ufsv,uf,nf) 
      call cmult (ufsv,2.0,nf)         ! times 2.
c     
c     P(eriodic), do nothing 
c     call add2  (dvgf,bf1,nf)         ! P, implied in next step, a- + a+(bc)
      call gs_op (dg_hndl,uf,1,1,0)    ! 1 ==> +      , a- + a+(itr)
      call chsign(uf,nf)               ! - ( a + b ) 
      call add2  (ufsv,uf,nf)          !  
c   . area, eta, he 
      k = 0
      do e=1,nelt
      do f=1,nfaces
      do i=1,nxz
         k=k+1
c        if(if3d) then
            huf(k) = ufsv(k)*area(i,1,f,e)*eta(k)/he(k) 
c        else ! 2D, 
c        endif
      enddo
      enddo
      enddo
c   . back to volume 
      call face2full(hu,huf)           !  
c     write(6,*) 'doen in ihu'

      return
      end
c-----------------------------------------------------------------------
c----| Central flux, with aux. solve  |---------------------------------
c-----------------------------------------------------------------------
      subroutine comp_dg_rhs_ht_ctr(rhs,tm) ! rhs for htflx 1 
c
c    rhs =  volume + surface
c     Input: 
c         . tmtp 
c     Output: 
c         . rhs - first (neg) volume res, then vol + surf 
c     Interm: 
c         . sht - auxiliary var. 
c         . rsf - surf 
c
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      real     rhs(1), tm
      real     rsf(lt)
c
      n  = nx1*ny1*nz1*nelt
      nf = nx1*nz1*2*ndim*nelt
c 
      call rzero(rsf,n) 
c     write(6,*) '1' 
c     call outfldz1(tmtp,'tmtp     1',1)  ! 

c   Auxiliary solve , boundary routine not complete yet 
      call ht_aux_slv  (sht)       ! sht: aux var. for heat 
c     write(6,*) '2' 
c     call outfldz1(sht(1,1),'sht 1 at 2',1)  
c     call outfldz1(sht(1,2),'sht 2 at 2',1)  
c     stop 

c   Volume 
      call ht_vol_eval             ! evaluate components, for vol int. 
c     write(6,*) '3, nothing changes' 
c     call outfldz1(sht(1,1),'sht 1 at 3',1)  
c     call outfldz1(sht(1,2),'sht 2 at 3',1)  

      call ht_vol_res  (rhs)       ! evaluate volume residue, with - 
c     write(6,*) '4'   
c     call outfldz1(rhs,'rhs      4',1)  

c   Surface 
      call ht_surf_res (rsf)       ! evaluate surface residue 
c     write(6,*) '5' 
c     call outfldz1(rsf,'rsf      5',1)  

      call add2        (rhs,rsf,n) ! DONE 
c     write(6,*) '6' 
c     call outfldz1(rhs,'rhs fin  6',1)  
c     stop 
c 
      return
      end
c-----------------------------------------------------------------------
      subroutine ht_aux_slv(tmq) ! 
c    Solve for aux variable 
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      real     tmq(lt,1) 
      real     vua(lt,ldim) 
      real     uaf(lf,ldim) 

      n  = nx1*ny1*nz1*nelt
      nf = nx1*nz1*2*ndim*nelt

c    Vol integral 
      call ht_aux_vol (vua)             ! 
c     write(6,*) '2.1' 
c     call outfldz1(vua(1,1),'vua1 in  2',1)  
c     call outfldz1(vua(1,2),'vua2 in  2',1)  

c    Surf integral 
      call ht_aux_srf (uaf)             ! 
c     write(6,*) '2.2' 
c     call outfacez1(uaf(1,1),'uaf1 in  2',1)  
c     call outfacez1(uaf(1,2),'uaf2 in  2',1)  

c    Solve for aux 
      call ht_aux_stp (tmq, vua, uaf)   ! 
c     write(6,*) '2.3' 
c     call outfldz1(tmq(1,1),'tmq1 in  2',1)  
c     call outfldz1(tmq(1,2),'tmq2 in  2',1)  
c     call outfldz1(bm1,  'bm1       ',1)  
c     call outfldz1(jacm1,'jacm1     ',1)  
c     call outfldz1_e1(w3m1,'w3m1      ')  
c     write(6,*) wxm1
c     write(6,*) wym1
c     stop 

      return
      end
c-----------------------------------------------------------------------
      subroutine ht_aux_vol(volq) ! 
c    
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      parameter(le=lx1*ly1*lz1)
      real     volq(lt,1)  ! second index for dimension 
      real     mtp(le,lelt) 
      integer  e 

      n  = nx1*ny1*nz1*nelt
      nf = nx1*nz1*2*ndim*nelt

c    Vol integral 
c    Temperature in tmtp field 
      do e=1,nelt 
        call col3 (mtp(1,e),tmtp(1,1,1,e),w3m1,le) ! multi mass matrix 
      enddo 

c     write(6,*) '2.1.1' 
c     call outfldz1(mtp,'mtp  in   ',1)  
      if(if3d) then  ! using routine from diffusion part ? 
        call grad_t_u (volq(1,1),volq(1,2),volq(1,3),mtp) ! transpose of grad
      else 
        call grad_t_u2(volq(1,1),volq(1,2),mtp)
      endif

      return
      end
c-----------------------------------------------------------------------
      subroutine ht_aux_srf(srfq) ! 
c    
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      real     srfq(lf,ldim) 
      integer  e, f, i, k
      real     ctrq(lf), difq(lf)

      n  = nx1*ny1*nz1*nelt
      nf = nx1*nz1*2*ndim*nelt
      nfaces=2*ndim
      nxz=nx1*nz1

c    Surf integral 
      call ht_aux_ctr(ctrq)     !  central 
c     call ht_aux_dif(difq)     !  difference, can wait a bit 
c     write(6,*) '2.2.1, should be the same across surf ' 
c     call outfacez1(ctrq,'ctrq in  2',1)  
c     yep they are the same 

      k = 0
      do e=1,nelt
      do f=1,nfaces
      do i=1,nxz
         k=k+1
         if(if3d) then
            srfq(k,1) = ctrq(k)*unx(i,1,f,e)*area(i,1,f,e)
            srfq(k,2) = ctrq(k)*uny(i,1,f,e)*area(i,1,f,e)
            srfq(k,3) = ctrq(k)*unz(i,1,f,e)*area(i,1,f,e)
         else ! 2D, just central flux  
            srfq(k,1) = ctrq(k)*unx(i,1,f,e)*area(i,1,f,e)
            srfq(k,2) = ctrq(k)*uny(i,1,f,e)*area(i,1,f,e)
         endif
      enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine ht_aux_stp ( tmq, vua, uaf)  ! ! 
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      integer  i, f, e, k, id 
      real     tmq(lt,1)  ! def. on volume 
      real     vua(lt,1)
      real     uaf(lf,1)
      real     tmp(lt)
c
      n=nx1*ny1*nz1*nelt
c
      do id=1,ndim  ! 3, 2 
        call face2full (tmp,uaf(1,id)) 
        call sub2      (tmp,vua(1,id),n)  ! surf - vol 
c    
        call invbf     (tmq(1,id), tmp)   ! inv bm1 
      enddo 
c
      return
      end
c-----------------------------------------------------------------------
      subroutine ht_aux_ctr(ctrq) !  ! central 
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      real     btf(lf), tpf (lf) 
      real     ctrq(1)

      n=nx1*ny1*nz1*nelt
      nf=nx1*nz1*2*ndim*nelt !total number of points on faces

      call full2face  (tpf, tmtp) 
      call rzero      (btf,nf) 
      call htbc_aux   (btf)           ! boundary for heat 

      call add2  (tpf,btf,nf)         ! first add bc , a- + a+(bc) 
      call gs_op (dg_hndl,tpf,1,1,0)  ! 1 ==> +      , a- + a+(itr)
      call cmult (tpf,0.5,nf)         ! times 1/2, (a- + a+)/2 

      call copy  (ctrq,tpf,nf)        ! just a copy 

      return
      end
c-----------------------------------------------------------------------
      subroutine htbc_aux(btf) ! 3+2
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      common /nekcb/ cb
      character cb*3
      integer  i, f, e, k
      real     btf(lf)
      n   =nx1*ny1*nz1*nelt
      nf  =nx1*nz1*2*ndim*nelt 
      nxz =nx1*nz1
      nfaces=2*ndim
      k=0

c     try to fix the thermal boundary thing by setting ifield
c     seems ifield = 2 will work 
      ifield=2 

      do e=1,nelt
      do f=1,nfaces
         ieg=lglel(e)
         cb =cbc(f,e,ifield)
c        write(6,*)'ifield', ifield 
c        write(6,*)'in the aux sovle, cb = ', cb
         if(cb.eq.'t  ' .or. cb.eq.'T  ' .or.  ! temp 
     $      cb.eq.'f  ' .or. cb.eq.'F  '       ! flux 
     $      ) then !  
            ia=0
            call facind(kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,f)
            do iz=kz1,kz2
            do iy=ky1,ky2
            do ix=kx1,kx2
               ia = ia + 1
               k  = ia + (f-1)*nxz + (e-1)*nfaces*nxz
               if(cb.eq.'t  '.or. cb .eq.'T  ') then  ! Thermal dirich
c     According to the one dimensional thing, on Dir, up = 2*ud - um 
                   call nekasgn(ix,iy,iz,e)
                   call userbc_t(ix,iy,iz,f,e,ieg) ! temperature 
                   btf(k) = 2.*tem - tmtp(ix,iy,iz,e)  !! 
               else if(cb.eq.'f  ' .or. cb.eq.'F  ' ) then ! ! solid wall 
c     According to the one dimensional thing, on Neu, up =  um 
                   call nekasgn(ix,iy,iz,e)
                   call userbc_f(ix,iy,iz,f,e,ieg) ! rho, rux, ruy..
                   btf(k) = tmtp(ix,iy,iz,e)           !! 
               endif
            enddo 
            enddo 
            enddo 
         endif
      enddo 
      enddo 

      return
      end
c-----------------------------------------------------------------------
c----- Bottom of Boundary 
c-----------------------------------------------------------------------
      subroutine ht_vol_eval ! 3+2? 
c
c --\
c    - Volume eval. Multi. nv to aux q := grad T 
c --/
c
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      parameter(le=lx1*ly1*lz1)
      integer  e, id 

      n=nx1*ny1*nz1*nelt
      do id=1,ndim  ! 3, 2, sht def. in DGUSE 
         call cmult(sht(1,id), nu, n) 
      enddo 
      
      return
      end
c-----------------------------------------------------------------------
      subroutine ht_surf_res(rsf) ! ! 3+2? 
c --\
c    - Surface residue. Central now? simplest one 
c --/
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      real     rsf(1)
      real     flx(lf)
      real     tm1(lt)

      n=nx1*ny1*nz1*nelt
      ne=nx1*ny1*nz1
      
      call ht_flx_srf (flx)      ! save into flx array  
c     write(6,*) '4.1' 
c     call outfacez1(flx,'flx    4-5',1)  

      call face2full  (tm1,flx)  ! 
      call invbf      (rsf,tm1)  ! all elements 
c     write(6,*) '4.2' 
c     call outfldz1(rsf,'rsf    4.2',1)  

      return
      end
c-----------------------------------------------------------------------
c----- Surface flx routines 
c-----------------------------------------------------------------------
      subroutine ht_flx_srf(flx) ! 3+2
c     Lax-Friedrichs flux, only weak form now 
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      integer  i, f, e, k
      real     flx(1)      ! defined on surface 
      real     ctrh(lf,3)  ! defined on surface, ndim components 
     $       , dift(lf)    ! defined on surface
     $       , tau (lf)    ! defined on surface

      n=nx1*ny1*nz1*nelt
      nxz=nx1*nz1
      nfaces=2*ndim

c    Surf integral 
      call ht_flx_ctr(ctrh) ! central flux of h 
c     write(6,*) '4.1.1' 
c     call outfacez1(ctrh(1,1),'ctrh1     ',1)  
c     call outfacez1(ctrh(1,2),'ctrh2     ',1)  

      call ht_flx_dif(dift) ! difference in temp
      call ht_flx_tau(tau ) ! eval tau 

      k = 0
      do e=1,nelt
      do f=1,nfaces
      do i=1,nxz
         k=k+1
         if(if3d) then
            flx(k) =  ctrh(k,1)*unx(i,1,f,e)
     $              + ctrh(k,2)*uny(i,1,f,e)
     $              + ctrh(k,3)*unz(i,1,f,e)
     $              - tau(k)*dift(k)/2. 
            flx(k) = flx(k)*area(i,1,f,e)
         else ! 2D 
            flx(k) =  ctrh(k,1)*unx(i,1,f,e)
     $              + ctrh(k,2)*uny(i,1,f,e)
     $              - tau(k)*dift(k)/2. 
            flx(k) = flx(k)*area(i,1,f,e)
         endif
      enddo
      enddo
      enddo

      return
      end
c-----------------------------------------------------------------------
      subroutine ht_flx_tau(tau) ! eval tau 
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      real     tau(1)   ! defined on surface 

      n=nx1*ny1*nz1*nelt
      nf=nx1*nz1*2*ndim*nelt !

      call rone       (tau, nf)
      call cmult      (tau, 0.5 , nf)      ! Hack to be 0.5 for now 

      return
      end
c-----------------------------------------------------------------------
      subroutine ht_flx_dif(dift) !  ! difference in tem
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      real     dift(1)
      real     ftsv(lf), btf  (lf) 

      n=nx1*ny1*nz1*nelt
      nf=nx1*nz1*2*ndim*nelt !

      call full2face  (dift, tmtp)          !

      call rzero      (btf, nf)
      call htbc_aux   (btf)                 ! 3+2

      call copy       (ftsv, dift, nf)      !
      call cmult      (ftsv, 2. , nf)       !

      call add2       (dift,btf,nf)         ! first add bc 
      call gs_op      (dg_hndl,dift,1,1,0)  ! 1 ==> + , a + b 
      call chsign     (dift,nf)             ! - ( a + b ) 
      call add2       (dift,ftsv,nf)        ! 2 a - ( a + b ) = a - b

      return
      end
c-----------------------------------------------------------------------
      subroutine ht_flx_ctr(ctrh) !  ! central 
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      real     btf(lf),  btg(lf), bth(lf), tpf(lf)
      real     ctrh(lf,1)

      n=nx1*ny1*nz1*nelt
      nf=nx1*nz1*2*ndim*nelt !total number of points on faces

      call copy    (fd, sht(1,1), n) 
      call copy    (gd, sht(1,2), n) 
      call rzero   (hd, n) 
c     
c     write(6,*) '4.1.1.1 ' 
c     call outfldz1(fd,'fd        ',1)  
c     call outfldz1(gd,'gd        ',1)  

      call full2face  (fdf, fd) 
      call full2face  (gdf, gd) 
      if(if3d)   call full2face  (hdf, hd) 

      call htbc_tmp   (btf,btg,bth)   ! ! boundary on aux

      call copy  (fcd,fdf,nf)         ! copy 
c     write(6,*) '4.1.1.2'
c     call outfacez1(fcd,'fcd       ',1)  
      call add2  (fcd,btf,nf)         ! first add bc , a- + a+(bc) 
c     write(6,*) '4.1.1.3'
c     call outfacez1(btf,'btf       ',1)  
c     call outfacez1(fcd,'fcd       ',1)  
      call gs_op (dg_hndl,fcd,1,1,0)  ! 1 ==> +      , a- + a+(itr)
c     write(6,*) '4.1.1.4'
c     call outfacez1(fcd,'fcd       ',1)  
      call cmult (fcd,0.5,nf)         ! times 1/2, (a- + a+)/2 
c     write(6,*) '4.1.1.5'
c     call outfacez1(fcd,'fcd       ',1)  

      call copy  (gcd,gdf,nf)         ! copy 
      call add2  (gcd,btg,nf)         ! first add bc , a- + a+(bc) 
      call gs_op (dg_hndl,gcd,1,1,0)  ! 1 ==> +      , a- + a+(itr)
      call cmult (gcd,0.5,nf)         ! times 1/2, (a- + a+)/2 

      if(if3d) then 
        call copy  (hcd,hdf,nf)         ! copy 
        call add2  (hcd,bth,nf)         ! first add bc , a- + a+(bc) 
        call gs_op (dg_hndl,hcd,1,1,0)  ! 1 ==> +      , a- + a+(itr)
        call cmult (hcd,0.5,nf)         ! times 1/2, (a- + a+)/2 
      endif

c ! ! ! 
      call copy    (ctrh(1,1), fcd, nf) 
      call copy    (ctrh(1,2), gcd, nf) 
      if(if3d) call copy    (ctrh(1,3), hcd, nf) 

      return
      end
c-----------------------------------------------------------------------
      subroutine htbc_tmp(bf1,bg1,bh1) ! 3 + 2
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      include 'DG'      ! dg_face is stored
      common /nekcb/ cb
      character cb*3
      integer  i, f, e, k
      real     bf1(1), bg1(1), bh1(1) 

      n  = nx1*ny1*nz1*nelt
      nf = nx1*nz1*2*ndim*nelt 
      nxz= nx1*nz1
      nfaces=2*ndim
c
      ifield=2 ! so that it reads temp bc 
c
      call rzero(bf1,nf)
      call rzero(bg1,nf)
      call rzero(bh1,nf)
c     
      k=0
      do e=1,nelt
      do f=1,nfaces
         ieg=lglel(e)
         cb =cbc(f,e,ifield)
c        write(6,*)'ifield', ifield 
c        write(6,*)'r u freaking kidding me, cb = ', cb
         if(cb.eq.'f  ' .or. cb.eq.'F  ' .or.  ! Flux 
     $      cb.eq.'t  ' .or. cb.eq.'T  '       ! Temp 
     $      ) then
            ia=0
            call facind(kx1,kx2,ky1,ky2,kz1,kz2,nx1,ny1,nz1,f)
            do iz=kz1,kz2
            do iy=ky1,ky2
            do ix=kx1,kx2
               ia = ia + 1
               k  = ia + (f-1)*nxz + (e-1)*nfaces*nxz

               if(cb.eq.'t  ' .or. cb.eq.'T  ') then ! 
                 call nekasgn(ix,iy,iz,e)
c                call userbc_f(ix,iy,iz,f,e,ieg) 
                 bf1(k) = fdf(k) 
                 bg1(k) = gdf(k)  ! BC for F_v
                 if(if3d) bh1(k) = hdf(k) 
               else if(cb.eq.'f  '.or. cb.eq.'F  ') then
c if no boundary conditions on \nabla u \cdot n 
c set bf1 = fdf 
                 call nekasgn(ix,iy,iz,e)
                 call userbc_f(ix,iy,iz,f,e,ieg) 
                 bf1(k) = - fdf(k) 
                 bg1(k) = - gdf(k)  ! BC for F_v
                 if(if3d) bh1(k) = - hdf(k) 
               endif
            enddo 
            enddo 
            enddo 
         endif
      enddo 
      enddo 

      return
      end
c-----------------------------------------------------------------------
      subroutine ht_vol_res(rhs) ! ! 3+2? 
c
c  weak form: 1. times mass matrix 
c             2. transpose of derivatives 
c             3. times mass matrix inverse 
c
      include 'SIZE'
      include 'TOTAL'
      include 'DGUSE'
      parameter(le=lx1*ly1*lz1)
      integer  e
      real     rhs(le,1), tm2(le,lelt)

      n=nx1*ny1*nz1*nelt
      ne=nx1*ny1*nz1
      
      if(if3d) then
        call rzero(tm2,n)
        write(6,*) 'ht_vol_res :: 3d not ther e' 
        call exitt 
      else  ! 
        call rzero      (tm2,n)
        call rzero      (hd,n)
        call grad_bm1_t (tm2,sht(1,1),sht(1,2),hd) ! no de-aliasing now 

c       write(6,*) '3.1'   
c       call outfldz1(tm2,'tm2      3',1)  

      endif

      call invbf        (rhs,tm2) ! all elements 
c     write(6,*) '3.2'   
c     call outfldz1(bm1,  'bm1       ',1)  
c     call outfldz1(rhs,'rhs      3',1)  

      call chsign       (rhs,n)   ! all elements 
c     write(6,*) '3.3, change sign'   
c     call outfldz1(rhs,'rhs      3',1)  
c 
      return
      end
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
c+---------------------------------------------------------------------+
c|                    +-----+    +-    +-+    +--_                     |
c|                    | |---|    | \   | |    |   \                    |
c|                    | +-|      | |\  | |    | |  |                   |
c|                    | +-|      | | \ | |    | |  |                   |
c|                    | |---|    | |  \| |    |   /                    |
c|                    +-----+    +-+    -+    +---                     |
c+---------------------------------------------------------------------+
