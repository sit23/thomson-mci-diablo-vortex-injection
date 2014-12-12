c********************************************************************************|
c
c qgbt_diablo.f --> Barotropic QG model based on Diablo DNS solver.
c
c This is a Fortran 77 code that computes barotropic QG flow.
c
C
C The original code was developed as a joint project in MAE 223 (CFD), taught by
C T. Bewley, at UC San Diego (spring 2001), and was provided to Stephen Thomson by 
C Andy Thompson and Emma Boland in October 2011.
C The present version of the code encorperates a vortex injection scheme developed by
C Stephen Thomson and Michael McIntyre, and further information can be found on their websites:
C http://www.damtp.cam.ac.uk/user/sit23/ 
C http://www.atm.damtp.cam.ac.uk/mcintyre/
C Any enquires should be sent to stephen.i.thomson@gmail.com
C Formal documentation discussing the vortex injection scheme will be available in due course,
C although details are discussed in our submitted paper, available on our websites.
C 
C******************************************************************************|
C
C This code is free software; you can redistribute it and/or modify it
C under the terms of the GNU General Public License as published by the
C Free Software Foundation; either version 2 of the License, or (at your
C option) any later version. This code is distributed in the hope that it
C will be useful, but WITHOUT ANY WARRANTY; without even the implied
C warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
C GNU General Public License for more details. You should have received a
C copy of the GNU General Public License along with this code; if not,
C write to the Free Software Foundation, Inc., 59 Temple Place - Suite
C 330, Boston, MA 02111-1307, USA.
C
C******************************************************************************|

c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      program qgbt_diablo
      include 'qgbt_diablo.h'

      write(6,*)
      write(6,*) ' ***** Welcome to QG Diablo (BT version) *****'
      write(6,*)

      call initialize
      
      do time_step = time_step+1, time_step + n_time_steps
         

         call AB_PER
c         write(6,*) 'time = ', time

         time = time+delta_t
         if (mod(time_step,save_stats_int).eq.0) then
            call save_stats
         end if
         if (mod(time_step,save_prof_int).eq.0) then
            call save_prof
         end if
         if (mod(time_step,save_phys_int).eq.0) then
            call save_phys
         end if
         if (mod(time_step,save_flow_int).eq.0) then
            call save_flow(.false.)
         end if
      
      end do
      
      write(6,*)
      write(6,*) '* Your simulation has now come to a close,' 
      write(6,*) '  thank you for flying QG Diablo *'
      write(6,*)

      end

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine initialize
      include 'qgbt_diablo.h'

      real    version, current_version
      logical create_new_flow, reset_time
      integer i, j, k

      namelist /inpt/ version,beta,kappa,nu,hx,hy,lx,ly,
     &                csx,csy,n_time_steps,delta_t,reset_time,
     &                num_per_dir,time_ad_meth,verbosity,
     &                save_flow_int,save_phys_int,save_prof_int,
     &                save_stats_int, create_new_flow,ic_type,
     &                energy,ic_1,ic_2,bc1_type,bc2_type,hamp,
     &                kcut,filter_exp,ldsqd,
     &		    kt,sr,cmc,a,b,to,xo,psilim,tcool,
     &		    bias_factor,time_bet_storms,
     &	          q_thresh,use_new_storms


c Read input file.
c   (Note--if you change the following section of code, update the 
c    current_version number to make obsolete previous input files!)

      write(6,*) 'We are now reading inputs and initializing the flow.'
      
      current_version = 1.0
      read(5,inpt)
      print inpt

c Check version

      if (version .ne. current_version) stop 'Wrong input format.'

c Initialize grid
      write(6,*) 'Grid size: nx =',nx,', ny =',ny,'.'
      nxm=nx-1
      nym=ny-1

      if (num_per_dir .gt. 0) then
         write(6,*) 'Fourier in X'
         do i=0,nx
            gx(i)=(i*lx)/nx
            if (verbosity .gt. 3) write(6,*) 'gx(',i,') = ',gx(i)
         end do
      else
         write(6,*) 'Finite-difference in X'
         do i=0,nx
            gxf(i)=(lx/2.0)*tanh(csx*((2.0*i)/nx-1.0))/tanh(csx)
            if (verbosity .gt. 3) write(6,*) 'gxf(',i,') = ',gxf(i)
         end do
         do i=1,nx
            gx(i)=(gxf(i)+gxf(i-1))/2.0
            if (verbosity .gt. 3) write(6,*) 'gx(',i,') = ',gx(i)
            dx(i)=(gxf(i)-gxf(i-1))
         end do
         gx(0)=gxf(0)
         gx(nx+1)=gxf(nx)
         do i=0,nx
            dxf(i)=(gx(i+1)-gx(i))
         end do
      end if

      if (num_per_dir .gt. 1) then
         write(6,*) 'Fourier in Y'
         do j=0,ny
            gy(j)=(j*ly)/ny
            if (verbosity .gt. 3) write(6,*) 'gy(',j,') =',gy(j)
         end do
      else
         write(6,*) 'Finite-difference in Y'
         do j=0,ny
            gyf(j)=(ly/2.0)*tanh(csy*((2.0*j)/ny-1.0))/tanh(csy)
            if (verbosity .gt. 3) write(6,*) 'gyf(',j,') = ',gyf(j)
         end do
         do j=1,ny
            gy(j)=(gyf(j)+gyf(j-1))/2.0
            if (verbosity .gt. 3) write(6,*) 'gy(',j,') = ',gy(j)
            dy(j)=(gyf(j)-gyf(j-1))
         end do
         gy(0)=gyf(0)
         gy(ny+1)=gyf(ny)
         do j=0,ny
            dyf(j)=(gy(j+1)-gy(j))
         end do
      end if

c Initialize store arrays
      
      do k=0,2
         do j=0,ny+1
            do i=0,nx+1
               UG(i,j,k)=0.
               VG(i,j,k)=0.
               QX(i,j,k)=0.
               QY(i,j,k)=0.
            end do
         end do
      end do

c Initialize FFT package (includes defining the wavenumber vectors).

      call init_fft
      if (verbosity .eq. 100) call test_fft

c Initialize RKW3 parameters.

c      alpha(1)=8.0/15.0
c      alpha(2)=2.0/15.0
c      alpha(3)=5.0/15.0
c      beta(1)=8.0/15.0
c      beta(2)=5.0/12.0
c      beta(3)=3.0/4.0
c      zeta(1)=0.0
c      zeta(2)=-17.0/60.0
c      zeta(3)=-5.0/12.0

c Initialize Adams-Bashforth parameters.
      
c      abp(1)=23.0/12.0
c      abp(2)=-16.0/12.0
c      abp(3)=5.0/12.0
      new = 2
      now = 1
      old = 0

c Initialize the Roberts filter parameter.

      rob = 0.005

c Initialize case-specific packages.

      if (num_per_dir.eq.2) then
         call init_per
      elseif (num_per_dir.eq.1) then
c         call init_chan
      elseif (num_per_dir.eq.0) then
c         call init_basin
      end if

c Initialize flow.

      if (reset_time .or. create_new_flow) then
         previous_time_step = 0
         time_step = 0
         time = 0
         write(6,*) 'The time has been reset.'
      end if

      if (create_new_flow) then
         if (num_per_dir.eq.2) then
            call create_flow_per
         elseif (num_per_dir.eq.1) then
c            call create_flow_chan
         elseif (num_per_dir.eq.0) then
c            call create_flow_basin
         end if
         call save_stats
         call save_phys
c         call save_flow(.false.)
      else
         if (num_per_dir.eq.2) then
            call read_flow_per
         elseif (num_per_dir.eq.1) then
c            call read_flow_chan
         elseif (num_per_dir.eq.0) then
c            call read_flow_basin
         end if
      end if

c Initialize filter

      call init_filter

c Initialize the topography

      call init_topo

c initialise some fields

	energy_new=0.
	energy_now=0.
	energy_old=0.

      return
      end

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine read_flow
      include 'qgbt_diablo.h'

      integer i, j, nx_t, ny_t, num_per_dir_t
      
c      fname = 'qg_diablo.res'
      write(6,*) 'Reading streamfunction from fort.40'

c      open(unit=10, file=fname, status="old", form="unformatted")
      read(40,*) nx_t, ny_t, num_per_dir_t, time, time_step,psibar2(new)

      if ((nx .ne. nx_t) .or. (ny .ne. ny_t)) 
     *   STOP 'Old flowfield wrong dimensions'
      if (num_per_dir .ne. num_per_dir_t) 
     *   STOP 'Old flowfield wrong per directions'

      if (num_per_dir.eq.2) then
         read(40,*) ((CQ(i,j,new),i=0,nkx),j=0,tnky)
      elseif (num_per_dir.eq.1) then
         read(40,*) ((CQ(i,j,new),i=0,nkx),j=1,ny)
      elseif (num_per_dir.eq.0) then
         read(40,*) ((CQ(i,j,new),i=1,nx),j=1,ny)
      end if
!       psibar2(new)=0.0
      return
      end

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine save_flow(final)
      include 'qgbt_diablo.h'
      
      integer i,j
      logical final
      
      write(6,*) 'Writing streamfunctions to fort.30'
      write(6,*) 'This is record ', time_step/save_flow_int,
     *           '.  The time is', time,'.'

      write(30,*) nx, ny, num_per_dir, time, time_step,psibar2(new)
      write(30,*) ((CQ(i,j,new),i=0,nkx),j=0,tnky)

      return
      end


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine save_phys
      include 'qgbt_diablo.h'

      real*8 W(0:nx+1,0:ny+1,0:1), W2(0:nx+1,0:ny+1),jjj
      complex*16 CW(0:nx/2,0:ny+1,0:1)
      equivalence(W,CW)
      integer i,j,jj,ii

      jjj=0.0
      write(6,*) 'Writing the physical fields to file at time = ',time


c ********************************************************
c     Adding this bit to calculate u and v profiles

      do j=0,tnky 
         do i=0,nkx
            CUG2(i,j,new) =-ciky(j)*CP(i,j,new)  ! u
            CVG2(i,j,new) = cikx(i)*CP(i,j,new)  ! v
            CUYY(i,j,new) = ciky(j)*ksqd(0,j)*CP(i,j,new)
            CUSQD(i,j,new) =-ciky(j)*CP(i,j,new)*ciky(j)*
     *      CONJG(CP(i,j,new)) ! usqd
            CVSQD(i,j,new) = -cikx(i)*CP(i,j,new)*cikx(i)*
     *      CONJG(CP(i,j,new))  ! vsqd
            CPSQD(i,j,new) = CP(i,j,new)*CONJG(CP(i,j,new))  ! psisqd

         end do
      end do
      cpsqd(0,0,new) = CMPLX(psibar2(new)*psibar2(new),0.0)

      call fft_xy_to_physical(CUG2,UG2,new)
      call fft_xy_to_physical(CVG2,VG2,new)
      call fft_xy_to_physical(CUYY,UYY,new)



c **********************************************************

      do j=0,tnky
         do i=0,nkx
            CW(i,j,0) = CP(i,j,new)
            CW(i,j,1) = CQ(i,j,new)-(1./ldsqd)*CP(i,j,new)
         end do
      end do

      call fft_xy_to_physical(CW,W,0)
      call fft_xy_to_physical(CW,W,1)

      do jj=0,nym	
		if(jj.gt.0) then 
		jjj=jjj+1.0
		end if
			do ii=0,nxm
      		W2(ii,jj)=W(ii,jj,1)+beta*(ly/NY)*jjj+H(ii,jj,0)
			end do
      end do

      write(12,*) ((W(i,j,0),i=0,nxm),j=0,nym) !Streamfunction
      write(14,*) ((W(i,j,1),i=0,nxm),j=0,nym) !PV without beta*y+psi_2/ldsqd
      write(16,*) ((W2(i,j),i=0,nxm),j=0,nym) !PV with beta*y+psi_2/ldsqd
      write(18,*) ((UG2(i,j,new),i=0,nxm),j=0,nym) !g/s u
      write(20,*) ((VG2(i,j,new),i=0,nxm),j=0,nym) !g/s v
      write(22,*) ((UYY(i,j,new),i=0,nxm),j=0,nym)
      write(25,*) ((CUSQD(i,j,new),i=0,nkx),j=0,tnky)
      write(26,*) ((CVSQD(i,j,new),i=0,nkx),j=0,tnky)
      write(27,*) ((CPSQD(i,j,new),i=0,nkx),j=0,tnky)
      write(29,*) ((CP(i,j,new),i=0,nkx),j=0,tnky) !psi FFT
      write(19,*) (upvpbar(j),j=0,nym)
      write(89,*) e_rate_average(0),e_rate_average_u(0),
     *		e_rate_average_v(0),e_rate_average_psi(0),
     *		e_rate_average(1),e_rate_average_u(1),
     *		e_rate_average_v(1),e_rate_average_psi(1),
     *		e_rate_average(2),e_rate_average_u(2),
     *		e_rate_average_v(2),e_rate_average_psi(2)
      flush(19)
      flush(89)
      
      do j=0,nym
      upvpbar(j)=0. !Resets upvpbar sum. Need to divide outputted upvpbar/save_phys_int, as then it is average over each window of that length.
      end do
	e_rate_average(0)=0.
	e_rate_average_u(0)=0.
	e_rate_average_v(0)=0.
	e_rate_average_psi(0)=0.
	e_rate_average(1)=0.
	e_rate_average_u(1)=0.
	e_rate_average_v(1)=0.
	e_rate_average_psi(1)=0.
	e_rate_average(2)=0.
	e_rate_average_u(2)=0.
	e_rate_average_v(2)=0.
	e_rate_average_psi(2)=0.


!      write(28,*) (cikx(i),i=0,nkx)
!      write(28,*) (ciky(i),i=0,tnky)
!      write(28,*) nkx,tnky


      return
      end

C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine save_prof
      include 'qgbt_diablo.h'

      real*8 W(0:nx+1,0:ny+1,0:1)
      complex*16 CW(0:nx/2,0:ny+1,0:1)
      equivalence(W,CW)
      integer i,j

      write(6,*) 'Writing the zonal profiles to file at time = ',time

      CW = 0.
      do j=0,tnky
         CW(0,j,0) = CP(0,j,new)
         CW(0,j,1) = CQ(0,j,new)
      end do
      
      call fft_xy_to_physical(CW,W,0)
      call fft_xy_to_physical(CW,W,1)
      write(31,*) (W(0,j,0),j=0,nym)
      write(33,*) (W(0,j,1),j=0,nym)

      CW = 0.
      do j=0,tnky
         do i=0,nkx
            CW(i,j,0) = cikx(i)*CP(i,j,new)
         end do
      end do

      call fft_xy_to_physical(CW,W,0)
      do j=0,nym
         do i=0,nxm
            W(i,j,0) = W(i,j,0)*sin(gx(i))
         end do
      end do
      call fft_xy_to_fourier(W,CW,0)
      CW(:,:,1) = 0.
      do j=0,tnky
         CW(0,j,1) = CW(0,j,0)
      end do
      call fft_xy_to_physical(CW,W,1)
      write(35,*) (W(0,j,1),j=0,nym)

      CW = 0.
      do j=0,tnky
         do i=0,nkx
            CW(i,j,0) = cikx(i)*CP(i,j,new)
            CW(i,j,1) = CQ(i,j,new)
         end do
      end do

      call fft_xy_to_physical(CW,W,0)
      call fft_xy_to_physical(CW,W,1)
      do j=0,nym
         do i=0,nxm
            W(i,j,0) = W(i,j,0)*W(i,j,1)
         end do
      end do
      call fft_xy_to_fourier(W,CW,0)
      CW(:,:,1) = 0.
      do j=0,tnky
         CW(0,j,1) = CW(0,j,0)
      end do
      call fft_xy_to_physical(CW,W,1)
      write(36,*) (W(0,j,1),j=0,nym)

      CW = 0.
      do j=0,tnky
         do i=0,nkx
            CW(i,j,0) =  cikx(i)*CP(i,j,new)
            CW(i,j,1) = -ciky(j)*CP(i,j,new)
         end do
      end do

      call fft_xy_to_physical(CW,W,0)
      call fft_xy_to_physical(CW,W,1)
      do j=0,nym
         do i=0,nxm
            W(i,j,0) = W(i,j,0)*W(i,j,1)
         end do
      end do
      call fft_xy_to_fourier(W,CW,0)
      CW(:,:,1) = 0.
      do j=0,tnky
         CW(0,j,1) = CW(0,j,0)
      end do
      call fft_xy_to_physical(CW,W,1)
      write(38,*) (W(0,j,1),j=0,nym)


      return
      end


C----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine save_stats
      include 'qgbt_diablo.h'
      
      real psisq, gpsisq, nabpsisq, nabgpsisq, nabgpsibar, 
     *     psiysq, psixsq, psi2ysq, psi2xsq, psibar, gpsibar, 
     *     nabpsibar, psixf,psisq_2,psiysq_2,psixsq_2,
     *     psiybarsq,psixbarsq
     
      real temp1
      integer i,j

      write(6,*) 'Saving flow statistics at time t = ',time

      psisq = 0.
      psisq_2= 0.
      gpsisq = 0.
      nabpsisq = 0.
      nabgpsisq = 0.
      psiysq = 0.
      psixsq = 0.
      psixf = 0.

      psibar = 0.
      gpsibar = 0.
      nabpsibar = 0.
      nabgpsibar = 0.
      psiybarsq = 0.
      psixbarsq = 0.

      do j=0,tnky
         do i=1,nkx
            temp1 = ksqd(i,j)
            psisq = psisq + dble( CP(i,j,new)
     *             * conjg(CP(i,j,new)) )
            gpsisq = gpsisq + temp1*dble( CP(i,j,new)
     *              * conjg(CP(i,j,new)) )
            nabpsisq = nabpsisq + temp1*temp1
     *                * dble( CP(i,j,new)*conjg(CP(i,j,new)) )
            nabgpsisq = nabgpsisq + temp1*temp1*temp1
     *                 * dble( CP(i,j,new)*conjg(CP(i,j,new)) )
            psiysq = psiysq + ky2(j)*dble( CP(i,j,new) 
     *              * conjg(CP(i,j,new)) )
            psixsq = psixsq + kx2(i)*dble( CP(i,j,new)
     *              * conjg(CP(i,j,new)) )
         end do
      end do

      do j=0,nky
         i=0
            temp1 = ksqd(i,j)
            psisq = psisq + 0.25*dble( CP(i,j,new)
     *             * conjg(CP(i,j,new)) )
            gpsisq = gpsisq + 0.25*temp1*dble( CP(i,j,new)
     *              * conjg(CP(i,j,new)) )
            nabpsisq = nabpsisq + 0.25*temp1*temp1
     *                * dble( CP(i,j,new)*conjg(CP(i,j,new)) )
            nabgpsisq = nabgpsisq + 0.25*temp1*temp1*temp1
     *                 * dble( CP(i,j,new)*conjg(CP(i,j,new)) )
            psiysq = psiysq + 0.25*ky2(j)*dble( CP(i,j,new) 
     *              * conjg(CP(i,j,new)) )
            psixsq = psixsq + 0.25*kx2(i)*dble( CP(i,j,new)
     *              * conjg(CP(i,j,new)) )
            psixf = psixf + 0.25*imag(cikx(32)*CP(32,j,new))
         
      end do
	
	do j=0,tnky
      do i=0,nkx
         psisq_2 = psisq_2 + dble( CP(i,j,new)
     *             * conjg(CP(i,j,new)) )
         psiysq_2 = psiysq_2 + ky2(j)*dble( CP(i,j,new) 
     *              * conjg(CP(i,j,new)) )
         psixsq_2 = psixsq_2 + kx2(i)*dble( CP(i,j,new)
     *              * conjg(CP(i,j,new)) )
	enddo
	enddo
	

      i=0
      do j=1,nky
         psibar = psibar + dble( CP(i,j,new)*conjg(CP(i,j,new)) )
         gpsibar = gpsibar + ky2(j)*dble( CP(i,j,new)
     *            * conjg(CP(i,j,new)) )
         nabpsibar = nabpsibar + ky2(j)*ky2(j)*dble( CP(i,j,new)
     *              * conjg(CP(i,j,new)) )
         nabgpsibar = nabgpsibar + ky2(j)*ky2(j)*ky2(j)
     *               * dble( CP(i,j,new)*conjg(CP(i,j,new)) )
            psiybarsq = psiybarsq + ky2(j)*dble( CP(i,j,new) 
     *              * conjg(CP(i,j,new)) )
            psixbarsq = psixbarsq + kx2(i)*dble( CP(i,j,new)
     *              * conjg(CP(i,j,new)) )
      end do

c Multiply by two for Fourier sums!
      
      psisq = 2.*psisq
      psisq_2 = 2.*psisq_2
      gpsisq = 2.*gpsisq
      nabpsisq = 2.*nabpsisq
      nabgpsisq = 2.*nabgpsisq
      psiysq = 2.*psiysq
      psixsq = 2.*psixsq
      psixf = 2.*psixf

      psibar = 2.*psibar
      gpsibar = 2.*gpsibar
      nabpsibar = 2.*nabpsibar
      nabgpsibar = 2.*nabgpsibar
      psiybarsq= 2.*psiybarsq
      psixbarsq = 2.*psixbarsq
      

      write(21,*) time, psisq, gpsisq, nabpsisq, nabgpsisq, 
     *            psiysq, psixsq, psixf, psibar, gpsibar, 
     *            nabpsibar, nabgpsibar,psisq_2,psixsq_2,
     *		psiysq_2,psiybarsq,psixbarsq

      write(24,*) psibar2(new)
      flush(21)
      flush(24)


      write(6,*) 'GradPsiSq = ',gpsisq,'Psi_xF = ',psixf
      write(6,*) 'UBarSq = ',gpsibar
      write(6,*) ' '

      return
      end
      
