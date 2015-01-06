c******************************************************************************|
c periodic.f, the fully periodic barotropic QG solvers for diablo.  VERSION 5.1 of Diablo
c Version 1.0 of released vortex injection code.
c
c These solvers were written by AFT (March 2009).
c The vortex injection code was written, and the original code modified by 
c Stephen Thomson (Dec 2014).
c******************************************************************************|

c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      SUBROUTINE AB_PER
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c Main time-stepping algorithm for the fully periodic case
c using a leap-frog time step with a weak Roberts filter.  Also
c an exponential wavenumber filters a la Smith.
c 
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      include 'qgbt_diablo.h'
      integer i, j,initial

      new = mod(new+1,3)
      now = mod(now+1,3)
      old = mod(old+1,3)
      rob = 0.005
	initial=0
c Time step the potential vorticities using leap-frog time step
c with a weak Roberts filter to eliminate the computational mode.

      do j=0,tnky
         do i=0,nkx
         

            CQ(i,j,new) = filter(i,j)*( CQ(i,j,old)
     *                       + 2.*delta_t*CR(i,j,now) )
		
          end do
      end do
      

      psibar2(new)=psibar2(old) + real(2.*delta_t*CR2(now))

	if(injection_flag.ne.1) then !Only apply the Roberts filter if no injection is currently taking place. This is so as to not allow the filter to EXCITE a computational mode during injections. Note that the injection flag is SET during the calculation of CR. The flag is therefore set during the calculation of CR which has just been used above in the timestepping. Therefore it is correct to not apply the filter when forcing from an injection is being used in CR, which is now, if this injection_flag is 1.
      do j=0,tnky
         do i=0,nkx
            CQ(i,j,now) = (1.-2.*rob)*CQ(i,j,now) +
     *                     rob*(CQ(i,j,new) + CQ(i,j,old))
         end do
      end do
      endif

      call q_to_psi(new)
      call compute_rhs(new,initial)
      call compute_upvpbar(new)
      call compute_energy_rate(new)

      return
      end


c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine compute_rhs(k,initial)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      include 'qgbt_diablo.h'

      integer i, j, k, iii,jjj,juplimit,jlowlimit, iuplimit,ilowlimit,
     *		ii,jj,initial,deltai,p,nstorms,ij,ilowlimitradius,
     *          iuplimitradius,jlowlimitradius,juplimitradius,
     *      iradius,jradius,exactlocation,iuptemp,
     *      ilowtemp,jjk,jjjk,iik,imid,jmid,stormpair,stormtype,
     *      time_step_read,i_read,
     *      j_read,i_read_str,j_read_str,time_step_read_str,
     *	k_read,i_read_fil,j_read_fil,mvcpflag_read,deepflag,
     *	i_fif,j_fif

      real*8 temp1, temp2,tbegin,W(0:nx+1,0:ny+1,0:1),
     *	 W2(0:nx+1,0:ny+1),ranf,signe,loopmeasure,strengthfactor,
     *	rsqd,maxdepth,
     *      W_read_current,W2_read_current,W2_read_new,pairflag_read,
     *      oriflag_read,icentral_read,strength_read,ro_read,
     *	depth_read,t_thr_storm_read, icentralnow_read,
     *	depthnow_read,qp_fil,nextstormi_fil,nextstormj_fil,
     *	qp_new_fil,time_step_fil,pairflag_fil,maxdepth_fil,
     *	time_step_fif,zero_fif,pairflag_fif,oriflag_fif,
     *	depthnow_fif,icentralnow_fif,delta_q,delta_q_acyc,
     *	delta_q_cyc,new_strength,req_delta_q,scale_factor,
     *	strength_acyc,strength_cyc,new_strength_acyc,fudge,
     *	thresh_check,prev_thresh_check,storm_present_flag,
     *	overlap_flag

      complex*16 CW(0:nx/2,0:ny+1,0:1)
      equivalence(W,CW)




c Begin with the non-linear terms and the forcing term

      do j=0,tnky
         do i=0,nkx
            CUG(i,j,k) =-ciky(j)*CP(i,j,k)  ! u
            CVG(i,j,k) = cikx(i)*CP(i,j,k)  ! v
            CQX(i,j,k) = cikx(i)*CQ(i,j,k)  ! q_x
            CQY(i,j,k) = ciky(j)*CQ(i,j,k)  ! q_y
            CW(i,j,0) = CP(i,j,k) !Psi is required at this stage to allow the PV field to be calculated, in order to apply condition (3.8) later.
         end do
      end do
      CW(0,0,0)=0.0
      call fft_xy_to_physical(CUG,UG,k)
      call fft_xy_to_physical(CVG,VG,k)
      call fft_xy_to_physical(CQX,QX,k)
      call fft_xy_to_physical(CQY,QY,k)
      call fft_xy_to_physical(CW,W,0)


c **************************  **********************************
c**********  *************

	    do j=0,tnky
		  do i=1,nkx ! Want all of the field except the zonal mean
		  CW(i,j,1) = CQ(i,j,k)-(1./ldsqd)*CP(i,j,k) ! Calculate PV for later on.
		  end do
		  CW(0,j,1)=0.0 !Don't forget to set previously un-set parameters to zero!
	    end do

	    call fft_xy_to_physical(CW,W,1)

c ********************************************************************************************

 	if(initial.eq.1) nstorms=0.0
 
 	if((initial.ge.1).and.(k.eq.1)) then

	!Storm initialisation

	if(use_new_storms) then
        oriflag=0
       nextstormi=128 !This defines where the first storm is.
      nextstormj=128
      nextstorm=time_step+3
      pairflag=1 !We then have to say that the partner to this first storm must be on the right, otherwise it won't work! :) 
      storm_id=1.
 	write(6,*) nstorms
	else

! If we are going to be using old storms from a previous simulation, we have to wind files on to point we're restarting from.

			time_step_read=0

                  if(time_step_read.ne.time_step) then

                  do while(time_step_read.lt.time_step)
                  read(74,*) time_step_read,i_read,j_read,
     *            nextstorm,nextstormi,nextstormj,
     *            pairflag_read,W_read_current,W2_read_current,
     *            W2_read_new,oriflag_read,icentral_read !The problem here is that we want nextstorm and nextstormi and nextstormj to come from the record BEFORE the one that this ticks over to. So, for the initial step, let's take nextstorm = time_step_read etc...
                  end do
		backspace(74)
		nextstorm=time_step_read
		nextstormi=i_read
		nextstormj=j_read
		oriflag=oriflag_read
		icentralnow=icentral_read
		pairflag=pairflag_read !Have to initialise the pairflag, otherwise it doesn't work!

		write(6,*) time_step_read,time_step,i_read,i,j_read,j
     *	,'initialising storms'


! New initialising subsection, so that we can read the injection strengths from the fort.90 file correctly.

		time_step_read_str=0
		
     		do while(time_step_read_str.lt.time_step)

		read(75,*) i_read_str,j_read_str,strength_read,ro_read,
     *	time_step_read_str,k_read
     *		,depth_read,t_thr_storm_read,
     *            icentralnow_read,depthnow_read
		end do
		backspace(75)
		
		time_step_fil=0
		
     		do while(time_step_fil.lt.time_step)
     		
	    read(76,*) i_read_fil,j_read_fil,qp_fil,
     *	nextstormi_fil,nextstormj_fil,
     *      qp_new_fil,mvcpflag_read,
     *      time_step_fil,pairflag_fil,maxdepth_fil
		end do
		backspace(76)

		
		time_step_fif=0
		
		do while(time_step_fif.lt.time_step)

            read(77,*) time_step_fif,i_fif,j_fif,zero_fif,
     *	pairflag_fif,oriflag_fif,
     *      depthnow_fif,icentralnow_fif

		end do
		backspace(77)


                  end if

	end if

	end if


!******************************************************************************************************************
!********************************End storm initilisation, and begin actual calulation of u.grad q. ****************
!******************************************************************************************************************

      do j=0,nym !Beware - this do loop is about 200 lines before you actually find the 'end do'. I will set it out better in a future version!
         do i=0,nxm
		!This is where we do the u.grad (q_rel + psi_2/ldsqd) part. H is the topography = psi_deep/ldsqd. As explained briefly in the README.txt file, this setup is equivalent to doing u.grad (q_total), but only once we re-scale the wavenumbers when CR is calulated later on. Note that the variables we need (psi for injection code in this do loop) has already been caluclated, the below is to calculated the tendancy (CR) for the next timestep. Beta * v will also be added later.
             QX(i,j,k) = -UG(i,j,k)*(QX(i,j,k)+H(i,j,1))
	     QX(i,j,k) = QX(i,j,k)-VG(i,j,k)*(QY(i,j,k)+H(i,j,2))

c *****************************************************************
c !********************Beginning of the storm injection code********

		injection_flag=0
	if((abs(location2(i,j,0)).eq.0.0) !Only enters this part of the code if i and j are equal to nextstormi and nextstormj, and time_step = nextstorm, i.e. that the time and place for the next storm is now! :) We will very quickly choose a new 'future storm' pair. This 'end if' lasts for a long time - down to end of the do loops above.
     *	.and.(initial.ne.1).and.(time_step.eq.nextstorm)
     *  .and.(i.eq.nextstormi).and.(j.eq.nextstormj)) then
                
                
          if(use_new_storms) then
		
		call storm_location_chooser(i,j) !Runs the storm location chooser. It is simply a subroutine that uses ALL GLOBAL VARIABLES, apart from i and j, which it is fed. The code works by always having a current injection pair, and a future injection pair. This subroutine chooses the future injection pair, as i and j is the location for the current storm, part of the existing pair. Note that each storm is dealt with individually. So if the above rotuine chooses to start a new pair, then the next time it is run, the routine will make the 2nd storm in the pair. Once this second storm in the pair is the current storm (at i,j) then the routine will choose to start a new pair, etc. 


            else ! If we're using old storms, then do this.
            read(74,*) time_step_read,i_read,j_read,
     *      nextstorm,nextstormi,nextstormj,
     *      pairflag_read,W_read_current,W2_read_current,
     *      W2_read_new,oriflag_read,icentral_read
		icentralnow=icentral_read
            oriflag=oriflag_read
		
            end if

            do while((abs(location2(nint(nextstormi),nint(nextstormj), !Checks if the suggested future storm is already the current storm (could happen by chance), and if this is the case, then chooses a new location for the new storm.
     *      0)).eq.1.0).or.((i.eq.nextstormi).and.(j.eq.nextstormj)))

            write(58,*) time_step,nextstormi,nextstormj,
     *      W2(nint(nextstormi),nint(nextstormj)),
     *      abs(location2(nint(nextstormi),nint(nextstormj),
     *      0)),i,j

            flush(58)
            if(use_new_storms) then !Another extra bit! This will only ever occur in the second half of a pair, and so the next storm we're choosing is already the correct pairflag to continue.

                  nextstormi=nint(NX*rand())
                  nextstormj=nint(NY*rand())
                  nextstormi=mod(nextstormi+REAL(NX),REAL(NX))
			nextstormj=mod(nextstormj+REAL(NY),REAL(NY))
			storm_edgeflag=0.

! BIT COPIED FROM STORM_LOCATION_CHOOSER, as the above bit of code chooses a new left hand storm, but we have to be sure it'll work! This will be explained more in the comments in the 'storm_location_chooser' subroutine, which is near the bottom!
			if((nextstormi-nint((2.*xo*NX/LX)+1.)).lt.0) then
			storm_edgeflag=1.
                  nextstormi=nextstormi!+nint((2.*xo*NX/LX)+1.)
                  icentral=nint(nextstormi+nint((xo*NX/LX)+1.))
                  icentral=nint(icentral-(NX-1.)*
     *            (floor(icentral/(NX-1.))))


			else

                  nextstormi=nextstormi-nint((2.*xo*NX/LX)+1.) !Go for the Left hand storm.
                  nextstormi=mod(nextstormi+REAL(NX),REAL(NX))
                  icentral=nint(nextstormi+nint((xo*NX/LX)+1.))
                  icentral=nint(icentral-(NX-1.)*
     *            (floor(icentral/(NX-1.))))
			end if
! END COPIED BIT!!!!******************************************************************************************************************
            else ! If we're not using new storms.
            read(74,*) time_step_read,i_read,j_read,
     *      nextstorm,nextstormi,nextstormj,
     *      pairflag_read,W_read_current,W2_read_current,
     *      W2_read_new,oriflag_read,icentral_read
            oriflag=oriflag_read
            icentralnow=icentral_read

            end if
            end do


            write(73,*) time_step,i,j,nextstorm,nextstormi,nextstormj, !So this is one of the main outputs for storms that are PROPOSED, BEFORE the condition 3.8 is checked and applied, which is done in 'storm_strength_chooser'.
     *      pairflag,W(i,j,0),W2(i,j),
     *      W2(nint(nextstormi),nint(nextstormj)),oriflag,icentral
		flush(73)
                signe=-signflag
           	location2(i,j,0)=signe !This means, confusingly, that signe is +1 for a/cycs and -1 for cycs.

		if(pairflag.eq.1.) then
		location2(i,j,3)=storm_id
		else
		location2(i,j,3)=storm_id
		storm_id=storm_id+1.!Gives each injected pair a unique number.
		end if
		

!********************************************************************************
!***************SET STRENGTH BASED ON PV BIASED MODEL****************************
!********************************************************************************

!This works by checking both storms in the pair at once. So it will run the check and either allow both storms, or cancel both.

      if(pairflag.eq.1) then !Only runs the check once per pair. 
            depthnow=W(icentralnow,j,0) !Psi value at the point in the centre of the proposed storm pair.

	   mvcpflag=0 !MVCP is the code name for condition 3.8, and stands for the 'monster vortex censorship program'. :) If mvcpflag=0, it doesn't need to be applied. This is the initialisation step.
	   maxdepth=0
	   if(use_new_storms) then !Because of drift, cannot trust filters to work properly after a restart, so we only ask the model to calculate the new strength if we're using new storms. If we're using old storms, we'll simply use the strength it ended up with last time.

		call storm_strength_chooser(i,j,
     *	new_strength_acyc,W) ! Subroutine that applies condition 3.8. (i,j) is first storm in pair, and nextstormi,nextstormj is second.
	
		else !If we're using old storms
! !     What happens when we don't explicely use the filter, is we have to read from the old fort.79 file.

	    	read(76,*) i_read_fil,j_read_fil,qp_fil,
     *	nextstormi_fil,nextstormj_fil,
     *      qp_new_fil,mvcpflag_read,
     *      time_step_fil,pairflag_fil,maxdepth_fil 

		if((i_read_fil.eq.i).and.(j_read_fil.eq.j).and.
     * 	(time_step_fil.eq.time_step)) then
	 		mvcpflag=mvcpflag_read
		else
			backspace(76)
			mvcpflag=0
		end if
	
         end if
	end if


!What the above has done is pick the new strengths for the storm and it's partner if 3.8 needed applying. If it did, then mvcpflag=1. If it didn't, then the strengths have not yet been set, and mvcpflag=0. Now we define the strengths that will be used in either case.



	  if(use_new_storms) then
		
            if((depthnow-(H(icentralnow,j,0))*ldsqd)+
     *	psilim/(2.*bias_factor)
     *      .lt.0.) then !This is checking what zeta is, although the model uses a variant on zeta. Here, if this quantity is less than zero, then zeta>-1/2, and so an injection is allowed. The quantity is really (zeta+1/2)*psilim/bias_factor, where here psilim<0, whereas in the paper it's >0.
		
			if(mvcpflag.eq.0.) then 
			location2(i,j,0)=strengthfactor(location2(i,j,0),depthnow,
     *     		i,j)*location2(i,j,0) !If no 3.8 needed, choose the strengths as normal. 
			else !If 3.8 was needed, apply the actual strengths.
				if(location2(i,j,0).gt.0.) then
				location2(i,j,0)=new_strength_acyc*
     *			location2(i,j,0)
				else
				location2(i,j,0)=new_strength_cyc* !new_strength_cyc is a global variable, so it will have got to this point alright!
     *			location2(i,j,0)
				end if
			end if
     
     
            else !If zeta is less than -1/2, then the injection is forbidden.
			location2(i,j,0)=0.
			write(59,*) time_step,i,j,location2(i,j,0),pairflag,oriflag,
     *	      depthnow,icentralnow !Always gets rid of pairs, which is good! :) 
			flush(59)
            end if
	   else !If we're using old storms, it's a bit more of a faff.

         read(77,*) time_step_fif,i_fif,j_fif,zero_fif,
     *   pairflag_fif,oriflag_fif,
     *   depthnow_fif,icentralnow_fif
     
		if((i_fif.eq.i).and.(j_fif.eq.j).and.
     *	 (time_step_fif.eq.time_step)) then
		deepflag=1
		else
		backspace(77)
		deepflag=0
		end if

            if((deepflag.eq.0).and.(mvcpflag.eq.0.)) then

		do while((i.ne.i_read_str).or.(j.ne.j_read_str))

		read(75,*) i_read_str,j_read_str,strength_read,ro_read,
     *	time_step_read_str,k_read
     *		,depth_read,t_thr_storm_read,
     *            icentralnow_read,depthnow_read
		end do
		
		location2(i,j,0)=strength_read
		else
            location2(i,j,0)=0.
            write(59,*) time_step,i,j,location2(i,j,0),pairflag,oriflag,
     *      depthnow,icentralnow
            flush(59)
            end if
		

	  end if
		


!********************************************************************************
!***************END SET STRENGTH BASED ON PV BIASED MODEL************************
!********************************************************************************

            pairflag=-pairflag
	     location2(i,j,1)=0.0 !This is how far through the storm injection we are (time-wise)
            location2(i,j,2)=xo !Storm radius.

            ilowlimitradius=nint(i-(xo*NX/LX))
            if(ilowlimitradius.lt.0.0) ilowlimitradius=
     *      ilowlimitradius+NX

            iuplimitradius=nint(i+(xo*NX/LX))
            if(iuplimitradius.gt.NX) iuplimitradius=
     *      iuplimitradius-NX

            jlowlimitradius=nint(j-(xo*NX/LX))
            if(jlowlimitradius.lt.0.0) jlowlimitradius=
     *      jlowlimitradius+NY


            juplimitradius=nint(j+(xo*NX/LX))
            if(juplimitradius.gt.NY) juplimitradius=
     *      juplimitradius-NY

	loopmeasure=0.0



		end if !End of 'if' that started at the end do's that finish below.
c		
c *****************************************************************
         end do
      end do !Where the first do loop over spatial variables ends.

	tbegin=0.0 !Now we're actually applying the injection. 
      do j=0,nym
         do i=0,nxm !Loop over space

	   if((abs(location2(i,j,0)).ne.0.0).and.( !Look for storm location
     *location2(i,j,1).ne.tbegin)) then !Note also that it will look for where the time through the injection is non-zero. But, it is set up to be applied for exactly t0/delta_t timesteps.

	   	jlowlimit=nint(j-((location2(i,j,2))*NYM/LY)) !Works out square that encompasses the circular storm.
		if(jlowlimit.lt.(j-((location2(i,j,2))*NYM/LY)))
     *             jlowlimit=jlowlimit+1

		juplimit=nint(j+((location2(i,j,2))*NYM/LY))
		if(juplimit.gt.(j+((location2(i,j,2))*NYM/LY)))
     *             juplimit=juplimit-1

	   	ilowlimit=nint(i-((location2(i,j,2))*NXM/LX))
		if(ilowlimit.lt.(i-((location2(i,j,2))*NXM/LX)))
     *             ilowlimit=ilowlimit+1
		
	   	iuplimit=nint(i+((location2(i,j,2))*NXM/LX))
		if(iuplimit.gt.(i+((location2(i,j,2))*NXM/LX)))
     *             iuplimit=iuplimit-1

		


		do jj=jlowlimit,juplimit !Now runs the loop over the square surrounding the storm, and applies forcing within the circle.
	deltai=nint((location2(i,j,2)*NXM/LX)-((nxm*location2(i,j,2))
     *            /(lx))*(1.-((LY*(jj-j))/
     *            (location2(i,j,2)*NYM))**2.0)**0.5)

                                do ii=ilowlimit+deltai,iuplimit-deltai
				if(ii.lt.0) then
				iii=ii+nx
				elseif(ii.gt.nxm) then
				iii=ii-nx
				else
				iii=ii
				end if

				if(jj.lt.0) then
				jjj=jj+ny
				elseif(jj.gt.nym) then
				jjj=jj-ny
				else
				jjj=jj
				end if


                QX(iii,jjj,k)=QX(iii,jjj,k) !Now we're adding the forcing. 
     *			-location2(i,j,0)*(Sr/cmc)*
     *                  (1.)*
     *			(2.-(2./((location2(i,j,2)*(NXM/LX))**2.))*((ii-i)**2 !Parabolic in space.
     *			+(jj-j)**2))
		if((ii.eq.i).and.(jj.eq.j)) then 
		injection_flag=1 !Set the injection flag, so that the robert filter can be turned off during injection, in order to allow for injection with no computational mode.
		end if
			end do
		end do
		

		if(location2(i,j,1).gt.(to-delta_t).and.
     *		   location2(i,j,1).lt.(to+delta_t)) then
			location2(i,j,0) = 0.0 
			location2(i,j,1) = 0.0 !Designed to cutoff injection after it has been applied when location2(i,j,0) eq t0.
		end if

    	location2(i,j,1)=location2(i,j,1)+delta_t


	end if

	if((abs(location2(i,j,0)).ne.0.0).and.
     *      (location2(i,j,1).eq.tbegin)) then
	location2(i,j,1)=location2(i,j,1)+delta_t !Gives new storm initial time.
	end if


      if((abs(location2(i,j,0)).ne.0.0).and.(initial.le.2)) then
      write(90,*) i,j,location2(i,j,0),location2(i,j,2),time_step,k !Writes full output for each storm after it is injected.
     *		,W(i,j,0),location2(i,j,1),
     *            icentralnow,depthnow,location2(i,j,3),
     *		REALPART(CP(0,1,k)),REALPART(CP(1,0,k)),
     *		REALPART(CP(1,1,k)),
     *		IMAGPART(CP(0,1,k)),IMAGPART(CP(1,0,k)),
     *		IMAGPART(CP(1,1,k))
	flush(90)
	end if

        end do
      end do

 !************************************************************
 !****** Begin calculation of CR *****************************

      call fft_xy_to_fourier(QX,CQX,k) !QX now contains both (u.grad (q_rel +psi_2/ldsqd)) and the storm injections.
      do j=0,tnky
         do i=0,nkx
            CR(i,j,k) = CQX(i,j,k)
         end do
      end do

      CR2(k)=-LDSQD*CR(0,0,k) !This is just the domain average depth, so has no effect anywhere else..


      do j=0,tnky
         do i=0,nkx
	!CR is calulated here. Once this operation is done, we have effectively made it such that we are timestepping PV, and not just relative vorticity.
            CR(i,j,k) = ( CR(i,j,k) - beta*cikx(i)*CP(i,j,k) !note that beta v is included here. 
     *                 -(hy*cikx(i) - hx*ciky(j))*CP(i,j,k) 
     *                 - kappa*CQ(i,j,k) + (((ksqd(i,j))**2.0)*
     *                  CP(i,j,k)*nu)
     *                  )/(1.0+(1.0)/(LDSQD*ksqd(i,j))) !Dividing by this is what means we're effectively timestepping PV. 


         end do
      end do

      CR(0,0,k)=0.0 !Sets domain average depth change to zero, as described in the paper.

      return
      end

c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine psi_to_q(k)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      include 'qgbt_diablo.h'

      integer i, j, k

      do j=0,tnky
         do i=0,nkx
            CQ(i,j,k) = -ksqd(i,j)*CP(i,j,k)  
         end do
      end do

      return
      end


c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine q_to_psi(k)
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      include 'qgbt_diablo.h'
      
      integer i, j, k

      do j=0,tnky
         do i=0,nkx
            CP(i,j,k) = CQ(i,j,k)/(-ksqd(i,j))
         end do
      end do

      CP(0,0,k) = 0.

      return
      end


c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine create_flow_per !Set up intial conditions.
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      include 'qgbt_diablo.h'
      
      real*8 total_en, ekin, aa, alp, phase, ranf
      real*8 temp1, temp2, temp3
      integer i, j, k,initial
      
	initial=1 !Don't add MC events during this time

      write(6,*) 'Creating new flow from scratch.'

      if (ic_type.eq.1) then



      total_en = 0.0

      ekin = 0.0

      CP(0,1,0)=560.!Initialise the streamfunction such that the upper jets are sinusoidal and peak at 35ms-1.
      alp = dble(CP(0,1,0)*conjg(CP(0,1,0)))
      ekin = ekin+ksqd(0,1)*alp

      total_en = ekin

c  Factor of two for Fourier sums

      total_en = total_en*2.
      print*, 'Total Energy = ', total_en


          aa=1.0

      do j=1,tnky
         do i=0,nkx
            CQ(i,j,0) = -ksqd(i,j)*CP(i,j,0) 
         enddo
      enddo
      psibar2(0)=0.0 !Start with the domain averaged layer depth contribution being zero.

      else
         
         write(6,*) 'Beginning a synthetic initial condition'

         CQ(10,0,0) = 0.001*cmplx(1.0,0.)
         
      end if

c Fill in the other initial time steps using Euler

      call compute_rhs(0,initial)
      do j=0,tnky
         do i=0,nkx
            CQ(i,j,1) = CQ(i,j,0)+delta_t*CR(i,j,0)
         end do
      end do

      call q_to_psi(1)
      call compute_rhs(1,initial+2) !The addition to initial here is to stop the 1st time step (when it is done 3 times to set up the leap-frog) doesn't record injections in all 3 initial timesteps.
      do j=0,tnky
         do i=0,nkx
            CQ(i,j,2) = CQ(i,j,1)+delta_t*CR(i,j,1)
         end do
      end do
      call q_to_psi(2)
      call compute_rhs(2,initial+2)!The addition to initial here is to stop the 1st time step (when it is done 3 times to set up the leap-frog) doesn't record injections in all 3 initial timesteps.

      return
      end



c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine read_flow_per
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      include 'qgbt_diablo.h'

      integer i, j, k, nx_t, ny_t, num_per_dir_t, no_rec,initial
      real k1, k2, k3, k4
	initial = 2 ! Don't add moist convection events during this time

      write(6,*) 'Reading PVs from fort.40'

      read(40,*) nx_t, ny_t, num_per_dir_t, time, time_step

c Check that old parameters agree with current parameters

      if ((nx .ne. nx_t) .or. (ny .ne. ny_t))
     *   STOP 'Old flowfield wrong dimensions'
c      if (num_per_dir .ne. num_per_dir_t)
c         STOP 'Old flowfield wrong per directions'
         
      read(40,*) ((CQ(i,j,0),i=0,nkx),j=0,tnky)

      call q_to_psi(0)
      call compute_rhs(0,initial)

c Use simple Euler step to fill in the prev and next fields

      do k=1,2
         do j=0,tnky
            do i=0,nkx
               CQ(i,j,k) = CQ(i,j,k-1) + delta_t*CR(i,j,k-1)
            end do
         end do
         call q_to_psi(k)
         call compute_rhs(k,initial+1)!The addition to initial here is to stop the 1st time step (when it is done 3 times to set up the leap-frog) doesn't record injections in all 3 initial timesteps.
      end do

      return
      end


c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine init_filter
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      include 'qgbt_diablo.h'
      integer i,j
      real*8 kmax, kcut2, kmask, dk, dl, filtdec, temp1

c This is a subroutine to make a filter that both creates
c an isotropic mask to eliminate aliasing as well as introduces
c a exponential cutoff filter to mimic a high order diffusion.
c The scheme is based on Smith et al. (2002).

      write(6,*) 'Initializing the filter.'

c First introduce the aliasing mask.

      kmax = nkx
      kcut2 = kcut*kcut
      kmask = 8./9.*(kmax+1)**2 ! They have done [8/9]*(nkx+1)**2, as the number in the round bracket should be the max number of wavemodes in x, which is nkx+1, as all of the fourier sums in x are done from 0 to nkx. 
      dk = kx(2)-kx(1)

      do j=0,tnky
         do i=0,nkx
            temp1 = ksqd(i,j)/dk/dk ! Note that all of these calculations for the filter are done with the length dimensions removed. This is worrying if you have L_x \ne L_y, as then you'll think, well i^2 and j^2 don't represent the same wavemodes. But, the good news is that because the wavemodes are calculated by ksqd(i,j)/(dk^2), you will end up with temp1=i^2 + (j^2)*(L_x/L_y)^2. So, even though we are dealing with non-dimensional numbers for calculating distances in wavenumber space, it is still alright, as dividing by dkx^2 means we have re-scaled j^2 correctly.
            
            if (temp1.ge.kmask) then
               filter(i,j) = 0.
            else
               filter(i,j) = 1.
            end if
         end do
      end do

c Now introduce the exponential cut-off

      filtdec = -log(1.+2.*pi/kmax)/((kmax-kcut)**filter_exp)
      do j=0,tnky
         do i=0,nkx
            temp1 = ksqd(i,j)/dk/dk
            if (temp1.gt.kcut2) then 
               if (temp1.le.kmask) then
               filter(i,j) = 
     *         exp(filtdec*(sqrt(temp1)-kcut)**filter_exp)
               end if
            end if
         end do
      end do

      return
      end


c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine init_topo
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      include 'qgbt_diablo.h'
      integer i, j


      write(6,*) 'Initializing the topography.'

c Define topography:  H(:,:,0) holds the topography
c                     H(:,:,1) holds the x-derivative of topography
c                     H(:,:,2) holds the y-derivative of topography

      do j=0,nym
         do i=0,nxm
            H(i,j,0) =  (Hamp*cos((j*KT*LY)/ny)-
     *      hamp/(kt**2)*cos((j*KT*LY)/ny)*(ksqd(1,0)))*2.34545 !Gives jets that are 35ms-1, as the upper ones are.
            H(i,j,1) = 0.
            H(i,j,2) =  (-Hamp*KT*sin((j*KT*LY)/ny)+hamp/kt*
     *                  sin((j*KT*LY)/ny)*(ksqd(1,0)))*2.34545
         end do
      end do

c Calculate derivaties and save in physical space.
c dH/dx is held in H(:,:,1) and dH/dy is in H(:,:,2).

c      call fft_xy_to_fourier(H,CH,0)
c      do j=0,tnky
c         do i=0,nkx
c            CH(i,j,1) = cikx(i)*CH(i,j,0)
c         end do
c      end do

c      call fft_xy_to_physical(CH,H,1)

c      do j=0,tnky
c         do i=0,nkx
c            CH(i,j,2) = ciky(j)*CH(i,j,0)
c         end do
c      end do
c      call fft_xy_to_physical(CH,H,2)


      write(6,*) 'Writing topography to fort.50s'
      write(50,*) ((H(i,j,0),i=0,nxm),j=0,nym)
      write(51,*) ((H(i,j,1),i=0,nxm),j=0,nym)
      write(52,*) ((H(i,j,2),i=0,nxm),j=0,nym)

      
      return
      end


c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      subroutine init_per
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      return
      end

c----*|--.---------.---------.---------.---------.---------.---------.-|-------|



!********************************************************************************
!***************calculate STRENGTH BASED ON PV BIASED MODEL**********************
!********************************************************************************
      FUNCTION strengthfactor(pairflagpassed,psibac,ipass,jpass)
      include 'qgbt_diablo.h'
      real*8 psibac,strengthfactor,pairflagpassed,psinew
      integer ipass,jpass

      psinew=psibac-(H(ipass,jpass,0))*
     *      ldsqd+psilim/(2.*bias_factor) !Now calculating zeta field (effectively). It's really (zeta+1/2)*(psilim/bias_factor). But when you end up multiplying it by the appropriate factors, you are really using zeta. As psilim here is negative, the signs are a bit confusing. 
     !Note that psibac is calculated using the central i position between the injections, so BT is the same for each pair. ipass may not be equal to icentral, but as it's only subscripting H, and H is uniform in i, it doesn't matter.

      if(pairflagpassed.eq.1) then !If it's the anticyclone's strength that we need to decide on. 

		if(psinew.gt.(1.*psilim/(bias_factor))) then !(effectively, if zeta is <1/2, note that we've already tested zeta>-1/2 before we enter the subrotuine)
            strengthfactor=(psinew)/psilim !(straight line)
		else
            strengthfactor=(1.*psilim/bias_factor)/psilim !(if zeta>1/2, then we have a constant strength).
		end if

      else !If it's the cyclone strength we need to decide on.

		if(psinew.gt.(1.*psilim/(bias_factor))) then !(if zeta < 1/2)
                  strengthfactor=(psinew*(psilim-psinew))/
     *            (psilim**2.)
		else !(if zeta > 1/2)
                  strengthfactor=((1.*psilim/bias_factor)* 
     *		(psilim-(1.*psilim/bias_factor)))/
     *            (psilim**2.)
		end if

      end if
      return
      end


!********************************************************************************
!***************End calculate STRENGTH BASED ON PV BIASED MODEL******************
!********************************************************************************


! ! !******** upvp calculator**************************
      subroutine compute_upvpbar(k)
! ! 
      include 'qgbt_diablo.h'
      
      real*8 upvptemp(0:nx+1,0:ny+1),up(0:nx+1,0:ny+1),vp(0:nx+1,0:ny+1)
     *      ,upvpbartemp(0:ny+1),ubartemp(0:ny+1),vbartemp(0:ny+1)
     *      ,ubartemp_tick(0:ny+1),utemp,vbartemp_tick(0:ny+1),
     *      vtemp,upvp_tick,upvpbar_old(0:ny+1)
      
      integer k,i,j


      do j=0,nym
         utemp=0.
         vtemp=0.
         do i=0,nxm
         ubartemp_tick(j)=utemp+UG(i,j,k)
         utemp=ubartemp_tick(j)
         vbartemp_tick(j)=vtemp+VG(i,j,k)
         vtemp=vbartemp_tick(j)
         end do
      end do

      ubartemp=ubartemp_tick/nx
      vbartemp=vbartemp_tick/nx


      do j=0,nym
         do i=0,nxm
            up(i,j)=UG(i,j,k)-ubartemp(j)
            vp(i,j)=VG(i,j,k)-vbartemp(j)
            upvptemp(i,j)=up(i,j)*vp(i,j)
         end do
      end do


      do j=0,nym
      upvp_tick=0.
         do i=0,nxm
      upvpbartemp(j)=upvp_tick+upvptemp(i,j)
      upvp_tick=upvpbartemp(j)
         end do
      end do

      upvpbar_old=upvpbar

      upvpbar=upvpbartemp+upvpbar_old !This is essentially just taking an average over all of the timesteps in each 'save_phys_int' window, as that's when upvpbar is reset. Nowhere is it divided by save_phys_int, so have to do that in post-processing.

      return
      end
      
!!!************************************************************
! ! !******** storm location chooser **************************
!!!************************************************************
      subroutine storm_location_chooser(i,j)
! ! 
      include 'qgbt_diablo.h'

	integer i,j
	
	   nextstorm=time_step+nint(0.5*(1.-pairflag))*  !What we now do is always inject the left storm first, and then the RH one at the same timestep, so they always start at the same time.
     *		nint(4.+time_bet_storms*rand())         !Time to next storm that is not part of this pair is = (time_bet_storms*rand() + 4.). Average of rand() =1/2.
            icentralnow=icentral
            if(pairflag.eq.-1) then 	!So we're thinking of the current storm at i,j being the second in the pair (which is always on the RIGHT), so we had better start a new pair for nextstorm, which invloves making the nextstorm the left-hand storm in the new pair.

	
				
                  if(oriflag.eq.0) then !There are 2 possible orientations of a pair along the EW directions. Oriflag=0 means A/cyc on the left.
			signflag=1. !RH Storm at i,j, is a cyclone.Signflag is the sign of q' that the i,j storm will be (signe=-signflag is what is assigned to the storm that is about to be injected at i,j.)
			else
			signflag=-1. !RH Storm at i,j, is an anticyclone.
			end if

			storm_edgeflag=0.
			
			nextstormi=nint((NX-1.)*rand()) !Pick random location.
                  nextstormj=nint((NY-1.)*rand())

			nextstormi=mod(nextstormi+REAL(NX),REAL(NX))
			nextstormj=mod(nextstormj+REAL(NY),REAL(NY))

			if((nextstormi-nint((2.*xo*NX/LX)+1.)).lt.0) then !asking the question ' is that location (nextstormi - (2xo+1)) OK to be the left storm?'. If no, then do this:
			storm_edgeflag=1.
                  nextstormi=nextstormi !Means the left hand storm is at the location nextstormi, as chosen above.
                  icentral=nint(nextstormi+nint((xo*NX/LX)+1.)) !Choose the midpoint of the pair, which is where zeta is measured.
                  icentral=nint(icentral-(NX-1.)*
     *            (floor(icentral/(NX-1.))))

			else !If (nextstormi - (2xo+1)) is OK to be the left hand storm, then

                  nextstormi=nextstormi-nint((2.*xo*NX/LX)+1.) !Now that we've tested above that the position here is in the domain, put the left hand most storm there.
                  nextstormi=mod(nextstormi+REAL(NX),REAL(NX))
                  icentral=nint(nextstormi+nint((xo*NX/LX)+1.))
                  icentral=nint(icentral-(NX-1.)*
     *            (floor(icentral/(NX-1.))))
			end if


			
                  oriflag=MOD(oriflag+1.,2.) !Cycle the oriflag

            else !If the current storm is the left hand one, then we need to add its neighbour, which is to its right, and on the same y gridpoint.
      
                  nextstormj=j
                  nextstormi=i+nint((2.*xo*NX/LX)+1.) !Always what should be done, as the LH storm is always what's come first!

                  if(nextstormi.gt.(REAL(NX)-1.)) then !This shoudn't happen, as the previously chosen storm should always be to the left of the domain edge, so the right hand storm should always fit. But, just in case, swap it around.
                  nextstormi=i-nint((2.*xo*NX/LX)+1.)
                  write(6,*) 'outrageous', i,j,time_step,nextstormi,
     *				(i+nint((2.*xo*NX/LX)+1.))
                  end if


                  nextstormi=mod(nextstormi+REAL(NX),REAL(NX)) !Make sure they fit in the domain.
                  nextstormj=mod(nextstormj+REAL(NY),REAL(NY))


                       if(oriflag.eq.0) then !There are 2 possible orientations of a pair along the E-W directions. Oriflag=0 means A/cyc on the left.
				signflag=-1. !LH Storm at i,j, is an anticyclone. This is consistent with other if statement like this one above.
				else
				signflag=1.!LH Storm at i,j, is a cyclone. So storm on the right will be a anticyclone.
				end if


            end if


	return
	end

!!!************************************************************
! ! !******** storm strength chooser **************************
!!!************************************************************
!This is the subroutine where condition 3.8 is tested, and applied if necessary.
!This is only applied once per pair. So it's when LH storm at i,j, and the RH storm is proposed at nextstormi,nextstormj.
      subroutine storm_strength_chooser(i,j,
     *		new_strength_acyc,W)
      
      include 'qgbt_diablo.h'

      real*8 strength_acyc,strength_cyc,delta_q_acyc,
     *	delta_q_cyc,new_strength_acyc,
     *	prev_thresh_check,thresh_check,overlap_flag,
     *	delta_q,rsqd,req_delta_q,new_strength,
     *	maxdepth,scale_factor,strengthfactor,
     *	W(0:nx+1,0:ny+1,0:1),scale_factor_used,
     *	stormtype

      
      integer i,j,stormpair,jmid,imid,ii,jj,ii_use,jj_use


            if((depthnow-(H(icentralnow,j,0))*ldsqd)+
     *	(psilim/(2.*bias_factor)) ! This is so that we can check if it's worth bothering checking the injection site first, or not! (so checking that zeta>-1/2)
     *      .lt.0.) then
     
     
		new_strength_acyc=0.
		new_strength_cyc=0.
		prev_thresh_check=0.!-q_thresh !setting this up for censorship later.
		thresh_check=0.
		overlap_flag=0.
		scale_factor_used=1.

	   do stormpair=0,1 !Once for each storm in the pair - idea is to find out, and implement, 3.8 by finding the gridpoint (which is part of the two storms) that requires it most.

	   if(stormpair.eq.0) then
	   jmid=j
	   imid=i
	   stormtype=-(location2(i,j,0)) !stormtype is the same sign as q'.
		if(stormtype.lt.0) then
	      strength_acyc=strengthfactor(-stormtype,depthnow,i,j) !Find out the proposed cyclone and anticyclone strengths.
	      strength_cyc=strengthfactor(stormtype,depthnow,
     *	nint(nextstormi),nint(nextstormj))
	      else
	      strength_acyc=strengthfactor(stormtype,depthnow,
     *	nint(nextstormi),nint(nextstormj)) !Find out the proposed cyclone and anticyclone strengths.
	      strength_cyc=strengthfactor(-stormtype,depthnow,i,j)


		end if
		delta_q_acyc=2.*strength_acyc*(SR/CMC)*to !work out their associated delta_q values for the centre of the vortices.
		delta_q_cyc=2.*strength_cyc*(SR/CMC)*to
	   else
	   jmid=nextstormj
	   imid=nextstormi
	   stormtype=-stormtype
	   end if

	   if(stormtype.lt.0) then 
	   delta_q=delta_q_acyc
	   else
	   delta_q=delta_q_cyc
	   end if


! 	First bit is to check that no gridpoints are going to clash with any others from a storm that is currently being injected.
	    do jj=(jmid-CEILING(2.*NY*Xo/LY)), 
     *		(jmid+CEILING(2.*NY*Xo/LY)) !Go over a square surrounding centre of the vortex.
	       do ii=(imid-CEILING(2.*NX*Xo/LX)),
     *		(imid+CEILING(2.*NX*Xo/LX))

		rsqd=floor(10000.*((((jj-jmid)**2.+(ii-imid)**2.))*
     *		((LX/NX)**2.)))/10000. !Compute radius from centre, to see if gridpoint is inside circle.

            jj_use=mod(jj+NY,NY) !We want jj and ii to be distance measurers, so it's fine if they go over the domain size, so define jj_use for indexing.
            ii_use=mod(ii+NX,NX)

		if(rsqd.lt.(4.*(xo**2.))) then !If we're inside the circle of the storm refered to by stormpair, then:
		
		
		if(((.not.((ii_use.eq.i).and.(jj_use.eq.j))).and. !Firstly check that there aren't currently any other storms (other than the ones in this pair) currently being injected within the circle. And if there is another storm there, then stop this pair being injected.
     *	(.not.((ii_use.eq.nextstormi).and.(jj_use.eq.nextstormj))))
     *	.and.(abs(location2(ii_use,jj_use,0)
     *	).ne.0.)) then
		new_strength_acyc=0.
		new_strength_cyc=0.
		overlap_flag=1. !Means that we don't check any of the other points once we know storm will be rejected
		mvcpflag=1. !Means that we will apply these new strengths (of 0) when they come to be injected.
		nextstormi=nint(NX*rand())
            nextstormj=nint(NY*rand())
		nextstormi=mod(nextstormi+REAL(NX),REAL(NX))
		nextstormj=mod(nextstormj+REAL(NY),REAL(NY))
		nextstorm=time_step+2
		pairflag=-1.
		storm_edgeflag=0.

! START COPIED BIT!!!!******************************************************************************************************************
! BIT COPIED FROM STORM_LOCATION_CHOOSER, as the above bit of code chooses a new left hand storm, but we have to be sure it'll work!
			if((nextstormi-nint((2.*xo*NX/LX)+1.)).lt.0) then
			storm_edgeflag=1.
                  nextstormi=nextstormi
                  icentral=nint(nextstormi+nint((xo*NX/LX)+1.))
                  icentral=nint(icentral-(NX-1.)*
     *            (floor(icentral/(NX-1.))))
			else

                  nextstormi=nextstormi-nint((2.*xo*NX/LX)+1.) !Go for the Left hand storm.
                  nextstormi=mod(nextstormi+REAL(NX),REAL(NX))
                  icentral=nint(nextstormi+nint((xo*NX/LX)+1.))
                  icentral=nint(icentral-(NX-1.)*
     *            (floor(icentral/(NX-1.))))
			end if
! END COPIED BIT!!!!******************************************************************************************************************
		
		write(80,*) ii_use,i,jj_use,j,time_step,mvcpflag,
     *	nextstormi,nextstormj,pairflag,icentral !If this happens, log it.
		end if
		end if
		
		if((rsqd.lt.(xo**2.)).and.(overlap_flag.eq.0.)) then !Now check whether the max resulting abs(q') is less than threshold, and if not, where needs capping?

		thresh_check=(abs(W(ii_use,jj_use,1)+(delta_q*
     *			(1.-rsqd/(xo**2.))*stormtype))) ! What is the proposed new q'? 

		  if((thresh_check.le.q_thresh)) then
				if(mvcpflag.ne.1) mvcpflag=0 ! Do nothing to exisitng proposed vortex based on the value at this point. Only make mvcpflag=0 if it is not already been made to =1.

		  else !If condition 3.8 does need to be applied, then proceed as follows.


                  mvcpflag=1 !Set the flag. 
		req_delta_q = (q_thresh-abs(W(ii_use,jj_use,1)))/stormtype ! Find out how much delta q we need to make it up to the threshold.

		new_strength=req_delta_q/(2.*(SR/CMC)*to*stormtype*
     *			(1.-rsqd/(xo**2.))) ! Find out what that corresponds to in terms of non-dimensionalisation.
		
		if(new_strength.lt.0.) then !This only occurs if req_delta_q is less than zero, which means the background itself has already exceeded the threshold. Therefore, we stop the injection here.
		new_strength=0. 
		write(6,*) new_strength,q_thresh,W(ii_use,jj_use,1),
     *	rsqd,i,j,ii_use,jj_use,'prevented injection'
     		nextstormi=nint(NX*rand())
            nextstormj=nint(NY*rand())
		nextstormi=mod(nextstormi+REAL(NX),REAL(NX))
		nextstormj=mod(nextstormj+REAL(NY),REAL(NY))
		nextstorm=time_step+2
		pairflag=-1.
		storm_edgeflag=0.
! START COPIED BIT!!!!******************************************************************************************************************
! BIT COPIED FROM STORM_LOCATION_CHOOSER, as the above bit of code chooses a new left hand storm, but we have to be sure it'll work!
			if((nextstormi-nint((2.*xo*NX/LX)+1.)).lt.0) then
			storm_edgeflag=1.
                  nextstormi=nextstormi
                  icentral=nint(nextstormi+nint((xo*NX/LX)+1.))
                  icentral=nint(icentral-(NX-1.)*
     *            (floor(icentral/(NX-1.))))
			else

                  nextstormi=nextstormi-nint((2.*xo*NX/LX)+1.) !Go for the Left hand storm.
                  nextstormi=mod(nextstormi+REAL(NX),REAL(NX))
                  icentral=nint(nextstormi+nint((xo*NX/LX)+1.))
                  icentral=nint(icentral-(NX-1.)*
     *            (floor(icentral/(NX-1.))))
			end if
! END COPIED BIT!!!!******************************************************************************************************************
		end if


		if(stormtype.lt.0) then !If the cap that needs to be applied is due to the anticyclone, then we will need to make the anticyclone equal to new_strength, and scale the cyclone by the same amount, so we keep the fractional bias the same. So define a scale factor. 
		scale_factor=new_strength/strength_acyc
		else 
		scale_factor=new_strength/strength_cyc !If it's the other way around, and the cyclone has caused 3.8 to be applied, make a scale factor.
		end if

		if(scale_factor.lt.scale_factor_used) then  !For each storm, keep track of what the largest scale factor is, as that's what needs to be applied.
		scale_factor_used=scale_factor
		end if
		




		if(stormtype.lt.0) then !Apply the new scaings. 
		new_strength_acyc=strength_acyc*scale_factor_used
		new_strength_cyc=strength_cyc*scale_factor_used 
		else
		new_strength_cyc=strength_cyc*scale_factor_used
		new_strength_acyc=strength_acyc*scale_factor_used 
		end if
		end if
	      if(ABS(W(ii_use,jj_use,1)).gt.maxdepth) maxdepth=
     *		ABS(W(ii_use,jj_use,1)) !This is just a diagnostic.

            end if

	    end do
	    end do
	    end do
	    write(79,*) i,j,(W(i,j,1)),nextstormi,nextstormj, !Each time we run the check, output results.
     *      (W(nint(nextstormi),nint(nextstormj),1)),mvcpflag,
     *      time_step,pairflag,maxdepth,strength_acyc,strength_cyc,
     *	new_strength_acyc,new_strength_cyc,overlap_flag,nextstorm
		flush(79)
	    end if


	return
	end

! ! !******** energy calculator************************** !Calculates energy rate of change in real time via a few different methods.
      subroutine compute_energy_rate(k)
! ! 
      include 'qgbt_diablo.h'
      
      real*8 e_rate_temp_now(0:2),temp1,psisq(0:2),psiysq(0:2),
     *	psixsq(0:2),
     *	e_rate_temp_now_u(0:2),e_rate_temp_now_v(0:2),
     *	e_rate_temp_now_psi(0:2)
      
      integer k,i,j,type_run
	psisq(0)=0.
	psiysq(0)=0.
	psixsq(0)=0.
	psisq(1)=0.
	psiysq(1)=0.
	psixsq(1)=0.
	psisq(2)=0.
	psiysq(2)=0.
	psixsq(2)=0.

	do type_run=0,1
      do j=0,tnky
      do i=type_run,nkx
            temp1 = ksqd(i,j)
            psisq(type_run) = psisq(type_run) + dble( CP(i,j,k)
     *             * conjg(CP(i,j,k)) )
            psiysq(type_run) = psiysq(type_run) + ky2(j)*dble( CP(i,j,k) 
     *              * conjg(CP(i,j,k)) )
            psixsq(type_run) = psixsq(type_run) + kx2(i)*dble( CP(i,j,k)
     *              * conjg(CP(i,j,k)) )
	end do
	end do
	end do

      do j=0,nky
         i=0
            temp1 = ksqd(i,j)
            psisq(2) = psisq(2) + dble( CP(i,j,k)
     *             * conjg(CP(i,j,k)) )
            psiysq(2) = psiysq(2) + ky2(j)*dble( CP(i,j,k) 
     *              * conjg(CP(i,j,k)) )
            psixsq(2) = psixsq(2) + kx2(i)*dble( CP(i,j,k)
     *              * conjg(CP(i,j,k)) )
	end do


	do type_run=0,2
! 	Total
	energy_old(type_run)=energy_now(type_run)
	energy_now(type_run)=energy_new(type_run)
	energy_new(type_run)=(psiysq(type_run)+psixsq(type_run)+
     *(psisq(type_run)/LDSQD))


	if(energy_old(type_run).ne.0.) then
	e_rate_temp_now(type_run)=(energy_new(type_run)-
     * energy_old(type_run))/(2.*delta_t)

	e_rate_average(type_run)=e_rate_average(type_run)+
     * e_rate_temp_now(type_run)
	end if
! 	ZONAL
	energy_u_old(type_run)=energy_u_now(type_run)
	energy_u_now(type_run)=energy_u_new(type_run)
	energy_u_new(type_run)=(psiysq(type_run))


	if(energy_u_old(type_run).ne.0.) then
	e_rate_temp_now_u(type_run)=(energy_u_new(type_run)-
     *energy_u_old(type_run))/(2.*delta_t)

	e_rate_average_u(type_run)=e_rate_average_u(type_run)+
     *e_rate_temp_now_u(type_run)
	end if
! 	Meridional

	energy_v_old(type_run)=energy_v_now(type_run)
	energy_v_now(type_run)=energy_v_new(type_run)
	energy_v_new(type_run)=(psixsq(type_run))


	if(energy_v_old(type_run).ne.0.) then
	e_rate_temp_now_v(type_run)=(energy_v_new(type_run)-
     * energy_v_old(type_run))/(2.*delta_t)

	e_rate_average_v(type_run)=e_rate_average_v(type_run)+
     * e_rate_temp_now_v(type_run)
	end if
	! 	potential

	energy_psi_old(type_run)=energy_psi_now(type_run)
	energy_psi_now(type_run)=energy_psi_new(type_run)
	energy_psi_new(type_run)=psisq(type_run)


	if(energy_psi_old(type_run).ne.0.) then
	e_rate_temp_now_psi(type_run)=(energy_psi_new(type_run)-
     *energy_psi_old(type_run))/(2.*delta_t)

	e_rate_average_psi(type_run)=e_rate_average_psi(type_run)+
     * e_rate_temp_now_psi(type_run)
	end if
	
	enddo
	
	
	
	
	

      return
      end
      

!*********************************************************************************
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      function ranf()
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      double precision g05caf
      ranf=dble(rand())
      return
      end


