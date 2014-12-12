c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c This is the header file for qgbt_diablo.               version 1.0
c This file contains definitions of global parameters and global variables.
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      
      IMPLICIT NONE

c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c Spatial resolution parameters
c (We hardwire these into the code so that the compiler may perform
c  optimizations based on the grid size at compile time).
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
      
      integer nx, ny
      include 'grid_def'

c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c Input parameters and runtime variables
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      real*8 kappa, nu, beta, hx, hy, lx, ly, csx, csy, delta,
     *       bc1_1, bc1_2, bc1_3, bc1_4,
     *       bc2_1, bc2_2, bc2_3, bc2_4,
     *       delta_t, ic_1, ic_2, energy, hamp,ldsqd,
     *	 kt,sr,cmc,a,b,to,xo,psilim,tcool,bias_factor,
     *	 time_bet_storms,q_thresh

      integer n_time_steps, num_per_dir, time_ad_meth, verbosity,
     *    save_flow_int, save_phys_int, save_prof_int, save_stats_int, 
     *    ic_type, bc1_type, bc2_type, previous_time_step,
     *    kcut, filter_exp

      logical use_new_storms

      common /input/
     *        kappa, nu, beta, hx, hy, lx, ly, csx, csy, 
     *        delta_t, ic_1, ic_2,
     *        bc1_1, bc1_2, bc1_3, bc1_4,
     *        bc2_1, bc2_2, bc2_3, bc2_4, energy, hamp,
     *        n_time_steps, num_per_dir, time_ad_meth, verbosity,
     *        save_flow_int, save_phys_int, save_prof_int, 
     *        save_stats_int, ic_type, bc1_type, bc2_type, 
     *        previous_time_step, kcut, filter_exp,ldsqd,
     *	  kt,sr,cmc,a,b,to,xo,psilim,tcool,bias_factor,
     *	  time_bet_storms,q_thresh,use_new_storms

      real*8 time
      integer time_step, advec_step, nxm, nym, tnky
      common /runtime_vars/ time, time_step, advec_step, nxm, nym, tnky

c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c Adams-Bashforth parameters and depth parameters
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
  
      real*8 abp(3), delt1, delt2, rob
      integer new, now, old
      common /advec_params/ abp, rob, new, now, old
 
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c FFT parameters
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      real*8 kx(0:nx/2), ky(0:2*(ny/2)), kx2(0:nx/2), ky2(0:2*(ny/2)),
     *       ksqd(0:nx/2,0:2*(ny/2)), pi, eps, rnx, rny


      complex*16 cikx(0:nx/2), ciky(0:2*(ny/2)), ci,
     *           cyx_plane(0:ny,0:nx/2)

      
      integer*8 fftw_x_to_p_plan, fftw_x_to_f_plan,
     *          fftw_y_to_p_plan, fftw_y_to_f_plan,
     *          fftw_z_to_p_plan, fftw_z_to_f_plan

      integer nkx, nky

      common /fft_params/
     *        kx, ky, kx2, ky2, ksqd,
     *        pi, eps, rnx, rny, cikx, ciky, ci, 
     *        fftw_x_to_p_plan, fftw_x_to_f_plan, 
     *        fftw_y_to_p_plan, fftw_y_to_f_plan, 
     *        fftw_z_to_p_plan, fftw_z_to_f_plan, 
     *        nkx, nky

c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c Grid parameters
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      real*8 gx(0:nx+1), gy(0:ny+1), dx(0:nx+1), dy(0:nx+1),
     *       gxf(0:nx) , gyf(0:ny) , dxf(0:nx) , dyf(0:ny)

      common /grid_params/ gx, gy, dx, dy, gxf, gyf, dxf, dyf
  
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|
c Global variables
c----*|--.---------.---------.---------.---------.---------.---------.-|-------|

      real*8     UG (0:nx+1,0:ny+1,0:2), VG (0:nx+1,0:ny+1,0:2),
     *           QX (0:nx+1,0:ny+1,0:2), QY (0:nx+1,0:ny+1,0:2),
     *           H  (0:nx+1,0:ny+1,0:2), filter(0:nx/2,0:ny+1),
     * 	     PSI(0:nx+1,0:ny+1,0:2), 
     *	     location2(0:nx+1,0:ny+1,0:3),
     *       UG2 (0:nx+1,0:ny+1,0:2), VG2 (0:nx+1,0:ny+1,0:2),
     *       UYY(0:nx+1,0:ny+1,0:2),psibar2(0:2),nextstorm
     *       ,nextstormi,nextstormj,pairflag,oriflag,
     *       upvpbar(0:ny+1),depthnow,signflag,storm_edgeflag,
     *	 storm_id,new_strength_cyc,e_rate_average(0:2),
     *	 energy_new(0:2),
     *	 energy_now(0:2),energy_old(0:2),e_rate_average_u(0:2),
     *	 energy_u_new(0:2),
     *	 energy_u_now(0:2),energy_u_old(0:2),
     *	 e_rate_average_v(0:2),energy_v_new(0:2),
     *	 energy_v_now(0:2),energy_v_old(0:2),
     *	 e_rate_average_psi(0:2),
     *	 energy_psi_new(0:2),energy_psi_now(0:2),
     *	 energy_psi_old(0:2)

      integer icentral,icentralnow,mvcpflag,injection_flag
      
      complex*16 CP (0:nx/2,0:ny+1,0:2),
     *           CQ (0:nx/2,0:ny+1,0:2), 
     *           CUG(0:nx/2,0:ny+1,0:2), CVG(0:nx/2,0:ny+1,0:2),
     *           CQX(0:nx/2,0:ny+1,0:2), CQY(0:nx/2,0:ny+1,0:2),
     *           CR (0:nx/2,0:ny+1,0:2),
     *           CH (0:nx/2,0:ny+1,0:2),
     *           CUG2(0:nx/2,0:ny+1,0:2), CVG2(0:nx/2,0:ny+1,0:2),
     *           CUYY(0:nx/2,0:ny+1,0:2), CR2(0:2),
     *           CUSQD(0:nx/2,0:ny+1,0:2),CVSQD(0:nx/2,0:ny+1,0:2),
     *           CPSQD(0:nx/2,0:ny+1,0:2)


      equivalence (UG,CUG), (VG,CVG),
     *            (QX,CQX), (QY,CQY), (H,CH),(UG2,CUG2),(VG2,CVG2),
     *            (UYY,CUYY)

      common /global_vars/ UG, VG, QX, QY, CP, CQ,
     *                     CR, H, filter, PSI, location2,CR2,psibar2,
     *			   nextstorm,nextstormi,nextstormj,pairflag,
     *                     oriflag,upvpbar,icentral,icentralnow,
     *                     depthnow,new_strength_cyc,signflag,
     *			   storm_edgeflag,storm_id,
     *			   injection_flag,mvcpflag,e_rate_average,
     *		         energy_now,energy_old,energy_new,
     *			   e_rate_average_u,energy_u_new,
     *	 		   energy_u_now,energy_u_old,
     *			   e_rate_average_v,energy_v_new,
     *	 		   energy_v_now,energy_v_old,
     *			   e_rate_average_psi,
     *	 		   energy_psi_new,energy_psi_now,
     *			   energy_psi_old
