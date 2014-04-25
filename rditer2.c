/**
 * @file rditer2.c
 * @brief Calculate the disc scale radius via an iterative process
 * @author Tomas E. Tecce
 */
#include "allvars.h"
#include "proto.h"
#include "readparameters.h"

#define RADMIN -4.0
#define RADMAX 1.0
#define RDLOGSTEP 5e-3
#define RDSTEPS (int)ceil((RADMAX-RADMIN)/RDLOGSTEP)
#define RFMAXITER 100                  /* Max iterations for root finder */
#define RFTOL 1e-4                        /* Root finder error tolerance */
#define RDDIVERGE 1e6                               /* Stability control */
#define RDMAXITER 100                 /* Max iterations for scale radius */
#define RDTOL 1e-2                             /* Scale radius tolerance */


/**
 * @brief Calculate the scale radius for a galaxy using an iterative method
 * @param chalo Concentration of host halo (assuming an NFW profile)
 * @param md Ratio of disc mass to halo mass
 * @param lambda Spin parameter of host halo
 * @param mb Ratio of bulge mass to halo mass
 * @param zcurr Current redshift
 * @param guess Initial guess value for root finder
 * @param *rdresult Pointer to store the scale radius obtained
 * @param *vcresult Pointer to store the circular velocity obtained
 * @param *roptresult Pointer to store the optical radius obtained
 */
void rditer(double chalo, double md, double lambda, double mb,
	    double zcurr, double guess, double *rdresult, double *vcresult, 
	    double *roptresult)
{
  int i, j, count, check, indexropt;
  unsigned long k;
  double rhalo;
  double grav;
  double rd, fr, fc, delta0, rdprev, m0nfw;
  double ropt;
  /*double mopt;*/
  double xp, xn, xnn, fxn, fxp;
  double vc2dm, vc2bulge;
  double trap, f_err;
  double y;
  double *r, *rinit, *md_of_r, *mf_of_r, *vcdisc, *vc;
  double *fxi;
  double retry;

  /* 
   * System of units:
   *   - Unit of mass is the halo virial mass corresponding to the given
   *     concentration
   *   - Unit of density: universal critical density for collapse
   *   - Unit of velocity: halo virial velocity G M_200 / r_200
   * In this system, the unit of length is (M_200 / rho_c)^1/3 and the
   * virial radius r_200 is equal to the gravitational constant G
   * (numerically).
   */
  
  grav = rhalo = pow(3.0/(800*M_PI), 1.0/3.0);
  
  if (BulgeAngularMomentumOn == 1)
    lambda *= 1.0 + mb/md;

#ifdef TESTMODE
  printf("c200   = %g\n"
	 "md     = %g\n"
	 "mb     = %g\n"
	 "lambda = %g\n"
	 "z      = %g\n",
	 chalo, md, mb, lambda, zcurr);
  fflush(stdout);
#endif

  r       = dvector(1, RDSTEPS);
  rinit   = dvector(1, RDSTEPS);
  md_of_r = dvector(1, RDSTEPS);
  mf_of_r = dvector(1, RDSTEPS);
  vcdisc  = dvector(1, RDSTEPS);
  vc      = dvector(1, RDSTEPS);
  fxi     = dvector(1, NPOINTSGL);

  /* Points where the rotation curve will be calculated */
  for (i=1; i<=RDSTEPS; i++)
    r[i] = rhalo*pow(10.0, RADMIN+i*RDLOGSTEP);
  
    fr = guess;                       /* Initial guess for the iteration */
    fc = 2.0/3.0 + pow(chalo/21.5, 0.7);        /* MMW 98 approx. for fc */
    rd = (1.0/sqrt(2))*lambda*rhalo*fr/sqrt(fc);
    rdprev = rd/2;
    delta0 = (200.0/3.0)*chalo*chalo*chalo/(log(1+chalo) -
					    chalo/(1+chalo));
    m0nfw = 3*delta0/(200*chalo*chalo*chalo);

    count = 0;
#ifdef TESTMODE
    printf("Initial guess for rd: %g (fr = %g)\n\n", rd, fr);
    printf("-------------------------------------\n"
	   "| Iteration | rd       | diff       |\n"
	   "-------------------------------------\n"
	   "| %3d %15g  %12g |\n",
	   count, rd, fabs(rd-rdprev)/rdprev);
#endif

    /* Start iterating until rd converges */
    do {
      for (i=1; i<=RDSTEPS; i++) {
	/* First, find disc mass with the guessed rd */
	md_of_r[i] = md*(1.0 - (1+r[i]/rd)*exp(-r[i]/rd));

	/* Find rinit(r) with a nonlinear root finder */
	xp = pow(10.0, RADMIN+(i-0.5)*RDLOGSTEP);
	xn = pow(10.0, RADMIN+(i+0.5)*RDLOGSTEP);

	/* We need to make sure that a root is bracketed by xp and xn */
	check = dzbrac(&xp, &xn, m0nfw, rhalo, chalo, r[i], md,
		       md_of_r[i], mb);

	if (check == 0)
	  rinit[i] = r[i];

	j = 0;
	do {
	  fxn = m_of_rinit(xn, m0nfw, rhalo, chalo, r[i], md, md_of_r[i],
			   mb);
	  fxp = m_of_rinit(xp, m0nfw, rhalo, chalo, r[i], md, md_of_r[i],
			   mb);
	  xnn = xn - fxn*(xn-xp)/(fxn-fxp);
	  xp = xn;
	  xn = xnn;
	  j++;	  
	} while(fabs(fxn-fxp)/fxp >= RFTOL && j < RFMAXITER);
	rinit[i] = xn*rhalo;
	if (isnan(rinit[i])) rinit[i] = r[i];    /* PATCH */
      }

      /* Once rinit is determined, calculate mass profile after adiabatic
	 contraction */
      for (i=1; i<=RDSTEPS; i++)
	mf_of_r[i] = md_of_r[i] + mb + mass_nfw(rinit[i]/rhalo,
						m0nfw,
						chalo)*(1-md-mb);

      /* Determine the rotation curve */
      /* Note: values in vectors vcdisc and vc are actually in units of
	 velocity^2 */
      for (i=1; i<=RDSTEPS; i++) {
	/* Contribution from the exponential disc (Binney & Tremaine 1st
	   ed. p. 77). Note: for large y, the product I*K -> 0 */
	y = r[i]/2.0/rd;
	if (y < 91) {
	  vcdisc[i] = (2.0*grav*md/rd)*y*y*
	    ((dbessi0(y)*dbessk0(y)) - (dbessi1(y)*dbessk1(y)));
	}
	else {
	  vcdisc[i] = 0.0;
	}

	/* Add the contribution of the DM and bulge in quadrature to get
	   the full rotation curve */
	vc2dm = grav*mass_nfw(rinit[i]/rhalo, m0nfw, chalo)/r[i];
	vc2bulge = grav*mb/r[i];
	vc[i] = vc2dm + vcdisc[i] + vc2bulge;
	if (vc[i] < 0) vc[i] = 0.0;
      }

      /* Integration of the rotation curve:
	 Integrate exp(-u)*u^2*vc(rd*u) by Gauss-Laguerre quadrature. Need
	 to interpolate to find f(xi) */
      trap = 0;
      for (i=1; i<=NPOINTSGL; i++) {
	dlocate(r, RDSTEPS, rd*Xgl[i], &k);

	if (k >= RDSTEPS-2)
	  fxi[i] = vc[RDSTEPS];
	else if (k <= 2)
	  fxi[i] = vc[1];
	else
	  dpolint(&r[k-2], &vc[k-2], 4, rd*Xgl[i], &fxi[i], &f_err);
	trap += Wgl[i]*sqrt(fxi[i]);
      }
	
      fr = 2.0/trap;

      /* With the new value of fr, recalculate the scale radius */
      rdprev = rd;
      rd = (1.0/sqrt(2))*lambda*rhalo*fr/sqrt(fc);
      count++;

#ifdef TESTMODE
      printf("| %3d %15g  %12g |\n", count, rd, fabs(rd-rdprev)/rdprev);
#endif
    } while (fabs(rd-rdprev)/rdprev >= RDTOL && count < RDMAXITER &&
	     fabs(rd-rdprev) < RDDIVERGE);

    if (count >= RDMAXITER) {
      /* 2010-12-15 I have found that in a few cases, for a certain 
	 combination of input parameters, the iterative method oscillates
	 between two values and never achieves convergence. Nudging the 
	 value of one parameter a little seems to solve the problem.
	 Call rditer recursively, shifting the value of lambda by a small 
	 amount, until the result converges */
#ifdef TESTMODE
      printf("Maximum number of iterations reached without "
	     "achieving convergence for rd.\n\n");
      printf("rd = %g, fabs(rd-rdprev) = %g\n", rd, fabs(rd-rdprev));
#endif
#ifdef DEBUG
      printf("    Warning: max iterations %d reached for: \n"
	     "      z      = %g\n"
	     "      c200   = %g\n"
	     "      md     = %g\n"
	     "      mb     = %g\n"
	     "      lambda = %g\n"
	     "      rdprev = %g\n"
	     "      rd     = %g\n\n", 
	     RDMAXITER, zcurr, chalo, md, mb, lambda, rdprev, rd);
#endif
      Count_itmax++;
      /* Safety check */
      /*if (Count_itmax > 10000) {
	fprintf(stderr, "Error in rditer: too many retries to fix a max "
		"iterations problem - Exit (A)\n"); fflush(stderr);
	exit(EXIT_FAILURE);
	}*/
      /* *rdresult = RDUNDEFINED;
      *vcresult = RDUNDEFINED;
      *roptresult = RDUNDEFINED;*/
      
      retry = lambda*1.01;
      if (retry/lambda > 1.1) {
	fprintf(stderr, "Error in rditer: too many retries to fix a max "
		"iterations problem - Exit\n"); fflush(stderr);
	exit(EXIT_FAILURE);
      }
      rditer(chalo, md, retry, mb, zcurr, guess, rdresult, vcresult, 
	     roptresult);
#ifdef DEBUG
      printf("Corrected value: rd = %g for lambda = %g\n\n",
	     *rdresult, retry); fflush(stdout);
#endif
    }                                  /* End if max iterations exceeded */
    else {
      if (fabs(rd-rdprev) < RDDIVERGE) {
#ifdef TESTMODE
	printf("-------------------------------------\n");
	printf("Convergence for rd reached after %d iterations.\n\n",
	       count);
      
	printf("*******\n"
	       " Values found (in internal calculation units)\n"
	       "*******\n");
	printf("Scale radius: %g (f_R = %g)\n", rd, fr);

#endif
	/* Find the disc's optical radius (user-defined as a percentage of
	       total baryonic mass) */
	/*mopt = MassFractionForVc*(md+mb);
	md_of_r[1] = 0;
	for (i=1; i<=RDSTEPS; i++) {
	  md_of_r[i] = md*(1.0 - (1+r[i]/rd)*exp(-r[i]/rd));
	  if (md_of_r[i]+mb > mopt && md_of_r[i-1]+mb <= mopt) {
	    ropt = r[i];
	    indexropt = i;
	  }
	  }*/
	/* TOMAS 2010-12-29: if the total baryonic mass is used, it may
	   happen that the bulge mass alone is MassFractionForVc or greater.
	   When this happens, vc=0. An alternative is to consider the
	   velocity at a given radius (in units of the scale radius) */
	indexropt = RDSTEPS;   /* Initialise using the largest
				  possible value */
	ropt = OpticalRadius*rd;
	for (i=1; i<=RDSTEPS; i++) {
	  if (r[i] > ropt && r[i-1] <= ropt) {
	    indexropt = i;
	  }
	}

	*rdresult = rd;
	*vcresult = sqrt(vc[indexropt]);
	*roptresult = ropt;

#ifdef TESTMODE
	/*printf("Radius for %g per cent of mass: %g \n",
	  MassFractionForVc*100, *roptresult);
	printf("Circular velocity at %g per cent of mass: %g \n\n",
	MassFractionForVc*100, *vcresult);*/
	printf("Selected radius for flat rotation curve: %g\n", ropt);
	printf("Circular velocity: %g \n\n", *vcresult);
#endif
      }  /* END if result does not diverge */
      else {
	printf("WARNING: divergence\n\n");
	Count_div++;
	*rdresult = RDUNDEFINED;
	*vcresult = RDUNDEFINED;
	*roptresult = RDUNDEFINED;
      }
      
    }                         /* END if max iterations were not exceeded */

    free_dvector(r, 1, RDSTEPS);
    free_dvector(rinit, 1, RDSTEPS);
    free_dvector(md_of_r, 1, RDSTEPS);
    free_dvector(mf_of_r, 1, RDSTEPS);
    free_dvector(vcdisc, 1, RDSTEPS);
    free_dvector(vc, 1, RDSTEPS);
    free_dvector(fxi, 1, NPOINTSGL);
    
    return;
}


/*
 * Mass profile at the radius x where the mass currently at r was
 * before adiabatic contraction.
 * a = 4*PI*delta0*Rhocrit*rs*rs*rs
 * b = rhalo
 * c = chalo
 * d = r (where the mass is now)
 * e = md
 * f = Md_of_r
 * g = mb
 */
double m_of_rinit(double x, double a, double b, double c, double d, double e,
		  double f, double g)
{
  return a*(1.0/(1+c*x) - 1 + log(1+c*x))*((x*b/d) - 1 + e + g) - f - g;
}


/*
 * Mass within x = r/r200 in a NFW halo of concentration c
 */
double mass_nfw(double x, double a, double c)
{
  return a*(1.0/(1+c*x) - 1 + log(1+c*x));
}


/*
 * Calculate scale radius using the approximations by MMW
 */
double rdaprox(double rvir, double lambda, double c, double mhalo, double md)
{
  double fc, fr;

  fc = 2.0/3.0 + pow(c/21.7, 0.7);

  fr = pow(lambda/0.1, -0.06+2.71*md+0.0047/lambda) * (1 - 3*md + 5.2*md*md)
    * (1 - 0.019*c + 0.00025*c*c + 0.52/c);

  return rvir*lambda*fr/sqrt(2*fc);
}


/**
 * @brief Find values that bracket a root
 * 
 * Given a function func and an initial guessed range x1 to x2, the
 * routine expands the range geometrically until a root is bracketed 
 * by the returned values x1 and x2 (in which case zbrac returns 1)
 * or until the range becomes unacceptably large (in which case 
 * zbrac returns 0).
 * This is the double precision version of routine zbrac from NR,
 * modified to use the function m_of_rinit
 */
int dzbrac(double *x1, double *x2, double a, double b, double c, 
	   double d, double e, double f, double g)
{
  int j;
  double f1, f2;

  if (*x1 == *x2) nrerror("Bad initial range in dzbrac");
  f1 = m_of_rinit(*x1, a, b, c, d, e, f, g);
  f2 = m_of_rinit(*x2, a, b, c, d, e, f, g);
  //printf("    In dzbrac: f1, f2 %g, %g\n", f1, f2);
  for (j=1; j<=ZBRAC_NTRY; j++) {
    if (f1*f2 < 0.0) return 1;
    if (fabs(f1) < fabs(f2)) {
      *x1 += ZBRAC_FACTOR*(*x1-*x2);
      //if (*x1 < 0) *x1 = 0;
      f1 = m_of_rinit(*x1, a, b, c, d, e, f, g);
    }
    else {
      *x2 += ZBRAC_FACTOR*(*x2-*x1);
      //if (*x2 < 0) *x2 = 0;
      f2 = m_of_rinit(*x2, a, b, c, d, e, f, g);
    }
    //printf("              correcting x1, x2, f1, f2 %g, %g, %g, %g\n", *x1, *x2, f1, f2);
  }
  return 0;
}


/*
 * To test performance of the program: user inputs halo concentration, md,
 * mb, lambda and percentage of disc mass corresponding to optical radius.
 * Program returns scale radius, optical radius, circular velocity at
 * optical radius and exits.
 * 2011-05-19: if the total baryonic mass is used, it may
 * happen that the bulge mass alone is MassFractionForVc or greater.
 * When this happens, vc=0. An alternative is to consider the
 * velocity at a given radius (in units of the scale radius) 
 */
#ifdef TESTMODE
void test_mode(void)
{
  double mhalo, chalo, md, mb, lambda, rscale, ropt, vcopt, zcurr;
  double hubble_of_z, rhocrit, grav, rap, rhalo;
  
  printf("Redshift?\n");
  scanf("%le",&zcurr);
  printf("Halo concentration?\n");
  scanf("%le", &chalo);
  printf("Disc mass fraction?\n");
  scanf("%le", &md);
  printf("Spin parameter?\n");
  scanf("%le", &lambda);
  printf("Bulge mass fraction?\n");
  scanf("%le", &mb);
  /*printf("Mass fraction at optical radius?\n");
    scanf("%le", &fopt);*/
  printf("\n");

  /* Universal critical density */
  hubble_of_z = HUBBLE*HUBBLE*(OmegaLambda +
			       ((1.0-OmegaLambda-Omega0) +
				Omega0*(1+zcurr))*(1+zcurr)*(1+zcurr));
  rhocrit = 3.0*hubble_of_z/(8.0*M_PI*GGRAV);      /* Now in Msun Mpc^-3 */
  rhocrit *= pow(UnitLength_in_kpc/1.0e3,3)/UnitMass_in_Msun;

  grav = GGRAV*1e3*UnitMass_in_Msun/                  /* G in code units */
    UnitLength_in_kpc/
    (UnitVelocity_in_km_per_s*UnitVelocity_in_km_per_s);
  
  /* Use the concentration tables to find the halo mass M_200 that
     corresponds to the given concentration. This function returns the halo
     mass in Msun */
  get_mass_from_concentration(chalo, zcurr, &mhalo);
  mhalo /= UnitMass_in_Msun;
  rditer(chalo, md, lambda, mb, zcurr, 1.0, &rscale, &vcopt, &ropt);

  /* Convert the values returned by the rditer function 
     to physical units */
  rhalo = pow(3.0*mhalo/(800*M_PI*rhocrit), 1.0/3.0);
  rscale *= pow(mhalo/rhocrit, 1.0/3.0);
  ropt   *= pow(mhalo/rhocrit, 1.0/3.0);
  vcopt  *= sqrt(grav*mhalo/rhalo);

  printf("\n*******\n"
	 "Characteristics of the resultant disc galaxy\n"
	 "*******\n");
  printf("Assuming a halo mass of %g M_Sun/h for c = %g at z = %g\n",
	 mhalo*UnitMass_in_Msun, chalo, zcurr);
  printf("Virial radius and velocity: %g kpc/h, %g km/s\n\n",
	 rhalo, sqrt(grav*mhalo/rhalo));
  printf("Scale radius: %g kpc/h \n", rscale*UnitLength_in_kpc);
  /*printf("Radius for %g per cent of mass: %g kpc/h\n",
    MassFractionForVc*100, ropt*UnitLength_in_kpc);
    printf("Circular velocity at %g per cent of mass: %g km/s\n\n",
    MassFractionForVc*100, sqrt(vcopt)*UnitVelocity_in_km_per_s);*/
  printf("Radius for flat rotation curve: %g kpc/h \n", 
	 ropt*UnitLength_in_kpc);
  printf("Circular velocity: %g km/s \n\n", vcopt);

  /* Compare with the MMW 98 approximation */
  rap =  rdaprox(rhalo, lambda, chalo, mhalo, md);
  printf("Using the MMW 98 approximation: scale radius %g kpc/h\n",
	 rap*UnitLength_in_kpc);
  printf("Relative difference: %g per cent\n",
  fabs(rscale-rap)*100/rscale);

  printf("\nTest mode complete\n");
  fflush(stdout);

  free_arrays();
  free_concentration_tables();

  exit(EXIT_SUCCESS);
}
#endif

