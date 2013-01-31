//
//  DFDyn_constants.h
//  DigitalFluid
//
//  Created by G5 G5 on Tue May 10 2005.
//  Copyright (c) 2005 Ivano Ras. All rights reserved.
//

#import "DFDynamics.h"


// Macro definitions
// =================

// ---------------------------------------------------------------------------

// data_grid dim init
/*                               */ #define   Width       60
/*                               */ #define   Height      60
/*                               */ #define   Depth       1
/*                               */ #define   GridCells   Width*Height*Depth
/*                               */ #define   WH          Width*Height

/* delta_t init					    */ #define   dt         0.01

// ---------------------------------------------------------------------------

// Physical constants
// ••••••••••••••••••

/*							       [m]   */ #define   L		 1.0
/*                                 [m/s] */ #define   uinf   3
/* Water density                         */ #define   rinf   1000
//   (T= 4°C) = 1000    Kg/m^3

/* Water viscosity                       */ #define   mu	 0.001
//   (T=10°C) = 0.0013  Ns/m^2
//   (T=20°C) = 0.001   Ns/m^2 = 1 cP (centiPoise)
//   (T=25°C) = 0.00089 Ns/m^2
//   (T=30°C) = 0.0008  Ns/m^2 

/* dimensionless Reynolds number         */ #define   Re	 2
//   rinf * uinf * L / mu;

/* dimensionless x gravity [ not m/s^2 ] */ #define   gx	 0.0*L/(uinf*uinf)
/* dimensionless y gravity [ not m/s^2 ] */ #define   gy	 9.8*L/(uinf*uinf)
/* dt* = uinf*dt/L; dt*=dimensionless dt */ #define   c1_dt  1/dt

/* p* = (p-pinf)/(rinf*uinf^2);  dimensionless p */
/* t* = uinf*t/L;                dimensionless t */

/* ambient temperature [K]				*/ #define   Tamb    293

/* 1st constant of Fbuoy ** to decide!  */ #define   alpha   0.2
/* 2nd constant of Fbuoy ** to decide!  */ #define   beta    0.0

/*									    */ #define   DeltaX    L/(Width+1)
/*									    */ #define   DeltaY    DeltaX
/*									    */ #define   DeltaZ    DeltaX
//  DeltaY..Z must be equal to DeltaX (!) else DXDX must be changed...

// ---------------------------------------------------------------------------

// SOR constants  (Successive Over Relaxation)
// •••••••••••••

/*	max error accepted in SOR		 */ #define    eps      0.0001
/*	egative epsilon.				 */ #define   neps      -eps

/*	omega parameter of SOR method	 */ #define   omega     1.89
/*	quarter of omega				 */ #define   omq       omega/(-4)
/*									 */ #define   cDXDX_dt  DeltaX*DeltaX/dt

// ---------------------------------------------------------------------------

// Upwind constants
// ••••••••••••••••

/*									 */ #define   gamma     0.01
// gamma  >=  max(|u(i,j)*dt/DX|,|v(i,j)*dt/DY|);  for each i,j
// [According to Hirt et al., 1975]

// ---------------------------------------------------------------------------

// Derivative constants
// ••••••••••••••••••••

/*								*/ #define   c1_DXDXRe    1/(DeltaX*DeltaX*Re)
/*								*/ #define   c1_DYDYRe    1/(DeltaY*DeltaY*Re)
/*								*/ #define   c1_4DX       1/(4*DeltaX)
/*								*/ #define   c_gam_4DX    gamma/(4*DeltaX)
/*								*/ #define   c1_4DY       1/(4*DeltaY)
/*								*/ #define   c_gam_4DY    gamma/(4*DeltaY)
/*								*/ #define   c1_DX        1/(DeltaX)
/*								*/ #define   c1_DY        1/(DeltaY)
/*								*/ #define   c1_DZ        1/(DeltaZ)

/*								*/ #define   cdt_DX       dt/(DeltaX)
/*								*/ #define   cdt_DY       dt/(DeltaY)

/*								*/ #define   c1_2DX       1/(2*DeltaX)
/*								*/ #define   c1_2DY       1/(2*DeltaY)

/*								*/ #define   DXDX		  DeltaX*DeltaX

/*								*/ #define   c1_4DXDY     1/(4*DeltaX*DeltaY)

// ---------------------------------------------------------------------------

// Misc constants
// ••••••••••••••

/*								*/ #define   Wminus1       Width-1
/*								*/ #define   Wminus2       Width-2
/*								*/ #define   Hminus1       Height-1
/*								*/ #define   Hminus2       Height-2
/*								*/ #define   Dminus1       Depth-1

/*								*/ #define   BytePerColor  1

// ---------------------------------------------------------------------------

// SemiLagrangian constants
// ••••••••••••••••••••••••


// --->  AddForce constants

/*                               */  #define   dtgx         dt*gx
/*                               */  #define   dtgy         dt*gy


// --->  Advect constants

/*                               */  #define   halfdt       dt*0.5
/*                               */  #define   DeltaY_2     DeltaY/(2)


// --->  Diffuse constants

/* see above the water viscosity */  #define   visc			0.010

/*  omega1 = 1 => Gauss-Seidel   */  #define   omega1		1.00

/* a useful constant             */  #define   aa	        visc*dt/(DXDX)
/* another useful constant       */  #define   cc           1+4*aa
/* another useful constant       */  #define   om_cc        omega1/(cc)

/*                               */  #define   epsd         0.0000001
/*                               */  #define   nepsd		-epsd


// --->  PressureProject constants

/*                               */  #define   DeltaX_2     DeltaX/(2)


// ---------------------------------------------------------------------------

// Marker Particles constants
// ••••••••••••••••••••••••••

/*                               */  #define   numMarkers   1000



// ---------------------------------------------------------------------------

// Coniugated Gradients constants
// ••••••••••••••••••••••••••••••

/*  absolute tolerance  10 ^ -8  */  #define   epsCG       0.00000001

//  maximum number of iterations used to reach the
//  stopping criteria of sqrt(rho) <= (epsCG*bnorm)
/*                               */  #define   itermax     1000




// ---------------------------------------------------------------------------





@interface DFDynamics (DFDyn_constants)   // categories

- (void)constant_init;

@end
