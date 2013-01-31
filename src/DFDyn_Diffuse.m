//
//  DFDyn_Diffuse.m
//  DigitalFluid
//
//  Created by G5 G5 on Sun Jun 05 2005.
//  Copyright (c) 2005 Ivano Ras. All rights reserved.
//

#import "DFDyn_Diffuse.h"
#import "DFDyn_constants.h"
#import "DFDyn_Set_BC.h"


@implementation DFDynamics (DFDyn_Diffuse)


//  Diffusing Velocities
//  --------------------
//
//					    __2
//  (1)	 du/dt  =  visc \/ u  =  visc ( d2u/dx2 + d2u/dy2 )
//					    __2
//  (2)  dv/dt  =  visc \/ v  =  visc ( d2v/dx2 + d2v/dy2 )
//                                          laplacian
//
//  The implicit form for the 1st equation is:
//
//  (3)  u(n+1) i,j - visc dt( u(n+1) i-1,j + u(n+1) i+1,j +
//                             u(n+1) i,j-1 + u(n+1) i,j+1 -
//						     4 u(n+1) i,j) / DXDX  =  u(n) i,j
//
//  that handled become:
//
//  (4)  u(n+1) i,j * (4 + DXDX / visc dt) - u(n+1) i-1,j -
//             u(n+1) i+1,j - u(n+1) i,j-1 - u(n+1) i,j+1 = 
//							( u(n) i,j * DXDX / visc dt)
//
//  same thing for v(n+1).




- (void)Diffuse
{

struct Fdata        *q;        // point to a Fluid grid cell
unsigned long int   i,j;       // counters
double              un, und;   // u at step n;   und = u(n+1) - u(n);
double              vn, vnd;   // v at step n;   vnd = v(n+1) - v(n);
unsigned long int   max;       // max = boolean flag to exit from SOR loop



//  u  Diffusion
//  ------------

//  computation of rhs  =  ( u(n) i,j * DXDX / visc dt)
//  ееееееееееееееееее
//
// 	i=1,..,imax;   j=1,..,jmax;

q = Q + Width;					   // point to the start of the 2nd row

for (j=1; j<Hminus1; j++) {        // *** j-th row

	q++;                           // jump to the 2nd cell (1st time)
	for (i=1; i<Wminus2; i++) {

		(q)->rhs = (q)->u;

		q++;    };                 // next cell of the row
	q++;
	q++;     };                    // next row




//  Gauss-Seidel solver
//  еееееееееееееееееее

uiter = 0;

do {

	max = 0;
	uiter++;
	q = Q + Width;		// point to the start of the 2nd row (bottom-left)

	for (j=1; j<Hminus1; j++) {         // *** j-th row

		q++;
		for (i=1; i<Wminus2; i++) {

			un = (q)->u;							    // u at time n.	
			(q)->u = /*un +*/ om_cc * ( (q)->rhs //- cc * un
									+ aa * (
									(q+1)->u + (q-1)->u + (q+Width)->u +
								    (q-Width)->u) );    // u at time n+1.
			und = (q)->u - un;
			if ((und > epsd) || (und < nepsd)) max = 1;   // Error Check

			q++;
		};

		q++;
		q++;
	};


    //  u no_slip BCs
    //  еееееееееееее
    //
    [self u_BC];


} while (max == 1);



//  ######################################################################


//  v  Diffusion
//  ------------

//  computation of rhs  =  ( v(n) i,j * DXDX / visc dt)
//  ееееееееееееееееее
//
// 	i=1,..,imax;   j=1,..,jmax;

q = Q + Width;					   // point to the start of the 2nd row

for (j=1; j<Hminus2; j++) {         // *** j-th row

	q++;                           // jump to the 2nd cell (1st time)
	for (i=1; i<Wminus1; i++) {

		(q)->rhs = (q)->v;

		q++;    };                 // next cell of the row
	q++;     };                    // next row




//  Gauss-Seidel solver
//  еееееееееееееееееее

viter = 0;

do {

	max = 0;
	viter++;
	q = Q + Width;		// point to the start of the 2nd row (bottom-left)

	for (j=1; j<Hminus2; j++) {         // *** j-th row

		q++;
		for (i=1; i<Wminus1; i++) {

			vn = (q)->v;							    // v at time n.
			(q)->v = /*vn +*/ om_cc * ( (q)->rhs //- cc * vn
							+ aa * ((q+1)->v + (q-1)->v + (q+Width)->v +
									(q-Width)->v) );    // v at time n+1.

			vnd = (q)->v - vn;
			if ((vnd > epsd) || (vnd < nepsd)) max = 1;    // Error Check

			q++;
		};

		q++;
	};


    //  v no_slip BCs
    //  еееееееееееее
    //
    [self v_BC];


} while (max == 1);


}

@end
