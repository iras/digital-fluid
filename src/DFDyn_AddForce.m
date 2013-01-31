//
//  DFDyn_AddForce.m
//  DigitalFluid
//
//  Created by G5 G5 on Sun Jun 05 2005.
//  Copyright (c) 2005 Ivano Ras. All rights reserved.
//

#import "DFDyn_AddForce.h"
#import "DFDyn_constants.h"
#import "DFDyn_Set_BC.h"


//         +---+---+---+---+---+---+
//         | • | • | • | • | • | • |  ----> [i,jmax+1]   last row
//         +---+---+---+---+---+---+
//         | • |   |   |   |   | • |  ----> [i,jmax]     last-but-one row
//         +---+---+---+---+---+---+
//         | • |   |   |   |   | • | 
//         +---+---+---+---+---+---+
//         | • |   |   |   |   | • |
//		   +---+---+---+---+---+---+
//		   | • |   |   |   |   | • |  ----> [i,1]   2nd row
//         +---+---+---+---+---+---+
//         | • | • | • | • | • | • |  ----> [i,0]   1st row
//         +---+---+---+---+---+---+
//
//           |   |			 |   |
//          1st 2nd          |   +------> [imax+1,j]   last column
//	         |   |           +----------> [imax,j]     last-but-one column
//           |   |
//           |   +-----> [1,j]   second column
//           +---------> [0,j]   first  column



@implementation DFDynamics (DFDyn_AddForce)


- (void)AddForce
{

struct Fdata       *q;     //point to a Fluid grid cell
unsigned long int  i,j;    //counters

double xM,yM,uM,vM;




//  Add Force
//  •••••••••
//
//      u(n+1),v(n+1):     i=1,..,imax;  j=1,..,jmax;

q = Q + Width;                // point to the start of the 2nd row

for (j=1; j<Hminus1; j++) {    // *** j-th row

	q++;					   // jump to the 2nd cell of the row : [1,1]

	for (i=1; i<Wminus1; i++) {

			xM = (double)(xMouse/5);
			yM = (double)(yMouse/5);
			uM = (double)(uMouse/25);
			vM = (double)(vMouse/25);

		    (q)->u += uM*exp(-0.035*((i+0.5-xM)*(i+0.5-xM) + (j-yM)*(j-yM)));
			(q)->v += vM*exp(-0.035*((i-xM)*(i-xM) + (j+0.5-yM)*(j+0.5-yM)));

/*
			// located Pulse of Force
			if (t < 2.5) {

			(q)->u += 0;    //dtgx;
			(q)->v += exp(-0.03*((i-20.3)*(i-20.3) + (j-31.4)*(j-31.4)));
            }
*/

		q++;     };                // next cell of the row
	q++;    };					   // next row


/*

//  located horiz impulse of force
//  ••••••••••••••••••••••••••••••

q = Q + (Width/4) + Height*(Height/4);

if (t < 2.1) {

	for (i=0; i<10; i++) {
		(q + i)->u  =  0;
		(q + i)->v  =  50;
	}
}

*/

//  u,v no_slip Boundary Conditions
//  •••••••••••••••••••••••••••••••
//
[self u_BC];
[self v_BC];


}


@end
