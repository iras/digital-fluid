//
//  DFDyn_Advect.m
//  DigitalFluid
//
//  Created by G5 G5 on Sun Jun 05 2005.
//  Copyright (c) 2005 Ivano Ras. All rights reserved.
//

#import "DFDyn_Advect.h"
#import "DFDyn_constants.h"
#import "DFDyn_Set_BC.h"


@implementation DFDynamics (DFDyn_Advect)




- (void)Advect3
{
//  Advect  (Semi-Lagrangian Advection)
//  ••••••
//
//      u(n+1), v(n+1), w(n+1)     for each fluid cell;

struct Fdata             *q;   // point to a Fluid grid cell
unsigned long int     i,j,k;   // counters

// TraceBackward section
double	a, b, c, a1, b1, c1;   // float-coords of the cell to advect.
double  ai, bi, ci;            // intermediate float-coords
double  idb, jdb, kdb;         // i,j,k "doubled"
double  u_avg, v_avg, w_avg;   // estimations on the edge of a cell.
double     uip1, vip1, wip1;   // intermediate velocities



q = Q;                         // point to the start of the GRID

for (k=0; k<Depth; k++) {
    for (j=0; j<Height; j++) {
        for (i=0; i<Width; i++) {

            // It's implied, the container is coated with solid cells.

            if ( (q)->cell_type > 1 ) {      // fluid_cell Checking
                
                idb = (double)(i);
                jdb = (double)(j);
                kdb = (double)(k);

        //  TraceBackward
        //  •••••••••••••
        //
        //  this tracing must be done THREE TIMES for each cell,
        //  because u,v and w have different positions, in the cell !
        //
        //               .--------- i,j+Dy/2,k
        //               :
        //            +==:=====+ ,---(w)---- i,j,k+Dz/2
        //          /   (v)  / |/
        //        +========+   |
        //        |      • |(u)------- i+Dx/2,j,k
        //        |   i,j,k|   +
        //        |        | /
        //        +========+




        //  1st tracing:  u value
        //
        //  Estimate (Linear but it'd be better Cubic Monotonic) of
        //  the vertical and trasversal component at i+0.5 position,
        //  where u is placed. That is the Average, respectively, of
        //  the 4 adjacent vertical v-values and of the 4 planar w-values.
        //
        //      +---(v)---+---(v)---+
        //      |         |         |           ___(w)___ ___(w)___
        //      |    •   [u]        |         /    •   [u]        /
        //      |         |         |       /___(w)___/___(w)___/
        //      +---(v)---+---(v)---+
        //
        //  After this, proceed with the backward RK2 tracing.
        //                                       (RK2 = Runge-Kutta order 2)
        //  e.g., (really it's tridimensional)
        //                   +---------+---------+
        //                   |         |         |
        //                   |    •  ,[u]        |
        //                   |      /  |         |     Pay attention:
        //                   +-----|---+---------+     In RK2 approximation
        //                   |     |   |         |     the curve is a
        //                   |      \  |         |     2-straight-line curve.
        //                   |       ''--(U)     |
        //                   +---------+---------+


        //  RK2 Particle Tracer  (u value)
        //  •••••••••••••••••••
        //
        //  preliminary (linear) tracing     C H A N G E  the average !!!
        //                                   ----------------------------
        v_avg = 0.25 * ( (q)->v + (q+1)->v + (q-Width)->v + (q-Width+1)->v);
        w_avg = 0.25 * ( (q)->w + (q+1)->w + (q-WH)->w + (q-WH+1)->w );
        ai =   - dt * (q)->u;
        bi =   - dt * v_avg;
        ci =   - dt * w_avg;


        //  (intermediate point)
        a = idb + 0.5 + ai*0.5;
        b = jdb       + bi*0.5;
        c = kdb       + ci*0.5;
        // intermediate path tracking
        stop_go_forward = 0;       // reset
        [self LineTracking_x0:(idb+0.5)  y0:jdb  z0:kdb  x1:a  y1:b z1:c];

        if (stop_go_forward == 1) {
            // do not continue the RK2,
            // just interpolate where it is arrived (al,bl,cl)  !
        };


        a = idb + 0.5 + ai;
        b = jdb       + bi;
        c = kdb       + ci;

        //  intermediate velocities
        [self Interp3_U_x:a y:b z:c];
        [self Interp3_V_x:a y:b z:c];
        [self Interp3_W_x:a y:b z:c];
        uip1 = uI;
        vip1 = vI;
        wip1 = wI;

        //  RK2 tracing
        a1 = idb + 0.5 - halfdt * ((q)->u + uip1);
        b1 = jdb       - halfdt * ( v_avg + vip1);
        c1 = kdb       - halfdt * ( w_avg + wip1);

        //  final velocity (order 2)
        [self Interp3_U_x:a1 y:b1 z:c1];
        (q)->u1 = uI;



        //  2nd tracing:  v-value
        //
        //  Estimate (Linear but it'd be better Cubic Monotonic) of
        //  the horizontal and trasversal component at j+0.5 position,
        //  where v is placed. That is the Average, respectively, of
        //  the 4 adjacent horizontal u-values and of the 4 planar w-values.
        //                                    +
        //       +---------+               /  |
        //       |         |            +    (w)
        //      (u)       (u)           |     |
        //       |         |           (w)    +
        //       +---[v]---+            | [v] |
        //       |         |            +    (w)
        //      (u)   •   (u)           |  •  |
        //       |         |           (w)    +
        //       +---------+            |  /
        //                              +
        //
        //  After this, proceed with the backward RK2 tracing.
        //
        //  e.g., (really it's tridimensional)
        //
        //              +---[v]---+---------+
        //              |  /      |         |
        //              | |  •    |         |
        //              | |       |         |
        //              +-|-------+---------+
        //              |  \      |         |
        //              |   \     |         |
        //              |    '(V) |         |
        //              +---------+---------+


        //  RK2 Particle Tracer  (v value)
        //  •••••••••••••••••••
        //
        //  preliminary (linear) tracing
        u_avg = 0.25 * ((q-1)->u + (q)->u+ (q+Width-1)->u + (q+Width)->u);
        w_avg = 0.25 * ((q-WH)->w+ (q)->w+ (q+Width-WH)->w+ (q+Width)->w);
        ai = - dt * u_avg;
        bi = - dt * (q)->v;
        ci = - dt * w_avg;
        
        //  (intermediate point)
        a = idb       + ai*0.5;
        b = jdb + 0.5 + bi*0.5;
        c = kdb       + ci*0.5;
        // ---------------------> ready for the intermediate path tracking
        
        a = idb       + ai;
        b = jdb + 0.5 + bi;
        c = kdb       + ci;

        //  intermediate velocities
        [self Interp3_U_x:a y:b z:c];
        [self Interp3_V_x:a y:b z:c];
        [self Interp3_W_x:a y:b z:c];
        uip1 = uI;
        vip1 = vI;
        wip1 = wI;

        //  RK2 tracing
        a1 = idb       - halfdt * ( u_avg + uip1);
        b1 = jdb + 0.5 - halfdt * ((q)->v + vip1);
        c1 = kdb       - halfdt * ( w_avg + wip1);

        //  final velocity (order 2)
        [self Interp3_V_x:a1 y:b1 z:c1];
        (q)->v1 = vI;



        //  3rd tracing:  w-value
        //
        //  Estimate (Linear but it'd be better Cubic Monotonic) of
        //  the horizontal and the vertical component at k+0.5 position,
        //  where w is placed. That is the Average, respectively, of the
        //  4 adjacent horizontal u-values and of the 4 vertical v-values.
        //
        //                                            +
        //                                        (v) |
        //               _________              +     |
        //            (u)      (u)          (v) |     |
        //           /___[w]___/          +    [w]    +
        //        (u)   •   (u)           |  •  | (v)
        //       /_________/              |     +
        //                                | (v)
        //                                +
        //
        //  After this, proceed with the backward RK2 tracing.


        
        //  RK2 Particle Tracer  (w value)
        //  •••••••••••••••••••
        //
        //  preliminary (linear) tracing
        u_avg = 0.25 * ((q-1)->u + (q)->u + (q+WH-1)->u + (q+WH)->u);
        v_avg = 0.25 * ((q-Width)->v+ (q)->v+ (q-Width+WH)->v+ (q+WH)->v);
        ai = - dt * u_avg;
        bi = - dt * v_avg;
        ci = - dt * (q)->w;
        
        //  (intermediate point)
        a = idb       + ai*0.5;
        b = jdb       + bi*0.5;
        c = kdb + 0.5 + ci*0.5;
        // ---------------------> ready for the intermediate path tracking
        
        a = idb       + ai;
        b = jdb       + bi;
        c = kdb + 0.5 + ci;

        //  intermediate velocities
        [self Interp3_U_x:a y:b z:c];
        [self Interp3_V_x:a y:b z:c];
        [self Interp3_W_x:a y:b z:c];
        uip1 = uI;
        vip1 = vI;
        wip1 = wI;

        //  RK2 tracing
        a1 = idb       - halfdt * ( u_avg + uip1);
        b1 = jdb       - halfdt * ( v_avg + vip1);
        c1 = kdb + 0.5 - halfdt * ((q)->w + wip1);

        //  final velocity (order 2)
        [self Interp3_W_x:a1 y:b1 z:c1];
        (q)->w1 = wI;


            };

            q++;  //  next grid cell

        };
    };
};




//  copy  u1,v1,w1  in  u,v,w
//  •••••••••••••••••••••••••••••

q = Q;                      // point to the start of the GRID

for (k=0; k<Depth; k++) {
    for (j=0; j<Height; j++) {
        for (i=0; i<Width; i++) {
            
            if ( (q)->cell_type > 1 ) {      // fluid_cell Checking

                (q)->u = (q)->u1;
                (q)->v = (q)->v1;
                (q)->w = (q)->w1;

            };
            
            q++;  //  next grid cell
            
        };
    };
};




//  u,v,w  no_slip BCs  +  u,v,w  Divergence_free Conditions
//  ••••••••••••••••••••••••••••••••••••••••••••••••••••••••

[self uvw_BCs_1];
[self uvw_DivFreeCs];
[self uvw_BCs_2];


}







//  This method track a straight line and check the
//  possible obstacles or empty cell on that path.

- (void)LineTracking_x0: (double) a0 y0: (double) b0 z0: (double) c0
                     x1: (double) a1 y1: (double) b1 z1: (double) c1
{
    
struct Fdata         *q;  // point to a Fluid grid cell

double x, y, z;           //  x=x(t);  y=y(t);  z=z(t);
double xprv, yprv, zprv;  //  previous good values
double dx, dy, dz;        //  differences
double increment;         //  increment to the next sample
double tt;                //  parameter to move along the line: (0,1]
unsigned long int xint, yint, zint;


//  init
tt = 0;                   //  parameter on the starting point (a0,b0,c0)
x = a0;  xprv = x;
y = b0;  yprv = y;
z = c0;  zprv = z;

//  compute the differences
dx = a1 - a0;
dy = b1 - b0;
dz = c1 - c0;

//  compute the increment  ( less the numeric factor, more accurate
//                           and more expensive is the estimation )
//  (a better increment, but more slow, is by square root (dx^2+dy^2+dz^2))
increment = 0.3 / ( abs(dx) + abs(dy) + abs(dz) );

//  MAIN LOOP
//
while (tt < 1) {
    
    //  if  x,y,z >= 0  the following casting works like the Floor function
    xint = (unsigned long int)(x);
    yint = (unsigned long int)(y);
    zint = (unsigned long int)(z);
    q = Q + xint + Width*yint + WH*zint;
    
    if ( ((q)->cell_type < 2) ||
         (x < 0) || (x > Wminus1) ||
         (y < 0) || (y > Hminus1) ||
         (z < 0) || (z > Dminus1) )  {  tt = 1;   // exit by the loop
                                        stop_go_forward = 1;  // flag
    } else {

        //  save the last good values.
        xprv = x;  yprv = y;  zprv = z;
        
        //  compute the new values by the parametric straight-line eq.
        x = dx*tt + a0;
        y = dy*tt + b0;
        z = dz*tt + c0;
        tt += increment;
    };

};

// export the last good values !
al = xprv;
bl = yprv;
cl = zprv;

}









- (void)Interp3_U_x:(double)a y:(double)b z:(double)c
{
}
- (void)Interp3_V_x:(double)a y:(double)b z:(double)c
{
}
- (void)Interp3_W_x:(double)a y:(double)b z:(double)c
{
}






















//  ------O-N-L-Y---2-D---V-E-R-S-I-O-N-------------------------
//
- (void)Advect
{

//  Advect  (Semi-Lagrangian Advection)
//  ••••••
//
//      u(n+1), v(n+1):     i=1,..,imax;  j=1,..,jmax;



struct Fdata       *q;        // point to a Fluid grid cell
unsigned long int  i,j;       // counters

// TraceBackward section
double			   a, b;      // floating-coords of the cell to advect.
double         v_avg, u_avg;  // estimations on the edge of a cell.
double          uip1, vip1;   // intermediate velocities



q  = Q + Width;               // point to the start of the 2nd row

for (j=1; j<Hminus1; j++) {   // *** j-th row

	q++;					  // jump to the 2nd cell of the row: [1,1]

	for (i=1; i<Wminus1; i++) {


		//  TraceBackward
		//  •••••••••••••
        //
		//  this tracing must be done TWICE for each cell, because u and v
		//  have different positions, in the cell !
		//
		//        i,j+Dy/2
		//    +---[v]---+
		//    |         |
		//    |    •   [u]  i+Dx/2,j
		//    |   i,j   |
	    //    +---------+




		//  1st tracing:  u value
		//
		//  Estimate (Linear but it would be better Cubic Monotonic) of
		//  the vertical component at i+0.5 position, where is placed u.
		//  That is the Average of the 4 adiacent vertical values.
		//
		//    +---(v)---+---(v)---+
		//    |         |         |
		//    |    •   [u]        |
	    //    |         |         |
		//    +---(v)---+---(v)---+
		//
		//  After this, proceed with the backward RK2 tracing.
		//  e.g.:                        (RK2 is for Runge-Kutta order 2)
		//
		//    +---------+---------+
		//    |         |         |
		//    |    •  ,[u]        |
	    //    |      /  |         |
        //    +-----|---+---------+
	    //    |     |   |         |
	    //    |      \  |         |
	    //    |       ''--(U)     |
	    //    +---------+---------+



		//  RK2 Particle Tracer  (u value)
		//  •••••••••••••••••••
		//
		//  preliminary (linear) tracing
		v_avg = 0.25 * ( (q)->v + (q+1)->v + (q-Width)->v + (q-Width+1)->v);
		a = (double)(i) + 0.5 - dt * (q)->u;
		b = (double)(j)       - dt * v_avg;

		//  intermediate velocities
        [self Interp_U_x:a y:b];
        [self Interp_V_x:a y:b];
        uip1 = uI;
        vip1 = vI;

		//  RK2 tracing
		a = (double)(i) + 0.5 - halfdt * ((q)->u + uip1);
		b = (double)(j)       - halfdt * ( v_avg + vip1);

		//  final velocity (order 2)
        [self Interp_U_x:a y:b];
		(q)->u1 = uI;



		//  2nd tracing:  v-value
		//
		//  Estimate (Linear but it would be better Cubic Monotonic) of
		//  the horizontal component at j+0.5 position, where is placed v.
		//  That is the Average of the 4 adiacent horizontal values.
		//
		//    +---------+
		//    |         |
		//   (u)       (u)
		//    |         |
		//    +---[v]---+
		//    |         |
		//   (u)   •   (u)
		//    |         |
        //    +---------+
		//
		//  After this, proceed with the backward RK2 tracing.
		//  e.g.:
		//
		//    +---[v]---+---------+
		//    |  /      |         |
		//    | |  •    |         |
	    //    | |       |         |
        //    +-|-------+---------+
	    //    |  \      |         |
	    //    |   \     |         |
	    //    |    '(V) |         |
        //    +---------+---------+



		//  RK2 Particle Tracer  (v value)
		//  •••••••••••••••••••
		//
		//  preliminary (linear) tracing
		u_avg = 0.25 * ( (q-1)->u + (q)->u + (q+Width-1)->u + (q+Width)->u);
		a = (double)(i)       - dt * u_avg;
		b = (double)(j) + 0.5 - dt * (q)->v;

		//  intermediate velocities
        [self Interp_U_x:a y:b];
        [self Interp_V_x:a y:b];
        uip1 = uI;
        vip1 = vI;

		//  RK2 tracing
		a = (double)(i)       - halfdt * ( u_avg + uip1);
		b = (double)(j) + 0.5 - halfdt * ((q)->v + vip1);

		//  final velocity (order 2)
        [self Interp_V_x:a y:b];
		(q)->v1 = vI;



		q++;     };            // next cell of the row
	q++;    };				   // next row




//  transfer  u1,v1  in  u,v
//  ••••••••••••••••••••••••

q  = Q + Width;               // point to the start of the 2nd row

for (j=1; j<Hminus1; j++) {    // *** j-th row

	q++;					  // jump to the 2nd cell of the row: [1,1]

	for (i=1; i<Wminus1; i++) {

		(q)->u = (q)->u1;
		(q)->v = (q)->v1;

		q++;      };		  // next cell of the row
	q++;     };				  // next row




//  u,v no_slip Boundary Conditions
//  •••••••••••••••••••••••••••••••

[self u_BC];
[self v_BC];

}











- (void)Interp_U_x:(double)a y:(double)b
{

struct Fdata       *f;            // point to a Fluid grid cell
unsigned long int  ci,cj;         // integer-coords of the cell to advect.
double			   xdist,ydist;          // see Figure.
double			   onemxdist,onemydist;  // auxiliary


//  Interpolate_U
//  •••••••••••••
//
//  the sign 'x' represent the a,b coords passed to this method.
//  Floor_function of a,b return the integer pair ci,cj.
//
//                +---------+---------+
//                |         |         |
//                |    •    |    •    |
//              ^ |         |  x      |
//        ydist | +---------+---------+
//              | |         |         |
//              |_|    o    |    •    |
//                |  ci,cj  |         |
//                +---------+---------+
//                     |------->
//                        xdist


//  Check it don't advect outside the grid.
if (a < 1.0)     { a = 1.0; }
if (b < 1.0)     { b = 1.0; }
if (a > Wminus1) { a = Wminus1; }
if (b > Hminus1) { b = Hminus1; }

//  for (a,b >= 0) the following casting works like the Floor_function.
ci = (unsigned long int)(a);
cj = (unsigned long int)(b);

//  xdist,ydist belongs to the interval [0,1).
xdist = a - (double)(ci);
ydist = b - (double)(cj);

//  If xdist >= 0.5 then it will be used the left approximation
//  to find the interpolation of u, else will be used the right one.
//
//       +---------+---------+           +---------+---------+
//       |         |         |           |         |         |
//       |    •   (•)   •   (•)         (•)   •   (•)   •    |
//       |         |  x      |           |         |         |
//       +---------+---------+           +---------+---------+
//       |         |         |           |       x |         |
//       |    o   (•)   •   (•)         (•)   •   (•)   •    |
//       |         |         |           |         |         |
//       +---------+---------+           +---------+---------+
//            |------->                       |-->
//           xdist >= 0.5                    xdist < 0.5

if (xdist >= 0.5) {	 xdist -= 0.5;         }
else              {  xdist += 0.5;  ci--;  }

onemxdist = 1 - xdist;
onemydist = 1 - ydist;

//  point to the (ci,cj) cell
f  = Q + ci + cj*Width;

//  Weighted Average to obtain the value of u
//  at the sign 'x', to copy in the cell (i,j).
uI = (     xdist * (ydist*(f+1+Width)->u + onemydist*(f+1)->u) + 
	   onemxdist * (ydist*(f+Width)->u   + onemydist*(f)->u)  );

}







- (void)Interp_V_x:(double)a y:(double)b
{

struct Fdata       *f;            // point to a Fluid grid cell
unsigned long int  ci,cj;         // integer-coords of the cell to advect.
double			   xdist,ydist;          // see Figure.
double			   onemxdist,onemydist;  // auxiliary


//  Interpolate_V
//  •••••••••••••
//
//  the sign 'x' represent the a,b coords passed to this method.
//  Floor_function of a,b return the integer pair ci,cj.
//
//                +---------+---------+
//                |         |         |
//                |    •    |    •    |
//              ^ |         |  x      |
//        ydist | +---------+---------+
//              | |         |         |
//              |_|    o    |    •    |
//                |  ci,cj  |         |
//                +---------+---------+
//                     |------->
//                       xdist


//  Check it don't advect outside the grid.
if (a < 1.0)     { a = 1.0; }
if (b < 1.0)     { b = 1.0; }
if (a > Wminus1) { a = Wminus1; }
if (b > Hminus1) { b = Hminus1; }

//  for (a,b >= 0) the following casting works like the Floor_function.
ci = (unsigned long int)(a);
cj = (unsigned long int)(b);

//  xdist,ydist belongs to the interval [0,1).
xdist = a - (double)(ci);
ydist = b - (double)(cj);

//  If ydist >= 0.5 then it will be used the left approximation
//  to find the interpolation of v, else will be used the right one.
//
//            +---(•)---+---(•)---+            +---------+---------+
//            |         |         |            |         |         |
//            |    •    |    •    |            |    •    |    •    |
//          ^ |         |  x      |            |         |         |
//   ydist  | +---(•)---+---(•)---+    ydist   +---(•)---+---(•)---+
//   >= 0.5 | |         |         |    < 0.5 ^ |       x |         |
//          |_|    o    |    •    |          |_|    •    |    •    |
//            |         |         |            |         |         |
//            +---------+---------+            +---(•)---+---(•)---+


if (ydist >= 0.5) {	 ydist -= 0.5;         }
else              {  ydist += 0.5;  cj--;  }

onemxdist = 1 - xdist;
onemydist = 1 - ydist;

//  point to the (ci,cj) cell
f  = Q + ci + cj*Width;

//  Weighted Average to obtain the value of v
//  at the sign 'x', to copy in the cell (i,j).
vI = (     ydist * (xdist*(f+1+Width)->v + onemxdist*(f+Width)->v) + 
       onemydist * (xdist*(f+1)->v       + onemxdist*(f)->v)  );

}



@end