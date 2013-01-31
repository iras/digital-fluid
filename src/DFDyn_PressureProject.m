//
//  DFDyn_PressureProject.m
//  DigitalFluid
//
//  Created by G5 G5 on Sun Jun 05 2005.
//  Copyright (c) 2005 Ivano Ras. All rights reserved.


#import "DFDyn_PressureProject.h"
#import "DFDyn_constants.h"
#import "DFDyn_Set_BC.h"


@implementation DFDynamics (DFDyn_PressureProject)


//  Pressure Projection  (Conservation of Mass)
//  -------------------
//  ending the Hodge-decomposition
//
//    __2      __
//    \/ p  =  \/ U  =  ( du/dx + dv/dy )     scalar Poisson equation
//							divergence
//
//  the numeric form is:
//
//    p i+1,j + p i-1,j + p i,j+1 + p i,j-1 - 4 p i,j   = 
//    DeltaX/2 (u i+1,j - u i-1,j  +  v i,j+1 - v i,j-1)

//  --------------------------------------------------------------------

//  The above p in reality is NOT the real pressure, but
//  this p is equal to (pressure * delta_t / density).
//  Calculate the exact pressure is useless because later be needed
//  re-multiply the 2 constants to make the subtraction of the gradient.

//  The above equation derive from the following reasoning:
//  (Stam 99 article  +  Carlson 04 PhD Thesis)
//                      __
//   (w4 - w3) / dt = - \/pr / ro;     (operator splitting)
//                  __
//   w4 = w3 - dt * \/pr / ro;     next apply the gradient to both sides 
//       __          __2
//   0 = \/w3 - dt * \/pr / ro;  lhs=0 because w4 must be divergence-free
//
//  where with    U = w3    and    p = dt*pressure/ro,
//  we have the initial (Poisson) equation (above).

//  --------------------------------------------------------------------






//  Assemble_Matrix
//  ---------------
//  extract the fluid-filled cells from the CONTAINER
//  according to the   i-j-k order  (see below) to
//  build the matrix for the Pressure projection.

//  The CONTAINER contain additionally the border wall-cells
//  and the possible internal obstacle wall-cells.

//  N = number of such fluid-filled cells that correspond to
//      the side of the matrix.

//  diag   = value of the diagonal entry, equal to:  6 - # of wall_faces.
//           (Neumann conditions: dp/dn = 0)

//  numOff = the number of off diagonal (-1) terms:  6 - # of empty_faces.

//  offIndices[6] = array index of the off diagonal terms.
//                  for each empty-face no one "-1"s must be
//                  inserted in the row.
//                  (Dirichlet conditions: p = 0 for empty-cell)

//  e.g. : 2x2 CONTAINER
//  ---
//                     +=======+=======+
//                   / :  6  / :  7  / |
//                 +=======+=======+   |          j (2nd)
//               / :     / :     / |  -+             ^
//        4    +=======+=======+   | / |             |   k (3rd)
//        :....|     2 |     3 |  -+ 5 |             | /
//             |.'     |.'     | / |  -+             +------>
//             +=======+=======+   | /                        i (1st)
//             |     0 |     1 |  -+
//             |.'     |.'     | /
//             +=======+=======+


- (long int)Assemble_Matrix
{

long int            N;          //  number of the fluid cells (or rather
                                //      the length of all pressure_vectors)
unsigned long int   i,j,k;      //  GRID counters
struct MatrixRow    *r;         //  point to a MatrixRow struct
struct Fdata        *q;         //  point to a Fluid grid cell

int     temp;        //  main_scan temp variable
int     index;       //  main_scan temp variable


//  preliminary scan
//  ----------------
//  it allows to obtain the # of the fluid_cells (to
//  allocate the memory) and to number each fluid_cell.
//  Recall: cell_type 0 & 1  are respectively empty and solid cells.

N = 0;     //  reset of the counter

for (i=0; i<GridCells; i++) {
	if ( (Q+i)->cell_type > 1 ) {
                                   N++;
                                   (Q+i)->number = N;
	};
};

//  Memory Allocation
//  -----------------
//  Start pointer in Arow.
//    •   sizeof(double) = 8 bytes:   x, b, pp, qq, rr
//    ••  sizeof(int)    = 4 bytes:   diag, numOff, offIndices[6]
Arow  =  malloc ( N * (sizeof(int) * 5 + sizeof(int) * 8) );



//  Main Scan
//  ---------
//  it allows to assemble the matrix
r = Arow;                   // point to the start of the 1st row
q = Q;                      // point to the start of the GRID

for (k=0; k<Depth; k++) {
   for (j=0; j<Height; j++) {
      for (i=0; i<Width; i++) {

	      if ( (q)->cell_type > 1 ) {      // fluid_cell Checking

              (r)->x = (q)->p;   // old pressure as initial condition in CG

              index = 0;         // index reset
              (r)->diag   = 6;
              (r)->numOff = 6;

              //  NOW, IDENTIFY EACH OF THE SIX ADIACENT CELL

              //   ==========================(  q+1  )==============
              if ( i < (Width-1) ) {

                  temp = (q+1)->cell_type;
                  if ( temp == 0 ) {  (r)->numOff--;  };   //  empty-cell
                  if ( temp == 1 ) {  (r)->diag--;    };   //  solid-cell
                  if ( temp  > 1 ) {                       //  fluid-cell
                      (r)->offIndices[index] = (q+1)->number;
                      index++;
                  };
               } else {
                   // pretend to find a solid-cell (temporary measure)
                   (r)->diag--;
               };


              //   ==========================(  q-1  )==============
              if ( i > 0 ) {

                  temp = (q-1)->cell_type;
                  if ( temp == 0 ) {  (r)->numOff--;  };   //  empty-cell
                  if ( temp == 1 ) {  (r)->diag--;    };   //  solid-cell
                  if ( temp  > 1 ) {                       //  fluid-cell
                      (r)->offIndices[index] = (q-1)->number;
                      index++;
                  };
              } else {
                  // pretend to find a solid-cell (temporary measure)
                  (r)->diag--;
              };


              //   ==========================(  q+Width  )==============
              if ( j < (Height-1) ) {

                  temp = (q+Width)->cell_type;
                  if ( temp == 0 ) {  (r)->numOff--;  };   //  empty-cell
                  if ( temp == 1 ) {  (r)->diag--;    };   //  solid-cell
                  if ( temp  > 1 ) {                       //  fluid-cell
                      (r)->offIndices[index] = (q+Width)->number;
                      index++;
                  };
              } else {
                  // pretend to find a solid-cell (temporary measure)
                  (r)->diag--;
              };


              //   ==========================(  q-Width  )==============
              if ( j > 0 ) {

                  temp = (q-Width)->cell_type;
                  if ( temp == 0 ) {  (r)->numOff--;  };   //  empty-cell
                  if ( temp == 1 ) {  (r)->diag--;    };   //  solid-cell
                  if ( temp  > 1 ) {                       //  fluid-cell
                      (r)->offIndices[index] = (q-Width)->number;
                      index++;
                  };
              } else {
                  // pretend to find a solid-cell (temporary measure)
                  (r)->diag--;
              };


              //   ==========================(  q+WH  )==============
              if ( k < (Depth-1) ) {

                  temp = (q+WH)->cell_type;
                  if ( temp == 0 ) {  (r)->numOff--;  };   //  empty-cell
                  if ( temp == 1 ) {  (r)->diag--;    };   //  solid-cell
                  if ( temp  > 1 ) {                       //  fluid-cell
                      (r)->offIndices[index] = (q+WH)->number;
                      index++;
                  };
              } else {
                  // pretend to find a solid-cell (temporary measure)
                  (r)->diag--;
              };


              //   ==========================(  q-WH  )==============
              if ( k > 0 ) {

                  temp = (q-WH)->cell_type;
                  if ( temp == 0 ) {  (r)->numOff--;  };   //  empty-cell
                  if ( temp == 1 ) {  (r)->diag--;    };   //  solid-cell
                  if ( temp  > 1 ) {                       //  fluid-cell
                      (r)->offIndices[index] = (q-WH)->number;
                      index++;
                  };
              } else {
                  // pretend to find a solid-cell (temporary measure)
                  (r)->diag--;
              };


              r++;   //  next matrix cell

          };

          q++;  //  next grid cell

      };
   };
};

return(N);   //  return the number of the fluid cells (or rather
             //                the length of all pressure_vectors).
}








//  fill the b_vector with the known values.
- (void)Assemble_b
{
//   ••••••••••••••••••••••••••••••••••
//   ••••••••••••••••••••••••••••••••••
//   ••••••••••••••••••••••••••••••••••
//   ••••••••••••••••••••••••••••••••••
}









//  Recall:   The method of conjugate gradients (CG) is an iterative
//  technique 4 solving symmetric positive-definite linear systems (Ax=b).

- (void)PressureProject_CG
{

long int   N;              //  number of the fluid cells (or rather
                           //     the length of all pressure_vectors)
long int   c1,c2;          //  counters.
long int   i,j,k;          //  GRID counters.
struct   MatrixRow    *r;  //  point to a Matrix Row struct.
struct   Fdata        *q;  //  point to a Fluid grid cell.
//  scalars
double   alphaCG, betaCG;
double   rho, rho_old, bnorm;

//  The epsCG, itermax are constants  (see DFDyn_constants.h)


//  MatrixRow variables resume :
// ----------------------------------------+
//  A is an vector of type MatrixRow, and  |
//  include the following variables:       |
//                                         |
//  b is for the vector of known values.   |
//                                         |
//  The x vector initially holds a guess;  |
//  each loop x is updated with the new    |
//  iterate, and at the end it holds the   |
//  solution.                              |
//                                         |
//  The rr vector holds the residuals, and |
//  pp and qq store intermediate values.   |
// ----------------------------------------+


//  Assemble the Matrix and return the length N  (of any pressure_vector).
N = [self Assemble_Matrix];
//  fill the b_vector with the known values.
[self Assemble_b];

rho = 0.0;

//  compute: bnorm = √ (b·b)  ( vector_vector multiply )
r = Arow;                      //  point to the start of the 1st row
bnorm = 0;                     //  bnorm reset.
for (c1 = 0; c1 < N; c1++) {
       bnorm += (r)->b * (r)->b;
       r++;
};
bnorm = sqrt(bnorm);

//  compute: rr = b - Ax   ( matrix_vector multiply + vector_vector subtr)
r = Arow;                 //  point to the start of the 1st row
for (c1 = 0; c1 < N; c1++)  {
    (r)->rr = (r)->diag * (r)->x;
    for (c2 = 0; c2 < (r)->numOff; c2++) {
        (r)->rr += - ( Arow + (r)->offIndices[c2] )->x;   //  (-1) * ...
    };
    r++;
};
r = Arow;                 //  point to the start of the 1st row
for (c1 = 0; c1 < N; c1++)  {
    (r)->rr = (r)->b - (r)->rr;
    r++;
};

//  M A I N   L O O P
//  -----------------
for (i = 0; i < itermax; i++)  {

     rho_old = rho;

     //  compute: rho = rr · rr    ( vector_vector multiply )
     r = Arow;                      //  point to the start of the 1st row
     rho = 0;                       //  rho reset.
     for (c1 = 0; c1 < N; c1++) {
         rho += (r)->rr * (r)->rr;
         r++;
     };

     // Check for convergence.
     if ( (rho != 0.0) && (sqrt(rho) > (epsCG*bnorm)) ) {

        if (i == 0) {
                    //  compute: pp = rr  ( vector copying )
                    r = Arow;       //  point to the start of the 1st row
                    for (c1 = 0; c1 < N; c1++) {
                        (r)->pp = (r)->rr;
                        r++;
                    };
        } else {
                    betaCG = rho/rho_old;

                    //  compute: pp *= betaCG  ( constant_vector multiply )
                    r = Arow;       //  point to the start of the 1st row
                    for (c1 = 0; c1 < N; c1++) {
                        (r)->pp *= betaCG;
                        r++;
                    };
                    //  compute: pp += rr      ( vector_vector add )
                    r = Arow;       //  point to the start of the 1st row
                    for (c1 = 0; c1 < N; c1++) {
                        (r)->pp += (r)->rr;
                        r++;
                    };
        }

        //  compute: qq = App   ( matrix_vector multiply )
        r = Arow;                   //  point to the start of the 1st row
        for (c1 = 0; c1 < N; c1++)  {
            (r)->qq = (r)->diag * (r)->pp;
            for (c2 = 0; c2 < (r)->numOff; c2++) {
                (r)->qq += - ( Arow + (r)->offIndices[c2] )->pp;  //(-1)*...
            };
            r++;
        };

        //  compute: alphaCG = rho/(pp·qq)    ( vector_vector multiply )
        r = Arow;                   //  point to the start of the 1st row
        alphaCG = 0;                //  alphaCG reset.
        for (c1 = 0; c1 < N; c1++) {
            alphaCG += (r)->pp * (r)->qq;
            r++;
        };
        alphaCG = rho / alphaCG;

        //  compute: x  += alphaCG*pp   ( constant_vector multiply )
        //  compute: rr -= alphaCG*qq   ( constant_vector multiply )
        r = Arow;                   //  point to the start of the 1st row
        for (c1 = 0; c1 < N; c1++) {
            (r)->x  += alphaCG * (r)->pp;
            (r)->rr -= alphaCG * (r)->qq;
            r++;
        };

     };
};


//  Now CG is ended, x hold the NEW_pressure !
//  The remaining operations are:


//   (i)  replace the NEW_pressure in the CONTAINER cells.
//
r = Arow;                       // point to the start of the 1st row
q = Q;                          // point to the start of the GRID
for (i=0; i<GridCells; i++) {

    if ( (q)->cell_type > 1 ) {      // fluid_cell Checking

        (q)->p = (r)->x;     // NEW_pressure in the CONTAINER
        r++;
    };
    q++;

};


//   (ii)  compute the NEW_velocities for each fluid_cell.
//
//   subtraction of gradient
//   -----------------------
//	                       __
//          Unew  =  U  -  \/ p
//
//  in scalar notation:
//
//          unew  =  u  -  dp/dx;
//		    vnew  =  v  -  dp/dy;
//		    wnew  =  w  -  dp/dw;
//
q = Q;                              // point to the start of the GRID
for (k=0; k<Depth; k++) {
    for (j=0; j<Height; j++) {
        for (i=0; i<Width; i++) {

            if ( (q)->cell_type > 1 ) {    // fluid_cell Checking


                //   ==========================(  i  )==============
                if ( i < (Width-1) ) {
                    //  if next is a solid-cell => set the p-Neumann condition
                    if ( (q+1)->cell_type == 1 )  {  (q+1)->p = (q)->p;  };
                    //  (q)->u must be computed on the right edge of each cell.
                    (q)->u -= c1_DX * ( (q+1)->p  -  (q)->p );

                } else {
                    // pretend to find a solid-cell (temporary measure)
                    (q)->u = 0;
                };


                //   ==========================(  j  )==============
                if ( j < (Height-1) ) {
                    //  if next is a solid-cell => set the p-Neumann condition
                    if ( (q+Width)->cell_type == 1 ) { (q+Width)->p=(q)->p; };
                    //  (q)->v must be computed on the upper edge of each cell.
                    (q)->v -= c1_DY * ( (q+Width)->p  -  (q)->p );

                } else {
                    // pretend to find a solid-cell (temporary measure)
                    (q)->v = 0;
                };


                //   ==========================(  k  )==============
                if ( k < (Depth-1) ) {
                    //  if next is a solid-cell => set the p-Neumann condition
                    if ( (q+WH)->cell_type == 1 )  {  (q+WH)->p = (q)->p;  };
                    //  (q)->w must be computed on the inner edge of each cell.
                    (q)->w -= c1_DZ * ( (q+WH)->p  -  (q)->p );
                } else {
                    // pretend to find a solid-cell (temporary measure)
                    (q)->w = 0;
                };

            };


            if ( (q)->cell_type == 0 ) {    //  empty_cell Checking
                (q)->p = 0;                 //  reset to zero of the Pressure
            };


            q++;   // next CONTAINER cell

        };
    };
};


//   (iii)  set the   u,v,w no_slip BCs   +   u,v,w Div_free Cs
//
//  u,v,w  no_slip Boundary  Conditions
//  -----------------------------------
//  They deal the velocities where a fluid-cell is adjacent to a solid-cell.
//  This velocity-handling concern either the fluid-cell and/or the solid-cell,
//  depending on the case.
//  The  u,v,w no_slip BCs 'uvw_BCs_1' must be set "before" the u,v,w Div_free Cs.
//  The  u,v,w no_slip BCs 'uvw_BCs_2' must be set "after" the u,v,w Div_free Cs.
[self uvw_BCs_1];
//
//  u,v,w  Divergence_free   Conditions
//  -----------------------------------
//  They deals the velocities where a fluid-cell is adjacent to an empty-cell.
//  This velocity-handling concern either the fluid-cell and/or the empty-cell,
//  depending on the case.
[self uvw_DivFreeCs];

[self uvw_BCs_2];



//   (iv)  release the memory previously allocated by 'Assemble_Matrix'.
//
free (Arow);


}






























//  ------O-N-L-Y---2-D---V-E-R-S-I-O-N-------------------------
//
- (void)PressureProject
{
struct Fdata        *q;        // point to a Fluid grid cell
unsigned long int   i,j;       // counters
double              pn, pnd;   // p at step n;   pnd = p(n+1) - p(n);
unsigned long int   max;       // max = boolean flag to exit from SOR loop




//  computation of rhs
//  ••••••••••••••••••
//
//		rhs  =  (u i+1,j - u i-1,j  +  v i,j+1 - v i,j-1) * DX/2
//
// 	i=1,..,imax;   j=1,..,jmax;

q = Q + Width;					   // point to the start of the 2nd row

for (j=1; j<Hminus1; j++) {         // *** j-th row

	q++;                           // jump to the 2nd cell (1st time)
	for (i=1; i<Wminus1; i++) {

		//  (q)->rhs must be computed on the center of each cell
		//  to follow the position of the pressure in the cell.
		//
		(q)->rhs = DeltaX * ( (q)->u - (q-1)->u   +
							  (q)->v - (q-Width)->v );

		/* // centered cell version
		(q)->rhs = DeltaX_2 * ( (q+1)->u   - (q-1)->u   +
								(q+Width)->v - (q-Width)->v );
		*/

		q++;    };                 // next cell of the row
	q++;     };                    // next row




//  SOR solver
//  ••••••••••

iter = 0;

do {

	max = 0;
	iter++;
	q = Q + Width;     // point to the start of the 2nd row (bottom-left)

	for (j=1; j<Hminus1; j++) {         // *** j-th row

		q++;
		for (i=1; i<Wminus1; i++) {

			pn = (q)->p;			                     // p at time n.
			(q)->p = pn + omq * ( (q)->rhs + pn+pn+pn+pn -
								  (q+1)->p - (q-1)->p - (q+Width)->p -
								  (q-Width)->p );        // p at time n+1.
			pnd = (q)->p - pn;
			if ((pnd > eps) || (pnd < neps)) max = 1;    // Error Check

			q++;
		};

		q++;
	};


    //  p-gradient Boundary Conditions
    //  ••••••••••••••••••••••••••••••
    //
    [self p_BC];


} while (max == 1);






//  subtraction of gradient
//  •••••••••••••••••••••••
//	                    __
//       Unew  =  U  -  \/ p
//
//  in scalar notation:
//
//       unew  =  u  -  dp/dx;
//		 vnew  =  v  -  dp/dy;
//
// 	i=1,..,imax;   j=1,..,jmax;

q = Q + Width;                      // point to the start of the 2nd row

for (j=1; j<Hminus1; j++)  {        //  *** j-th row

	q++;                            // jump to the 2nd cell (first time)
	for (i=1; i<Wminus1; i++)  {

		//  (q)->u must be computed on the right edge of each cell.
		(q)->u -= c1_DX * (  (q+1)->p   - (q)->p );
		//  (q)->v must be computed on the upper edge of each cell.
		(q)->v -= c1_DY * ((q+Width)->p - (q)->p );

		/* // centered cell version
		(q)->u -= c1_2DX * (  (q+1)->p   - (q-1)->p );
		(q)->v -= c1_2DY * ((q+Width)->p - (q-Width)->p );
		*/

		q++;     };                 // next cell of the row
	q++;     };                     // next row




//  u,v no_slip Boundary Conditions
//  •••••••••••••••••••••••••••••••

[self u_BC];
[self v_BC];



}

@end
