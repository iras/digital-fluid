//
//  DFDynamics.h
//  DigitalFluid
//
//  Created by G5 G5 on Thu May 05 2005.
//  Copyright (c) 2005 Ivano Ras. All rights reserved.
//

#import <Foundation/Foundation.h>


@interface DFDynamics : NSObject
{
@public     // each instance-var is accessible everywhere.

	// mouse position and motion-direction
	int xMouse, yMouse;
	int uMouse, vMouse;

    // time elapsed
    double t;

	// SOR-iterations to converge (pressure projection)
	unsigned int  iter;

	// GaussSeidel-iterations to converge (diffuse velocity)
	unsigned int  uiter, viter, witer;

	// Interpolation velocity results
	double  uI, vI, wI;
    
    // Line-tracking positon resuls + flag
    double  al, bl, cl;
    int     stop_go_forward;

    // Marker particle data + array declaration
	struct Mdata {
		double x;
		double y;
		double z;
	}  MP[1000];      // declare MP[ ] as an array of Mdata type

    // Basic staggered cell (Fluid or not) definition + pointer declaration
    struct Fdata     {
		int     cell_type;  // identify the type of the fluid-cell  (see below)
		int     number;     // # of the cell in the seq (PressureProject matrix)
		double  p;          // pressure
		double  u;          // x velocity (step n)
		double  v;          // y velocity (step n)
		double  w;          // z velocity (step n)
		double  u1;         //   x velocity (step n+1)
		double  v1;			//   y velocity (step n+1)
		double  w1;			//   z velocity (step n+1)
		double  rhs;        // right hand side (needed by the iterative solvers)
	} *Q;              // declare Q as a pointer variable to Fdata type

    // Row of the matrix A (+ other stuff) for pressure projectn + pointer decl
    struct MatrixRow  {
        double x,b;         //  variables used in  Ax = b  (CG loop)
        double pp,qq,rr;    //  auxiliary variables + residual (CG loop)
		int    diag;           //  value of the diagonal entry (A row)
		int    numOff;         //  number of off diagonal terms (A row)
		int    offIndices[6];  //  array index of the off diagonal terms (A row)
    } *Arow;           // declare Arow as a pointer variable to MatrixRow type

}

- (void)instance_vars_init;
- (void)dynamics;

@end




//    cell (i,j,k)
//    ============

//    4 values :   p    at the center of the cell.
//               u,v,w  onto the center of the three orthogonal faces.
//
//
//                     +======================+
//                   / :                     /|
//                 /   '                   /  |
//               /           X       ..../....|....  w(i,j,k)  <= hidden face
//             /         v(i,j,k)  .'  /      |
//           /                    x  /        |
//         +=======================+          |
//         |                       |          |
//         |      o          *     |     X •--+----  u(i,j,k)
//         |             p(i,j,k)  |          |
//         |                       |        --+
//         |           O           |         /
//         |                       |       /          j
//         |                 o     |     /             ^
//         |  ,                    |   /               |   k
//         |.'                     | /                 | /
//         +=======================+                   +------->
//                                                              i


//   cell_type
//   ===+=====
//      |
//      +----> 0 : empty
//      |
//      +----> 1 : solid (obstacle/wall)
//      |
//      +----> 2 : fluid







//   OLD guess cell_type Tree
//   ========================
//      |
//      +----> 0 : empty
//      |
//      +----> 1 : solid (obstacle/wall)
//      |
//      |
//      +-- fluid-filled -+-> 2 : full-fluid  (Navier-Stokes eqs)
//                        |
//                        |
//                        +- interface -+-> e*100  : fluid-empty (div-free vel)
//                                      |
//                                      +-> s*1000 : fluid-solid (BCs)
//                                      |
//                                      +-> e*100 + s*1000 : fluid-empty-solid
//                                                           (div-free + BCs)

