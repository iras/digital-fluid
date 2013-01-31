//
//  DFDynamics.m
//  DigitalFluid
//
//  Created by G5 G5 on Thu May 05 2005.
//  Copyright (c) 2005 Ivano Ras. All rights reserved.
//

#import "DFDynamics.h"
#import "DFDyn_constants.h"
#import "DFDyn_AddForce.h"
#import "DFDyn_Advect.h"
// #import "DFDyn_Diffuse.h"
#import "DFDyn_PressureProject.h"
#import "DFDyn_Set_BC.h"
#import "DFDyn_UpdateMarkers.h"


@implementation DFDynamics

- (void)instance_vars_init
{
	// iter reset
	iter = 0;

    //-----------------------------------

    // counter
    long int n;

    // Allocate mem_block to store the Fluid grid.
	// Start pointer in Q.
	//   •     sizeof(double) = 8 bytes.
	//   ••    8 = # of "double"s in the struct.
	//   •••   sizeof(int)    = 4 bytes.
    Q  =  malloc (GridCells * ( sizeof(double) * 8 + sizeof(int) * 2 ) );
                                                   
	// reset of the grid values
    for (n=0; n<GridCells; n++) {
		(Q+n)->p = (Q+n)->u  = (Q+n)->v   =
                   (Q+n)->u1 = (Q+n)->v1  =
							   (Q+n)->rhs = 0.0;
	}
	
    //-----------------------------------

	t = 0.0;
}



- (void)dynamics    // Explicit Euler scheme
{


	[self AddForce];
//  [self PressureProject];  // Stam suggest to use still this step
	[self Advect];
	//[self Diffuse];
	[self PressureProject];

}

@end
