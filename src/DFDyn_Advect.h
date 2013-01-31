//
//  DFDyn_Advect.h
//  DigitalFluid
//
//  Created by G5 G5 on Sun Jun 05 2005.
//  Copyright (c) 2005 Ivano Ras. All rights reserved.
//

#import "DFDynamics.h"


@interface DFDynamics (DFDyn_Advect)

- (void)Advect3;
- (void)LineTracking_x0: (double) a0 y0: (double) b0 z0: (double) c0
                     x1: (double) a1 y1: (double) b1 z1: (double) c1;
- (void)Interp3_U_x: (double) a y: (double) b z: (double) c;
- (void)Interp3_V_x: (double) a y: (double) b z: (double) c;
- (void)Interp3_W_x: (double) a y: (double) b z: (double) c;

//  ------O-N-L-Y---2-D---V-E-R-S-I-O-N-------------------------

- (void)Advect;
- (void)Interp_U_x: (double) a y: (double) b;
- (void)Interp_V_x: (double) a y: (double) b;


@end
