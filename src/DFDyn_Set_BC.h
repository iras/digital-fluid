//
//  DFDyn_Set_BC.h
//  DigitalFluid
//
//  Created by G5 G5 on Mon Jun 27 2005.
//  Copyright (c) 2005 Ivano Ras. All rights reserved.
//

#import "DFDynamics.h"


@interface DFDynamics (DFDyn_Set_BC)   // categories

- (void)uvw_BCs_1;
- (void)uvw_BCs_2;
- (void)uvw_DivFreeCs;



//  ------O-N-L-Y---2-D---V-E-R-S-I-O-N-------------------------

- (void)u_BC;
- (void)v_BC;
- (void)p_BC;

@end
