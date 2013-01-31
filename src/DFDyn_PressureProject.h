//
//  DFDyn_PressureProject.h
//  DigitalFluid
//
//  Created by G5 G5 on Sun Jun 05 2005.
//  Copyright (c) 2005 Ivano Ras. All rights reserved.
//

#import "DFDynamics.h"


@interface DFDynamics (DFDyn_PressureProject)

- (void)PressureProject_CG;
- (long int)Assemble_Matrix;
- (void)Assemble_b;

//  ------O-N-L-Y---2-D---V-E-R-S-I-O-N-------------------------

- (void)PressureProject;

@end
