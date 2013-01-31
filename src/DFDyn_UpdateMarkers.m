//
//  DFDyn_UpdateMarkers.m
//  DigitalFluid
//
//  Created by iMac on 04/09/05.
//  Copyright 2005 Ivano Ras. All rights reserved.
//

#import "DFDyn_UpdateMarkers.h"


@implementation DFDynamics (DFDyn_UpdateMarkers)     // categories


- (void)Update_Markers
{

//  RK2 if possible, else Extrapolation !

}



- (void)Identify_Cells
{


//  Does the cell contain at least one Marker-particle ?
//  ----------------------------------------------------
//
//     •   NO   :  empty-cell or solid-cell  ->  cell_type = 0 or 1;
//
//     •  YES  :  fluid-filled cell



//  If YES, which cell type ?
//  -------------------------
//
//     •  full-fluid  :  6 faces are shared with others full-fluid cells
//                       (governed by the Navier-Stokes eqs)
//
//     •  fluid-empty :  at least one face is shared with an empty cell.
//        interface      (divergence-free velocity conditions)
//
//     •  fluid-wall  :  at least one face is shared with a wall cell.
//        interface      (BCs : no-slip  &&/||  free-slip conditions)
//
//     •  fluid-      :  at least one face is shared with a wall cell
//        empty-wall     and at least one face is shared with a wall
//        interface      cell.
//                       (both div-free velocity + boundary conditions)


//     •  obstacle    :  empty cell with at least one surface adiacent
//                       to a fluid cell. (wall-cell & obstacle-cell)
//                       (BCs : no-slip  &&/||  free-slip conditions)


}


@end
