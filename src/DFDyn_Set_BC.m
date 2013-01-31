//
//  DFDyn_Set_BC.m
//  DigitalFluid
//
//  Created by G5 G5 on Mon Jun 27 2005.
//  Copyright (c) 2005 Ivano Ras. All rights reserved.
//

#import "DFDyn_Set_BC.h"
#import "DFDyn_constants.h"


@implementation DFDynamics (DFDyn_Set_BC)    // categories



//  This method keeps in line the velocities according to
//  the Divergence_free condition (continuity eq.) when a
//  fluid_cell is adjacent to, at least, one empty_cell.
//          ==>  63  possible cell configs!
//
- (void)uvw_DivFreeCs
{

struct     Fdata      *q;  //  point to a Fluid grid cell.
long int   i,j,k;          //  GRID counters.
int        empty_counter;  //  number of adjacent empty-faces.
int    re,le,ue,de,ie,oe;  //  empty_flags (see below: stylized_cube)
double     temp;           //  temporary variable.


//                                +=======+  ..... ie = inner-empty
//            up-empty = ue     /   ue  / |.'
//                            +=======+   |
//         left-empty = le ...|       | re+     re = right-empty
//                            |  oe   | /
//       outer-empty = oe     +=======+
//                                  :...... de = down-empty
//
//        re = (q)->u = u i,j,k;     le = (q-1)->u = u i-1,j,k;
//        ue = (q)->v = v i,j,k;     de = (q-Width)->v = v i,j-1,k;
//        ie = (q)->w = w i,j,k;     oe = (q-WH)->w = w i,j,k-1;
//
//                      continuity  equation:
//                      --------------------
//   u i,j,k - u i-1,j,k + v i,j,k - v i,j-1,k + w i,j,k - w i,j,k-1 = 0


//  scan through the fluid_cells to search how many and which faces
//  are adjacent to an empty cell, thus apply the suitable DivFree Cs.
//
q = Q;                              // point to the start of the GRID
for (k=0; k<Depth; k++) {
    for (j=0; j<Height; j++) {
        for (i=0; i<Width; i++) {

            if ( (q)->cell_type > 1 ) {    // cell_type Checking

                empty_counter = 0;         //  reset to zero.
                re=le=ue=de=ie=oe = 0;     //  reset all flags to zero.


                //  empty-face COUNTING and POSITION-DETECTING
                //  ------------------------------------------

                if ( i < (Width-1) ) {  // ============(  q+1  )=====
                    if ( (q+1)->cell_type == 0 ) {   empty_counter++;
                                                     re = 1;    };
                } else {
                    // pretend to find a solid-cell (temporary measure)
                };

                if ( i > 0 ) {  // ====================(  q-1  )=====
                    if ( (q-1)->cell_type == 0 ) {   empty_counter++;
                                                     le = 1;    };
                } else {
                    // pretend to find a solid-cell (temporary measure)
                };

                if ( j < (Height-1) ) {  // =========(  q+Width  )======
                    if ( (q+Width)->cell_type == 0 ) {  empty_counter++;
                                                        ue = 1;   };
                } else {
                    // pretend to find a solid-cell (temporary measure)
                };

                if ( j > 0 ) {  // ==================(  q-Width  )======
                    if ( (q-Width)->cell_type == 0 ) {  empty_counter++;
                                                        de = 1;   };
                } else {
                    // pretend to find a solid-cell (temporary measure)
                };

                if ( k < (Depth-1) ) {  // =============(  q+WH  )=====
                    if ( (q+WH)->cell_type == 0 ) {    empty_counter++;
                                                       ie = 1;   };
                } else {
                    // pretend to find a solid-cell (temporary measure)
                };

                if ( k > 0 ) {  // ====================(  q-WH  )======
                    if ( (q-WH)->cell_type == 0 ) {    empty_counter++;
                                                       oe = 1;   };
                } else {
                    // pretend to find a solid-cell (temporary measure)
                };


                //  CASE STUDY
                //  ----------

                //  one empty-face     (6 configs)
                //  ----> only in this case the solution is unique.
                if (empty_counter == 1) {

                    if (re == 1) {    
                        (q)->u = (q-1)->u - (q)->v + (q-Width)->v -
                                               (q)->w + (q-WH)->w;     };
                    if (le == 1) {
                        (q-1)->u = (q)->u + (q)->v - (q-Width)->v +
                                               (q)->w - (q-WH)->w;     };
                    if (ue == 1) {
                        (q)->v = - (q)->u + (q-1)->u + (q-Width)->v -
                                               (q)->w + (q-WH)->w;     };
                    if (de == 1) {
                        (q-Width)->v = (q)->u - (q-1)->u + (q)->v +
                                               (q)->w - (q-WH)->w;     };
                    if (ie == 1) {
                        (q)->w = - (q)->u + (q-1)->u - (q)->v +
                                         (q-Width)->v + (q-WH)->w;     };
                    if (oe == 1) {
                        (q-WH)->w = (q)->u - (q-1)->u + (q)->v -
                                            (q-Width)->v + (q)->w;     };
                };


                //  two empty-faces    (15 configs -> 3 + 12)
                //  ---
                if (empty_counter == 2) {

                    //  opposite empty-faces (3 sub_configs)
                    if ( ((re == 1) && (le == 1)) ||
                         ((ue == 1) && (de == 1)) ||
                         ((ie == 1) && (oe == 1)) )   {
                        //  we do nothing. It might results not
                        //  div-free, but the tries say OK.

                    //  not opposite empty-faces (12 sub_configs)
                    } else {
                        //  copy the velocities on the empty-face from
                        //  the opposite side face and add (or '-') half
                        //  the difference of the remaining two faces to
                        //  each surface-face, to get a div-free cell.

                        if ((re == 1) && (oe == 1)) {
                                    temp = 0.5 * ((q)->v - (q-Width)->v);
                                    (q)->u    = (q-1)->u - temp;
                                    (q-WH)->w = (q)->w + temp;     };

                        if ((re == 1) && (ie == 1)) {
                                    temp = 0.5 * ((q)->v - (q-Width)->v);
                                    (q)->u = (q-1)->u - temp;
                                    (q)->w = (q-WH)->w - temp;     };

                        if ((re == 1) && (ue == 1)) {
                                    temp = 0.5 * ((q)->w - (q-WH)->w);
                                    (q)->u = (q-1)->u - temp;
                                    (q)->v = (q-Width)->v - temp;   };

                        if ((re == 1) && (de == 1)) {
                                    temp = 0.5 * ((q)->w - (q-WH)->w);
                                    (q)->u       = (q-1)->u - temp;
                                    (q-Width)->v = (q)->v + temp;    };


                        if ((le == 1) && (oe == 1)) {
                                    temp = 0.5 * ((q)->v - (q-Width)->v);
                                    (q-1)->u  = (q)->u + temp;
                                    (q-WH)->w = (q)->w + temp;     };

                        if ((le == 1) && (ie == 1)) {
                                    temp = 0.5 * ((q)->v - (q-Width)->v);
                                    (q-1)->u = (q)->u + temp;
                                    (q)->w   = (q-WH)->w - temp;     };

                        if ((le == 1) && (ue == 1)) {
                                    temp = 0.5 * ((q)->w - (q-WH)->w);
                                    (q-1)->u = (q)->u + temp;
                                    (q)->v   = (q-Width)->v - temp;  };

                        if ((le == 1) && (de == 1)) {
                                    temp = 0.5 * ((q)->w - (q-WH)->w);
                                    (q-1)->u     = (q)->u + temp;
                                    (q-Width)->v = (q)->v + temp;    };


                        if ((ue == 1) && (oe == 1)) {
                                    temp = 0.5 * ((q)->u - (q-1)->u);
                                    (q)->v    = (q-Width)->v - temp;
                                    (q-WH)->w = (q)->w + temp;       };
                        
                        if ((ue == 1) && (ie == 1)) {
                                    temp = 0.5 * ((q)->u - (q-1)->u);
                                    (q)->v = (q-Width)->v - temp;
                                    (q)->w = (q-WH)->w - temp;       };
                        
                        if ((de == 1) && (oe == 1)) {
                                    temp = 0.5 * ((q)->u - (q-1)->u);
                                    (q-Width)->v = (q)->v + temp;
                                    (q-WH)->w = (q)->w + temp;       };
                        
                        if ((de == 1) && (ie == 1)) {
                                    temp = 0.5 * ((q)->u - (q-1)->u);
                                    (q-Width)->v = (q)->v + temp;
                                    (q)->w       = (q-WH)->w - temp; };
                    };
                };

                //  three empty-faces    (20 configs -> 8 + 12)
                //  -----
                if (empty_counter == 3) {

                    //  none opposite empty-face    (8 sub_configs)
                    if ( ((oe == 1) && (re == 1) && (ue == 1)) ||
                         ((oe == 1) && (le == 1) && (ue == 1)) ||
                         ((oe == 1) && (re == 1) && (de == 1)) ||
                         ((oe == 1) && (le == 1) && (de == 1)) ||
                         ((ie == 1) && (re == 1) && (ue == 1)) ||
                         ((ie == 1) && (le == 1) && (ue == 1)) ||
                         ((ie == 1) && (re == 1) && (de == 1)) ||
                         ((ie == 1) && (le == 1) && (de == 1)) ) {
                        //  only copy from the opposite faces.

                        if ((oe == 1) && (re == 1) && (ue == 1)) {
                            (q-WH)->w = (q)->w;
                            (q)->u    = (q-1)->u;
                            (q)->v    = (q-Width)->v;   };
 
                        if ((oe == 1) && (le == 1) && (ue == 1)) {
                            (q-WH)->w = (q)->w;
                            (q-1)->u  = (q)->u;
                            (q)->v    = (q-Width)->v;   };

                        if ((oe == 1) && (re == 1) && (de == 1)) {
                            (q-WH)->w    = (q)->w;
                            (q)->u       = (q-1)->u;
                            (q-Width)->v = (q)->v;      };

                        if ((oe == 1) && (le == 1) && (de == 1)) {
                            (q-WH)->w = (q)->w;
                            (q-1)->u  = (q)->u;
                            (q-Width)->v = (q)->v;      };


                        if ((ie == 1) && (re == 1) && (ue == 1)) {
                            (q)->w = (q-WH)->w;
                            (q)->u    = (q-1)->u;
                            (q)->v    = (q-Width)->v;   };

                        if ((ie == 1) && (le == 1) && (ue == 1)) {
                            (q)->w = (q-WH)->w;
                            (q-1)->u  = (q)->u;
                            (q)->v    = (q-Width)->v;   };

                        if ((ie == 1) && (re == 1) && (de == 1)) {
                            (q)->w = (q-WH)->w;
                            (q)->u    = (q-1)->u;
                            (q-Width)->v = (q)->v;      };

                        if ((ie == 1) && (le == 1) && (de == 1)) {
                            (q)->w = (q-WH)->w;
                            (q-1)->u  = (q)->u;
                            (q-Width)->v = (q)->v;      };

                        //  two opposite empty-faces and one not.
                        //                       (12 sub_configs)
                    } else {
                        //  Here we solve for the velocity for the empty
                        //  face which is opposite a non-empty face the
                        //  same way we do for the one empty-face case,
                        //  ignoring the fact that two of the other faces
                        //  are by air.

                        if ((oe == 1) && (ie == 1) && (re == 1)) {
                            (q)->u = (q-1)->u - (q)->v + (q-Width)->v -
                            (q)->w + (q-WH)->w;     };

                        if ((oe == 1) && (ie == 1) && (le == 1)) {
                            (q-1)->u = (q)->u + (q)->v - (q-Width)->v +
                            (q)->w - (q-WH)->w;     };

                        if ((oe == 1) && (ie == 1) && (ue == 1)) {
                            (q)->v = - (q)->u + (q-1)->u + (q-Width)->v -
                            (q)->w + (q-WH)->w;     };

                        if ((oe == 1) && (ie == 1) && (de == 1)) {
                            (q-Width)->v = (q)->u - (q-1)->u + (q)->v +
                            (q)->w - (q-WH)->w;     };


                        if ((ue == 1) && (de == 1) && (re == 1)) {
                            (q)->u = (q-1)->u - (q)->v + (q-Width)->v -
                            (q)->w + (q-WH)->w;     };

                        if ((ue == 1) && (de == 1) && (le == 1)) {
                            (q-1)->u = (q)->u + (q)->v - (q-Width)->v +
                            (q)->w - (q-WH)->w;     };

                        if ((ue == 1) && (de == 1) && (ie == 1)) {
                            (q)->w = - (q)->u + (q-1)->u - (q)->v +
                            (q-Width)->v + (q-WH)->w;     };

                        if ((ue == 1) && (de == 1) && (oe == 1)) {
                            (q-WH)->w = (q)->u - (q-1)->u + (q)->v -
                            (q-Width)->v + (q)->w;     };


                        if ((re == 1) && (le == 1) && (ue == 1)) {
                            (q)->v = - (q)->u + (q-1)->u + (q-Width)->v -
                            (q)->w + (q-WH)->w;     };

                        if ((re == 1) && (le == 1) && (de == 1)) {
                            (q-Width)->v = (q)->u - (q-1)->u + (q)->v +
                            (q)->w - (q-WH)->w;     };

                        if ((re == 1) && (le == 1) && (ie == 1)) {
                            (q)->w = - (q)->u + (q-1)->u - (q)->v +
                            (q-Width)->v + (q-WH)->w;     };

                        if ((re == 1) && (le == 1) && (oe == 1)) {
                            (q-WH)->w = (q)->u - (q-1)->u + (q)->v -
                            (q-Width)->v + (q)->w;     };
                    };
                };

                //  four empty-faces    (15 configs -> 3 + 12)
                //  ----
                if (empty_counter == 4) {

                    //  The four empty-face have an opposite empty face.
                    //                          (3 sub_configs)
                    if ( ((ue == 0) && (de == 0)) ||
                         ((re == 0) && (le == 0)) ||
                         ((ie == 0) && (oe == 0)) )   {
                        //  compute the div of the cell (ignoring the
                        //  fact that four of the other faces are by air)
                        //  and add (or '-') a quarter of that Div
                        //  equally to each empty-face.

                        if ((ue == 0) && (de == 0))  {
                            temp = 0.25 * ( (q)->u - (q-1)->u +
                                            (q)->v - (q-Width)->v +
                                            (q)->w - (q-WH)->w);
                            (q)->u -= temp;    (q-1)->u  += temp;
                            (q)->w -= temp;    (q-WH)->w += temp;    };

                        if ((re == 0) && (le == 0))  {
                            temp = 0.25 * ( (q)->u - (q-1)->u +
                                            (q)->v - (q-Width)->v +
                                            (q)->w - (q-WH)->w);
                            (q)->v -= temp;    (q-Width)->v += temp;
                            (q)->w -= temp;    (q-WH)->w    += temp;  };

                        if ((ie == 0) && (oe == 0))  {
                            temp = 0.25 * ( (q)->u - (q-1)->u +
                                            (q)->v - (q-Width)->v +
                                            (q)->w - (q-WH)->w);

                            (q)->u -= temp;    (q-1)->u += temp;
                            (q)->v -= temp;    (q-Width)->v = temp;  };

                        //  two are opposite and two are not.
                        //                      (12 sub_configs)
                    } else {
                        //  copy from the opposite side face those not
                        //  opposite and again add them the half Div
                        //  of the other two empty-faces.
                        
                        if ((oe == 0) && (ue == 0))  {
                            temp = 0.5 * ( (q)->u - (q-1)->u );
                            (q)->w = (q-WH)->w - temp;
                            (q-Width)->v = (q)->v + temp;   };

                        if ((oe == 0) && (de == 0))  {
                            temp = 0.5 * ( (q)->u - (q-1)->u );
                            (q)->w = (q-WH)->w - temp;
                            (q)->v = (q-Width)->v - temp;   };

                        if ((oe == 0) && (re == 0))  {
                            temp = 0.5 * ( (q)->v - (q-Width)->v );
                            (q)->w = (q-WH)->w - temp;
                            (q-1)->u = (q)->u + temp;   };

                        if ((oe == 0) && (le == 0))  {
                            temp = 0.5 * ( (q)->v - (q-Width)->v );
                            (q)->w = (q-WH)->w - temp;
                            (q)->u = (q-1)->u - temp;   };


                        if ((ie == 0) && (ue == 0))  {
                            temp = 0.5 * ( (q)->u - (q-1)->u );
                            (q-WH)->w = (q)->w + temp;
                            (q-Width)->v = (q)->v + temp;   };

                        if ((ie == 0) && (de == 0))  {
                            temp = 0.5 * ( (q)->u - (q-1)->u );
                            (q-WH)->w = (q)->w + temp;
                            (q)->v = (q-Width)->v - temp;   };

                        if ((ie == 0) && (re == 0))  {
                            temp = 0.5 * ( (q)->v - (q-Width)->v );
                            (q-WH)->w = (q)->w + temp;
                            (q-1)->u = (q)->u + temp;   };

                        if ((ie == 0) && (le == 0))  {
                            temp = 0.5 * ( (q)->v - (q-Width)->v );
                            (q-WH)->w = (q)->w + temp;
                            (q)->u = (q-1)->u - temp;   };


                        if ((re == 0) && (ue == 0))  {
                            temp = 0.5 * ( (q)->w - (q-WH)->w );
                            (q-1)->u = (q)->u + temp;
                            (q-Width)->v = (q)->v + temp;   };

                        if ((re == 0) && (de == 0))  {
                            temp = 0.5 * ( (q)->w - (q-WH)->w );
                            (q-1)->u = (q)->u + temp;
                            (q)->v = (q-Width)->v - temp;   };

                        if ((le == 0) && (ue == 0))  {
                            temp = 0.5 * ( (q)->w - (q-WH)->w );
                            (q)->u = (q-1)->u - temp;
                            (q-Width)->v = (q)->v + temp;   };

                        if ((le == 0) && (de == 0))  {
                            temp = 0.5 * ( (q)->w - (q-WH)->w );
                            (q)->u = (q-1)->u - temp;
                            (q)->v = (q-Width)->v - temp;   };
                    };
                };

                //  five empty-faces       (6 configs )
                //  ----
                if (empty_counter == 5) {
                    
                    //  We treat this case the same as
                    //  the one empty-face case.
                    if (re == 0) {    
                        (q)->u = (q-1)->u - (q)->v + (q-Width)->v -
                        (q)->w + (q-WH)->w;     };
                    if (le == 0) {
                        (q-1)->u = (q)->u + (q)->v - (q-Width)->v +
                        (q)->w - (q-WH)->w;     };
                    if (ue == 0) {
                        (q)->v = - (q)->u + (q-1)->u + (q-Width)->v -
                        (q)->w + (q-WH)->w;     };
                    if (de == 0) {
                        (q-Width)->v = (q)->u - (q-1)->u + (q)->v +
                        (q)->w - (q-WH)->w;     };
                    if (ie == 0) {
                        (q)->w = - (q)->u + (q-1)->u - (q)->v +
                        (q-Width)->v + (q-WH)->w;     };
                    if (oe == 0) {
                        (q-WH)->w = (q)->u - (q-1)->u + (q)->v -
                        (q-Width)->v + (q)->w;   };
                };


            };


            q++;   // next CONTAINER cell

        };
    };
};


}











//  F O R   N O W ,  O N L Y   N O T   M O V I N G   O B S T A C L E S  !
//                   A N D   A T   L E A S T   2 - p i x e l   T H I C K
//                   I N   A N Y   D I R E C T I O N .

//  This 2 methods take care of the solid-fluid interface (no-slip BCs),

//  • 1st : takes care of the velocity on the face shared by both
//          the solid cell and the fluid cell, the velocity must be
//          equal to zero (no-slip).
//                                               |
//              (6 conditions)            solid '>' fluid
//                                               |
//
//          It is needs to be called "before" the 'uvw_DivFreeCs' method.

//  • 2nd : takes care of the 2 remaining velocities (orthogonal vels),
//          these require to be opposite to those in the fluid.
//
//                                        solid  |  fluid
//             (18 conditions)          ---'^'---+--- v ---
//                                        solid  |  fluid
//
//          It is needs to be called "after" the 'uvw_DivFreeCs' method,
//          when the velocity to be copied is well-balanced according to
//          the Divergence Free property.

- (void)uvw_BCs_1
{

struct     Fdata      *q;  //  point to a Fluid grid cell.
long int   i,j,k;          //  GRID counters.

q = Q;                              // point to the start of the GRID
for (k=0; k<Depth; k++) {
    for (j=0; j<Height; j++) {
        for (i=0; i<Width; i++) {

            if ( (q)->cell_type > 1 ) {    // fluid_cell Checking

                //              |
                //      FLUID  '0' solid       for each i,j,k
                //              |
                //
                if ( i < (Width-1) ) {
                    if ( (q+1)->cell_type == 1 ) {   (q)->u = 0;   };
                } else {
                         // pretend to find a solid-cell (temp measure)
                         (q)->u = 0;   };

                if ( j < (Height-1) ) {
                    if ( (q+Width)->cell_type == 1 ) {   (q)->v = 0;   };
                } else {
                         // pretend to find a solid-cell (temp measure)
                         (q)->v = 0;   };

                if ( k < (Depth-1) ) {
                    if ( (q+WH)->cell_type == 1 ) {   (q)->w = 0;   };
                } else {
                         // pretend to find a solid-cell (temp measure)
                         (q)->w = 0;   };
            };


            if ( (q)->cell_type == 1 ) {    // solid_cell Checking

                //              |
                //      SOLID  '0' fluid       for each i,j,k
                //              |
                //
                if ( i < (Width-1) )  {
                    if ( (q+1)->cell_type     > 1 )  {  (q)->u = 0;  };
                };
                if ( j < (Height-1) ) {
                    if ( (q+Width)->cell_type > 1 )  {  (q)->v = 0;  };
                };
                if ( k < (Depth-1) )  {
                    if ( (q+WH)->cell_type    > 1 )  {  (q)->w = 0;  };
                };
            };


            q++;   // next CONTAINER cell

        };
    };
};

}





//  pay attention :   If the obstacles (wall not included) were only
//                    1-pixel thick then the following method has to
//                    be modified !

- (void)uvw_BCs_2
{

struct     Fdata      *q;  //  point to a Fluid grid cell.
long int   i,j,k;          //  GRID counters.

q = Q;                              // point to the start of the GRID
for (k=0; k<Depth; k++) {
    for (j=0; j<Height; j++) {
        for (i=0; i<Width; i++) {

            if ( (q)->cell_type == 1 ) {        //  solid_cell Checking 

                //   e.g.             ---+
                //                 |    /|
                //                 |--+ •----- another solid/empty cell
                //     solid/empty | /|
                //    ------'^'----+---- v ----    Check for each i,j,k
                //        SOLID    |   fluid
                //
                if ( (i<(Width-1)) && (j<(Height-1)) && (k<(Depth-1)) ) {

                    //  (i)           a = is for 'any';      se a 
                    //                                        S f
                    //                                         \__'se'
                    //                                          behind S
                    if ( ((q+1)->cell_type > 1) &&
                         ((q+WH)->cell_type < 2) &&
                         ((q+Width)->cell_type < 2)  )  {

                        (q)->v = -(q+1)->v;              };
                    //                                       se a
                    //                                        S se
                    //                                         \__'f'
                    if ( ((q+1)->cell_type < 2) &&
                         ((q+WH)->cell_type > 1) &&
                         ((q+Width)->cell_type < 2)  )  {

                        (q)->v = -(q+WH)->v;             };
                    //                                       se a
                    //                                        S f
                    //                                         \__'f'
                    if ( ((q+1)->cell_type > 1) &&
                         ((q+WH)->cell_type > 1) &&
                         ((q+Width)->cell_type < 2)  )  {

                        (q)->v = -0.5*( (q+WH)->v + (q+1)->v);   };

                    //   (ii)                                 f a
                    //                                        S se
                    //                                         \__'se'
                    //                                          behind S
                    if ( ((q+Width)->cell_type > 1) &&
                         ((q+WH)->cell_type < 2) &&
                         ((q+1)->cell_type < 2)  )      {

                        (q)->u = -(q+Width)->u;          };
                    //                                       se a
                    //                                        S se
                    //                                         \__'f'
                    if ( ((q+Width)->cell_type < 2) &&
                         ((q+WH)->cell_type > 1) &&
                         ((q+1)->cell_type < 2)  )      {

                        (q)->u = -(q+WH)->u;          };
                    //                                        f a
                    //                                        S se
                    //                                         \__'f'
                    if ( ((q+Width)->cell_type > 1) &&
                         ((q+WH)->cell_type > 1) &&
                         ((q+1)->cell_type < 2)  )      {

                        (q)->u = -0.5*( (q+WH)->u + (q+Width)->u );  };

                    //   (iii)                               f
                    //                                'a'__/ S se
                    //                                        \__'se'
                    //                                         behind S
                    if ( ((q+Width)->cell_type > 1) &&
                         ((q+1)->cell_type  < 2) &&
                         ((q+WH)->cell_type < 2)  )  {

                        (q)->w = -(q+Width)->w;          };
                    //                                      se
                    //                                'a'__/ S f
                    //                                        \__'se'
                    if ( ((q+Width)->cell_type < 2) &&
                         ((q+1)->cell_type  > 1) &&
                         ((q+WH)->cell_type < 2)  )  {

                        (q)->w = -(q+1)->w;              };
                    //                                      se
                    //                                'a'__/ S f
                    //                                        \__'f'
                    if ( ((q+Width)->cell_type > 1) &&
                         ((q+1)->cell_type  > 1) &&
                         ((q+WH)->cell_type < 2)  )  {

                        (q)->w = -0.5*( (q+1)->w + (q+Width)->w );  };
                };

                //   e.g.               another   ---+
                //                 |   solid/empty  /|
                //                 |          ----+  |
                //          any    |solid/empty  /|
                //    ------ ^ ----+----'v'----+  |   for each i,j,k
                //         fluid   |   SOLID   |
                //
                if ( (i>0) && (j>0) && (k>0) ) {

                    //   (iv)                                 a se
                    //                                        f S
                    //                                   'se'__/
                    //                                    in front of S
                    if ( ((q-1)->cell_type > 1) &&
                         ((q-WH)->cell_type < 2) &&
                         ((q+Width)->cell_type < 2)  )  {

                        (q)->v = -(q-1)->v;              };
                    //                                        a se
                    //                                       se S
                    //                                    'f'__/
                    //                                    in front of S
                    if ( ((q-1)->cell_type < 2) &&
                         ((q-WH)->cell_type > 1) &&
                         ((q+Width)->cell_type < 2)  )  {

                        (q)->v = -(q-WH)->v;             };
                    //                                        a se
                    //                                        f S
                    //                                    'f'__/
                    //                                    in front of S
                    if ( ((q-1)->cell_type > 1) &&
                         ((q-WH)->cell_type > 1) &&
                         ((q+Width)->cell_type < 2)  )  {

                        (q)->v = -0.5*( (q-1)->v + (q-WH)->v );   };
                    
                    //   (v)                                  S se
                    //                                'se'__/ f  a
                    //                            front of S
                    if ( ((q-Width)->cell_type > 1) &&
                         ((q-WH)->cell_type < 2) &&
                         ((q+1)->cell_type < 2)  )      {

                        (q)->u = -(q-Width)->u;    };
                    //                                        S se
                    //                                 'f'__/ se a
                    //                            front of S
                    if ( ((q-Width)->cell_type < 2) &&
                         ((q-WH)->cell_type > 1) &&
                         ((q+1)->cell_type < 2)  )      {

                        (q)->u = -(q-WH)->u;       };
                    //                                        S se
                    //                                 'f'__/ f  a
                    //                            front of S
                    if ( ((q-Width)->cell_type > 1) &&
                         ((q-WH)->cell_type > 1) &&
                         ((q+1)->cell_type < 2)  )      {
                        
                        (q)->u = -0.5*( (q-WH)->u + (q-Width)->u);  };

                    //   (vi)                              se S
                    //                                 'se'__/f
                    //                                         \__'a'
                    //                             (no error!) behind S
                    if ( (k<(Depth-1)) &&
                         ((q-Width)->cell_type > 1) &&
                         ((q-1)->cell_type < 2) &&
                         ((q+WH)->cell_type < 2)  )  {

                        (q)->w = -(q-Width)->w;      };
                    //                                      f S
                    //                                 'se'__/se
                    //                                         \__'a'
                    //                                       behind S
                    if ( (k<(Depth-1)) &&
                         ((q-Width)->cell_type < 2) &&
                         ((q-1)->cell_type > 1) &&
                         ((q+WH)->cell_type < 2)  )  {
                        
                        (q)->w = -(q-1)->w;      };
                    //                                      f S
                    //                                 'se'__/f
                    //                                         \__'a'
                    //                                       behind S
                    if ( (k<(Depth-1)) &&
                         ((q-Width)->cell_type > 1) &&
                         ((q-1)->cell_type > 1) &&
                         ((q+WH)->cell_type < 2)  )  {
                        
                        (q)->w = -0.5*( (q-Width)->w + (q-1)->w);  };
                };

            };


            q++;   // next CONTAINER cell

        };
    };
};

}




















//  ------O-N-L-Y---2-D---V-E-R-S-I-O-N-------------------------



//     +---+---+---+---+---+---+---+
//     | • | • | • | • | • | • | • |  ---> boundary strip (outside domain)
//     +---+---+---+---+---+---+---+  ---> boundary domain (linear)
//     | • |   |   | x |   |   | • |
//       0   1       |     imax imax+1
//                   |
//                   +-------------------> inside domain


- (void)u_BC
{

//  u no_slip Boundary Conditions
//  •••••••••••••••••••••••••••••

// adjust the  u-Velocity on the Boundary strip to make 
// the velocity = 0 on the Boundary domain (no-slip BC) :
//
//              U(i,jmax+1) = -u(i,jmax);
//                     +--------+
//                     | +----+ |
//        U(0,j) = 0;  | |    | |  U(imax,j) = 0;
//                     | +----+ |
//                     +--------+
//                  U(i,0) = -u(i,1);
//
// where:  U = u(uiter+1);  u = u(uiter);


struct Fdata        *q;        // point to a Fluid grid cell
unsigned long int   i,j;       // counters


// first column (left of the grid)     ==> column [0,j];     j = 1,..,jmax
q = Q;                              // point to the start of the grid
for (j=1; j<Hminus1; j++)    {
	q += Width;                     // jump to the 2nd row (first time)
	(q)->u = 0;	 };                 // U(0,j) = 0

// first row    (bottom of the grid)   ==> row [i,0];        i = 1,..,imax
q = Q;                              // point to the start of the grid
for (i=1; i<Wminus1; i++)     {
	q++;                            // jump to 2nd cell of the row (1st t)
	(q)->u = -(q+Width)->u; };      // U(i,0) = -u(i,1)

// last row     (top of the grid)      ==> row [i,jmax+1];   i = 1,..,imax
q = Q + GridCells - Width ;         // point to the last row
for (i=1; i<Wminus1; i++)     {
	q++;                            // jump to 2nd cell of the row (1st t)
//  (q)->u  = 10;                   // ••• Lid_Driven Test
	(q)->u = -(q-Width)->u; };      // U(i,jmax+1) = -u(i,jmax)

// last but one row column  (right)    ==> column [imax,j];  j = 1,..,jmax
q = Q + Width - 2;                  // point to the last but one column
for (j=1; j<Hminus1; j++)    {
	q += Width;                     // jump to the 2nd row (first time)
	(q)->u = 0;	 };                 // U(imax,j) = 0

}





- (void)v_BC
{

//  v no_slip Boundary Conditions
//  •••••••••••••••••••••••••••••

// adjust the  v-Velocity on the Boundary strip to make 
// the velocity = 0 on the Boundary domain (no-slip BC) :
//
//                   V(i,jmax) = 0);
//                     +--------+
//                     | +----+ |
//   V(0,j) = -v(1,j); | |    | | V(imax+1,j) = -v(imax,j);
//                     | +----+ |
//                     +--------+
//                    V(i,0) = 0;
//
// where:  V = v(viter+1);  v = v(viter);


struct Fdata        *q;        // point to a Fluid grid cell
unsigned long int   i,j;       // counters


// first column (left of the grid)     ==> column [0,j];     j = 1,..,jmax
q = Q;                              // point to the start of the grid
for (j=1; j<Hminus1; j++)    {
	q += Width;                     // jump to the 2nd row (first time)
	(q)->v = -(q+1)->v;	 };         // V(0,j) = -v(1,j)

// first row    (bottom of the grid)   ==> row [i,0];        i = 1,..,imax
q = Q;                              // point to the start of the grid
for (i=1; i<Wminus1; i++)     {
	q++;                            // jump to 2nd cell of the row (1st t)
	(q)->v = 0;          };         // V(i,0) = 0

// last but one row (top of the grid)  ==> row [i,jmax];     i = 1,..,imax
q = Q + GridCells - Width - Width;  // point to the last but one row
for (i=1; i<Wminus1; i++)     {
	q++;                            // jump to 2nd cell of the row (1st t)
	(q)->v = 0;          };         // V(i,jmax) = 0

// last column  (right of the grid)    ==> column [imax+1,j];j = 1,..,jmax
q = Q + Width - 1;                  // point to the last column
for (j=1; j<Hminus1; j++)    {
	q += Width;                     // jump to the 2nd row (first time)
	(q)->v = -(q-1)->v;	 };         // V(imax+1,j) = -v(imax,j)

}





- (void)p_BC
{

//  p-gradient Boundary Conditions
//  ••••••••••••••••••••••••••••••

//  copy of the Pressure on the Boundary strip   (see Griebel p. 38)
//  (to make the pressure gradient = 0 on the Boundary domain)
//
//               P(i,jmax+1) = p(i,jmax);
//                     +--------+
//                     | +----+ |
//    P(0,j) = p(1,j); | |    | | P(imax+1,j) = p(imax,j);
//                     | +----+ |
//                     +--------+
//                  P(i,0) = p(i,1);
//
//  where:  P = p (iter+1);  p = p(iter);


struct Fdata        *q;        // point to a Fluid grid cell
unsigned long int   i,j;       // counters


// first column (left of the grid)     ==> column [0,j];     j = 1,..,jmax
q = Q;                              // point to the start of the grid
for (j=1; j<Hminus1; j++)    {
	q += Width;                     // jump to the 2nd row (first time)
	(q)->p = (q+1)->p;       };     // P(0,j) = p(1,j)

// first row    (bottom of the grid)   ==> row [i,0];        i = 1,..,imax
q = Q;                              // point to the start of the grid
for (i=1; i<Wminus1; i++)    {
	q++;                            // jump to 2nd cell of the row (1st t)
	(q)->p = (q+Width)->p;   };     // P(i,0) = p(i,1)

// last row     (top of the grid)      ==> row [i,jmax+1];   i = 1,..,imax
q = Q + GridCells - Width ;         // point to the last row
for (i=1; i<Wminus1; i++)    {
	q++;                            // jump to 2nd cell of the row (1st t)
	(q)->p = (q-Width)->p;   };     // P(i,jmax+1) = p(i,jmax)

// last column  (right of the grid)    ==> column [imax+1,j];j = 1,..,jmax
q = Q + Width - 1;                  // point to the last column
for (j=1; j<Hminus1; j++)    {
	q += Width;                     // jump to the 2nd row (first time)
	(q)->p = (q-1)->p;	     };     // P(imax+1,j) = p(imax,j)

}



@end
