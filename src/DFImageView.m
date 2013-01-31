//
//  DFImageView.m
//  DigitalFluid
//
//  Created by G5 G5 on Sun May 01 2005.
//  Copyright (c) 2005 Ivano Ras. All rights reserved.
//

#import "DFImageView.h"
#import "DFDynamics.h"
#import "DFDyn_constants.h"


@implementation DFImageView

- (void)awakeFromNib
{
	unsigned int w = Width  * 5;
	unsigned int h = Height * 5;
	
    // init the Destination Image
	destImage = [[NSImage alloc]  initWithSize:NSMakeSize(w,h)];
	[destImage retain];
	// init the Bitmap Representation of the Dest Image
	destImageRep = [[[NSBitmapImageRep alloc]
                      initWithBitmapDataPlanes:NULL
								    pixelsWide:w
								    pixelsHigh:h
								 bitsPerSample:8
							   samplesPerPixel:BytePerColor
									  hasAlpha:NO
									  isPlanar:NO
							    colorSpaceName:NSCalibratedWhiteColorSpace
								   bytesPerRow:w
								  bitsPerPixel:8] autorelease];
	[destImageRep retain];

	// add the Bitmap representation to the NSImage destImage
	[destImage addRepresentation:destImageRep];

	// pointer to bitmap data of the Dest Image
	destData = [destImageRep bitmapData];

	//--------------------------------

	// DFDynamics object instancing
	compute = [[DFDynamics alloc] init];
    // instance variables init
	[compute instance_vars_init];
	// ensure the existence of this object
	[compute retain];

	//--------------------------------

	[self displayParticlesInit];
	DisplayModalityFlag = 2;
}



- (void)drawRect:(NSRect)rect
{
    [super drawRect:rect];
}



- (void)stepAnim
{
   	// compute one step of dynamics
	[compute dynamics];	

    // update the Window with actual time and iters needed.
	compute->t += dt;
    [timeField setFloatValue:compute->t];
    [iterField setFloatValue:compute->iter];
	[uiterField setFloatValue:compute->uiter];
	[viterField setFloatValue:compute->viter];

	// force to displayParticles
    DisplayModalityFlag = 2;

	// visualize the new data
	if (DisplayModalityFlag == 2) {
	    NSImage *newImage = [self displayParticles];
	    [self setImage:newImage];
	}
	if (DisplayModalityFlag == 1) {
		NSImage *newImage = [self displayField];
		[self setImage:newImage];
	}
	[self setNeedsDisplay:YES];
}



- (void)startTimer
{
	// (0.04 =  25 fps)
	// (0,01 = 100 fps)
	timer = [NSTimer scheduledTimerWithTimeInterval:0.02
										     target:self
                                           selector:@selector(stepAnim)
										   userInfo:nil
                                            repeats:YES];
}


- (void)stopTimer
{
	[timer invalidate];
}



- (void)resetTimer
{
	// reset time and iter
	compute->t = 0;
    [timeField setFloatValue:compute->t];
    [iterField setFloatValue:0];
	[uiterField setFloatValue:0];
	[viterField setFloatValue:0];

	// reset Fluid grid
	[compute instance_vars_init];

	// save the value of the DisplayModalityFlag before
	// the calling of the displayParticlesInit.
    unsigned int flag = DisplayModalityFlag;

	// reset displayParticles grid
	[self displayParticlesInit];

	// restore the value in DisplayModalityFlag.
    DisplayModalityFlag = flag;

	// visualize the new data
	if (DisplayModalityFlag == 2) {
	    NSImage *newImage = [self displayParticles];
	    [self setImage:newImage];
	}
	if (DisplayModalityFlag == 1) {
		NSImage *newImage = [self displayField];
		[self setImage:newImage];
	}
	[self setNeedsDisplay:YES];
}






- (void)mouseDown:(NSEvent *)theEvent
{
	// get the actual position of the mouse.
    NSPoint loc = [theEvent locationInWindow];
	loc.x -= [self frame].origin.x;
	loc.y -= [self frame].origin.y;

	mouseActive = 1;

	// save the first position
	compute->xMouse = loc.x;
	compute->yMouse = loc.y;
	
	// reset the velocity
	compute->uMouse = 0;
	compute->vMouse = 0;
}


- (void)mouseDragged:(NSEvent *)theEvent
{
	// get the actual position of the mouse.
    NSPoint loc = [theEvent locationInWindow];
	loc.x -= [self frame].origin.x;
	loc.y -= [self frame].origin.y;
	
	// get the actual velocity of the mouse.
	compute->uMouse = (loc.x - compute->xMouse)/dt;
	compute->vMouse = (loc.y - compute->yMouse)/dt;
	
	// store the actual position of the mouse.
	compute->xMouse = loc.x;
	compute->yMouse = loc.y;
}


- (void)mouseUp:(NSEvent *)theEvent
{	
	// reset the velocity
	compute->uMouse = 0;
	compute->vMouse = 0;
}






- (NSImage *)displayField
{

//	        (top)      X                   Y
//         +-------->                   A
//  (left) |				transf		|
//         |	image        ===>       |  fluid-grid
//         |	Frame                   |   Frame
//         v                            +------->
//      Y							                X
//
//  the code is written to visualize the data in the fluid-grid way.

	
	unsigned char *p = destData;         // point to image_data
	struct Fdata  *q;                    // point to fluid_grid

	unsigned long int i,j,n;             // counters

    double u2,v2,u_v,v_u,v_ux05,u_vx05;
	unsigned int val;

    unsigned int  uflag, vflag;			 // boolean value
	unsigned int  dir;					 // direction flag

	unsigned long int nextR1 = Width * 1 * 5 * BytePerColor;
	unsigned long int nextR2 = Width * 2 * 5 * BytePerColor;
	unsigned long int nextR3 = Width * 3 * 5 * BytePerColor;

	unsigned long int lenghtGrid = GridCells * 5 * 5 * BytePerColor;


	// Clear Screen
	for (n=0; n<lenghtGrid; n++) { *(p) = 255;  p++; }

	// Write Screen
	p = destData;						// rewind the image pointer
	q = compute->Q + GridCells - Width; // point to last row (fluid_grid)

    for (j=0; j<Height; j++) {
		for (i=0; i<Width; i++) {

			u2 = (q)->u*(q)->u;         // ^2 serve to avoid more 
			v2 = (q)->v*(q)->v;         //    complex sign_checks.

			if ((u2==0.0) && (v2==0.0)) { 
				p++; p++; *(p+nextR2) = 0; p++; p++; p++;}

			else {

				uflag = vflag = 0;      // reset of the boolean values
				dir = 0;

				if (u2>0) { 
					v_u = (q)->v / (q)->u;
					if ((v_u <= 1) && (v_u >= -1)) {
						vflag = 1;
						if ((q)->u < 0) {dir = 1;}
					}
				}
				if (v2>0) {
					u_v = (q)->u / (q)->v;
					if ((u_v <  1) && (u_v >  -1)) {
						uflag = 1;
						if ((q)->v > 0) {dir = 1;}
					}
				}

				if (vflag) {
					v_ux05 = 0.5 * v_u;
					// val = ( 1.5 * v_u + 2) * Wdt12;
					val = (unsigned int)( v_u + v_ux05 + 2) * nextR1;
					*(p + val) = 160*(1-dir);  p++;
					// val = ( 0.5 * v_u + 2) * Wdt12;
					val = (unsigned int)(       v_ux05 + 2) * nextR1;
					*(p + val) = 160;		   p++;
					// val = (-0.5 * v_u + 2) * Wdt12;
					val = (unsigned int)(      -v_ux05 + 2) * nextR1;
					*(p + val) = 160;          p++;
					// val = (-1.5 * v_u + 2) * Wdt12;
					val = (unsigned int)(-v_u - v_ux05 + 2) * nextR1;
					*(p + val) = 160*dir;      p++; p++;
					}

				else {
					u_vx05 = 0.5 * u_v;
					// val = ( 1.5 * u_v + 2);
					val = (unsigned int)( u_v + u_vx05 + 2);
					*(p + val) =          160*(1-dir);
					// val = ( 0.5 * u_v + 2);
					val = (unsigned int)(		u_vx05 + 2);
					*(p + val + nextR1) = 160;
					// val = (-0.5 * u_v + 2);
					val = (unsigned int)(     - u_vx05 + 2);
					*(p + val + nextR2) = 160;
					// val = (-1.5 * u_v + 2);
					val = (unsigned int)(-u_v - u_vx05 + 2);
					*(p + val + nextR3) = 160*dir;

					p++; p++; p++; p++; p++; 
				}

			}

			q++;
		}

        q += - Width - Width;

		p += 20 * Width;    //   20   =    5    *   4    *    1 
						    //   --      side   remaining   bytes x 
	}                       //          pixels    rows      color

	return destImage;
}




- (NSImage *)displayParticles
{

//	        (top)      X                   Y
//         +-------->                   A
//  (left) |				transf		|
//         |	image        ===>       |  fluid-grid
//         |	Frame                   |   Frame
//         v                            +-------->
//      Y							                X
//
//  the code is written to visualize the data in the fluid-grid way.


	unsigned char *p;                    // point to image_data
	struct Fdata  *q;                    // point to fluid_grid

	unsigned long int i,j,n;             // counters
	unsigned long int lenghtGrid = GridCells * 5 * 5 * BytePerColor;
	unsigned long int oldX, oldY;
	unsigned long int Width5 = Width * 5 * BytePerColor;
	
	double currP_x, currP_y;			 // current particles positions


	// clear screen
	p = destData;
	for (n=0; n<lenghtGrid; n++) { *(p) = 255;  p++; }

    // position computation + display
    for (j=1; j<=Height; j++) {
		for (i=1; i<=Width; i++) {

			oldX = (unsigned long int)(particle[i][j].x);
			oldY = (unsigned long int)(particle[i][j].y);

			q = compute->Q + (oldX-1) + (oldY-1)*Width;

			// update x,y for each particle.
			particle[i][j].x += (q)->u * dt;
			particle[i][j].y += (q)->v * dt;

			// display the particle..
			currP_x = particle[i][j].x;
			currP_y = particle[i][j].y;

			if ((currP_x >= 0) && (currP_x <= Width) &&
				(currP_y >= 0) && (currP_y <= Height))  {

			       p = destData + 
				       (unsigned long int)((currP_x - 1) *5 ) +
			           (unsigned long int)((Height - currP_y) *5) *Width5;
                   // ..as a black dot
			       *(p) = 0;

				   //  point DOUBLED on the Boundary Strip
				   if ((currP_x < 2) || (currP_x > (Wminus1)) ||
					   (currP_y < 2) || (currP_y > (Hminus1)))  { *(p+1) = 0; }
				}

		}
	}

	return destImage;
}




- (void)displayParticlesInit
{
	unsigned long int i,j;             // counters

	// particles placement.
    for (j=1; j<=Height; j++) {
		for (i=1; i<=Width; i++) {
			particle[i][j].x = (double)(i);
			particle[i][j].y = (double)(j);
		}   }

}





//=========================================================================
//=========================================================================


- (NSImage *)testFilter:(NSImage *)srcImage
{
	// create the Bitmap Representation of the Source Image
	NSBitmapImageRep *srcImageRep = [NSBitmapImageRep
                     imageRepWithData:[srcImage TIFFRepresentation]];

	int w1 = [srcImageRep pixelsWide];
	int h1 = [srcImageRep pixelsHigh];
	// create the Destination Image
	NSImage *destImage1 = [[NSImage alloc] initWithSize:NSMakeSize(w1,h1)];
	// create the Bitmap Representation of the Dest Image
	NSBitmapImageRep *destImageRep1 = [[[NSBitmapImageRep alloc]
                     initWithBitmapDataPlanes:NULL
								   pixelsWide:w1
								   pixelsHigh:h1
								bitsPerSample:8
							  samplesPerPixel:1
									 hasAlpha:NO
									 isPlanar:NO
							   colorSpaceName:NSCalibratedWhiteColorSpace
								  bytesPerRow:w1
								 bitsPerPixel:8] autorelease];

    // pointer to bitmap data of the Source Image
	unsigned char *srcData = [srcImageRep bitmapData];

	// pointer to bitmap data of the Dest Image
	unsigned char *destData1 = [destImageRep1 bitmapData];

	// pointers to the pixel (x,y) in Source & Dest
	unsigned char *p1, *p2;
	// number of bytes per pixel of the Source Image
	int n = [srcImageRep bitsPerPixel] / 8;
	// counters
	int x, y;

	// main loop
	for ( y = 0; y < h1; y++ ) {
		for ( x = 0; x < w1; x++ ) {
			p1 = srcData + n * (y * w1 + x);
			p2 = destData1 + y * w1 + x;
			// Average of the Source RGB colors -> one byte of the Dest
			p2[0] = (unsigned char)rint((p1[0] + p1[1] + p1[2]) / 3);
		}
	}

	// add the Bitmap representation to the NSImage destImage
	[destImage1 addRepresentation:destImageRep1];

	return destImage1;
}


@end
