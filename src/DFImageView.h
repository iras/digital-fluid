//
//  DFImageView.h
//  DigitalFluid
//
//  Created by G5 G5 on Sun May 01 2005.
//  Copyright (c) 2005 Ivano Ras. All rights reserved.
//

#import "DFDyn_constants.h"
#import "DFDynamics.h"
#import <AppKit/AppKit.h>


@interface DFImageView : NSImageView
{
	NSImage *destImage;
	NSBitmapImageRep *destImageRep;
	unsigned char *destData;

	IBOutlet id timeField;
	IBOutlet id iterField;
	IBOutlet id uiterField;
	IBOutlet id viterField;

	NSTimer *timer;

	DFDynamics *compute;    // declare a DFDynamics object

	struct Pdata {
		double x;
		double y;
	}  particle[Width][Height];  // declare particle[ ][ ] as a pointer to a Pdata type


@public
	unsigned int DisplayModalityFlag;
    unsigned int mouseActive;
}


- (NSImage *)displayField;
- (NSImage *)displayParticles;
- (void)displayParticlesInit;

- (void)stepAnim;
- (void)startTimer;
- (void)stopTimer;
- (void)resetTimer;

- (NSImage *)testFilter:(NSImage *)srcImage;


@end
