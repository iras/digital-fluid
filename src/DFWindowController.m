//
//  DFWindowController.m
//  DigitalFluid
//
//  Created by G5 G5 on Sun May 01 2005.
//  Copyright (c) 2005 Ivano Ras. All rights reserved.
//

#import "DFWindowController.h"
#import "DFImageView.h"
#import <MyDocument.h>


@implementation DFWindowController

- (void)windowDidLoad
{
	NSImage *image = [[self document] activeImage];
	[view setImage:image];
}


- (IBAction)Go:(id)sender
{
	[view startTimer];
}


- (IBAction)Alt:(id)sender
{
	[view stopTimer];
}


- (IBAction)Step:(id)sender
{
	[view stepAnim];
}


- (IBAction)Reset:(id)sender
{
	[view resetTimer];
}


- (IBAction)DisplayP:(id)sender
{
	view->DisplayModalityFlag = 2;
}


- (IBAction)DisplayF:(id)sender
{
	view->DisplayModalityFlag = 1;
}


- (IBAction)Test:(id)sender
{
	NSImage *image = [[self document] activeImage];
	NSImage *newImage = [view testFilter:image];
	[view setImage:newImage];
}

@end