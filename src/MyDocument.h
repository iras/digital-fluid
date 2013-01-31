//
//  MyDocument.h
//  DigitalFluid
//
//  Created by G5 G5 on Sun May 01 2005.
//  Copyright (c) 2005 Ivano Ras. All rights reserved.
//


#import <Cocoa/Cocoa.h>
#import "DFWindowController.h"

@interface MyDocument : NSDocument
{
	DFWindowController *windowController;
	NSImage *activeImage;
}

-(NSImage *)activeImage;

@end
