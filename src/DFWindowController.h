//
//  DFWindowController.h
//  DigitalFluid
//
//  Created by G5 G5 on Sun May 01 2005.
//  Copyright (c) 2005 Ivano Ras. All rights reserved.
//

#import <AppKit/AppKit.h>
#import "DFImageView.h"


@interface DFWindowController : NSWindowController
{
	IBOutlet DFImageView *view;
}

- (IBAction)Alt:(id)sender;
- (IBAction)Go:(id)sender;
- (IBAction)Step:(id)sender;
- (IBAction)Reset:(id)sender;

- (IBAction)DisplayP:(id)sender;
- (IBAction)DisplayF:(id)sender;

- (IBAction)Test:(id)sender;

@end
