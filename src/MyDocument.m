//
//  MyDocument.m
//  DigitalFluid
//
//  Created by G5 G5 on Sun May 01 2005.
//  Copyright (c) 2005 Ivano Ras. All rights reserved.
//

#import "MyDocument.h"

@implementation MyDocument

- (id)init
{
    self = [super init];
    if (self) {
    
        // Add your subclass-specific initialization here.
        // If an error occurs here, send a [self release] message and return nil.
    }
    return self;
}


- (void)makeWindowControllers
{
    //	 windowController = [[DFWindowController alloc] initWithWindowNibName:@"DFWindow"];
    //	 [self addWindowController:windowController];
    // the previous two rows are grouped together into the next row to avoid a warning !
	[self addWindowController:[[DFWindowController alloc] initWithWindowNibName:@"DFWindow"]];
}


- (void)windowControllerDidLoadNib:(NSWindowController *) aController
{
    [super windowControllerDidLoadNib:aController];
    // Add any code here that needs to be executed once the windowController
	// has loaded the document's window.
}


- (NSData *)dataRepresentationOfType:(NSString *)aType
{
    // Insert code here to write your document from the given data.
	// You can also choose to override -fileWrapperRepresentationOfType:
	// or -writeToFile:ofType: instead.
    return nil;
}


- (BOOL)loadDataRepresentation:(NSData *)data ofType:(NSString *)aType
{
	activeImage = [[NSImage alloc] initWithData:data];
	return (activeImage != nil);
}


- (NSImage *)activeImage
{
	return activeImage;
}


@end
