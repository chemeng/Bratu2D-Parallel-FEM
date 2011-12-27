//
//  pause.c
//  Bratu-C
//
//  Created by Tim Ioannidis on 12/19/11.
//  Copyright 2011 __MyCompanyName__. All rights reserved.
//

#include <stdio.h>
#include "methods_decl.h"

void mypause()
{
	char pause;
	printf("Press [Enter] to continue. . .");
	scanf("%c",&pause);
}

