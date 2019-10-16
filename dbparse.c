/* Philip T.L.C. Clausen Jan 2017 plan@dtu.dk */

/*
 * Copyright (c) 2017, Philip Clausen, Technical University of Denmark
 * All rights reserved.
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 *		http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
*/

#include <stdio.h>
#include <stdlib.h>
#include "dbparse.h"
#include "pherror.h"
#include "qseqs.h"

char * nameLoad(Qseqs *name, FILE *infile) {
	
	int c, size, len;
	unsigned char *ptr;
	
	ptr = name->seq;
	size = name->size - 1;
	len = 0;
	while((c = fgetc(infile)) != '\n' && c != EOF) {
		*ptr++ = c;
		++len;
		if(--size == 0) {
			size = name->size - 1;
			name->seq = realloc(name->seq, (name->size <<= 1));
			if(!name->seq) {
				ERROR();
			}
			ptr = name->seq + size;
		}
	}
	*ptr = 0;
	name->len = len;
	if(c == EOF) {
		return 0;
	}
	
	return (char *) name->seq;
}
