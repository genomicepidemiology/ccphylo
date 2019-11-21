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

#include "meth.h"
#include "pherror.h"

/* here */
/* remember Makefile */

MethMotif * newMethMotif(int len) {
	
	int size;
	MethMotif *dest;
	
	size = (len << 5) + (len & 31);
	dest = smalloc(sizeof(MethMotif) + size * (sizeof(long unsigned) + sizeof(unsigned)));
	dest->len = len;
	dest->size = size;
	dest->motif = (long unsigned *)(dest + 1);
	dest->mask = (unsigned *)(dest->motif + size);
	dest->next = 0;
	
	return dest;
}

void destroyMethMotifs(MethMotif *src) {
	
	MethMotif *src_next;
	
	while(src) {
		src_next = src->next;
		free(src);
		src = src_next;
	}
}

long unsigned matchMotif32(long unsigned seq, long unsigned *motif, unsigned num) {
	
	/* matches a 32 bp sequence to a 32 bp motif */
	
	/* num = # of alternate motifs, e.g:
	N included -> num = 4
	Y -> num = 2
	*/
	
	long unsigned match, mer;
	
	/* mark unmatched parts */
	mer = *motif ^ seq;
	/* mirror nibbles, and mask them out */
	match = mer | ((mer << 1) & 0xAAAAAAAAAAAAAAAA) | ((mer >> 1) & 0x5555555555555555);
	while(--num && match) {
		mer = *++motif ^ seq;
		//mer |= (((mer << 1) & 0xAAAAAAAAAAAAAAAA) | ((mer >> 1) & 0x5555555555555555));
		match &= mer | ((mer << 1) & 0xAAAAAAAAAAAAAAAA) | ((mer >> 1) & 0x5555555555555555);
	}
	
	return match;
}
