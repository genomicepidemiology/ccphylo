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

#include <stdlib.h>
#include "pherror.h"
#include "ulist.h"

uList * uList_init(unsigned size) {
	
	uList *dest;
	
	dest = smalloc(sizeof(uList));
	dest->n = 0;
	dest->size = size;
	dest->list = smalloc(size * sizeof(unsigned));
	
	return dest;
}

void uList_realloc(uList *src, unsigned newsize) {
	
	if(!(src->list = realloc(src->list, newsize * sizeof(unsigned)))) {
		ERROR();
	}
	src->size = newsize;
}

void uList_destroy(uList *src) {
	
	free(src->list);
	free(src);
}

void uList_push(uList *src, unsigned num) {
	
	if(++src->n == src->size) {
		uList_realloc(src, src->size << 1);
	}
	
	src->list[src->n - 1] = num;
}
