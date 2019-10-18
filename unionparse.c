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

#include "filebuff.h"
#include "pherror.h"
#include "unionparse.h"

UnionEntry * UnionEntry_init(unsigned size, unsigned num) {
	
	UnionEntry *dest;
	
	dest = smalloc(sizeof(UnionEntry));
	dest->len = 0;
	dest->size = size;
	dest->num = 0;
	dest->fSize = num;
	dest->target = smalloc(size);
	dest->filenames = smalloc(num * sizeof(unsigned));
	
	return dest;
}

void UnionEntry_destroy(UnionEntry *src) {
	
	free(src->target);
	free(src->filenames);
	free(src);
}

char ** UnionEntry_getHeader(FileBuff *src, unsigned *n) {
	
	unsigned char *buff, *name, **filenames;
	int size, avail, num, namesize, (*buffFileBuff)(FileBuff *);
	
	/* init */
	avail = src->bytes;
	buff = src->next;
	buffFileBuff = src->buffFileBuff;
	if(avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	
	/* get num */
	num = 0;
	while((size = *buff++) != '\t') {
		num = num * 10 + (size - '0');
		/* rebuff */
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				return 0;
			}
			buff = src->buffer;
		}
	}
	
	/* rebuff */
	if(--avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	
	/* get filenames */
	*n = num;
	filenames = smalloc(num * sizeof(unsigned char *));
	*filenames = smalloc((namesize = 32));
	num = 0;
	size = namesize;
	name = *filenames - 1;
	while((*++name = *buff++) != '\n') {
		/* rebuff */
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				return 0;
			}
			buff = src->buffer;
		}
		/* realloc */
		if(--size == 0) {
			size = namesize;
			namesize <<= 1;
			filenames[num] = realloc(filenames[num], namesize);
			if(!filenames[num]) {
				ERROR();
			}
			name = (unsigned char *)(filenames[num]) + size - 1;
		}
		
		/* end of filename */
		if(*name == '\t') {
			*name = 0;
			/* realloc for file suffixes */
			if(size < 6 && !(filenames[num] = realloc(filenames[num], namesize + 8))) {
				ERROR();
			}
			filenames[++num] = smalloc((namesize = 32));
			size = namesize;
			name = filenames[num] - 1;
		}
	}
	*name = 0;
	/* realloc for file suffixes */
	if(size < 6 && !(filenames[num] = realloc(filenames[num], namesize + 8))) {
		ERROR();
	}
	
	/* set filebuff */
	src->bytes = avail - 1;
	src->next = buff;
	
	return (char **) filenames;
}

int UnionEntry_get(FileBuff *src, UnionEntry *dest) {
	
	unsigned char *buff, *name, c;
	int size, avail, num, (*buffFileBuff)(FileBuff *);
	
	/* init */
	avail = src->bytes;
	buff = src->next;
	buffFileBuff = src->buffFileBuff;
	if(avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	
	/* get target */
	size = dest->size;
	name = (unsigned char *) dest->target - 1;
	while((*++name = *buff++) != '\t') {
		/* rebuff */
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				return 0;
			}
			buff = src->buffer;
		}
		/* realloc */
		if(--size == 0) {
			size = dest->size;
			dest->size <<= 1;
			dest->target = realloc(dest->target, dest->size);
			if(!dest->target) {
				ERROR();
			}
			name = (unsigned char *)(dest->target) + size - 1;
		}
	}
	*name = 0;
	dest->len = dest->size - size + 1;
	if(--avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	
	/* get num */
	num = 0;
	while((size = *buff++) != '\t') {
		num = num * 10 + (size - '0');
		/* rebuff */
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				return 0;
			}
			buff = src->buffer;
		}
	}
	dest->num = num;
	
	/* rebuff */
	if(--avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	
	/* get filenames */
	num = 0;
	size = 0;
	while((c = *buff++) != '\n') {
		/* rebuff */
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				return 0;
			}
			buff = src->buffer;
		}
		
		if(c != '\t') {
			size = size * 10 + (c - '0');
		} else {
			dest->filenames[num++] = size;
			size = 0;
		}
	}
	dest->filenames[num] = size;
	
	/* set filebuff */
	src->bytes = avail - 1;
	src->next = buff;
	return 1;
}
