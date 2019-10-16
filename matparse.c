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
#include <string.h>
#include "filebuff.h"
#include "matparse.h"
#include "pherror.h"
#include "qseqs.h"

NucCount * initNucCount(unsigned size) {
	
	NucCount *dest;
	
	dest = smalloc(sizeof(NucCount));
	dest->name = smalloc(128);
	dest->size = 128;
	
	return dest;
}

int FileBuffGetRow(FileBuff *src, NucCount *dest) {
	
	unsigned char *buff, *name;
	int size, avail, num;
	short unsigned *count;
	int (*buffFileBuff)(FileBuff *);
	
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
	
	/* get ref */
	dest->ref = *buff++;
	if(--avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	
	/* evaluate row */
	if(dest->ref == '\n') {
		/* end of entry */
		src->bytes = avail;
		src->next = buff;
		dest->ref = 0;
		return 1;
	} else if(dest->ref == '#') {
		/* new template */
		name = dest->name;
		size = dest->size;
		while((*name++ = *buff++) != '\n') {
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
				dest->name = realloc(dest->name, dest->size);
				if(!dest->name) {
					ERROR();
				}
				name = dest->name + size;
			}
		}
		name[-1] = 0;
		src->bytes = avail - 1;
		src->next = buff;
		dest->ref = 0;
		return 1;
	}
	
	/* get counts */
	count = ((short unsigned *) (dest->counts)) - 1;
	dest->total = (num = 0);
	while((size = *buff++) != '\n') {
		/* rebuff */
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				return 0;
			}
			buff = src->buffer;
		}
		
		/* get count */
		if(size == '\t') {
			if(num) {
				*count = num;
				dest->total += num;
			}
			*++count = (num = 0);
		} else {
			num = 10 * num + (size - '0');
		}
	}
	if(num) {
		*count = num;
		dest->total += num;
	}
	src->bytes = avail - 1;
	src->next = buff;
	return 1;
}

int FileBuffSkipTemplate(FileBuff *src, NucCount *dest) {
	
	unsigned char *buff;
	int avail;
	int (*buffFileBuff)(FileBuff *);
	
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
	
	/* skip entry */
	while(*buff++ != '#') {
		/* rebuff */
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				return 0;
			}
			buff = src->buffer;
		}
	}
	
	/* get new name */
	src->bytes = avail;
	src->next = buff - 1;
	
	return FileBuffGetRow(src, dest);
}

MatrixCounts * initMat(unsigned matSize, unsigned nameSize) {
	
	MatrixCounts *dest;
	
	dest = smalloc(sizeof(MatrixCounts));
	dest->name = setQseqs(nameSize);
	dest->len = 0;
	dest->size = matSize;
	dest->refs = smalloc(matSize);
	dest->counts = smalloc(7 * matSize * sizeof(short unsigned));
	
	return dest;
}

void setMatName(MatrixCounts *dest, NucCount *src) {
	
	Qseqs *name;
	
	name = dest->name;
	if(name->size < src->size) {
		name->size = src->size;
		free(name->seq);
		name->seq = smalloc(name->size);
	}
	//name->len = sprintf((char *) name->seq, "%s", src->name);
	name->len = strxfrm((char *) name->seq, (char *) src->name, src->size);
}

int FileBuffLoadMat(MatrixCounts *dest, FileBuff *src) {
	
	unsigned char *buff, *refs, c;
	unsigned size, avail, num, total, len;
	short unsigned *counts;
	int (*buffFileBuff)(FileBuff *);
	
	/* init */
	avail = src->bytes;
	buffFileBuff = src->buffFileBuff;
	if(avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
	}
	buff = src->next;
	size = dest->size;
	refs = dest->refs;
	counts = dest->counts - 1;
	len = 0;
	total = (num = 0);
	*refs++ = *buff++;
	if(--avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	while((c = *buff++) != '#') {
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				return 0;
			}
			buff = src->buffer;
		}
		
		if(c == '\n') {
			/* finalize row */
			total += num;
			*counts++ = num;
			*counts = total;
			++len;
			num = 0;
			total = 0;
			/* get ref */
			if((*refs++ = *buff++) == '\n') {
				*--refs = 0;
				dest->len = len;
				src->bytes = avail - 1;
				src->next = buff;
				return 1;
			}
			if(--avail == 0) {
				if((avail = buffFileBuff(src)) == 0) {
					dest->len = len;
					return 0;
				}
				buff = src->buffer;
			}
			if(--size == 0) {
				size = dest->size;
				dest->size <<= 1;
				dest->counts = realloc(dest->counts, 7 * dest->size * sizeof(short unsigned));
				dest->refs = realloc(dest->refs, dest->size);
				if(!dest->counts || !dest->refs) {
					ERROR();
				}
				refs = dest->refs + size;
				counts = dest->counts + 7 * size - 1;
			}
		} else if(c == '\t') {
			if(num) {
				*counts = num;
				total += num;
			}
			*++counts = (num = 0);
		} else {
			num = 10 * num + (c - '0');
		}
	}
	
	dest->len = len;
	src->bytes = avail;
	src->next = buff - 1;
	
	return 1;
}