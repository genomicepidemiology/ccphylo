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

#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include "filebuff.h"
#include "pherror.h"
#include "qseqs.h"
#include "seqparse.h"

int FileBuffgetFsa(FileBuff *src, Qseqs *header, Qseqs *qseq, unsigned char *trans) {
	
	unsigned char *buff, *seq;
	unsigned size, avail;
	int (*buffFileBuff)(FileBuff *);
	
	/* init */
	avail = src->bytes;
	buff = src->next;
	buffFileBuff = src->buffFileBuff;
	header->len = 0;
	*(header->seq) = 0;
	qseq->len = 0;
	*(qseq->seq) = 0;
	if(avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	
	/* get header */
	seq = header->seq;
	size = header->size;
	while((*seq++ = *buff++) != '\n') {
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				return 0;
			}
			buff = src->buffer;
		}
		if(--size == 0) {
			size = header->size;
			header->size <<= 1;
			header->seq = realloc(header->seq, header->size);
			if(!header->seq) {
				ERROR();
			}
			seq = header->seq + size;
		}
	}
	if(--avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	/* chomp header */
	while(isspace(*--seq)) {
		++size;
	}
	*++seq = 0;
	header->len = header->size - size + 1;
	/* get qseq */
	seq = qseq->seq;
	size = qseq->size;
	while(*buff != '>') {
		*seq = trans[*buff++];
		if(((*seq) >> 3) == 0) {
			if(--size == 0) {
				size = qseq->size;
				qseq->size <<= 1;
				qseq->seq = realloc(qseq->seq, qseq->size);
				if(!qseq->seq) {
					ERROR();
				}
				seq = qseq->seq + size;
			} else {
				++seq;
			}
		}
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				/* chomp header */
				while(*--seq == 8) {
					++size;
				}
				*++seq = 0;
				qseq->len = qseq->size - size;
				
				src->bytes = 0;
				src->next = buff;
				return 1;
			}
			buff = src->buffer;
		}
	}
	/* chomp header */
	while(*--seq == 8) {
		++size;
	}
	*++seq = 0;
	qseq->len = qseq->size - size;
	
	src->bytes = avail;
	src->next = buff;
	
	return 1;
}

int FileBuffgetFsaHeader(FileBuff *src, Qseqs *header) {
	
	unsigned char *buff, *seq;
	unsigned size, avail;
	int (*buffFileBuff)(FileBuff *);
	
	/* init */
	avail = src->bytes;
	buff = src->next;
	buffFileBuff = src->buffFileBuff;
	header->len = 0;
	*(header->seq) = 0;
	if(avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	
	/* find next header */
	while(*buff++ != '>') {
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				return 0;
			}
			buff = src->buffer;
		}
	}
	if(--avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	
	/* get header */
	seq = header->seq;
	size = header->size;
	while((*seq++ = *buff++) != '\n') {
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				return 0;
			}
			buff = src->buffer;
		}
		if(--size == 0) {
			size = header->size;
			header->size <<= 1;
			header->seq = realloc(header->seq, header->size);
			if(!header->seq) {
				ERROR();
			}
			seq = header->seq + size;
		}
	}
	while(isspace(*--seq)) {
		++size;
	}
	*++seq = 0;
	header->len = header->size - size + 1;
	
	src->bytes = avail - 1;
	src->next = buff;
	
	return 1;
}

int FileBuffgetFsaSeq(FileBuff *src, Qseqs *qseq, unsigned char *trans) {
	
	unsigned char *buff, *seq;
	unsigned size, avail;
	int (*buffFileBuff)(FileBuff *);
	
	/* init */
	avail = src->bytes;
	buff = src->next;
	buffFileBuff = src->buffFileBuff;
	qseq->len = 0;
	*(qseq->seq) = 0;
	if(avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	
	/* get qseq */
	seq = qseq->seq;
	size = qseq->size;
	while(*buff != '>') {
		*seq = trans[*buff++];
		if(*seq < 32) {
			if(--size == 0) {
				size = qseq->size;
				qseq->size <<= 1;
				qseq->seq = realloc(qseq->seq, qseq->size);
				if(!qseq->seq) {
					ERROR();
				}
				seq = qseq->seq + size;
			} else {
				++seq;
			}
		}
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				qseq->len = qseq->size - size;
				src->bytes = 0;
				src->next = buff;
				return 1;
			}
			buff = src->buffer;
		}
	}
	qseq->len = qseq->size - size;
	
	src->bytes = avail;
	src->next = buff;
	
	return 1;
}
