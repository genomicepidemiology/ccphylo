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
#include "filebuff.h"
#include "meth.h"
#include "methparse.h"
#include "pherror.h"
#include "qseqs.h"

unsigned char * getMethBitTable() {
	
	int i;
	unsigned char *table;
	
	table = smalloc(384); /* 128 * 3 = 384 -> OS independent */
	i = 385;
	--table;
	while(--i) {
		*++table = 32;
	}
	table -= 255;
	table['\n'] = 16;
	
	table['a'] = 0;
	table['c'] = 1;
	table['g'] = 2;
	table['t'] = 3;
	table['u'] = 3;
	/* &8 -> meth site */
	table['A'] = 8;
	table['C'] = 9;
	table['G'] = 10;
	table['T'] = 11;
	table['U'] = 11;
	
	/* unknowns */
	table['n'] = 4;
	table['N'] = 4;
	table['-'] = 4;
	table['R'] = 4;
	table['Y'] = 4;
	table['S'] = 4;
	table['W'] = 4;
	table['K'] = 4;
	table['M'] = 4;
	table['B'] = 4;
	table['D'] = 4;
	table['H'] = 4;
	table['V'] = 4;
	table['X'] = 4;
	table['r'] = 4;
	table['y'] = 4;
	table['s'] = 4;
	table['w'] = 4;
	table['k'] = 4;
	table['m'] = 4;
	table['b'] = 4;
	table['d'] = 4;
	table['h'] = 4;
	table['v'] = 4;
	table['x'] = 4;
	
	return table;
}

int FileBuffgetFsaMethSeq(FileBuff *src, Qseqs *qseq, unsigned char *trans) {
	
	unsigned char *buff, *seq;
	int size, avail, (*buffFileBuff)(FileBuff *);
	
	/* init */
	avail = src->bytes;
	buff = src->next;
	buffFileBuff = src->buffFileBuff;
	qseq->len = 0;
	if(avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	
	/* skip header */
	if(*buff == '>') {
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				return 0;
			}
			buff = src->buffer;
		}
		while(*++buff != '\n') {
			if(--avail == 0) {
				if((avail = buffFileBuff(src)) == 0) {
					return 0;
				}
				buff = src->buffer;
			}
		}
		++buff;
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				return 0;
			}
			buff = src->buffer;
		}
	}
	
	/* get qseq */
	seq = qseq->seq;
	size = qseq->size;
	while(*buff != '>') {
		*seq = trans[*buff++];
		if(((*seq) >> 4) == 0) {
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

MethMotif * qseq2methMotif(Qseqs *qseq) {
	
	unsigned i, len, *mask;
	long unsigned *motif;
	unsigned char *seq;
	MethMotif *dest;
	
	/* init */
	dest = newMethMotif(qseq->len);
	
	/* get motif and mask */
	motif = dest->motif - 1;
	mask = dest->mask - 1;
	seq = qseq->seq;
	len = qseq->len;
	for(i = 0; i < len; ++i) {
		if((i & 31) == 0) {
			*++motif = 0;
			*++mask = 0;
		}
		
		/* update motif and mask */
		*motif <<= 1;
		*motif |= (*seq & 3);
		*mask <<= 1;
		if(*seq & 4) {
			*mask |= 1;
		}		
	}
	
	if(qseq->len & 31) {
		i = 32 - (qseq->len & 31);
		*mask <<= i;
		*motif <<= (i << 1);
	}
	
	return dest;
}

MethMotif * getMethMotifs(FileBuff *infile, Qseqs *qseq) {
	
	unsigned char *trans;
	MethMotif *dest, *node;
	
	/* init */
	dest = 0;
	trans = getMethBitTable();
	
	/* get motifs */
	while(FileBuffgetFsaMethSeq(infile, qseq, trans)) {
		/* convert motif */
		node = qseq2methMotif(qseq);
		
		/* link new motif */
		node->next = dest;
		dest = node;
	}
	
	/* clean */
	free(trans - 128);
	
	return dest;
}
