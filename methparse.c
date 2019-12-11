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
		*++table = 64;
	}
	table -= 255;
	table['\n'] = 32;
	table['-'] = 32;
	table['.'] = 32;
	
	/* standard bases */
	table['a'] = 1;
	table['c'] = 2;
	table['g'] = 4;
	table['t'] = 8;
	table['u'] = 8;
	table['r'] = 5;
	table['y'] = 10;
	table['s'] = 6;
	table['w'] = 9;
	table['k'] = 12;
	table['m'] = 3;
	table['b'] = 14;
	table['d'] = 13;
	table['h'] = 11;
	table['v'] = 7;
	table['x'] = 15;
	table['n'] = 15;
	
	/* &16 -> meth site */
	table['A'] = 17;
	table['C'] = 18;
	table['G'] = 20;
	table['T'] = 24;
	table['U'] = 24;
	table['R'] = 21;
	table['Y'] = 26;
	table['S'] = 22;
	table['W'] = 25;
	table['K'] = 28;
	table['M'] = 19;
	table['B'] = 30;
	table['D'] = 29;
	table['H'] = 27;
	table['V'] = 23;
	table['X'] = 31;
	table['N'] = 31;
	
	return table;
}

void strrcMeth(unsigned char *qseq, int q_len) {
	
	int seqlen;
	unsigned char carry, *rseq;
	unsigned char comp[] = {0, 8, 4, 12, 2, 10, 6, 14, 1, 9, 5, 13, 3, 11, 7, 15, 16, 24, 20, 28, 18, 26, 22, 30, 17, 25, 21, 29, 19, 27, 23, 31};
	
	seqlen = (q_len >> 1) + 1;
	rseq = qseq + q_len;
	--qseq;
	while(--seqlen) {
		carry = comp[*++qseq];
		*qseq = comp[*--rseq];
		*rseq = carry;
	}
	if(q_len & 1) {
		*qseq = comp[*qseq];
	}
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
		if((*seq = trans[*buff++]) < 32) {
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
	
	unsigned i, len, num, *mask;
	long unsigned *motif, *mPtr;
	unsigned char base, mBase, *seq;
	unsigned char nums[] = {0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4};
	unsigned char bases[] = {0, 0, 1, 0, 2, 0, 1, 0, 3, 0, 1, 0, 2, 0, 1, 0};
	MethMotif *dest;
	
	/* init */
	len = qseq->len;
	i = len;
	seq = qseq->seq;
	num = nums[*seq];
	while(--i) {
		if(num < nums[*++seq]) {
			num = nums[*seq];
		}
	}
	dest = newMethMotif(num, len);
	
	/* get motif and mask */
	motif = dest->motif - dest->num;
	mask = dest->mask - 1;
	seq = qseq->seq - 1;
	for(i = 0; i < len; ++i) {
		if((i & 31) == 0) {
			motif += dest->num;
			*++mask = 0;
		}
		
		/* update mask */
		base = *++seq;
		if(base & 16) {
			base ^= 16;
			*mask = (*mask << 1) | 1;
		} else {
			*mask <<= 1;
		}
		
		/* update motif */
		num = nums[base] + 1;
		mPtr = motif - 1;
		while(--num) {
			*++mPtr <<= 2;
			mBase = bases[base];
			*mPtr |= mBase;
			base ^= (1 << mBase);
		}
		/* update remainder */
		base = *seq & 31;
		num = dest->num - nums[base] + 1;
		mBase = bases[base];
		while(--num) {
			*++mPtr <<= 2;
			*mPtr |= mBase;
		}
	}
	
	/* fence post */
	if(len & 31) {
		i = 32 - (len & 31);
		*mask <<= i;
		i <<= 1;
		num = dest->num + 1;
		mPtr = motif - 1;
		while(--num) {
			*++mPtr <<= i;
		}
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
		
		/* rc */
		strrcMeth(qseq->seq, qseq->len);
		node = qseq2methMotif(qseq);
		node->next = dest;
		dest = node;
	}
	
	/* clean */
	free(trans - 128);
	
	return dest;
}
