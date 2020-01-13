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
#include "meth.h"
#include "pherror.h"

MethMotif * newMethMotif(int num, int len) {
	
	int size;
	MethMotif *dest;
	
	size = num * ((len << 5) + (len & 31));
	dest = smalloc(sizeof(MethMotif) + size * (sizeof(long unsigned) + sizeof(unsigned)));
	dest->num = num;
	dest->len = len;
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
	long unsigned match, mer;
	
	/* mark unmatched parts */
	mer = *motif ^ seq;
	/* mirror nibbles, and mask them out */
	match = mer | ((mer << 1) & 0xAAAAAAAAAAAAAAAA) | ((mer >> 1) & 0x5555555555555555);
	while(--num && match) {
		mer = *++motif ^ seq;
		match &= mer | ((mer << 1) & 0xAAAAAAAAAAAAAAAA) | ((mer >> 1) & 0x5555555555555555);
	}
	
	return match;
}

long unsigned getMmer(long unsigned *seq, unsigned cPos, const unsigned shifter, const long unsigned mask) {
	
	unsigned iPos = (cPos & 31) << 1;
	cPos >>= 5;
	
	return (iPos <= shifter) ? ((seq[cPos] << iPos) & mask) : (((seq[cPos] << iPos) | (seq[cPos + 1] >> (64-iPos))) & mask);
}

int matchMotif(long unsigned *seq, int pos, int seqlen, MethMotif *motif) {
	
	int sPos, mPos, len, num, shifter;
	long unsigned kmer, mask, match, *motifP;
	
	len = motif->len;
	num = motif->num;
	shifter = len <= 32 ? 64 - (len << 1) : 0;
	mask = 0xFFFFFFFFFFFFFFFF << shifter;
	motifP = motif->motif;
	seqlen -= len - 1;
	--pos;
	
	while(++pos < seqlen) {
		kmer = getMmer(seq, pos, shifter, mask);
		if(!matchMotif32(kmer, motifP, num)) {
			motifP += num;
			sPos = pos + 32;
			mPos = len <= 32 ? 0 : 32;
			match = 0;
			while(!match && mPos) {
				mPos += 32;
				if(len <= mPos) {
					shifter = 64 - (((mPos - len) + 32) << 1);
					mask = 0xFFFFFFFFFFFFFFFF << shifter;
					mPos = 0;
				}
				kmer = getMmer(seq, sPos, shifter, mask);
				match = matchMotif32(kmer, motifP, num);
				sPos += 32;
				motifP += num;
			}
			
			if(match) {
				shifter = len <= 32 ? 64 - (len << 1) : 0;
				mask = 0xFFFFFFFFFFFFFFFF << shifter;
				motifP = motif->motif;
			} else {
				return pos;
			}
		}
	}
	
	return -1;
}

void maskMotif(unsigned *include, int pos, unsigned *mask, int len) {
	
	unsigned mShifter, rShifter;
	
	mShifter = pos & 31;
	rShifter = 32 - mShifter;
	include += pos >> 5;
	--mask;
	while(0 < len) {
		*include &= 0xFFFFFFFF ^ (*++mask >> mShifter);
		if(mShifter && rShifter < len) {
			*++include &= 0xFFFFFFFF ^ ((*mask << rShifter));
		} else {
			++include;
		}
		len -= 32;
	}
}

int maskMotifs(long unsigned *seq, unsigned *include, int seqlen, MethMotif *motif) {
	
	/* returns number of found methylation sites */
	int pos, n;
	
	/* search motifs */
	n = 0;
	while(motif) {
		/* search seq */
		pos = -1;
		while(0 <= (pos = matchMotif(seq, pos + 1, seqlen, motif))) {
			maskMotif(include, pos, motif->mask, motif->len);
			++n;
		}
		motif = motif->next;
	}
	
	return n;
}
