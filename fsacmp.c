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

#include "fsacmp.h"
#include "pherror.h"
#include "qseqs.h"

unsigned char * get2BitTable(unsigned flag) {
	
	int i;
	unsigned char *to2Bit;
	
	to2Bit = smalloc(384); /* 128 * 3 = 384 -> OS independent */
	i = 385;
	--to2Bit;
	while(--i) {
		*++to2Bit = 8;
	}
	to2Bit -= 256;
	to2Bit['\n'] = 16;
	to2Bit['A'] = 0;
	to2Bit['C'] = 1;
	to2Bit['G'] = 2;
	to2Bit['T'] = 3;
	to2Bit['U'] = 3;
	to2Bit['N'] = 4;
	to2Bit['-'] = 4;
	/* include insignificant bases */
	if(flag & 1) {
		to2Bit['a'] = 0;
		to2Bit['c'] = 1;
		to2Bit['g'] = 2;
		to2Bit['t'] = 3;
		to2Bit['u'] = 3;
		to2Bit['n'] = 4;
	} else {
		to2Bit['a'] = 4;
		to2Bit['c'] = 4;
		to2Bit['g'] = 4;
		to2Bit['t'] = 4;
		to2Bit['u'] = 4;
		to2Bit['n'] = 4;
	}
	to2Bit['R'] = 4;
	to2Bit['Y'] = 4;
	to2Bit['S'] = 4;
	to2Bit['W'] = 4;
	to2Bit['K'] = 4;
	to2Bit['M'] = 4;
	to2Bit['B'] = 4;
	to2Bit['D'] = 4;
	to2Bit['H'] = 4;
	to2Bit['V'] = 4;
	to2Bit['X'] = 4;
	to2Bit['r'] = 4;
	to2Bit['y'] = 4;
	to2Bit['s'] = 4;
	to2Bit['w'] = 4;
	to2Bit['k'] = 4;
	to2Bit['m'] = 4;
	to2Bit['b'] = 4;
	to2Bit['d'] = 4;
	to2Bit['h'] = 4;
	to2Bit['v'] = 4;
	to2Bit['x'] = 4;
	
	return to2Bit;
}

void getMethPos(short unsigned *include, Qseqs *ref) {
	
	/* here */
	
	/* ask Malte how he did it */
	
	
	
}

void getIncPos(unsigned *include, Qseqs *seq, Qseqs *ref, unsigned proxi) {
	
	int i, lastSNP;
	short unsigned *iPtr, inc;
	unsigned char *cPtr, *rPtr, c, r;
	
	lastSNP = -1;
	cPtr = seq->seq - 1;
	rPtr = ref->seq - 1;
	for(i = 0; i < len; ++i) {
		c = ++*cPtr;
		r = ++*rPtr;
		
		/* here */
		/* incorporate proxi */
		
		/* mask position */
		if(c != r || c == 4) {
			include[i >> 5] |= 1 << i & 31;
		}
	}
	
	
	
	
	
}

void getIncPos_old(unsigned *include, Qseqs *seq, Qseqs *ref, unsigned proxi) {
	
	int i, inc, n, lastSNP, nProxi, pos;
	unsigned char *ptr, *refPtr, *proxiPtr, c, r;
	
	/* fill in proximity exclusions */
	lastSNP = -1;
	n = 0;
	--include;
	ptr = seq->seq - 1;
	refPtr = ref->seq - 1;
	i = seq->len + 1;
	while(--i) {
		if((c = *++ptr) == 4) {
			inc = 0;
		}
		r = *++refPtr;
		if(c != r) {
			/* check proximity */
			if(proxiPtr) {
				/* already in proximity */
				pos = seq->len - i;
				if(proxi < i) { /* check proxy exceeds end of seq */
					nProxi += (pos - lastSNP);
				} else {
					/* break the loop */
					nProxi = i + 1;
					i = 1;
				}
				lastSNP = pos;
			} else if(seq->len - i - lastSNP < proxi) {
				/* new proximity */
				pos = lastSNP - proxi;
				if(pos < 0) {
					nProxi = proxi + pos;
					pos = 0;
				} else {
					nProxi = proxi;
				}
				proxiPtr = seq->seq + pos;
				*proxiPtr = 0;
				lastSNP = seq->len - i;
				if(proxi < i) { /* check proxy exceeds end of seq */
					nProxi += proxi + 1;
				} else {
					/* break the loop */
					nProxi = i + 1;
					i = 1;
				}
			}
			lastSNP = seq->len - i;
		}
		if() {
			*++include = inc;
			inc = 0;
		}
		/* clear proximity */
		if(proxiPtr) {
			if(--nProxi == 0) {
				proxiPtr = 0;
			} else {
				*++proxiPtr = 0;
			}
		}
	}
	
	/* clear remaining proxi */
	if(proxiPtr) {
		while(--nProxi) {
			*++proxiPtr = 0;
		}
	}
}

int getNpos(unsigned *include, int len) {
	
	unsigned n, inc;
	
	n = 0;
	--inlcude;
	len = (len >> 5) + 2; /* convert length to compressed length */
	while(--len) {
		if((inc = *++include)) {
			while(inc) {
				n += inc & 1;
				inc >>= 1;
			}
		}
	}
	
	return n;
}

unsigned fsacmp(long unsigned *seq1, long unsigned *seq2, unsigned *include, int len) {
	
	unsigned dist, inc;
	long unsigned kmer1, kmer2;
	
	dist = 0;
	--include;
	--seq1;
	--seq2;
	len = (len >> 5) + 2; /* convert length to compressed length */
	while(--len) {
		kmer1 = *++seq1;
		kmer2 = *++seq2;
		/* at least one difference */
		if(*++include && kmer1 != kmer2) {
			inc = *include;
			while(inc) {
				if((inc & 1) && (kmer1 & 3) != (kmer2 & 3)) {
					++dist;
				}
				kmer1 >>= 2;
				kmer2 >>= 2;
				inc >>= 1;
			}
		} else {
			++include;
		}
	}
	
	return dist;
}
