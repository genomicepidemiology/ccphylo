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

#define _XOPEN_SOURCE 600
#include <limits.h>
#include <stdio.h>
#include "bytescale.h"
#include "fsacmp.h"
#include "matrix.h"
#include "pherror.h"
#include "qseqs.h"
#include "threader.h"

void (*getIncPosPtr)(unsigned*, Qseqs*, Qseqs*, unsigned) = &getIncPos;

unsigned char * get2BitTable(unsigned flag) {
	
	int i;
	unsigned char *to2Bit;
	
	to2Bit = smalloc(384); /* 128 * 3 = 384 -> OS independent */
	i = 385;
	--to2Bit;
	while(--i) {
		*++to2Bit = 32;
	}
	to2Bit -= 255;
	to2Bit['A'] = 0;
	to2Bit['C'] = 1;
	to2Bit['G'] = 2;
	to2Bit['T'] = 3;
	to2Bit['U'] = 3;
	to2Bit['N'] = 4;
	to2Bit['-'] = 4;
	/* include insignificant bases */
	if(flag & 8) {
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

unsigned char * getIupacBitTable(unsigned flag) {
	
	int i;
	unsigned char *to2Bit;
	
	to2Bit = smalloc(384); /* 128 * 3 = 384 -> OS independent */
	i = 385;
	--to2Bit;
	while(--i) {
		*++to2Bit = 32;
	}
	to2Bit -= 255;
	to2Bit['A'] = 0;
	to2Bit['C'] = 1;
	to2Bit['G'] = 2;
	to2Bit['T'] = 3;
	to2Bit['U'] = 3;
	to2Bit['N'] = 4;
	to2Bit['-'] = 5;
	to2Bit['R'] = 6;
	to2Bit['Y'] = 7;
	to2Bit['S'] = 8;
	to2Bit['W'] = 9;
	to2Bit['K'] = 10;
	to2Bit['M'] = 11;
	to2Bit['B'] = 12;
	to2Bit['D'] = 13;
	to2Bit['H'] = 14;
	to2Bit['V'] = 15;
	to2Bit['X'] = 4;
	/* include insignificant bases */
	if(flag & 1) {
		to2Bit['a'] = 4;
		to2Bit['c'] = 4;
		to2Bit['g'] = 4;
		to2Bit['t'] = 4;
		to2Bit['u'] = 4;
		to2Bit['n'] = 4;
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
	} else {
		to2Bit['a'] = 0 | 16;
		to2Bit['c'] = 1 | 16;
		to2Bit['g'] = 2 | 16;
		to2Bit['t'] = 3 | 16;
		to2Bit['u'] = 3 | 16;
		to2Bit['n'] = 4;
		to2Bit['r'] = 6 | 16;
		to2Bit['y'] = 7 | 16;
		to2Bit['s'] = 8 | 16;
		to2Bit['w'] = 9 | 16;
		to2Bit['k'] = 10 | 16;
		to2Bit['m'] = 11 | 16;
		to2Bit['b'] = 12 | 16;
		to2Bit['d'] = 13 | 16;
		to2Bit['h'] = 14 | 16;
		to2Bit['v'] = 15 | 16;
	}
	to2Bit['x'] = 4;
	
	return to2Bit;
}

void initIncPos(unsigned *include, int len) {
	
	int complen;
	
	--include;
	/* convert length to compressed length */
	if(len & 31) {
		complen = (len >> 5) + 2;
	} else {
		complen = (len >> 5) + 1;
	}
	while(--complen) {
		*++include = UINT_MAX;
	}
	*include <<= (32 - (len & 31));
}

void getIncPos(unsigned *include, Qseqs *seq, Qseqs *ref, unsigned proxi) {
	
	int lastSNP;
	unsigned i, end, len, mask, topBit, *includePtr;
	unsigned char *cPtr, *rPtr, c, r;
	
	topBit = UINT_MAX ^ (UINT_MAX >> 1);
	len = seq->len;
	lastSNP = -1;
	cPtr = seq->seq - 1;
	rPtr = ref->seq - 1;
	for(i = 0; i < len; ++i) {
		c = *++cPtr;
		r = *++rPtr;
		/* mask position */
		if(c != r || c == 4 || (c & 16)) {
			/* unknown base */
			if(c == 4 || r == 4) {
				/* mask position */
				include[i >> 5] &= (UINT_MAX ^ (1 << (31 - (i & 31))));
			} else if((c & 16) || (r & 16)) {
				/* mask position */
				include[i >> 5] &= (UINT_MAX ^ (1 << (31 - (i & 31))));
				*cPtr &= 15;
				*rPtr &= 15;
			}
			
			/* check proximity */
			if(i - lastSNP <= proxi) {
				/*
				lastSNP = (lastSNP - proxi) < 0 ? 0 : (lastSNP - proxi);
				end = len < i + proxi + 1 ? len : i + proxi + 1;
				*/
				end = i + 1;
				mask = UINT_MAX ^ (1 << (31 - (lastSNP & 31)));
				includePtr = include + (lastSNP >> 5);
				while(lastSNP < end) {
					if((*includePtr &= mask)) {
						if(mask & 1) {
							mask = (mask >> 1) | topBit;
						} else {
							++includePtr;
							mask = UINT_MAX >> 1;
						}
						++lastSNP;
					} else {
						lastSNP = ((lastSNP >> 5) + 1) << 5;
						++includePtr;
						mask = UINT_MAX >> 1;
					}
				}
			}
			
			/* mark position as the last seen SNP */
			lastSNP = i;
		}
	}
}

void getIncPosInsigPrune(unsigned *include, Qseqs *seq, Qseqs *ref, unsigned proxi) {
	
	int lastSNP;
	unsigned i, end, len, mask, topBit, *includePtr;
	unsigned char *cPtr, *rPtr, c, r;
	
	topBit = UINT_MAX ^ (UINT_MAX >> 1);
	len = seq->len;
	lastSNP = -1;
	cPtr = seq->seq - 1;
	rPtr = ref->seq - 1;
	for(i = 0; i < len; ++i) {
		c = *++cPtr;
		r = *++rPtr;
		/* mask position */
		if(c == 4 || r == 4) {
			/* mask position */
			include[i >> 5] &= (UINT_MAX ^ (1 << (31 - (i & 31))));
		} else if((c & 16) || (r & 16)) {
			/* mask position */
			include[i >> 5] &= (UINT_MAX ^ (1 << (31 - (i & 31))));
			*cPtr &= 15;
			*rPtr &= 15;
		}else if(c != r) {
			/* check proximity */
			if(i - lastSNP <= proxi) {
				/*
				lastSNP = (lastSNP - proxi) < 0 ? 0 : (lastSNP - proxi);
				end = len < i + proxi + 1 ? len : i + proxi + 1;
				*/
				end = i + 1;
				mask = UINT_MAX ^ (1 << (31 - (lastSNP & 31)));
				includePtr = include + (lastSNP >> 5);
				while(lastSNP < end) {
					if((*includePtr &= mask)) {
						if(mask & 1) {
							mask = (mask >> 1) | topBit;
						} else {
							++includePtr;
							mask = UINT_MAX >> 1;
						}
						++lastSNP;
					} else {
						lastSNP = ((lastSNP >> 5) + 1) << 5;
						++includePtr;
						mask = UINT_MAX >> 1;
					}
				}
			}
			
			/* mark position as the last seen SNP */
			lastSNP = i;
		}
	}
}

void getIncPosInsig(unsigned *include, Qseqs *seq, Qseqs *ref, unsigned proxi) {
	
	int lastSNP;
	unsigned i, end, len, mask, topBit, *includePtr;
	unsigned char *cPtr, *rPtr, c, r;
	
	topBit = UINT_MAX ^ (UINT_MAX >> 1);
	len = seq->len;
	lastSNP = -1;
	cPtr = seq->seq - 1;
	rPtr = ref->seq - 1;
	for(i = 0; i < len; ++i) {
		c = *++cPtr;
		r = *++rPtr;
		if(c == 4 || r == 4) {
			/* mask unknown position */
			include[i >> 5] &= (0xFFFFFFFF ^ (1 << (31 - (i & 31))));
		} else if(c != r) {
			/* check proximity */
			if(i - lastSNP <= proxi) {
				/*
				lastSNP = (lastSNP - proxi) < 0 ? 0 : (lastSNP - proxi);
				end = len < i + proxi + 1 ? len : i + proxi + 1;
				*/
				end = i + 1;
				mask = 0xFFFFFFFF ^ (1 << (31 - (lastSNP & 31)));
				includePtr = include + (lastSNP >> 5);
				while(lastSNP < end) {
					if((*includePtr &= mask)) {
						if(mask & 1) {
							mask = (mask >> 1) | topBit;
						} else {
							++includePtr;
							mask = 0xFFFFFFFF >> 1;
						}
						++lastSNP;
					} else {
						lastSNP = ((lastSNP >> 5) + 1) << 5;
						++includePtr;
						mask = 0xFFFFFFFF >> 1;
					}
				}
				/*
				this implementation does not scale well
				while(lastSNP < end) {
					if(include[lastSNP >> 5] & (1 << (31 - (lastSNP & 31)))) {
						include[lastSNP >> 5] ^= (1 << (31 - (lastSNP & 31)));
					}
					++lastSNP;
				}
				*/
			}
			
			/* mark position as the last seen SNP */
			lastSNP = i;
		}
	}
}

void maskProxi(unsigned *include, unsigned *include1, unsigned *include2, long unsigned *seq1, long unsigned *seq2, unsigned len, unsigned proxi) {
	
	int lastSNP, seqlen, proxiMasking;
	unsigned i, j, end, mask, topBit, inc, *includePtr, *includeOrg;
	long unsigned kmer1, kmer2;
	
	/* init */
	proxiMasking = 0;
	topBit = 0xFFFFFFFF ^ (0xFFFFFFFF >> 1);
	includeOrg = include;
	lastSNP = len + proxi;
	seqlen = len;
	/* convert length to compressed length */
	if(len & 31) {
		len = (len >> 5) + 1;
	} else {
		len = (len >> 5);
	}
	include += len;
	include1 += len;
	include2 += len;
	seq1 += len;
	seq2 += len;
	if((i = seqlen) & 31) {
		i = ((i >> 5) + 1) << 5;
	}
	++len;
	while(--len) {
		kmer1 = *--seq1;
		kmer2 = *--seq2;
		inc = *--include1 & *--include2;
		*--include = inc;
		
		/* possiblity for difference(s) */
		if(proxi && inc && kmer1 != kmer2) {
			while(inc) {
				if((inc & 1) && (kmer1 & 3) != (kmer2 & 3)) {
					/* check proximity */
					if(lastSNP - i <= proxi) {
						/* mask prev bases */
						mask = 0xFFFFFFFF ^ (1 << (31 - (i & 31)));
						includePtr = includeOrg + (i >> 5);
						end = lastSNP + 1;
						/* uncomment if past proxi masking is wanted */
						/*
						end = seqlen < lastSNP + proxi ? seqlen : lastSNP + proxi;
						*/
						
						j = i;
						while(j < end) {
							*includePtr &= mask;
							if(*includePtr) {
								if(mask & 1) {
									mask = (mask >> 1) | topBit;
								} else {
									++includePtr;
									mask = 0xFFFFFFFF >> 1;
								}
								++j;
							} else {
								j = ((j >> 5) + 1) << 5;
								++includePtr;
								mask = 0xFFFFFFFF >> 1;
							}
						}
						
						/* uncomment if prior proxi masking is wanted */
						/*
						proxiMasking = 1;
						*/
					} else if(proxiMasking) {
						/* mask bases in front of last SNP */
						j = lastSNP - proxi < 0 ? 0 : lastSNP - proxi;
						mask = 0xFFFFFFFF ^ (1 << (31 - (j & 31)));
						includePtr = includeOrg + (j >> 5);
						while(j < lastSNP) {
							*includePtr &= mask;
							if(*includePtr) {
								if(mask & 1) {
									mask = (mask >> 1) | topBit;
								} else {
									++includePtr;
									mask = 0xFFFFFFFF >> 1;
								}
								++j;
							} else {
								j = ((j >> 5) + 1) << 5;
								++includePtr;
								mask = 0xFFFFFFFF >> 1;
							}
						}
						proxiMasking = 0;
					}
					
					/* mark position as the last seen SNP */
					lastSNP = i;
				}
				
				kmer1 >>= 2;
				kmer2 >>= 2;
				inc >>= 1;
				--i;
			}
			i = (i >> 5) << 5;
		} else {
			i -= 32;
		}
	}
	if(proxiMasking) {
		/* mask bases in front of last SNP */
		j = lastSNP - proxi < 0 ? 0 : lastSNP - proxi;
		mask = 0xFFFFFFFF ^ (1 << (31 - (j & 31)));
		includePtr = includeOrg + (j >> 5);
		while(j < lastSNP) {
			*includePtr &= mask;
			if(*includePtr) {
				if(mask & 1) {
					mask = (mask >> 1) | topBit;
				} else {
					++includePtr;
					mask = 0xFFFFFFFF >> 1;
				}
				++j;
			} else {
				j = ((j >> 5) + 1) << 5;
				++includePtr;
				mask = 0xFFFFFFFF >> 1;
			}
		}
	}
}

int getNpos(unsigned *include, int len) {
	
	unsigned n, inc;
	
	n = 0;
	--include;
	len = (len >> 5) + ((len & 31) ? 2 : 1);
	while(--len) {
		inc = *++include;
		while(inc) {
			n += inc & 1;
			inc >>= 1;
		}
	}
	
	return n;
}

void pseudoAlnPrune(unsigned *include, unsigned char **seqs, int len, int n) {
	
	int i, shifter;
	unsigned char *ref, *refptr, *seqptr;
	unsigned *consensus, *conptr;
	
	if(!len || !n) {
		return;
	} else if(!(consensus = calloc((len >> 5) + ((len & 31) ? 1 : 0), sizeof(unsigned)))) {
		ERROR();
	}
	
	
	/* go over alignments */
	while(n && !(ref = *seqs)) {
		++seqs;
		--n;
	}
	while(--n) {
		conptr = consensus;
		refptr = ref - 1;
		if((seqptr = *++seqs)) {
			--seqptr;
			/* go through alignment */
			for(i = 0, shifter = 31; i < len; ++i, --shifter) {
				if(*++refptr != *++seqptr) {
					*conptr |= (1 << shifter);
				}
				if(!shifter) {
					shifter = 32;
					++conptr;
				}
			}
		}
	}
	
	/* mask inclusion array */
	i = (len >> 5) + ((len & 31) ? 1 : 0) + 1;
	--include;
	conptr = consensus - 1;
	while(--i) {
		*++include &= *++conptr;
	}
	
	/* clean up */
	free(consensus);
}

unsigned fsacmp(long unsigned *seq1, long unsigned *seq2, unsigned *include, int len) {
	
	unsigned dist, inc;
	long unsigned kmer1, kmer2;
	
	dist = 0;
	--include;
	--seq1;
	--seq2;
	/* convert length to compressed length */
	if(len & 31) {
		len = (len >> 5) + 2;
	} else {
		len = (len >> 5) + 1;
	}
	while(--len) {
		kmer1 = *++seq1;
		kmer2 = *++seq2;
		inc = *++include;
		/* at least one difference */
		if(inc && kmer1 != kmer2) {
			while(inc) {
				if((inc & 1) && (kmer1 & 3) != (kmer2 & 3)) {
					++dist;
				}
				kmer1 >>= 2;
				kmer2 >>= 2;
				inc >>= 1;
			}
		}
	}
	
	return dist;
}

long unsigned fsacmpair(long unsigned *seq1, long unsigned *seq2, unsigned *include, int len) {
	
	unsigned dist, n, inc;
	long unsigned kmer1, kmer2;
	
	n = 0;
	dist = 0;
	--include;
	--seq1;
	--seq2;
	/* convert length to compressed length */
	if(len & 31) {
		len = (len >> 5) + 2;
	} else {
		len = (len >> 5) + 1;
	}
	
	while(--len) {
		kmer1 = *++seq1;
		kmer2 = *++seq2;
		inc = *++include;
		/* possiblity for difference(s) */
		if(inc && kmer1 != kmer2) {
			while(inc) {
				if((inc & 1)) {
					++n;
					dist += ((kmer1 & 3) != (kmer2 & 3));
				}
				kmer1 >>= 2;
				kmer2 >>= 2;
				inc >>= 1;
			}
		} else {
			while(inc) {
				n += inc & 1;
				inc >>= 1;
			}
		}
	}
	
	/* return distance in upper 4 bytes, and n in the lower */
	kmer1 = dist;
	kmer1 <<= 32;
	kmer1 |= n;
	
	return kmer1;
}

void printDiff(FILE *outfile, int samplei, int samplej, int nuci, int pos, int nucj) {
	
	static volatile int Lock = 0;
	volatile int *lock = &Lock;
	char bases[4] = "ACGT";
	
	lock(lock);
	fprintf(outfile, "(%d, %d)\t%c%d%c\n", samplei, samplej, bases[nuci], pos, bases[nucj]);
	unlock(lock);
}

unsigned fsacmprint(FILE *outfile, int samplei, int samplej, long unsigned *seq1, long unsigned *seq2, unsigned *include, int len) {
	
	unsigned dist, inc, pos;
	long unsigned kmer1, kmer2;
	
	dist = 0;
	--include;
	--seq1;
	--seq2;
	/* convert length to compressed length */
	if(len & 31) {
		len = (len >> 5) + 2;
	} else {
		len = (len >> 5) + 1;
	}
	pos = 1;
	while(--len) {
		kmer1 = *++seq1;
		kmer2 = *++seq2;
		inc = *++include;
		/* at least one difference */
		if(inc && kmer1 != kmer2) {
			while(inc) {
				if((inc & 1) && (kmer1 & 3) != (kmer2 & 3)) {
					printDiff(outfile, samplei, samplej, (kmer1 & 3), pos, (kmer2 & 3));
					++dist;
				}
				kmer1 >>= 2;
				kmer2 >>= 2;
				inc >>= 1;
				++pos;
			}
		} else {
			pos += 32;
		}
	}
	
	return dist;
}

long unsigned fsacmpairint(FILE *outfile, int samplei, int samplej, long unsigned *seq1, long unsigned *seq2, unsigned *include, int len) {
	
	unsigned dist, n, inc, pos;
	long unsigned kmer1, kmer2;
	
	n = 0;
	dist = 0;
	--include;
	--seq1;
	--seq2;
	/* convert length to compressed length */
	if(len & 31) {
		len = (len >> 5) + 2;
	} else {
		len = (len >> 5) + 1;
	}
	pos = 1;
	while(--len) {
		kmer1 = *++seq1;
		kmer2 = *++seq2;
		inc = *++include;
		/* possiblity for difference(s) */
		if(inc && kmer1 != kmer2) {
			while(inc) {
				if((inc & 1)) {
					++n;
					if(((kmer1 & 3) != (kmer2 & 3))) {
						printDiff(outfile, samplei, samplej, (kmer1 & 3), pos, (kmer2 & 3));
						++dist;
					}
				}
				kmer1 >>= 2;
				kmer2 >>= 2;
				inc >>= 1;
				++pos;
			}
		} else {
			while(inc) {
				n += inc & 1;
				inc >>= 1;
			}
			pos += 32;
		}
	}
	
	/* return distance in upper 4 bytes, and n in the lower */
	kmer1 = dist;
	kmer1 <<= 32;
	kmer1 |= n;
	
	return kmer1;
}

unsigned cmpFsa(Matrix *D, int n, int len, long unsigned **seqs, unsigned char *include, unsigned **includes, unsigned norm, FILE *diffile) {
	
	unsigned i, j, inc, *includesPtr;
	long unsigned **seqi, **seqj, dist;
	unsigned char *includei, *includej;
	double *Dptr, nFactor;
	float *Dfptr;
	short unsigned *Dsptr;
	unsigned char *Dbptr;
	
	/* init */
	includesPtr = *includes;
	inc = getNpos(includesPtr, len);
	if(norm) {
		nFactor = norm;
		nFactor /= inc;
	} else {
		nFactor = 1.0;
	}
	
	/* seek first included sample */
	while(*include == 0) {
		++seqs;
		++include;
		if(!--n) {
			return 0;
		}
	}
	seqi = seqs;
	Dptr = 0;
	Dfptr = 0;
	Dsptr = 0;
	Dbptr = 0;
	if(D->mat) {
		Dptr = *(D->mat) - 1;
	} else if(D->fmat) {
		Dfptr = *(D->fmat) - 1;
	} else if(D->smat) {
		Dsptr = *(D->smat) - 1;
	} else {
		Dbptr = *(D->bmat) - 1;
	}
	includei = include;
	
	/* get distances */
	if(diffile) {
		D->n = 1;
		for(i = 0; i < n; ++i) {
			++seqi;
			if(*++includei) {
				++D->n;
				seqj = seqs - 1;
				includej = include - 1;
				for(j = 0; j < i; ++j) {
					if(*++includej) {
						dist = fsacmprint(diffile, i, j, *seqi, *++seqj, includesPtr, len);
						if(Dptr) {
							*++Dptr = nFactor * dist;
						} else if(Dfptr) {
							*++Dfptr = nFactor * dist;
						} else if(Dsptr) {
							*++Dsptr = dtouc(nFactor * dist, 0.5);
						} else {
							*++Dbptr = dtouc(nFactor * dist, 0.5);
						}
					} else {
						++seqj;
					}
				}
			}
		}
	} else {
		i = n;
		n = 1;
		while(--i) {
			++seqi;
			if(*++includei) {
				seqj = seqs - 1;
				includej = include - 1;
				j = ++n;
				while(--j) {
					if(*++includej) {
						dist = fsacmp(*seqi, *++seqj, includesPtr, len);
						if(Dptr) {
							*++Dptr = nFactor * dist;
						} else if(Dfptr) {
							*++Dfptr = nFactor * dist;
						} else if(Dsptr) {
							*++Dsptr = dtouc(nFactor * dist, 0.5);
						} else {
							*++Dbptr = dtouc(nFactor * dist, 0.5);
						}
					} else {
						++j;
						++seqj;
					}
				}
			}
		}
		D->n = n;
	}
	
	return inc;
}

void cmpairFsa(Matrix *D, Matrix *N, int n, int len, long unsigned **seqs, unsigned char *include, unsigned **includes, unsigned norm, unsigned minLength, double minCov, unsigned proxi, FILE *diffile) {
	
	unsigned i, j, inc, *includePair, **includesi, **includesj;
	long unsigned **seqi, **seqj, dist;
	unsigned char *includei, *includej;
	double *Dptr, *Nptr;
	float *Dfptr, *Nfptr;
	short unsigned *Dsptr, *Nsptr;
	unsigned char *Dbptr, *Nbptr;
	
	/* init */
	minLength = minLength < minCov * len ? minCov * len : minLength;
	Dptr = 0;
	Nptr = 0;
	Dfptr = 0;
	Nfptr = 0;
	Dsptr = 0;
	Nsptr = 0;
	Dbptr = 0;
	Nbptr = 0;
	if(D->mat) {
		Dptr = *(D->mat) - 1;
		Nptr = N ? *(N->mat) - 1 : 0;
	} else if(D->fmat) {
		Dfptr = *(D->fmat) - 1;
		Nfptr = N ? *(N->fmat) - 1 : 0;
	} else if(D->smat) {
		Dsptr = *(D->smat) - 1;
		Nsptr = N ? *(N->smat) - 1 : 0;
	} else {
		Dbptr = *(D->bmat) - 1;
		Nbptr = N ? *(N->bmat) - 1 : 0;
	}
	
	/* seek first incuded sample */
	while(n && *include == 0) {
		++seqs;
		++includes;
		++include;
		--n;
	}
	seqi = seqs;
	includesi = includes;
	includei = include;
	includePair = smalloc((len / 32 + 1) * sizeof(unsigned));
	
	/* get distances */
	if(diffile) {
		D->n = 1;
		for(i = 0; i < n; ++i) {
			++seqi;
			++includesi;
			if(*++includei) {
				++D->n;
				seqj = seqs - 1;
				includesj = includes - 1;
				includej = include - 1;
				for(j = 0; j < i; ++j) {
					if(*++includej) {
						/* mask out proximity SNPs */
						maskProxi(includePair, *includesi, *++includesj, *seqi, *++seqj, len, proxi);
						
						/* get distance */
						dist = fsacmpairint(diffile, i, j, *seqi, *seqj, includePair, len);
						
						/* separate distance and included bases */
						if(Dptr) {
							if(minLength <= (inc = dist & UINT_MAX)) {
								if(norm) {
									*++Dptr = (dist >> 32);
								} else {
									*++Dptr = (dist >> 32) * norm;
									*Dptr /= inc;
								}
							} else {
								*++Dptr = dtouc(-1.0, 0);
							}
							if(N) {
								*++Nptr = inc;
							}
						} else if(Dfptr) {
							if(minLength <= (inc = dist & UINT_MAX)) {
								if(norm) {
									*++Dfptr = (dist >> 32);
								} else {
									*++Dfptr = (dist >> 32) * norm;
									*Dfptr /= inc;
								}
							} else {
								*++Dfptr = dtouc(-1.0, 0);
							}
							if(N) {
								*++Nfptr = inc;
							}
						} else if(Dsptr) {
							if(minLength <= (inc = dist & UINT_MAX)) {
								if(norm) {
									*++Dsptr = dtouc((dist >> 32), 0.5);
								} else {
									*++Dsptr = (dtouc((dist >> 32) * norm, 0.5)) / inc;
								}
							} else {
								*++Dsptr = 65535;
							}
							if(N) {
								*++Nsptr = dtouc(inc, 0.5);
							}
						} else {
							if(minLength <= (inc = dist & UINT_MAX)) {
								if(norm) {
									*++Dbptr = dtouc((dist >> 32), 0.5);
								} else {
									*++Dbptr = (dtouc((dist >> 32) * norm, 0.5)) / inc;
								}
							} else {
								*++Dbptr = 255;
							}
							if(N) {
								*++Nbptr = dtouc(inc, 0.5);
							}
						}
					} else {
						++seqj;
						++includesj;
					}
				}
			}
		}
	} else {
		i = n;
		n = 1;
		while(--i) {
			++seqi;
			++includesi;
			if(*++includei) {
				seqj = seqs - 1;
				includesj = includes - 1;
				includej = include - 1;
				j = ++n;
				while(--j) {
					if(*++includej) {
						/* mask out proximity SNPs */
						maskProxi(includePair, *includesi, *++includesj, *seqi, *++seqj, len, proxi);
						
						/* get distance */
						dist = fsacmpair(*seqi, *seqj, includePair, len);
						
						/* separate distance and included bases */
						if(Dptr) {
							if(minLength <= (inc = dist & UINT_MAX)) {
								if(norm) {
									*++Dptr = (dist >> 32);
								} else {
									*++Dptr = (dist >> 32) * norm;
									*Dptr /= inc;
								}
							} else {
								*++Dptr = -1.0;
							}
							if(N) {
								*++Nptr = inc;
							}
						} else if(Dfptr) {
							if(minLength <= (inc = dist & UINT_MAX)) {
								if(norm) {
									*++Dfptr = (dist >> 32);
								} else {
									*++Dfptr = (dist >> 32) * norm;
									*Dfptr /= inc;
								}
							} else {
								*++Dfptr = -1.0;
							}
							if(N) {
								*++Nfptr = inc;
							}
						} else if(Dsptr) {
							if(minLength <= (inc = dist & UINT_MAX)) {
								if(norm) {
									*++Dsptr = dtouc((dist >> 32), 0.5);
								} else {
									*++Dsptr = (dtouc((dist >> 32) * norm, 0.5)) / inc;
								}
							} else {
								*++Dsptr = dtouc(-1.0, 0);
							}
							if(N) {
								*++Nsptr = dtouc(inc, 0.5);
							}
						} else {
							if(minLength <= (inc = dist & UINT_MAX)) {
								if(norm) {
									*++Dbptr = dtouc((dist >> 32), 0.5);
								} else {
									*++Dbptr = (dtouc((dist >> 32) * norm, 0.5)) / inc;
								}
							} else {
								*++Dbptr = dtouc(-1.0, 0);
							}
							if(N) {
								*++Nbptr = dtouc(inc, 0.5);
							}
						}
					} else {
						++seqj;
						++includesj;
						++j;
					}
				}
			}
		}
		D->n = n;
	}
	
	if(N) {
		N->n = D->n;
	}
	free(includePair);
}
