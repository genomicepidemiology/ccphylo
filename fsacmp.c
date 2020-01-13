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
#include "fsacmp.h"
#include "matrix.h"
#include "pherror.h"
#include "qseqs.h"
#include "threader.h"

unsigned char * get2BitTable(unsigned flag) {
	
	int i;
	unsigned char *to2Bit;
	
	to2Bit = smalloc(384); /* 128 * 3 = 384 -> OS independent */
	i = 385;
	--to2Bit;
	while(--i) {
		*++to2Bit = 8;
	}
	to2Bit -= 255;
	to2Bit['\n'] = 16;
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
		if(c != r || c == 4) {
			/* unknown base */
			if(c == 4 || r == 4) {
				/* mask position */
				include[i >> 5] &= (UINT_MAX ^ (1 << (31 - (i & 31))));
			}
			
			/* check proximity */
			if(i - lastSNP < proxi) {
				lastSNP -= proxi;
				if(lastSNP < 0) {
					lastSNP = 0;
				}
				end = i + proxi;
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
	
	int lastSNP;
	unsigned i, j, end, mask, topBit, inc, *includePtr;
	long unsigned kmer1, kmer2;
	
	/* init */
	topBit = UINT_MAX ^ (UINT_MAX >> 1);
	lastSNP = -1;
	i = 0;
	--include;
	--include1;
	--include2;
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
		inc = *++include1 & *++include2;
		*++include = inc;
		
		/* possiblity for difference(s) */
		if(inc && kmer1 != kmer2 && proxi) {
			j = i;
			while(inc) {
				if((inc & 1) && (kmer1 & 3) != (kmer2 & 3)) {
					/* check proximity */
					if(j - lastSNP < proxi) {
						lastSNP -= proxi;
						if(lastSNP < 0) {
							lastSNP = 0;
						}
						end = j + proxi;
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
					lastSNP = j;
				}
				kmer1 >>= 2;
				kmer2 >>= 2;
				inc >>= 1;
				++j;
			}
		}
		i += 32;
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
	
	static volatile int lock[1] = {0};
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
	Dptr = *(D->mat) - 1;
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
						*++Dptr = nFactor * dist;
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
						*++Dptr = nFactor * dist;
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
	
	/* init */
	minLength = minLength < minCov * len ? minCov * len : minLength;
	Dptr = *(D->mat) - 1;
	Nptr = N ? *(N->mat) - 1 : 0;
	
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
