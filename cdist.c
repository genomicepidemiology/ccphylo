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

#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cdist.h"
#include "fbseek.h"
#include "filebuff.h"
#include "fsacmp.h"
#include "fsacmpthrd.h"
#include "matrix.h"
#include "meth.h"
#include "pherror.h"
#include "phy.h"
#include "seqparse.h"
#define exchange(src1, src2, tmp) tmp = src1; src1 = src2; src2 = tmp;

int ltdFsaMatrix_get(Matrix *D, Matrix *N, int numFile, long unsigned **seqs, int cSize, FileBuff *infile, TimeStamp **targetStamps, unsigned char *include, unsigned **includes, char *targetTemplate, char **filenames, unsigned char *trans, Qseqs *ref, Qseqs *seq, Qseqs *header, unsigned norm, unsigned minLength, double minCov, unsigned flag, unsigned proxi, MethMotif *motif, FILE *diffile, int tnum) {
	
	unsigned i, j, pair, len, inludeN, **includesPtr;
	long unsigned **seqsPtr;
	unsigned char *includePtr, *tmpseq;
	TimeStamp **targetStamp;
	
	/* init */
	pair = flag & 2 ? 1 : 0;
	len = 0;
	ref->len = 0;
	inludeN = numFile;
	includesPtr = includes;
	includePtr = include - 1;
	seqsPtr = seqs - 1;
	targetStamp = targetStamps - 1;
	--filenames;
	
	/* load sequences, get included positions and compress */
	i = numFile + 1;
	while(--i) {
		if(*++includePtr) {
			/* open and validate file */
			if(openAndDetermine(infile, *++filenames) != '>') {
				fprintf(stderr, "\"%s\" is not fasta.\n", *filenames);
				exit(1);
			}
			
			/* find the correct entry */
			if(*++targetStamp) {
				seekFileBiff(infile, *targetStamp);
			}
			while(FileBuffgetFsaHeader(infile, header) && strcmp((char *) header->seq, targetTemplate) != 0);
			
			/* make timestamp */
			*targetStamp = timeStampFileBuff(infile, *targetStamp);
			
			/* get sequence */
			if(!header->len) {
				fprintf(stderr, "Missing template entry (\"%s\") in file:\t%s\n", targetTemplate, *filenames);
				*includePtr = 0;
				++seqsPtr;
				--inludeN;
			} else if(FileBuffgetFsaSeq(infile, seq, trans)) {
				if(ref->len) {
					if(seq->len != ref->len) {
						fprintf(stderr, "Sequences does not match: %s\n", *filenames);
						exit(1);
					}
					
					/* make / update inclusion array(s) */
					if(pair) {
						initIncPos(*includesPtr, len);
						qseq2nibble(seq, *++seqsPtr);
						maskMotifs(*seqsPtr, *includesPtr, len, motif);
						getIncPosPtr(*includesPtr, seq, seq, proxi);
						if(getNpos(*includesPtr, len) < minLength) {
							fprintf(stderr, "Template (\"%s\") did not exceed threshold for inclusion:\t%s\n", targetTemplate, *filenames);
							*includePtr = 0;
							--inludeN;
						}
					} else {
						qseq2nibble(seq, *++seqsPtr);
						maskMotifs(*seqsPtr, *includes, len, motif);
						getIncPosPtr(*includes, seq, ref, proxi);
					}
				} else {
					len = seq->len;
					minLength = minLength < minCov * len ? minCov * len : minLength;
					if(cSize < len / 32 + 1) {
						cSize = len / 32 + 1;
						j = numFile;
						while(--j) {
							if(!(seqs[j] = realloc(seqs[j], cSize * sizeof(long unsigned)))) {
								ERROR();
							}
							if(pair && (includes[j] = realloc(includes[j], cSize * sizeof(unsigned))) == 0) {
								ERROR();
							}
						}
						if(!(*seqs = realloc(*seqs, cSize * sizeof(long unsigned)))) {
							ERROR();
						}
						if(!(*includes = realloc(*includes, cSize * sizeof(unsigned)))) {
							ERROR();
						}
					}
					initIncPos(*includesPtr, len);
					qseq2nibble(seq, *++seqsPtr);
					maskMotifs(*seqsPtr, *includesPtr, len, motif);
					getIncPosPtr(*includesPtr, seq, seq, proxi);
					
					if(getNpos(*includesPtr, len) < minLength) {
						fprintf(stderr, "Template (\"%s\") did not exceed threshold for inclusion:\t%s\n", targetTemplate, *filenames);
						*includePtr = 0;
						ref->len = 0;
						--inludeN;
					} else {
						/* swap seq and ref */
						exchange(seq->size, ref->size, len);
						exchange(seq->len, ref->len, len);
						exchange(seq->seq, ref->seq, tmpseq);
					}
				}
			} else {
				fprintf(stderr, "Missing template sequence (\"%s\") in file:\t%s\n", targetTemplate, *filenames);
				*includePtr = 0;
				++seqsPtr;
				--inludeN;
			}
			closeFileBuff(infile);
		} else {
			++filenames;
			++targetStamp;
			++seqsPtr;
			--inludeN;
		}
		includesPtr += pair;
	}
	
	/* adjust number of threads */
	if(numFile * (numFile - 1) / 2 < tnum) {
		fprintf(stderr, "Adjustning number of nodes to %d, to conform with the matrix size.\n", (tnum = numFile * (numFile - 1) / 2));
	}
	
	/* make ltd matrix */
	if(!inludeN) {
		D->n = 0;
		errno = 1;
		fprintf(stderr, "All sequences were trimmed away.\n");
	} else if(pair) {
		fsaCmpThreadOut(tnum, &cmpairFsaThrd, D, N, numFile, len, seqs, include, includes, norm, minLength, minCov, diffile, targetTemplate, ref, 0, proxi);
		//cmpairFsa(D, N, numFile, len, seqs, include, includes, norm, minLength, minCov, proxi, diffile);
	} else if(minLength <= getNpos(*includes, len)) {
		fsaCmpThreadOut(tnum, &cmpFsaThrd, D, N, numFile, len, seqs, include, includes, norm, minLength, minCov, diffile, targetTemplate, ref, 0, proxi);
		/* i = cmpFsa(D, numFile, len, seqs, include, includes, norm, diffile);
		fprintf(stderr, "# %d / %d bases included in distance matrix.\n", i, len); */
	} else {
		D->n = 0;
		errno = 1;
		fprintf(stderr, "No sufficient overlap was found.\n");
	}
	
	return cSize;
}

int ltdMsaMatrix_get(FileBuff *infile, FILE *outfile, FILE *noutfile, Qseqs *ref, Qseqs *seq, Qseqs *header, unsigned char *trans, unsigned norm, unsigned minLength, double minCov, unsigned flag, unsigned proxi, MethMotif *motif, FILE *diffile, int tnum) {
	
	int i, n, size, cSize, len, pair;
	unsigned **includes, **includesPtr;
	long unsigned **seqs, **seqsPtr;
	char **headers, **headersPtr;
	unsigned char *include;
	Matrix *D, *N;
	
	/* init */
	pair = flag & 2 ? 1 : 0;
	len = 0;
	cSize = 0;
	n = 0;
	size = 16;
	includes = pair ? smalloc(size * sizeof(unsigned *)) : smalloc(sizeof(unsigned *));
	seqs = smalloc(size * sizeof(long unsigned *));
	headers = smalloc(size * sizeof(char *));
	includesPtr = includes;
	seqsPtr = seqs;
	headersPtr = headers;
	
	/* load sequences from msa */
	while(FileBuffgetFsa(infile, header, seq, trans)) {
		/* realloc */
		if(n == size) {
			size <<= 1;
			includes = pair ? realloc(includes, size * sizeof(unsigned *)) : includes;
			seqs = realloc(seqs, size * sizeof(long unsigned *));
			headers = realloc(headers, size * sizeof(char *));
			if(!includes || !seqs || !headers) {
				ERROR();
			}
			for(i = n; i < size; ++i) {
				seqs[i] = smalloc(cSize * sizeof(long unsigned));
				if(pair) {
					includes[i] = smalloc(cSize * sizeof(unsigned));
				}
			}
			
			/* set pointers */
			includesPtr = pair ? includes + n : includes;
			seqsPtr = seqs + n;
			headersPtr = headers + n;
		}
		++n;
		
		/* add seq */
		if(ref->len) {
			if(seq->len != ref->len) {
				fprintf(stderr, "Sequences does not match: %s\n", header->seq);
				exit(1);
			}
			
			/* make / update inclusion array(s) */
			if(pair) {
				initIncPos(*includesPtr, len);
				qseq2nibble(seq, *seqsPtr);
				maskMotifs(*seqsPtr, *includesPtr, len, motif);
				getIncPosPtr(*includesPtr, seq, seq, proxi);
				*headersPtr = smalloc(header->len); 
				strcpy(*headersPtr, (char *)(header->seq + 1));
				if(getNpos(*includesPtr, len) < minLength) {
					fprintf(stderr, "Template (\"%s\") did not exceed threshold for inclusion\n", header->seq);
					--n;
				} else {
					++seqsPtr;
					++headersPtr;
					includesPtr += pair;
				}
			} else {
				qseq2nibble(seq, *seqsPtr);
				maskMotifs(*seqsPtr, *includes, len, motif);
				getIncPosPtr(*includes, seq, ref, proxi);
				*headersPtr = smalloc(header->len); 
				strcpy(*headersPtr, (char *)(header->seq + 1));
				++seqsPtr;
				++headersPtr;
				includesPtr += pair;
			}
		} else { /* use first valid sequence as validater */
			len = seq->len;
			minLength = minLength < minCov * len ? minCov * len : minLength;
			cSize = len / 32 + 1;
			for(i = 0; i < size; ++i) {
				seqs[i] = smalloc(cSize * sizeof(long unsigned));
				if(pair) {
					includes[i] = smalloc(cSize * sizeof(unsigned));
				}
			}
			if(!pair) {
				*includes = smalloc(cSize * sizeof(unsigned));
			}
			initIncPos(*includesPtr, len);
			qseq2nibble(seq, *seqsPtr);
			maskMotifs(*seqsPtr, *includesPtr, len, motif);
			getIncPosPtr(*includesPtr, seq, seq, proxi);
			
			if(getNpos(*includesPtr, len) < minLength) {
				fprintf(stderr, "Template (\"%s\") did not exceed threshold for inclusion\n", header->seq);
				ref->len = 0;
				--n;
			} else {
				*headersPtr = smalloc(header->len); 
				strcpy(*headersPtr, (char *)(header->seq + 1));
				/* swap seq and ref */
				exchange(seq->size, ref->size, len);
				exchange(seq->len, ref->len, len);
				exchange(seq->seq, ref->seq, include);
				++seqsPtr;
				++headersPtr;
				includesPtr += pair;
			}
		}
	}
	
	/* adjust size of arrays */
	includes = pair ? realloc(includes, n * sizeof(unsigned *)) : includes;
	seqs = realloc(seqs, n * sizeof(long unsigned *));
	headers = realloc(headers, n * sizeof(char *));
	if(!includes || !seqs || !headers) {
		ERROR();
	}
	include = smalloc(n);
	memset(include, 1, n);
	D = ltdMatrix_init(n);
	if(noutfile) {
		N = ltdMatrix_init(n);
	} else {
		N = 0;
	}
	
	/* adjust number of threads */
	if(n * (n - 1) / 2 < tnum) {
		fprintf(stderr, "Adjustning number of nodes to %d, to conform with the matrix size.\n", (tnum = n * (n - 1) / 2));
	}
	
	/* make ltd matrix */
	if(!n) {
		D->n = 0;
		errno = 1;
		fprintf(stderr, "All sequences were trimmed away.\n");
	} else if(pair) {
		fsaCmpThreadOut(tnum, &cmpairFsaThrd, D, N, n, len, seqs, include, includes, norm, minLength, minCov, diffile, 0, ref, 0, proxi);
		//cmpairFsa(D, N, numFile, len, seqs, include, includes, norm, minLength, minCov, proxi, diffile);
	} else if(minLength <= getNpos(*includes, len)) {
		fsaCmpThreadOut(tnum, &cmpFsaThrd, D, N, n, len, seqs, include, includes, norm, minLength, minCov, diffile, 0, ref, 0, proxi);
		/* i = cmpFsa(D, numFile, len, seqs, include, includes, norm, diffile);
		fprintf(stderr, "# %d / %d bases included in distance matrix.\n", i, len); */
	} else {
		D->n = 0;
		errno = 1;
		fprintf(stderr, "No sufficient overlap was found.\n");
	}
	
	/* print matrix */
	if(1 < D->n) {
		printphy(outfile, D, headers, include, 0, flag);
		if(N && 1 < N->n) {
			printphy(outfile, N, headers, include, 0, flag);
		}
	}
	
	/* clean up */
	for(i = 0; i < n; ++i) {
		free(seqs[i]);
		free(headers[i]);
		if(pair) {
			free(includes[i]);
		}
	}
	if(!pair && n) {
		free(*includes);
	}
	free(seqs);
	free(headers);
	free(include);
	free(includes);
	Matrix_destroy(D);
	Matrix_destroy(N);
	
	return D->n;
}

int ltdFsaRow_get(double *D, double *N, FileBuff *infile, char *targetTemplate, char *addfilename, char *diffilename, Qseqs **filenames, int n, unsigned norm, unsigned minLength, double minCov, unsigned flag, unsigned proxi) {
	/* doesn't work with floats, but its deprecated and unused anyway */
	char *filename;
	unsigned char *trans;
	unsigned *includeseq, *includeadd;
	unsigned i, j, len, inc;
	long unsigned *addL, *seqL, dist;
	double *Dptr, *Nptr;
	FILE *diffile;
	Qseqs *header, *seq;
	
	/* init  */
	len = 0;
	minLength = minLength < minCov * len ? minCov * len : minLength;
	trans = get2BitTable(flag);
	header = setQseqs(32);
	seq = setQseqs(1048576);
	
	/* load new sample into memory */
	/* open file, and find target */
	openAndDetermine(infile, addfilename);
	
	/* find the correct entry */
	while(FileBuffgetFsaHeader(infile, header) && strcmp((char *) header->seq, targetTemplate) != 0);
	
	/* get sequence */
	if(!header->len) {
		fprintf(stderr, "Missing template entry (\"%s\") in file:\t%s\n", targetTemplate, addfilename);
		exit(1);
	} else if(FileBuffgetFsaSeq(infile, seq, trans)) {
		closeFileBuff(infile);
	} else {
		fprintf(stderr, "\"%s\" is not fasta.\n", addfilename);
		exit(1);
	}
	if(len != 0 && len != seq->len) {
		fprintf(stderr, "New sequence does not match the existing sequences.\n");
		exit(1);
	} else {
		len = seq->len;
	}
	
	/* validate new seq */
	includeseq = smalloc(((seq->size >> 5) + 1) * sizeof(unsigned));
	includeadd = smalloc(((seq->size >> 5) + 1) * sizeof(unsigned));
	addL = smalloc(((seq->size >> 5) + 1) * sizeof(long unsigned));
	seqL = smalloc(((seq->size >> 5) + 1) * sizeof(long unsigned));
	initIncPos(includeadd, len);
	getIncPosPtr(includeadd, seq, seq, proxi);
	if(getNpos(includeadd, len) < minLength) {
		fprintf(stderr, "Template (\"%s\") did not exceed threshold for inclusion:\t%s\n", targetTemplate, addfilename);
		/* clean */
		destroyTable(trans);
		destroyQseqs(header);
		destroyQseqs(seq);
		free(includeseq);
		free(includeadd);
		free(addL);
		free(seqL);
		return 1;
	}  else if(diffilename) {
		if(*diffilename == '-' && diffilename[1] == 0) {
			diffile = stdout;
		} else {
			diffile = sfopen(diffilename, "ab");
		}
	} else {
		diffile = 0;
	}
	
	/* convert seq to nibbles */
	qseq2nibble(seq, addL);
	
	/* calculate distances */
	Dptr = D - 1;
	Nptr = N ? N - 1 : 0;
	--filenames;
	i = ++n;
	j = 0;
	while(--n) {
		filename = (char *) (*++filenames)->seq;
		/* open matrix file, and find target */
		openAndDetermine(infile, filename);
		while(FileBuffgetFsaHeader(infile, header) && strcmp((char *) header->seq, targetTemplate) != 0);
		FileBuffgetFsaSeq(infile, seq, trans);
		closeFileBuff(infile);
		
		/* get included positions */
		memcpy(includeseq, includeadd, (len / 32 + 1) * sizeof(unsigned));
		getIncPosPtr(includeseq, seq, seq, proxi);
		
		/* convert seq to nibbles */
		qseq2nibble(seq, seqL);
		
		/* get distance */
		if(diffile) {
			dist = fsacmpairint(diffile, i, ++j, addL, seqL, includeseq, len);
		} else {
			dist = fsacmpair(addL, seqL, includeseq, len);
		}
		
		/* separate distance and included bases */
		if(minLength <= (inc = dist & UINT_MAX)) {
			*++Dptr = (dist >> 32) * norm;
			*Dptr /= inc;
		} else {
			*++Dptr = -1.0;
			fprintf(stderr, "No sufficient overlap between samples:\t%s, %s\n", filename, addfilename);
		}
		if(N) {
			*++Nptr = inc;
		}
	}
	if(diffile && diffile != stdout) {
		fclose(diffile);
	}
	
	/* clean */
	destroyTable(trans);
	destroyQseqs(header);
	destroyQseqs(seq);
	free(includeseq);
	free(includeadd);
	free(addL);
	free(seqL);
	
	return 0;
}
