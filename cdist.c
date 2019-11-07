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
#include "pherror.h"
#include "phy.h"
#include "seqparse.h"

int ltdFsaMatrix_get(Matrix *D, Matrix *N, int numFile, long unsigned **seqs, int cSize, FileBuff *infile, TimeStamp **targetStamps, unsigned char *include, unsigned **includes, char *targetTemplate, char **filenames, unsigned char *trans, Qseqs *ref, Qseqs *seq, Qseqs *header, unsigned norm, unsigned minLength, double minCov, unsigned flag, unsigned proxi, FILE *diffile, int tnum) {
	
	unsigned i, pair, len, **includesPtr;
	long unsigned **seqsPtr;
	unsigned char *includePtr;
	TimeStamp **targetStamp;
	
	/* Validate reference file */
	if(openAndDetermine(infile, *filenames) == '>') {
		/* find the correct entry */
		if(*(targetStamp = targetStamps)) {
			seekFileBiff(infile, *targetStamp);
		}
		while(FileBuffgetFsaHeader(infile, header) && strcmp((char *) header->seq, targetTemplate) != 0);
		
		/* make timestamp */
		*targetStamp = timeStampFileBuff(infile, *targetStamp);
		
		/* get sequence */
		if(!header->len) {
			fprintf(stderr, "Missing template entry (\"%s\") in file:\t%s\n", targetTemplate, *filenames);
			exit(1);
		} else if(FileBuffgetFsaSeq(infile, ref, trans)) {
			/* make initial inclusion array */
			len = ref->len;
			minLength = minLength < minCov * len ? minCov * len : minLength;
			pair = flag & 2;
			if(cSize < len / 32 + 1) {
				cSize = len / 32 + 1;
				i = numFile;
				while(--i) {
					if(!(seqs[i] = realloc(seqs[i], cSize * sizeof(long unsigned)))) {
						ERROR();
					}
					if(pair && (includes[i] = realloc(includes[i], cSize * sizeof(unsigned))) == 0) {
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
			initIncPos(*includes, len);
			getIncPos(*includes, ref, ref, proxi);
			
			if(getNpos(*includes, len) < minLength) {
				fprintf(stderr, "Template (\"%s\") did not exceed threshold for inclusion:\t%s\n", targetTemplate, *filenames);
				exit(1);
			}
			if((*include = (flag & 16) ? 1 : 0)) {
				qseq2nibble(ref, *seqs);
			}
		} else {
			fprintf(stderr, "Missing template sequence (\"%s\") in file:\t%s\n", targetTemplate, *filenames);
			exit(1);
		}
		
		closeFileBuff(infile);
	} else {
		fprintf(stderr, "\"%s\" is not fasta.\n", *filenames);
		exit(1);
	}
	
	/* init */
	includesPtr = includes;
	includePtr = include;
	seqsPtr = seqs;
	targetStamp = targetStamps;
	
	/* load sequences, get included positions and compress */
	i = numFile;
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
				++includesPtr;
				++seqsPtr;
			} else if(FileBuffgetFsaSeq(infile, seq, trans)) {
				if(seq->len != ref->len) {
					fprintf(stderr, "Sequences does not match: %s\n", *filenames);
					exit(1);
				}
				
				/* make / update inclusion array(s) */
				if(pair) {
					memcpy(*++includesPtr, *includes, ((len >> 5) + 1) * sizeof(unsigned));
					getIncPos(*includesPtr, seq, ref, proxi);
					if(getNpos(*includesPtr, len) < minLength) {
						fprintf(stderr, "Template (\"%s\") did not exceed threshold for inclusion:\t%s\n", targetTemplate, *filenames);
						*includePtr = 0;
					} else {
						*includePtr = 1;
					}
				} else {
					getIncPos(*includes, seq, ref, proxi);
					*includePtr = 1;
				}
				if(*includePtr) {
					/* convert seq to nibbles */
					qseq2nibble(seq, *++seqsPtr);
				} else {
					++seqsPtr;
				}
			} else {
				fprintf(stderr, "Missing template sequence (\"%s\") in file:\t%s\n", targetTemplate, *filenames);
				*includePtr = 0;
				++includesPtr;
				++seqsPtr;
			}
			
			closeFileBuff(infile);
		} else {
			++filenames;
			++targetStamp;
			++includesPtr;
			++seqsPtr;
		}
	}
	
	/* make ltd matrix */
	if(pair) {
		fsaCmpThreadOut(tnum, &cmpairFsaThrd, D, N, numFile, len, seqs, include, includes, norm, minLength, minCov, diffile, targetTemplate, ref, 0, proxi);
		//cmpairFsa(D, N, numFile, len, seqs, include, includes, norm, minLength, minCov, diffile);
	} else if(minLength <= getNpos(*includes, len)) {
		fsaCmpThreadOut(tnum, &cmpFsaThrd, D, N, numFile, len, seqs, include, includes, norm, minLength, minCov, diffile, targetTemplate, ref, 0, proxi);
		/* i = cmpFsa(D, numFile, len, seqs, include, includes, norm, diffile);
		fprintf(stderr, "# %d / %d bases included in distance matrix.\n", i, len); */
	} else {
		D->n = 0;
		fprintf(stderr, "No sufficient overlap was found.\n");
	}
	
	return cSize;
}

int ltdFsaRow_get(double *D, double *N, FileBuff *infile, char *targetTemplate, char *addfilename, char *diffilename, Qseqs **filenames, int n, unsigned norm, unsigned minLength, double minCov, unsigned flag, unsigned proxi) {
	
	char *filename;
	unsigned char *trans;
	unsigned *includeseq, *includeref, *includeadd, i, j, len, inc;
	long unsigned *addL, *seqL, dist;
	double *Dptr, *Nptr;
	FILE *diffile;
	Qseqs *header, *seq, *ref;
	
	/* init  */
	len = 0;
	minLength = minLength < minCov * len ? minCov * len : minLength;
	trans = get2BitTable(flag);
	header = setQseqs(32);
	ref = setQseqs(1048576);
	
	/* load reference */
	if(flag & 16) {
		openAndDetermine(infile, (char *) (*filenames)->seq);
		while(FileBuffgetFsaHeader(infile, header) && strcmp((char *) header->seq, targetTemplate) != 0);
		
		/* get sequence */
		if(!header->len) {
			fprintf(stderr, "Missing template entry (\"%s\") in file:\t%s\n", targetTemplate, (*filenames)->seq);
			exit(1);
		} else if(FileBuffgetFsaSeq(infile, ref, trans)) {
			closeFileBuff(infile);
		} else {
			fprintf(stderr, "\"%s\" is not fasta.\n", (*filenames)->seq);
			exit(1);
		}
		includeref = smalloc(((ref->size >> 5) + 1) * sizeof(unsigned));
		len = ref->len;
		initIncPos(includeref, len);
		getIncPos(includeref, ref, ref, proxi);
		closeFileBuff(infile);
	} else {
		includeref = smalloc(((ref->size >> 5) + 1) * sizeof(unsigned));
		initIncPos(includeref, ref->size);
	}
	seq = setQseqs(ref->size);
	includeseq = smalloc(((seq->size >> 5) + 1) * sizeof(unsigned));
	includeadd = smalloc(((seq->size >> 5) + 1) * sizeof(unsigned));
	addL = smalloc(((seq->size >> 5) + 1) * sizeof(long unsigned));
	seqL = smalloc(((seq->size >> 5) + 1) * sizeof(long unsigned));
	
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
	memcpy(includeadd, includeref, ((seq->len >> 5) + 1) * sizeof(unsigned));
	getIncPos(includeadd, seq, ref, proxi);
	if(getNpos(includeadd, len) < minLength) {
		fprintf(stderr, "Template (\"%s\") did not exceed threshold for inclusion:\t%s\n", targetTemplate, addfilename);
		return 1;
	} else if(diffilename) {
		diffile = fopen(diffilename, "ab");
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
		memcpy(includeseq, includeref, ((len >> 5) + 1) * sizeof(unsigned));
		getIncPos(includeseq, seq, ref, proxi);
		
		/* convert seq to nibbles */
		qseq2nibble(seq, seqL);
		
		/* get distance */
		if(diffile) {
			dist = fsacmpairint(diffile, i, ++j, addL, seqL, includeadd, includeseq, len);
		} else {
			dist = fsacmpair(addL, seqL, includeadd, includeseq, len);
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
	fclose(diffile);
	
	/* clean */
	free(trans);
	destroyQseqs(header);
	destroyQseqs(ref);
	destroyQseqs(seq);
	free(includeref);
	free(includeseq);
	free(includeadd);
	free(addL);
	free(seqL);
	
	return 0;
}
