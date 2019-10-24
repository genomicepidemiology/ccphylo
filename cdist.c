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
#include "cdist.h"
#include "fbseek.h"
#include "matrix.h"
#include "pherror.h"
#include "phy.h"
#include "str.h"
#include "unionparse.h"

void ltd_get(Matrix *dest, Matrix *nDest, int numFile, int len, long unsigned **seqs, FileBuff *infile, TimeStamp **targetStamps, unsigned char *include, char *targetTemplate, char **filenames, unsigned norm, unsigned minLength, double minCov, unsigned flag, unsigned proxi) {
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
}

void makeMatrix_fsa(unsigned numFile, char **filenames, char *outputfilename, char *noutputfilename, char *targetTemplate, double minCov, unsigned norm, unsigned minLength, unsigned format) {
	
	unsigned n, pos, len;
	unsigned char *include;
	long unsigned **seqs;
	FILE *outfile, *noutfile;
	FileBuff *infile, *unionfile;
	Matrix *distMat, *nDest;
	TimeStamp **targetStamps;
	UnionEntry *entry;
	
	/* open output */
	if(*outputfilename == '-' && outputfilename[1] == '-' && outputfilename[2] == 0) {
		outfile = stdout;
	} else {
		outfile = sfopen(outputfilename, "wb");
	}
	if(noutputfilename) {
		if(strcmp(noutputfilename, outputfilename) == 0) {
			noutfile = outfile;
		} else if(*noutputfilename == '-' && noutputfilename[1] == '-' && noutputfilename[2] == 0) {
			noutfile = stdout;
		} else {
			noutfile = sfopen(noutputfilename, "wb");
		}
	} else {
		noutfile = 0;
	}
	
	/* init */
	infile = setFileBuff(1048576);
	
	if(targetTemplate) {
		/* init */
		distMat = ltdMatrix_init(numFile);
		nDest = ltdMatrix_init(numFile);
		if(!(targetStamps = calloc(numFile, sizeof(TimeStamp *)))) {
			ERROR();
		}
		include = smalloc(numFile);
		memset(include, 1, numFile);
		
		/* load seqs */
		seqs = smalloc(numFile * sizeof(long unsigned *));
		
		
		/* make ltd matrix */
		/* here */
		/* rewrite ltdMatrix_get
		ltdMatrix_get(distMat, nDest, mat1, mat2, infile, targetStamps, include, targetTemplate, filenames, numFile, norm, minDepth, minLength, minCov, veccmp);
		*/
		
		/* print ltd matrix */
		if(1 < distMat->n) {
			printphy(outfile, distMat, filenames, include, format);
			if(noutputfilename) {
				printphy(noutfile, nDest, filenames, include, format);
			}
		}
	} else if(numFile == 1) {
		/* create many matrices */
		/* requires *.union */
		unionfile = setFileBuff(524288);
		if(filenames) {
			openAndDetermine(unionfile, *filenames);
		} else {
			openAndDetermine(unionfile, "--");
		}
		/* get samplenames */
		if(!(filenames = UnionEntry_getHeader(unionfile, &numFile))) {
			fprintf(stderr, "Malformed union input.\n");
			exit(1);
		}
		/* init */
		entry = UnionEntry_init(32, numFile);
		distMat = ltdMatrix_init(numFile);
		nDest = ltdMatrix_init(numFile);
		if(!(targetStamps = calloc(numFile, sizeof(TimeStamp *)))) {
			ERROR();
		}
		include = smalloc(numFile);
		
		/* get filenames */
		n = numFile;
		while(n--) {
			/* truncate */
			len = strlen(filenames[n]);
			if((pos = rstrpos(filenames[n], '.', len)) != -1) {
				filenames[n][pos] = 0;
				len = pos;
			}
			
			/* add .mat.gz */
			strcpy(filenames[n] + len, ".fsa");
			len += 4;
		}
		
		/* parse candidates */
		while(UnionEntry_get(unionfile, entry)) {
			/* get included templates */
			memset(include, 0, numFile);
			n = entry->num;
			while(n--) {
				include[entry->filenames[n]] = 1;
			}
			
			/* get matrix */
			/* here */
			/* rewrite ltdMatrix_get
			ltdMatrix_get(distMat, nDest, mat1, mat2, infile, targetStamps, include, entry->target, filenames, numFile, norm, minDepth, minLength, minCov, veccmp);
			*/
			
			/* print ltd matrix */
			if(1 < distMat->n) {
				printphy(outfile, distMat, filenames, include, format);
				if(noutputfilename) {
					printphy(noutfile, nDest, filenames, include, format);
				}
			}
		}
		closeFileBuff(unionfile);
		destroyFileBuff(unionfile);
		UnionEntry_destroy(entry);
	} else {
		fprintf(stderr, "Missing input.\n");
		exit(1);
	}
	
	Matrix_destroy(distMat);
	Matrix_destroy(nDest);
	destroyFileBuff(infile);
	free(targetStamps);
	free(include);
	
	/* close */
	fclose(outfile);
	if(noutputfilename && strcmp(noutputfilename, outputfilename) != 0) {
		fclose(noutfile);
	}
}

int add2Matrix_fsa(char *path, char *addfilename, char *outputfilename, char *noutputfilename, char *targetTemplate, double minCov, unsigned norm, unsigned minLength, unsigned format) {
	
	int n, pos;
	double *D, *N;
	char *ptr;
	FILE *outfile, *noutfile;
	FileBuff *infile;
	Qseqs **names;
	
	/* init */
	infile = setFileBuff(1048576);
	openAndDetermine(infile, outputfilename);
	
	/* convert path to dir */
	pos = -1;
	n = 0;
	ptr = path - 1;
	while(*++ptr) {
		++n;
		if(*ptr == '/') {
			pos = n;
		}
	}
	if(0 <= pos) {
		path[pos] = 0;
	}
	
	/* get n */
	n = getSizePhy(infile);
	D = smalloc(n * sizeof(double));
	N = smalloc(n * sizeof(double));
	
	/* get names */
	if(!(names = getFilenamesPhy(path, n, infile))) {
		ERROR();
	} else if(infile->bytes) {
		fprintf(stderr, "Cannot update a multi distance phylip file.\n");
		return 1;
	}
	
	/* close phyfile */
	closeFileBuff(infile);
	
	/* open output and init new row(s) */
	outfile = sfopen(outputfilename, "rb+");
	if(noutputfilename) {
		noutfile = sfopen(noutputfilename, "rb+");
	} else {
		noutfile = 0;
	}
	
	/* calculate new row */
	/* here */
	/* rewrite ltdRow_get
	if(ltdRow_get(D, N, mat1, mat2, infile, targetTemplate, addfilename, names, n, norm, minDepth, minLength, minCov, veccmp)) {
		fprintf(stderr, "Distance measures failed and thus the matrix was not updated.\n");
		return 1;
	}
	*/
	
	/* print updated matrix */
	printphyUpdate(outfile, ++n, addfilename, D, format);
	fclose(outfile);
	if(noutfile) {
		printphyUpdate(noutfile, n, addfilename, N, format);
		fclose(noutfile);
	}
	
	/* clean */
	pos = n - 1;
	while(pos--) {
		destroyQseqs(names[pos]);
	}
	free(names);
	free(D);
	free(N);
	destroyFileBuff(infile);
	
	return 0;
}
