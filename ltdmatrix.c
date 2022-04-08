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
#include "bytescale.h"
#include "fbseek.h"
#include "filebuff.h"
#include "matcmp.h"
#include "matparse.h"
#include "matrix.h"
#include "pherror.h"
#define strcmp2(str1, str2) (*(str1) == *(str2) && strcmp(str1 + 1, str2 + 1) == 0)

void ltdMatrix_get(Matrix *dest, Matrix *nDest, MatrixCounts *mat1, NucCount *mat2, FileBuff *infile, TimeStamp **targetStamps, unsigned char *include, char *targetTemplate, char **filenames, int numFile, unsigned norm, unsigned minDepth, unsigned minLength, double minCov, double (*veccmp)(short unsigned*, short unsigned*, int, int)) {
	
	int i, j, n, srtd;
	char **filename;
	double dist, *mat, *nMat;
	float *fmat, *nfMat;
	short unsigned *smat, *nsMat;
	unsigned char *bmat, *nbMat;
	TimeStamp **targetStamp;
	
	/* get distances */
	srtd = 1;
	dest->n = 0;
	mat = 0;
	nMat = 0;
	fmat = 0;
	nfMat = 0;
	smat = 0;
	nsMat = 0;
	bmat = 0;
	nbMat = 0;
	if(dest->mat) {
		mat = *(dest->mat);
		nMat = *(nDest->mat);
	} else if(dest->fmat) {
		fmat = *(dest->fmat);
		nfMat = *(nDest->fmat);
	} else if(dest->smat) {
		smat = *(dest->smat);
		nsMat = *(nDest->smat);
	} else {
		bmat = *(dest->bmat);
		nbMat = *(nDest->bmat);
	}
	
	for(i = 1; i < numFile; ++i) {
		if(include[i]) {
			/* open matrix file, and find target */
			openAndDetermine(infile, filenames[i]);
			if(*(targetStamp = targetStamps + i) && srtd) {
				seekFileBiff(infile, *targetStamp);
				/* first time -> seek to new template */
				while(FileBuffSkipTemplate(infile, mat2) && !(strcmp2(targetTemplate, (char *) mat2->name)));
			} else {
				while(FileBuffSkipTemplate(infile, mat2) && !(strcmp2(targetTemplate, (char *) mat2->name)));
			}
			
			/* make timestamp */
			*targetStamp = timeStampFileBuff(infile, *targetStamp);
			include[i] = 2;
			
			/* initialize mat1 */
			setMatName(mat1, mat2);
			if(FileBuffLoadMat(mat1, infile, minDepth) == 0) {
				fprintf(stderr, "Input is not DB sorted.\n");
				/* reset strm */
				sfseek(infile->file, 0, SEEK_SET);
				while(FileBuffSkipTemplate(infile, mat2) && !(strcmp2(targetTemplate, (char *) mat2->name)));
				
				if(strcmp2(targetTemplate, (char *) mat2->name)) {
					/* try load again */
					setMatName(mat1, mat2);
					if(FileBuffLoadMat(mat1, infile, minDepth) == 0) {
						fprintf(stderr, "Malformed matrix in:\t%s\n", filenames[i]);
						exit(1);
					}
				} else {
					fprintf(stderr, "Template (\"%s\") was not found in sample:\t%s\n", targetTemplate, filenames[i]);
					mat1->nNucs = 0;
					include[i] = 0;
				}
				/* mark run as unsorted */
				srtd = 0;
			}
			closeFileBuff(infile);
			
			/* validate matrix */
			if(mat1->nNucs < minLength || mat1->nNucs < minCov * mat1->len) {
				fprintf(stderr, "Template (\"%s\") did not exceed threshold for inclusion:\t%s\n", targetTemplate, filenames[i]);
				include[i] = 0;
			} else {
				/* strip matrix for insersions */
				stripMat(mat1);
			}
		}
		
		if(include[i]) {
			targetStamp = targetStamps;
			/* calculate distances */
			filename = filenames;
			n = 0;
			j = i;
			while(j--) {
				if(include[n]) {
					/* open matrix file, and find target */
					openAndDetermine(infile, *filename);
					if(*targetStamp && srtd) {
						seekFileBiff(infile, *targetStamp);
						if(include[n] == 1) {
							while(FileBuffSkipTemplate(infile, mat2) && !(strcmp2(targetTemplate, (char *) mat2->name)));
							if(!(strcmp2(targetTemplate, (char *) mat2->name))) {
								srtd = 0;
								sfseek(infile->file, 0, SEEK_SET);
								while(FileBuffSkipTemplate(infile, mat2) && !(strcmp2(targetTemplate, (char *) mat2->name)));
							}
						} else {
							/* set name */
							memcpy(mat2->name, mat1->name->seq, mat1->name->len);
						}
					} else {
						while(FileBuffSkipTemplate(infile, mat2) && !(strcmp2(targetTemplate, (char *) mat2->name)));
					}
					
					/* make timestamp */
					if(include[n] == 1) {
						*targetStamp = timeStampFileBuff(infile, *targetStamp);
						include[n] = 2;
					}
					
					/* get distance between the matrices */
					dist = cmpMats(mat1, mat2, infile, norm, minDepth, minLength, minCov, veccmp);
					if(dist < 0) {
						if(dist == -1.0) {
							fprintf(stderr, "No sufficient overlap between samples:\t%s, %s\n", filenames[i], *filename);
						} else if(dist == -2.0) {
							fprintf(stderr, "Template (\"%s\") did not exceed threshold for inclusion:\t%s\n", targetTemplate, *filename);
						} else {
							fprintf(stderr, "Failed to produce a distance metric between samples:\t%s, %s\n", filenames[i], *filename);
						}
					}
					
					if(-1.0 <= dist) {
						if(mat) {
							*mat++ = dist;
							*nMat++ = mat2->total;
						} else if(fmat) {
							*fmat++ = dist;
							*nfMat++ = mat2->total;
						} else if(smat) {
							*smat++ = dtouc(dist, 0.5);
							*nsMat++ = dtouc(mat2->total, 0.5);
						} else {
							*bmat++ = dtouc(dist, 0.5);
							*nbMat++ = dtouc(mat2->total, 0.5);
						}
					} else {
						include[n] = 0;
					}
					
					/* close mtrix file */
					closeFileBuff(infile);
				}
				++targetStamp;
				++filename;
				++n;
			}
		}
	}
	
	/* get number of include templates */
	++numFile;
	--include;
	n = 0;
	while(--numFile) {
		if(*++include) {
			++n;
		}
	}
	dest->n = n;
	nDest->n = n;
	
}

int ltdRow_get(double *D, double *N, char *targetTemplate, char *addfilename, Qseqs **filenames, int n, unsigned norm, unsigned minDepth, unsigned minLength, double minCov, double (*veccmp)(short unsigned*, short unsigned*, int, int)) {
	
	char *filename;
	double dist, *mat, *nMat;
	FileBuff *infile;
	MatrixCounts *mat1;
	NucCount *mat2;
	
	/* init */
	infile = setFileBuff(1048576);
	mat1 = initMat(1048576, 128);
	mat2 = initNucCount(128);
	
	/* load new sample matrix into memory */
	/* open matrix file, and find target */
	openAndDetermine(infile, addfilename);
	while(FileBuffSkipTemplate(infile, mat2) && !(strcmp2(targetTemplate, (char *) mat2->name)));
	
	/* initialize mat1 */
	setMatName(mat1, mat2);
	if(FileBuffLoadMat(mat1, infile, minDepth) == 0) {
		fprintf(stderr, "Malformed matrix in:\t%s\n", addfilename);
		exit(1);
	}
	closeFileBuff(infile);
	
	/* validate matrix */
	if(mat1->nNucs < minLength || mat1->nNucs < minCov * mat1->len) {
		fprintf(stderr, "Template (\"%s\") did not exceed threshold for inclusion:\t%s\n", targetTemplate, addfilename);
		return 1;
	}
	
	/* strip matrix for insersions */
	stripMat(mat1);
	
	/* calculate distances */
	mat = D;
	nMat = N;
	--filenames;
	++n;
	while(--n) {
		filename = (char *) (*++filenames)->seq;
		/* open matrix file, and find target */
		openAndDetermine(infile, filename);
		while(FileBuffSkipTemplate(infile, mat2) && !(strcmp2(targetTemplate, (char *) mat2->name)));
		
		/* get distance between the matrices */
		dist = cmpMats(mat1, mat2, infile, norm, minDepth, minLength, minCov, veccmp);
		if(dist < 0) {
			if(dist == -1.0) {
				fprintf(stderr, "No sufficient overlap between samples:\t%s, %s\n", filename, addfilename);
			} else if(dist == -2.0) {
				fprintf(stderr, "Template (\"%s\") did not exceed threshold for inclusion:\t%s\n", targetTemplate, filename);
				return 1;
			} else {
				fprintf(stderr, "Failed to produce a distance metric between samples:\t%s, %s\n", filename, addfilename);
				return 1;
			}
		}
		
		*mat++ = dist;
		if(N) {
			*nMat++ = mat2->total;
		}
		
		/* close mtrix file */
		closeFileBuff(infile);
	}
	
	/* clean */
	destroyFileBuff(infile);
	destroyMat(mat1);
	destroyNucCount(mat2);
	
	return 0;
}
