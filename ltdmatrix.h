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

#include "fbseek.h"
#include "filebuff.h"
#include "matrix.h"
#include "matparse.h"
#include "matrix.h"

void ltdMatrix_get(Matrix *dest, Matrix *nDest, MatrixCounts *mat1, NucCount *mat2, FileBuff *infile, TimeStamp **targetStamps, unsigned char *include, char *targetTemplate, char **filenames, int numFile, unsigned norm, unsigned minDepth, unsigned minLength, double minCov, double (*veccmp)(short unsigned*, short unsigned*, int, int));
int ltdRow_get(double *D, double *N, char *targetTemplate, char *addfilename, Qseqs **filenames, int n, unsigned norm, unsigned minDepth, unsigned minLength, double minCov, double (*veccmp)(short unsigned*, short unsigned*, int, int));
