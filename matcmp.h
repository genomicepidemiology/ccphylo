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

#include "filebuff.h"
#include "matparse.h"

void stripMat(MatrixCounts *mat1);
double nl1cmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2);
double nl2cmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2);
double nlncmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2);
double nlinfcmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2);
double l1cmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2);
double l2cmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2);
double lncmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2);
double linfcmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2);
double nbccmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2);
double bccmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2);
double nccmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2);
double ccmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2);
double zcmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2);
double pcmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2);
double npcmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2);
double chi2cmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2);
double nchi2cmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2);
double coscmp(short unsigned *counts1, short unsigned *counts2, int tot1, int tot2);
double cmpMats(MatrixCounts *mat1, NucCount *mat2, FileBuff *infile, unsigned norm, unsigned minDepth, unsigned minLength, double minCov, double (*veccmp)(short unsigned*, short unsigned*, int, int));
