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
#include "meth.h"
#include "qseqs.h"

int ltdFsaMatrix_get(Matrix *D, Matrix *N, int numFile, long unsigned **seqs, int cSize, FileBuff *infile, TimeStamp **targetStamps, unsigned char *include, unsigned **includes, char *targetTemplate, char **filenames, unsigned char *trans, Qseqs *ref, Qseqs *seq, Qseqs *header, unsigned norm, unsigned minLength, double minCov, unsigned flag, unsigned proxi, MethMotif *motif, FILE *diffile, int tnum);
int ltdMsaMatrix_get(FileBuff *infile, FILE *outfile, FILE *noutfile, Qseqs *ref, Qseqs *seq, Qseqs *header, unsigned char *trans, unsigned norm, unsigned minLength, double minCov, unsigned flag, unsigned proxi, MethMotif *motif, FILE *diffile, int tnum);
int ltdFsaRow_get(double *D, double *N, FileBuff *infile, char *targetTemplate, char *addfilename, char *diffilename, Qseqs **filenames, int n, unsigned norm, unsigned minLength, double minCov, unsigned flag, unsigned proxi);
