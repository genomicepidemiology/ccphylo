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

#ifndef METH
typedef struct methMotif MethMotif;
struct methMotif {
	int num;
	int len;
	long unsigned *motif;
	unsigned *mask;
	struct methMotif *next;
};
#define METH 1
#endif

MethMotif * newMethMotif(int num, int len);
void destroyMethMotifs(MethMotif *src);
long unsigned matchMotif32(long unsigned seq, long unsigned *motif, unsigned num);
int matchMotif(long unsigned *seq, int pos, int seqlen, MethMotif *motif);
void maskMotif(unsigned *include, int pos, unsigned *mask, int len);
int maskMotifs(long unsigned *seq, unsigned *include, int seqlen, MethMotif *motif);
