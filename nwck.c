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

#include "nwck.h"
#include "qseqs.h"
#include "pherror.h"
#define exchange(src1, src2, tmp) tmp = src1; src1 = src2; src2 = tmp;

void formNode(Qseqs *node1, Qseqs *node2, double L1, double L2) {
	
	unsigned char *seq;
	unsigned newsize;
	double L;
	
	/* move largest qseq down */
	if(node1->size < node2->size) {
		exchange(node1->seq, node2->seq, seq);
		exchange(node1->size, node2->size, newsize);
		exchange(node1->len, node2->len, newsize);
		exchange(L1, L2, L);
	}
	
	/* alloc */
	newsize = node1->len + node2->len + 32;
	if(node1->size < newsize) {
		if(!(node1->seq = realloc(node1->seq, newsize))) {
			ERROR();
		}
		node1->size = newsize;
	}
	
	/* form node */
	node1->seq[node1->len] = 0;
	if(L1 < 0 && L2 < 0) {
		node1->len = sprintf((char *) node1->seq, "%c%s,%s)", *node1->seq, (char *) node1->seq, (char *) node2->seq);
	} else {
		node1->len = sprintf((char *) node1->seq, "%c%s:%.2f,%s:%.2f)", *node1->seq, (char *) node1->seq, L1, (char *) node2->seq, L2);
	}
	*node1->seq = '(';
}

void formLastNode(Qseqs *node1, Qseqs *node2, double L) {
	
	unsigned char *seq;
	unsigned newsize;
	
	/* move largest qseq down */
	if(node1->size < node2->size) {
		exchange(node1->seq, node2->seq, seq);
		exchange(node1->size, node2->size, newsize);
		exchange(node1->len, node2->len, newsize);
	}
	
	/* alloc */
	newsize = node1->len + node2->len + 32;
	if(node1->size < newsize) {
		if(!(node1->seq = realloc(node1->seq, newsize))) {
			ERROR();
		}
		node1->size = newsize;
	}
	
	/* truncate node */
	node1->seq[--node1->len] = 0;
	
	/* form node */
	if(L < 0) {
		node1->len += sprintf((char *) node1->seq + node1->len, ",%s)", (char *) node2->seq);
	} else {
		node1->len += sprintf((char *) node1->seq + node1->len, ",%s:%.2f)", (char *) node2->seq, L);
	}
}
