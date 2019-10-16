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

#include <stdlib.h>
#include <string.h>
#include "filebuff.h"
#include "pherror.h"
#include "qseqs.h"
#include "resparse.h"

ResEntry * ResEntry_init(unsigned size) {
	
	ResEntry *dest;
	
	dest = smalloc(sizeof(ResEntry));
	dest->Template = setQseqs(128);
	
	return dest;
}

int FileBuffValidateHeader(FileBuff *src) {
	
	char header[] = "#Template\tScore\tExpected\tTemplate_length\tTemplate_Identity\tTemplate_Coverage\tQuery_Identity\tQuery_Coverage\tDepth\tq_value\tp_value\n";
	
	if(src->bytes < 129 || strncmp((char *) src->next, header, 129) != 0) {
		return 1;
	}
	src->bytes -= 129;
	src->next += 129;
	
	return 0;
}

int FileBuffGetEntry(FileBuff *src, ResEntry *dest) {
	
	unsigned char *buff, *seq, c, stop;
	char *msg, strbuff[32];
	int i, size, avail, (*buffFileBuff)(FileBuff *);
	long unsigned num, numbuff[3];
	double doublebuff[7];
	Qseqs *Template;
	
	/* init */
	avail = src->bytes;
	buff = src->next;
	buffFileBuff = src->buffFileBuff;
	if(avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	
	/* get template */
	Template = dest->Template;
	seq = Template->seq;
	size = Template->size;
	while((*seq++ = *buff++) != '\t') {
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				fprintf(stderr, "Malformatted res file\n");
				errno |= 1;
				return 0;
			}
			buff = src->buffer;
		}
		if(--size == 0) {
			size = Template->size;
			Template->size <<= 1;
			Template->seq = realloc(Template->seq, Template->size);
			if(!Template->seq) {
				ERROR();
			}
			seq = Template->seq + size;
		}
	}
	if(--avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			fprintf(stderr, "Malformatted res file\n");
			errno |= 1;
			return 0;
		}
		buff = src->buffer;
	}
	Template->len = Template->size - size;
	*--seq = 0;
	
	/* get whole numbers */
	i = 3;
	while(i--) {
		num = 0;
		while((c = *buff++) != '\t') {
			if(--avail == 0) {
				if((avail = buffFileBuff(src)) == 0) {
					return 0;
				}
				buff = src->buffer;
			}
			if(c != ' ') {
				num = num * 10 + (c - '0');
			}
		}
		numbuff[i] = num;
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				return 0;
			}
			buff = src->buffer;
		}
	}
	dest->Template_length = numbuff[0];
	dest->Expected = numbuff[1];
	dest->Score = numbuff[2];
	
	/* get doubles */
	i = 7;
	while(i--) {
		stop = i != 0 ? '\t' : '\n';
		seq = (unsigned char *) strbuff;
		while((c = *buff++) != stop) {
			if(--avail == 0) {
				if((avail = buffFileBuff(src)) == 0) {
					return 0;
				}
				buff = src->buffer;
			}
			if(c != ' ') {
				*seq++ = c;
			}
		}
		*seq = 0;
		
		doublebuff[i] = strtod(strbuff, &msg);
		if(*msg != 0) {
			fprintf(stderr, "Malformatted res file\n");
			errno |= 1;
			return 0;
		} else if(--avail == 0 && stop != '\n') {
			if((avail = buffFileBuff(src)) == 0) {
				return 0;
			}
			buff = src->buffer;
		}
	}
	dest->p_value = doublebuff[0];
	dest->q_value = doublebuff[1];
	dest->Depth = doublebuff[2];
	dest->Query_Coverage = doublebuff[3];
	dest->Query_Identity = doublebuff[4];
	dest->Template_Coverage = doublebuff[5];
	dest->Template_Identity = doublebuff[6];
	
	src->bytes = avail;
	src->next = buff;
	return 1;
}

int FileBuffGetTemplate(FileBuff *src, ResEntry *dest) {
	
	unsigned char *buff, *seq;
	int size, avail;
	Qseqs *Template;
	int (*buffFileBuff)(FileBuff *);
	
	/* init */
	avail = src->bytes;
	buff = src->next;
	buffFileBuff = src->buffFileBuff;
	if(avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			return 0;
		}
		buff = src->buffer;
	}
	
	/* get template */
	Template = dest->Template;
	seq = Template->seq;
	size = Template->size;
	while((*seq++ = *buff++) != '\t') {
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				return 0;
			}
			buff = src->buffer;
		}
		if(--size == 0) {
			size = Template->size;
			Template->size <<= 1;
			Template->seq = realloc(Template->seq, Template->size);
			if(!Template->seq) {
				ERROR();
			}
			seq = Template->seq + size;
		}
	}
	if(--avail == 0) {
		if((avail = buffFileBuff(src)) == 0) {
			fprintf(stderr, "Malformatted res file\n");
			errno |= 1;
			return 0;
		}
		buff = src->buffer;
	}
	Template->len = Template->size - size;
	*--seq = 0;
	
	/* skip the rest */
	while(*buff++ != '\n') {
		if(--avail == 0) {
			if((avail = buffFileBuff(src)) == 0) {
				fprintf(stderr, "Malformatted res file\n");
				errno |= 1;
				return 0;
			}
			buff = src->buffer;
		}
	}
	
	src->bytes = avail - 1;
	src->next = buff;
	return 1;
}
