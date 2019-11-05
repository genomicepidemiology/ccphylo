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
#include "hashmapstr.h"
#include "pherror.h"

int minimalStandard(int rand) {
	
	rand = 16807 * (rand % 127773) - 2836 * (rand / 127773);
	if (rand <= 0) {
			rand += 0x7fffffff;
	}
	
	return rand;
}

long unsigned djb2(unsigned char *str) {
	
	long unsigned hash;
	unsigned char c;
	
	hash = 5381;
	--str;
	while((c = *++str)) {
		hash = ((hash << 5) + hash) + c;
	}
	
	return minimalStandard(hash);
}

long unsigned next2power(long unsigned num) {
	
	long unsigned newNum;
	
	newNum = 1;
	while(newNum < num) {
		newNum <<= 1;
	}
	
	return newNum;
}

unsigned char * ustrdup(unsigned char *src) {
	
	unsigned char *dest;
	int len;
	
	len = strlen((char *) src);
	dest = smalloc(len + 1);
	
	if((*dest = *src) == 0) {
		return dest;
	}
	
	while((*++dest = *++src));
	
	return dest - len;
}

HashMapStr * HashMapStr_init(long unsigned size) {
	
	HashMapStr *dest;
	
	dest = smalloc(sizeof(HashMapStr));
	dest->n = 0;
	dest->size = next2power(size);
	if(!(dest->table = calloc(size, sizeof(BucketStr)))) {
		ERROR();
	}
	/* convert to masking problem */
	--dest->size;
	
	return dest;
}

void HashMapStr_grow(HashMapStr *src) {
	
	unsigned size, pos;
	BucketStr *node, *node_next, **table;
	
	/* reallocate hashtable */
	size = ++src->size;
	src->size <<= 1;
	src->table = realloc(src->table, src->size * sizeof(BucketStr));
	if(!src->table) {
		ERROR();
	}
	--src->size;
	
	/* move buckets top down */
	table = src->table + size;
	++size;
	while(--size) {
		node = *--table;
		*table = 0;
		while(node) {
			node_next = node->next;
			pos = node->hash & src->size;
			node->next = src->table[pos];
			src->table[pos] = node;
			node = node_next;
		}
	}
}

int HashMapStr_add(HashMapStr *src, unsigned char *str, unsigned n) {
	
	unsigned hash, pos;
	BucketStr *node;
	
	hash = djb2(str);
	pos = hash & src->size;
	
	/* search hashmap */
	for(node = src->table[pos]; node; node = node->next) {
		if(hash == node->hash && strcmp((char *) str, (char *) node->str) == 0) {
			if(!(node->uList = realloc(node->uList, (++node->num + 1) * sizeof(unsigned)))) {
				ERROR();
			}
			node->uList[node->num] = n;
			return node->num;
		}
	}
	
	if(++src->n == src->size) {
		HashMapStr_grow(src);
	}
	node = smalloc(sizeof(BucketStr));
	node->str = ustrdup(str);
	node->hash = hash;
	node->num = 0;
	node->uList = smalloc(8 * sizeof(unsigned));
	*(node->uList) = n;
	node->next = src->table[pos];
	src->table[pos] = node;
	
	return 0;
}

BucketStr * HashMapStr_get(HashMapStr *src, unsigned char *str) {
	
	unsigned hash, pos;
	BucketStr *node, *prev;
	
	hash = djb2(str);
	pos = hash & src->size;
	
	/* search hashmap */
	prev = 0;
	node = src->table[pos];
	while(node) {
		if(hash == node->hash && strcmp((char *) str, (char *) node->str) == 0) {
			if(prev) {
				prev->next = node->next;
			} else {
				src->table[pos] = node->next;
			}
			--src->n;
			return node;
		}
		prev = node;
		node = node->next;
	}
	
	return 0;
}

int HashMapStr_print(HashMapStr *src, FILE *out) {
	
	unsigned size, nc, num, *ptr;
	BucketStr *node, **table;
	
	nc = 0;
	size = src->size + 2;
	table = src->table - 1;
	while(--size) {
		for(node = *++table; node; node = node->next) {
			if((num = node->num)) {
				nc += fprintf(out, "%s\t%d", node->str, ++num);
				++num;
				ptr = node->uList - 1;
				while(--num) {
					nc += fprintf(out, "\t%u", *++ptr);
				}
				nc += fprintf(out, "\n");
			}
		}
	}
	
	return nc;
}

void HashMapStr_destroy(HashMapStr *src) {
	
	unsigned size;
	BucketStr *node, *node_next, **table;
	
	size = ++src->size + 1;
	table = src->table - 1;
	while(--size && src->n) {
		for(node = *++table; node; node = node_next) {
			node_next = node->next;
			free(node->str);
			free(node);
			--src->n;
		}
	}
	free(src->table);
	free(src);
}
