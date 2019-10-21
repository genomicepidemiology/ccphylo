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

#ifndef HASHMAPSTR
typedef struct hashMapStr HashMapStr;
typedef struct bucketStr BucketStr;
struct bucketStr {
	unsigned char *str;
	unsigned hash;
	unsigned num;
	unsigned *uList;
	struct bucketStr *next;
};
struct hashMapStr {
	unsigned n;
	unsigned size;
	struct bucketStr **table;
};
#define HASHMAPSTR 1
#endif

int minimalStandard(int rand);
long unsigned djb2(unsigned char *str);
long unsigned next2power(long unsigned num);
unsigned char * ustrdup(unsigned char *src);
HashMapStr * HashMapStr_init(long unsigned size);
void HashMapStr_grow(HashMapStr *src);
int HashMapStr_add(HashMapStr *src, unsigned char *str, unsigned n);
BucketStr * HashMapStr_get(HashMapStr *src, unsigned char *str);
int HashMapStr_print(HashMapStr *src, FILE *out);
void HashMapStr_destroy(HashMapStr *src);
