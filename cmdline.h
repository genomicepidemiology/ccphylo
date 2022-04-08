/* Philip T.L.C. Clausen Mar 2022 plan@dtu.dk */

/*
 * Copyright (c) 2022, Philip Clausen, Technical University of Denmark
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

#define missArg(opt) fprintf(stderr, "Missing argument at %s.\n", opt); exit(1);
#define invaArg(opt) fprintf(stderr, "Invalid value parsed at %s.\n", opt); exit(1);
#define unknArg(arg) fprintf(stderr, "Unknown argument or option: \"%s\"\n", arg); exit(1);
#define nonOptError() fprintf(stderr, "Too many non-option arguments handed.\n"); exit(1);

int getOptArg(const char *Arg);
char * getArg(char ***Arg, int *argc, const int len);
char * getArgDie(char ***Arg, int *argc, const int len, const char *opt);
int getcArg(char ***Arg, int *argc, int len);
int getcArgDie(char ***Arg, int *argc, int len, const char *opt);
char * getDefArg(char ***Arg, int *argc, const int len, char *def);
long getNumArg(char ***Arg, int *argc, int len, const char *opt);
long getNumDefArg(char ***Arg, int *argc, int len, int def, const char *opt);
double getdArg(char ***Arg, int *argc, int len, const char *opt);
double getdDefArg(char ***Arg, int *argc, int len, double def, const char *opt);
char ** getArgList(char ***Arg, int *argc, const int len);
char ** getArgListDie(char ***Arg, int *argc, const int len, const char *opt);
int getArgListLen(char ***Arg, int *argc);
