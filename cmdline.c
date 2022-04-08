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

#include <stdlib.h>
#include "cmdline.h"

int getOptArg(const char *Arg) {
	
	char *arg;
	
	arg = (char *)(Arg);
	while(*arg && *arg != '=') {
		++arg;
	}
	
	return (int)(arg - Arg);
}

char * getArg(char ***Arg, int *argc, const int len) {
	
	char *arg;
	
	arg = (char *)(**Arg + len);
	if(*arg == 0) {
		if(argc) {
			arg = *(++*Arg);
			--*argc;
		} else {
			arg = 0;
		}
	}
	
	return arg;
}

char * getArgDie(char ***Arg, int *argc, const int len, const char *opt) {
	
	char *arg;
	
	if(!(arg = getArg(Arg, argc, len))) {
		missArg(opt);
	}
	
	return arg;
}

int getcArg(char ***Arg, int *argc, int len) {
	
	int c;
	char *arg;
	
	arg = getArg(Arg, argc, len);
	if(arg == 0) {
		c = 0;
	} else if(*arg == 0 || arg[1] != 0) {
		c = -1;
	} else {
		c = *arg;
	}
	
	return c;
}

int getcArgDie(char ***Arg, int *argc, int len, const char *opt) {
	
	int c;
	
	c = getcArg(Arg, argc, len);
	if(!c) {
		missArg(opt);
	} else if(c == -1) {
		invaArg(opt);
	}
	
	return c;
}

char * getDefArg(char ***Arg, int *argc, const int len, char *def) {
	
	char *arg;
	
	arg = (char *)(**Arg + len);
	if(*arg == 0) {
		if(argc) {
			arg = *(++*Arg);
			if(*arg == '-') {
				--*Arg;
				arg = 0;
			} else {
				--*argc;
			}
		} else {
			arg = 0;
		}
	}
	
	return arg;
}

long getNumArg(char ***Arg, int *argc, int len, const char *opt) {
	
	char *arg, *errorMsg;
	
	arg = getArg(Arg, argc, len);
	if(!arg) {
		missArg(opt);
	}
	len = strtol(arg, &errorMsg, 10);
	if(*errorMsg != 0) {
		invaArg(opt);
	}
	
	return len;
}

long getNumDefArg(char ***Arg, int *argc, int len, int def, const char *opt) {
	
	char *arg, *errorMsg;
	
	arg = getArg(Arg, argc, len);
	if(arg && *arg != '-') {
		len = strtol(arg, &errorMsg, 10);
		if(*errorMsg != 0) {
			invaArg(opt);
		}
	} else {
		--*Arg;
		++*argc;
		len = def;
	}
	
	return len;
}

double getdArg(char ***Arg, int *argc, int len, const char *opt) {
	
	char *arg, *errorMsg;
	double num;
	
	arg = getArg(Arg, argc, len);
	num = strtod(arg, &errorMsg);
	if(*errorMsg != 0) {
		invaArg(opt);
	}
	
	return num;
}

double getdDefArg(char ***Arg, int *argc, int len, double def, const char *opt) {
	
	char *arg, *errorMsg;
	double num;
	
	arg = getArg(Arg, argc, len);
	if(arg && *arg != '-') {
		num = strtod(arg, &errorMsg);
		if(*errorMsg != 0) {
			invaArg(opt);
		}
	} else {
		--*Arg;
		++*argc;
		num = def;
	}
	
	return num;
}

char ** getArgList(char ***Arg, int *argc, const int len) {
	
	char *arg, **argList;
	
	arg = (char *)(**Arg + len);
	if(*arg == 0) {
		if(argc) {
			argList = ++*Arg;
			--*argc;
		} else {
			argList = 0;
		}
	} else {
		argList = *Arg;
		*argList = arg;
	}
	
	return argList;
}

char ** getArgListDie(char ***Arg, int *argc, const int len, const char *opt) {
	
	char **argList;
	
	if(!(argList = getArgList(Arg, argc, len))) {
		missArg(opt);
	}
	
	return argList;
}

int getArgListLen(char ***Arg, int *argc) {
	
	int len;
	char *arg;
	
	len = 0;
	arg = **Arg;
	while(*argc && (*arg != '-' || arg[1] == 0)) {
		++len;
		--*argc;
		arg = *++*Arg;
	}
	--*Arg;
	++*argc;
	
	return len;
}
