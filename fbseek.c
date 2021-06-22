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

#include <zlib.h>
#include <stdlib.h>
#include <string.h>
#include "fbseek.h"
#include "filebuff.h"
#include "pherror.h"

TimeStamp * timeStampFileBuff(FileBuff *src, TimeStamp *dest) {
	
	int status;
	z_stream *strm;
	
	if(src->buffFileBuff == &BuffgzFileBuff) {
		/* deflateCopy does not work -> opt out for gz */
		return 0;
	} else if(!dest) {
		dest = smalloc(sizeof(TimeStamp));
		dest->strm = 0;
		dest->size = 0;
		dest->buffer = 0;
	}
	
	dest->filePos = ftell(src->file);
	if(src->buffFileBuff == &BuffgzFileBuff) {
		/* gz file */
		if((strm = dest->strm)) {
			/* free state */
			deflateEnd(strm);
		} else {
			/* make new stream */
			dest->strm = (strm = smalloc(sizeof(z_stream)));
		}
		
		/* copy stream state */
		status = deflateCopy(strm, src->strm);
		strm->avail_in = 0;
		if(status != Z_OK) {
			fprintf(stderr, "Gzip error %d\n", status);
			errno |= status;
			exit(errno);
		}
		
		/* move filePos to next compressed chunk */
		dest->filePos -= strm->avail_in;
		
		/* copy remaining data */
		if(dest->size < src->bytes) {
			dest->buffer = smalloc((dest->size = src->bytes));
		}
		memcpy(dest->buffer, src->next, (dest->bytes = src->bytes));
	} else {
		/* regular file */
		dest->filePos -= src->bytes;
	}
	
	return dest;
}

int seekFileBiff(FileBuff *src, TimeStamp *dest) {
	
	sfseek(src->file, dest->filePos, SEEK_SET);
	if(src->buffFileBuff == &BuffgzFileBuff) {
		/* gz file */
		deflateEnd(src->strm);
		src->strm = dest->strm;
		
		/* get uncompressed leftover */
		memcpy(src->buffer, dest->buffer, (src->bytes = dest->bytes));
		src->next = src->buffer;
	} else {
		/* regular file */
		return src->buffFileBuff(src);
	}
	
	return src->bytes;
}

void freeTimeStamp(TimeStamp *dest) {
	
	if(dest->strm) {
		/* free state */
		deflateEnd(dest->strm);
		free(dest->strm);
		free(dest->buffer);
	}
	free(dest);
}
