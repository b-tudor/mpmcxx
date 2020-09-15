#pragma once
#ifndef SAFEOPS_H
#define SAFEOPS_H

#include <stdio.h>
#include <stdint.h>
#include <cstdlib>

#include "constants.h"
#include "Output.h"

class SafeOps {
	SafeOps() {}
	~SafeOps() {}

public:

	static void   closeFile(FILE *fp );
	static FILE *  openFile( const char *inFilename, const char *mode, int line, const char *file ); 
	// Attempts to open files and throws exception upon failure

	

	 // string-to-data converters
	///////////////////////////////////////////////////////////////////////////////

	//convert string *s to a double and store at *d.
	static bool atod( const char * s, double &d );
	
	//convert string *s to an int and store at *i.
	static bool atoi( const char * s, int &i );
	
	//convert string *s to an unsigned int and stores at *u.
	static bool atou( const char * s, uint32_t &u );
	
	//converts string *s to a long and stores at *l.
	static bool atol( const char * s, long &l );
	
	// compares first n characters of A and B (ignoring case) and returns true if they match
	static bool strncasecmp( const char *A, const char *B, size_t n );
	
	// returns true if all elements of string are identical, compared case-insensitively
	static bool iequals( const char *A, const char *B );
	
	// returns true if a file exists
	static bool file_exists(const char* filename) {

		if( FILE *file = fopen(filename, "r")) {
			fclose(file);
			return true;
		}

		return false;
	}

	// The following memory functions were written to hopefully deal with some memory  
	// issues that have been plaguing very large runs on PCN61. I believe that we may  
	// be running into memory fragmentation. I am adding a function that will add an
	// option defrag the thole_matricies every X steps.
	//     Also, as memory becomes scarce, it is neccessary to check the return values of
	// malloc, calloc and realloc to make sure we're not getting NULL's, particularly
	// for larger requests.


	// Allocates memory and throws an exception upon failure
	// Eg (to allocate an array of 100 ints):
	// int *myInt; 
	// SafeOps::calloc( myInt, 100, sizeof(int), __LINE__, __FILE__ );
	// __LINE__ and __FILE__ are preprocessor macros that are replaced with the line number
	// and file name of the line/file in which they appear. These will help us track down 
	// memory errors when they occur. 
	template<typename T>
	static void calloc( T &ptr, size_t qty, size_t size, int line, const char *file )
	{
		char linebuf[maxLine];

		// written to hopefully deal with some memory issues that have been plaguing very large
		// runs on PCN61. I believe that we may be running into memory fragmentation.
		// I am adding a function that will add an option defrag the thole_matricies every
		// X steps.
		// Also, as memory becomes scarce, it is neccessary to check the return values of
		// malloc, calloc and realloc to make sure we're not getting NULL's, particularly
		// for larger requests.
	
		
		if( (qty * size) <= 0 ) {	
			sprintf(linebuf, "[internal] Requested %d bytes: [%d] %s", (uint32_t)(qty * size), line, file );
			Output::err(linebuf);
			throw memory_request_invalid;
		}

		// Attempt to allocate the memory 
		ptr = (T) std::calloc(qty, size);
		
		// Check the allocation
		if (ptr == nullptr) {
			sprintf(linebuf, "[runtime system] Failed to allocate %lu bytes: [%d] %s", (long)(size * qty), line, file);
			Output::err(linebuf);
			throw memory_request_fail;
		}
	}



	// Allocates memory and throws an exception upon failure
	// Eg. (to allocate an array of 100 ints):
	// int *myInt; 
	// SafeOps::malloc( myInt, sizeof(int)*100, __LINE__, __FILE__ );
	// __LINE__ and __FILE__ are preprocessor macros that are replaced with the line number
	// and file name of the line/file in which they appear. These will help us track down 
	// memory errors when they occur. 
	template<typename T>
	static void malloc( T &ptr, size_t qty, int line, const char *file )
	{
		char linebuf[maxLine];
		

		if( qty <= 0 ) {	
			sprintf(linebuf, "[internal] Requested %ld bytes: [%d] %s", (long) qty, line, file );
			Output::err(linebuf);
			throw memory_request_invalid;
		}

		// Attempt to allocate the memory 
		ptr = (T) std::malloc( qty );
		
		// Check the allocation
		if( ptr == nullptr ) {
			sprintf(linebuf, "[runtime system] Failed to allocate %lu bytes: [%d] %s", (long) qty, line, file );
			Output::err(linebuf);
			throw memory_request_fail;
		}
	}


	// Shrinks/grows a memory block and throws an exception upon failure. 
	// Original block may be relocated and copied in some instances. 
	// Eg. (to resize an array to 100 ints):
	// int *myInt; 
	// SafeOps::realloc( myInt, sizeof(int)*100, __LINE__, __FILE__ );
	// __LINE__ and __FILE__ are preprocessor macros that are replaced with the line number
	// and file name of the line/file in which they appear. These will help us track down 
	// memory errors when they occur. 
	template<typename T>
	static void realloc( T &ptr, size_t qty, int line, const char *file )
	{
		char linebuf[maxLine];
		
		if( qty < 0 ) {	
			sprintf(linebuf, "[internal] Requested %ld bytes: [%d] %s", (long) qty, line, file );
			Output::err(linebuf);
			throw memory_request_invalid;
		}

		// Attempt to allocate the memory 
		T tempPtr = nullptr;
		tempPtr = (T) std::realloc( ptr, qty );
		
		// Check the allocation
		if(  (tempPtr == nullptr)   &&   qty  ) {
			sprintf(linebuf, "[runtime system] Failed to allocate %lu bytes: [%d] %s", (unsigned long) qty, line, file );
			Output::err(linebuf);
			throw memory_request_fail;
		}

		ptr = tempPtr;
	}

	// Frees memory if pointer is not null, and then assigns nullptr to said pointer. 
	template<typename T>
	static void free( T &pointer ) {
		if (pointer)
			std::free(pointer);
		pointer = nullptr;
	}

};


#endif // SAFEOPS_H