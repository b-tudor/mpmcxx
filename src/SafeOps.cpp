#include "SafeOps.h"

#include <chrono>
#include <cstring>
#include <string>
#include <thread>



// Opens files and throws exception upon failure
FILE * SafeOps::openFile( const char *inFilename, const char *mode, int line, const char *file ) {

	FILE *fp;
	char linebuf[maxLine];

	if (mode[0] != 'r' && mode[0] != 'w' && mode[0] != 'a') {
		sprintf(linebuf, "[internal] Invalid file mode specified (\"%s\" in %s (line %d)", mode, file, line );
		Output::err(linebuf);
		throw internal_error;
	}

	fp = fopen(inFilename, mode);
	if (fp == nullptr) {
		std::this_thread::sleep_for(std::chrono::seconds(1));
		fp = fopen(inFilename, mode);
	}

	if( fp == nullptr ) {

		int errorCode = 0;

		switch (mode[0]) {
		case 'r':
			sprintf(linebuf, "Failed to open file (for read): %s\n", inFilename );
			errorCode = fopen_fail_read;
			break;
		case 'w':
			sprintf(linebuf, "Failed to open file (for write): %s\n", inFilename );
			errorCode = fopen_fail_write;
			break;
		case 'a':
			sprintf(linebuf, "Failed to open file (for append): %s\n", inFilename);
			errorCode = fopen_fail_append;
			break;
		default:
			sprintf(linebuf, "Failed to open file: \"%s\"\n", inFilename);
			errorCode = fopen_fail_unknown;
		}

		Output::err(linebuf);
		throw errorCode;
	}

	return fp;
}



void SafeOps::closeFile(FILE *fp) {
	if( fp )
		fclose( fp );
	fp = nullptr;
}


 // Safe string-to-data converters
///////////////////////////////////////////////////////////////////////////////

//convert string *a to a double and store at *d.
bool SafeOps::atod( const char * s, double &d) {
	size_t idx;
	std::string myDoubleString(s);
	try {
		d = std::stod(myDoubleString, &idx);
		if(idx != myDoubleString.length())     // Causes program to reject input like 123.4m1 or 100.o
			return fail;
	}
	catch (...) {
		return fail;
	}
	return ok;
}

//convert string *a to an int and store at *i.
bool SafeOps::atoi( const char * s, int &i) {
	size_t idx;
	std::string myIntString(s);
	try {
		i = std::stoi(myIntString, &idx);
		if(idx != myIntString.length())     // Causes program to reject input like 1000m00 or 5b
			return fail;
	}
	catch (...) {
		return fail;
	}
	return ok;
}

//convert string *a to an unsigned int and stores at *u.
bool SafeOps::atou( const char * s, uint32_t &u) {
	long l;
	size_t idx;
	std::string myUInt32String(s);
	try {
		l = std::stol(myUInt32String, &idx);
		if(idx != myUInt32String.length())     // Causes program to reject input like 1000m00 or 5b
			return fail;
	}
	catch (...) {
		return fail;
	}
	if( l>=0 && l<=UINT32_MAX )
		u = (uint32_t) l;
	return ok;

}

//converts string *a to a long and stores at *l.
bool SafeOps::atol( const char * s, long &l) {
	size_t idx;
	std::string myLongString(s);
	try {
		l = std::stol(myLongString, &idx );
		if(idx != myLongString.length())     // Causes program to reject input like 1000m00 or 5b
			return fail;
	}
	catch (...) {
		return fail;
	}
	return ok;
}

bool SafeOps::strncasecmp(const char *A, const char *B, size_t n ) {

	char *a, *b;
	bool returnVal;

	size_t ALen = strlen(A);
	size_t BLen = strlen(B);
	SafeOps::calloc( a, n, sizeof(char), __LINE__, __FILE__ );
	SafeOps::calloc( b, n, sizeof(char), __LINE__, __FILE__ );

	for( size_t i = 0; i < n; i++ ) {
		if( i < ALen ) 
			a[i] = (char) tolower(A[i]);
		else
			a[i] = (char) 0;

		if( i < BLen ) 
			b[i] = (char) tolower(B[i]);
		else
			b[i] = (char) 0;
	}

	returnVal = ( strcmp(a,b) == 0 );

	free(a);
	free(b);

	return returnVal;
}



bool SafeOps::iequals(const char *A, const char *B) {

	int A_length = 0,
	    B_length = 0;

	// Get length of string A, set to 0 if pointer to string is null
	if( nullptr == A )
		A_length = 0;
	 else
		A_length = (int) strlen(A);

	// Get length of string B, set to 0 if pointer to string is null
	if( nullptr == B )
		B_length = 0;
	else 
		B_length = (int) strlen(B);
	
	// If both strings are non-existent, we say they are equal
	if( A_length==0 && B_length==0 )
		return true;
	
	// If they have different lengths, we don't bother to compare individual characters
	if( A_length != B_length )
		return false;
	
	for( int i=0; i<A_length; i++ ) {
		// If any letter is different (disregarding case), return false
		if( toupper(A[i]) != toupper(B[i]) )
			return false;
	}

	return true;
}