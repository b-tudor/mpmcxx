#pragma once
#ifndef OUTPUT_H
#define OUTPUT_H

using uint = size_t;

class Output
{
public:
	Output();
	~Output();

	// Immediately displays a message to stdout
	static void err( const char *msg);

	// Immediately displays a message to stderr
	static void out1( const char *msg );
	static void out(  const char *msg );
	
	// Forms a numbered filename, based on a common basename
	static char * make_filename( const char * basename, int fileno );
	
	static double calctimediff(struct timeval a, struct timeval b);

	static int GetTimeOfDay( struct timeval *tv );
	
};



#endif // OUTPUT_H