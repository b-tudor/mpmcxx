#include <stdio.h>
#include <cstring>

// (OS Dependent) Timing Includes
#if defined(WIN32) || defined(_WIN32) || defined(__WIN32)
#	include <time.h>
#	include <windows.h>
#else
#	include <sys/time.h>
#endif

#include "Output.h"
#include "constants.h"
#include "SafeOps.h"


extern int  rank;
extern bool mpi;




Output::Output() { }
Output::~Output() { }

void Output::err( const char *msg ) 
{

	if (!rank) fprintf(stderr, "(ERROR) %s", msg);
	fflush(stderr);

}

void Output::out1(const char *msg) {
	if (mpi && rank)
		return;
	else
		out(msg);
}
void Output::out( const char *msg )
{
	printf("%s", msg);
	fflush(stdout);
}

char * Output::make_filename( const char * basename, int fileno ) 
{
	size_t len = std::strlen(basename);
	size_t outlen;
	char * rval = nullptr;
	char START[maxLine], STOP[maxLine];

	if(  !strncmp("/dev/null", basename, 9)  ) {
		SafeOps::calloc(rval, 10, sizeof(char), __LINE__, __FILE__ );
		sprintf(rval,"/dev/null");
		return rval;
	}
	else 
	{
	//check if the file has a three character extension
		if ( len - 4 > 0 ) 
		{
			if ( basename[len-4] == '.' ) //only check this if previous check passes so we don't have a memory fault
			{
				size_t i;
				size_t j;
				for ( i=0; i<(len-4); i++ ) //format pre-extension
					START[i] = basename[i];
				START[i]='\0';
				for ( i = (len-4), j=0; i<len; i++, j++ ) //format extension
					STOP[j] = basename[i];
				STOP[j]='\0';
				//set string length
				outlen = std::strlen(START) + std::strlen(STOP)+7;
				SafeOps::calloc( rval, outlen, sizeof(char), __LINE__, __FILE__ );
				//make filename
				sprintf(rval,"%s-%04d%s", START, fileno, STOP);
			}
		}
	}

	//if rval is still NULL, then it's neither /dev/null nor has a proper 3-character file extension
	//just add the number to the end
	if( rval == nullptr ) {
		outlen = len + 7;
		SafeOps::calloc( rval, outlen, sizeof(char), __LINE__, __FILE__ );
		//generate filename
		sprintf(rval,"%s-%04d", basename, fileno);
	}
	
	return rval;	
}


double Output::calctimediff(struct timeval a, struct timeval b) {
	return (size_t)a.tv_sec - b.tv_sec + 1.0e-6 * ((size_t)a.tv_usec - b.tv_usec);
}



#if defined(WIN32) || defined(_WIN32) || defined(__WIN32)
#	include <time.h>
#	include <windows.h>
#	if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
#		define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
#	else
#		define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
#	endif
	
	

	int Output::GetTimeOfDay(struct timeval *tv )
	{
		// Thank you: https://gist.github.com/ugovaretto/5875385#file-win-gettimeofday-c-L17

		FILETIME         ft;
		unsigned __int64 tmpres = 0;
		static int       tzflag = 0;

		if( tv != nullptr )
		{
			GetSystemTimeAsFileTime(&ft);

			tmpres |= ft.dwHighDateTime;
			tmpres <<= 32;
			tmpres |= ft.dwLowDateTime;

			tmpres /= 10;  //convert into microseconds
			//converting file time to unix epoch
			tmpres -= DELTA_EPOCH_IN_MICROSECS; 
			tv->tv_sec = (long)(tmpres / 1000000UL);
			tv->tv_usec = (long)(tmpres % 1000000UL);
		}
		return 0;
	}
#else
#	include <sys/time.h>
	int Output::GetTimeOfDay(struct timeval *tv) {
		return gettimeofday( tv, nullptr );
	}
#endif