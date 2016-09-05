// ==========================================================================
// $Id:$
// insert a T-junction into the IO pipeline
// ==========================================================================
// License: Internal use at UBC only! External use is a copyright violation!
// ==========================================================================
// (C)opyright:
//
// 2007-, UBC
//
// Creator: atcheson (Bradley Atcheson)
// Email:   atcheson@cs.ubc.ca
// ==========================================================================

#include <stdlib.h>
#include <stdio.h>

#include <iostream>

#if defined (_WIN32) || defined (_WIN64)
#include <io.h>
#else
#include <unistd.h>
#endif

#define BUFSIZE 2048

using namespace std;


int
main(int argc, char *argv[])
{
	register FILE* subpipeline = NULL;

	// check arguments
	if (argc == 2)
	{
		if (*argv[1])
		{
			// create child process
#if defined (_WIN32) || defined (_WIN64)

			if (NULL == (subpipeline = _popen(argv[1], "w")))
#else
			if (NULL == (subpipeline = popen(argv[1], "w")))
#endif
			{
				cerr << argv[0] << " cannot create subpipeline " << argv[1] << endl;
				exit(1);
			}
		}
	}
	else if (argc > 2)
	{
		cerr << "Usage: " << argv[0] << " [pipeline]" << endl;
		exit(2);
	}

	// pass stdin to stdout, and copy to child process
	char buf[BUFSIZE];
	register unsigned int n; // num bytes read
#if defined (_WIN32) || defined (_WIN64)
	while (0 < (n = _read(0, buf, BUFSIZE)))
	{
		_write(1, buf, n); // write to stdout

		if (subpipeline)
		{
			_write((int)_fileno(subpipeline), buf, n);
		}
	}

#else
	while (0 < (n = read(0, buf, BUFSIZE)))
	{
		write(1, buf, n); // write to stdout

		if (subpipeline)
		{
			write((int)fileno(subpipeline), buf, n);
		}
	}
#endif

	if (subpipeline)
	{
#if defined (_WIN32) || defined (_WIN64)
		_pclose(subpipeline);
#else
		pclose(subpipeline);
#endif
	}

	return 0;
}

