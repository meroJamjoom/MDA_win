
#include <sstream>
#include <math.h>

#include <time.h>

#include "MDA/Array/Array.hh"


#define BITS 8
#define NUM_PATTERNS (1<<BITS)
#define PATTERN_LEN (1<<BITS)


using namespace MDA;
using namespace std;


void
makeOnes( bool *p )
{
  for( int i= 0 ; i< NUM_PATTERNS ; i++ )
    for( int j= 0 ; j< PATTERN_LEN ; j++ )
      p[i*PATTERN_LEN+j]= true;
}


void
makeLinear( bool *p )
{
  int i, j, k, l;
  
  for( i= 0 ; i< NUM_PATTERNS ; i++ )
    for( j= k= 0 ; k< BITS ; k++ )
      for( l= 0 ; l< PATTERN_LEN/BITS; j++, l++ )
	p[i*PATTERN_LEN+j]= ((i&(1<<k)) > 0);
}

void
makeLogLSB( bool *p )
{
  int i, j, k, l;
  
  for( i= 0 ; i< NUM_PATTERNS ; i++ )
  {
    for( j= k= 0 ; k< BITS ; k++ )
      for( l= 0 ; l< (1<<k)*PATTERN_LEN/(1<<BITS) && j< PATTERN_LEN; j++, l++ )
	p[i*PATTERN_LEN+j]= ((i&(1<<k)) > 0);
    // log patterns are actually only 2^BITS-1 entries long, so copy last entry
    p[i*PATTERN_LEN+PATTERN_LEN-1]= p[i*PATTERN_LEN+PATTERN_LEN-2];
  }
}

void
makeLogMSB( bool *p )
{
  int i, j, k, l;
  
  for( i= 0 ; i< NUM_PATTERNS ; i++ )
  {
    for( j= k= 0 ; k< BITS ; k++ )
      for( l= 0 ; l< (1<<k)*PATTERN_LEN/(1<<BITS) && j< PATTERN_LEN; j++, l++ )
	p[i*PATTERN_LEN+(PATTERN_LEN-2-j)]= ((i&(1<<k)) > 0);
    // log patterns are actually only 2^BITS-1 entries long, so copy last entry
    p[i*PATTERN_LEN+PATTERN_LEN-1]= p[i*PATTERN_LEN+PATTERN_LEN-2];
  }
}


void
makeRandom( bool *p )
{
  int i, j, k, l;
  bool tmp;
#if defined (_WIN32) || defined (_WIN64)
#define srand48(num) srand(num)
#define lrand48() rand() 

#else
  srand48( (long)clock() );
  for( i= 0 ; i< NUM_PATTERNS ; i++ )
  {
    for( j= 0 ; j< PATTERN_LEN ; j++ )
      p[i*PATTERN_LEN+j]= (j<i);
    
    // random shuffle
    for( j= 0 ; j< PATTERN_LEN*10 ; j++ )
    {
      k= (int)(drand48() * PATTERN_LEN -.5);
      l= (int)(drand48() * PATTERN_LEN -.5);
      
      tmp= p[i*PATTERN_LEN+k];
      p[i*PATTERN_LEN+k]= p[i*PATTERN_LEN+l];
      p[i*PATTERN_LEN+l]= tmp;
    }
  }
#endif
}


int
main( int argc, char *argv[] )
{
  CoordinateVector dim( 2, NUM_PATTERNS );
  Array<float> out(dim);
  out.addChannel();
  
  bool pattern1[NUM_PATTERNS*PATTERN_LEN];
  bool pattern2[NUM_PATTERNS*PATTERN_LEN];
  
  makeRandom( pattern1 );
  makeRandom( pattern2 );
  
  float *data= &((*out[0])[0]);
  for( int i= 0 ; i < NUM_PATTERNS ; i++ )
    for( int j= 0 ; j < NUM_PATTERNS ; j++, data++ )
    {
      *data= 0.0;
      for( int k= 0 ; k < PATTERN_LEN ; k++ )
	*data+= (pattern1[i*PATTERN_LEN+k] && pattern2[j*PATTERN_LEN+k]);
      *data/= (double)PATTERN_LEN;
      //      *data= pattern1[i*PATTERN_LEN+j];
    }
  
  out.write();
}
