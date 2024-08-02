/********************************************************************
***                       Monte Carlo Routines                    ***
***                         V0.11 18/05/2000                      ***
***                                                               ***
***                                                               ***
********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h> 
#include <string.h>
#include <time.h>  


#define MC_VERSION  "0.11"
#define FNORM   (2.3283064365e-10)
#define RANDOM  ( (ira[ip++]=ira[ip1++]+ira[ip2++]) ^ira[ip3++] )        
#define FRANDOM (FNORM * RANDOM)

/********************************************************************
***                      Variable Declarations                    ***
********************************************************************/

unsigned long randomize, seed;                    /* Random number */
static clock_t   begin_cpu_time;                 /* Time variables */
static time_t    begin_real_time;
unsigned zseed, ira[256];
unsigned char ip,ip1,ip2,ip3;

/********************************************************************
*            Random Number Generator by Parisi & Rapuano            *
*                  Last Modified: 18/05/2000                        *
*                                                                   *
*  First, the function start_randomic() should be called to create  *
*  the seed (odd number):                                           *
*                 seed = start_randomic();                          *
*                                                                   *
*  Return: start_randomic -> unsigned long int                      *
*          FRANDOM -> double in [0,1)                               *
********************************************************************/
unsigned rand4init(void)
{
     unsigned long long y;
   
     y = (zseed*16807LL);
     zseed = (y&0x7fffffff) + (y>>31);
     if (zseed&0x80000000)
         zseed = (zseed&0x7fffffff) + 1;
     return zseed;
}

void Init_Random(void)
{
     int i;
   
     ip=128;
     ip1=ip-24;
     ip2=ip-55;
     ip3=ip-61;
   
     for (i=ip3; i<ip; i++)
         ira[i] = rand4init();
}

/****************************************************************
*                          Random Seed                          *
*                  Last Modified: 18/05/2000                    *
*                                                               *
* When debugging, use always the same seed, otherwise, take it  *
* as the numbers of seconds ellapsed since ?/?/1970.            *
****************************************************************/ 
unsigned long start_randomic(void)

{
  unsigned long semente;
   
  semente = (unsigned long) time(NULL);
  if (semente%2==0) ++semente; /*odd number*/ /* random seed */
#ifdef DEBUG
  semente = 1206544417 /* 1189007899 */ /*2067511789*//*2067512563*//*1852686583*/;
#endif
  zseed = semente;
  Init_Random();
  return semente;
}

/****************************************************************
*                          Random Seed                          *
*                  Last Modified: 18/05/2000                    *
*                                                               *
* When debugging, use always the same seed, otherwise, take it  *
* as the numbers of seconds ellapsed since ?/?/1970.            *
****************************************************************/ 
unsigned long start_randomic_seed_ext(unsigned long semente)

{
  //unsigned long semente;
   
  //semente = (unsigned long) time(NULL);
  if (semente%2==0) ++semente; /*odd number*/ /* random seed */
#ifdef DEBUG
  semente = 1206544417 /* 1189007899 */ /*2067511789*//*2067512563*//*1852686583*/;
#endif
  zseed = semente;
  Init_Random();
  return semente;
}


/****************************************************************
*                          Random Seed                          *
*                  Last Modified: 01/04/2008                    *
*                                                               *
* When debugging, use always the same seed, otherwise, take it  *
* as the numbers of seconds ellapsed since ?/?/1970.            *
****************************************************************/ 
void start_randomic_parallel(unsigned long semente)

{
  zseed = semente;
  Init_Random();
  return;
}

/*****************************************************************
***                          2D Neighbours                     ***
***                    Last Modified: 24/01/1999               ***
***                                                            ***
***  In a 2D, with periodic boundary conditions, create the    ***
***  matrixes with the neighbours of all sites.                ***
***                                                            ***
***  Input: right,left,up & down -> matrixes with neighb. sites***
***         l_size -> linear dimension (l_size x l_size)       ***
***                                                            ***
**   From right to left, top to bottom, front to back.         ***
*****************************************************************/
void neighbours_2d(unsigned long *right,unsigned long *left,
		   unsigned long *up, unsigned long *down,
		   unsigned long lsize)

{
unsigned long i,l2;
   
l2 = lsize*lsize;
for (i=0; i<l2; ++i)     
    {
     if (i % lsize==lsize-1) *(right+i) = i-lsize+1;/* last col. */
                        else *(right+i) = i+1; 
     if (i % lsize==0) *(left+i) = i+lsize-1;       /* first col.*/
                  else *(left+i) = i-1;
     if (i<lsize) *(up+i) = l2-lsize+i;             /* first row */
             else *(up+i) = i - lsize;
     if (i>=l2-lsize) *(down+i) = (i % lsize);      /* last row  */
                 else *(down+i) = i + lsize;
    }  
return;
}
/*****************************************************************
***                    2D Hexagonal Neighbours                 ***
***                    Last Modified: 13/08/2012               ***
***                                                            ***
***  In a 2D, with periodic boundary conditions, create the    ***
***  matrixes with the neighbours of all sites. (Honeycomb)    ***
***                                                            ***
***  Input: right,left,up1,up2,down1 & down2 -> matrixes with  ***
***     neighb. sites                                          ***  
***         l_size -> linear dimension (l_size x l_size)       ***
***                                                            ***
**   From right to left, top to bottom,                        ***
*****************************************************************/
void neighbours_hexagonal_2d(unsigned long *right,unsigned long *left,
		   unsigned long *up1, unsigned long *down1, 
		   unsigned long *up2, unsigned long *down2,
		   unsigned long lsize)

{
unsigned long i,l2;

 if (lsize%2==1) 
   {
     printf("For honeycomb lattice, LSIZE must be even.\n");
     exit(8);
   }
 else
   {
     unsigned long *up,*down;
     
     l2 = lsize*lsize;
     up=(unsigned long *)calloc(l2,sizeof(unsigned long));
     down=(unsigned long *)calloc(l2,sizeof(unsigned long));
     if ((up!=NULL) && (down!=NULL))
       {
	 /*Regular Square lattice*/
	 for (i=0; i<l2; ++i)     
	   {
	     if (i % lsize==lsize-1) *(right+i) = i-lsize+1;/* last col. */
	     else *(right+i) = i+1; 
	     if (i % lsize==0) *(left+i) = i+lsize-1;       /* first col.*/
	     else *(left+i) = i-1;
	     if (i<lsize) *(up+i) = l2-lsize+i;             /* first row */
	     else *(up+i) = i - lsize;
	     if (i>=l2-lsize) *(down+i) = (i % lsize);      /* last row  */
	     else *(down+i) = i + lsize;
	   }  
	  /*Modification to hexagonal lattice*/
	  for (i=0; i<l2; ++i)     
	    {
	      if ((i/lsize)%2==1) 
		{
		  *(up1+i)=*(up+*(left+i));
		  *(up2+i)=*(up+i);
		  *(down1+i)=*(down+*(left+i));
		  *(down2+i)=*(down+i);
		}
	      else
		{
		  *(up1+i)=*(up+i);
		  *(up2+i)=*(up+*(right+i));
		  *(down1+i)=*(down+i);
		  *(down2+i)=*(down+*(right+i));
		}
	    }
	 free(up);
	 free(down);
       }
     else 
       {
	 fprintf(stderr,"Error in routine neighbour_hexagonal_2d\n");
	 exit(8);
       }
   }
 return;
}
/*****************************************************************
***                          3D Neighbours                     ***
***                    Last Modified: 23/03/1999               ***
***                                                            ***
***  In a 3D, with periodic boundary conditions, create the    ***
***  matrixes with the neighbours of all sites.                ***
***                                                            ***
***  Input: right,left,up,down,front, back -> neighb. sites    ***
***         lsize -> linear dimension (lsize^3)                ***
*****************************************************************/
void neighbours_3d(unsigned long *right,unsigned long *left,
		   unsigned long *up, unsigned long *down,
		   unsigned long *front, unsigned long *back,
		   unsigned long lsize)

{
unsigned long i,l2,l3;

l2 = lsize*lsize;
l3 = l2*lsize;
for (i=0; i<l3; ++i) 
    {
     if (i % lsize==lsize-1) *(right+i) = i - lsize + 1; /* righmost plane */
                        else *(right+i) = i+1;
     if (i % lsize==0) *(left+i) = i + lsize - 1;        /* leftmost plane */
                  else *(left+i) = i-1;
     if (i % l2 < lsize) *(up+i) = l2 - lsize + i;       /* top plane */
                    else *(up+i) = i - lsize;
     if (i % l2 >= (l2-lsize)) *(down+i) = i + lsize - l2;  /* bottom plane */
                          else *(down+i) = i + lsize;
     if (i < l2) *(front+i) = l3 - l2 + i;             /* frontal plane */
            else *(front+i) = i - l2;
     if (i >= (l3-l2)) *(back+i) = i % l2;
                  else *(back+i) = i + l2;             /* back plane */
    }
return;
}         
/*****************************************************************
***                          3D Neighbours -INT                ***
***                    Last Modified: 23/04/2014               ***
***                                                            ***
***  In a 3D, with periodic boundary conditions, create the    ***
***  matrixes with the neighbours of all sites.                ***
***                                                            ***
***  Input: right,left,up,down,front, back -> neighb. sites    ***
***         lsize -> linear dimension (lsize^3)                ***
*****************************************************************/
void neighbours_3d_uint(unsigned int *right,unsigned int *left,
		   unsigned int *up, unsigned int *down,
		   unsigned int *front, unsigned int *back,
		   unsigned int lsize)

{
unsigned int i,l2,l3;

l2 = lsize*lsize;
l3 = l2*lsize;
for (i=0; i<l3; ++i) 
    {
     if (i % lsize==lsize-1) *(right+i) = i - lsize + 1; /* righmost plane */
                        else *(right+i) = i+1;
     if (i % lsize==0) *(left+i) = i + lsize - 1;        /* leftmost plane */
                  else *(left+i) = i-1;
     if (i % l2 < lsize) *(up+i) = l2 - lsize + i;       /* top plane */
                    else *(up+i) = i - lsize;
     if (i % l2 >= (l2-lsize)) *(down+i) = i + lsize - l2;  /* bottom plane */
                          else *(down+i) = i + lsize;
     if (i < l2) *(front+i) = l3 - l2 + i;             /* frontal plane */
            else *(front+i) = i - l2;
     if (i >= (l3-l2)) *(back+i) = i % l2;
                  else *(back+i) = i + l2;             /* back plane */
    }
return;
}         
/*****************************************************************
***                          4D Neighbours                     ***
***                    Last Modified: 09/05/2002               ***
***                                                            ***
***  In a 4D, with periodic boundary conditions, create the    ***
***  matrixes with the neighbours of all sites.                ***
***                                                            ***
***  Input: right,left,up,down,front,back,                     ***
***             front4d,back4d -> neighb. sites                ***
***         lsize -> linear dimension (lsize^4)                ***
*****************************************************************/
void neighbours_4d(unsigned long *right,unsigned long *left,
		   unsigned long *up, unsigned long *down,
		   unsigned long *front, unsigned long *back,
		   unsigned long *front4d, unsigned long *back4d,
		   unsigned long lsize)

{
  unsigned long i,l2,l3,l4;

  l2 = lsize*lsize;
  l3 = l2*lsize;
  l4 = l2*l2;
  for (i=0; i<l4; ++i) 
    {
      if (i % lsize==lsize-1) *(right+i) = i - lsize + 1; /* righmost plane */
                         else *(right+i) = i+1;
      if (i % lsize==0) *(left+i) = i + lsize - 1;        /* leftmost plane */
                   else *(left+i) = i-1;
      if (i % l2 < lsize) *(up+i) = l2 - lsize + i;       /* top plane */
                     else *(up+i) = i - lsize;
      if (i % l2 >= (l2-lsize)) *(down+i) = i + lsize - l2;  /* bottom plane */
                           else *(down+i) = i + lsize;
      if (i % l3 < l2) *(front+i) = l3 - l2 + i;           /* frontal plane */
                  else *(front+i) = i - l2;
      if (i % l3 >= (l3-l2)) *(back+i) = i - (lsize-1)*l2;
                   else *(back+i) = i + l2;             /* back plane */
      if (i < l3) *(front4d+i) = l4 - l3 + i;          /* 'front' cube */
             else *(front4d+i) = i - l3;
      if (i >=(l4-l3)) *(back4d+i) = i % l3;           /* 'back' cube */
                  else *(back4d+i) = i + l3;
    }
return;
}         
/*******************************************************************************
*                   Instantaneous Percolation Measurements                     *
*                     Last modified:  13/06/2005                               *
*                                                                              *
* Remember that a cluster is identified by the last site in the chain of       *
* pointers.                                                                    *
*******************************************************************************/
void unionfind(unsigned long i,unsigned long j, unsigned long *lab)
          
{
unsigned long i1,j1;

i1 = lab[i];                            /* check where i points to          */
while (lab[i1] != i1)                   /* while it doesn't point to itself */
      i1 = lab[i1];

j1 = lab[j];                            /* check where j points to          */
while (lab[j1] != j1)                   /* while it doesn't point to itself */
      j1 = lab[j1];

if (lab[i1] > lab[j1]) lab[i1] = lab[j1];
                  else lab[j1] = lab[i1];
return;
}

/**********************************************************************
***                         Hamming Distance                        ***
***                   Last modified: 24/01/1999                     ***
***    Returns the number of different sites between two vectors.   ***
**********************************************************************/ 
int hamming_distance(int *s1, int *s2, int lsize)

{
unsigned long i,temp=0;
  
for (i=0; i<lsize; ++i)
     if (s1[i] != s2[i]) ++temp;
return temp;
}

/*********************************************************************
***                          Magnetization                         ***
***                   Last Modified: 30/05/1999                    ***
***                                                                ***
*** Return: LSIZE*magnetization                                    ***
*********************************************************************/
int magnetization(int *s, int lsize)

{
unsigned long i,temp=0;
   
for (i=0; i<lsize; ++i)
    temp += s[i];
return temp;
}

int magnetization_diluted(int *s, int *n, int lsize)

{
unsigned long i, temp=0;
   
for (i=0; i<lsize; ++i)
    temp += (s[i]*n[i]);
return temp;
}

/*********************************************************************
***                              Overlap                           ***
***                    Last Modified: 30/05/1999                   ***
***                                                                ***
*** Return: LSIZE*overlap                                          ***
*********************************************************************/
int overlap(int *s1, int *s2, int lsize)

{
unsigned long i,temp=0;
   
for (i=0; i<lsize; ++i)
    temp += s1[i]*s2[i];
return temp;
}                       

int overlap_diluted(int *s1, int *n1, int *s2, int *n2, int lsize)

{
unsigned long i,temp=0;
   
for (i=0; i<lsize; ++i)
    temp += s1[i]*n1[i]*s2[i]*n2[i];
return temp;
}          

/********************************************************************
***                   Gaussian Random Number Generator            ***
***                     Last Modified: 18/05/2000                 ***
********************************************************************/ 
float ngaussian(void)

{
static int iset=0;
static float gset;
float fac,r,v1,v2;
   
if (iset==0) {
	      do {
	  	  v1=2.0*FRANDOM-1.0;
		  v2=2.0*FRANDOM-1.0;
		  r=v1*v1+v2*v2;
	         } 
	      while (r>=1 || r==0.0);
	      fac=sqrt(-2.0*log(r)/r);
	      gset=v1*fac;
	      iset=1;
	      return v2*fac;
             }
        else {
	      iset=0;
	      return gset;
             }
}

/*********************************************************************
***                Create Initial Bond Configuration               ***
***                    Last modified: 18/05/2000                   ***
***                                                                ***
***  The (symmetric) bonds are created depending on the model      ***
***  up to dimension 3:                                            ***
***    0. Ferromagnet: all bonds = 1                               ***
***    1. Spin Glass: same probability of being +1 and -1          ***
***    2. Spin Glass: gaussian distributed                         ***
***    3. Fully Frustrated: horizontal connections = +1 and        ***
***       vertical ones alternating signals such that all minimal  ***
***       loops are frustrated (LSIZE must be EVEN!)               ***
***       2D ONLY!!                                                ***
***                                                                ***
***    4. Fully Frustrated (Federico)                              ***
***                                                                ***
***  The first index follows the rule:                             ***
***    0. Horizontal coupling;                                     ***
***    1. Vertical coupling;                                       ***
***    2. Back coupling;                                           ***
***                                                                ***
***  Use 3 for the FF-2D!                                          ***
*********************************************************************/
void create_bonds(int model, int coordination, int lsize, 
		  double *connect)

{
  unsigned long i,k,ll,x,y,z; //l2;

ll = lsize;
// l2 = lsize*lsize;
switch (coordination)
       {
   	case 2: ll *= lsize;  /* 2D */
	        break;
        case 3: ll *= lsize*lsize; /* 3D */
	        break;
       default: fprintf(stdout,"Dimension should be either 2 or 3!\n");
	        exit(1);
	        break;
      }

for (k=0; k<coordination; ++k)
     for (i=0; i<ll; ++i)
       switch (model){
	      case 0: *(connect+k*ll+i) = 1;
	              break;
	      case 1: if (FRANDOM<0.5) *(connect+k*ll+i) = 1;
	                             else *(connect+k*ll+i) = -1;
		      break;
	      case 2: *(connect+k*ll+i) = ngaussian();
	              break;
	      case 3: *(connect+i) = 1; /* connect[0][i] */
	              if ((i % lsize) % 2 == 0) *(connect+ll+i) = -1;
	                                   else *(connect+ll+i) = 1;
     		      break;
              case 4: x = i % lsize; /* From Federico */
                      y = ((i - x) / lsize) % lsize;
                      z = i / (lsize*lsize);
                      *(connect+i) = 1;
                      *(connect+ll+i) = 1;
                      *(connect+2*ll+i) = 1;
                      if ((y % 2 == 0) && (z % 2 == 1))
                         *(connect+i) = -1;
                      if ((x % 2 == 0) && (z % 2 == 0))
                         *(connect+ll+i) = -1;
                      if ((x % 2 == 1) && (y % 2 == 1))
                         *(connect+2*ll+i) = -1;
		      break;
	     default: fprintf(stdout,"Bond set not defined!\n");
		      exit(1);
	              break;
       }  
return;
}


void create_bonds_a(double aa, int coordination, int lsize, double *connect)

{
unsigned long i,k,ll;

ll = lsize*lsize;
if (coordination==3) ll *= lsize; 
   else if (coordination!=2)
           fprintf(stdout,"Dimension should be either 2 or 3!\n");

for (k=0; k<coordination; ++k)
     for (i=0; i<ll; ++i)
       {
	  if (FRANDOM<aa) *(connect+k*ll+i) = 1;
	             else *(connect+k*ll+i) = -1;
       }  
return;
}

/*********************************************************************
***                     Create Initial Configuration               ***
***                      Last modified: 25/04/2000                 ***
***                                                                ***
*** Create at random an initial configuration for both S_i and n_i.***
*********************************************************************/
void initial_configuration(float mag, float dens, int ll, 
			   int *s, int *n)

{
unsigned long i;
float rm;
   
rm = 0.5*(1.0 - mag);  /* fraction of -1 sites */
if (dens!=0.0)
     for (i=0; i<ll; ++i)
       {
	  if (s != NULL)
	    {
	       if (FRANDOM<rm) s[i] = -1;
	       else s[i] =  1;
	    }
	  if (n != NULL)
	    {
	       if (FRANDOM<dens) n[i] = 1;
	       else n[i] = 0;
	    }
       }
   return;
}

/*********************************************************************
***                    Measuring Performance (time)                ***
***                     Last modified: 24/03/1999                  ***
***                                                                ***
*** Use:                                                           ***
***  cpu_time = (clock() - begin_clock) / (double) CLOCKS_PER_SEC; ***
***  real_time = difftime(time(NULL), begin_time);                 ***
*********************************************************************/ 
void start_stopwatches(void)

{
begin_cpu_time = clock();
begin_real_time = time(NULL);
return;
}

/*********************************************************************
***                     Simple Matrix Operations                   ***
***                    Last Modified: 18/05/2000                   ***
***                                                                ***
*********************************************************************/
void copymatrix(int *s1, int *s2, int lsize)
{
unsigned long i;
   
for (i=0; i<lsize; ++i)
    s1[i]=s2[i];
return;
}               

void copymatrix_f(float *s1, float *s2, int lsize)
{
  unsigned long i;
  
  for (i=0; i<lsize; ++i)
    s1[i]=s2[i];
  return;
}

void zeromatrix_i(int *s1, int lsize)
{
  
  unsigned long i;
  
  for (i=0; i<lsize; ++i)
    s1[i]=0;
  return;
}

void zeromatrix_l(long *s1, int lsize)
{
   
unsigned long i;
   
for (i=0; i<lsize; ++i)
    s1[i]=0;
return;
}

void zeromatrix_ul(unsigned long *s1, int lsize)
{
   
  unsigned long i;
  
  for (i=0; i<lsize; ++i)
    s1[i]=0;
  return;
}


void zeromatrix_sl(signed long *s1, int lsize)
{
unsigned long i;
   
for (i=0; i<lsize; ++i)
    s1[i]=0;
return;
}

void zeromatrix_f(float *s1, int lsize)
{
  unsigned long i;
  
  for (i=0; i<lsize; ++i)
    s1[i]=0.0;
  return;
}


void zeromatrix_d(double *s1, int lsize)
{
unsigned long i;
   
for (i=0; i<lsize; ++i)
    s1[i]=0.0;
return;
}

/**********************************************************************
***                          Create Time Table                      ***
***                    Last Modified: 04/08/1999                    ***
**********************************************************************/
void create_time_table(long *t1,long *t2,int power, int base)

{
long i,j;
   
t1[0] = 1;
for (i=1; i<power; i++)
     t1[i] = t1[i-1]*base;
if (t2 != NULL)
   {
     for (i=0; i<power; i++)
       for (j=0; j<power; j++)
	 {
	   *(t2+i*power+j) = t1[i] + t1[j];
	   if (*(t2+i*power+j) == base*t1[power-1]) *(t2+i*power+j) =base*t1[power-1];
	   if (*(t2+i*power+j) > base*t1[power-1]) *(t2+i*power+j) =-1;
	   
	 }
   }
return;
}                    

/*********************************************************************
***                         Time Table 2                           ***
***                    Last Modified: 08/01/2000                   ***
***  The total number of time steps and the number of measures are ***
***  specified.                                                    ***
*********************************************************************/
void create_time_table_2(long *t1, int total_time, int measures)

{
unsigned long i,k;
double temp;
   
t1[0] = 0;
temp = pow((double) total_time,1.0/(measures-1));
k=0;
for (i=1; i<measures; ++i)
    {
     t1[i] = (int) pow(temp, (double) i);
     if (t1[i]<=k) t1[i]=k+1;
     k = t1[i];
    }
return;
}


/*********************************************************************
***                         Time Table                             ***
*********************************************************************/
void create_time_table_fa(float *t1, float *t2, int decades)
  
{
  
  unsigned long j,k=1;
  signed int i=-3;
  
  t1[0] = 0;
  while (i <= decades)
    {
      
      t1[k] = pow(10,i-1);
      ++k;
      for (j=2; j<10; ++j)
	{
	  t1[k] = t1[k-1] + pow(10,i-1);
	  ++k;
	}
      ++i;
    }
  t1[k] = pow(10,i-1);
  
  if (t2 != NULL)
    {
      
      for (i=0; i<decades+5; i++)
        for (j=0; j<(decades+4)*9+2; j++)
	  {
	    if (i==0) *(t2+i*((decades+4)*9+2)+j) = t1[j];
	    else 
	      *(t2+i*((decades+4)*9+2)+j) = t1[9*i+1] + t1[j];
	    /*if (*(t2+i*10+j) == base*t1[power-1]) *(t2+i*power+j) =base*t1[power-1];
	      if (*(t2+i*10+j) > base*t1[power-1]) *(t2+i*power+j) =-1;*/
	  }
    }


  return;
}

