#include <stdio.h>
#include <math.h>

#define PI acos(-1)

/*****************************************************************
***                          2D Neighbours                     ***
***                    Last Modified: 29/07/2024               ***
***                      Modified by: Gustavo                  ***
***                                                            ***
***  In a 2D, with periodic boundary conditions, create the    ***
***  matrixes with the neighbours of all sites.                ***
***                                                            ***
*****************************************************************/
void neighbours_2d(unsigned long *right,unsigned long *left,
		   unsigned long *up, unsigned long *down,
		   unsigned long lsize)

{
unsigned long i, l2;

l2 = lsize * lsize;
for (i = 0; i < l2; ++i)
    {
    if (i % lsize == lsize-1) *(right+i) = i-lsize+1;/* last col. */
    else *(right+i) = i+1;

    if (i % lsize == 0) *(left+i) = i+lsize-1;       /* first col.*/
    else *(left+i) = i-1;

    if (i < lsize) *(up+i) = l2-lsize+i;             /* first row */
    else *(up+i) = i - lsize;

    if (i >= l2-lsize) *(down+i) = (i % lsize);      /* last row  */
    else *(down+i) = i + lsize;
    }
return;
}

void generate_matrix(long *s1, int lsize)
{
unsigned long i;

for (i = 0; i < lsize; ++i)
    s1[i] = (float)rand()/(float)(RAND_MAX/360);
return;
}

float simulation(){
    int length         = 32;
    int total_mc_steps = 50000;
    // 2pi
    float theta        = 2 * PI * (330 / 360);

    float total_samples_in_step = length * length;

    return generate_matrix(length, length);
}

int main()
{
    float temperature = 0.4;
    printf(simulation(temperature));
}
