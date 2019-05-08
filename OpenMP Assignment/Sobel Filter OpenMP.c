#include<math.h>
#include<omp.h>
#include<stdio.h>
#include<string.h> 
#include<time.h>
#include<stdlib.h>



#define rand 1111
typedef double ttype;

ttype tdiff(struct timespec a, struct timespec b)
{
  ttype dt = (( b.tv_sec - a.tv_sec ) + ( b.tv_nsec - a.tv_nsec ) / 1E9);    //Finding the time difference (copied from sample file on e-learning)
  return dt;
}
struct timespec now()
{
  struct timespec t;
  clock_gettime(clock, &t);
  return t;
}

int main(int argc, char *argv[])
	{
		int i,j,k,l;
		int p = 3840 ,q = 2160;
		srand(rand);
		struct timespec begin, end;
        double time_spent;
  		int num_threads = atof(argv[1]);
		int r = p +4;
		int s = q+4;
	
	float Arr1[p][q];
	float x[p+2][q+2], y[p+2][q+2];
	float out[p][q];
	float edge[p][q];
	memset(out,0,(p+2)*(q+2)*sizeof(int));
	memset(edge,0,(p+2)*(q+2)*sizeof(int));
	memset(x,0,(p+2)*(q+2)*sizeof(int));
	memset(y,0,(p+2)*(q+2)*sizeof(int));
	int ker_x[3][3] = { {1,0,-1},{2,0,-2},{1,0,-1}};
	int ker_y[3][3] = { {-1,-2,-1},{0,0,0},{1,2,1}};
	int iker_x[3][3] = {0};
	int iker_y[3][3] = {0};
	float B[r][s];
	memset(B,0,r*s*sizeof(int));
    omp_set_num_threads(num_threads); 
    printf("Number of threads:%d \n", num_threads);
	for(i=0;i<p;i++)
	{
		for(j=0;j<q;j++)
		{
			Arr1[i][j] = (rand() % (256));      	
		}
	}
  	for(i=0;i<p;i++)
	  {
		for(j=0;j<q;j++)
		{
			B[i+2][j+2]=Arr1[i][j];
		}
	}
	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			iker_x[i][j] = ker_x[2-i][2-j];
			iker_y[i][j] = ker_y[2-i][2-j];	
		}
	}
	printf("before x \n");
    begin = now();
	for(i=0;i<p+2;i++){
		 #pragma omp parallel for 
		for(j=0;j<q+2;j++){
    
			for(k=0;k<3;k++){
				for(l=0;l<3;l++){
					x[i][j]+= (B[i+k][j+l]*iker_x[k][l]);
          y[i][j]+= (B[i+k][j+l]*iker_y[k][l]);
          out[i][j] = sqrt(powf(x[i][j],2) + powf(y[i][j],2));
				}
			}
		}
	}
 printf("before end \n");
	 end = now();
        time_spent = tdiff(begin, end);
        printf("Total execution time is %.8f s\n", time_spent);
}
	
	
