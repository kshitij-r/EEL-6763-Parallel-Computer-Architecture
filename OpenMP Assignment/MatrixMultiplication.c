#include<omp.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
typedef double ttype;
#define rand 1111
ttype tdiff(struct timespec a, struct timespec b)                    

{
  ttype dt = (( b.tv_sec - a.tv_sec ) + ( b.tv_nsec - a.tv_nsec ) / 1E9);   //Finding the time difference (copied from sample file on e-learning)
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
    int i,j,k;
    srand(rand);
    struct timespec begin, end;
    double time_spent;
    int a1,a2,a3,a4;;
	int num_threads = atof(argv[1]);
  a1 = atof(argv[2]);
  a2 = atof(argv[3]);
  a3 = atof(argv[4]);
  a4 = atof(argv[5]);
 
	int **A =(int **)malloc(a1 * sizeof(int*));   //dynamically allocating space for matrix A 
	for(i=0;i<a1;i++)
    {
		A[i] = (int *)malloc(a2 * sizeof(int));
	}
	for(i=0;i<a1;i++)
    {
		for(j=0;j<a2;j++)
        {
			A[i][j]= rand();                     // Initializing matrix A with random values
		}
	}
	

	int **B = (int **)malloc(a3 * sizeof(int*));  //dynamically allocating space for matrix B
	for(i=0;i<a3;i++)
    {
		B[i]=(int *)malloc(a4*sizeof(int));
	}
	for(i=0;i<a3;i++)
    {
		for(j=0;j<a4;j++)
        {
			B[i][j]=rand();                       // Initializing matrix A with random values
		}
	}
	
	
	int **C=(int **)malloc(a1*sizeof(int*));     //dynamically allocating space for matrix C
	for(i=0;i<a1;i++)
    {
		C[i]=(int *)malloc(a4*sizeof(int));
	}
	
	for(i=0;i<a1;i++)
    {
		for(j=0;j<a4;j++)
        {
			C[i][j] = 0;
		}
	}

  
    omp_set_num_threads(num_threads);  
   	begin = now();                        //start calculating the time
	if(a2 != a3)
    {         
		printf("Cannot Multipy matrices \n");     //checking for order of matrices to tbe multiplied.
 	}
	else 
    {
   
		for(i=0;i<a1;i++)
        {
            //Parallel execution starts
     #pragma omp parallel for shared(a1,a4,a2) private(j,k) firstprivate(A,B)
			for(k=0;k<a4;k++)
            {
				for(j=0;j<a2;j++)
                {
					C[i][k]+= A[i][j]*B[j][k];
				}
			}
		}
	}
	
    end = now();
    time_spent = tdiff(begin, end);
    printf("Total time in execution is %.8f sec\n", time_spent);
return 0;
}
