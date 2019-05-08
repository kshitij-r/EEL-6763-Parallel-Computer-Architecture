//Kshitij Raj
#include<mpi.h>
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>

	enum
	{
		MASTER=0,
		FROM_MASTER=1,
		FROM_WORKER=2
	};

double estimate_g(double lower_bound, double upper_bound, long long int N);
void collect_results(double *addr);
double initialize_data();
double result;
double SEED;


	int main(int argc,char* argv[]){
		
		int num_tasks; //No. of nodes//
    int i;
    
		MPI_Init(&argc, &argv);
     double lower_bound = atof(argv[1]); /* a */
     printf("a= %lf \n", lower_bound);
		double upper_bound = atof(argv[2]); /* b */
     printf("b= %lf \n", upper_bound);
		long long int N = atof(argv[3]); /* n */
   printf("N= %lli \n", N);
   
	
		initialize_data();
    result= estimate_g(lower_bound, upper_bound, N);
		collect_results(&result);
		
		MPI_Finalize();
		return 0;
	}

double estimate_g(double a, double b, long long int N){

	static	int         num_tasks;     /* number of tasks     */
	      	int         task_id;       /* number of processes  */
  static	int         num_workers;   /* number of worker tasks */
		      int         source;        /* rank of sender       */
          int         dest;          /* rank of receiver     */
		      int         mtype;         /* message type */
		      int         tag = 0;       /* tag for messages     */
		      char        message[100];  /* storage for message  */
		      MPI_Status  status;        /* return status for receive */
		      int         samples;
		      int         offset; 
		
   //printf("entered estimate_g \n");
   
		MPI_Comm_rank(MPI_COMM_WORLD, &task_id);
		MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
		
		num_workers =num_tasks-1;
 
   
		if(task_id == MASTER) {
			long long int average_samples = N/ num_workers;
		printf("No.of average samples = %lli \n", average_samples);
			long long int extra_samples = N % num_workers;
		printf("No.of extra samples = %lli \n", extra_samples);
	  
	   
			// Send data to workers
			int mtype = 1;
			long long int offset = 0;
			long long int samples = 0;
			
			for (dest = 1; dest < num_tasks; dest++) {
				
				samples = (dest <= extra_samples) ? average_samples+1 : average_samples; //Deciding number of samples to send
				printf( " Sending %d samples to worker %d \n ",samples, dest);
				
				MPI_Send(&offset, 1, MPI_LONG_LONG, dest, mtype, MPI_COMM_WORLD);
				MPI_Send(&samples, 1 , MPI_LONG_LONG, dest, mtype, MPI_COMM_WORLD);
				
				offset+= samples;  
			}
			printf("Sent all the sample numbers \n");
		}
		else {
			mtype = FROM_MASTER;
			
			MPI_Recv(&offset, 1, MPI_LONG_LONG, MASTER, mtype, MPI_COMM_WORLD, &status);
			MPI_Recv(&samples, 1, MPI_LONG_LONG, MASTER, mtype, MPI_COMM_WORLD, &status);
			
			//Compute f(xi) value
				double x =0, num, den, fx[samples];
				long long int i;
     
         srand(SEED);
			for(i =0; i<samples ; i++){  
				x= ((double)rand()/(double)RAND_MAX )*(b-a)+a;
				num = 8* (sqrt(2*M_PI));
				den = pow((M_E),(pow((2*x),2)));
        fx[i] = ((b-a)/N)*(num/den);
			  result = result + fx[i];
      }
	 
}
 
	return result;
}
	
void collect_results(double *addr) {
		int         num_tasks;     /* number of tasks     */
		int         task_id;       /* number of processes  */
		int         num_workers;   /* number of worker tasks */
		int         source;        /* rank of sender       */
		int         dest;          /* rank of receiver     */
		int         mtype;         /* message type */
		int         tag = 0;       /* tag for messages     */
		char        message[100];  /* storage for message  */
		MPI_Status  status;        /* return status for receive */
		int         samples;
		int         offset; 

    double h;
    
   // printf("entered collect_results \n");
    MPI_Comm_rank(MPI_COMM_WORLD, &task_id);
  	MPI_Comm_size(MPI_COMM_WORLD, &num_tasks);
   int i;
		if(task_id == MASTER){
      int j;
      mtype = FROM_WORKER;
        
			for (i=1; i<num_tasks; i++){
				source = i;
				MPI_Recv(&offset, 1, MPI_LONG_LONG, source, mtype, MPI_COMM_WORLD, &status);
		
        MPI_Recv(&samples, 1, MPI_LONG_LONG, source, mtype, MPI_COMM_WORLD, &status);
		
        MPI_Recv(&result, 1, MPI_DOUBLE, source, mtype, MPI_COMM_WORLD, &status);
						h+= result;		
			}	
			printf("estimate_g = %lf \n", h);
		}		
	else {
			mtype =FROM_WORKER;
				MPI_Send( &offset, 1, MPI_LONG_LONG, MASTER, mtype, MPI_COMM_WORLD);
				MPI_Send(&samples, 1, MPI_LONG_LONG, MASTER, mtype, MPI_COMM_WORLD);
				MPI_Send(&result, 1, MPI_DOUBLE, MASTER, mtype, MPI_COMM_WORLD);				
		}	
}

double initialize_data(){
    
    int num_taks, task_id;
    MPI_Comm_rank(MPI_COMM_WORLD, &task_id);
    SEED = (double)rand()*task_id;  
    return SEED;
}
