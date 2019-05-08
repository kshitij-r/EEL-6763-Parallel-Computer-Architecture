#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include<time.h>
#include<math.h>

#define DEBUG 2

/*Function Prototypes */
void initialize_data(int *,int *,int );
void scatter_data(int *, int);
void mask_operation(int *, int, int *);
void gather_results(int *, int);

/*Global Variables */
int rank, size;
int master = 0 ;
int rows_available;

int main(int argc, char *argv[])
{
    int *A; 
    int *Ap;
	int *B;
    int i,j,k;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);  
    MPI_Comm_size (MPI_COMM_WORLD, &size);  
    
    if (rank ==0) {
        printf("Change DEBUG level to see matrix values\n");
    }
    int N = atoi(argv[1]);
	int M= N+4;
    int L = N+2;

    A = (int *)malloc(N*N*sizeof(int));
    Ap = (int *)malloc(L*L*sizeof(int));
	B= (int *)malloc(M*M*sizeof(int));
    initialize_data(A,B, N);
    scatter_data(A, N);
    mask_operation(A, N, Ap);
    gather_results(Ap, N);
    MPI_Finalize();
    free(A);
    free(Ap);
    return 0; 
} 

void initialize_data(int *A,int *B, int N)
{
    int i, j, k = 0;
    int M= N+4;
 int L = N+2;
    srand((unsigned)time(NULL));
    if(rank == 0){
        for(i = 0; i < N; i++){
            for(j = 0; j < N; j++){
                A[j+N*i] = i+j;
            }
        }
		for (i=0; i< 2*M+2;i++){
			B[i]= 0;
		}
		for(i=0;i<N;i++){
			for(j=0;j<N;j++){
				B[(j+2)+((i+2)*M)]= A[j+N*i];
			}
		}
		for(i=((N+3)+((N+3)*M));i<M;i++){
			B[i]=0;
		}
  }

    
    if(DEBUG > 0 && rank == 0){
        printf("\n(%d)Initial: \n",rank);
        for(i = 0; i < N; i++){
            for(j = 0; j < N; j++){
                printf ("%d ",A[j + (N * i)]);
            }
            printf("\n");
        }
	
    }
}

// Send data to even workers first and then to odd workers .
void scatter_data(int *B, int N)
{
    int M= N+4;
   int L = N+2;
    int i , j,k;
    // send count and displacement array for EVEN workers
    int  sendcnts_even[size], displs_even[size];
    
    // send count and displacement array for ODD workers
    int sendcnts_odd[size],displs_odd[size];
    
    // temporary variables for calculating the send count and displacement array
    int sendcnts_present, displs_present, displs_previous = 0, sendcnts_previous = 0;
    
    // minimum number of rows sent to any of the worker
    int min_row = L / (size);

    for(i=0 ;i<size ; i++){
        
        if(L % size >= i+1 ){  // adjusting for extra rows
            sendcnts_present = (min_row +1)*L;
        }
        else{
            sendcnts_present = min_row*L;
        }
        
        // finding displacement wrt to the start of matrix buffer
        displs_present = displs_previous + sendcnts_previous ;
        
        // Adjusting the sendcount and displacement according to the overlap.
        if(i%2 == 0){
            //printf("EVEN %d",i);
            if( i == size-1 || i == 0 )
                sendcnts_even[i] = sendcnts_present + L;
            else
                sendcnts_even[i] = sendcnts_present + 2*L;
            
            displs_even[i] = displs_present - L ;
            
            /*
             MPI standard allows send count value to be zero but the
             Mpich and OpenMPI do not allow it so sending single element
            */
            sendcnts_odd[i] = 1;
            
            displs_odd[i] = displs_present+ L ;
            
        }else{
            //printf("ODD %d",i);
            if( i == size-1 || i == 0)
                sendcnts_odd[i] = sendcnts_present + L;
            else
                sendcnts_odd[i] = sendcnts_present+ 2*L;
            
            displs_odd[i] = displs_present - L ; 
            
            sendcnts_even[i] = 1;
            displs_even[i] = displs_present+L ;
        }
        
        //settin the intial displacement to zero
        if(i == 0){
            displs_even[i] = 0;
            displs_odd [i] = 0;
        }
        
        displs_previous =displs_present ;
        sendcnts_previous = sendcnts_present;
        
        if(DEBUG > 1 ){
            if(rank== master){
                printf("\nNo adjustment rank = %d , displacement = %d , count = %d , end = %d \n" , i , displs_present, sendcnts_present , displs_present+sendcnts_present);
            }
            if(rank== master){
                printf("Adjusted Even rank = %d , displacement = %d , count= %d ,end = %d \n" , i , displs_even[i], sendcnts_even[i], displs_even[i]+sendcnts_even[i]);
            }
            
            if(rank== master){
                printf("Adjusted Odd  rank = %d , displacement = %d , count = %d ,end = %d \n" , i , displs_odd[i], sendcnts_odd[i], displs_odd[i]+sendcnts_odd[i]);
            }
        }
        
    }
    
    //Maximum number of elements received by any of the worker node.
    int max_recv_count =(min_row+3)*L;
    
    // Using the A buffer as receive buffer
    MPI_Scatterv(B,sendcnts_even,displs_even ,MPI_INT,B, max_recv_count, MPI_INT,master,MPI_COMM_WORLD);
    
    //storing the first element of receive buffer
    int temp = B[0];
    
    MPI_Scatterv(B,sendcnts_odd,displs_odd ,MPI_INT,B, max_recv_count, MPI_INT,master,MPI_COMM_WORLD);
    
    //adjusting for the one element which overwrote in the receive buffer
    // Read above why can't we send zero count
    if(rank%2 == 0){
        B[0] = temp;
    }
    
    if(DEBUG > 1 ){
        
        for(k= 0 ; k < size ; k++){
            if(k == rank){
                printf("\n(%d)Scattered Data :\n",rank);
                for(i = 0; i < min_row+3; i++){
                    for(j = 0; j < M; j++){
                        printf ("%d ",B[j + (N * i)]);
                    }
                    printf("\n");
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
}


void mask_operation(int * B, int N, int *Ap)
{
    int i,j,k,l;
    int M= N+4;
 int L = N+2;
    int *magx;
    int *magy;
    
    magx = (int *)malloc(l*L*sizeof(int));
    magy = (int *)malloc(L*L*sizeof(int));
    
    int ker_x[3][3] = { {1,0,-1},{2,0,-2},{1,0,-1}};
		int ker_y[3][3] = { {-1,-2,-1},{0,0,0},{1,2,1}};
		
		int iker_x[3][3] = {0};
		int iker_y[3][3] = {0};
		
		for(i=0;i<3;i++){
		for(j=0;j<3;j++){
			iker_x[i][j] = ker_x[2-i][2-j];
			iker_y[i][j] = ker_y[2-i][2-j];
			
		}
	}
    /*
     rows_available is calculated which is nothing but the number of rows available with the worker.
     */
    int rows_available = L / (size);
    if(L % size >= rank+1 ){
        // adjusting for extra rows
        rows_available++;
    }
    //adjusting for overlap of rows
    if( rank == size-1 || rank == 0 )
        rows_available++;
    else
        rows_available= rows_available+2;
    
    
    
    //Performing the mask operation
    for( i = 0; i < rows_available; i++)
    {
        for(i=0;i<L;i++){
		for(j=0;j<L;j++){
     // printf("magx[%d][%d] perfomed by thread %d \n",i,j,omp_get_thread_num());
			for(k=0;k<3;k++){
				for(l=0;l<3;l++){
					magx[i*(N+2)+j]+= (B[(i+k)*M+(j+l)]*iker_x[k][l]);
          magy[i*(N+2)+j]+= (B[(i+k)*M+(j+l)]*iker_y[k][l]);
          Ap[i*(N+2)+j] = sqrt(powf(magx[i*(N+2)+j],2)+powf(magy[i*(N+2)+j],2));
				}
			}
		}
	}
        
    }
    
    if(DEBUG > 1 ){ //check if the matrix is initialized properly or not
        MPI_Barrier(MPI_COMM_WORLD);
        for(k= 0 ; k < size ; k++){
            if(k == rank){
                printf("\n--------------------(%d)Processed Data------------------\n",rank);
                for(i = 0; i < rows_available; i++){
                    for(j = 0; j < L; j++){
                        printf ("%d ",Ap[j + (L * i)]);
                    }
                    printf("\n");
                }
            }
            MPI_Barrier(MPI_COMM_WORLD);
        }
    }
}



void gather_results(int *Ap, int N) {
    int i,j,k;
    int M= N+4;
 int L = N+2;
    int send_row_count = M / (size);
    int * sendBuf;
    if(N % size >= rank+1 ){
        // adjusting for extra rows
        send_row_count++;
    }
    
    
    int sendcnts[size];
    int disp[size];
    
    if (rank == master) {
        
        for( i = 0; i < size; i++)
        {
            int row_count = L / (size);
            if(L % size >= i+1 ){
                // adjusting for extra rows
                row_count++;
            }
            
            sendcnts[i] = row_count*M;
            if (i==0) {
                disp[i] = 0;
            }else{
                disp[i] = sendcnts[i-1]+disp[i-1];
            }
        }
    }
    
    if(rank == 0 ){
        sendBuf = &Ap[0];
    }else{
        sendBuf = &Ap[M];
    }
    
    MPI_Gatherv(sendBuf, send_row_count*L, MPI_INT, Ap, sendcnts, disp, MPI_INT, master, MPI_COMM_WORLD);
    
    if(DEBUG > 0 ){ 
        MPI_Barrier(MPI_COMM_WORLD);
        if(rank == master){
            printf("\n(%d)RESULT :\n",rank);
            for(i = 0; i < L; i++){
                for(j = 0; j < L; j++){
                    printf ("%d ",Ap[j + (N+2 * i)]);
                }
                printf("\n");
            }
        }
    }
    
}







