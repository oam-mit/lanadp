
#include <stdio.h>
#include <cuda.h>

__constant__ int powers[3][3]={
    {8,4,2},
    {16,0,1},
    {32,64,128}
};



//dichom, dcmtk
//arr to unsigned short int
__global__ void TLBAP (short int *arr, unsigned char *out,unsigned int rows,unsigned int cols,float threshold_1) {
    
    
    int global_row = threadIdx.y + blockDim.y * blockIdx.y;
    int global_col = threadIdx.x + blockDim.x * blockIdx.x;

    if(global_row <rows && global_col<cols){
         if(global_row == 0 ||global_col ==0 || global_col==cols-1 || global_row==rows-1 ){
        //do nothing
        }

          else {
              
              short int max = arr[(global_row-1)*cols+global_col];
              
              for(int i=global_row-1;i<=global_row+1;i++) {
                  for(int i1=global_col-1;i1<=global_col+1;i1++){
                      if(i!=global_row || i1!=global_col){
                          if(arr[i*cols+i1]>max){
                              max=arr[i*cols+i1];
                          }
                      }
                  }
              }
            
              float threshold_2 = threshold_1 * (float)max;

              short int ele = arr[global_row*cols + global_col];
              //printf("row:%d\tcol:%d\tThres:%f\nMax:%d\n",global_row,global_col,threshold_2,max);

              int ans=0;
              int power_i=0;

              for(int i=global_row-1;i<=global_row+1;i++){
                  int power_i1=0;
                  for(int i1=global_col-1; i1<=global_col+1 ;i1++) {
                      if(i!=global_row || i1!=global_col) {
                          if(arr[i*cols+i1]>ele && (float)arr[i*cols+i1]>threshold_2) {
                              ans+=powers[power_i][power_i1];
                          }
                      }
                      power_i1++;
                  }
                  power_i++;
              }

                out[(global_row-1)*(cols-2)+(global_col-1)] = ans;
                
        }
    }

   
  
  
}


int main() {

    int row=5,column=5;
    short int arr[row][column]={
      {3,5,4,1,1},
      {10,5,3,1,1},
      {12,2,6,1,1},
      {1,2,3,4,5},
      {1,2,3,4,5} 
    };

    unsigned char output[row-2][column-2];
 
    int output_rows=row-2,output_colums=column-2;

    int size_arr = row*column*sizeof(short int);
    int size_out = output_rows*output_colums*sizeof(unsigned char);

    //memory variables
    short int *d_arr; 
    unsigned char *d_output;

    //allocating memory
    cudaMalloc((void **)&d_arr,size_arr);
    cudaMalloc((void **)&d_output,size_out);
 
    //copying data
    cudaMemcpy(d_arr,arr,size_arr,cudaMemcpyHostToDevice); 


    dim3 numberOfBlocks(ceil((float)row/16.0),ceil((float)column/16.0),1);
    dim3 numberOfThreads(16,16,1);
 
    float elapsed=0;
    cudaEvent_t start, stop;

    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    cudaEventRecord(start, 0);


    TLBAP<<<numberOfBlocks,numberOfThreads>>>(d_arr,d_output,row,column,0.9);
 
    cudaEventRecord(stop, 0);
    cudaEventSynchronize (stop);

    cudaEventElapsedTime(&elapsed, start, stop);
 
    printf("Time Elapsed: %f\n",elapsed);

    cudaMemcpy(output,d_output,size_out,cudaMemcpyDeviceToHost);

    for(int i=0;i<row-2;i++) {
        
        for(int i1=0;i1<column-2;i1++) {
            
            printf("%d  ",output[i][i1]);
        }
        printf("\n");
    }




}