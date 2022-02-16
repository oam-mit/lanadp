#include <stdio.h>
#include <cuda.h>
#include "dcmtk/dcmimgle/dcmimage.h"

__constant__ int powers[3][3]={
    {8,4,2},
    {16,0,1},
    {32,64,128}
};


#define NUMBER_OF_THREADS 16

__device__ int getGlobalID(int row,int col,int cols){
    return row*cols+col;
}

__device__ int getCircularIndex(int i){
    
    if (i<1) {
        return 8 + i;
    }
    else if (i>8 )
      return (i%9 )+1;
    else 
      return i;
  
}

//dichom, dcmtk
__global__ void TLBAP (unsigned int *arr, unsigned char *out,unsigned int rows,unsigned int cols,float threshold_1) {
    
    
    int global_row = threadIdx.y + blockDim.y * blockIdx.y;
    int global_col = threadIdx.x + blockDim.x * blockIdx.x;

    global_row++;
    global_col++;


    int shared_no_of_rows = NUMBER_OF_THREADS+2;
    int shared_no_of_cols = NUMBER_OF_THREADS+2;

    __shared__ short int shared_data[NUMBER_OF_THREADS+2][NUMBER_OF_THREADS+2];

    

    int local_row = threadIdx.y;
    int local_col = threadIdx.x;

    int shared_row = local_row+1;
    int shared_col = local_col+1;

    //printf("global_row: %d\nglobal_col:%d\nlocal_row:%d\nlocal_col:%d\nsahred_row: %d\nshared_col:%d\n",global_row,global_col,local_row,local_col,shared_row,shared_col);

    shared_data[shared_row][shared_col]=arr[getGlobalID(global_row,global_col,cols)];

    if (local_row==0){
        //if first_row, load data above itself compusorily
        shared_data[shared_row-1][shared_col]=arr[getGlobalID(global_row-1,global_col,cols)];

        //if first column, load diagonally above
        if(local_col==0){
            shared_data[shared_row-1][shared_col-1] = arr[getGlobalID(global_row-1,global_col-1,cols)];
        }

       //if last column load digonally above 
       if(local_col == gridDim.x-1){
            shared_data[shared_row-1][shared_col+1] = arr[getGlobalID(global_row-1,global_col+1,cols)];
        }
    }

    if (local_row==gridDim.y-1) {

        //if last row load elements below it
        shared_data[shared_row+1][shared_col] = arr[getGlobalID(global_row+1,global_col,cols)];

        //if first column, load diagonally below 
        if(local_col==0){
            shared_data[shared_row+1][shared_col-1] = arr[getGlobalID(global_row+1,global_col-1,cols)];
        }

        //if last column, load diagonally below
        if(local_col == gridDim.x-1){
            shared_data[shared_row+1][shared_col+1] = arr[getGlobalID(global_row+1,global_col+1,cols)];
        }
    }

    if (local_col==0){
        
        //if first column load elements to left
        shared_data[shared_row][shared_col-1] = arr[getGlobalID(global_row,global_col-1,cols)];
    
    }

    if(local_col==gridDim.x-1) {
        
        //if last column load elements to right
        shared_data[shared_row][shared_col+1] = arr[getGlobalID(global_row,global_col+1,cols)];
    }
    global_row--;
    global_col--;

    __syncthreads();

  /*for (int i=0;i<3;i++){
      for(int i1=0;i1<3;i1++){
          printf("%d  ",shared_data[i][i1]);
      }
      printf("\n");
  }*/

  if(global_row <rows && global_col<cols){
    if(global_row == 0 ||global_col ==0 || global_col==cols-1 || global_row==rows-1 ){
   //do nothing
   }

     else {
         
         short int max = shared_data[(shared_row-1)][shared_col];
         
         for(int i=shared_row-1;i<=shared_row+1;i++) {
             for(int i1=shared_col-1;i1<=shared_col+1;i1++){
                 if(i!=shared_row || i1!=shared_col){
                     if(shared_data[i][i1]>max){
                         max=shared_data[i][i1];
                     }
                 }
             }
         }
       
         float threshold_2 = threshold_1 * (float)max;

         short int ele = shared_data[shared_row][shared_col];
         //printf("row:%d\tcol:%d\tThres:%f\nMax:%d\n",global_row,global_col,threshold_2,max);

         int ans=0;
         int power_i=0;

         for(int i=shared_row-1;i<=shared_row+1;i++){
             int power_i1=0;
             for(int i1=shared_col-1; i1<=shared_col+1 ;i1++) {
                 if(i!=shared_row || i1!=shared_col) {
                     if(shared_data[i][i1]>ele && (float)shared_data[i][i1]>threshold_2) {
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


int main(int argc, char *argv[]) {

    int row,column;
    printf("before readd\n");
    
    DicomImage *img=new DicomImage(argv[1]);
    

    if(img != NULL && img->getStatus()==EIS_Normal) {

        if(img->isMonochrome()) {

            img->setMinMaxWindow();
			
			row=img->getHeight();
			column=img->getWidth();

           

            unsigned int *arr = new unsigned int [row*column];

            printf("Rows: %d\tColumns:%d\n",row,column);
            printf("Image depth %d\n",img->getDepth());

            img->getOutputData(arr,row*column*sizeof(unsigned int));

            printf("after 4reading img\n");

            unsigned char output[row-2][column-2];

            int output_rows=row-2,output_colums=column-2;

            int size_arr = row*column*sizeof(unsigned int);
            int size_out = output_rows*output_colums*sizeof(unsigned char);

            //memory variables
            unsigned int *d_arr; 
            unsigned char *d_output;

            //allocating memory
            cudaMalloc((void **)&d_arr,size_arr);
            cudaMalloc((void **)&d_output,size_out);
        
            //copying data
            cudaMemcpy(d_arr,arr,size_arr,cudaMemcpyHostToDevice); 


            dim3 numberOfBlocks(ceil((float)(row-2)/16.0),ceil((float)(column-2)/16.0),1);
            dim3 numberOfThreads(NUMBER_OF_THREADS,NUMBER_OF_THREADS,1);
        
            float elapsed=0;
            cudaEvent_t start, stop;

            cudaEventCreate(&start);
            cudaEventCreate(&stop);

            cudaEventRecord(start, 0);

            for(int i=1;i<=50;i++)
                TLBAP<<<numberOfBlocks,numberOfThreads>>>(d_arr,d_output,row,column,0.9);
        
            cudaEventRecord(stop, 0);
            cudaEventSynchronize (stop);

            cudaEventElapsedTime(&elapsed, start, stop);
        
            printf("Time Elapsed: %f\n",elapsed);

            cudaMemcpy(output,d_output,size_out,cudaMemcpyDeviceToHost);
        
            
            // for(int i=0;i<row;i++)
            // {
            //     for(int i1=0;i1<column;i1++)
            //     {
            //         printf("%d ",output[i][i1]);
            //     }
            //     printf("\n");
            // }

            std::ofstream outfile;

            outfile.open("tlbap_shared.txt", std::ios_base::app); // append instead of overwrite
            outfile <<argv[1]<<": " <<elapsed<<"\n"; 



        }
    }








    




}