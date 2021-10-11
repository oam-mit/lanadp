#include <stdio.h>
#include <cuda.h>
#include "dcmtk/dcmimgle/dcmimage.h"


__constant__ int powers[9]={0,1,2,4,8,16,32,64,128};

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
__global__ void LANADP (int *arr, unsigned char *out,unsigned int rows,unsigned int cols) {
    
    
    int global_row = threadIdx.y + blockDim.y * blockIdx.y;
    int global_col = threadIdx.x + blockDim.x * blockIdx.x;

    if(global_row <rows && global_col<cols){
         if(global_row == 0 ||global_col ==0 || global_col==cols-1 || global_row==rows-1 ){
            //    if(global_row==0 && global_col==0){
            //         for(int i=0;i<3;i++){
            //             for(int i1=0;i1<3;i1++){
            //                 printf("%d  ",arr[i*cols+i1]);
            //             }

            //             printf("\n");
            //         }
            //    }
        }

          else {
              
              int mapping[9];
              mapping[1]=arr[getGlobalID(global_row,global_col+1,cols)];
              mapping[2]=arr[getGlobalID(global_row-1,global_col+1,cols)];
              mapping[3]=arr[getGlobalID(global_row-1,global_col,cols)];
              mapping[4]=arr[getGlobalID(global_row-1,global_col-1,cols)];
              mapping[5]=arr[getGlobalID(global_row,global_col-1,cols)];
              mapping[6]=arr[getGlobalID(global_row+1,global_col-1,cols)];
              mapping[7]=arr[getGlobalID(global_row+1,global_col,cols)];
              mapping[8]=arr[getGlobalID(global_row+1,global_col+1,cols)];


              int ans=0;
              for (int i=1;i<=8;i++) {
                  
                  float avg1 = ((float)(mapping[getCircularIndex(i+1)] + mapping[getCircularIndex(i+2)]))/2.0;
                  float avg2 = ((float)(mapping[getCircularIndex(i-1)] + mapping[getCircularIndex(i-2)]))/2.0;

                  //printf("Number: %d, Avg1: %f, Avg2: %f\n",mapping[i],avg1,avg2);

                  if (avg1>=(float)mapping[i] && avg2>=(float)mapping[i] || (avg1<=(float)mapping[i] && avg2<=(float)mapping[i]) ){
                      
                      //printf("Adding:%d\n",powers[i]);
                      ans+=powers[i];
                  }
              }

              

               out[(global_row-1)*(cols-2)+(global_col-1)] = ans;


              
            
                
        }
    
    }

   
  
  
}


int main() {

    int row,column;
    
    DicomImage *img=new DicomImage("512x512CT.dcm");
    

    if(img != NULL && img->getStatus()==EIS_Normal) {

        if(img->isMonochrome()) {

            img->setMinMaxWindow();
			
			row=img->getHeight();
			column=img->getWidth();

           

            int *arr = new int [row*column];

            printf("Rows: %d\tColumns:%d\n",row,column);

            img->getOutputData(arr,row*column*sizeof(int));

            unsigned char output[row-2][column-2];

            int output_rows=row-2,output_colums=column-2;

            int size_arr = row*column*sizeof(int);
            int size_out = output_rows*output_colums*sizeof(unsigned char);

            //memory variables
            int *d_arr; 
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


            LANADP<<<numberOfBlocks,numberOfThreads>>>(d_arr,d_output,row,column);
        
            cudaEventRecord(stop, 0);
            cudaEventSynchronize (stop);

            cudaEventElapsedTime(&elapsed, start, stop);
        
            printf("Time Elapsed: %f\n",elapsed);

            cudaMemcpy(output,d_output,size_out,cudaMemcpyDeviceToHost);
              

            // for(int i=0;i<row-2;i++) {
                
            //     for(int i1=0;i1<column-2;i1++) {
                    
            //         printf("%d  ",output[i][i1]);
            //     }
            //     printf("\n");
            // }


        }
    }

    
 
   

   

  



}