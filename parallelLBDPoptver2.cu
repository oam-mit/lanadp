#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include<iostream>
#include<stdlib.h>
#include "dcmtk/dcmimgle/dcmimage.h"
#include<algorithm>
#include<stdio.h>


#define TILE_WIDTH 16
#define TILE_HEIGHT 16

typedef  unsigned short int USint;

using namespace std;

__constant__ USint multiple2[16];

__constant__ char x[8];

__constant__ char y[8];


__global__ void getLBDP(USint *pdata,int h, int w, unsigned char *out) {

	int tx=blockIdx.x*blockDim.x+threadIdx.x;
	int ty=blockIdx.y*blockDim.y+threadIdx.y;

	//printf("inside ty=%d tx=%d  w=%d  h=%d",ty,tx,w,h);	
	//printf("x[0]=%d",x[0]);
	

	__shared__ USint bpdata[TILE_HEIGHT+2][TILE_WIDTH+2];

	bpdata[threadIdx.y+1][threadIdx.x+1]=pdata[ty*w+tx];

	//printf("YESAYESUYSE");
	USint bitp[16];


	//Since notion of threadIdx.x is changed to threadIdx.y
	if(threadIdx.x==0)
	{
		if(tx != 0)
		{
			//printf("A");
			bpdata[threadIdx.y+1][0]=pdata[ty*w+tx-1];
			if(threadIdx.y==0 && ty!=0)
				bpdata[0][0]=pdata[(ty-1)*w+tx-1];
		}
		
	}
	if(threadIdx.x==TILE_WIDTH-1)
	{
		//printf("B");
		if(tx != w-1) {
			bpdata[threadIdx.y+1][threadIdx.x+2]=pdata[ty*w+tx+1];
			if(threadIdx.y==TILE_HEIGHT-1 && ty!=h-1)
				bpdata[TILE_HEIGHT+1][TILE_WIDTH+1]=pdata[(ty+1)*w+tx+1];
		}
	}

	if(threadIdx.y==0)
	{
		if(ty != 0)
		{
			//printf("C");
			bpdata[0][threadIdx.x+1]=pdata[(ty-1)*w+tx];
			if(threadIdx.x == TILE_WIDTH-1 && tx!=w-1)
				bpdata[0][TILE_WIDTH+1]=pdata[(ty-1)*w+tx+1];
		}
	}
	if(threadIdx.y==TILE_HEIGHT-1)
	{
		if(ty!=h-1)
		{
			//printf("D");
			bpdata[threadIdx.y+2][threadIdx.x+1]=pdata[(ty+1)*w+tx];
			if(threadIdx.x==0 && tx!=0)
				bpdata[TILE_HEIGHT+1][0]=pdata[(ty+1)*w+tx-1];
		}
	} 


			
	/*if(threadIdx.y==0 && ty !=0) {
		bpdata[0][threadIdx.x]=pdata[(ty-1)*w+tx];
		if(threadIdx.x == TILE_WIDTH-1 && tx!=w-1)
			bpdata[0][TILE_WIDTH+1]=pdata[(ty-1)*w+tx+1];
	}	
	if(threadIdx.y==TILE_HEIGHT-1 && ty != h-1)
	{
		bpdata[TILE_HEIGHT+1][threadIdx.x+1]=pdata[(ty+1)*w+tx];
		if(threadIdx.x == TILE_WIDTH-1 && tx !=0)
			bpdata[TILE_HEIGHT+1][0]=pdata[(ty+1)*w+tx-1];
	}		
	if(threadIdx.x==0 && tx !=0)
	{
		bpdata[threadIdx.y+1][0]=pdata[ty*w+tx-1];
		if(threadIdx.y==0 && ty!=0)
			bpdata[0][0]=pdata[(ty-1)*w+tx-1];
	}
	if(threadIdx.x==TILE_WIDTH-1 && tx != w-1)
	{
		bpdata[threadIdx.y+1][TILE_WIDTH+1]=pdata[ty*w+tx+1];
		if(threadIdx.y == TILE_HEIGHT-1 && ty!=h-1)
			bpdata[TILE_HEIGHT+1][TILE_WIDTH+1]=pdata[(ty+1)*w+tx+1];
	}*/


	if(ty==0 || ty==h-1 || tx==0 || tx==w-1)
		return;

	__syncthreads();
	
	/*int m,n;
	if(ty==1 && tx==15)
	{
		for(m=0;m<16+2;m++)
		{
			for(n=0;n<16+2;n++)
				printf("%d   ",bpdata[m][n]);
			printf("\n");
		}
	}*/
		
	//char *out;
	//out=(char *)calloc((h-1)*(w-1),sizeof(char));

	//printf("LKJLKJLK");
			
	USint centre;
	centre=bpdata[threadIdx.y+1][threadIdx.x+1];
	//printf("CENTER=%d",centre);
	USint cur;
	USint pow;
	USint btot;
	//printf("BBBBBB");		
	btot=0;

	//printf("before LOOP1");
	for(int l=0;l<16;l++)
		bitp[l]=0;
        for(int l=0;l<16;l++)
	{
		pow=1;

		//printf("cur = %d\n",cur);
		for(int k=0;k<8;k++)
		{
			cur=bpdata[threadIdx.y+1+y[k]][threadIdx.x+1+x[k]];
			//printf("cur = %d bitp[%d]=%d\n",cur,l,bitp[l]);
			if(cur & multiple2[l])
				bitp[l] = bitp[l] + pow;

			//printf("%d \n",pow);	
			pow = pow*2;
			//printf("pow =%d\n",pow);
			//if(ty==1 && tx==15 && l==0)
			//	printf("%d ",cur);
		}
	}
	//printf("before LOOP2");
	for(int l=15;l>=0;l--)
	{	
		//printf("bitp[%d]=%d\n",l,bitp[l]);
		if(bitp[l] > centre)
			btot = btot+multiple2[15-l];
	}
	//printf("xxxxx centre=%d btot=%d",centre,btot);
	out[(ty-1)*(w-2)+(tx-1)]=(unsigned char)(btot>>8);

	//if(ty==1 && tx==15)
	//	printf("%d ",out[(ty-1)*(w-2)+(tx-1)]);
}

					


int main(int argc, char *argv[]) {
	
	cudaEvent_t pstart,pstop;

	float elapsedTime;
	cudaEventCreate(&pstart);
	cudaEventCreate(&pstop);

	cudaEventRecord(pstart,0);

	DicomImage *img=new DicomImage("../images/512x512CT.dcm");

	if(img != NULL && img->getStatus()==EIS_Normal) {
		if(img->isMonochrome()) {
			img->setMinMaxWindow();
			
			int h=img->getHeight();
			int w=img->getWidth();

			/*int h=3;
                        int w=3;

                        USint pdata[]={10,99,70,1,25,82,55,71,35}; 
			correct output = 238

				*/

			//cout<<"size = "<<img->getOutputDataSize(8)<<" "<<sizeof(Uint)<<" ";
			USint *pdata = (USint *)new USint[h*w];

			//cout<<"  "<<sizeof(pdata)<<" ";
			//(Uint *)img->getOutputData(8);

			//cout<<"Depth = "<<img->getDepth();

			img->getOutputData(pdata,w*h*sizeof(USint));


			USint *d_pdata;
			cudaMalloc((void **)&d_pdata,w*h*sizeof(USint));

			cudaMemcpy(d_pdata,pdata,w*h*sizeof(USint),cudaMemcpyHostToDevice);

			//double min;
			//double max;
			//img->getMinMaxValues(min,max);

			//cout<<"Min = "<<min<<" max = "<<max;

			//DiPixel *indata = (DiPixel *)img->getInterData();
			
			//for(int i=0;i<8;i++)
			//	cout<<pdata[i]<<"  ";
			
			/*if(pdata != NULL) {
				
				//cout<<w<<"   "<<h;

				for(int i=0;i<h;i++) {
					for(int j=0;j<w;j++) {
						if(pdata[i*w+j])
							cout<<pdata[i*w+j]<<" ";
					}
					//cout<<endl;
				}
			*/
			/*ofstream ofs;
			ofs.open("MR.pgm");
			ofs<<"P2"<<endl<<w<<endl<<h<<endl;
			ofs<<65535<<endl;

			if(pdata != NULL) {
				
				//cout<<w<<"   "<<h;

				for(int i=0;i<h;i++) {
					for(int j=0;j</*int h=3;
                        int w=3;*/

                        //USint pdata[]={15,0,56,0,50,0,160,0,100};w;j++) {
			/*			ofs<<pdata[i*w+j]<<endl;
					}
				}
			}*/


                     	unsigned char *out;

			out=(unsigned char *)calloc((h-2)*(w-2),sizeof(unsigned char));
			
			unsigned char *d_out;
			cudaMalloc((void **)&d_out,((h-2)*(w-2)*sizeof(unsigned char)));

			USint *h_multiple2;
			h_multiple2=(USint *)calloc(sizeof(USint)*8,sizeof(USint));

			h_multiple2[0]=1;
			for(int l=1;l<16;l++) 
			{
				h_multiple2[l]=h_multiple2[l-1]<<1;
				//printf("%d   \n",multiple2[l]);		
			}
			char h_x[8]={1,1,0,-1,-1,-1,0,1};
			char h_y[8]={0,-1,-1,-1,0,1,1,1};

			cudaMemcpyToSymbol(multiple2,h_multiple2,sizeof(USint)*8*sizeof(USint),0,cudaMemcpyHostToDevice);

			cudaMemcpyToSymbol(x,h_x,sizeof(char)*8,0,cudaMemcpyHostToDevice);
			cudaMemcpyToSymbol(y,h_y,sizeof(char)*8,0,cudaMemcpyHostToDevice);


			cudaEvent_t kstart, kstop;

			dim3 blk(std::max(w/TILE_WIDTH,1),std::max(h/TILE_HEIGHT,1));
			dim3 thrd(std::min(w,TILE_WIDTH),std::min(h,TILE_HEIGHT));
			//dim3 thrd(3,3);

			cudaEventCreate(&kstart);
			cudaEventCreate(&kstop); 

			cudaEventRecord(kstart, 0);
			getLBDP<<<blk,thrd>>>(d_pdata,h,w,d_out);
			
			cudaEventRecord(kstop, 0);
			cudaEventSynchronize(kstop);
			float elapsedTime;
			cudaEventElapsedTime(&elapsedTime, kstart, kstop);			
			
			printf("Kernel Elapsed time =%f\n",elapsedTime);

			cudaMemcpy(out,d_out,((h-2)*(w-2)*sizeof(unsigned char)),cudaMemcpyDeviceToHost);

			/*
			for(int i=1;i<h-1;i++)
			{
				for(int j=1;j<w-2;j++)
				{
					/*	
					unsigned char  ndig=0;
					unsigned char pdig=1;
					unsigned char btot=0;
					for(int k=-1;k<2;k++)
					{
						for(int l=-1;l<2;l++)
						{
								if(!(k|l))
									continue;
								if(pdata[i*w+j]>pdata[(i+k)*w+(j+l)])
									btot+=pdig;
								ndig++;
								pdig*=2;
						}
					}
					if(btot !=out[(i-1)*(w-1)+(j-1)])
						printf("Error i=%d,j=%d btot=%d out=%d\n",i,j,btot,out[(i-1)*(w-1)+(j-1)]);
					*/
			/*		cout<<"h="<<i<<"  w="<<j<<"  "<<((USint)out[(i-1)*(w-2)+(j-1)])<<endl;
				}
				cout<<endl;
					;
			}*/

			//cout<<sizeof(pdata);
			//cout<<" "<<sizeof(indata);
			//ofs.close();
			delete out;
			//delete pdata;
			//delete indata;
		}
	}

	cudaEventRecord(pstop, 0);
	cudaEventSynchronize(pstop);
	cudaEventElapsedTime(&elapsedTime, pstart, pstop);
	printf("Program Elapsed time =%f\n",elapsedTime);

	return 0;
}	
