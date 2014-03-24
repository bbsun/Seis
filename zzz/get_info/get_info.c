#include"stdio.h"
#include"stdlib.h"
#include "mysu.h"
#include "string.h"
char info1[1024]="    CSP(SP)     NTR       NT        DT     XMINALL   XMAXALL   YMINALL   YMAXALL";
char info2[1024]="        SP       NTR     BEGIN       END        SX       SY       XMIN      XMAX      YMIN      YMAX";
long getfilesize(FILE *fp);
void getinfo(char *infile,char *outfile,char* outfile2,char* dir,char *key); 
int  main( int argc, char *argv[])
{
  if(argc != 6){
    printf(" \n  get_info -- statistics header information\n\n");
    printf("  Usage:\n\n  get_info  [infile] [infofile] [key] \n\n");
    printf("  Parameters: \n");
    printf("  infile  :  inputfile in SU format\n");
    printf("  infofile:  file to store statistics information\n");
    printf("  coodfile:  file to store the coordinates\n");
    printf("  dir     :   working directory\n");
    printf("  key=fldr/ep/cdp/sx    keyword describing the gathers\n\n");
    return 0;
  }
  printf("input file is %s\n",argv[1]);
  printf("info file  is %s\n",argv[2]);
  printf("coord file is %s\n",argv[3]);
  printf("dir        is %s\n",argv[4]);
  printf("keyword    is %s\n",argv[5]);
  getinfo(argv[1],argv[2],argv[3],argv[4],argv[5]);
  return 0;
}

void getinfo(char *infile,char *outfile,char*outfile_coord,char*dir,char *key)
{
  int nt,nseek;
  int dt,shotnumber,k,olds,traceall;
  long POS=0;
  int *shot,ikey;
  double  scalco=0;
  su tr;
  FILE* outfp;
  FILE* infp;
  FILE* outcord;
  infp=fopen(infile,"r");
  outfp=fopen(outfile,"w");
  outcord=fopen(outfile_coord,"w");
  if(! strcmp(key,"sx"))
    ikey=18;
  else if(! strcmp(key,"ep"))
    ikey=4;
  else if(! strcmp(key,"cdp"))
    ikey=5;
  else 
    ikey=2;
    //fprintf(stderr,"ikey=%d\n",ikey);
  shot=&tr.tracl+ikey;        
  k=1;shotnumber=0;
  fread(&tr,1,240,infp);//读第一道
  dt=tr.dt;//tr.dt:采样间隔
  nt=tr.ns;//nt:采样点数 tr.ns：采样点数
  olds=*shot;//olds:原先炮号
  nseek=nt*4+240;//nseek:一道的数据大小
  char buf[nseek-240];
  traceall=getfilesize(infp)/nseek;
  int ntr=1;
  int maxntr=1;
  int starttrace=1;
  int sx,sy,xmin,xmax,ymin,ymax;
  int sz,gx,gy,gz;
  int xminall,xmaxall,yminall,ymaxall;
  float sxf,szf;
  float gxf,gzf;
  // specify the sz and gz coordinate
  sz=0;
  gz=0;
  if(tr.gx<tr.sx){
    xmin=tr.gx;
    xmax=tr.sx;
  }
  else{
    xmin=tr.sx;
    xmax=tr.gx;
  }
  if(tr.gy<tr.sy){
    ymin=tr.gy;
    ymax=tr.sy;
  }
  else{
    ymin=tr.sy;
    ymax=tr.gy;
  }
  sx=tr.sx;
  sy=tr.sy;
  scalco = tr.scalco>0?tr.scalco:1.0/abs(tr.scalco); // normally this will not change 
  xminall=xmin;
  xmaxall=xmax;
  yminall=ymin;
  ymaxall=ymax;
  fprintf(outfp,"%s\n",info1);
  fprintf(outfp,"%10d%10d%10d%10d%10d%10d%10d%10d\n",0,0,0,0,0,0,0,0);
  fprintf(outfp,"%s\n",info2);
  char filename[1024];
  int ishot = 0;
  FILE * CSGF=0;
  sprintf(filename,"%sCSG%d.dat",dir,ishot++);
  CSGF=fopen(filename,"w");
  int itt=1;
  sxf = tr.sx * scalco; // set sx
  szf = 0.0f;            // set sz
  gzf = 0.0f;
  fprintf(outcord,"%12.2f %12.2f %12.2f %12.2f\n",sxf,szf,tr.gx*scalco,gzf);
  int size_trace_data = nt*4;
  fseek(infp,POS,0);//指针指向文件开始
  fread(&tr,1,240,infp);//读取240道头信息到结构体tr中
  fread(buf,1,size_trace_data,infp);//readin trace data
  fwrite(buf,1,size_trace_data,CSGF);//writeout trace data;
  
  //读su文件
  do {
    POS+=nseek;//指针指向下一道道头位置
    fseek(infp,POS,0);//指针指向文件开始
    fread(&tr,1,240,infp);//读取240道头信息到结构体tr中
    if(k<=traceall-1){
    if(*shot==olds ){
    fread(buf,1,size_trace_data,infp);//readin trace data
    fwrite(buf,1,size_trace_data,CSGF);//writeout trace data;
    }
    else
      {
	fclose(CSGF);
	char filename2[1024];
	sprintf(filename2,"%sCSG%d.dat",dir,ishot++);
	CSGF=fopen(filename2,"w");
	fread(buf,1,size_trace_data,infp);//readin trace data
	fwrite(buf,1,size_trace_data,CSGF);//writeout trace data;
      }
    }
    ntr++;//ntr：某炮的道数
    itt = *shot!=olds?1:itt;

    if(k<=traceall-1)
      if((*shot)==olds)
	fprintf(outcord,"%12.2f %12.2f %12.2f %12.2f\n",sxf,szf,tr.gx*scalco,gzf);
    if(*shot!=olds)//判断是否到了下一炮
      {
        shotnumber++;//统计总炮数
	olds=*shot;
	sxf = (tr.sx) * scalco; // set sx
	szf = 0.0;          //   set sz
	gzf = 0.0;
	fprintf(outcord,"%12.2f %12.2f %12.2f %12.2f\n",sxf,szf,tr.gx*scalco,gzf);
	if(maxntr<ntr-1)
	  maxntr=ntr-1;
	//输出第一炮的信息
	fprintf(outfp,"%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d\n",
		shotnumber,ntr-1,starttrace,starttrace+ntr-2,sx,sy,xmin,xmax,ymin,ymax);
	starttrace=k+1;
	sx=tr.sx;
	sy=tr.sy;
	ntr=1;
	if(tr.gx<tr.sx){
	  xmin=tr.gx;
	  xmax=tr.sx;
	}
	else{
	  xmin=tr.sx;
	  xmax=tr.gx;
	}
	if(tr.gy<tr.sy){
	  ymin=tr.gy;
	  ymax=tr.sy;
	}
	else{
	  ymin=tr.sy;
	  ymax=tr.gy;
	}
      }
    k++;
    if(tr.gx<xmin)
      xmin=tr.gx;
    if(tr.gx>xmax)
      xmax=tr.gx;
    if(tr.gy<ymin)
      ymin=tr.gy;
    if(tr.gy>ymax)
      ymax=tr.gy;
    if(xminall>xmin)
      xminall=xmin;
    if(xmaxall<xmax)
      xmaxall=xmax;
    if(yminall>ymin)
      yminall=ymin;
    if(ymaxall<ymax)
      ymaxall=ymax;
    if(k%1000==0) fprintf(stderr,"trcae=%-10d      finished(%) %5.1f\n",k,100.*k/traceall);
  }while(k<=traceall);//traceall是总道数，循环里面是读每一道的信息
  shotnumber++;
  if(CSGF) fclose(CSGF);
  if(maxntr<ntr-1)
    maxntr=ntr-1;
  fprintf(outfp,"%10d%10d%10d%10d%10d%10d%10d%10d%10d%10d\n",
	  shotnumber,ntr-1,starttrace,starttrace+ntr-2,sx,sy,xmin,xmax,ymin,ymax);
  
  rewind(outfp);
  fprintf(outfp,"%s\n",info1);
  fprintf(outfp,"%10d%10d%10d%10d%10d%10d%10d%10d\n",shotnumber,maxntr,nt,dt,xminall,xmaxall,yminall,ymaxall);
  fclose(infp);
  fclose(outfp);
  fclose(outcord);
}

long getfilesize(FILE *fp)
{
  if(fp==NULL) return -1;
  fseek(fp,0L,SEEK_END);
  return ftell(fp);
}

