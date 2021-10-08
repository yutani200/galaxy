#include <stdio.h>
#include <stdlib.h>
//void output(int n,double T,double m[],double x[][3],double vel[][3],double acc[][3],double H,double dt,double eps2);
void output(int n,double T,double dH,double dt,Istar *star){
  char filename[100];
  FILE *fp;
  sprintf(filename,"./data/%s.%03d.txt","star",(int)T);
  fp = fopen(filename,"w");
  for(int i = 0;i<n;i++){
    fprintf(fp,"%lf %lf %lf %lf ",T,dt,star[i].eps,star[i].m);
    fprintf(fp,"%lf %lf %lf ",star[i].x[0],star[i].x[1],star[i].x[2]);
    fprintf(fp,"%lf %lf %lf ",star[i].v[0],star[i].v[1],star[i].v[2]);
    fprintf(fp,"%lf %lf %lf ",star[i].a[0],star[i].a[1],star[i].a[2]);
    fprintf(fp,"%.16lf\n",dH);
  }
  fclose(fp);
  return;
}
