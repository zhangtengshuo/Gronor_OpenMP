#define _GNU_SOURCE
#include <time.h>
#include <sys/types.h>
#include <string.h>
#include <sched.h>
double linux_cputime_(void)
{
  return ((double) clock()) / ((double) CLOCKS_PER_SEC);
}
double linux_walltime_(void)
{
    struct timeval time;
    int gettimeofday ( struct timeval *tp , int *tz );
    if (gettimeofday(&time,NULL)){return 0;}
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}
void linux_date_(char adate[10])
/* Return date as 'YY/MM/DD  ' */
{time_t t=time((time_t *) 0);
 char *tmp = ctime(&t);
 char mm[3];
 int i;
 adate[0]=*(tmp+22);
 adate[1]=*(tmp+23);
 adate[2]='/';
 for (i=0; i<3; i++) mm[i]=*(tmp+i+4);
 if(mm[0]=='J'&&mm[1]=='a') {adate[3]='0'; adate[4]='1';};
 if(mm[0]=='F')             {adate[3]='0'; adate[4]='2';};
 if(mm[0]=='M'&&mm[2]=='r') {adate[3]='0'; adate[4]='3';};
 if(mm[0]=='A'&&mm[1]=='p') {adate[3]='0'; adate[4]='4';};
 if(mm[0]=='M'&&mm[2]=='y') {adate[3]='0'; adate[4]='5';};
 if(mm[0]=='J'&&mm[2]=='n') {adate[3]='0'; adate[4]='6';};
 if(mm[0]=='J'&&mm[2]=='l') {adate[3]='0'; adate[4]='7';};
 if(mm[0]=='A'&&mm[1]=='u') {adate[3]='0'; adate[4]='8';};
 if(mm[0]=='S')             {adate[3]='0'; adate[4]='9';};
 if(mm[0]=='O')             {adate[3]='1'; adate[4]='0';};
 if(mm[0]=='N')             {adate[3]='1'; adate[4]='1';};
 if(mm[0]=='D')             {adate[3]='1'; adate[4]='2';};
 adate[5]='/';
 adate[6]=*(tmp+8); if(adate[6]==' ') adate[6]='0';
 adate[7]=*(tmp+9);
 adate[8]=' ';
 adate[9]=' ';}
void linux_time_(char atime[10])
/* Return time as 'HH:MM:SS  ' */
{time_t t=time((time_t *) 0);
 char *tmp = ctime(&t);
 int i;
 for(i=0; i<8; i++) atime[i]=*(tmp+i+11); }
int getcpucount_()
{cpu_set_t cs; int i; CPU_ZERO(&cs);
 sched_getaffinity(0,sizeof(cs),&cs);
 int count=0;
 for ( i=0; i<CPU_COUNT(&cs); i++){if(CPU_ISSET(i,&cs)) count++;}
 return count;
}
void setfirstonly_() {cpu_set_t cs; CPU_ZERO(&cs); CPU_SET(0,&cs); sched_setaffinity(0,sizeof(cs),&cs);}
void skipfirstonly_()
{cpu_set_t cs; int i; CPU_ZERO(&cs); for ( i=1; i<CPU_COUNT(&cs); i++){CPU_SET(i,&cs);}; sched_setaffinity(0,sizeof(cs),&cs); }
