#include "cirsym.h"

extern FILE *out,*inpa,*set;

extern char *b,*c,sp;
extern unsigned long leng,lengl,lepr;
extern int flag_fi3,flag_fi4,flag_fi5,flag_fsn,flag_nul,extract,
       flag_e,cirf,noequ,flag_pln,nuz,flag_sp,flag_cL,
       flag_mir,flag_opt,flag_mat,flag_two,flag_thr,
	   flag_2,flag_3,flag_4,flag_5,fl_g,fl_rg;
extern float range2,range3,range4,range5,
             range2p,range3p,range4p,range5p;

/* нахождение определителя схемы (при отсутствии приемников откликов) */
 void detanm(int act_n,SOURCE *act)
{
 printf("\n detan\n");
 allntr999(&act_n,act); /* нейтрализация источников воздействия */
 parallm(&act_n,act);  
 c[0]='\0';
 b=c;
 lengl=0;
 fputs("\ndetan= \n",out);
 gggfm(act_n,act);
 free(act);
 fputs(c,out);
 leng+=(b-c);
 lengl+=(b-c);
 fputs(";",out);
}
/* Проверка на связность */
 void bond1f(int n,GRAPH *a1,char s1,char s2,int *k)
{
 int i;
 char s4;

 for(i=0;i<n;i++)
  {
   if(a1[i].v1==s1)
    {
     if(a1[i].v2==s2)
      {
       a1[i].v1='z'; (*k)++;
       continue;
      }
     s4=a1[i].v2;
     a1[i].v1='z'; (*k)++;
     bond1f(n,a1,s4,s2,&*k);
    }
   if (a1[i].v1 != 'z' && a1[i].v2==s1)
    {
     if (a1[i].v1==s2)
      {
       a1[i].v1='z'; (*k)++;
       continue;
      }
     s4=a1[i].v1;
     a1[i].v1='z'; (*k)++;
     bond1f(n,a1,s4,s2,&*k);
    }
  }
}
/***********************************************************/
 void verstr(int n,GRAPH *matr,int *p,char *str)
{
 int i;
 char *st=str;

 (*p)=0;
 str[0]='\0';
 for (i=0;i<n;i++)
  {
   if (strchr(str,matr[i].v1)==NULL)
    {
     *st=matr[i].v1;
     *++st='\0';
     (*p)++;
    }
   if (strchr(str,matr[i].v2)==NULL)
    {
     *st=matr[i].v2;
     *++st='\0';
     (*p)++;
    }
  }
}
/***********************************************************/
void ver_n(int n,GRAPH *matr,int *p)
{
 char *str;

 str=(char*) malloc(2*n+1);
 verstr(n,matr,p,str);
 free(str);
}
/*********************************************************/
void node(int act_n,SOURCE *act,int *p)
{
 char *str;

 str=(char*) malloc(4*act_n+1);
 nodes(act_n,act,p,str);
 free(str);
}
/*************************************************************/
void nodesm(int act_n,SOURCE *act,int *p,char *str)
{
 int i;

 GRAPH *a2=( GRAPH *) malloc((act_n+1)*sizeof(GRAPH));

 for (i=0;i<act_n;i++)
  {
   a2[i].v1=act[i].v1;
   a2[i].v2=act[i].v3;
  }
 verstr(act_n,a2,p,str);
 free(a2);
}
/*************************************************************/
void nodes(int act_n,SOURCE *act,int *p,char *str)
{
 int j,na;

 GRAPH *a2=( GRAPH *) malloc((2*act_n+1)*sizeof(GRAPH));

 copyactm(act_n,act,a2);
 na = 2*act_n;
 verstr(na,a2,p,str);
 free(a2);
}
/*******************************/
void ster(unsigned long k)
{
 while(leng+strlen(c)>k)
  {
   b--;
   b[0]='\0';
   if (c>b)
     {
      int j,j1;
      fseek(out,-72,SEEK_CUR);
      for (j=0;j<70;j++)
	 {
	  c[j]=fgetc(out);
	  if (c[j]==' ') {j1=j-1; break;}
	  j1=j;
	 }
      c[j1+1]='\0';
      fseek(out,-j1-2,SEEK_CUR);
      fprintf(out,"%d",feof(out));
      fseek(out,-1,SEEK_CUR);
      b=c+j1+1;
      leng-=j1+1;
     }
  }
}
/*********************************/
void kontrol(void)
{
 char *a;
 int i,i0,j;

 while((b-c)>70)
{
 j=0;
 for (i=1;i<70;i++) if (*(c+70-i)=='(') {j=1; break;}
 if (!j) 
 for (i=1;i<70;i++) if (*(c+70-i)=='*' || *(c+70-i)=='+') break;
 i=70-i;
 a=c+i;
 for (j=0;j<i;j++) {fputc(c[j],out);}
 i0=0;
 for (j=i;j<70;j++) {i0++; fputc(' ',out);}
 leng+=70-i0;
 fputs("\n",out);
 strcpy(c,a);
 b-=70-i0;
 if (lepr > 50000) {printf("\n %ld ",leng); lepr=0;}
 lepr++;
}
}
/******************************************************/
void printa(char *t,int numb,int sign)
{
 int i,j;

 if (t[0] == '-' && 
    (abs(numb) == 1 || abs(numb) == 888 || abs(numb) == 555))
   {
    i=strlen(t);
    for (j=0;j<i;j++) t[j]=t[j+1];
    if (sign) sign=0; else sign=1;
   }
 if (abs(numb) == 1 || abs(numb) == 888 || abs(numb) == 555)
   {
    if (!sign) {strcat(b,t); b+=strlen(t);}
    else
      {
       if (*(b-1)=='+')
	 {
	  b--; b[0]='\0';
          strcat(b++,"-");
	 }
       else
       if (*(b-1)=='-')
	 {
	  b--; b[0]='\0';
          strcat(b++,"+");
	 }
       else strcat(b++,"-");
       strcat(b,t); b+=strlen(t);
      }
   }
 if (abs(numb) > 1 && abs(numb) != 888 && abs(numb) != 555)
   {
    if (!sign)
      {
       strcat(b++,"(");
       strcat(b,t); b+=strlen(t);
       strcat(b++,")");
      }
    else
      {
       if (*(b-1)=='+') {b--; b[0]='\0';}
       strcat(b,"-("); b+=2;
       strcat(b,t); b+=strlen(t);
       strcat(b++,")");
      }
   }
 kontrol();
}
/* Освобождение j-го элемента в структуре act */
 void fract(int j,int *act_n,SOURCE *act)
{
 int i;

 free(act[j].reb);
 (*act_n)--; for (i=j;i<*act_n;i++) act[i]=act[i+1];
}
/* Освобождение строк в структуре act */
 void freact(int act_n, SOURCE *act)
{
 int i;

 for (i=0;i<act_n;i++) free(act[i].reb);
}
/******************************************************/
void printi(char **t,char *a,int numb)
{
 if (abs(numb) == 1 || abs(numb) == 888 || abs(numb) == 555)
   {
    strcat(*t,a); (*t)+=strlen(a);
   }
 else
   {
    strcat((*t)++,"(");
    strcat((*t),a); (*t)+=strlen(a);
    strcat((*t)++,")");
   }
 strcat((*t)++,"*");
}
/*************************************************************/
void actcopy(int act_n,SOURCE *act,SOURCE *act1)
{
 int j;

 for (j=0;j<act_n;j++)
    {
     act1[j].v1=act[j].v1; act1[j].v2=act[j].v2;
     act1[j].v3=act[j].v3; act1[j].v4=act[j].v4;
     act1[j].kol=act[j].kol;
     act1[j].reb=strdup(act[j].reb);
    }
}
/*************************************************************/
void copyactm(int act_n,SOURCE *act,GRAPH *a2)
{
 int j,n=-1;
 for (j=0;j<act_n;j++)
    {
     a2[++n].v1=act[j].v1; a2[n].v2=act[j].v2;
     a2[++n].v1=act[j].v3; a2[n].v2=act[j].v4;
    }
}
/***************************************************************/
void uniact(char s1,char s2,int *act_n,SOURCE *act)
{
 int i,j;
 char sa,sb;

 if (s1 != '0') {sa=s1; sb=s2;}
 if (s1 == '0') {sa=s2; sb=s1;}
 for (i=0;i<*act_n;i++)
    {
     if (act[i].v1 == sa) act[i].v1=sb;
     if (act[i].v2 == sa) act[i].v2=sb;
     if (act[i].v3 == sa) act[i].v3=sb;
     if (act[i].v4 == sa) act[i].v4=sb;
     if (act[i].v1 == act[i].v2 || act[i].v3 == act[i].v4)
       {
	free(act[i].reb);
	(*act_n)--; for (j=i;j<*act_n;j++) act[j]=act[j+1];
	i--;
       }
    }
}
