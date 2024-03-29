﻿/* Вспомогательные функции и подпрограммы */

#include "cirsym.h"

extern FILE *out,*inpa,*set;

extern char *b,*c,sp;
extern unsigned long leng,lengl,lepr;
extern int flag_fi3,flag_fi4,flag_fi5,flag_fsn,extract,
       flag_e,cirf,noequ,flag_pln,nuz,flag_sp,flag_cL,st_pln;
extern float range2,range3,range4,range5;

/* Проверка наличия и выбор оптимального для разложения УИ */
/* int noideal(int act_n,SOURCE *act,int *first)
{
 int i;

 *first=-1;

 for (i=0;i<act_n;i++)
    if (act[i].kol > 0) {*first=i; return(1);}
	 
 for (i=0;i<act_n;i++)
    {
     if ((act[i].reb[0] == 'E' || act[i].reb[0] == 'J') && act[i].kol > 0
	 && (   act[i].v1==act[i].v3 || act[i].v1==act[i].v4
	     || act[i].v2==act[i].v3 || act[i].v2==act[i].v4))
       {*first=i; return(1);}
    }
 for (i=0;i<act_n;i++)
    {
     if ((act[i].reb[0] == 'E' || act[i].reb[0] == 'J') && act[i].kol > 0)
       {*first=i; return(1);}
    }
 for (i=0;i<act_n;i++)
    {
     if (act[i].reb[0] == 'B' && act[i].kol > 0
	 && (   act[i].v1==act[i].v3 || act[i].v1==act[i].v4
	     || act[i].v2==act[i].v3 || act[i].v2==act[i].v4))
       {*first=i; return(1);}
    }
 for (i=0;i<act_n;i++)
    {
     if (act[i].reb[0] == 'B' && act[i].kol > 0)
       {*first=i; return(1);}
    }
 for (i=0;i<act_n;i++)
    {
     if (act[i].reb[0] == 'K' && act[i].kol > 0
	 && (   act[i].v1==act[i].v3 || act[i].v1==act[i].v4
	     || act[i].v2==act[i].v3 || act[i].v2==act[i].v4))
       {*first=i; return(1);}
    }
 for (i=0;i<act_n;i++)
    {
     if (act[i].reb[0] == 'K' && act[i].kol > 0)
       {*first=i; return(1);}
    }
 for (i=0;i<act_n;i++)
    {
     if (act[i].reb[0] == 'H' && act[i].kol > 0
	 && (   act[i].v1==act[i].v3 || act[i].v1==act[i].v4
	     || act[i].v2==act[i].v3 || act[i].v2==act[i].v4))
       {*first=i; return(1);}
    }
 for (i=0;i<act_n;i++)
    {
     if (act[i].reb[0] == 'H' && act[i].kol > 0)
       {*first=i; return(1);}
    }
 for (i=0;i<act_n;i++)
    {
     if (act[i].reb[0] == 'G' && act[i].kol > 0
	 && (   act[i].v1==act[i].v3 || act[i].v1==act[i].v4
	     || act[i].v2==act[i].v3 || act[i].v2==act[i].v4))
       {*first=i; return(1);}
    }
 for (i=0;i<act_n;i++)
    {
     if (act[i].reb[0] == 'G' && act[i].kol > 0)
       {*first=i; return(1);}
    }
 for (i=0;i<act_n;i++)
    {
     if (act[i].kol > 0)
       {*first=i; return(1);}
    }
 return(0);
}*/
/* Функция проверки наличия в схеме проводимостей g, y, Y */
 int yescL_n(int n,PASSIVE *matr)
{
 int i,k=0;

 for (i=0;i<n;i++)
     if (   matr[i].reb[0]=='c' || matr[i].reb[0]=='L') k++;
 return(k);
}
 void nodestr(int n,PASSIVE *matr,int act_n,SOURCE *act,int *p,char *str)
{
 GRAPH *a2;
 int i,na;

 if (n+act_n==0) {*p=0; return;}
 a2=( GRAPH *) malloc((n+4*act_n+1)*sizeof(GRAPH));
 copypas(n,matr,a2);
 copyact(n,act_n,act,&na,a2);
 verstr(na,a2,p,str);
 if (yesmqt(act_n,act))
   {
    for (i=0;i<*p;i++)
       if (str[i]=='0') {free(a2); return;}
    str[*p]='0';
    (*p)++;
   }
 free(a2);
}
 void node_n(int n,PASSIVE *matr,int act_n,SOURCE *act,int *p)
{
 char *str;

 str=(char*) malloc(2*n+8*act_n+1);
 nodestr(n,matr,act_n,act,p,str);
 free(str);
}
/* void nodes(int act_n,SOURCE *act,int *p,char *str)
{
 GRAPH *a2;
 int i,na;

 if (!act_n) {*p=0; return;}
 a2=( GRAPH *) malloc(4*act_n*sizeof(GRAPH));
 copyact(0,act_n,act,&na,a2);
 verstr(na,a2,p,str);
 if (yesmqt(act_n,act))
   {
    for (i=0;i<*p;i++)
       if (str[i]=='0') {free(a2); return;}
    str[*p]='0';
    (*p)++;
   }
 free(a2);
}*/
/* Копирование структуры двухполюсников matr d в структуру matr1 */
 void pascopy(int n,PASSIVE *matr,PASSIVE *matr1)
{
 int j;

 for (j=0;j<n;j++)
    {
     matr1[j].v1=matr[j].v1; matr1[j].v2=matr[j].v2;
     matr1[j].kol=matr[j].kol;
     matr1[j].reb=strdup(matr[j].reb);
    }
}
/* Копирование структуры управляемых источников act в структуру act1 */
 
/* Копирование структуры двухполюсников matr d в структуру "граф" */
 void copypas(int n,PASSIVE *matr,GRAPH *a1)
{
 int i;

 for (i=0;i<n;i++)
    {
     a1[i].v1=matr[i].v1; a1[i].v2=matr[i].v2;
    }
}
/* Копирование структуры управляемых источников act в структуру "граф" */
 void copyact(int n,int act_n,SOURCE *act,int *na,GRAPH *a2)
{
 int j;

 n--;
 for (j=0;j<act_n;j++)
    {
     if (act[j].reb[0] == 'M')
       {
	a2[++n].v1=act[j].v1; a2[n].v2='0';
	a2[++n].v1=act[j].v2; a2[n].v2='0';
	a2[++n].v1=act[j].v3; a2[n].v2='0';
	a2[++n].v1=act[j].v4; a2[n].v2='0';
       }
     else
       if (act[j].reb[0] == 'Q')
	 {
	  a2[++n].v1=act[j].v1; a2[n].v2=act[j].v2;
	  a2[++n].v1=act[j].v3; a2[n].v2='0';
	  a2[++n].v1=act[j].v4; a2[n].v2='0';
	 }
       else
	 if (act[j].reb[0] == 'T')
	   {
	    a2[++n].v1=act[j].v1; a2[n].v2='0';
	    a2[++n].v1=act[j].v2; a2[n].v2='0';
	    a2[++n].v1=act[j].v3; a2[n].v2=act[j].v4;
	   }
	 else
	   {
	    a2[++n].v1=act[j].v1; a2[n].v2=act[j].v2;
	    a2[++n].v1=act[j].v3; a2[n].v2=act[j].v4;
	   }
    }
 *na=n+1;
}
/* Занесение в структуру act нового управляемого источника */
 void actall(char s1,char s2,char s3,char s4,
	     int kol,char *str,int *act_n,SOURCE *act)
{
 act[*act_n].v1=s1;
 act[*act_n].v2=s2;
 act[*act_n].v3=s3;
 act[*act_n].v4=s4;
 act[*act_n].kol=kol;
 act[*act_n].reb=strdup(str);
 (*act_n)++;
}
/* Занесение в структуру matr нового двухполюсного элемента */
 void matrall(char s1,char s2,int kol,char *str,int *n,PASSIVE *matr)
{
 matr[*n].v1=s1;
 matr[*n].v2=s2;
 matr[*n].kol=kol;
 matr[*n].reb=strdup(str);
 (*n)++;
}
/* Занесение в структуру act нового неудаляемого управляемого источника */
 void nuiall(char s1,char s2,char s3,char s4,int *act_n,SOURCE *act)
{
 char *w,w0[2];

 w0[0]='\0'; w=w0;
 strcat(w++,"1");
 act[*act_n].v1=s1;
 act[*act_n].v2=s2;
 act[*act_n].v3=s3;
 act[*act_n].v4=s4;
 act[*act_n].kol=-1;
 act[*act_n].reb=strdup(w0);
 (*act_n)++;
}
/* Занесение в структуру act нового вырожденного НУИ - "закоротки" */
 void nui(char s1,char s2,int *act_n,SOURCE *act)
{
 char *w,w0[2];

 w0[0]='\0'; w=w0;
 strcat(w++,"1");
 act[*act_n].v1=s1;
 act[*act_n].v2=s2;
 act[*act_n].v3=s1;
 act[*act_n].v4=s2;
 act[*act_n].kol=-1;
 act[*act_n].reb=strdup(w0);
 (*act_n)++;
}
/* Занесение в структуру act "закоротки" вместо k-го элемента  */
 void nuintr(char s1,char s2,int k,SOURCE *act)
{
 char *w,w0[2];

 w0[0]='\0'; w=w0;
 strcat(w++,"1");
 free(act[k].reb);
 act[k].v1=s1;
 act[k].v2=s2;
 act[k].v3=s1;
 act[k].v4=s2;
 act[k].kol=-1;
 act[k].reb=strdup(w0);
}
/* Освобождение структуры matr */
 void freematr(int n, PASSIVE *matr)
{
 int i;

 if (n)
   {
    for (i=0;i<n;i++) free(matr[i].reb); free(matr);
   }
}
/* Освобождение структуры act */
 void freeact(int act_n, SOURCE *act)
{
 int i;

 if (act_n)
   {
    for (i=0;i<act_n;i++) free(act[i].reb); free(act);
   }
}
/* Освобождение строк в структуре matr */
 void frematr(int n, PASSIVE *matr)
{
 int i;

 for (i=0;i<n;i++) free(matr[i].reb);
}
/* Разбор строки исходного cir-файла */
 void extra(int el_n,char *el,int be,int *ee,char *sr,int *k)
{
 *k=0;
  for (*ee=be+1;*ee<el_n;(*ee)++)
     {
      if ( (el[*ee] == ' ' || el[*ee] == '\t'  
	      || el[*ee]=='\n') && *ee == 1) break;
      if (el[*ee] != ' ' && el[*ee] != '\t')
        {
	 sr[(*k)++]=el[*ee];
         if (   el[*ee+1]==' ' || el[*ee+1]=='\t' 
		     || el[*ee+1]=='\n') break;
	}
     }
 sr[*k]='\0';
}
/******************************************************/
/*void printa(char *t,int numb,int sign)
{
 int i,j;

 if (t[0] == '-' && 
    (abs(numb) == 1 || abs(numb) == 888 || abs(numb) == 555))
   {
    i=strlen(t);
    for (j=0;j<i;j++) t[j]=t[j+1];
    if (sign) sign=0; else sign=1;
   }
 else  
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
}*/
 /******************************************************/
/*void printi(char **t,char *a,int numb)
{
 if (   abs(numb) == 1 || abs(numb) == 888
     || abs(numb) == 555
    )
   {
    strcat(*t,a); (*t)+=strlen(a);
   }
 else if (abs(numb) > 1 && abs(numb) != 888 && abs(numb) != 555)
   {
    strcat((*t)++,"(");
    strcat((*t),a); (*t)+=strlen(a);
    strcat((*t)++,")");
   }
 strcat((*t)++,"*");
}*/
/********************************************/
 void unipasp2(char s1,char s2,int *n,int n0,int *n3,PASSIVE *matr)
{
 int i,j;

 for (i=n0;i<*n;i++)
    {
     if (matr[i].v1 == s1) matr[i].v1=s2;
	if (matr[i].v2 == s1) matr[i].v2=s2;
	
//	 for (i=n0;i<*n;i++)
     if (   matr[i].v1 == matr[i].v2
         && (   matr[i].reb[0]=='g' || matr[i].reb[0]=='c'
	     || matr[i].reb[0]=='y' || matr[i].reb[0]=='Y'))
       {
		frmatr(i,n3,matr); i--; (*n)--;
	  }
 }
}
 void uniactp2(char s1,char s2,
                int *act_n,int act0_n,int *act3_n,SOURCE *act)
{
 int i,j;

 for (i=act0_n;i<*act_n;i++)
    {
     if (act[i].v1 == s1) act[i].v1=s2;
     if (act[i].v2 == s1) act[i].v2=s2;
     if (act[i].v3 == s1) act[i].v3=s2;
	 if (act[i].v4 == s1) act[i].v4=s2;
 // for (i=act0_n;i<*act_n;i++)
     if (   (act[i].v1 == act[i].v2 || act[i].v3 == act[i].v4)
         && act[i].kol > 0 && act[i].reb[0]=='G')      
       {
	  fract(i,act3_n,act); i--; (*act_n)--;
       }
 }
}
 void unipasp(char s1,char s2,int *n,int n0,PASSIVE *matr)
{
 int i,j;

 for (i=n0;i<*n;i++)
    {
     if (matr[i].v1 == s1) matr[i].v1=s2;
	if (matr[i].v2 == s1) matr[i].v2=s2;
	
	//i=n0-1;
//rep:	i++;
//	 for (i=n0;i<*n;i++)
     if (   matr[i].v1 == matr[i].v2
         && (   matr[i].reb[0]=='g' || matr[i].reb[0]=='c'
	     || matr[i].reb[0]=='y' || matr[i].reb[0]=='Y'))
       {
	 frmatr(i,n,matr); i--;
       }
	  // if (i<*n) goto rep;
 }
}
 void uniactp(char s1,char s2,int *act_n,int act0_n,SOURCE *act)
{
 int i,j;

 for (i=act0_n;i<*act_n;i++)
    {
     if (act[i].v1 == s1) act[i].v1=s2;
     if (act[i].v2 == s1) act[i].v2=s2;
     if (act[i].v3 == s1) act[i].v3=s2;
	 if (act[i].v4 == s1) act[i].v4=s2;
//i=act0_n-1;
//rep:	i++;
// for (i=act0_n;i<*act_n;i++)
	 if (   (act[i].v1 == act[i].v2 || act[i].v3 == act[i].v4)
         && act[i].kol > 0 && act[i].reb[0]=='G')  
		 {fract(i,act_n,act); i--;
       }
	 //  if (i<*act_n) goto rep;
 }
}
void unipas(char s1,char s2,int *n,PASSIVE *matr)
{
 int i,j;

 for (i=0;i<*n;i++)
    {
     if (matr[i].v1 == s1) matr[i].v1=s2;
	if (matr[i].v2 == s1) matr[i].v2=s2;
    // for (i=0;i<*n;i++)
    //{
	 if (   matr[i].v1 == matr[i].v2
         && (   matr[i].reb[0]=='g' || matr[i].reb[0]=='c'
	     || matr[i].reb[0]=='y' || matr[i].reb[0]=='Y'))
       {
	frmatr(i,n,matr); i--;
       }
    //}
 }	  
} 
/* Освобождение j-го элемента в структуре matr */
 void frmatr(int j,int *n,PASSIVE *matr)
{
 int i;

 free(matr[j].reb);
 (*n)--; for (i=j;i<*n;i++) matr[i]=matr[i+1];
}
/* Функция проверки наличия в схеме УИ */
 int yesd(int act_n,SOURCE *act)
{
 int i;

 for (i=0;i<act_n;i++)
 {
     if (act[i].kol < 0) continue; 
     // && act[i].kol != 888) 
	 return(1);
 }
 return(0);
}
/* Функция проверки наличия в схеме проводимостей g, y, Y */
 int yesc(int n,PASSIVE *matr)
{
 int i;

 for (i=0;i<n;i++)
     if (matr[i].reb[0]=='c') return(1);
 return(0);
}
/* Функция проверки наличия в схеме сопротивлений R, z, Z */
 int yesL(int n,PASSIVE *matr)
{
 int i;

 for (i=0;i<n;i++)
     if (matr[i].reb[0] == 'L') return(1);
 return(0);
}
/* Функция проверки наличия в схеме проводимостей g, y, Y */
 int yescL(int n,PASSIVE *matr)
{
 int i;

 for (i=0;i<n;i++)
     if (   matr[i].reb[0]=='c' || matr[i].reb[0]=='L') return(1);
 return(0);
}
/* Функция проверки наличия в схеме проводимостей g, y, Y */
 int yesg(int n,PASSIVE *matr)
{
 int i;

 for (i=0;i<n;i++)
     if (   matr[i].reb[0]=='g'
	 || matr[i].reb[0]=='y' || matr[i].reb[0]=='Y') return(1);
 return(0);
}
/*  Упрощение третьего в треугольнике c двумя генераторами или приемниками */
 int yeslm(int n,PASSIVE *matr,int act_n,SOURCE *act,
	    int *first,int *fm,char *sm1,char *sm2, int *smf)
{
 int i,j,k,m,pr=0;
 char s1,s2,t1,t2,x1,x2;

beg:
 for (i=0;i<act_n;i++)
    {
     if (   act[i].reb[0] == '1' && act[i].kol < 0
	 || !pr && act[i].reb[0] == 'Q' && act[i].kol < 0
	 ||  pr && act[i].reb[0] == 'T' && act[i].kol < 0)
       {
	if (!pr) {s1=act[i].v1; s2=act[i].v2;}
	else {s1=act[i].v3; s2=act[i].v4;}
	// if (s1=='0' || s2=='0') continue;
	for (j=0;j<act_n;j++)
	   {
	    if (i != j && act[j].reb[0] == 'M' && act[j].kol < 0)
	      {
	       if (!pr) {t1=act[j].v1; t2=act[j].v2;}
	       else {t1=act[j].v3; t2=act[j].v4;}
	       // if (t1=='0' || t2=='0') continue;
	       if (t1==s1) {x1=t2; x2=s2; goto sear;}
	       if (t1==s2) {x1=t2; x2=s1; goto sear;}
	       if (t2==s1) {x1=t1; x2=s2; goto sear;}
	       if (t2==s2) {x1=t1; x2=s1; goto sear;}
	      }
	    continue;
sear:
	    for (k=0;k<n;k++)
	       {
		if (   (matr[k].reb[0] == 'g' || matr[k].reb[0] == 'y'
		    || matr[k].reb[0] == 'Y' || matr[k].reb[0] == 'c')
		    && (   matr[k].v1==x1 && matr[k].v2==x2
			|| matr[k].v1==x2 && matr[k].v2==x1))
		  {
		   if (pr) {t1=act[j].v1; t2=act[j].v2;}
		   else {t1=act[j].v3; t2=act[j].v4;}
		   *first=k; *fm=j; *sm1=t1; *sm2=t2; *smf=pr;
		   return(1);
		  }
	       }
	    }
	}
    }
 pr++;
 if (pr==1) goto beg;
 return(0);
}
/*  Упрощение третьего в треугольнике c двумя генераторами или приемниками */
 int yeslq(int n,PASSIVE *matr,int act_n,SOURCE *act,
	    int *first,int *fm,char *sm1,char *sm2, int *smf)
{
 int i,j,k,m,pr=0;
 char s1,s2,t1,t2,x1,x2;

 for (i=0;i<act_n;i++)
    {
     if (   act[i].reb[0] == '1' && act[i].kol < 0
	 || act[i].reb[0] == 'T' && act[i].kol < 0)
       {
	s1=act[i].v3; s2=act[i].v4;
	if (s1=='0' || s2=='0') continue;
	for (j=0;j<act_n;j++)
	   {
	    if (i != j && act[j].reb[0] == 'Q' && act[j].kol < 0)
	      {
	       t1=act[j].v3; t2=act[j].v4;
	       if (t1=='0' || t2=='0') continue;
	       if (t1==s1) {x1=t2; x2=s2; goto sear;}
	       if (t1==s2) {x1=t2; x2=s1; goto sear;}
	       if (t2==s1) {x1=t1; x2=s2; goto sear;}
	       if (t2==s2) {x1=t1; x2=s1; goto sear;}
	      }
	    continue;
sear:
	    for (k=0;k<n;k++)
	       {
		if (   (matr[k].reb[0] == 'g' || matr[k].reb[0] == 'y'
		    || matr[k].reb[0] == 'Y' || matr[k].reb[0] == 'c')
		    && (   matr[k].v1==x1 && matr[k].v2==x2
			|| matr[k].v1==x2 && matr[k].v2==x1))
		  {
		   *first=k; *fm=j; *sm1=t1; *sm2=t2; *smf=pr;
		   return(1);
		  }
	       }
	   }
	}
    }
 return(0);
}
/*  Упрощение третьего в треугольнике c двумя генераторами или приемниками */
 int yeslt(int n,PASSIVE *matr,int act_n,SOURCE *act,
	    int *first,int *fm,char *sm1,char *sm2, int *smf)
{
 int i,j,k,m,pr=0;
 char s1,s2,t1,t2,x1,x2;

 for (i=0;i<act_n;i++)
    {
     if (   act[i].reb[0] == '1' && act[i].kol < 0
	 || act[i].reb[0] == 'Q' && act[i].kol < 0)
       {
	s1=act[i].v1; s2=act[i].v2;
	if (s1=='0' || s2=='0') continue;
	for (j=0;j<act_n;j++)
	   {
	    if (i != j && act[j].reb[0] == 'T' && act[j].kol < 0)
	      {
	       t1=act[j].v1; t2=act[j].v2;
	       if (t1=='0' || t2=='0') continue;
	       if (t1==s1) {x1=t2; x2=s2; goto sear;}
	       if (t1==s2) {x1=t2; x2=s1; goto sear;}
	       if (t2==s1) {x1=t1; x2=s2; goto sear;}
	       if (t2==s2) {x1=t1; x2=s1; goto sear;}
	      }
	    continue;
sear:
	    for (k=0;k<n;k++)
	       {
		if (   (matr[k].reb[0] == 'g' || matr[k].reb[0] == 'y'
		    || matr[k].reb[0] == 'Y' || matr[k].reb[0] == 'c')
		    && (   matr[k].v1==x1 && matr[k].v2==x2
			|| matr[k].v1==x2 && matr[k].v2==x1))
		  {
		   *first=k; *fm=j; *sm1=t1; *sm2=t2; *smf=pr;
		   return(1);
		  }
	       }
	   }
       }
    }
 return(0);
}
 int yesgm(int n,PASSIVE *matr,int act_n,SOURCE *act,
	   int *first,int *fm,char *sm1,char *sm2, int *smf)
{
 int i,j;
 char s1,s2,s3,s4;

 for (j=0;j<act_n;j++)
    {
     if (act[j].reb[0] == 'M' && act[j].kol==-1)
       {
	s1=act[j].v1; s2=act[j].v2; s3=act[j].v3; s4=act[j].v4;
	for (i=0;i<n;i++)
	   {
	    if (   matr[i].v1==s1 && matr[i].v2==s2
		|| matr[i].v1==s2 && matr[i].v2==s1)
	      {
	       *first=i; *fm=j; *sm1=s3; *sm2=s4; *smf=1; return(1);
	      }
	    if (   matr[i].v1==s3 && matr[i].v2==s4
		|| matr[i].v1==s4 && matr[i].v2==s3)
	      {
	       *first=i; *fm=j; *sm1=s1; *sm2=s2; *smf=0; return(1);
	      }
	   }
       }
    }
 return(0);
}
 int yesgq(int n,PASSIVE *matr,int act_n,SOURCE *act,
	   int *first,int *fm)
{
 int i,j;
 char s3,s4;

 for (j=0;j<act_n;j++)
    {
     if (act[j].reb[0] == 'Q' && act[j].kol==-1)
       {
	s3=act[j].v3; s4=act[j].v4;
	for (i=0;i<n;i++)
	   {
	    if (   matr[i].v1==s3 && matr[i].v2==s4
		|| matr[i].v1==s4 && matr[i].v2==s3)
	      {
	       *first=i; *fm=j; return(1);
	      }
	   }
       }
    }
 return(0);
}
 int yesgt(int n,PASSIVE *matr,int act_n,SOURCE *act,
	   int *first,int *fm)
{
 int i,j;
 char s1,s2;

 for (j=0;j<act_n;j++)
    {
     if (act[j].reb[0] == 'T' && act[j].kol==-1)
       {
	s1=act[j].v1; s2=act[j].v2;
	for (i=0;i<n;i++)
	   {
	    if (   matr[i].v1==s1 && matr[i].v2==s2
		|| matr[i].v1==s2 && matr[i].v2==s1)
	      {
	       *first=i; *fm=j; return(1);
	      }
	   }
       }
    }
 return(0);
}
/* Функция проверки наличия в схеме сопротивлений R, z, Z */
 int yesr(int n,PASSIVE *matr)
{
 int i;

 for (i=0;i<n;i++)
     if (   matr[i].reb[0] == 'R' 
     || matr[i].reb[0] == 'z'
	 || matr[i].reb[0] == 'Z') return(1);
 return(0);
}
/* Функция проверки наличия в схеме сопротивлений R, z, Z */
 int yesrL(int n,PASSIVE *matr)
{
 int i;

 for (i=0;i<n;i++)
     if (   matr[i].reb[0] == 'R' || matr[i].reb[0] == 'L' 
     || matr[i].reb[0] == 'z'
	 || matr[i].reb[0] == 'Z') return(1);
 return(0);
}
/* Функция проверки наличия в схеме ИТУН */
 int yess(int act_n,SOURCE *act)
{
 int i;

 for (i=0;i<act_n;i++)
     if (act[i].reb[0] == 'G') return(1);
 return(0);
}
/* Функция проверки наличия в схеме НУИ */
 int yesn(int act_n,SOURCE *act)
{
 int i;

 for (i=0;i<act_n;i++)
     if (act[i].kol < 0 &&
	     (   act[i].reb[0] == '1'
	      || (act[i].reb[0] != 'M' || act[i].reb[0] != 'Q'
	      || act[i].reb[0] != 'T'))
	) return(1);
 return(0);
}
 int yesmqt(int act_n,SOURCE *act)
{
 int i;

 for (i=0;i<act_n;i++)
      if (   act[i].reb[0] == 'M' || act[i].reb[0] == 'Q'
	  || act[i].reb[0] == 'T') return(1);
 return(0);
}
 int yesm(int act_n,SOURCE *act,int *first)
{
 int i;

 for (i=0;i<act_n;i++)
      if (act[i].reb[0] == 'M' && act[i].kol < 0)
	{
	 *first=i;
	 return(1);
	}
 return(0);
}
 int yesq(int act_n,SOURCE *act,int *first)
{
 int i;

 for (i=0;i<act_n;i++)
      if (act[i].reb[0] == 'Q' && act[i].kol < 0)
	{
	 *first=i;
	 return(1);
	}
 return(0);
}
 int yest(int act_n,SOURCE *act,int *first)
{
 int i;

 for (i=0;i<act_n;i++)
      if (act[i].reb[0] == 'T' && act[i].kol < 0)
	{
	 *first=i;
	 return(1);
	}
 return(0);
}
/* Нахождение оптимальной замкнутой на генераторы и приемники УИ y-ветви */
 void hangpas_p(int n,PASSIVE *matr,int act_n,SOURCE *act,int *first)
{
 int j,i,k,x,p,send,rec,activ=-1,f=0,k1,k2;
 char *str;
 GRAPH *a2;

 if (!n) return;
 a2=( GRAPH *) malloc(n*sizeof(GRAPH));
 str = (char*) malloc(2*n+4*act_n+1);
 copypas(n,matr,a2);
 verstr(n,a2,&p,str); free(a2);
 *first=-1;
 for (i=0;i<p;i++)
    {
     for (j=0,k=0;j<n;j++)
	 if (matr[j].v1 == str[i] || matr[j].v2 == str[i]) {x=j;k++;}
     if (k==1 && matr[x].reb[0] == 'c')
       {
	*first=x;
	send=0; rec=0;
	for (j=0;j<act_n;j++)
	   {
	    if (str[i]==act[j].v1 || str[i]==act[j].v2) {send++; k1=j;}
	    if (str[i]==act[j].v3 || str[i]==act[j].v4) {rec++; k2=j;}
	   }
	if (send==1 && rec==1) break;
	if (   send == 1
	    && (act[k1].reb[0] == 'K'
		|| act[k1].reb[0] == 'E' || act[k1].reb[0] == 'U'
		|| act[k1].reb[0] == 'H')) break;
	if (   rec == 1
	    && (act[k2].reb[0] == 'B'
	    || act[k2].reb[0] == 'J' || act[k2].reb[0] == 'I'
	    || act[k2].reb[0] == 'H')) break;
	if (abs(send-rec) >= activ) {activ=abs(send-rec);f=*first;}
	  else  *first=f;
       }
    }
 free(str);
}
/* Нахождение оптимальной замкнутой на генераторы и приемники УИ y-ветви */
 void hangpas(int n,PASSIVE *matr,int act_n,SOURCE *act,int *first)
{
 int j,i,k,x,p,send,rec,activ=-1,f=0,k1,k2;
 char *str;
 GRAPH *a2;

 if (!n) return;
 a2=( GRAPH *) malloc(n*sizeof(GRAPH));
 str = (char*) malloc(2*n+4*act_n+1);
 copypas(n,matr,a2);
 verstr(n,a2,&p,str); free(a2);
 *first=-1;
 for (i=0;i<p;i++)
    {
     for (j=0,k=0;j<n;j++)
	if (matr[j].v1 == str[i] || matr[j].v2 == str[i]) {x=j;k++;}
     if (k==1 && (   matr[x].reb[0] == 'g' || matr[x].reb[0] == 'y'
	     || matr[x].reb[0] == 'Y'))
	   {
	*first=x;
	send=0; rec=0;
	for (j=0;j<act_n;j++)
	   {
	    if (str[i]==act[j].v1 || str[i]==act[j].v2) {send++; k1=j;}
	    if (str[i]==act[j].v3 || str[i]==act[j].v4) {rec++; k2=j;}
	   }
	if (send==1 && rec==1) break;
	if (   send == 1
	    && (act[k1].reb[0] == 'K'
		|| act[k1].reb[0] == 'E' || act[k1].reb[0] == 'U'
		|| act[k1].reb[0] == 'H')) break;
	if (   rec == 1
	    && (act[k2].reb[0] == 'B'
	    || act[k2].reb[0] == 'J' || act[k2].reb[0] == 'I'
	    || act[k2].reb[0] == 'H')) break;
	if (abs(send-rec) >= activ) {activ=abs(send-rec);f=*first;}
	  else  *first=f;
       }
    }
 free(str);
}
/* Выбор оптимальной для разложения y-ветви (проводимости) first */
 void choiceg_p(int n,PASSIVE *matr,int act_n,SOURCE *act,int *first)
{
 int j,i,s1,s2,pow1,pow2,pow0=100,pow12=100,t1,t2,k,number=0,
     send1,send2,rec1,rec2,f1;

*first=0;
 for (i=0;i<n;i++)
    {
     if (matr[i].reb[0]!='c') continue;
    s1=matr[i].v1; s2=matr[i].v2;
/*    send1=0; send2=0; rec1=0; rec2=0;
    for (j=0;j<act_n;j++)
       { 
	if (   act[j].reb[0]=='1' // || act[j].reb[0]=='M'
	    // || act[j].reb[0]=='Q' || act[j].reb[0]=='T'
		)
	  {
	   if (s1 == act[j].v1 || s1 == act[j].v2) send1++;
	   if (s1 == act[j].v3 || s1 == act[j].v4) rec1++;
	   if (s2 == act[j].v1 || s2 == act[j].v2) send2++;
	   if (s2 == act[j].v3 || s2 == act[j].v4) rec2++;
	  }
       }*/ 
	   /* rec1, rec2 - число приемников, инцидентных узлам ветви */
    pow1=0; pow2=0;
    for (k=0;k<n;k++)
       { /* pow1, pow2 - степени узлов рассматриваемой y-ветви */
	if (k==i) continue;
	t1=matr[k].v1; t2=matr[k].v2;
	if (s1==t1 || s1==t2) pow1++;
	if (s2==t1 || s2==t2) pow2++;
       }
 //   if (send1 && rec1) pow1 = pow1+2;
 //   if (send2 && rec2) pow2 = pow2+2;
 //   if (!send1 && rec1 || send1 && !rec1) pow1 = pow1 - 2;
 //   if (!send2 && rec2 || send2 && !rec2) pow2 = pow2 - 2;

    if (   pow1 < pow0 || pow2 < pow0)
      {
       if (pow1 < pow0) {pow0 = pow1; *first=i;}
       if (pow2 < pow0) {pow0 = pow2; *first=i;}
       if (pow1+pow2 < pow12) {pow12=pow1+pow2; *first=i;}
       if (pow1+pow2 <= pow12 && matr[i].kol > number)
	 {number=matr[i].kol; *first=i;}
      }
   }
 }
/* Выбор оптимальной для разложения y-ветви (проводимости) first */
 void choiceg(int n,PASSIVE *matr,int act_n,SOURCE *act,int *first)
{
 int j,i,s1,s2,pow1,pow2,pow0=100,pow12=100,t1,t2,k,number=0,
     send1,send2,rec1,rec2,f1;

*first=0;
 for (i=0;i<n;i++)
    {
     if (   matr[i].reb[0]=='R' || matr[i].reb[0]=='z'
	 || matr[i].reb[0]=='Z' || matr[i].reb[0]=='L'
	 || matr[i].reb[0]=='c') continue;
    s1=matr[i].v1; s2=matr[i].v2;
/*    send1=0; send2=0; rec1=0; rec2=0;
    for (j=0;j<act_n;j++)
       { 
	if (   act[j].reb[0]=='1' // || act[j].reb[0]=='M'
	    // || act[j].reb[0]=='Q' || act[j].reb[0]=='T'
		)
	  {
	   if (s1 == act[j].v1 || s1 == act[j].v2) send1++;
	   if (s1 == act[j].v3 || s1 == act[j].v4) rec1++;
	   if (s2 == act[j].v1 || s2 == act[j].v2) send2++;
	   if (s2 == act[j].v3 || s2 == act[j].v4) rec2++;
	  }
       } */
    pow1=0; pow2=0;
    for (k=0;k<n;k++)
       { /* pow1, pow2 - степени узлов рассматриваемой y-ветви */
	if (k==i) continue;
	t1=matr[k].v1; t2=matr[k].v2;
	if (s1==t1 || s1==t2) pow1++;
	if (s2==t1 || s2==t2) pow2++;
       }
//    if (send1 && rec1) pow1 = pow1+2;
//    if (send2 && rec2) pow2 = pow2+2;
//    if (!send1 && rec1 || send1 && !rec1) pow1 = pow1 - 2;
//    if (!send2 && rec2 || send2 && !rec2) pow2 = pow2 - 2;

    if (   pow1 < pow0 || pow2 < pow0)
      {
       if (pow1 < pow0) {pow0 = pow1; *first=i;}
       if (pow2 < pow0) {pow0 = pow2; *first=i;}
       if (pow1+pow2 < pow12) {pow12=pow1+pow2; *first=i;}
       if (pow1+pow2 <= pow12 && matr[i].kol > number)
	 {number=matr[i].kol; *first=i;}
      }
   }
}
/* Выбор минимального числа из двух целых чисел */
 int minimum(int pow1,int pow2)
{
 if (pow1>pow2) return(pow2); else return(pow1);
}
/* Выбор максимального числа из двух целых чисел */
 int maximum(int activ1,int activ2)
{
 if (activ1>activ2) return(activ1); else return(activ2);
}
/* Выбор оптимальной для разложения z-ветви (сопротивления) - first */
 void choicer_p(int n,PASSIVE *matr,int act_n,SOURCE *act,int *first)
{
 int j,i,s1,s2,pow1,pow2,pow10=0,pow20=0,t1,t2,k,number=0,
     send1,send2,rec1,rec2,activ1,activ2,activ10=0,activ20=0; //,ind=0;

 *first=0;
 for (i=0;i<n;i++)
    { /* begin branch */
    if (matr[i].reb[0]!='L') continue;
    s1=matr[i].v1; s2=matr[i].v2;
    send1=0; send2=0; rec1=0; rec2=0;
    for (j=0;j<act_n;j++)
       { /* send1, send2 - число генераторов, инцидентных узлам ветви */
	if (s1 == act[j].v1 || s1 == act[j].v2) send1++;
	if (s1 == act[j].v3 || s1 == act[j].v4) rec1++;
	if (s2 == act[j].v1 || s2 == act[j].v2) send2++;
	if (s2 == act[j].v3 || s2 == act[j].v4) rec2++;
       } /* rec1, rec2 - число приемников, инцидентных узлам ветви */
    activ1=abs(send1-rec1);
    activ2=abs(send2-rec2);
    pow1=0; pow2=0;
    for (k=0;k<n;k++)
       { /* pow1, pow2 - степени узлов рассматриваемой z-ветви */
	t1=matr[k].v1; t2=matr[k].v2;
	if (s1==t1 || s1==t2) pow1++;
	if (s2==t1 || s2==t2) pow2++;
       }
    if (maximum(pow1,pow2) > maximum(pow10,pow20))
      { /* если рассматриваемая ветвь предпочтительнее ранее выбранной */
       pow10=pow1; pow20=pow2;
       number=matr[i].kol;
       activ10=activ1; activ20=activ2;
       *first=i;
      // ind=1;
      }
    else
      if (maximum(pow1,pow2) == maximum(pow10,pow20)
	  && pow1+pow2 > pow10+pow20 )
	{ /* если рассматриваемая ветвь предпочтительнее ранее выбранной */
	 pow10=pow1; pow20=pow2;
	 number=matr[i].kol;
	 activ10=activ1; activ20=activ2;
	 *first=i;
	// ind=1;
	}
    else
      if (maximum(pow1,pow2) == maximum(pow10,pow20)
	  && pow1+pow2 == pow10+pow20
	  && matr[i].kol > number )
	 { /* если рассматриваемая ветвь предпочтительнее ранее выбранной */
	  pow10=pow1; pow20=pow2;
	  number=matr[i].kol;
	  activ10=activ1; activ20=activ2;
	  *first=i;
	 // ind=1;
	 }
    else
      if (maximum(pow1,pow2) == maximum(pow10,pow20)
	  && pow1+pow2 == pow10+pow20
	  && matr[i].kol == number
	  && maximum(activ1,activ2) > maximum(activ10,activ20) )
	{ /* если рассматриваемая ветвь предпочтительнее ранее выбранной */
	 pow10=pow1; pow20=pow2;
	 number=matr[i].kol;
	 activ10=activ1; activ20=activ2;
	 *first=i;
	 // ind=1;
	}
    } /* end branch */
/* if (ind) return;
 for (i=0;i<n;i++)
    {
    if (   matr[i].reb[0]=='g' || matr[i].reb[0]=='c'
	|| matr[i].reb[0]=='y' || matr[i].reb[0]=='Y'
       ) continue;
    s1=matr[i].v1; s2=matr[i].v2;
    send1=0; send2=0; rec1=0; rec2=0;
    for (j=0;j<act_n;j++)
       {
	if (s1 == act[j].v1 || s1 == act[j].v2) send1++;
	if (s1 == act[j].v3 || s1 == act[j].v4) rec1++;
	if (s2 == act[j].v1 || s2 == act[j].v2) send2++;
	if (s2 == act[j].v3 || s2 == act[j].v4) rec2++;
       }
    activ1=abs(send1-rec1);
    activ2=abs(send2-rec2);
    pow1=0; pow2=0;
    for (k=0;k<n;k++)
       { 
	t1=matr[k].v1; t2=matr[k].v2;
	if (s1==t1 || s1==t2) pow1++;
	if (s2==t1 || s2==t2) pow2++;
       }
    if (maximum(pow1,pow2) > maximum(pow10,pow20))
      { 
       pow10=pow1; pow20=pow2;
       number=matr[i].kol;
       activ10=activ1; activ20=activ2;
       *first=i;
      }
    else
      if (maximum(pow1,pow2) == maximum(pow10,pow20)
	  && pow1+pow2 > pow10+pow20 )
	{ 
	 pow10=pow1; pow20=pow2;
	 number=matr[i].kol;
	 activ10=activ1; activ20=activ2;
	 *first=i;
	}
    else
      if (maximum(pow1,pow2) == maximum(pow10,pow20)
	  && pow1+pow2 == pow10+pow20
	  && matr[i].kol > number )
	 { 
      pow10=pow1; pow20=pow2;
	  number=matr[i].kol;
	  activ10=activ1; activ20=activ2;
	  *first=i;
	 }
    else
      if (maximum(pow1,pow2) == maximum(pow10,pow20)
	  && pow1+pow2 == pow10+pow20
	  && matr[i].kol == number
	  && maximum(activ1,activ2) > maximum(activ10,activ20) )
	{ 
	 pow10=pow1; pow20=pow2;
	 number=matr[i].kol;
	 activ10=activ1; activ20=activ2;
	 *first=i;
	}
    } */
 }
/* Выбор оптимальной для разложения z-ветви (сопротивления) - first */
 void choicer(int n,PASSIVE *matr,int act_n,SOURCE *act,int *first)
{
 int j,i,s1,s2,pow1,pow2,pow10=0,pow20=0,t1,t2,k,number=0,
     send1,send2,rec1,rec2,activ1,activ2,activ10=0,activ20=0,ind=0;

 *first=0;
 for (i=0;i<n;i++)
    { /* begin branch */
    if (   matr[i].reb[0]=='g' || matr[i].reb[0]=='c'
	|| matr[i].reb[0]=='y' || matr[i].reb[0]=='Y'
	|| matr[i].reb[0]=='L'
       ) continue;
    s1=matr[i].v1; s2=matr[i].v2;
    send1=0; send2=0; rec1=0; rec2=0;
    for (j=0;j<act_n;j++)
       { /* send1, send2 - число генераторов, инцидентных узлам ветви */
	if (s1 == act[j].v1 || s1 == act[j].v2) send1++;
	if (s1 == act[j].v3 || s1 == act[j].v4) rec1++;
	if (s2 == act[j].v1 || s2 == act[j].v2) send2++;
	if (s2 == act[j].v3 || s2 == act[j].v4) rec2++;
       } /* rec1, rec2 - число приемников, инцидентных узлам ветви */
    activ1=abs(send1-rec1);
    activ2=abs(send2-rec2);
    pow1=0; pow2=0;
    for (k=0;k<n;k++)
       { /* pow1, pow2 - степени узлов рассматриваемой z-ветви */
	t1=matr[k].v1; t2=matr[k].v2;
	if (s1==t1 || s1==t2) pow1++;
	if (s2==t1 || s2==t2) pow2++;
       }
    if (maximum(pow1,pow2) > maximum(pow10,pow20))
      { /* если рассматриваемая ветвь предпочтительнее ранее выбранной */
       pow10=pow1; pow20=pow2;
       number=matr[i].kol;
       activ10=activ1; activ20=activ2;
       *first=i;
       ind=1;
      }
    else
      if (maximum(pow1,pow2) == maximum(pow10,pow20)
	  && pow1+pow2 > pow10+pow20 )
	{ /* если рассматриваемая ветвь предпочтительнее ранее выбранной */
	 pow10=pow1; pow20=pow2;
	 number=matr[i].kol;
	 activ10=activ1; activ20=activ2;
	 *first=i;
	 ind=1;
	}
    else
      if (maximum(pow1,pow2) == maximum(pow10,pow20)
	  && pow1+pow2 == pow10+pow20
	  && matr[i].kol > number )
	 { /* если рассматриваемая ветвь предпочтительнее ранее выбранной */
	  pow10=pow1; pow20=pow2;
	  number=matr[i].kol;
	  activ10=activ1; activ20=activ2;
	  *first=i;
	  ind=1;
	 }
    else
      if (maximum(pow1,pow2) == maximum(pow10,pow20)
	  && pow1+pow2 == pow10+pow20
	  && matr[i].kol == number
	  && maximum(activ1,activ2) > maximum(activ10,activ20) )
	{ /* если рассматриваемая ветвь предпочтительнее ранее выбранной */
	 pow10=pow1; pow20=pow2;
	 number=matr[i].kol;
	 activ10=activ1; activ20=activ2;
	 *first=i;
	 ind=1;
	}
    } /* end branch */
 if (ind) return;
 for (i=0;i<n;i++)
    { /* begin branch */
    if (   matr[i].reb[0]=='g' || matr[i].reb[0]=='c'
	|| matr[i].reb[0]=='y' || matr[i].reb[0]=='Y'
       ) continue;
    s1=matr[i].v1; s2=matr[i].v2;
    send1=0; send2=0; rec1=0; rec2=0;
    for (j=0;j<act_n;j++)
       { /* send1, send2 - число генераторов, инцидентных узлам ветви */
	if (s1 == act[j].v1 || s1 == act[j].v2) send1++;
	if (s1 == act[j].v3 || s1 == act[j].v4) rec1++;
	if (s2 == act[j].v1 || s2 == act[j].v2) send2++;
	if (s2 == act[j].v3 || s2 == act[j].v4) rec2++;
       } /* rec1, rec2 - число приемников, инцидентных узлам ветви */
    activ1=abs(send1-rec1);
    activ2=abs(send2-rec2);
    pow1=0; pow2=0;
    for (k=0;k<n;k++)
       { /* pow1, pow2 - степени узлов рассматриваемой z-ветви */
	t1=matr[k].v1; t2=matr[k].v2;
	if (s1==t1 || s1==t2) pow1++;
	if (s2==t1 || s2==t2) pow2++;
       }
    if (maximum(pow1,pow2) > maximum(pow10,pow20))
      { /* если рассматриваемая ветвь предпочтительнее ранее выбранной */
       pow10=pow1; pow20=pow2;
       number=matr[i].kol;
       activ10=activ1; activ20=activ2;
       *first=i;
      }
    else
      if (maximum(pow1,pow2) == maximum(pow10,pow20)
	  && pow1+pow2 > pow10+pow20 )
	{ /* если рассматриваемая ветвь предпочтительнее ранее выбранной */
	 pow10=pow1; pow20=pow2;
	 number=matr[i].kol;
	 activ10=activ1; activ20=activ2;
	 *first=i;
	}
    else
      if (maximum(pow1,pow2) == maximum(pow10,pow20)
	  && pow1+pow2 == pow10+pow20
	  && matr[i].kol > number )
	 { /* если рассматриваемая ветвь предпочтительнее ранее выбранной */
	  pow10=pow1; pow20=pow2;
	  number=matr[i].kol;
	  activ10=activ1; activ20=activ2;
	  *first=i;
	 }
    else
      if (maximum(pow1,pow2) == maximum(pow10,pow20)
	  && pow1+pow2 == pow10+pow20
	  && matr[i].kol == number
	  && maximum(activ1,activ2) > maximum(activ10,activ20) )
	{ /* если рассматриваемая ветвь предпочтительнее ранее выбранной */
	 pow10=pow1; pow20=pow2;
	 number=matr[i].kol;
	 activ10=activ1; activ20=activ2;
	 *first=i;
	}
    } /* end branch */
 }