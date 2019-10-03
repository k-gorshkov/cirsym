/* Управляющий модуль: ввод схемы, анализ и диагностика */

#include "cirsym.h"

extern FILE *out,*inpa,*set;

extern char *b, *c, sp;
extern unsigned char powr[100];
extern unsigned long leng,lengl,lepr;
extern int flag_fi3,flag_fi4,flag_fi5,flag_fsn,extract,
           fl_pow,fl_rg,fl_sp,flag_e,cirf,noequ,flag_pln,
		   nuz,flag_sp,flag_cL,st_pln,flag_out,fl_ran,ranmax;
extern float range2,range3,range4,range5;

/**************************************************************/
 void detanp(int n,PASSIVE *matr,int act_n,SOURCE *act,int ncL)
{
 int j,k,k1;
 char *str;
 PASSIVE *matr1;
 SOURCE *act1;

 printf("\n detan\n");
 c[0]='\0';
 b=c;
 if (flag_e)
  { /* простейшие упрощения и вывод подформул */
    fputs("\n",out);
    reallyy(&n,matr);
    seqzz(&n,matr,act_n,act);
	red2zz(&n,matr);
	if (fl_rg)
	  {red2yy(&n,matr,act_n,act); reallyy(&n,matr);
       reallyz(&n,matr); seqzy(&n,matr,act_n,act);
	   seqzz(&n,matr,act_n,act);}
    strcat(b++,"\n");
    kontrol();
    fputs(c,out);
   }
 else {
 reallf(&n,matr);
 parall(&act_n,act); }
 allntr999(&act_n,act); /* нейтрализация источников воздействия */
 fputs("\ndetan=",out);
  if (fl_ran == -1 && ranmax == -1)
 {
  for (st_pln=0;st_pln<=ncL;st_pln++)
   {
     c[0]='\0';
     b=c;
     lengl=0;
     for (j=0;j<st_pln;j++) fprintf(out,"%c*",sp);
     printf("\n%c**%d*\n",sp,st_pln);
     fputs("\n(\n",out);
     matr1=(PASSIVE *) malloc (n*sizeof(PASSIVE));
     act1=(SOURCE *) malloc (act_n*sizeof(SOURCE));
     pascopy(n,matr,matr1);
     actcopy(act_n,act,act1);
     gggp(n,matr1,act_n,act1,0,st_pln);
     free(matr1); free(act1);
     fputs(c,out);
     leng+=(b-c);
     lengl+=(b-c);
     if (!(b-c)) fputs(" 0",out),leng++,lengl++;
     fputs("\n)",out);
     if (st_pln != ncL) fputs("\n+",out);
    }
     freematr(n,matr); freeact(act_n,act);
	}
  else if (fl_ran != -1 && ranmax != -1)
 {
  for (st_pln=fl_ran;st_pln<=ranmax;st_pln++)
   {
     c[0]='\0';
     b=c;
     lengl=0;
     for (j=0;j<st_pln;j++) fprintf(out,"%c*",sp);
     printf("\n%c**%d*\n",sp,st_pln);
     fputs("\n(\n",out);
     matr1=(PASSIVE *) malloc (n*sizeof(PASSIVE));
     act1=(SOURCE *) malloc (act_n*sizeof(SOURCE));
     pascopy(n,matr,matr1);
     actcopy(act_n,act,act1);
     gggp(n,matr1,act_n,act1,0,st_pln);
     free(matr1); free(act1);
     fputs(c,out);
     leng+=(b-c);
     lengl+=(b-c);
     if (!(b-c)) fputs(" 0",out),leng++,lengl++;
     fputs("\n)",out);
     if (st_pln != ranmax) fputs("\n+",out);
    }
     freematr(n,matr); freeact(act_n,act);
	}
 else if (   fl_ran == -1 && ranmax != -1
          || fl_ran != -1 && ranmax == -1)
{
 k1=strlen(powr)-1;
 for (k=0;k<k1;k++)
	{
     if (powr[k]==' ') continue;
	 st_pln=powr[k]-48;
     printf("\n        degree =  %d\n",st_pln);
     c[0]='\0';
     b=c;
     lengl=0;
	 for (j=0;j<st_pln;j++) fprintf(out,"%c*",sp);
     printf("\n%c**%d*\n",sp,st_pln);
     fputs("\n(\n",out);
     matr1=(PASSIVE *) malloc (n*sizeof(PASSIVE));
     act1=(SOURCE *) malloc (act_n*sizeof(SOURCE));
     pascopy(n,matr,matr1);
     actcopy(act_n,act,act1);
     gggp(n,matr1,act_n,act1,0,st_pln);
     free(matr1); free(act1);
     fputs(c,out);
     leng+=(b-c);
     lengl+=(b-c);
     if (!(b-c)) fputs(" 0",out),leng++,lengl++;
     fputs("\n)",out);
 	}
  freematr(n,matr); freeact(act_n,act);
 }
  fputs("\n;",out);
 }
/* нахождение определителя схемы (при отсутствии приемников откликов) */
 void detan(int n,PASSIVE *matr,int act_n,SOURCE *act)
{
 printf("\n detan\n");
 allntr999(&act_n,act); /* нейтрализация источников воздействия */
 c[0]='\0';
 b=c;
  if (flag_e)
  { /* простейшие упрощения и вывод подформул */
    fputs("\n",out);
    reallyy(&n,matr);
    parallm(&act_n,act);
    seqzz(&n,matr,act_n,act);
	red2zz(&n,matr);
	if (fl_rg)
	  {red2yy(&n,matr,act_n,act); reallyy(&n,matr);
       reallyz(&n,matr); seqzy(&n,matr,act_n,act);
	   seqzz(&n,matr,act_n,act);}
    strcat(b++,"\n");
    kontrol();
    fputs(c,out);
   }
 else {
 reallf(&n,matr);
 parall(&act_n,act); }
 c[0]='\0';
 b=c;
 lengl=0;
 fputs("\ndetan= \n",out);
 gggf(n,matr,act_n,act);
 free(matr); free(act);
 fputs(c,out);
 leng+=(b-c);
 lengl+=(b-c);
 fputs(";",out);
}
/***************************************************************************/
/*            Формирование схемной функции без единого сечения             */
/*                      для  числителя и знаменателя                       */
 void cirfunst(int n,PASSIVE *matr,int act_n,SOURCE *act,PASSIVE *diagui)
{
 int j,ej;
 char o1,o2,*w,w0[5];
 PASSIVE *matr0; SOURCE *act0;

 c[0]='\0';
 b=c;
 if (flag_e)
   { /* простейшие упрощения и вывод подформул */
    fputs("\n",out);
    reallyy(&n,matr);
    parallm(&act_n,act);
    seqzz(&n,matr,act_n,act);
	red2zz(&n,matr);
	if (fl_rg)
	  {red2yy(&n,matr,act_n,act); reallyy(&n,matr);
       reallyz(&n,matr); seqzy(&n,matr,act_n,act);
	   seqzz(&n,matr,act_n,act);}
    strcat(b++,"\n");
    kontrol();
    fputs(c,out);
   }
 else
   {reallf(&n,matr); parall(&act_n,act);}
 for (j=0;j<act_n;j++) /* ej - номер УИ, задающего вход */
     if (act[j].reb[0]=='E' || act[j].reb[0]=='J')
       {ej=j; break;}
 o1=diagui[0].v1; o2=diagui[0].v2; /* узлы искомого отклика */
 act[ej].v3=o2; /* перенос управления на искомый отклик */
 act[ej].v4=o1;
 matr0=(PASSIVE *) malloc (n*sizeof(PASSIVE));
 act0=(SOURCE *) malloc (act_n*sizeof(SOURCE));
 pascopy(n,matr,matr0);
 actcopy(act_n,act,act0);
 if (flag_out)
   {
    fprintf(out,"\n%% Input: %s (%d,%d), output: %s (%d,%d)\n",
            act[ej].reb,act[ej].v1-48,act[ej].v2-48,
	    diagui[0].reb,diagui[0].v1-48,diagui[0].v2-48);
    printf("\n Input: %s (%d,%d), output: %s (%d,%d)\n",
           act[ej].reb,act[ej].v1-48,act[ej].v2-48,
	   diagui[0].reb,diagui[0].v1-48,diagui[0].v2-48);
   }
 c[0]='\0';
 b=c;
 fprintf(out,"\n%s = %s*\n(\n",diagui[0].reb,act0[ej].reb);
 printf("\n%s num\n",diagui[0].reb);
 lengl=0;
 /* преобразование УИ и взвешенного НУИ в единичный НУИ */
 act0[ej].kol=-1;
 w0[0]='\0'; w=w0;
 strcat(w++,"1");
 free(act0[ej].reb);
 act0[ej].reb=strdup(w0);
 gggf(n,matr0,act_n,act0);
 free(matr0); free(act0);
 fputs(c,out);
 leng+=(b-c);
 lengl+=(b-c);
 if (!(b-c)) fputs(" 0",out),leng++,lengl++;
 fputs("\n)\n/\n(\n",out);
 printf("\n%s den\n",diagui[0].reb);
 allntr999(&act_n,act); /* нейтрализация многомерного НУИ */
 c[0]='\0';
 b=c;
 lengl=0;
 gggf(n,matr,act_n,act);
 free(matr); free(act);
 fputs(c,out);
 leng+=(b-c);
 lengl+=(b-c);
 if (!(b-c)) fputs(" 0",out),leng++,lengl++;
 fputs("\n);",out);
}
 /***************************************************************************/
/*            Формирование полиномиальной схемной функции                  */
/*                                                                         */
 void cirfunpa(int n,PASSIVE *matr,int act_n,SOURCE *act,
	       PASSIVE *diagui,int ncL)
{
 int i,j,k,ej,max,min;
 char v1,v2,v3,v4,*w,w0[100];
 long unsigned dl;
 PASSIVE *matr1,*matr2;
 SOURCE *act1,*act2;

 c[0]='\0';
 b=c;
 v3=diagui[0].v1; v4=diagui[0].v2; /* узлы искомого отклика */
 for (j=0;j<act_n;j++)
    {
     if (act[j].kol==999)
       {
	act[j].v3=v3; /* перенос управления на искомый отклик */
	act[j].v4=v4;
	ej=j; break;
       }
    }
 if (flag_out) fprintf(out,"\n%% Output: %s (%d,%d)\n",
	               diagui[0].reb,diagui[0].v1-48,diagui[0].v2-48);
 fprintf(out,"\n%s = \n",diagui[0].reb);
 v1=act[ej].v1;
 v2=act[ej].v2;
 if (n) matr2=(PASSIVE *) malloc (n*sizeof(PASSIVE));
 if (act_n) act2=(SOURCE *) malloc (act_n*sizeof(SOURCE));
 pascopy(n,matr,matr2);
 actcopy(act_n,act,act2);
 if (flag_out) printf("\n Input: %s (%d,%d), output: %s (%d,%d)\n",
	              act2[ej].reb,act2[ej].v1-48,act2[ej].v2-48,
	              diagui[0].reb,diagui[0].v1-48,diagui[0].v2-48);
 fprintf(out,"%s*\n(",act2[ej].reb);
 // printf("\n%s num\n",diagui[0].reb);
 w0[0]='\0'; w=w0;
 strcat(w++,"1");
 free(act2[ej].reb);
 for (i=ej;i<act_n-1;i++) act2[i]=act2[i+1];
 act2[act_n-1].v1=v2;
 act2[act_n-1].v2=v1;
 act2[act_n-1].v3=v3;
 act2[act_n-1].v4=v4;
 act2[act_n-1].kol=-1;
 act2[act_n-1].reb=strdup(w0);
 // ej=act_n-1;
 if (flag_e)
   { /* простейшие упрощения и вывод подформул */
 //   fputs("\n",out);
    reallyy(&n,matr);
    kontrol();
    seqzz(&n,matr,act_n,act);
    kontrol();
    red2yy(&n,matr,act_n,act);
    kontrol();
    red2zz(&n,matr);
    kontrol();
    reallyz(&n,matr);
    seqzy(&n,matr,act_n,act);
 //   strcat(b++,"\n");
    kontrol();
    fputs(c,out);
   }
 else
   {
    reallf(&n,matr);
    parall(&act_n,act);
   }
 max=ncL; min=0;
// for (st_pln=max;st_pln>=min;st_pln--)
 for (st_pln=min;st_pln<=max;st_pln++)
    {
     fprintf(out,"\n");
     if (st_pln)
       {
	for (i=0;i<st_pln;i++) fprintf(out,"%c*",sp);
       }
     else
       {
	fprintf(out,"1*");
       }
     fprintf(out," (");
     printf("\n%c**%d*\n",sp,st_pln);
     if (n) matr1=(PASSIVE *) malloc (n*sizeof(PASSIVE));
     if (act_n) act1=(SOURCE *) malloc (act_n*sizeof(SOURCE));
     pascopy(n,matr2,matr1);
     actcopy(act_n,act2,act1);
     c[0]='\0';
     b=c;
     dl=leng+strlen(c);
     gggp(n,matr1,act_n,act1,0,st_pln);
     if (n) free(matr1); if (act_n) free(act1);
     if (dl>=leng+strlen(c))
       {
	if (st_pln) fseek(out,-2*st_pln-4,SEEK_CUR);
	else fseek(out,-6,SEEK_CUR);
       }
     else
       {
	fputs(c,out);
	leng+=(b-c); lengl+=(b-c);
	fputs(")",out);
//	if (!strlen(c)) break;
  //  else
		fprintf(out,"+\n");
	}
	}
 fseek(out,-3,SEEK_CUR);
 freematr(n,matr2);  freeact(act_n,act2);
 fputs("\n)\n/\n(",out);
 printf("\n%s den\n",diagui[0].reb);
 allntr999(&act_n,act); /* нейтрализация многомерного НУИ */
 if (diagui[0].reb[0]=='I')
   {
    act=(SOURCE *) realloc (act,(act_n+1)*sizeof(SOURCE));
    nui(v3,v4,&act_n,act);
   }
 // for (st_pln=max;st_pln>=min;st_pln--)
 for (st_pln=min;st_pln<=max;st_pln++)
	 {
     fprintf(out,"\n");
     if (st_pln)
       {
	for (i=0;i<st_pln;i++) fprintf(out,"%c*",sp);
       }
     else fprintf(out,"1*");
     fprintf(out," (");
     printf("\n%c**%d*\n",sp,st_pln);
     if (n) matr1=(PASSIVE *) malloc (n*sizeof(PASSIVE));
     if (act_n) act1=(SOURCE *) malloc (act_n*sizeof(SOURCE));
     pascopy(n,matr,matr1);
     actcopy(act_n,act,act1);
     c[0]='\0';
     b=c;
     dl=leng+strlen(c);
     gggp(n,matr1,act_n,act1,0,st_pln);

     if (dl>=leng+strlen(c))
       {
	if (st_pln) fseek(out,-2*st_pln-4,SEEK_CUR);
	else fseek(out,-6,SEEK_CUR);
       }
     else
       {
	fputs(c,out);
	leng+=(b-c); lengl+=(b-c);
	fputs(")",out);
	//if (!strlen(c)) break;
    //else
		fprintf(out,"+\n");
       }
    }
 if (n) free(matr); if (act_n) free(act);
 fseek(out,-3,SEEK_CUR);
 fprintf(out,"\n);");
}
/* анализ схемы - вывод выражений для искомых токов и напряжений */
 void analys(int n,PASSIVE *matri,int act_n,SOURCE *acti,
	     int dui,PASSIVE *diagui)
{
 int j,k,n0,act0_n;
 char o1,o2;
 PASSIVE *matr0,*matr;
 SOURCE *act0,*act;

 for (k=0;k<dui;k++)
    { /* перебор искомых откликов U и I, k - текущий отклик */
     n0=n;
     act0_n=act_n;
     if (n0) matr0=(PASSIVE *) malloc (n0*sizeof(PASSIVE));
     if (act0_n) act0=(SOURCE *) malloc (act0_n*sizeof(SOURCE));
     pascopy(n0,matri,matr0);
     actcopy(act0_n,acti,act0);
     for (j=0;j<dui;j++)
	{ /* стягивание откликов тока, отличных от текущего отклика I */
	 if (j==k) continue;
	 if (diagui[j].reb[0]=='I')
	   {
	    act0=(SOURCE *) realloc (act0,(act0_n+1)*sizeof(SOURCE));
	    nui(diagui[j].v1,diagui[j].v2,&act0_n,act0);
	   }
	}
     o1=diagui[k].v1; o2=diagui[k].v2; /* узлы текущего отклика */
     for (j=0;j<act0_n;j++)
	{ /* образование управляющей ветви многомерного НУИ */
	 if (act0[j].reb[0]=='E' || act0[j].reb[0]=='J')
	   {
	    act0[j].v3=o2;
	    act0[j].v4=o1;
	   }
	}
     c[0]='\0';
     b=c;
     fprintf(out,"\n");
     if (flag_e)
       { /* простейшие упрощения и вывод подформул */
	fputs("\n",out);
	reallyy(&n0,matr0);
	kontrol();
	seqzz(&n0,matr0,act0_n,act0);
	kontrol();
	red2yy(&n0,matr0,act0_n,act0);
	kontrol();
	red2zz(&n0,matr0);
	kontrol();
	reallyz(&n0,matr0);
	seqzy(&n0,matr0,act0_n,act0);
	strcat(b++,"\n");
	kontrol();
	fputs(c,out);
       }
     else
       {
	reallf(&n0,matr0);
	parall(&act0_n,act0);
       }
     if (n0) matr=(PASSIVE *) malloc (n0*sizeof(PASSIVE));
     if (act0_n) act=(SOURCE *) malloc (act0_n*sizeof(SOURCE));
     pascopy(n0,matr0,matr);
     actcopy(act0_n,act0,act);
     c[0]='\0';
     b=c;
     fprintf(out,"\n%s = \n(\n",diagui[k].reb);
     printf("\n%s num\n",diagui[k].reb);
     lengl=0;
     gggf(n0,matr0,act0_n,act0);
     if (n0) free(matr0); if (act0_n) free(act0);
     fputs(c,out);
     leng+=(b-c);
     lengl+=(b-c);
     if (!(b-c)) fputs(" 0",out),leng++,lengl++;
     fputs("\n)\n/\n(\n",out);
     printf("\n%s den\n",diagui[k].reb);
     allntr999(&act0_n,act); /* нейтрализация многомерного НУИ */
     if (diagui[k].reb[0]=='I')
       {
	act=(SOURCE *) realloc (act,(act0_n+1)*sizeof(SOURCE));
	nui(o1,o2,&act0_n,act);
       }
     c[0]='\0';
     b=c;
     lengl=0;
     gggf(n0,matr,act0_n,act);
     if (n0) free(matr); if (act0_n) free(act);
     fputs(c,out);
     leng+=(b-c);
     lengl+=(b-c);
     if (!(b-c)) fputs(" 0",out),leng++,lengl++;
     fputs("\n);",out);
    }
}
/* анализ схемы - вывод выражений для искомых токов и напряжений */
 void analysp(int n,PASSIVE *matri,int act_n,SOURCE *acti,
	     int dui,PASSIVE *diagui,int ncL)
{
 int i,j,k,n0,act0_n,max,min;
 char o1,o2;
 PASSIVE *matr0,*matr;
 SOURCE *act0,*act;
 long unsigned dl;
 PASSIVE *matr1;
 SOURCE *act1;

 for (k=0;k<dui;k++)
    { /* перебор искомых откликов U и I, k - текущий отклик */
     n0=n;
     act0_n=act_n;
     if (n0) matr0=(PASSIVE *) malloc (n0*sizeof(PASSIVE));
     if (act0_n) act0=(SOURCE *) malloc (act0_n*sizeof(SOURCE));
     pascopy(n0,matri,matr0);
     actcopy(act0_n,acti,act0);
     for (j=0;j<dui;j++)
	{ /* стягивание откликов тока, отличных от текущего отклика I */
	 if (j==k) continue;
	 if (diagui[j].reb[0]=='I')
	   {
	    act0=(SOURCE *) realloc (act0,(act0_n+1)*sizeof(SOURCE));
	    nui(diagui[j].v1,diagui[j].v2,&act0_n,act0);
	   }
	}
     o1=diagui[k].v1; o2=diagui[k].v2; /* узлы текущего отклика */
     for (j=0;j<act0_n;j++)
	{ /* образование управляющей ветви многомерного НУИ */
	 if (act0[j].reb[0]=='E' || act0[j].reb[0]=='J')
	   {
	    act0[j].v3=o2;
	    act0[j].v4=o1;
	   }
	}
     c[0]='\0';
     b=c;
     fprintf(out,"\n");
     if (flag_e)
       { /* простейшие упрощения и вывод подформул */
	fputs("\n",out);
	reallyy(&n0,matr0);
	kontrol();
	seqzz(&n0,matr0,act0_n,act0);
	kontrol();
	red2yy(&n0,matr0,act0_n,act0);
	kontrol();
	red2zz(&n0,matr0);
	kontrol();
	reallyz(&n0,matr0);
	seqzy(&n0,matr0,act0_n,act0);
	strcat(b++,"\n");
	kontrol();
	fputs(c,out);
       }
     else
       {
	reallf(&n0,matr0);
	parall(&act0_n,act0);
       }
     if (n0) matr=(PASSIVE *) malloc (n0*sizeof(PASSIVE));
     if (act0_n) act=(SOURCE *) malloc (act0_n*sizeof(SOURCE));
     pascopy(n0,matr0,matr);
     actcopy(act0_n,act0,act);
     fprintf(out,"\n%s = \n(\n",diagui[k].reb);
     printf("\n%s num\n",diagui[k].reb);


 max=ncL; min=0;
 for (st_pln=max;st_pln>=min;st_pln--)
    {
     fprintf(out,"\n");
     if (st_pln)
       {
	for (i=0;i<st_pln;i++) fprintf(out,"%c*",sp);
       }
     else
       {
	fprintf(out,"1*");
       }
     fprintf(out," (");
     printf("\n%c**%d*\n",sp,st_pln);
     if (n0) matr1=(PASSIVE *) malloc (n0*sizeof(PASSIVE));
     if (act0_n) act1=(SOURCE *) malloc (act0_n*sizeof(SOURCE));
     pascopy(n0,matr0,matr1);
     actcopy(act0_n,act0,act1);
     c[0]='\0';
     b=c;
     dl=leng+strlen(c);
     gggp(n0,matr1,act0_n,act1,0,st_pln);
     if (n0) free(matr1); if (act0_n) free(act1);
     if (dl>=leng+strlen(c))
       {
	if (st_pln) fseek(out,-2*st_pln-4,SEEK_CUR);
	else fseek(out,-6,SEEK_CUR);
       }
     else
       {
	fputs(c,out);
	leng+=(b-c); lengl+=(b-c);
	fputs(")",out);
	fprintf(out,"+\n");
       }
    }
 fseek(out,-3,SEEK_CUR);
 fputs("\n)\n/\n(",out);

 printf("\n%s den\n",diagui[k].reb);
 allntr999(&act0_n,act); /* нейтрализация многомерного НУИ */
 if (diagui[k].reb[0]=='I')
   {
    act=(SOURCE *) realloc (act,(act0_n+1)*sizeof(SOURCE));
    nui(o1,o2,&act0_n,act);
   }
 for (st_pln=max;st_pln>=min;st_pln--)
    {
     fprintf(out,"\n");
     if (st_pln)
       {
	for (i=0;i<st_pln;i++) fprintf(out,"%c*",sp);
       }
     else fprintf(out,"1*");
     fprintf(out," (");
     printf("\n%c**%d*\n",sp,st_pln);
     if (n0) matr1=(PASSIVE *) malloc (n0*sizeof(PASSIVE));
     if (act0_n) act1=(SOURCE *) malloc (act0_n*sizeof(SOURCE));
     pascopy(n0,matr,matr1);
     actcopy(act0_n,act,act1);
     c[0]='\0';
     b=c;
     dl=leng+strlen(c);
     gggp(n0,matr1,act0_n,act1,0,st_pln);

     if (dl>=leng+strlen(c))
       {
	if (st_pln) fseek(out,-2*st_pln-4,SEEK_CUR);
	else fseek(out,-6,SEEK_CUR);
       }
     else
       {
	fputs(c,out);
	leng+=(b-c); lengl+=(b-c);
	fputs(")",out);
	fprintf(out,"+\n");
       }
    }
 if (n) free(matr); if (act_n) free(act);
 fseek(out,-3,SEEK_CUR);
 fprintf(out,"\n);");
}
}
