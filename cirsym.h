#include <stdio.h>
#include <string.h>  
#include <dos.h>
#include <windows.h> 
#include <alloc.h>
#include <conio.h>              
#include <stdlib.h>
#include <math.h>

#define PASSIVE struct pas
#define SOURCE struct ist
#define GRAPH struct graph
#define DIAG struct diag

PASSIVE /* структура двухполюсных ветвей схемы */
{
 char v1; /* первый узел */
 char v2; /* второй узел */
 int kol; /* кратность - количество ветвей, входящих в двухполюсник */
          /* для y-двухполюсника - это количество параллельных ветвей */
          /* для z-двухполюсника - это количество последовательных ветвей */
 char *reb; /* идентификатор двухполюсника (отдельные сопротивление или */
            /* проводимость, сумма проводимостей или сопротивлений)     */
};
SOURCE /* структура управляемых источников схемы */
{
 char v1; /* начальный узел ветви генератора УИ */
 char v2; /* конечный узел ветви генератора УИ */
 char v3; /* начальный узел ветви приемника УИ */
 char v4; /* конечный узел ветви приемника УИ */
 int kol; /* кратность - количество отдельных УИ, образующих данный УИ */
          /* для ИТУН - это количество параллельных ИТУН               */
          /* для ИНУТ - то количество последовательных ИНУТ            */
 char *reb;
};
GRAPH /* структура графа, изоморфного схеме */
{
 char v1; /* первый узел */
 char v2; /* второй узел */
 int ess; /* тип ветви или ее номер (для графа схемы с УИ) */
};
DIAG
{
 int nu; /* номер диагностируемой ветви */
 char *reb; /* идентификатор диагностируемой ветви */
};

void tighten(char send1,char send2,char rec1,char rec2,
	     int *act_n,SOURCE *act);
int ideal(int *act_n,SOURCE *act,
	  char *send1,char *send2,char *rec1,char *rec2,
	  char **t,int *numb,int *sign);
void simplym(int *act_n,SOURCE *act,int *flag);
int smplhng(int act_n,SOURCE *act,int kand);
int hangedm(int act_n,SOURCE *act);
int connecm(int act_n,SOURCE *act);
void parallm(int *act_n,SOURCE *act);
void redumm(int *act_n,SOURCE *act,int *flag);
void trianglm(int *act_n,SOURCE *act,int *flag);
void freactm(int act_n, SOURCE *act);
void fractm(int j,int *act_n,SOURCE *act);
void copyactm(int act_n,SOURCE *act,GRAPH *a2);	  
void quartet(int act_n,SOURCE *act,int *fl);
void nodesm(int act_n,SOURCE *act,int *p,char *str);
void noid_opt(int act_n,SOURCE *act,int *first);
int maxsoum(int *act_n,SOURCE *act,
	    char *c1,char *c2,char **t,int *numb,int *sign);
int degenerm(int *act_n,SOURCE *act,
	  char *s1,char *s2,char **t,int *numb,int *sign);
int nodredmm(int *act_n,SOURCE *act,
	   char  *s1,char *s2, char **t,int *numb,int *sign);	  
void delet(int act1_n,SOURCE *act,int first,SOURCE *act1);
	  
void paslp(int *n,PASSIVE *matr);
void actlp(int *act_n, SOURCE *act,int *fl);
void main(int argc, char* argv[]);

int yesgm(int n,PASSIVE *matr,int act_n,SOURCE *act,
	   int *first,int *fm,char *sm1,char *sm2, int *smf);
int yesgq(int n,PASSIVE *matr,int act_n,SOURCE *act,
	   int *first,int *fm);
int yesgt(int n,PASSIVE *matr,int act_n,SOURCE *act,
	   int *first,int *fm);
int yesm(int act_n,SOURCE *act,int *first);
int yesq(int act_n,SOURCE *act,int *first);
int yest(int act_n,SOURCE *act,int *first);


void parse_command_line(int argc, char* argv[]);

/* Менеджер ввода схемы, проверки корректности, анализа и диагностики */
void input(int *n,PASSIVE *matr,int *act_n,SOURCE *act,
           int *nmd,PASSIVE *matrd,int *actd_n,SOURCE *actd,
           int *dui,PASSIVE *diagui,int *dact,int *diag,
	   int *dva,DIAG *diagva,int *nej,int *da,int *ncL);
void detanm(int act_n,SOURCE *act);
void detan(int n,PASSIVE *matr,int act_n,SOURCE *act);
void detanp(int n,PASSIVE *matr,int act_n,SOURCE *act,int ncL);
void cirfunpa(int n,PASSIVE *matr,int act_n,SOURCE *act,
	      PASSIVE *diagui,int ncL);
void cirfunst(int n,PASSIVE *matr,int act_n,SOURCE *act,
	      PASSIVE *diagui);
void analys(int n,PASSIVE *matr,int act_n,SOURCE *act,
	    int dui,PASSIVE *diagui);
void analysp(int n,PASSIVE *matri,int act_n,SOURCE *acti,
	     int dui,PASSIVE *diagui,int ncL);

/* Рекурсивное динамическое разложение схемного определителя */
void gggfm(int act_n,SOURCE *act);
void gggf(int n,PASSIVE *matr,int act_n,SOURCE *act);
void multipl(int *n,PASSIVE *matr,int *act_n,SOURCE *act,int *fl);
void gggp(int n,PASSIVE *matr,int act_n,SOURCE *act,int pln,int st_pln);
void multiplp(int *n,PASSIVE *matr,int *act_n,SOURCE *act,
              int *fl,int pln,int st_pln);
void newtral(int *n,PASSIVE *matr,int *act_n,SOURCE *act,int *fl);
void allntr999(int *act_n,SOURCE *act);
void autontr999(int k,int *act_n,SOURCE *act,int *fl);
void ntr999(int *act_n,SOURCE *act);
int is999(int j,int act_n,SOURCE *act);
void one999(int act_n,SOURCE *act);

/* Функции и подпрограммы для проверки вырождения и упрощения схемы */
void simply(int *n,PASSIVE *matr,int *act_n,SOURCE *act,int *flag);
void triangl(int *n,PASSIVE *matr,int *act_n,SOURCE *act,int *flag);
void redu(int *n,PASSIVE *matr,int *act_n,SOURCE *act);
void reallf(int *n,PASSIVE *a1);
void parall(int *act_n,SOURCE *act);
void seqzz0(int *n,PASSIVE *a1,int act_n,SOURCE *act);
int pasloop(int *n,PASSIVE *matr,
	     char *s1,char *s2,char **t,int *numb,int *sign);
int redur(int *n,PASSIVE *matr,int act_n,SOURCE *act,
	  char *ss1,char *ss2,char **t,int *numb,int *sign);
int degener(int *act_n,SOURCE *act,
	  char *s1,char *s2,char **t,int *numb,int *sign);
int pashang(int *n,PASSIVE *matr,int act_n,SOURCE *act,
	    char *s1,char *s2,char **t,int *numb,int *sign);
int nodred(int n,PASSIVE *matr,int *act_n,SOURCE *act,
	   char  *s1,char *s2, char **t,int *numb,int *sign);
int nodreda(int n,int *act_n,SOURCE *act,
	   char  *s1,char *s2, char **t,int *numb,int *sign,int *flm);
void ejuihang(int *n,PASSIVE *matr,int *act_n,SOURCE *act,int *fl);
void actloop(int *act_n, SOURCE *act,int *fl);
void redparm(int *act_n,SOURCE *act, int *fl);
void redpar(int *act_n,SOURCE *act, int *fl);
void redseq(int n,PASSIVE *matr,int *act_n,SOURCE *act,int *fl);
int hanged(int n,PASSIVE *matr,int act_n,SOURCE *act);
int connec(int n,PASSIVE *matr,int act_n,SOURCE *act);
void hangr(int *n,PASSIVE *matr,int *act_n,SOURCE *act,int *fl);

/* Функции вывода подформул */
void red2yy(int *n,PASSIVE *a1,int act_n,SOURCE *act);
void red2zz(int *n,PASSIVE *a1);
void seqzz(int *n,PASSIVE *a1,int act_n,SOURCE *act);
void reallyy(int *n,PASSIVE *a1);
void parallm(int *act_n,SOURCE *act);
void reallyz(int *n,PASSIVE *a1);
void seqzy(int *n,PASSIVE *a1,int act_n,SOURCE *act);

/* Вспомогательные функции и подпрограммы */
void extra(int el_n,char *el,int be,int *ee,char *sr,int *k);
void ster(unsigned long k);
void verstr(int n,GRAPH *matr,int *p,char *str);
void ver_n(int n,GRAPH *matr,int *p);
void nodestr(int n,PASSIVE *matr,int act_n,SOURCE *act,int *p,char *str);
void node_n(int n,PASSIVE *matr,int act_n,SOURCE *act,int *p);
void nodes(int act_n,SOURCE *act,int *p,char *str);
void node(int act_n,SOURCE *act,int *p);
void printa(char *t,int numb,int sign);
void printi(char **t,char *a,int numb);
int noideam(int act_n,SOURCE *act,int *first);
int noideal(int act_n,SOURCE *act,int *first);
void copypas(int n,PASSIVE *matr,GRAPH *a1);
void actcopy(int act_n,SOURCE *act,SOURCE *act1);
void pascopy(int n,PASSIVE *matr,PASSIVE *matr1);
void copyact(int n,int act_n,SOURCE *act,int *na,GRAPH *a2);
void uniact(char s1,char s2,int *act_n,SOURCE *act);
void uniactp(char s1,char s2,int *act_n,int act0_n,SOURCE *act);
void uniactp2(char s1,char s2,
              int *act_n,int act0_n,int *act3_n,SOURCE *act);
void unipas(char s1,char s2,int *n,PASSIVE *matr);
void unipasp(char s1,char s2,int *n,int n0,PASSIVE *matr);
void unipasp2(char s1,char s2,int *n,int n0,int *n3,PASSIVE *matr);
void kontrol(void);
int minimum(int pow1,int pow2);
int maximum(int activ1,int activ2);
void hangpas(int n,PASSIVE *matr,int act_n,SOURCE *act,int *first);
void hangpas_p(int n,PASSIVE *matr,int act_n,SOURCE *act,int *first);
void choiceg(int n,PASSIVE *matr,int act_n,SOURCE *act,int *first);
void choiceg_p(int n,PASSIVE *matr,int act_n,SOURCE *act,int *first);
void choicer(int n,PASSIVE *matr,int act_n,SOURCE *act,int *first);
void choicer_p(int n,PASSIVE *matr,int act_n,SOURCE *act,int *first);
void nui(char s1,char s2,int *act_n,SOURCE *act);
void nuiall(char s1,char s2,char s3,char s4,int *act_n,SOURCE *act);
void nuintr(char s1,char s2,int k,SOURCE *act);
void actall(char s1,char s2,char s3,char s4,
	    int kol,char *str,int *act_n,SOURCE *act);
void matrall(char s1,char s2,int kol,char *str,int *n,PASSIVE *matr);
void freematr(int n, PASSIVE *matr);
void freeact(int act_n, SOURCE *act);
void frematr(int n, PASSIVE *matr);
void freact(int act_n, SOURCE *act);
void fract(int j,int *act_n,SOURCE *act);
void frmatr(int j,int *n,PASSIVE *matr);
int yesd(int act_n,SOURCE *act);
int yesc(int n,PASSIVE *matr);
int yesL(int n,PASSIVE *matr);
int yescL(int n,PASSIVE *matr);
int yescL_n(int n,PASSIVE *matr);
int yesg(int n,PASSIVE *matr);
int yesr(int n,PASSIVE *matr);
int yesrL(int n,PASSIVE *matr);
int yess(int act_n,SOURCE *act);
int yesn(int act_n,SOURCE *act);
// float fabs(float x);
void bond1f(int n,GRAPH *a1,char s1,char s2,int *k);

/* Разложение схемного определителя по частям - бисекция схемы */
void binvec(int next, char *sn, char *on, 
            char nus[16][5],int *s);
void bond1(int n,GRAPH *a1,char s1,char s2,int *k);
void bond2(int n,GRAPH *a1,char s1,char s2,char s3,int *k);
void bond3(int n,GRAPH *a1,char s1,char s2,char s3,char s4,int *k);
void bond4(int n,GRAPH *a1,
	   char s1,char s2,char s3,char s4,char s5,int *k);
void bond5(int n,GRAPH *a1,
	   char s1,char s2,char s3,char s4,char s5,char s6,int *k);
int indep1(int nc,GRAPH *ac);
int indep2(int nc,GRAPH *ac,char *ac0,char s1,char s2,int *k);
int indep3(int nc,GRAPH *ac,
	   char *ac0,char s1,char s2,char s3,int *k);
int indep4(int nc,GRAPH *ac,
	   char *ac0,char s1,char s2,char s3,char s4,int *k);
int indep5(int nc,GRAPH *ac,
	    char *ac0,char s1,char s2,char s3,char s4,char s5,int *k);
void fortra(int num,char *ac,SOURCE *act,
	          SOURCE *act1,int key);
void fortrans(int num,char *ac,int n,PASSIVE *matr,SOURCE *act,
	          PASSIVE *a1,SOURCE *act1,int key);
int n_cL(int n,PASSIVE *matr,char *ac,int key);
void transf(int act_n,SOURCE *act,int *num,GRAPH *ac);
void transfor(int n,PASSIVE *matr,
	      int act_n,SOURCE *act,int *num,GRAPH *ac);
void forte(int num,char *ac,int *act1_n,int key);		  
void fortest(int num,int n,char *ac,int *n1,int *act1_n,int key);
void form1(int num,GRAPH *ac,int act_n,SOURCE *act);
void form2(int num,GRAPH *ac0,char s1,char s2,
	       int act_n,SOURCE *act);
void form1p(int num,char *ac,int n,PASSIVE *matr,
	   	   int act_n,SOURCE *act,int pln,int st_pln);
void form2p(int num,char *ac,char s1,char s2,
	 int n,PASSIVE *matr,int act_n,SOURCE *act,int pln,int st_pln);
void formm(int num,char *ac,int next,char *ext,
	       int act_n,SOURCE *act);
void form(int num,char *ac,int next,char *ext,
	   int n,PASSIVE *matr,int act_n,SOURCE *act);
void formp(int num,char *ac,int next,char *ext,
	   int n,PASSIVE *matr,int act_n,SOURCE *act,int pln,int st_pln);
void bisec1m(int act_n,SOURCE *act,int *flag);
void bisec1(int n,PASSIVE *matr,int act_n,SOURCE *act,int *flag);
void bisec1p(int n,PASSIVE *matr,int act_n,SOURCE *act,
             int *flag,int pln,int st_pln);
void bisec2m(int act_n,SOURCE *act,int *flag,float range);		 
void bisec2(int n,PASSIVE *matr,int act_n,SOURCE *act,int *flag,float range);
void bisec2p(int n,PASSIVE *matr,int act_n,SOURCE *act,
             int *flag,float range,int pln,int st_pln);
void bisec3m(int act_n,SOURCE *act,int *flag,float range);
void bisec3(int n,PASSIVE *matr,int act_n,SOURCE *act,int *flag,float range);
void bisec3p(int n,PASSIVE *matr,int act_n,SOURCE *act,
             int *flag,float range,int pln,int st_pln);
void bisec4m(int act_n,SOURCE *act,int *flag,float range);
void bisec4(int n,PASSIVE *matr,int act_n,SOURCE *act,int *flag,float range);
void bisec4p(int n,PASSIVE *matr,int act_n,SOURCE *act,
             int *flag,float range,int pln,int st_pln);
void bisec5m(int act_n,SOURCE *act,int *flag,float range);			 
void bisec5(int n,PASSIVE *matr,int act_n,SOURCE *act,int *flag,float range);
void bisec5p(int n,PASSIVE *matr,int act_n,SOURCE *act,
             int *flag,float range,int pln,int st_pln);

void mirnul(char s1,char s2,char s3,char s4,int *act_n,SOURCE *act,char type);
void mirall(char s1,char s2,char s3,char s4,int *p,char *str,
	    int *n,PASSIVE *matr,int *act_n,SOURCE *act);
void qirall(char s1,char s2,char s3,char s4,int *p,char *str,
	    int *n,PASSIVE *matr,int *act_n,SOURCE *act);
void tirall(char s1,char s2,char s3,char s4,int *p,char *str,
	    int *n,PASSIVE *matr,int *act_n,SOURCE *act);
void mloop(int *act_n, SOURCE *act,int *fl);
int yeslm(int n,PASSIVE *matr,int act_n,SOURCE *act,
	    int *first,int *fm,char *sm1,char *sm2, int *smf);
int yeslq(int n,PASSIVE *matr,int act_n,SOURCE *act,
	    int *first,int *fm,char *sm1,char *sm2, int *smf);
int yeslt(int n,PASSIVE *matr,int act_n,SOURCE *act,
	    int *first,int *fm,char *sm1,char *sm2, int *smf);
int mtria(int *n,PASSIVE *matr,int act_n,SOURCE *act,int *fl);
void mirror(int act_n,SOURCE *act);
int mqtloop(int act_n,SOURCE *act,
	     char *s1,char *s2,char **t,int *numb,int *sign);
int nodredm(int n,PASSIVE *matr,int *act_n,SOURCE *act,
	    char  *s1,char *s2, char **t,int *numb,int *sign,int *flm);
int nodredg(int n,PASSIVE *matr,int *act_n,SOURCE *act,
	    char  *s1,char *s2, char **t,int *numb,int *sign,int *flm);
int yesmqt(int act_n,SOURCE *act);
void redum(int act_n,SOURCE *act);