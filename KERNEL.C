/* Рекурсивное динамическое разложение схемного определителя */

#include "cirsym.h"

extern FILE *out,*inpa,*set;

extern char *b,*c,sp;
extern unsigned long leng,lengl,lepr;
extern int flag_fi3,flag_fi4,flag_fi5,
       flag_fsn1,flag_fsn2,flag_nul,extract,
       flag_e,cirf,noequ,flag_pln,nuz,flag_sp,flag_cL,
       flag_mir,flag_opt,flag_mat,flag_two,flag_thr,
	   flag_2,flag_3,flag_4,flag_5,fl_g,fl_mGc,fl_mGm;
extern float range2,range3,range4,range5,
             range2p,range3p,range4p,range5p;

/* Функция проверки существования многомерного НУИ */
 int is999(int j,int act_n,SOURCE *act)
{
 int i;

 for (i=0;i<act_n;i++)
    {
     if (i==j) continue;
     if (abs(act[i].kol) == 999) return(1);
    }
 return(0);
}
/* Нейтрализация многомерного НУИ */
 void allntr999(int *act_n,SOURCE *act)
{
 int i,act0_n;

 for (i=0;i<*act_n;i++)
    {
     act0_n=*act_n;
     if (act[i].kol==999)
       {
        if (act[i].reb[0] == 'E' || act[i].reb[0] == 'U')
          nuintr(act[i].v1,act[i].v2,i,act);
        else
          if (act[i].reb[0] == 'J' || act[i].reb[0] == 'I')
            {
	     fract(i,act_n,act);
	     if (*act_n != act0_n) i--;
	    }
       }
    }
}
/* Нейтрализация k-го элемента многомерного НУИ */
 void autontr999(int k,int *act_n,SOURCE *act,int *fl)
{
 *fl=2;
 if (act[k].kol < 0) return;
 if (act[k].kol != 999) {act[k].kol=0; *fl=3; return;}
 if (is999(k,*act_n,act))
   {
    *fl=1;
    if (act[k].reb[0] == 'E' || act[k].reb[0] == 'U')
      nuintr(act[k].v1,act[k].v2,k,act);
    else
    if (act[k].reb[0] == 'J' || act[k].reb[0] == 'I')
      fract(k,act_n,act);
   }
 one999(*act_n,act);
}
/* Нейтрализация оставшихся элементов многомерного НУИ */
 void ntr999(int *act_n,SOURCE *act)
{
 int i;

 for (i=0;i<*act_n;i++)
    {
     if (act[i].kol != 999) continue;
     if (act[i].reb[0] == 'E' || act[i].reb[0] == 'U')
       nuintr(act[i].v1,act[i].v2,i,act);
     else
     if (act[i].reb[0] == 'J' || act[i].reb[0] == 'I')
       fract(i,act_n,act);
    }
}
/* Преобразование многомерного НУИ (с одним элементом) в обычный НУИ*/
 void one999(int act_n,SOURCE *act)
{
 int i;

 for (i=0;i<act_n;i++)
    if (act[i].kol == 999) {act[i].kol=-999; break;}
}
/* Нейтрализация УИ (act[j].kol = 0) с проверкой на вырождение схемы */
 void newtral(int *n,PASSIVE *matr,int *act_n,SOURCE *act,int *fl)
{
 int j,p1,p2;
 char s1,s2,s3,s4,type;

nwtrl:
 for (j=0;j<*act_n;j++)
    {
     if (act[j].kol) continue;
     node_n(*n,matr,*act_n,act,&p1);
     type=act[j].reb[0];
     s1=act[j].v1; s2=act[j].v2; s3=act[j].v3; s4=act[j].v4;
     if (type == 'K')
       {
	nuintr(s1,s2,j,act);
	node_n(*n,matr,*act_n,act,&p2);
	if (p1 != p2) {*fl=2; goto ret;}
	else break;
       }
     else
     if (type == 'B')
       {
	nuintr(s3,s4,j,act);
	node_n(*n,matr,*act_n,act,&p2);
	if (p1 != p2) {*fl=2; goto ret;}
	else break;
       }
     else
       if (type == 'H')
	 {
                    if (s3=='0') {s3=s4; s4='0';}
	  nuintr(s1,s2,j,act);
	  unipas(s3,s4,n,matr);
	  uniact(s3,s4,act_n,act);
	 }
       else
       if (type == 'G')
	 {
	  fract(j,act_n,act);
	  node_n(*n,matr,*act_n,act,&p2);
	  if (p1 != p2) {*fl=2; goto ret;}
	  else break;
	 }
    }
 simply(n,matr,act_n,act,fl);
 if (*fl == 3) goto nwtrl;
 else return; /* *fl == 0 - det!=0 или *fl == 2 - det=0 */
ret:
 frematr(*n,matr); freact(*act_n,act);
}
/* Выделение элементов с максимальным участием (общих множителей) */
/*                  в полиномиальном определителе                 */
void multiplp(int *n,PASSIVE *matr,
              int *act_n,SOURCE *act,int *fl,int pln,int st_pln)
{
 int p2,sign,signs,numb,tnum,tn,sens,k=0,flm;
 char s1,s2,*t,*m,m0[200];
 unsigned long dl1,dl2;

 node_n(*n,matr,*act_n,act,&p2);
 if (!p2 && *(b-1)=='*' && pln==st_pln)
   {
    strcat(b++,"1");
    *fl=1;
    return;
   }
 else if (!p2 && pln != st_pln)
	{
	 *fl=2;
	 return;
	}
 newtral(n,matr,act_n,act,fl);
 if (*fl == 2) return;
 /* выявление общих множителей */
 flm=0;
 if (
 // !flag_nul && mqtloop(*act_n,act,&s1,&s2,&t,&numb,&sign)|| 
   pasloop(n,matr,&s1,&s2,&t,&numb,&sign) || 
   pashang(n,matr,*act_n,act,&s1,&s2,&t,&numb,&sign)
 || redur(n,matr,*act_n,act,&s1,&s2,&t,&numb,&sign)
 || degener(act_n,act,&s1,&s2,&t,&numb,&sign)
 || nodred(*n,matr,act_n,act,&s1,&s2,&t,&numb,&sign)
 || nodreda(*n,act_n,act,&s1,&s2,&t,&numb,&sign,&flm)
    )
    { /* t - очередной множитель; numb, sign - кратность и знак множителя */
    signs=0; /* знак строки сомножителей */
    tnum=0;
    m0[0]='\0'; m=m0;
rep:
    if (flm==1) {*fl=2; return;}
    if (sign && signs || !sign && !signs) signs=0; else signs=1;
    if (t[0]=='1') goto nopr;
    tn=tnum;
    if (abs(numb) == 1) tnum+=strlen(t)+1;
    else tnum+=strlen(t)+3;
    m0[tn]='\0'; m=m0+tn;
    if (   !strlen(m0) && !*n && !*act_n  && !signs
	&& ((*(b-1) == '(' || *(b-1) == '+' || *(b-1) == '-' /* || *(b-1) == '*' */)
	&& strlen(c)  || !lengl ))
      {
       numb=1;
      }
    if (strchr(t,'c')!=NULL || strchr(t,'L')!=NULL) k++;
    printi(&m,t,numb); /* вывод в строку t очередного сомножителя */
nopr:
    free(t);
    if (*n || *act_n)
      {
       if (s1 != ' ' && s1 != s2)
	 {
	  if (s1=='0') {s1=s2; s2='0';}
	  unipas(s1,s2,n,matr);
	  uniact(s1,s2,act_n,act);
	 }
      // mirror(*act_n,act);
       newtral(n,matr,act_n,act,fl);
       if (*fl == 2) return;
       node_n(*n,matr,*act_n,act,&p2);
       if (!p2 || p2==1 && !n) goto cont;
/* повторный проход */
       flm=0;
       if  (
  // !flag_nul && mqtloop(*act_n,act,&s1,&s2,&t,&numb,&sign) || 
   pasloop(n,matr,&s1,&s2,&t,&numb,&sign) || 
   pashang(n,matr,*act_n,act,&s1,&s2,&t,&numb,&sign)
 || redur(n,matr,*act_n,act,&s1,&s2,&t,&numb,&sign)
 || degener(act_n,act,&s1,&s2,&t,&numb,&sign)	
 || nodred(*n,matr,act_n,act,&s1,&s2,&t,&numb,&sign)
 || nodreda(*n,act_n,act,&s1,&s2,&t,&numb,&sign,&flm)	  
	  ) goto rep;
cont:
       dl2=leng+strlen(c);
       if (!strlen(m0)) strcat(m++,"1"), m[0]='\0'; /* множитель-единица */
       else {m--; m[0]='\0';} /* удаление лишнего знака "*" */
       if (k + pln > st_pln) {*fl=2; return;}
       else
       printa(m0,1,signs); /* вывод строки сомножителей в выходную строку с */
       if (strlen(m0)) strcat(b++,"*");
       if (   (*(b-3) == '(' || *(b-3) == '+'
	   || *(b-3) == '-' || strlen(c) == 2)
	   && *(b-2) == '1' && *(b-1) == '*')
	 {
	  b-=2; b[0]='\0'; /* удаление избыточного множителя "1*" */
	 }
       if (*(b-1)=='-') sens=1; else sens=0;
       if (*n + *act_n > 0 && (m0[0] != '1' || sens)) strcat(b++,"(");
       dl1=leng+strlen(c);
       gggp(*n,matr,*act_n,act,pln+k,st_pln);
       if (dl1>=leng+strlen(c))
	 {
	  ster(dl2); /* определитель схемы равен нулю */
	  *fl=2;
	  return;
	 }
       else
	 {
	  if (*n + *act_n > 0 && (m0[0] != '1' || sens)) strcat(b++,")");
	 }
      }
    else
      { /* определитель схемы является произведением сомножителей */
       if (k + pln != st_pln) {*fl=2; return;}
       if (!strlen(m0)) {strcat(m++,"1"); m[0]='\0';}
       else {m--; m[0]='\0';} /* удаление лишнего знака "*" */
       printa(m0,1,signs); /* вывод строки сомножителей в выходную строку с */
      }
    *fl=1; /* определитель раскрыт */
    return;
   }
 *fl=3; /* нет общих множителей */
 return;
}
/* Генератор формулы полиномиального схемного определителя */
 void gggp(int n,PASSIVE *matr,int act_n,SOURCE *act,int pln,int st_pln)
{
 int fl,p1,p2,first,fm,smf,n1,act1_n,extract0,q,i,k;
 char s1,s2,sm1,sm2;
 unsigned long dl,dl1,dl2;
 PASSIVE *a1;
 SOURCE *act1;

/* коррекция знака подвыражения */

if (yescL_n(n,matr) < st_pln-pln) {fl = 2; return;} 
beg:
 if (*(b-1)=='-')
   {
    b--; b[0]='\0';
    strcat(b++,"+");
   }   
 // выделение общих множителей
 // if (!flag_nul) mirror(act_n,act);
 multiplp(&n,matr,&act_n,act,&fl,pln,st_pln);
 if (fl == 1 || fl == 2) return; /* det схемы != 0 || det = 0 */
 node_n(n,matr,act_n,act,&p1);
 // степень полинома превышена - схема удаляется
 if (pln > st_pln) {frematr(n,matr); freact(act_n,act);}
 else
 if (pln==st_pln) {
	 if (yescL(n,matr))
       {  // удаление емкостей и индуктивностей
          // после достижения степени полинома
    i=0;
    while(i<n)
	{
	 if (matr[i].reb[0]=='c')
	   {frmatr(i,&n,matr); i--;}
	 else
	 if (matr[i].reb[0]=='L') {
		 if (!act_n) act=(SOURCE *) malloc (sizeof(SOURCE));
		 else act=(SOURCE *) realloc (act,(act_n+1)*sizeof(SOURCE));
    	 nui(matr[i].v1,matr[i].v2,&act_n,act);
	 	 frmatr(i,&n,matr); i--;
		}
	i++;
	}
	 node_n(n,matr,act_n,act,&p2);
	 if (   p1 != p2 || !hanged(n,matr,act_n,act)
     || !connec(n,matr,act_n,act))
	 {
	  frematr(n,matr); freact(act_n,act); return;}
    goto beg;
   }
 // на счет резистивной схемы
 /* деление на две подсхемы по одному узлу */
    bisec1(n,matr,act_n,act,&fl);
    if (fl) return;
 /* деление на две подсхемы по двум - пяти узлам */
    if (flag_fsn2 && p1 > flag_fsn2)
      {
       bisec2(n,matr,act_n,act,&fl,range2);
       if (fl) return;
      }
	if (flag_fi3 && p1 > flag_fi3)
      {
       bisec3(n,matr,act_n,act,&fl,range3);
       if (fl) return;
      }
	if (flag_fi4  &&  p1 > flag_fi4)
      {
       bisec4(n,matr,act_n,act,&fl,range4);
       if (fl) return;
      }
    if (flag_fi5 && p1 > flag_fi5)
      {
       bisec5(n,matr,act_n,act,&fl,range5);
       if (fl) return;
      }
 if (yesg(n,matr))
   {  /* выделение активной проводимости */
    hangpas(n,matr,act_n,act,&first);
    if (first==-1)
      choiceg(n,matr,act_n,act,&first);
    s1=matr[first].v1; s2=matr[first].v2;
    dl=leng+strlen(c);
    dl2=leng+strlen(c);
    printa(matr[first].reb,matr[first].kol,0);
    frmatr(first,&n,matr);
    strcat(b++,"*");
    n1=n; act1_n=act_n;
    if (n1)
      {
       a1=(PASSIVE *) malloc(n1*sizeof(PASSIVE));
       pascopy(n1,matr,a1);
      }
    act1=(SOURCE *) malloc ((act1_n+1)*sizeof(SOURCE));
    if (act1_n) actcopy(act1_n,act,act1);
    nui(s1,s2,&act1_n,act1);
    if (n1+act1_n > 1) strcat(b++,"(");
    dl1=leng+strlen(c);
    gggf(n1,a1,act1_n,act1);
    if (n1) free(a1); if (act1_n) free(act1);
    if (dl1>=leng+strlen(c)) ster(dl2);
    else
      {
	   if (*(b-3) == '*' && *(b-2) == '(' && *(b-1) == '1')
         {b-=3; b[0]='\0';}
       else	
	   if (n1+act1_n > 1) strcat(b++,")");
       dl=leng+strlen(c);
       strcat(b++,"+");
      }
two3:
    node_n(n,matr,act_n,act,&p2);
    dl1=leng+strlen(c);
    if (p2 == p1)
      {
       gggf(n,matr,act_n,act);
       if (dl1 >= leng+strlen(c)) ster(dl);
      }
    else
      {
       ster(dl);
       frematr(n,matr); freact(act_n,act);
      }
   }
 else
 if (yesr(n,matr))
   { /* выделение активного сопротивления */
    choicer(n,matr,act_n,act,&first);
    dl=leng+strlen(c);
    dl2=leng+strlen(c);
    s1=matr[first].v1; s2=matr[first].v2;
    printa(matr[first].reb,matr[first].kol,0);
    frmatr(first,&n,matr);
    strcat(b++,"*");
    n1=n; act1_n=act_n;
    node_n(n1,matr,act1_n,act,&p2);
    if ((n1 || act1_n) && p2 == p1)
      {
       if (n1)
	 {
	  a1=(PASSIVE *) malloc(n1*sizeof(PASSIVE));
	  pascopy(n1,matr,a1);
	 }
       if (act1_n)
	 {
	  act1=(SOURCE *) malloc (act1_n*sizeof(SOURCE));
	  actcopy(act1_n,act,act1);
	 }
       if (n1+act1_n > 1) strcat(b++,"(");
       dl1=leng+strlen(c);
       gggf(n1,a1,act1_n,act1);
       if (n1) free(a1); if (act1_n) free(act1);
       if (dl1>=leng+strlen(c)) ster(dl2);
       else
	 {
	  if (*(b-3) == '*' && *(b-2) == '(' && *(b-1) == '1')
         {b-=3; b[0]='\0';}
      else	
	  if (n1+act1_n > 1) strcat(b++,")");
	  dl=leng+strlen(c);
	  strcat(b++,"+");
	 }
      }
    else ster(dl2);
two4:
    act1_n=act_n;
    act1=(SOURCE*) malloc((act1_n+1)*sizeof(SOURCE));
    if (act1_n) actcopy(act1_n,act,act1);
    nui(s1,s2,&act1_n,act1);
    dl1=leng+strlen(c);
    gggf(n,matr,act1_n,act1);
    if (dl1 >= leng+strlen(c)) ster(dl);
    free(act1);
    freact(act_n,act);
   }
 else
 if (act_n)
   { /* выделение управляемого источника */
    char type;

	noideal(act_n,act,&first);
	if (first == -1) {
       ster(dl);
       freematr(n,a1); freeact(act1_n,act1);
	   return;}
    type=act[first].reb[0];
    if (n)
      {
       a1=(PASSIVE *) malloc(n*sizeof(PASSIVE));
       pascopy(n,matr,a1);
      }
    act1_n=act_n;
    act1=(SOURCE *) malloc (act1_n*sizeof(SOURCE));
    actcopy(act1_n,act,act1);
    act[first].kol=-act[first].kol;
    if (act[first].kol==-999) ntr999(&act_n,act);
    dl=leng+strlen(c);
    if (type == 'K' || type == 'B')
      {
       fl=0;
       if (!n && act_n==1) {strcat(b++,"("); fl=1;}
      }
    dl1=leng+strlen(c);
    gggf(n,matr,act_n,act);
    if (dl1 >= leng+strlen(c)) ster(dl);
    dl=leng+strlen(c);
    if (dl1 < leng+strlen(c)) strcat(b++,"+");
    if (type == 'K' || type == 'B')
      {
       if (fl)
	 {
	  strcat(b,"1)"); b+=2;
	  freematr(n,a1); freeact(act1_n,act1);
	  return;
	 }
      }
    autontr999(first,&act1_n,act1,&fl);
    if (fl == 1 || fl == 3)
      {
       dl1=leng+strlen(c);
       gggf(n,a1,act1_n,act1);
       if (n) free(a1); if (act1_n) free(act1);
       if (dl1>=leng+strlen(c)) ster(dl);
      }
    else
      {
       ster(dl);
       freematr(n,a1); freeact(act1_n,act1);
      }
   }   /* end active */
  }
 else     
 if (pln < st_pln) 
 {
  if (yescL(n,matr)) {
  // на счет реактивно-резистивной схемы
  if (p1 > 2) {
  bisec1p(n,matr,act_n,act,&fl,pln,st_pln);
  if (fl) return;}
  if (flag_2 && p1 > flag_2) {
  bisec2p(n,matr,act_n,act,&fl,range2p,pln,st_pln);
  if (fl) return;}
  if (flag_3 && p1 > flag_3) {
  bisec3p(n,matr,act_n,act,&fl,range3p,pln,st_pln);
  if (fl) return;}
  if (flag_4  &&  p1 > flag_4) {
  bisec4p(n,matr,act_n,act,&fl,range4p,pln,st_pln);
  if (fl) return;}
  if (flag_5 && p1 > flag_5) {
  bisec5p(n,matr,act_n,act,&fl,range5p,pln,st_pln);
  if (fl) return;}
  if (!fl_g && yesg(n,matr))
   {  
    hangpas(n,matr,act_n,act,&first);
    if (first==-1)
      choiceg(n,matr,act_n,act,&first);
    s1=matr[first].v1; s2=matr[first].v2;
    dl=leng+strlen(c);
    dl2=leng+strlen(c);
    printa(matr[first].reb,matr[first].kol,0);
    frmatr(first,&n,matr);
    strcat(b++,"*");
    n1=n; act1_n=act_n;
    if (n1)
      {
       a1=(PASSIVE *) malloc(n1*sizeof(PASSIVE));
       pascopy(n1,matr,a1);
      }
    act1=(SOURCE *) malloc ((act1_n+1)*sizeof(SOURCE));
    if (act1_n) actcopy(act1_n,act,act1);
    nui(s1,s2,&act1_n,act1);
    if (n1+act1_n > 1) strcat(b++,"(");
    dl1=leng+strlen(c);
	gggp(n1,a1,act1_n,act1,pln,st_pln);
    if (n1) free(a1); if (act1_n) free(act1);
    if (dl1>=leng+strlen(c)) ster(dl2);
    else
      {
	   if (*(b-3) == '*' && *(b-2) == '(' && *(b-1) == '1')
         {b-=3; b[0]='\0';}
       else	
       if (n1+act1_n > 1) strcat(b++,")");
       dl=leng+strlen(c);
       strcat(b++,"+");
      }
two33:
    node_n(n,matr,act_n,act,&p2);
    dl1=leng+strlen(c);
    if (p2 == p1)
      {
	   gggp(n,matr,act_n,act,pln,st_pln);
       if (dl1 >= leng+strlen(c)) ster(dl);
      }
    else
      {
       ster(dl);
       frematr(n,matr); freact(act_n,act);
      }
   }
 else   
  if (yesc(n,matr))
   { /* выделение емкости */
    hangpas_p(n,matr,act_n,act,&first);
	if (first==-1) choiceg_p(n,matr,act_n,act,&first);
    s1=matr[first].v1; s2=matr[first].v2;
    dl=leng+strlen(c);
    dl2=leng+strlen(c);
    printa(matr[first].reb,matr[first].kol,0);
    frmatr(first,&n,matr);
    strcat(b++,"*");
    n1=n; act1_n=act_n;
    if (n1)
      {
       a1=(PASSIVE *) malloc(n1*sizeof(PASSIVE));
       pascopy(n1,matr,a1);
      }
    act1=(SOURCE *) malloc ((act1_n+1)*sizeof(SOURCE));
    if (act1_n) actcopy(act1_n,act,act1);
    nui(s1,s2,&act1_n,act1);
    if (n1+act1_n > 1) strcat(b++,"(");
    dl1=leng+strlen(c);
    gggp(n1,a1,act1_n,act1,pln+1,st_pln);
    if (n1) free(a1); if (act1_n) free(act1);
    if (dl1>=leng+strlen(c)) ster(dl2);
    else
      {
       if (*(b-3) == '*' && *(b-2) == '(' && *(b-1) == '1')
         {b-=3; b[0]='\0';}
       else		   
       if (n1+act1_n > 1) strcat(b++,")");
       dl=leng+strlen(c);
       strcat(b++,"+");
      }
two1:
    node_n(n,matr,act_n,act,&p2);
    dl1=leng+strlen(c);
    if (p2 == p1)
      {
       gggp(n,matr,act_n,act,pln,st_pln);
       if (dl1 >= leng+strlen(c)) ster(dl);
      }
    else
      {
       ster(dl);
       frematr(n,matr); freact(act_n,act);
      }
   } /* end g-extraction */
 else
 if (yesL(n,matr))
   { /* выделение индуктивности */
    choicer_p(n,matr,act_n,act,&first);
    dl=leng+strlen(c);
    dl2=leng+strlen(c);
    s1=matr[first].v1; s2=matr[first].v2;
    printa(matr[first].reb,matr[first].kol,0);
    frmatr(first,&n,matr);
    strcat(b++,"*");
    n1=n; act1_n=act_n;
    node_n(n1,matr,act1_n,act,&p2);
    if ((n1 || act1_n) && p2 == p1)
      {
       if (n1)	 {
	  a1=(PASSIVE *) malloc(n1*sizeof(PASSIVE));
	  pascopy(n1,matr,a1); }
       if (act1_n)
	 {
	  act1=(SOURCE *) malloc (act1_n*sizeof(SOURCE));
	  actcopy(act1_n,act,act1);
	 }
       if (n1+act1_n > 1) strcat(b++,"(");
       dl1=leng+strlen(c);
       gggp(n1,a1,act1_n,act1,pln+1,st_pln);
       if (n1) free(a1); if (act1_n) free(act1);
       if (dl1>=leng+strlen(c)) ster(dl2);
       else
	 {
	  if (*(b-3) == '*' && *(b-2) == '(' && *(b-1) == '1')
        {b-=3; b[0]='\0';}
      else	
	  if (n1+act1_n > 1) strcat(b++,")");
	  dl=leng+strlen(c);
	  strcat(b++,"+");
	 }
      }
    else ster(dl2);
two2:
    act1_n=act_n;
    act1=(SOURCE*) malloc((act1_n+1)*sizeof(SOURCE));
    if (act1_n) actcopy(act1_n,act,act1);
    nui(s1,s2,&act1_n,act1);
    dl1=leng+strlen(c);
    gggp(n,matr,act1_n,act1,pln,st_pln);
    if (dl1 >= leng+strlen(c)) ster(dl);
    free(act1);
    freact(act_n,act);
   }  /* end r-extraction */
  }
  else {frematr(n,matr); freact(act_n,act);}
 } 
 }
 /* Выделение элементов с максимальным участием (общих множителей) */
 void multipl(int *n,PASSIVE *matr,int *act_n,SOURCE *act,int *fl)
{
 int p2,sign,signs,numb,tnum,tn,sens,flm;
 char s1,s2,*t,*m,m0[200];
 unsigned long dl1,dl2;

 newtral(n,matr,act_n,act,fl);
 if (*fl == 2) return;
 node_n(*n,matr,*act_n,act,&p2);
 if (!p2 && *(b-1)=='*') {strcat(b++,"1"); *fl=1; return;}
 flm=0;
    if (!flag_nul && mqtloop(*act_n,act,&s1,&s2,&t,&numb,&sign)
   || pasloop(n,matr,&s1,&s2,&t,&numb,&sign)
   || pashang(n,matr,*act_n,act,&s1,&s2,&t,&numb,&sign)
 || redur(n,matr,*act_n,act,&s1,&s2,&t,&numb,&sign)
 || degener(act_n,act,&s1,&s2,&t,&numb,&sign)
 || nodred(*n,matr,act_n,act,&s1,&s2,&t,&numb,&sign)
  || nodreda(*n,act_n,act,&s1,&s2,&t,&numb,&sign,&flm)
    )
   {
    signs=0;
    tnum=0;
    m0[0]='\0'; m=m0;
rep:
    if (flm==1) {*fl=2; goto ret;}
    if (sign && signs || !sign && !signs) signs=0; else signs=1;
    if (t[0]=='1') goto nopr;
    tn=tnum;
    if (abs(numb) == 1) tnum+=strlen(t)+1;
    else if (abs(numb) > 1) tnum+=strlen(t)+3;
    m0[tn]='\0'; m=m0+tn;

    printi(&m,t,numb);
nopr:
    free(t);
    if (*n || *act_n)
      {
       if (s1 != ' ' && s1 != s2)
	 {
	  if (s1=='0') {s1=s2; s2='0';}
	  unipas(s1,s2,n,matr); uniact(s1,s2,act_n,act);
	 }
       mirror(*act_n,act);
   	   newtral(n,matr,act_n,act,fl);
       if (*fl == 2) return;
       node_n(*n,matr,*act_n,act,&p2);
       if (!p2 || p2==1 && !(*n)  && !(*act_n)) goto cont;
       flm=0;
       if (!flag_nul && mqtloop(*act_n,act,&s1,&s2,&t,&numb,&sign)
   || pasloop(n,matr,&s1,&s2,&t,&numb,&sign)
   || pashang(n,matr,*act_n,act,&s1,&s2,&t,&numb,&sign)
 || redur(n,matr,*act_n,act,&s1,&s2,&t,&numb,&sign)
 || degener(act_n,act,&s1,&s2,&t,&numb,&sign)
 || nodred(*n,matr,act_n,act,&s1,&s2,&t,&numb,&sign)
 || nodreda(*n,act_n,act,&s1,&s2,&t,&numb,&sign,&flm)
	  ) goto rep;
cont:
       dl2=leng+strlen(c);
 if (!strlen(m0)) {strcat(m++,"1"); m[0]='\0';}
       else {m--; m[0]='\0';} 
       printa(m0,1,signs); 
       if (!p2 || p2==1) goto ret;
       if (strlen(m0)) strcat(b++,"*");
       if (   (*(b-3) == '(' || *(b-3) == '+'
	   || *(b-3) == '-' || strlen(c) == 2)
	   && *(b-2) == '1' && *(b-1) == '*')
	 {
	  b-=2; b[0]='\0'; 
	 }
       if (*(b-1)=='-') sens=1; else sens=0;
       if (*n + *act_n > 1 && (m0[0] != '1' || sens)) strcat(b++,"(");
       dl1=leng+strlen(c);
       gggf(*n,matr,*act_n,act);
       if (dl1>=leng+strlen(c)) {ster(dl2); *fl=2;  goto ret;}
       else
	 {
	  if (*n + *act_n > 1 && (m0[0] != '1' || sens)) strcat(b++,")");
	 }
      }
    else
      { 
        if (!strlen(m0)) {strcat(m++,"1"); m[0]='\0';}
        else {m--; m[0]='\0';} 
        printa(m0,1,signs); 
      }
    *fl=1; 
    goto ret;
   }
 *fl=3; 
 return;
ret: ; 
}
/* "ЇаRй?-Ёп ЇаЁ Ї а <<?<м-R┐ бR?¤Ё-?-ЁЁ ??-?а вRаRў Ё ЇаЁ?┐-ЁЄRў "? */
 void redpar(int *act_n,SOURCE *act, int *fl)
{
 int i,j;
 char typi,typj;

rep:
 for (i=0;i<*act_n;i++)
    {
     if (!act[i].kol) continue;
     for (j=0;j<*act_n;j++)
	{
	 if (i == j || !act[j].kol) continue;
	 typi=act[i].reb[0]; typj=act[j].reb[0];
	 if (   act[j].v1==act[i].v1 && act[j].v2==act[i].v2
	     || act[j].v1==act[i].v2 && act[j].v2==act[i].v1)
	   {
	    if (   (   typi == 'K' || typi == 'E' || typi == 'U'
		    || typi == 'H' || act[i].kol < 0)
		&& (   typj == 'K' || typj == 'E' || typj == 'U'
		    || typj == 'H' || act[j].kol < 0))
	      {
	       *fl=2; return;
	      }
	    else
	    if (   (   typi == 'K' || typi == 'E' || typi == 'U'
		    || typi == 'H' || act[i].kol < 0)
		&& ((   typj == 'G' || typj == 'B'
		     || typj == 'J' || typj == 'I') && act[j].kol > 0))
	     {
	      autontr999(j,act_n,act,fl);
	      if (*fl == 2 || *fl == 3) return;
	      else goto rep;
	     }
	   }
	 if (   act[j].v3==act[i].v3 && act[j].v4==act[i].v4
	     || act[j].v3==act[i].v4 && act[j].v4==act[i].v3)
	   {
	    if (   (typi == 'B' || typi == 'H' || act[i].kol < 0)
		&& (typj == 'B' || typj == 'H' || act[j].kol < 0))
	      {
	       *fl=2; return;
	      }
	    else
	      if (   (  (typi == 'G' || typi == 'K') && act[i].kol > 0)
		  && (  typj == 'B' || typj == 'H'
		   || 
act[j].kol < 0 || act[j].kol == 999))
		{
		 autontr999(i,act_n,act,fl);
		 if (*fl == 2 || *fl == 3) return;
		 else goto rep;
		}
	   }
	 if (   act[i].v1==act[j].v3 && act[i].v2==act[j].v4
	     || act[i].v1==act[j].v4 && act[i].v2==act[j].v3)
	   {
	    if (   act[i].kol < 0
		&& (   //typj == 'B'
		 //   || 
        typj == 'H') && act[j].kol > 0
		&& act[j].kol != 999)
	      {
	       act[j].kol = -act[j].kol;
	       goto rep;
	      }
	   }
	 if (   act[i].v3==act[j].v1 && act[i].v4==act[j].v2
	     || act[i].v3==act[j].v2 && act[i].v4==act[j].v1)
	   {
	    if (   act[i].kol < 0
		&& (   typj == 'H' || typj == 'K' 
		    || typj == 'E' || typj == 'U') && act[j].kol > 0)
	      {
	       act[j].kol = -act[j].kol;
	       if (act[j].kol==-999)
		 {
		  ntr999(act_n,act);
		  goto rep;
		 }
	      }
	   }
	}
    }
}
/* Генератор формулы схемного определителя */
 void gggf(int n,PASSIVE *matr,int act_n,SOURCE *act)
{
 int fl,p1,p2,p3,first,fm,smf,n1,act1_n,extract0,q,
     sign,numb,signs;
 char *t,*m,m0[100];
 char type,send1,send2,rec1,rec2,s1,s2,sm1,sm2;
 unsigned long dl,dl1,dl2;
 PASSIVE *a1;
 SOURCE *act1;

 if (*(b-1)=='-')
   {
    b--; b[0]='\0';
    strcat(b++,"+");
   }
 mirror(act_n,act);
 multipl(&n,matr,&act_n,act,&fl);
 if (fl == 1 || fl == 2 || !n && !act_n) return; 
 bisec1(n,matr,act_n,act,&fl);
 if (fl) return;
 node_n(n,matr,act_n,act,&p1);
 if (flag_fsn2 && p1 > flag_fsn2)
   {
    bisec2(n,matr,act_n,act,&fl,range2);
    if (fl) return;
   }
 if (flag_fi3 && p1 > flag_fi3)
   {
    bisec3(n,matr,act_n,act,&fl,range3);
    if (fl) return;
   }
 if (flag_fi4  &&  p1 > flag_fi4)
   {
    bisec4(n,matr,act_n,act,&fl,range4);
    if (fl) return;
   }
 if (flag_fi5 && p1 > flag_fi5)
   {
    bisec5(n,matr,act_n,act,&fl,range5);
    if (fl) return;
   }
 extract0=extract;
beg:
 m0[0]='\0'; m=m0;
 if (fl_mGc && ideal(&act_n,act,&send1,&send2,
	                        &rec1,&rec2,&t,&numb,&sign))
   {
    signs=0;
rephdi:
    if (sign && signs || !sign && !signs) signs=0; else signs=1;
    if (t[0]=='1') goto nopri;
    printi(&m,t,numb);
nopri:
    free(t);
    if (act_n)
      {
       tighten(send1,send2,rec1,rec2,&act_n,act);
       simply(&n,matr,&act_n,act,&fl);
       if (fl == 2) return;
	   //simplym(&act_n,act,&fl);
       //if (fl) return;
       node(act_n,act,&p1);
       if (ideal(&act_n,act,&send1,&send2,
		   &rec1,&rec2,&t,&numb,&sign)) goto rephdi;
     //  m--; m[0]='\0';
       dl2=leng+strlen(c);
       if (!strlen(m0)) strcat(m++,"1");
       else printa(m0,1,signs);
    //   if (strlen(m0)) strcat(b++,"*");
 /*      if (   (*(b-3) == '(' || strlen(c) == 2)
	   && *(b-2) == '1' && *(b-1) == '*')
	 {
	  b-=2; b[0]='\0';
	 }*/
	p3=0;
	if (*(b-1) != '(' && p1>2)
	  {
	   strcat(b++,"(");
	   p3=1;
	  }
       dl1=leng+strlen(c);
       gggf(n,matr,act_n,act);
       if (dl1>=leng+strlen(c)) ster(dl2);
       else
       if (p3) strcat(b++,")");
      }
    else
      {
       m--; m[0]='\0';
       if (!strlen(m0)) strcat(m0,"1");
       else printa(m0,1,signs);
      }
    return;
   }
 else  
 if (!extract0 && yesg(n,matr))
   {
    hangpas(n,matr,act_n,act,&first);
    if (first==-1)
      {if (!flag_opt) first=0;  
       else choiceg(n,matr,act_n,act,&first);}
    dl2=leng+strlen(c);
    printa(matr[first].reb,matr[first].kol,0);
    dl=leng+strlen(c);
    strcat(b++,"*");
    s1=matr[first].v1; s2=matr[first].v2;
    frmatr(first,&n,matr);
    n1=n; act1_n=act_n;
    if (n1)
      {
       a1=(PASSIVE *) malloc(n1*sizeof(PASSIVE));
       pascopy(n1,matr,a1);
      }
    act1=(SOURCE *) malloc ((act1_n+1)*sizeof(SOURCE));
    if (act1_n) actcopy(act1_n,act,act1);
    nui(s1,s2,&act1_n,act1);
    if (n1+act1_n > 1) strcat(b++,"(");
    dl1=leng+strlen(c);
    gggf(n1,a1,act1_n,act1);
    if (n1) free(a1); free(act1);
    if (dl1>=leng+strlen(c)) ster(dl2);
    else
      {
       if (n1+act1_n > 1) strcat(b++,")");
       dl=leng+strlen(c);
       strcat(b++,"+");
      }
    node_n(n,matr,act_n,act,&p2);
    dl1=leng+strlen(c);
    if (p2 == p1)
      {
       gggf(n,matr,act_n,act);
       if (dl1 >= leng+strlen(c)) ster(dl);
      }
    else
      {
       ster(dl);
       frematr(n,matr); freact(act_n,act);
      }
    return;
   }
 else
 if (!extract0 && yesr(n,matr))
   {
    choicer(n,matr,act_n,act,&first);
    dl2=leng+strlen(c);
    printa(matr[first].reb,matr[first].kol,0);
    dl=leng+strlen(c);
    strcat(b++,"*");
    s1=matr[first].v1; s2=matr[first].v2;
    frmatr(first,&n,matr);
    n1=n; act1_n=act_n;
       if (n1)
	 {
	  a1=(PASSIVE *) malloc(n1*sizeof(PASSIVE));
	  pascopy(n1,matr,a1);
	 }
       if (act1_n)
	 {
	  act1=(SOURCE *) malloc (act1_n*sizeof(SOURCE));
	  actcopy(act1_n,act,act1);
	 }
    node_n(n1,a1,act1_n,act1,&p2);
    if ((n1 || act1_n) && p2 == p1)
      {

       if (n1+act1_n > 1) strcat(b++,"(");
       dl1=leng+strlen(c);
       gggf(n1,a1,act1_n,act1);
       if (n1) free(a1); if (act1_n) free(act1);
       if (dl1>=leng+strlen(c)) ster(dl2);
       else
	 {
	  if (n1+act1_n > 1) strcat(b++,")");
	  dl=leng+strlen(c);
	  strcat(b++,"+");
	 }
      }
    else ster(dl2);
    act1=(SOURCE*) malloc((act_n+1)*sizeof(SOURCE));
    actcopy(act_n,act,act1);
    nui(s1,s2,&act_n,act1);
    dl1=leng+strlen(c);
    gggf(n,matr,act_n,act1);
    free(act1);
    if (dl1 >= leng+strlen(c)) ster(dl);
    freact(act1_n,act);
    return;
   }
 else
 if (extract0 && noideal(act_n,act,&first))    
   {
    type=act[first].reb[0];
    if (n)
      {
       a1=(PASSIVE *) malloc(n*sizeof(PASSIVE));
       pascopy(n,matr,a1);
      }
    act1_n=act_n;
    act1=(SOURCE *) malloc (act1_n*sizeof(SOURCE));
    actcopy(act1_n,act,act1);
    act[first].kol=-act[first].kol;
    if (act[first].kol==-999) ntr999(&act_n,act);
    dl=leng+strlen(c);
    if (type == 'K' || type == 'B')
      {
       fl=0;
       if (!n && act_n==1) {strcat(b++,"("); fl=1;}
      }
    dl1=leng+strlen(c);
    gggf(n,matr,act_n,act);
    if (dl1 >= leng+strlen(c)) ster(dl);
    dl=leng+strlen(c);
    if (dl1 < leng+strlen(c)) strcat(b++,"+");
    if (type == 'K' || type == 'B')
      {
       if (fl)
	 {
	  strcat(b,"1)"); b+=2;
	  freematr(n,a1); freeact(act1_n,act1);
	  return;
	 }
      }
    autontr999(first,&act1_n,act1,&fl);
    if (fl == 1 || fl == 3)
      {
       dl1=leng+strlen(c);
       gggf(n,a1,act1_n,act1);
       if (n) free(a1); if (act1_n) free(act1);
       if (dl1>=leng+strlen(c)) ster(dl);
      }
    else
      {
       ster(dl);
       freematr(n,a1); freeact(act1_n,act1);
      }
    return;
   }	 
 else
   {
    if (extract0 == 0) extract0=1;
    else if (extract0 == 1) extract0=0;
    if (!flag_nul) mirror(act_n,act);
      multipl(&n,matr,&act_n,act,&fl);
	  if (fl == 1 || fl == 2) return; // det схемы!=0 || !det 
     goto beg;}
}
/*************************************/
 void mirror(int act_n,SOURCE *act)
{
 int i;
 char s1,s2,s3,s4,s0;

 for (i=0;i<act_n;i++)
    {
     if (act[i].reb[0]=='M')
       {
	if (act[i].v3!='0' && act[i].v4!='0')
	  {
	   if (   act[i].v1=='0' && act[i].v2!='0'
	       || act[i].v1!='0' && act[i].v2=='0' )
	     {
	      act[i].reb[0]='Q';
	      if (act[i].v1=='0') {act[i].v1=act[i].v2; act[i].v2='0';}
	     }
	  }
	else
	if (act[i].v1!='0' && act[i].v2!='0')
	  {
	   if (   act[i].v3=='0' && act[i].v4!='0'
	       || act[i].v3!='0' && act[i].v4=='0' )
	     {
	      act[i].reb[0]='T';
	      if (act[i].v3=='0') {act[i].v3=act[i].v4; act[i].v4='0';}
	     }
	  }
	else
	if (   act[i].v1=='0' && act[i].v2!='0'
	    && act[i].v3=='0' && act[i].v4!='0'
	    ||
	       act[i].v1!='0' && act[i].v2=='0'
	    && act[i].v3!='0' && act[i].v4=='0'
	    ||
	       act[i].v1=='0' && act[i].v2!='0'
	    && act[i].v3!='0' && act[i].v4=='0'
	    ||
	       act[i].v1!='0' && act[i].v2=='0'
	    && act[i].v3=='0' && act[i].v4!='0')
	     {
	      act[i].reb[0]='1';
	      if (act[i].v1=='0') {act[i].v1=act[i].v2; act[i].v2='0';}
	      if (act[i].v3=='0') {act[i].v3=act[i].v4; act[i].v4='0';}
	     }
       }
     else
     if (act[i].reb[0]=='Q')
       {
	if (   act[i].v3!='0' && act[i].v4=='0'
	    || act[i].v3=='0' && act[i].v4!='0')
	  {
	   act[i].reb[0]='1';
	   if (act[i].v3=='0') {act[i].v3=act[i].v4; act[i].v4='0';}
	  }
       }
     else
     if (act[i].reb[0]=='T')
       {
	if (   act[i].v1!='0' && act[i].v2=='0'
	    || act[i].v1=='0' && act[i].v2!='0')
	  {
	   act[i].reb[0]='1';
	   if (act[i].v1=='0') {act[i].v1=act[i].v2; act[i].v2='0';}
	  }
       }
    }
}
/****************************************************/
void redum(int act_n,SOURCE *act)
{
 int i,j;
 char s1,s2,s3,s4,s0;

 for (i=0;i<act_n;i++)
    {
     if (   act[i].reb[0]=='M' || act[i].reb[0]=='Q')
       {
	for (j=0;j<act_n;j++)
	   {
	    if (i==j) continue;
	    if (   act[j].reb[0] == '1' || act[j].reb[0] == 'T')
	      {
	       if (   act[i].v3==act[j].v3 && act[j].v4=='0'
		   || act[i].v3==act[j].v4 && act[j].v3=='0')
		 {
		  if (act[i].reb[0]=='M') act[i].reb[0]='T';
		  else act[i].reb[0]='1';
		  act[i].v3=act[i].v4;
		  act[i].v4='0';
		  break;
		 }
	       else
		 if (   act[i].v4==act[j].v4 && act[j].v3=='0'
		     || act[i].v4==act[j].v3 && act[j].v4=='0')
		   {
		    if (act[i].reb[0]=='M') act[i].reb[0]='T';
		    else act[i].reb[0]='1';
		    act[i].v4='0';
		    break;
		   }
	      }
	   }
       }
    }
 for (i=0;i<act_n;i++)
    {
     if (   act[i].reb[0]=='M' || act[i].reb[0]=='T')
       {
	for (j=0;j<act_n;j++)
	   {
	    if (i==j) continue;
	    if (   act[j].reb[0] == '1' || act[j].reb[0] == 'Q')
	      {
	       if (   act[i].v1==act[j].v1 && act[j].v2=='0'
		   || act[i].v1==act[j].v2 && act[j].v1=='0')
		 {
		  if (act[i].reb[0]=='M') act[i].reb[0]='Q';
		  else act[i].reb[0]='1';
		  act[i].v1=act[i].v2;
		  act[i].v2='0';
		  break;
		 }
	       else
		 if (   act[i].v2==act[j].v2 && act[j].v1=='0'
		     || act[i].v2==act[j].v1 && act[j].v2=='0')
		   {
		    if (act[i].reb[0]=='M') act[i].reb[0]='Q';
		    else act[i].reb[0]='1';
		    act[i].v2='0';
		    break;
		   }
	      }
	   }
       }
    }
}
/******************** Выделение z-петли ********************/
 int mqtloop(int act_n,SOURCE *act,
	     char *s1,char *s2,char **t,int *numb,int *sign)
{
 int i,ind=0,ipr=0;
 char *w,w0[20];

 *sign=0;
 w0[0]='\0'; w=w0;
rep:
 ipr++;
 for (i=0;i<act_n;i++)
    {
     if (   (act[i].reb[0]=='M' && act[i].kol < 0
	 || act[i].reb[0]=='T' && act[i].kol < 0)
	 && act[i].v1 == act[i].v2)
       {
	printf(" r_l2");
	strcat(w++,"2");
	*t=strdup(w0);
	*numb=1;
	*s1=act[i].v1; *s2=act[i].v2;
	*sign=0;
	if (act[i].reb[0]=='M') act[i].reb[0]='Q';
	else act[i].reb[0]='N';
	act[i].v2='0';
	ind=1;
	goto con;
       }
    }
con:
 for (i=0;i<act_n;i++)
    {
     if (  (act[i].reb[0]=='M' && act[i].kol < 0
	 || act[i].reb[0]=='Q' && act[i].kol < 0)
	 && act[i].v3 == act[i].v4)
       {
	if (ind)
	  {
	   printf(" r_l4");
	   w0[0]='4';
	  }
	else
	  {
	   printf(" r_l2"); strcat(w++,"2");
	  }
	*t=strdup(w0);
	*numb=1;
	*s1=act[i].v3; *s2=act[i].v4;
	*sign=0;
	if (act[i].reb[0]=='M') act[i].reb[0]='T';
	else act[i].reb[0]='N';
	act[i].v4='0';
	ind=1;
       }
    }
 if (ipr<5) goto rep;
 if (ind) return(1); else return(0);
}
/*************************************/
int hangedm(int act_n,SOURCE *act)
{
 int j,i,p,send,rec;
 char *str;
 GRAPH *a3;

 if (!act_n) return(1);
 a3=(GRAPH *) malloc ((2*act_n+1)*sizeof(GRAPH));
 copyactm(act_n,act,a3);
 str=(char*) malloc(4*act_n+1);
 verstr(2*act_n,a3,&p,str); free(a3);
 for (i=0;i<p;i++)
    {
     send=0; rec=0;
     for (j=0;j<act_n;j++)
	{
	 if (str[i]==act[j].v1 || str[i]==act[j].v2) send=1;
	 if (str[i]==act[j].v3 || str[i]==act[j].v4) rec=1;
        }
     if (send && !rec || !send && rec)
       {
	free(str); return(0);
       }
    }
 free(str);
 return(1);
}
/***************************************************************/
int connecm(int act_n,SOURCE *act)
{
 int j,na,pr;
 char kr;
 GRAPH *a2;

 if (!act_n) return(1);
 a2=( GRAPH *) malloc ((2*act_n+1)*sizeof(GRAPH));
 copyactm(act_n,act,a2);
 na = 2*act_n;
 pr=0; kr=a2[0].v1; bond1f(na,a2,kr,kr,&pr); free(a2);
 if (na>pr) return(0); else return(1);
}
/************************************************************/
void trianglm(int *act_n,SOURCE *act,int *flag)
{
 int i,j,k,l,m,pr=0,p1,p2;
 char s1,s2,t1,t2,x1,x2;

 *flag=0;
 for (i=0;i<*act_n;i++) if (act[i].kol < 0)
    {
     node(*act_n,act,&p1);
     goto beg;
    }
 return;
beg:
 for (i=0;i<*act_n;i++)
    {
     if (act[i].kol > 0) continue;
     if (!pr) {s1=act[i].v1; s2=act[i].v2;}
     else {s1=act[i].v3; s2=act[i].v4;}
     for (j=0;j<*act_n;j++)
	{
	 if (act[j].kol > 0) continue;
	 if (j==i) continue;
	 if (!pr) {t1=act[j].v1; t2=act[j].v2;}
	 else {t1=act[j].v3; t2=act[j].v4;}
	 if (t1==s1) {x1=t2; x2=s2; goto sear;}
	 if (t1==s2) {x1=t2; x2=s1; goto sear;}
	 if (t2==s1) {x1=t1; x2=s2; goto sear;}
	 if (t2==s2) {x1=t1; x2=s1; goto sear;}
	 continue;
sear:
	 for (m=0;m<*act_n;m++)
	    {
	     if (m==i || m==j) continue;
	     if (!pr) {t1=act[m].v1; t2=act[m].v2;}
	       else {t1=act[m].v3; t2=act[m].v4;}
	     if (t1==x1 && t2==x2 || t1==x2 && t2==x1)
	       {
		if (act[m].kol < 0) {pr=3; goto ret;}
		free(act[m].reb);
		(*act_n)--;
		for (k=m;k<*act_n;k++) act[k]=act[k+1];
		if (i >= m) i--; if (j >= m) j--; m--;
	       }
	    }
	}
    }
 pr++;
 if (pr==1) goto beg;
ret:
 node(*act_n,act,&p2);
 if (pr==3 || p1!=p2) *flag=1;
}
/************************************************************/
void redumm(int *act_n,SOURCE *act,int *flag)
{
 int i,k,pr=0,m;
 char s1,s2,t1,t2;
 int p1,p2,j;

 *flag=0;
 if (!*act_n) return;
 for (i=0;i<*act_n;i++) if (act[i].kol < 0)
    {
     node(*act_n,act,&p1);
     goto beg;
    }
 return;
beg:
 for (i=0;i<*act_n;i++)
    {
     if (act[i].kol > 0) continue;
     if (!pr) {s1=act[i].v1; s2=act[i].v2;}
     else {s1=act[i].v3; s2=act[i].v4;}
     for (m=0;m<*act_n;m++)
	{
	 if (m==i) continue;
	 if (!pr) {t1=act[m].v1; t2=act[m].v2;}
	 else {t1=act[m].v3; t2=act[m].v4;}
	 if (t1==s1 && t2==s2 || t1==s2 && t2==s1)
	   {
	    if (act[m].kol < 0) {pr=3; goto ret;}
	    free(act[m].reb);
	    (*act_n)--; for (k=m;k<*act_n;k++) act[k]=act[k+1];
	    if (i >= m) i--; m--;
	   }
	}
    }
 pr++;
 if (pr==1) goto beg;
ret:
 node(*act_n,act,&p2);
 if (pr==3 || p1 != p2) *flag=1;
}
/**************************************************************/
int degenerm(int *act_n,SOURCE *act,
	  char *s1,char *s2,char **t,int *numb,int *sign)
{
 int i,j;

 for (i=0;i<*act_n;i++)
    {
     if (act[i].kol > 0) continue;
     for (j=0;j<*act_n;j++)
	{
	 if (act[j].kol > 0) continue;
	 *sign=-1;
	 if (act[i].v1 == act[j].v3 && act[i].v2 == act[j].v4) *sign=0;
	 if (act[i].v1 == act[j].v4 && act[i].v2 == act[j].v3) *sign=1;
	 if (*sign != -1)
	   {
	    if (i != j)
	      {
	       if (!*sign)
		 {
		  act[j].v3=act[i].v4;
		  act[j].v4=act[i].v3;
		 }
	       else
		 {
		  act[j].v3=act[i].v3;
		  act[j].v4=act[i].v4;
		  *sign=0;
		 }
	      }
	    *t=strdup(act[i].reb);
	    *numb=act[i].kol;
	    *s1=act[i].v1; *s2=act[i].v2;
	    free(act[i].reb);
	    (*act_n)--; for (j=i;j<*act_n;j++) act[j]=act[j+1];
	    return(1);
	   }
	}
    }
 return(0);
}
/*********************************************************/
 int smplhng(int act_n,SOURCE *act,int kand)
{
 int j,i,p,send,rec;
 char *str;
 GRAPH *a3;

 if (!act_n) return(1);
 a3=(GRAPH *) malloc ((2*act_n+1)*sizeof(GRAPH));
 copyactm(act_n,act,a3);
 str=(char*) malloc(4*act_n+1);
 verstr(2*act_n,a3,&p,str); free(a3);
 for (i=0;i<p;i++)
    {
     send=0; rec=0;
     for (j=0;j<act_n;j++)
	{
	 if (j == kand) continue;
	 if (str[i]==act[j].v1 || str[i]==act[j].v2) send++;
	 if (str[i]==act[j].v3 || str[i]==act[j].v4) rec++;
	}
     if (send && !rec || !send && rec)
       {
	free(str); return(0);
       }
    }
 free(str);
 return(1);
}
/* Упрощения при параллельном соединении генераторов и приемников УИ */
 void redparm(int *act_n,SOURCE *act, int *fl)
{
 int i,j;

 *fl=0;
rep:
 for (i=0;i<*act_n;i++)
    {
     if (!act[i].kol) continue;
     for (j=0;j<*act_n;j++)
	{
	 if (i == j || !act[j].kol) continue;
	 if (   act[j].v1==act[i].v1 && act[j].v2==act[i].v2
	     || act[j].v1==act[i].v2 && act[j].v2==act[i].v1)
	   {
	    if (act[i].kol < 0 && act[j].kol < 0)
	      {
	       *fl=1; return;
	      }
		else
		if (   act[i].kol > 0 && act[i].kol != 888
	        && act[j].kol == 888)
	      fract(i,act_n,act);
        else
		if (   act[i].kol == 888
	        && act[j].kol > 0 && act[j].kol != 888)
	      fract(j,act_n,act);
	   }
	 if (   act[j].v3==act[i].v3 && act[j].v4==act[i].v4
	     || act[j].v3==act[i].v4 && act[j].v4==act[i].v3)
	   {
	    if (act[i].kol < 0	&& act[j].kol < 0)
	      {
	       *fl=1; return;
	      }
		else
		if (   act[i].kol > 0 && act[i].kol != 888
	        && act[j].kol == 888)
	      fract(i,act_n,act);
        else
		if (   act[i].kol == 888
	        && act[j].kol > 0 && act[j].kol != 888)
	      fract(j,act_n,act);
	   }
	}
    }
}
/**********************************************************/
void simplym(int *act_n,SOURCE *act,int *flag)
{
 int fl,i;

 *flag=0;
 if (!*act_n) return;
 else
   {
    redparm(act_n,act,&fl);
    if (fl) goto ret;
    parall(act_n,act);
	for (i=0;i<*act_n;i++)
       {
	if (act[i].kol < 0) continue;
	if (!smplhng(*act_n,act,i))
	  {
	   act[i].kol=-act[i].kol;
	  }
       }
    trianglm(act_n,act,&fl);
    if (fl) goto ret;
    redumm(act_n,act,&fl);
    if (fl) goto ret;
   }
 if (!hangedm(*act_n,act) || !connecm(*act_n,act)) goto ret;
 return;
ret:
 *flag=1;
 freact(*act_n,act);
}
/**************************************************************/
int ideal(int *act_n,SOURCE *act,
	  char *send1,char *send2,char *rec1,char *rec2,
	  char **t,int *numb,int *sign)
{
 int i,first;
 char s1,s2,s3,s4;

 first=-1;
 for (i=0;i<*act_n;i++)
    {
     if (act[i].kol > 0) continue;
     first=i;
    }
 if (first== -1) return(0);
beg:
 s1=act[first].v1; s2=act[first].v2;
 s3=act[first].v3; s4=act[first].v4;
 //if (s1=='0') {s1=s2; s2='0';}
 //if (s3=='0') {s3=s4; s4='0';}
 if (s1==s3)
   {*send1=s1; *rec1=s3; *send2=s2; *rec2=s4; *sign=0;}
 else
 if (s2==s4)
   {*send1=s2; *rec1=s4; *send2=s1; *rec2=s3; *sign=0;}
 else
 if (s1==s4)
   {*send1=s1; *rec1=s4; *send2=s2; *rec2=s3; *sign=1;}
 else
 if (s2==s3)
   {*send1=s2; *rec1=s3; *send2=s1; *rec2=s4; *sign=1;}
 else
   {
      *send1=s1; *rec1=s4; *send2=s2; *rec2=s3; *sign=0;
/*      *send1=s2; *rec1=s3; *send2=s1; *rec2=s4; *sign=0;
      *send1=s1; *rec1=s3; *send2=s2; *rec2=s4; *sign=1;
      *send1=s2; *rec1=s4; *send2=s1; *rec2=s3; *sign=1; */
   }
 // *c1=' '; *c2='i';
 *t=strdup(act[first].reb);
 *numb=act[first].kol;
 free(act[first].reb);
 (*act_n)--; for (i=first;i<*act_n;i++) act[i]=act[i+1];
 return(1);
}
/**********************************************************/
 int maxsoum(int *act_n,SOURCE *act,
	    char *c1,char *c2,char **t,int *numb,int *sign)
{
 int na,j,i,k,p,send,rec,sou1,sou2,in,ou;
 char *str;
 GRAPH *a3;

 if (!act_n) return(0);
 a3=(GRAPH *) malloc ((2*(*act_n)+1)*sizeof(GRAPH));
 copyactm(*act_n,act,a3);
 na = 2*(*act_n);
 str=(char*) malloc(4*(*act_n)+1);
 verstr(na,a3,&p,str);
 free(a3);
 for (i=0;i<p;i++)
    {
     send=0; rec=0;
     for (j=0;j<*act_n;j++)
	{
	 if (str[i]==act[j].v1) {sou1=j; in=0; send++;}
	 if (str[i]==act[j].v2) {sou1=j; in=1; send++;}
	 if (str[i]==act[j].v3) {sou2=j; ou=0; rec++;}
	 if (str[i]==act[j].v4) {sou2=j; ou=1; rec++;}
	}
     if (send==1 && rec==1 && sou1!=sou2) goto beg;
       else continue;
beg:
     *c1=' '; *c2='m';
     if (in==0 && ou==0 || in==1 && ou==1) *sign=1;
     if (in==0 && ou==1 || in==1 && ou==0) *sign=0;
     if (act[sou1].kol>0) act[sou1].kol=-act[sou1].kol;
     if (act[sou2].kol>0) act[sou2].kol=-act[sou2].kol;
     *t=strdup(act[sou1].reb);
     *numb=act[sou1].kol;
     act[sou2].v3=act[sou1].v3;
     act[sou2].v4=act[sou1].v4;
     free(act[sou1].reb);
     (*act_n)--; for (i=sou1;i<*act_n;i++) act[i]=act[i+1];
     free(str); return(1);
    }
 free(str); return(0);
}
/*******************************************/
 void quartet(int act_n,SOURCE *act,int *fl)
{ 
 int j,p;
 unsigned long dl,dl1;
 SOURCE *act1;
 char *w,w0[2];

 *fl=0;
 w0[0]='\0'; w=w0;
 strcat(w++,"1");
 for (j=0;j<4;j++) 
 if (act_n==1 || act[j].kol != 888 || act[j].kol == 555) return;
 act1=(SOURCE *) malloc (4*sizeof(SOURCE));
 for (j=0;j<4;j++)
    {
     act1[j].v1=act[j].v1; act1[j].v2=act[j].v2;
     act1[j].v3=act[j].v3; act1[j].v4=act[j].v4;
     act1[j].kol=act[j].kol;
     act1[j].reb=strdup(act[j].reb);
	}
 for (j=0;j<4;j++)
    {
     act1[j].v2='0';
     act1[j].v4='0';
	}
 act1[0].v1='1';
 act1[0].v3='1';
 act1[1].v1='1';
 act1[1].v3='2';
 act1[2].v1='2';
 act1[2].v3='1';
 act1[3].v1='2';
 act1[3].v3='2';
 for (j=0;j<4;j++) {act1[j].kol = 555; act[4+j].kol = 555;}
 act[0].kol=-1; //act[0].kol;
 act[3].kol=-1;  //act[3].kol;
 free(act[0].reb); free(act[3].reb);
 act[0].reb=strdup(w0); act[3].reb=strdup(w0);
 strcat(b++,"(");
 dl1=leng+strlen(c);
 gggfm(4,act1);
 free(act1);
// dl=leng+strlen(c);
 if (dl1 < leng+strlen(c))
 {strcat(b++,")");
 strcat(b++,"*");
 }
 node(act_n,act,&p);
 dl1=leng+strlen(c);
 if (p)
   {
	strcat(b++,"(");
	gggfm(act_n,act);
    if (dl1<leng+strlen(c)) strcat(b++,")");
    else ster(dl1);
   }	
 *fl=1;  
} 
/*******************************************/
 int noideal(int act_n,SOURCE *act,int *first)
{
 int i,j,k,p,pow0[200],pow1[200],posi[200],posi_n=0;
 char *str,min,minj,mi0,mi0j,mi1,mi1j;

 *first = -1;
 for (i=0;i<act_n;i++)
    { 
     if (act[i].kol < 0) continue; 
	 posi[posi_n++]=i;
	}
 if (!posi_n) return(0);
 str=(char*) malloc(4*act_n+1);
 nodesm(act_n,act,&p,str);
 for (j=0;j<p;j++)
    {
     pow0[j]=0; pow1[j]=0;		
     for (k=0;k<posi_n;k++)
        {
         i=posi[k];
         if (   act[i].v1 == str[j] 
		     || act[i].v2 == str[j]) pow0[j]++;
		 if (   act[i].v3 == str[j] 
		     || act[i].v4 == str[j]) pow1[j]++;
        }
    }
 mi0=100; mi1=100; min=100;
 for (j=0;j<p;j++)
    {	 
     if (pow0[j] && pow0[j] < mi0) {mi0=pow0[j]; mi0j=j;}
     if (pow1[j] && pow1[j] < mi1) {mi1=pow1[j]; mi1j=j;}
     if (pow0[j] && pow1[j] && pow0[j]+pow1[j] < min) 
       {min=pow0[j]+pow1[j]; minj=j;} 
    }
 for (k=0;k<posi_n;k++)
    {
     i=posi[k];	
     if (   (act[i].v1 == str[mi0j] || act[i].v2 == str[mi0j]) 
         && (act[i].v3 == str[mi1j] || act[i].v4 == str[mi1j])) 
	   {*first=i; goto ret;}
    } 
 if (pow0[mi0j] < pow1[mi1j])
 for (k=0;k<posi_n;k++)
    {
     i=posi[k]; 
     if (   act[i].v1 == str[mi0j] || act[i].v2 == str[mi0j]) 
	   {*first=i; goto ret;}
    }	
 else
 for (k=0;k<posi_n;k++)
    {
     i=posi[k]; 
     if (act[i].v3 == str[mi1j] || act[i].v4 == str[mi1j]) 
	   {*first=i; goto ret;}
    }
 for (k=0;k<posi_n;k++)
    {
     i=posi[k]; 
     if (   (act[i].v1 == str[minj] || act[i].v2 == str[minj]) 
         && (act[i].v3 == str[minj] || act[i].v4 == str[minj])) 
	   {*first=i; goto ret;}
    }   
 if (pow0[minj] < pow1[minj])
 for (k=0;k<posi_n;k++)
    {
     i=posi[k]; 
     if (   act[i].v1 == str[minj] || act[i].v2 == str[minj]) 
	   {*first=i; break;}
    }	
 else
 for (k=0;k<posi_n;k++)
    {
     i=posi[k]; 
     if (act[i].v3 == str[minj] || act[i].v4 == str[minj]) 
	   {*first=i; break;}
    }
ret:	
 free(str);
 if (*first==-1) return(0); else return(1);
}
/**************************************************************/
 void noid_opt(int act_n,SOURCE *act,int *first)
{
 int i,j,m,p,p0,pow[200];
 char *str,*str0,min;

 *first = 0;
 str0=(char*) malloc(4*act_n+1);
 str=(char*) malloc(4*act_n+1);
 nodesm(act_n,act,&p,str0);
 p0=p;
 for (m=0;m<p;m++)
 {
 min=100; 
 for (i=0;i<p;i++)
    if (str0[i]<min) min=str0[i]; 
 str[m]=min;
 for (i=0;i<p0;i++)
    if (str0[i]==min) 
	 {p0--; for (j=i;j<p0;j++) str0[j]=str0[j+1]; break;}
 }
 for (j=0;j<p;j++) pow[j]=0;
 for (j=0;j<p;j++)
    {
     for (i=0;i<act_n;i++)
        {
         if (act[i].kol < 0) continue;
         if (act[i].v1 == str[j]) pow[j]++;
        }
    }
 free(str);
 for (m=0;m<p;m++) if (pow[m]) break;
 if (pow[m] > 2) *first=2;
 else if (pow[m] == 2)
 {
  for (i=m+1;i<p;i++) if (pow[i]) goto br;
  return;
br:
 if (pow[i] > 2) *first=4;
 else
 if (pow[i]==2)
   {
    for (j=m+2;j<p;j++)
	   if (pow[j] && act_n > 4 && act[4].kol > 0)
	     {*first=4; break;}
    act[0].kol=888; act[1].kol=888;
    act[2].kol=888; act[3].kol=888;
   }
  }
}
/**************************************************************/
 int noideam(int act_n,SOURCE *act,int *first)
{
 int i,j,jmin,p,pmin,pow[100];
 char *str;

 str=(char*) malloc(4*act_n);
 nodes(act_n,act,&p,str);
 for (j=0;j<p;j++) pow[j]=0;
 *first=-1;
 for (j=0;j<p;j++)
    {
//     if (str[j] == '0') continue;
     for (i=0;i<act_n;i++)
        {
         if (act[i].kol < 0) continue;
         if (act[i].v1 == str[j] || act[i].v2 == str[j])
           pow[j]++;
        }
    }
 for (j=0;j<p;j++)
    {
     if (str[j] == '0') continue;
     if (pow[j] < pmin) {pmin=pow[j]; jmin=j; break;}
    }
 for (i=0;i<act_n;i++)
    {
     if (act[i].v1 == str[jmin] || act[i].v2 == str[jmin])
       {*first=i; break;}
    }
 free(str);
 if (*first==-1) return(0); else return(1);
 }
/********************************************************************/
void gggfm(int act_n,SOURCE *act)
{
 int j,fl,flag=0,p1,p2,p3,first,act1_n,sign,numb,signs;
 char send1,send2,rec1,rec2,s1,s2,*t,*m,m0[200];
 unsigned long dl,dl1,dl2;
 SOURCE *act1;

 simplym(&act_n,act,&flag);
 if (flag) return;
 node(act_n,act,&p1);
 m0[0]='\0'; m=m0;
 if (   maxsoum(&act_n,act,&s1,&s2,&t,&numb,&sign) ||
        degenerm(&act_n,act,&s1,&s2,&t,&numb,&sign)
     || nodredmm(&act_n,act,&s1,&s2,&t,&numb,&sign)
    )
   {
    signs=0;

rephd:
    if (sign && signs || !sign && !signs) signs=0; else signs=1;
    if (t[0]=='1') goto nopr;
    printi(&m,t,numb);
nopr:
    free(t);
    if (act_n)
      {
       if (s1 != ' ') uniact(s1,s2,&act_n,act);
       simplym(&act_n,act,&flag);
       if (flag) return;
       if (s1 != ' ')
	 {
	  node(act_n,act,&p2);
	  if (p2 && p2 != p1-1)
             {for (j=0;j<act_n;j++) free(act[j].reb); return;}
	 }
       node(act_n,act,&p1);
       if (   maxsoum(&act_n,act,&s1,&s2,&t,&numb,&sign) ||
	          degenerm(&act_n,act,&s1,&s2,&t,&numb,&sign)
	       || nodredmm(&act_n,act,&s1,&s2,&t,&numb,&sign)
	      )
	 goto rephd;
       m--; m[0]='\0';
       dl2=leng+strlen(c);
       if (!strlen(m0)) strcat(m0,"1");
       printa(m0,1,signs);
       if (strlen(m0)) strcat(b++,"*");
       if (   (*(b-3) == '(' || strlen(c) == 2)
	   && *(b-2) == '1' && *(b-1) == '*')
	 {
	  b-=2; b[0]='\0';
	 }
	  if (   *(b-4) == '(' && *(b-3) == '-'
       	 && *(b-2) == '1' && *(b-1) == '*')
	 {
	  b-=2; b[0]='\0';
	 }
    p3=0;
	if (*(b-1) != '(' && p1>2)
	  {
	   strcat(b++,"(");
	   p3=1;
	  }
       dl1=leng+strlen(c);
       gggfm(act_n,act);
       if (dl1>=leng+strlen(c)) ster(dl2);
       else
       if (p3) strcat(b++,")");
      }
    else
      {
       m--; m[0]='\0';
       if (!strlen(m0)) strcat(m0,"1");
       printa(m0,1,signs);
      }
    return;
   }
 node(act_n,act,&p1);   
 if (p1 > 4 && flag_opt == 1) {
	 quartet(act_n,act,&flag);
     if (flag) return;} 
 else 
 if (flag_opt == -1)
 {	
  if (flag_fsn1 && p1 > flag_fsn1)
   {	 
    bisec1m(act_n,act,&flag);
    if (flag) return;
   }	
 if (flag_fsn2 && p1 > flag_fsn2)
   {
    bisec2m(act_n,act,&flag,range2);
    if (flag) return;
   }
 if (flag_fi3 && p1 > flag_fi3)
   {
    bisec3m(act_n,act,&flag,range3);
    if (flag) return;
   }
 if (flag_fi4  &&  p1 > flag_fi4)
   {
    bisec4m(act_n,act,&flag,range4);
    if (flag) return;
   }
 if (flag_fi5 && p1 > flag_fi5)
   {
    bisec5m(act_n,act,&flag,range5);
    if (flag) return;
   }
 }  
 m0[0]='\0'; m=m0;
 if (ideal(&act_n,act,&send1,&send2,
	      &rec1,&rec2,&t,&numb,&sign))
   {
    signs=0;

rephdi:
    if (sign && signs || !sign && !signs) signs=0; else signs=1;
    if (t[0]=='1') goto nopri;
    printi(&m,t,numb);
nopri:
    free(t);
    if (act_n)
      {
       tighten(send1,send2,rec1,rec2,&act_n,act);
       simplym(&act_n,act,&flag);
       if (flag) return;
       node(act_n,act,&p1);
       if (ideal(&act_n,act,&send1,&send2,
		    &rec1,&rec2,&t,&numb,&sign)) goto rephdi;
       dl2=leng+strlen(c);
       if (!strlen(m0)) strcat(m0,"1");
       else printa(m0,1,signs);
     /*  if (strlen(m0)) strcat(b++,"*");
       if (   (*(b-3) == '(' || strlen(c) == 2)
	   && *(b-2) == '1' && *(b-1) == '*')
	 {
	  b-=2; b[0]='\0';
	 }*/
	p3=0;
	if (*(b-1) != '(' && p1>2)
	  {
	   strcat(b++,"(");
	   p3=1;
	  }
       dl1=leng+strlen(c);
       gggfm(act_n,act);
       if (dl1>=leng+strlen(c)) ster(dl2);
       else
       if (p3) strcat(b++,")");
      }
    else
      {
       m--; m[0]='\0'; 
       if (!strlen(m0)) strcat(m0,"1");
       else printa(m0,1,signs);
      }
    return;
   }
 if (   flag_opt == -1 && noideal(act_n,act,&first)
	 || flag_opt != -1 && noideam(act_n,act,&first))
   { /* begin noideal */
    SOURCE *act1;

	if (flag_opt == 1 && p1 > 4) noid_opt(act_n,act,&first);
	act1_n=act_n-1;
    if (act1_n)
      {
       act1=(SOURCE *) malloc (act1_n*sizeof(SOURCE));
       delet(act1_n,act,first,act1);
      }
    act[first].kol=-act[first].kol;
    dl1=leng+strlen(c);
    gggfm(act_n,act);
    dl=leng+strlen(c);
    if (dl1 < leng+strlen(c)) strcat(b++,"+");
    node(act1_n,act1,&p2);
    dl1=leng+strlen(c);
    if (p1 == p2)
      {
       gggfm(act1_n,act1);
       if (act1_n) free(act1);
      }
    else
      {
       if (act1_n)
	 {
	  for (j=0;j<act1_n;j++) free(act1[j].reb); free(act1);
	 }
      }
    if (dl1>=leng+strlen(c)) ster(dl);
    return;
   } /* end noideal */
 else //freact(act_n,act);
   for (j=0;j<act_n;j++) free(act[j].reb);
}
/**************************************************************/
int nodredmm(int *act_n,SOURCE *act,
	   char  *s1,char *s2, char **t,int *numb,int *sign)
{
 char *str,kr,akr;
 GRAPH *a3;
 int i,j,k,na,p,sn1,sn2,ret=0,sou1,sou2;

 for (i=0;i<*act_n;i++) if (act[i].kol < 0) goto beg;
 return(0);
beg:
 a3=(GRAPH *) malloc ((2*(*act_n)+1)*sizeof(GRAPH));
 copyactm(*act_n,act,a3);
 na = 2*(*act_n);
 str=(char*) malloc(4*(*act_n)+1);
 verstr(na,a3,&p,str);
 free(a3);
 for (i=0;i<p;i++)
    {
     kr=str[i];
     sou1 = -1; sou2 = -1;
     for (j=0;j<*act_n;j++)
	{
	 if (act[j].kol > 0) continue;
	 if (kr==act[j].v1) {sou1=j; sn1=1;}
	 if (kr==act[j].v2) {sou1=j; sn1=0;}
	 if (kr==act[j].v3) {sou2=j; sn2=1;}
	 if (kr==act[j].v4) {sou2=j; sn2=0;}
	}
     if (sou1 == -1 || sou2 == -1) continue;
     if (sn1==sn2) *sign=0; else *sign=1;
     if (sou1 != sou2)
       {
	if (*sign) *sign=0; else *sign=1;
       }
     ret=1;
     break;
    }
 free(str);
 if (!ret) return(0);
 if (act[sou1].v1 == kr) akr=act[sou1].v2; else akr=act[sou1].v1;
 for (k=0;k<*act_n;k++)
    {
     if (k==sou1) continue;
     if (act[k].v1==kr) act[k].v1=akr;
     if (act[k].v2==kr) act[k].v2=akr;
    }
 if (act[sou2].v3 == kr) akr=act[sou2].v4; else akr=act[sou2].v3;
 for (k=0;k<*act_n;k++)
    {
     if (k==sou2) continue;
     if (act[k].v3==kr) act[k].v3=akr;
     if (act[k].v4==kr) act[k].v4=akr;
    }
 if (sou1!=sou2)
   {
    if (!*sign)
      {
       act[sou2].v3=act[sou1].v3;
       act[sou2].v4=act[sou1].v4;
      }
    else
      {
       act[sou2].v3=act[sou1].v4;
       act[sou2].v4=act[sou1].v3;
       *sign=0;
      }
   }
 *t=strdup(act[sou1].reb);
 *numb=act[sou1].kol;
 free(act[sou1].reb);
 (*act_n)--; for (k=sou1;k<*act_n;k++) act[k]=act[k+1];
 *s1=' '; *s2=' ';
 return(1);
}
/****************************************************************/
void tighten(char send1,char send2,char rec1,char rec2,
	     int *act_n,SOURCE *act)
{
 int i,j,k;

 for (i=0;i<*act_n;i++)
    {
     if (act[i].v1==send1) act[i].v1=send2;
     if (act[i].v2==send1) act[i].v2=send2;
     if (act[i].v3==rec1) act[i].v3=rec2;
     if (act[i].v4==rec1) act[i].v4=rec2;
    }
 if (send1==rec1) return;
 for (i=0;i<*act_n;i++)
    {
     if (act[i].v1==send1) act[i].v1=rec1;
     if (act[i].v2==send1) act[i].v2=rec1;
     if (act[i].v3==send1) act[i].v3=rec1;
     if (act[i].v4==send1) act[i].v4=rec1;
    }
}
/***************************************************************/
void delet(int act1_n,SOURCE *act,int first,SOURCE *act1)
{
 int i,i1;

 actcopy(first,act,act1);
 for (i=first;i<act1_n;i++)
    {
     i1=i+1;
     act1[i].v1=act[i1].v1; act1[i].v2=act[i1].v2;
     act1[i].v3=act[i1].v3; act1[i].v4=act[i1].v4;
     act1[i].kol=act[i1].kol;
     act1[i].reb=strdup(act[i1].reb);
    }
}
