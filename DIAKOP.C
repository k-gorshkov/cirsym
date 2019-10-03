/* Разложение схемного определителя по частям - бисекция схемы */

#include "cirsym.h"

extern FILE *out,*inpa,*set;

extern char *b,*c,sp;
extern unsigned long leng,lengl,lepr;
extern int lag_fi3,flag_fi4,flag_fi5,flag_fsn,extract,flag_e,
cirf,noequ,flag_pln,nuz,flag_sp,flag_cL,st_pln,flag_nul,
flag_mat,nnod,fl_2m,fl_ord;
extern float range2,range3,range4,range5,
             range2p,range3p,range4p,range5p;

void binvec(int next, char *sn, char *on, 
            char nus[16][5],int *s)
{
if (next==0) 
{
*s=0;
on[0]=0; 
} 
else
if (next==1)
{  //  0 1
// sn[0]=0; sn[1]=0;
on[0]=0; on[1]=1;
nus[1][0]=0;
*s=1;
}
else if (next==2)
{  //    00  01 10 11
sn[0]=0; sn[1]=0; sn[2]=1; sn[3]=1; 
on[0]=0; on[1]=1; on[2]=1; on[3]=2;
nus[1][0]=0;
nus[2][0]=1;
nus[3][0]=0; nus[3][1]=1;
*s=3;
}
else if (next==3)
{  // 000 001 010 011 100 101 110 111
 sn[0]=0; sn[1]=0; sn[2]=1; sn[3]=1; sn[4]=0; sn[5]=0; sn[6]=1; sn[7]=1; 
 on[0]=0; on[1]=1; on[2]=1; on[3]=2; on[4]=1; on[5]=2; on[6]=2; on[7]=3;
 nus[1][0]=0;
 nus[2][0]=1;
 nus[3][0]=0; nus[3][1]=1;
 nus[4][0]=2;
 nus[5][0]=0; nus[5][1]=2;
 nus[6][0]=1; nus[6][1]=2;
 nus[7][0]=0; nus[7][1]=1; nus[7][2]=2;
 *s=7;
}
else if (next==4)
{  // 0000 0001 0010 0011 0100 0101 0110 0111 
    // 1000 1001 1010 1011 1100 1101 1110 1111
sn[0]=0; sn[1]=0; sn[2]=1; sn[3]=1; sn[4]=0; sn[5]=0; sn[6]=1; sn[7]=1; 
sn[8]=1; sn[9]=1; sn[10]=0; sn[11]=0; sn[12]=1; sn[13]=1; sn[14]=0; sn[15]=0; 
on[0]=0; on[1]=1; on[2]=1; on[3]=2; on[4]=1; on[5]=2; on[6]=2; on[7]=3; on[8]=1;
on[9]=2; on[10]=2; on[11]=3; on[12]=2; on[13]=3; on[14]=3; on[15]=4;
nus[1][0]=0;
nus[2][0]=1;
nus[3][0]=0; nus[3][1]=1;
nus[4][0]=2;
nus[5][0]=0; nus[5][1]=2;
nus[6][0]=1; nus[6][1]=2;
nus[7][0]=0; nus[7][1]=1; nus[7][2]=2;
nus[8][0]=3;
nus[9][0]=0; nus[9][1]=3;
nus[10][0]=1; nus[10][1]=3;
nus[11][0]=0; nus[11][1]=1; nus[11][2]=3;
nus[12][0]=2; nus[12][1]=3;
nus[13][0]=0; nus[13][1]=2; nus[13][2]=3;
nus[14][0]=1; nus[14][1]=2; nus[14][2]=3;
nus[15][0]=0; nus[15][1]=1; nus[15][2]=2; nus[15][3]=3;
*s=15;
}
}
/***************************************************************/
void form(int num,char *ac,int next,char *ext,
                 int n,PASSIVE *matr,int act_n,SOURCE *act)
{
 int i,k,l,m,s,n1,act1_n,n2,act2_n,nn[2],na[2],na0,na1;
 unsigned long dl,dl1,dl2;
 char nus[16][5],sn[16],on[16],o,o1,o2;
 PASSIVE *a1;
 SOURCE *act1;
// 01
// 01
//    00  - 0      ─хё Єшўэюх яЁхфёЄртыхэшх ─┬ ЁрчьхЁэюёЄш n=2 (2n=4)
//    01  - 1
//    10  - 2
//    11  - 3
//    011223      ─хё Єшўэ√щ тхъЄюЁ яхЁтющ яюыютшэ√ ─┬
//    012123      ─хё Єшўэ√щ тхъЄюЁ тЄюЁющ яюыютшэ√ ─┬
//    00 11 12 21 22 33   ╧ю¤ыхьхэЄэюх юс·хфшэхэшх фхё Єшўэ√ї тхъЄюЁют яюфёїхь√ I
//    33 22 21 12 11 00   ─юяюыэ ■∙хх юЄюсЁрцхэшх яю ЇюЁьєых: u-1-i (u-1-j) фы  яюфёїхь√ II
 o=ext[next];
 binvec(next,sn,on,nus,&s);
 fortest(num,n,ac,&n1,&act1_n,0);
 fortest(num,n,ac,&n2,&act2_n,1);
 dl2=leng+strlen(c);
 nn[0]=n1; nn[1]=n2; na0=act1_n; na1=act2_n;
 for (k=0;k<=s;k++)
      for (l=0;l<=s;l++)
		if(on[k]==on[l]) 
		{
		na[0]=na0+on[k]; na[1]=na1+next-on[k];		
        m=0;
rep:
 a1=( PASSIVE *) malloc (nn[m]*sizeof(PASSIVE));
	act1=( SOURCE *) malloc (na[m]*sizeof(SOURCE));
     fortrans(num,ac,n,matr,act,a1,act1,m);	
    if (!m)
				{   act1_n=na0;
					dl=leng+strlen(c);
					if (dl != dl2) strcat(b++,"+");
	         if (k && next) {
				o1=ext[nus[k][0]]; o2=ext[nus[l][0]];
				if (next>1 && (sn[k]+sn[l])%2) 
				  nuiall(o,o1,o2,o,&act1_n,act1);
                else nuiall(o,o1,o,o2,&act1_n,act1);
                for (i=1;i<on[k];i++) 
	nuiall(o,ext[nus[k][i]],o,ext[nus[l][i]],&act1_n,act1);	
    				}   
                 }
				else 
				{act1_n=na1;
			     strcat(b++,"*");
	 		     for(i=0;i<next-on[k];i++)	
		nuiall(o,ext[nus[s-k][i]],o,ext[nus[s-l][i]],&act1_n,act1);	
	}
	    if (n1+act1_n > 1) strcat(b++,"(");
	    dl1=leng+strlen(c);
 	    gggf(nn[m],a1,na[m],act1);
	    free(a1); free(act1);
	    if (dl1>=leng+strlen(c)) {ster(dl); continue;}
	    else if (n1+act1_n > 1) strcat(b++,")");
        if (   *(b-4) == '*' && *(b-3) == '(' && *(b-2) == '1'
            && *(b-1) == ')' )
          {b-=4; b[0]='\0';}			
		if (!m) {m++; goto rep;} 
		kontrol();
	   }
 frematr(n,matr); freact(act_n,act);
}
/********************************************************/
 void formp(int num,char *ac,int next,char *ext,
	 int n,PASSIVE *matr,int act_n,SOURCE *act,int pln,int st_pln)
{
 char nus[16][5],sn[16],on[16],o,o1,o2;
 int i,k,l,m,n1,n2,nn[2],na[2],act1_n,act2_n,pln0,
     pln1,pln2,st,st0,st1,s,cL1,cL2,key[2],na0,na1;
 unsigned long dl,dl1,dl2;
 PASSIVE *a1;  
 SOURCE *act1;
 
 o=ext[next];
 st=st_pln-pln; 
 n1=n_cL(n,matr,ac,0);
 n2=n_cL(n,matr,ac,1);
 if (n1 < n2 && fl_ord==1 || n1 > n2 && fl_ord==2 || !fl_ord) 
 {key[0]=0; key[1]=1; cL1=n1; cL2=n2;}
 else 
 {key[0]=1; key[1]=0; cL1=n2; cL2=n1;} 
 fortest(num,n,ac,&n1,&act1_n,key[0]);
 fortest(num,n,ac,&n2,&act2_n,key[1]);
 nn[0]=n1; nn[1]=n2; na0=act1_n; na1=act2_n; 
 dl2=leng+strlen(c);
 binvec(next,sn,on,nus,&s);
 if (st<=cL2) st0=0; else st0=st-cL2;
 if (st<cL1) st1=st; else st1=cL1; 
 for (pln1=st0; pln1<=st1; pln1++)
   for (k=0;k<=s;k++)  
	for (l=0;l<=s;l++)
		if (on[k]==on[l])
		{
    na[0]=na0+on[k]; na[1]=na1+next-on[k]; 
	m=0;
rep:    a1=( PASSIVE *) malloc (nn[m]*sizeof(PASSIVE));
        act1=( SOURCE *) malloc (na[m]*sizeof(SOURCE));
	    fortrans(num,ac,n,matr,act,a1,act1,key[m]);				    
	              if (!m)
				{ act1_n=na0;
		pln0=pln1;
		dl=leng+strlen(c);
     if (dl != dl2) strcat(b++,"+");	
	    if (k && next) {
			o1=ext[nus[k][0]]; o2=ext[nus[l][0]];
			if (next>1 && (sn[k]+sn[l])%2) 
				  nuiall(o,o1,o2,o,&act1_n,act1);
                else nuiall(o,o1,o,o2,&act1_n,act1);
                for (i=1;i<on[k];i++) 
	nuiall(o,ext[nus[k][i]],o,ext[nus[l][i]],&act1_n,act1);	
			         	   	}
          	}
	else
    {	act1_n=na1;
     pln0=st-pln1; 
     strcat(b++,"*"); 	
	 for (i=0;i<next-on[k];i++)	
		nuiall(o,ext[nus[s-k][i]],o,ext[nus[s-l][i]],&act1_n,act1);
	}	
	if (n1+act1_n > 1) strcat(b++,"(");
     dl1=leng+strlen(c);
	 gggp(nn[m],a1,na[m],act1,0,pln0);
	 free(a1); free(act1);
	 if (dl1>=leng+strlen(c)) {ster(dl); continue;} 
     else if (n1+act1_n > 1) strcat(b++,")");
	 if (   *(b-4) == '*' && *(b-3) == '(' && *(b-2) == '1'
       && *(b-1) == ')' )
         {b-=4; b[0]='\0';}
	 kontrol();
	 if (!m) {m++; goto rep;} 
   }
 frematr(n,matr); freact(act_n,act);
}
/*  ********* ПОДПРОГРАММА УЧЕТА ШАРНИРОВ **********
    ************************************************ */
void bisec1p(int n,PASSIVE *matr,int act_n,SOURCE *act,
             int *flag,int pln,int st_pln)
{
 int i,k,p,num,next=0;
 char *str,*ac0,ext[1];
 GRAPH *ac=( GRAPH *) malloc ((n+2*act_n)*sizeof(GRAPH));

 str=(char*) malloc(2*n+4*act_n+1);
 ac0=(char*) malloc(n+2*act_n);
 nodestr(n,matr,act_n,act,&p,str);
 transfor(n,matr,act_n,act,&num,ac);
 for (k=0;k<num;k++) ac0[k]=ac[k].v1;
 *flag=0;
 for (i=0;i<p;i++)
    {
	// if (str[i]!='0') continue; 
     k=0;
     bond1(num,ac,str[i],str[i],&k);
     if (k == num || !indep1(num,ac)) goto rep;
 //    printf("*");
     *flag=1;
	 for (k=0;k<num;k++) ac0[k]=ac[k].v1;
	 formp(num,ac0,next,ext,n,matr,act_n,act,pln,st_pln);
     break;
rep:
     for (k=0;k<num;k++) ac[k].v1=ac0[k];
    }
 free(ac0); free(ac); free(str);
}
/*  ********* ПОДПРОГРАММА УЧЕТА ШАРНИРОВ **********
    ************************************************ */
void bisec1(int n,PASSIVE *matr,int act_n,SOURCE *act,
             int *flag)
{
 int i,k,p,num,next=0;
 char *str,*ac0,ext[1];
 GRAPH *ac=( GRAPH *) malloc ((n+2*act_n)*sizeof(GRAPH));

 str=(char*) malloc(2*n+4*act_n+1);
 ac0=(char*) malloc(n+2*act_n);
 nodestr(n,matr,act_n,act,&p,str);
 transfor(n,matr,act_n,act,&num,ac);
 for (k=0;k<num;k++) ac0[k]=ac[k].v1;
 *flag=0;
 for (i=0;i<p;i++)
    {
     k=0;
     bond1(num,ac,str[i],str[i],&k);
     if (k == num || !indep1(num,ac)) goto rep;
 //    printf("*");
     *flag=1;
	 for (k=0;k<num;k++) ac0[k]=ac[k].v1;
     form(num,ac0,next,ext,n,matr,act_n,act);
     break;
rep:
     for (k=0;k<num;k++) ac[k].v1=ac0[k];
    }
 free(ac0); free(ac); free(str);
}
/*************************************************************/
 void bisec2p(int n,PASSIVE *matr,int act_n,SOURCE *act,
             int *flag,float range,int pln,int st_pln)
{
 int i,j,l,i3,j3,h,m,k,p,num,fl11,fl12,fl21,fl22,mn,mk,next=1;
 //,kopt;
 char s1,s2,s1o,s2o,*str,*ac0,o1,o2,ext[2];
 float pr1,rel,nrel,drel;

 GRAPH *ac=( GRAPH *) malloc ((n+2*act_n)*sizeof(GRAPH));
 GRAPH *ace=( GRAPH *) malloc ((n+2*act_n)*sizeof(GRAPH));
 ac0=(char*) malloc(n+2*act_n);
 str=(char*) malloc(2*n+4*act_n+1);
 nodestr(n,matr,act_n,act,&p,str);
 pr1 = 1.;
 transfor(n,matr,act_n,act,&num,ac);
 for (k=0;k<num;k++) ac0[k]=ac[k].v1;
 *flag=0;
 for (i=0;i<p;i++)
    {
     s1=str[i];
     for (j=i+1;j<p;j++)
        {
         s2=str[j];
     if (s1 != '0' && s2 != '0') continue; 
	 for (k=0;k<act_n;k++)
	    {
	     if (act[k].kol != 999) continue;
	     o1=act[k].v3;
	     o2=act[k].v4;
	     if (s1==o1 && s2==o2 || s1==o2 && s2==o1) goto cnt;
	    }
         k=0;
         bond2(num,ac,s1,s1,s2,&k);
         if (   k==num || k==1
             || !indep2(num,ac,ac0,s1,s2,&k)) goto rep;
         nrel=num-k; drel=num; rel=nrel/drel;
         if (k==num-1) goto rep;
         if (rel <= 0.5-range || rel >= 0.5+range) goto rep;
         for (i3=0;i3<n;i3++)
            {
	     if (   ac0[i3]==s1 && ac[i3].v2==s2
	         || ac0[i3]==s2 && ac[i3].v2==s1)
	       {
	        mn=0; mk=0;
	        for (j3=0;j3<n;j3++)
	           {
	            if (j3==i3) continue;
                    mn=mn+matr[j3].kol;
	            if (ac[j3].v1=='я') mk=mk+matr[j3].kol;
	           }
	        h=0;
	        for (j3=n;j3<num;j3++)
	           {
	            if (ac[j3].ess < 0) continue;
	            mn=mn+abs(act[ac[j3].ess-1].kol);
	            if (ac[j3].v1=='я')
		      {h++; mk=mk+abs(act[ac[j3].ess-1].kol);}
	           }
	        if (k-h<n+act_n-k+h || k-h==n+act_n-k+h && mk<mn-mk)
	          {if (ac[i3].v1 != 'я') {ac[i3].v1='я'; k++;}}
                else
	          if (ac[i3].v1 == 'я') {ac[i3].v1=ac0[i3]; k--;}
	       }
            }
         if (k==1) goto rep;
         fl11=0;fl12=0;fl21=0;fl22=0;
         for (m=0;m<num;m++)
            {
	     if (ac[m].v1 =='я')
	       {
	        if (ac0[m]==s1 || ac[m].v2==s1) fl11=1;
	        if (ac0[m]==s2 || ac[m].v2==s2) fl12=1;
	       }
	     else
	       {
	        if (ac0[m]==s1 || ac[m].v2==s1) fl21=1;
	        if (ac0[m]==s2 || ac[m].v2==s2) fl22=1;
	       }
             if (fl11 && fl12 && fl21 && fl22) goto net;
            }
         goto rep;
     net:
         if (pr1 > fabs(rel-0.5))
	   {
	    *flag=1;
            for (l=0;l<num;l++) ace[l]=ac[l];
	    s1o=s1; s2o=s2; // kopt=k;
            pr1=fabs(rel-0.5);
           }
    rep:
	 for (l=0;l<num;l++) ac[l].v1=ac0[l];
    cnt: ;
        }
     }
 free(str); free(ac);
 if (*flag)
   {
	   for (l=0;l<num;l++) ac0[l]=ace[l].v1;
       if (s2o=='0') {ext[0]=s1o; ext[1]=s2o;}
       else {ext[0]=s2o; ext[1]=s1o;}
 	   formp(num,ac0,next,ext,n,matr,act_n,act,pln,st_pln);
   }	 
 free(ace); free(ac0);
}
/*************************************************************/
 void bisec2(int n,PASSIVE *matr,int act_n,SOURCE *act,
             int *flag,float range)
{
 int i,j,l,i3,j3,h,m,k,p,num,fl11,fl12,fl21,fl22,mn,mk,next=1;
 char s1,s2,s1o,s2o,*str,*ac0,o1,o2,ext[2];
 float pr1,rel,nrel,drel;

 GRAPH *ac=( GRAPH *) malloc ((n+2*act_n)*sizeof(GRAPH));
 GRAPH *ace=( GRAPH *) malloc ((n+2*act_n)*sizeof(GRAPH));
 ac0=(char*) malloc(n+2*act_n);
 str=(char*) malloc(2*n+4*act_n+1);
 nodestr(n,matr,act_n,act,&p,str);
 pr1 = 1.;
 transfor(n,matr,act_n,act,&num,ac);
 for (k=0;k<num;k++) ac0[k]=ac[k].v1;
 *flag=0;
 for (i=0;i<p;i++)
    {
     s1=str[i];
     for (j=i+1;j<p;j++)
        {
         s2=str[j];
    // if (s1 != '0' && s2 != '0') continue;
	 for (k=0;k<act_n;k++)
	    {
	     if (act[k].kol != 999) continue;
	     o1=act[k].v3;
	     o2=act[k].v4;
	     if (s1==o1 && s2==o2 || s1==o2 && s2==o1) goto cnt;
	    }
         k=0;
         bond2(num,ac,s1,s1,s2,&k);
         if (   k==num || k==1
             || !indep2(num,ac,ac0,s1,s2,&k)) goto rep;
         nrel=num-k; drel=num; rel=nrel/drel;
         if (k==num-1) goto rep;
         if (rel <= 0.5-range || rel >= 0.5+range) goto rep;
         for (i3=0;i3<n;i3++)
            {
	     if (   ac0[i3]==s1 && ac[i3].v2==s2
	         || ac0[i3]==s2 && ac[i3].v2==s1)
	       {
	        mn=0; mk=0;
	        for (j3=0;j3<n;j3++)
	           {
	            if (j3==i3) continue;
                    mn=mn+matr[j3].kol;
	            if (ac[j3].v1=='я') mk=mk+matr[j3].kol;
	           }
	        h=0;
	        for (j3=n;j3<num;j3++)
	           {
	            if (ac[j3].ess < 0) continue;
	            mn=mn+abs(act[ac[j3].ess-1].kol);
	            if (ac[j3].v1=='я')
		      {h++; mk=mk+abs(act[ac[j3].ess-1].kol);}
	           }
	        if (k-h<n+act_n-k+h || k-h==n+act_n-k+h && mk<mn-mk)
	          {if (ac[i3].v1 != 'я') {ac[i3].v1='я'; k++;}}
                else
	          if (ac[i3].v1 == 'я') {ac[i3].v1=ac0[i3]; k--;}
	       }
            }
         if (k==1) goto rep;
         fl11=0;fl12=0;fl21=0;fl22=0;
         for (m=0;m<num;m++)
            {
	     if (ac[m].v1 =='я')
	       {
	        if (ac0[m]==s1 || ac[m].v2==s1) fl11=1;
	        if (ac0[m]==s2 || ac[m].v2==s2) fl12=1;
	       }
	     else
	       {
	        if (ac0[m]==s1 || ac[m].v2==s1) fl21=1;
	        if (ac0[m]==s2 || ac[m].v2==s2) fl22=1;
	       }
             if (fl11 && fl12 && fl21 && fl22) goto net;
            }
         goto rep;
     net:
         if (pr1 > fabs(rel-0.5))
	   {
	    *flag=1;
            for (l=0;l<num;l++) ace[l]=ac[l];
	    s1o=s1; s2o=s2;
            pr1=fabs(rel-0.5);
           }
    rep:
	 for (l=0;l<num;l++) ac[l].v1=ac0[l];
    cnt: ;
        }
     }
 free(str); free(ac);
 if (*flag)
   {
	    for (l=0;l<num;l++) ac0[l]=ace[l].v1;
	   if (s2o=='0') {ext[0]=s1o; ext[1]=s2o;}
       else {ext[0]=s2o; ext[1]=s1o;}
	   form(num,ac0,next,ext,n,matr,act_n,act);
   }
 free(ace); free(ac0);
}
/**************************************************/
void transfor(int n,PASSIVE *matr,
	      int act_n,SOURCE *act,int *num,GRAPH *ac)
{
 int i;

 *num=0;
 for (i=0;i<n;i++)
    {
     ac[*num].v1=matr[*num].v1;
     ac[*num].v2=matr[*num].v2;
     ac[(*num)++].ess=0;
    }
 for (i=0;i<act_n;i++)
    {
     ac[*num].v1=act[i].v1;
     ac[*num].v2=act[i].v2;
     ac[*num].ess=i+1;
     ++(*num);
     ac[*num].v1=act[i].v3;
     ac[*num].v2=act[i].v4;
     ac[(*num)++].ess=-(i+1);
    }
}
/***************************************************************/
int n_cL(int n,PASSIVE *matr,char *ac,int key)
{
 int i,nc;
 
 nc=0;
 for (i=0;i<n;i++)
    if ((!key && ac[i] == 'я' || key && ac[i] != 'я')
		 && (matr[i].reb[0] == 'c' || matr[i].reb[0] == 'L')) nc++; 
 return(nc);
}
/***************************************************************/
void fortrans(int num,char *ac,int n,PASSIVE *matr,SOURCE *act,
	      PASSIVE *a1,SOURCE *act1,int key)
{
 int i,k=0,n1=0,act1_n=0;
 for (i=0;i<n;i++)
     if (!key && ac[i] == 'я' || key && ac[i] != 'я')
       {
	      a1[n1].v1=matr[i].v1;
		  a1[n1].v2=matr[i].v2;
		  a1[n1].kol=matr[i].kol;
		  a1[n1++].reb=strdup(matr[i].reb);
	    }
	for (i=n;i<num;i+=2)
	{
    	if (!key && ac[i] == 'я' || key && ac[i] != 'я')
	  {  
          act1[act1_n].v1=act[k].v1;
	      act1[act1_n].v2=act[k].v2;
	      act1[act1_n].v3=act[k].v3;
	      act1[act1_n].v4=act[k].v4;
	      act1[act1_n].kol=act[k].kol;
	      act1[act1_n++].reb=strdup(act[k].reb);
	  }
	 k++;
	}
}
/**********************************************************/
void fortest(int num,int n,char *ac,int *n1,int *act1_n,int key)
{
 int i;

 *n1=0, *act1_n=0;
 for (i=0;i<n;i++)
     if (!key && ac[i] == 'я' || key && ac[i] != 'я') (*n1)++;
 for (i=n;i<num;i+=2)
    if (!key && ac[i] == 'я' || key && ac[i] != 'я') (*act1_n)++;
}
/******************************************************/
void bond1(int n,GRAPH *a1,char s1,char s2,int *k)
{
 int i;
 char s4;

 for (i=0;i<n;i++)
    {
     if (*k ==n) break;
     if (a1[i].v1==s1)
       {
	if (a1[i].v2==s2)
	  {
	   a1[i].v1='я'; (*k)++;
	   continue;
	  }
	s4=a1[i].v2;
	a1[i].v1='я'; (*k)++;
	bond1(n,a1,s4,s2,&*k);
	if (s1==s2) break;
       }
     if (a1[i].v1!='я' && a1[i].v2==s1)
       {
	if (a1[i].v1==s2)
	  {
	   a1[i].v1='я'; (*k)++;
	   continue;
	  }
	s4=a1[i].v1;
	a1[i].v1='я'; (*k)++;
	bond1(n,a1,s4,s2,&*k);
	if (s1==s2) break;
       }
    }
}
/*********************************************************/
int indep1(int nc,GRAPH *ac)
{
 int i,j;

 for (i=0;i<nc;i++)
    {
     if (!ac[i].ess || ac[i].v1 == 'я') continue;
     for (j=0;j<nc;j++)
	{
	 if (!ac[j].ess || ac[j].v1 == 'я') continue;
	 if (ac[i].ess==-ac[j].ess) goto eq1;
	}
     return(0);
eq1: ;
    }
 return(1);
}
/***********************************************************/
void bond2(int n,GRAPH *a1,char s1,char s2,char s3,int *k)
{
 int i;
 char s4;

 for (i=0;i<n;i++)
    {
     if (*k > n-2)
       {
	*k=n;
	break;
       }
     if (a1[i].v1==s1)
       {
	if (a1[i].v2==s2 || a1[i].v2==s3)
	  {
	   a1[i].v1='я'; (*k)++;
	   continue;
	  }
	s4=a1[i].v2;
	a1[i].v1='я'; (*k)++;
	bond2(n,a1,s4,s2,s3,&*k);
	if(s1==s2) break;
       }
     if (a1[i].v1!='я' && a1[i].v2==s1)
       {
	if (a1[i].v1==s2 || a1[i].v1==s3)
	  {
	   a1[i].v1='я'; (*k)++;
	   continue;
	  }
	s4=a1[i].v1;
	a1[i].v1='я'; (*k)++;
	bond2(n,a1,s4,s2,s3,&*k);
	if(s1==s2) break;
       }
    }
}
/*********************************************************/
int indep2(int nc,GRAPH *ac,
	   char *ac0,char s1,char s2,int *k)
{
 int i,j;

 for (i=0;i<nc;i++)
    {
     if (!ac[i].ess || ac[i].v1 == 'я') continue;
     for (j=0;j<nc;j++)
	{
	 if (!ac[j].ess || ac[j].v1 == 'я') continue;
	 if (ac[i].ess==-ac[j].ess) goto eq1;
	}
     if (ac0[i]==s1 && ac[i].v2==s2 || ac[i].v2==s1 && ac0[i]==s2)
       {
	ac[i].v1='я';
	(*k)++;
       }
eq1: ;
    }
 for (i=0;i<nc;i++)
    {
     if (!ac[i].ess || ac[i].v1 != 'я') continue;
     for (j=0;j<nc;j++)
	{
	 if (!ac[j].ess || ac[j].v1 != 'я') continue;
	 if (ac[i].ess==-ac[j].ess) goto eq2;
	}
     if (ac0[i]==s1 && ac[i].v2==s2 || ac[i].v2==s1 && ac0[i]==s2)
       {
	ac[i].v1=ac0[i];
	(*k)--;
       }
     else return(0);
eq2: ;
    }
 return(1);
}
/****************************************************/
void bond3(int n,GRAPH *a1,char s1,char s2,char s3,char s4,int *k)
{
 int i;
 char s5;

 for(i=0;i<n;i++)
  {
   if (a1[i].v1==s1)
     {
      if (a1[i].v2==s2 || a1[i].v2==s3 || a1[i].v2==s4)
	{
	 a1[i].v1='я'; (*k)++;
	 continue;
	}
      s5=a1[i].v2;
      a1[i].v1='я'; (*k)++;
      bond3(n,a1,s5,s2,s3,s4,k);
      if(s1==s2) break;
     }
   if (a1[i].v1 != 'я' && a1[i].v2==s1)
     {
      if (a1[i].v1==s2 || a1[i].v1==s3 || a1[i].v1==s4)
	{
	 a1[i].v1='я'; (*k)++;
	 continue;
	}
      s5=a1[i].v1;
      a1[i].v1='я'; (*k)++;
      bond3(n,a1,s5,s2,s3,s4,k);
      if (s1==s2) break;
     }
  }
}
/*********************************************************/
int indep3(int nc,GRAPH *ac,
	   char *ac0,char s1,char s2,char s3,int *k)
{
 int i,j;

 for (i=0;i<nc;i++)
    {
     if (!ac[i].ess || ac[i].v1 == 'я') continue;
     for (j=0;j<nc;j++)
	{
	 if (!ac[j].ess || ac[j].v1 == 'я') continue;
	 if (ac[i].ess==-ac[j].ess) goto eq1;
	}
     if (   ac0[i]==s1 && ac[i].v2==s2 || ac[i].v2==s1 && ac0[i]==s2
	 || ac0[i]==s1 && ac[i].v2==s3 || ac[i].v2==s1 && ac0[i]==s3
	 || ac0[i]==s2 && ac[i].v2==s3 || ac[i].v2==s2 && ac0[i]==s3)
       {
	ac[i].v1='я';
	(*k)++;
       }
eq1: ;
    }
 for (i=0;i<nc;i++)
    {
     if (!ac[i].ess || ac[i].v1 != 'я') continue;
     for (j=0;j<nc;j++)
	{
	 if (!ac[j].ess || ac[j].v1 != 'я') continue;
	 if (ac[i].ess==-ac[j].ess) goto eq2;
	}
     if (   ac0[i]==s1 && ac[i].v2==s2 || ac[i].v2==s1 && ac0[i]==s2
	 || ac0[i]==s1 && ac[i].v2==s3 || ac[i].v2==s1 && ac0[i]==s3
	 || ac0[i]==s2 && ac[i].v2==s3 || ac[i].v2==s2 && ac0[i]==s3)
       {
	ac[i].v1=ac0[i];
	(*k)--;
       }
     else return(0);
eq2: ;
    }
 return(1);
}
/************************************************************/
void bisec3p(int n,PASSIVE *matr,int act_n,SOURCE *act,
             int *flag,float range,int pln,int st_pln)
{
 int i,j,h=0,k,k1,k2,k3,p,num,numa,next,r1,r2,r3,
     fl11,fl12,fl13,fl21,fl22,fl23; // ,kopt;
 char s1,s2,s3,s1o,s2o,*str,*ac0,ext[3],o1,o2;
 float pr1=1.,rel,nrel,drel;
 GRAPH *ac=( GRAPH *) malloc ((n+2*act_n+100)*sizeof(GRAPH));
 GRAPH *ace=( GRAPH *) malloc ((n+2*act_n+100)*sizeof(GRAPH));

 ac0=(char*) malloc(n+2*act_n+100);
 str=(char*) malloc(2*n+4*act_n+1);
 *flag=0;
 nodestr(n,matr,act_n,act,&p,str);
 transfor(n,matr,act_n,act,&num,ac);
 if (flag_mat) {
 i=0;
 for (k=num;k<num+p-1;k++)
    {
     ac[k].ess=0;
     ac[k].v1=str[i];
     ac[k].v2=str[i+1];
     i++;
    }
  numa=num+p-1;} else numa=num;
 for (i=0;i<p;i++)
    {
  if (str[i] != '0')
{
     k=0;
     for (j=0;j<num;j++)
	{
	 if (ac[j].v1 == str[i] || ac[j].v2 == str[i]) k++;
	}
     if (k > h) {h=k; s3=str[i];}  } else {s3=str[i]; break;}
    }
 if (!flag_nul && s3 != '0') goto ret;
 for (i=0;i<p;i++)
    {
     if (str[i]==s3) {for (j=i;j<p-1;j++) str[j]=str[j+1]; break;}
    }
 for (k=0;k<numa;k++) ac0[k]=ac[k].v1;
 for (i=0;i<p-1;i++)
    {
     s1=str[i];
     for (j=i+1;j<p-1;j++)
	{
                    s2=str[j];
	 k1=0;
	 bond3(numa,ac,s1,s1,s2,s3,&k1);
         for (h=0;h<numa;h++) ac[h].v1=ac0[h];
         k2=0;
         bond3(numa,ac,s2,s2,s1,s3,&k2);
         for (h=0;h<numa;h++) ac[h].v1=ac0[h];
         k3=0;
         bond3(numa,ac,s3,s3,s1,s2,&k3);
         r1=abs(num-2*k1);
         r2=abs(num-2*k2);
         r3=abs(num-2*k3);
         k=0;
         for (h=0;h<numa;h++) ac[h].v1=ac0[h];
         if (r1<=r2 && r1<=r3) bond3(numa,ac,s1,s1,s2,s3,&k);
         else
           if (r2<=r1 && r2<=r3) bond3(numa,ac,s2,s2,s1,s3,&k);
           else
             if (r3<=r1 && r3<=r2) bond3(numa,ac,s3,s3,s1,s2,&k);
         if (!indep3(num,ac,ac0,s1,s2,s3,&k)) goto rep;
         nrel=num-k; drel=num; rel=nrel/drel;
         if (rel < 0.5-range || rel > 0.5+range) goto rep;
         fl11=0;fl12=0;fl13=0;fl21=0;fl22=0;fl23=0;
         for (h=0;h<num;h++)
            {
	     if (ac[h].v1 =='я')
	       {
	        if (ac0[h]==s1 || ac[h].v2==s1) fl11=1;
	        if (ac0[h]==s2 || ac[h].v2==s2) fl12=1;
                if (ac0[h]==s3 || ac[h].v2==s3) fl13=1;
	       }
	     else
	       {
	        if (ac0[h]==s1 || ac[h].v2==s1) fl21=1;
	        if (ac0[h]==s2 || ac[h].v2==s2) fl22=1;
	        if (ac0[h]==s3 || ac[h].v2==s3) fl23=1;
	       }
             if (   fl11 && fl12 && fl13
                 && fl21 && fl22 && fl23) goto net;
            }
         goto rep;
     net:
         if (pr1 > fabs(rel-0.5))
           {
            s1o=s1; s2o=s2;
	    pr1=fabs(rel-0.5); // kopt=k;
            for (h=0;h<numa;h++) ace[h]=ac[h];
            *flag=1;
           }
     rep:
	 for (h=0;h<numa;h++) ac[h].v1=ac0[h];
     cnt: ;
	}
    }
 if (*flag)
   {
   ext[0]=s1o; ext[1]=s2o; ext[2]=s3; next=2;
 //   printf(" ' ");
    for (k=0;k<num;k++) ac0[k]=ace[k].v1;
    formp(num,ac0,next,ext,n,matr,act_n,act,pln,st_pln);
   }
ret:
 free(str); free(ac);
 free(ace); free(ac0);
}
/************************************************************/
void bisec3(int n,PASSIVE *matr,int act_n,SOURCE *act,
            int *flag,float range)
{
 int i,j,h=0,k,k1,k2,k3,p,num,numa,next,r1,r2,r3,
     fl11,fl12,fl13,fl21,fl22,fl23;
 char s1,s2,s3,s1o,s2o,*str,*ac0,ext[3],o1,o2;
 float pr1=1.,rel,nrel,drel;
 GRAPH *ac=( GRAPH *) malloc ((n+2*act_n+100)*sizeof(GRAPH));
 GRAPH *ace=( GRAPH *) malloc ((n+2*act_n+100)*sizeof(GRAPH));

 ac0=(char*) malloc(n+2*act_n+100);
 str=(char*) malloc(2*n+4*act_n+1);
 *flag=0;
 nodestr(n,matr,act_n,act,&p,str);
 transfor(n,matr,act_n,act,&num,ac);
 if (flag_mat) {
 i=0;
 for (k=num;k<num+p-1;k++)
    {
     ac[k].ess=0;
     ac[k].v1=str[i];
     ac[k].v2=str[i+1];
     i++;
    }
  numa=num+p-1;} else numa=num;
 for (i=0;i<p;i++)
    {
     if (str[i] != '0')
{
     k=0;
     for (j=0;j<num;j++)
	{
	 if (ac[j].v1 == str[i] || ac[j].v2 == str[i]) k++;
	}
     if (k > h) {h=k; s3=str[i];}  } else {s3=str[i]; break;}
    }
 if (!flag_nul && s3 != '0') goto ret;
 for (i=0;i<p;i++)
    {
     if (str[i]==s3) {for (j=i;j<p-1;j++) str[j]=str[j+1]; break;}
    }
 for (k=0;k<numa;k++) ac0[k]=ac[k].v1;
 for (i=0;i<p-1;i++)
    {
     s1=str[i];
     for (j=i+1;j<p-1;j++)
	{
                    s2=str[j];
	/*     for (k=0;k<act_n;k++)
	            {
		 if (act[k].kol != 999) continue;
		 o1=act[k].v3; o2=act[k].v4;
		 if (   s1==o1 && s2==o2 || s1==o2 && s2==o1
		     || s1==o1 && s3==o2 || s1==o2 && s3==o1
		     || s2==o1 && s3==o2 || s2==o2 && s3==o1) goto cnt;
		}*/
	 k1=0;
	 bond3(numa,ac,s1,s1,s2,s3,&k1);
         for (h=0;h<numa;h++) ac[h].v1=ac0[h];
         k2=0;
         bond3(numa,ac,s2,s2,s1,s3,&k2);
         for (h=0;h<numa;h++) ac[h].v1=ac0[h];
         k3=0;
         bond3(numa,ac,s3,s3,s1,s2,&k3);
         r1=abs(num-2*k1);
         r2=abs(num-2*k2);
         r3=abs(num-2*k3);
         k=0;
         for (h=0;h<numa;h++) ac[h].v1=ac0[h];
         if (r1<=r2 && r1<=r3) bond3(numa,ac,s1,s1,s2,s3,&k);
         else
           if (r2<=r1 && r2<=r3) bond3(numa,ac,s2,s2,s1,s3,&k);
           else
             if (r3<=r1 && r3<=r2) bond3(numa,ac,s3,s3,s1,s2,&k);
         if (!indep3(num,ac,ac0,s1,s2,s3,&k)) goto rep;
         nrel=num-k; drel=num; rel=nrel/drel;
         if (rel < 0.5-range || rel > 0.5+range) goto rep;
         fl11=0;fl12=0;fl13=0;fl21=0;fl22=0;fl23=0;
         for (h=0;h<num;h++)
            {
	     if (ac[h].v1 =='я')
	       {
	        if (ac0[h]==s1 || ac[h].v2==s1) fl11=1;
	        if (ac0[h]==s2 || ac[h].v2==s2) fl12=1;
                if (ac0[h]==s3 || ac[h].v2==s3) fl13=1;
	       }
	     else
	       {
	        if (ac0[h]==s1 || ac[h].v2==s1) fl21=1;
	        if (ac0[h]==s2 || ac[h].v2==s2) fl22=1;
	        if (ac0[h]==s3 || ac[h].v2==s3) fl23=1;
	       }
             if (   fl11 && fl12 && fl13
                 && fl21 && fl22 && fl23) goto net;
            }
         goto rep;
     net:
         if (pr1 > fabs(rel-0.5))
           {
            s1o=s1; s2o=s2;
	    pr1=fabs(rel-0.5);
            for (h=0;h<numa;h++) ace[h]=ac[h];
            *flag=1;
           }
     rep:
	 for (h=0;h<numa;h++) ac[h].v1=ac0[h];
     cnt: ;
	}
    }
 if (*flag)
   {
    ext[0]=s1o; ext[1]=s2o; ext[2]=s3; next=2;
 //   printf(" ' ");
    for (k=0;k<num;k++) ac0[k]=ace[k].v1;
    form(num,ac0,next,ext,n,matr,act_n,act);
   }
ret:
 free(str); free(ac);
 free(ace); free(ac0);
}
/*****************************************************************/
void bond4(int n,GRAPH *a1,
	   char s1,char s2,char s3,char s4,char s5,int *k)
{
 int i;
 char s6;

 for(i=0;i<n;i++)
  {
   if (a1[i].v1==s1)
     {
      if (a1[i].v2==s2 || a1[i].v2==s3 || a1[i].v2==s4 || a1[i].v2==s5)
	{
	 a1[i].v1='я'; (*k)++;
	 continue;
	}
      s6=a1[i].v2;
      a1[i].v1='я'; (*k)++;
      bond4(n,a1,s6,s2,s3,s4,s5,&*k);
      if(s1==s2) break;
     }
   if (a1[i].v1 != 'я' && a1[i].v2==s1)
     {
      if (a1[i].v1==s2 || a1[i].v1==s3 || a1[i].v1==s4 || a1[i].v1==s5)
	{
	 a1[i].v1='я'; (*k)++;
	 continue;
	}
      s6=a1[i].v1;
      a1[i].v1='я'; (*k)++;
      bond4(n,a1,s6,s2,s3,s4,s5,&*k);
      if (s1==s2) break;
     }
  }
}
/*********************************************************/
int indep4(int nc,GRAPH *ac,
	   char *ac0,char s1,char s2,char s3,char s4,int *k)
{
 int i,j;

 for (i=0;i<nc;i++)
    {
     if (!ac[i].ess || ac[i].v1 == 'я') continue;
     for (j=0;j<nc;j++)
	{
	 if (!ac[j].ess || ac[j].v1 == 'я') continue;
	 if (ac[i].ess==-ac[j].ess) goto eq1;
	}
     if (   ac0[i]==s1 && ac[i].v2==s2 || ac[i].v2==s1 && ac0[i]==s2
	 || ac0[i]==s1 && ac[i].v2==s3 || ac[i].v2==s1 && ac0[i]==s3
	 || ac0[i]==s1 && ac[i].v2==s4 || ac[i].v2==s1 && ac0[i]==s4
	 || ac0[i]==s2 && ac[i].v2==s3 || ac[i].v2==s2 && ac0[i]==s3
	 || ac0[i]==s2 && ac[i].v2==s4 || ac[i].v2==s2 && ac0[i]==s4
	 || ac0[i]==s3 && ac[i].v2==s4 || ac[i].v2==s3 && ac0[i]==s4)
       {
	ac[i].v1='я';
	(*k)++;
       }
eq1: continue;
    }
 for (i=0;i<nc;i++)
    {
     if (!ac[i].ess || ac[i].v1 != 'я') continue;
     for (j=0;j<nc;j++)
	{
	 if (!ac[j].ess || ac[j].v1 != 'я') continue;
	 if (ac[i].ess==-ac[j].ess) goto eq2;
	}
     if (   ac0[i]==s1 && ac[i].v2==s2 || ac[i].v2==s1 && ac0[i]==s2
	 || ac0[i]==s1 && ac[i].v2==s3 || ac[i].v2==s1 && ac0[i]==s3
	 || ac0[i]==s1 && ac[i].v2==s4 || ac[i].v2==s1 && ac0[i]==s4
	 || ac0[i]==s2 && ac[i].v2==s3 || ac[i].v2==s2 && ac0[i]==s3
	 || ac0[i]==s2 && ac[i].v2==s4 || ac[i].v2==s2 && ac0[i]==s4
	 || ac0[i]==s3 && ac[i].v2==s4 || ac[i].v2==s3 && ac0[i]==s4)
       {
	ac[i].v1=ac0[i];
	(*k)--;
       }
     else return(0);
eq2: continue;
    }
 return(1);
}
/**************************************************************************/
void bisec4p(int n,PASSIVE *matr,int act_n,SOURCE *act,
             int *flag,float range,int pln,int st_pln)
{
 int i,j,l,h=0,k,k1,k2,k3,k4,p,num,numa,next,r1,r2,r3,r4;
 int fl11,fl12,fl13,fl14,fl21,fl22,fl23,fl24; // kopt;
 char s1,s2,s3,s4,s1o,s2o,s3o,*str,*ac0,ext[4],o1,o2;
 float pr1=1.,rel,nrel,drel;
 GRAPH *ac=( GRAPH *) malloc ((n+2*act_n+100)*sizeof(GRAPH));
 GRAPH *ace=( GRAPH *) malloc ((n+2*act_n+100)*sizeof(GRAPH));

 ac0=(char*) malloc(n+2*act_n+100);
 str=(char*) malloc(2*n+4*act_n+1);
 *flag=0;
 nodestr(n,matr,act_n,act,&p,str);
 transfor(n,matr,act_n,act,&num,ac);
 if (flag_mat) {
 i=0;
 for (k=num;k<num+p-1;k++)
    {
     ac[k].ess=0;
     ac[k].v1=str[i];
     ac[k].v2=str[i+1];
     i++;
    }
 numa=num+p-1;} else numa=num;
for (i=0;i<p;i++)
    {
    if (str[i] != '0')
{
     k=0;
     for (j=0;j<num;j++)
	{
	 if (ac[j].v1 == str[i] || ac[j].v2 == str[i]) k++;
	}
     if (k > h) {h=k; s4=str[i];}  } else {s4=str[i]; break;}
    }
 if (!flag_nul && s4 != '0') goto ret;
 for (i=0;i<p;i++)
    {
     if (str[i]==s4) {for (j=i;j<p-1;j++) str[j]=str[j+1]; break;}
    }
  for (k=0;k<numa;k++) ac0[k]=ac[k].v1;
  for (l=0;l<p-1;l++)
    {
     s1=str[l];
     for (i=l+1;i<p-1;i++)
	{
         s2=str[i];
	 for (j=i+1;j<p-1;j++)
	    {
     //    if (str[j]=='G') goto rep;
			s3=str[j];
			 k1=0;
	     bond4(numa,ac,s1,s1,s2,s3,s4,&k1);
             for (h=0;h<numa;h++) ac[h].v1=ac0[h];
             k2=0;
             bond4(numa,ac,s2,s2,s1,s3,s4,&k2);
             for (h=0;h<numa;h++) ac[h].v1=ac0[h];
             k3=0;
             bond4(numa,ac,s3,s3,s1,s2,s4,&k3);
			 for (h=0;h<numa;h++) ac[h].v1=ac0[h];
             k4=0;
             bond4(numa,ac,s4,s4,s1,s2,s3,&k4);
			 for (h=0;h<numa;h++) ac[h].v1=ac0[h];
             r1=abs(num-2*k1);
             r2=abs(num-2*k2);
             r3=abs(num-2*k3);
             r4=abs(num-2*k4);
             k=0;
             if (r1<=r2 && r1<=r3 && r1<=r4) bond4(numa,ac,s1,s1,s2,s3,s4,&k);
             else
             if (r2<=r1 && r2<=r3 && r2<=r4) bond4(numa,ac,s2,s2,s1,s3,s4,&k);
             else
             if (r3<=r1 && r3<=r2 && r3<=r4) bond4(numa,ac,s3,s3,s1,s2,s4,&k);
             else
             if (r4<=r1 && r4<=r2 && r4<=r3) bond4(numa,ac,s4,s4,s1,s2,s3,&k);
             if (!indep4(num,ac,ac0,s1,s2,s3,s4,&k)) goto rep;
             nrel=num-k; drel=num; rel=nrel/drel;

             if (rel > 0.5-range && rel < 0.5+range)
{
             fl11=0;fl12=0;fl13=0;fl14=0;fl21=0;fl22=0;fl23=0;fl24=0;
             for (h=0;h<num;h++)
                {
	         if (ac[h].v1 =='я')
		   {
	            if (ac0[h]==s1 || ac[h].v2==s1) fl11=1;
	            if (ac0[h]==s2 || ac[h].v2==s2) fl12=1;
	            if (ac0[h]==s3 || ac[h].v2==s3) fl13=1;
	            if (ac0[h]==s4 || ac[h].v2==s4) fl14=1;
	           }
	         else
	           {
	            if (ac0[h]==s1 || ac[h].v2==s1) fl21=1;
	            if (ac0[h]==s2 || ac[h].v2==s2) fl22=1;
	            if (ac0[h]==s3 || ac[h].v2==s3) fl23=1;
	            if (ac0[h]==s4 || ac[h].v2==s4) fl24=1;
	           }
                 if (   fl11 && fl12 && fl13 && fl14
                     && fl21 && fl22 && fl23 && fl24) goto net;
                }
}
             goto rep;
         net:
             if (pr1 > fabs(rel-0.5))
            {
                s1o=s1; s2o=s2; s3o=s3;
		pr1=fabs(rel-0.5);  // kopt=k;
                for (h=0;h<num;h++) ace[h]=ac[h];
                *flag=1;   
               }
	 rep:
	     for (h=0;h<numa;h++) ac[h].v1=ac0[h];
         cnt: ;
	    }
     }
  }
	if (*flag)
   {
    printf(" : ");
    ext[0]=s1o; ext[1]=s2o; ext[2]=s3o; ext[3]=s4; next=3;
	for (k=0;k<num;k++) ac0[k]=ace[k].v1;
	formp(num,ac0,next,ext,n,matr,act_n,act,pln,st_pln);
   }
ret:
 free(str); free(ac);
 free(ace); free(ac0);
}
/**************************************************************************/
void bisec4(int n,PASSIVE *matr,int act_n,SOURCE *act,int *flag,float range)
{

 int i,j,l,h=0,k,k1,k2,k3,k4,p,num,numa,next,r1,r2,r3,r4;
 int fl11,fl12,fl13,fl14,fl21,fl22,fl23,fl24;
 char s1,s2,s3,s4,s1o,s2o,s3o,*str,*ac0,ext[4],o1,o2;
 float pr1=1.,rel,nrel,drel;
 GRAPH *ac=( GRAPH *) malloc ((n+2*act_n+100)*sizeof(GRAPH));
 GRAPH *ace=( GRAPH *) malloc ((n+2*act_n+100)*sizeof(GRAPH));

 ac0=(char*) malloc(n+2*act_n+100);
 str=(char*) malloc(2*n+4*act_n+1);
 *flag=0;
 nodestr(n,matr,act_n,act,&p,str);
 transfor(n,matr,act_n,act,&num,ac);
 if (flag_mat) {
 i=0;
 for (k=num;k<num+p-1;k++)
    {
     ac[k].ess=0;
     ac[k].v1=str[i];
     ac[k].v2=str[i+1];
     i++;
    }
 numa=num+p-1;} else numa=num;
for (i=0;i<p;i++)
    {
     if (str[i] != '0') 
{
     k=0; 
     for (j=0;j<num;j++)
	{
	 if (ac[j].v1 == str[i] || ac[j].v2 == str[i]) k++;
	}
     if (k > h) {h=k; s4=str[i];}  } else {s4=str[i]; break;}
    }
 if (!flag_nul && s4 != '0') goto ret;
 for (i=0;i<p;i++)
    {
     if (str[i]==s4) {for (j=i;j<p-1;j++) str[j]=str[j+1]; break;}
    }
  for (k=0;k<numa;k++) ac0[k]=ac[k].v1;
  for (l=0;l<p-1;l++)
    {
     s1=str[l];
     for (i=l+1;i<p-1;i++)
	{
         s2=str[i];
	 for (j=i+1;j<p-1;j++)
	    {
             s3=str[j];
	     k1=0;
	     bond4(numa,ac,s1,s1,s2,s3,s4,&k1);
             for (h=0;h<numa;h++) ac[h].v1=ac0[h];
             k2=0;
             bond4(numa,ac,s2,s2,s1,s3,s4,&k2);
             for (h=0;h<numa;h++) ac[h].v1=ac0[h];
             k3=0;
             bond4(numa,ac,s3,s3,s1,s2,s4,&k3);
			 for (h=0;h<numa;h++) ac[h].v1=ac0[h];
             k4=0;
             bond4(numa,ac,s4,s4,s1,s2,s3,&k4);
			 for (h=0;h<numa;h++) ac[h].v1=ac0[h];
             r1=abs(num-2*k1);
             r2=abs(num-2*k2);
             r3=abs(num-2*k3);
             r4=abs(num-2*k4);
             k=0;
             if (r1<=r2 && r1<=r3 && r1<=r4) bond4(numa,ac,s1,s1,s2,s3,s4,&k);
             else
             if (r2<=r1 && r2<=r3 && r2<=r4) bond4(numa,ac,s2,s2,s1,s3,s4,&k);
             else
             if (r3<=r1 && r3<=r2 && r3<=r4) bond4(numa,ac,s3,s3,s1,s2,s4,&k);
             else
             if (r4<=r1 && r4<=r2 && r4<=r3) bond4(numa,ac,s4,s4,s1,s2,s3,&k);
             if (!indep4(num,ac,ac0,s1,s2,s3,s4,&k)) goto rep;
             nrel=num-k; drel=num; rel=nrel/drel;

             if (rel > 0.5-range && rel < 0.5+range) 
{ 
             fl11=0;fl12=0;fl13=0;fl14=0;fl21=0;fl22=0;fl23=0;fl24=0;
             for (h=0;h<num;h++)
                {
	         if (ac[h].v1 =='я')
		   {
	            if (ac0[h]==s1 || ac[h].v2==s1) fl11=1;
	            if (ac0[h]==s2 || ac[h].v2==s2) fl12=1;
	            if (ac0[h]==s3 || ac[h].v2==s3) fl13=1;
	            if (ac0[h]==s4 || ac[h].v2==s4) fl14=1;
	           }
	         else
	           {
	            if (ac0[h]==s1 || ac[h].v2==s1) fl21=1;
	            if (ac0[h]==s2 || ac[h].v2==s2) fl22=1;
	            if (ac0[h]==s3 || ac[h].v2==s3) fl23=1;
	            if (ac0[h]==s4 || ac[h].v2==s4) fl24=1;
	           }
                 if (   fl11 && fl12 && fl13 && fl14
                     && fl21 && fl22 && fl23 && fl24) goto net;
                }
}
             goto rep;
         net:
             if (pr1 > fabs(rel-0.5))
            {
                s1o=s1; s2o=s2; s3o=s3;
		pr1=fabs(rel-0.5);
                for (h=0;h<num;h++) ace[h]=ac[h];
                *flag=1;
               }
	 rep:
	     for (h=0;h<numa;h++) ac[h].v1=ac0[h];
         cnt: ;
	    }
	}
    }
 if (*flag)
   {
    printf(" % ");
    ext[0]=s1o; ext[1]=s2o; ext[2]=s3o; ext[3]=s4; next=3;
	for (k=0;k<num;k++) ac0[k]=ace[k].v1;
    form(num,ac0,next,ext,n,matr,act_n,act);
   }
ret:
 free(str); free(ac);
 free(ace); free(ac0);
}
/*****************************************************************/
void bond5(int n,GRAPH *a1,
	   char s1,char s2,char s3,char s4,char s5,char s6,int *k)
{
 int i;
 char s7;

 for(i=0;i<n;i++)
  {
   if (a1[i].v1==s1)
     {
      if (   a1[i].v2==s2 || a1[i].v2==s3
	  || a1[i].v2==s4 || a1[i].v2==s5 || a1[i].v2==s6)
	{
	 a1[i].v1='я'; (*k)++;
	 continue;
	}
      s7=a1[i].v2;
      a1[i].v1='я'; (*k)++;
      bond5(n,a1,s7,s2,s3,s4,s5,s6,k);
      if (s1==s2) break;
     }
   if (a1[i].v1 != 'я' && a1[i].v2==s1)
     {
      if (   a1[i].v1==s2 || a1[i].v1==s3
	  || a1[i].v1==s4 || a1[i].v1==s5 || a1[i].v1==s6)
	{
	 a1[i].v1='я'; (*k)++;
	 continue;
	}
      s7=a1[i].v1;
      a1[i].v1='я'; (*k)++;
      bond5(n,a1,s7,s2,s3,s4,s5,s6,k);
      if (s1==s2)  break;
     }
  }
}
/*********************************************************/
int indep5(int nc,GRAPH *ac,
	    char *ac0,char s1,char s2,char s3,char s4,char s5,int *k)
{
 int i,j;

 for (i=0;i<nc;i++)
    {
     if (!ac[i].ess || ac[i].v1 == 'я') continue;
     for (j=0;j<nc;j++)
	{
	 if (!ac[j].ess || ac[j].v1 == 'я') continue;
	 if (ac[i].ess==-ac[j].ess) goto eq1;
	}
     if (   ac0[i]==s1 && ac[i].v2==s2 || ac[i].v2==s1 && ac0[i]==s2
	 || ac0[i]==s1 && ac[i].v2==s3 || ac[i].v2==s1 && ac0[i]==s3
	 || ac0[i]==s1 && ac[i].v2==s4 || ac[i].v2==s1 && ac0[i]==s4
	 || ac0[i]==s1 && ac[i].v2==s5 || ac[i].v2==s1 && ac0[i]==s5
	 || ac0[i]==s2 && ac[i].v2==s3 || ac[i].v2==s2 && ac0[i]==s3
	 || ac0[i]==s2 && ac[i].v2==s4 || ac[i].v2==s2 && ac0[i]==s4
	 || ac0[i]==s2 && ac[i].v2==s5 || ac[i].v2==s2 && ac0[i]==s5
	 || ac0[i]==s3 && ac[i].v2==s4 || ac[i].v2==s3 && ac0[i]==s4
	 || ac0[i]==s3 && ac[i].v2==s5 || ac[i].v2==s3 && ac0[i]==s5
	 || ac0[i]==s4 && ac[i].v2==s5 || ac[i].v2==s4 && ac0[i]==s5)
       {
	ac[i].v1='я';
	(*k)++;
       }
eq1: continue;
    }
 for (i=0;i<nc;i++)
    {
     if (!ac[i].ess || ac[i].v1 != 'я') continue;
     for (j=0;j<nc;j++)
	{
	 if (!ac[j].ess || ac[j].v1 != 'я') continue;
	 if (ac[i].ess==-ac[j].ess) goto eq2;
	}
     if (   ac0[i]==s1 && ac[i].v2==s2 || ac[i].v2==s1 && ac0[i]==s2
	 || ac0[i]==s1 && ac[i].v2==s3 || ac[i].v2==s1 && ac0[i]==s3
	 || ac0[i]==s1 && ac[i].v2==s4 || ac[i].v2==s1 && ac0[i]==s4
	 || ac0[i]==s1 && ac[i].v2==s5 || ac[i].v2==s1 && ac0[i]==s5
	 || ac0[i]==s2 && ac[i].v2==s3 || ac[i].v2==s2 && ac0[i]==s3
	 || ac0[i]==s2 && ac[i].v2==s4 || ac[i].v2==s2 && ac0[i]==s4
	 || ac0[i]==s2 && ac[i].v2==s5 || ac[i].v2==s2 && ac0[i]==s5
	 || ac0[i]==s3 && ac[i].v2==s4 || ac[i].v2==s3 && ac0[i]==s4
	 || ac0[i]==s3 && ac[i].v2==s5 || ac[i].v2==s3 && ac0[i]==s5
	 || ac0[i]==s4 && ac[i].v2==s5 || ac[i].v2==s4 && ac0[i]==s5)
       {
	ac[i].v1=ac0[i];
	(*k)--;
       }
     else return(0);
eq2: continue;
    }
 return(1);
}
/***********************************************************************/
void bisec5p(int n,PASSIVE *matr,int act_n,SOURCE *act,
            int *flag,float range,int pln,int st_pln)
{
 int i,j,l,h=0,k,k1,k2,k3,k4,k5,m,p,num,numa,next,r1,r2,r3,r4,r5;
 int fl11,fl12,fl13,fl14,fl15,fl21,fl22,fl23,fl24,fl25; //,kopt;
 char c,s1,s2,s3,s4,s5,s1o,s2o,s3o,s4o,*str,*ac0,ext[5],o1,o2;
 float pr1=1.,rel,nrel,drel;
 GRAPH *ac=( GRAPH *) malloc ((n+2*act_n+100)*sizeof(GRAPH));
 GRAPH *ace=( GRAPH *) malloc ((n+2*act_n+100)*sizeof(GRAPH));

 ac0=(char*) malloc(n+2*act_n+100);
 str=(char*) malloc(2*n+4*act_n+1);
 *flag=0;
 nodestr(n,matr,act_n,act,&p,str);
 transfor(n,matr,act_n,act,&num,ac);
 if (flag_mat) {
 i=0;
 for (k=num;k<num+p-2;k++)
    {
     ac[k].ess=0;
     ac[k].v1=str[i];
     ac[k].v2=str[i+1];
     i++;
    }
 numa=num+p-1;} else numa=num;
 for (i=0;i<p;i++)
    {
     if (str[i] != '0') 
{
     k=0;
     for (j=0;j<num;j++)
	{
	 if (ac[j].v1 == str[i] || ac[j].v2 == str[i]) k++;
	}
     if (k > h) {h=k; s5=str[i];}  } else {s5=str[i]; break;}
    }
 if (!flag_nul && s5 != '0') goto ret;
 for (i=0;i<p;i++)
    {
     if (str[i]==s5) {for (j=i;j<p-1;j++) str[j]=str[j+1]; break;}
    }
 for (k=0;k<numa;k++) ac0[k]=ac[k].v1;
 for (l=0;l<p-1;l++)
    {
     s1=str[l];
     for (i=l+1;i<p-1;i++)
	{
	 s2=str[i];
	 for (j=i+1;j<p-1;j++)
	    {
	     s3=str[j];
	 for (m=j+1;m<p-1;m++)
	    {
	     s4=str[m];
/*	     for (k=0;k<act_n;k++)
		{
		 if (act[k].kol== 999) {
		 o1=act[k].v3;
		 o2=act[k].v4;
		 if (   s1==o1 && s2==o2 || s1==o2 && s2==o1
		     || s1==o1 && s3==o2 || s1==o2 && s3==o1
		     || s1==o1 && s4==o2 || s1==o2 && s4==o1
		     || s1==o1 && s5==o2 || s1==o2 && s5==o1
		     || s2==o1 && s3==o2 || s2==o2 && s3==o1
		     || s2==o1 && s4==o2 || s2==o2 && s4==o1
		     || s2==o1 && s5==o2 || s2==o2 && s5==o1
		     || s3==o1 && s4==o2 || s3==o2 && s4==o1
		     || s3==o1 && s5==o2 || s3==o2 && s5==o1
	     || s4==o1 && s5==o2 || s4==o2 && s5==o1) goto cnt;
		}  */
	     k1=0;
	     bond5(numa,ac,s1,s1,s2,s3,s4,s5,&k1);
             for (h=0;h<numa;h++) ac[h].v1=ac0[h];
             k2=0;
             bond5(numa,ac,s2,s2,s1,s3,s4,s5,&k2);
             for (h=0;h<num;h++) ac[h].v1=ac0[h];
             k3=0;
             bond5(numa,ac,s3,s3,s1,s2,s4,s5,&k3);
             for (h=0;h<numa;h++) ac[h].v1=ac0[h];
             k4=0;
             bond5(numa,ac,s4,s4,s1,s2,s3,s5,&k4);
             for (h=0;h<numa;h++) ac[h].v1=ac0[h];
             k5=0;
             bond5(numa,ac,s5,s5,s1,s2,s3,s4,&k5);
             r1=abs(num-2*k1);
             r2=abs(num-2*k2);
             r3=abs(num-2*k3);
             r4=abs(num-2*k4);
             r5=abs(num-2*k5);
             k=0;
             for (h=0;h<numa;h++) ac[h].v1=ac0[h];
           if (r1<=r2 && r1<=r3 && r1<=r4 && r1<=r5)
               bond5(numa,ac,s1,s1,s2,s3,s4,s5,&k);
             else
             if (r2<=r1 && r2<=r3 && r2<=r4 && r2<=r5)
               bond5(numa,ac,s2,s2,s1,s3,s4,s5,&k);
             else
             if (r3<=r1 && r3<=r2 && r3<=r4 && r3<=r5)
               bond5(numa,ac,s3,s3,s1,s2,s4,s5,&k);
             else
             if (r4<=r1 && r4<=r2 && r4<=r3 && r4<=r5)
               bond5(numa,ac,s4,s4,s1,s2,s3,s5,&k);
             else
             if (r5<=r1 && r5<=r2 && r5<=r3 && r5<=r4) 
               bond5(numa,ac,s5,s5,s1,s2,s3,s4,&k);
             if (!indep5(num,ac,ac0,s1,s2,s3,s4,s5,&k)) goto rep;
	     nrel=num-k; drel=num; rel=nrel/drel;
	     if (rel < 0.5-range || rel > 0.5+range) goto rep;
             fl11=0;fl12=0;fl13=0;fl14=0;fl15=0;
             fl21=0;fl22=0;fl23=0;fl24=0;fl25=0;
             for (h=0;h<num;h++)
                {
	         if (ac[h].v1 =='я')
	           {
	            if (ac0[h]==s1 || ac[h].v2==s1) fl11=1;
	            if (ac0[h]==s2 || ac[h].v2==s2) fl12=1;
	            if (ac0[h]==s3 || ac[h].v2==s3) fl13=1;
	            if (ac0[h]==s4 || ac[h].v2==s4) fl14=1;
	            if (ac0[h]==s5 || ac[h].v2==s5) fl15=1;
	           }
	         else
	           {
	            if (ac0[h]==s1 || ac[h].v2==s1) fl21=1;
	            if (ac0[h]==s2 || ac[h].v2==s2) fl22=1;
	            if (ac0[h]==s3 || ac[h].v2==s3) fl23=1;
	            if (ac0[h]==s4 || ac[h].v2==s4) fl24=1;
	            if (ac0[h]==s5 || ac[h].v2==s5) fl25=1;
	           }
                 if (   fl11 && fl12 && fl13 && fl14 && fl15
                     && fl21 && fl22 && fl23 && fl24 && fl25) goto net;
                }
             goto rep;
         net:
             if (pr1 > fabs(rel-0.5))
               {
                s1o=s1; s2o=s2; s3o=s3; s4o=s4;
		pr1=fabs(rel-0.5);   // kopt=k;
                for (h=0;h<numa;h++) ace[h]=ac[h];
                *flag=1;    goto fl;
               }
	 rep:
	     for (h=0;h<numa;h++) ac[h].v1=ac0[h];
         cnt: ;
	    }
	    }
        }
    }
fl:
 if (*flag)
   {
    printf(" ^ ");
    ext[0]=s1o; ext[1]=s2o; ext[2]=s3o; ext[3]=s4o; ext[4]='0';//s5; 
    next=4;
    for (k=0;k<num;k++) ac0[k]=ace[k].v1;
	formp(num,ac0,next,ext,n,matr,act_n,act,pln,st_pln);
   }
ret:
 free(str); free(ac);
 free(ace); free(ac0);
}
/***********************************************************************/
void bisec5(int n,PASSIVE *matr,int act_n,SOURCE *act,
            int *flag,float range)
{
 int i,j,l,h=0,k,k1,k2,k3,k4,k5,m,p,num,numa,next,r1,r2,r3,r4,r5;
 int fl11,fl12,fl13,fl14,fl15,fl21,fl22,fl23,fl24,fl25;
 char c,s1,s2,s3,s4,s5,s1o,s2o,s3o,s4o,*str,*ac0,ext[5],o1,o2;
 float pr1=1.,rel,nrel,drel;
 GRAPH *ac=( GRAPH *) malloc ((n+2*act_n+100)*sizeof(GRAPH));
 GRAPH *ace=( GRAPH *) malloc ((n+2*act_n+100)*sizeof(GRAPH));

 ac0=(char*) malloc(n+2*act_n+100);
 str=(char*) malloc(2*n+4*act_n+1);
 *flag=0;
 nodestr(n,matr,act_n,act,&p,str);
 transfor(n,matr,act_n,act,&num,ac);
 if (flag_mat) {
 i=0;
 for (k=num;k<num+p-2;k++)
    {
     ac[k].ess=0;
     ac[k].v1=str[i];
     ac[k].v2=str[i+1];
     i++;
    }
 numa=num+p-1;} else numa=num;
 for (i=0;i<p;i++)
    {
     if (str[i] != '0') 
{
     k=0;
     for (j=0;j<num;j++)
	{
	 if (ac[j].v1 == str[i] || ac[j].v2 == str[i]) k++;
	}
     if (k > h) {h=k; s5=str[i];}  } else {s5=str[i]; break;}
    }
 if (!flag_nul && s5 != '0') goto ret;
 for (i=0;i<p;i++)
    {
     if (str[i]==s5) {for (j=i;j<p-1;j++) str[j]=str[j+1]; break;}
    }
 for (k=0;k<numa;k++) ac0[k]=ac[k].v1;
 for (l=0;l<p-1;l++)
    {
     s1=str[l];
     for (i=l+1;i<p-1;i++)
	{
	 s2=str[i];
	 for (j=i+1;j<p-1;j++)
	    {
	     s3=str[j];
	 for (m=j+1;m<p-1;m++)
	    {
	     s4=str[m];
/*	     for (k=0;k<act_n;k++)
		{
		 if (act[k].kol== 999) {
		 o1=act[k].v3;
		 o2=act[k].v4;
		 if (   s1==o1 && s2==o2 || s1==o2 && s2==o1
		     || s1==o1 && s3==o2 || s1==o2 && s3==o1
		     || s1==o1 && s4==o2 || s1==o2 && s4==o1
		     || s1==o1 && s5==o2 || s1==o2 && s5==o1
		     || s2==o1 && s3==o2 || s2==o2 && s3==o1
		     || s2==o1 && s4==o2 || s2==o2 && s4==o1
		     || s2==o1 && s5==o2 || s2==o2 && s5==o1
		     || s3==o1 && s4==o2 || s3==o2 && s4==o1
		     || s3==o1 && s5==o2 || s3==o2 && s5==o1
	     || s4==o1 && s5==o2 || s4==o2 && s5==o1) goto cnt;
		}  */
	     k1=0;
	     bond5(numa,ac,s1,s1,s2,s3,s4,s5,&k1);
             for (h=0;h<numa;h++) ac[h].v1=ac0[h];
             k2=0;
             bond5(numa,ac,s2,s2,s1,s3,s4,s5,&k2);
             for (h=0;h<num;h++) ac[h].v1=ac0[h];
             k3=0;
             bond5(numa,ac,s3,s3,s1,s2,s4,s5,&k3);
             for (h=0;h<numa;h++) ac[h].v1=ac0[h];
             k4=0;
             bond5(numa,ac,s4,s4,s1,s2,s3,s5,&k4);
             for (h=0;h<numa;h++) ac[h].v1=ac0[h];
             k5=0;
             bond5(numa,ac,s5,s5,s1,s2,s3,s4,&k5);
             r1=abs(num-2*k1);
             r2=abs(num-2*k2);
             r3=abs(num-2*k3);
             r4=abs(num-2*k4);
             r5=abs(num-2*k5);
             k=0;
             for (h=0;h<numa;h++) ac[h].v1=ac0[h];
           if (r1<=r2 && r1<=r3 && r1<=r4 && r1<=r5)
               bond5(numa,ac,s1,s1,s2,s3,s4,s5,&k);
             else
             if (r2<=r1 && r2<=r3 && r2<=r4 && r2<=r5)
               bond5(numa,ac,s2,s2,s1,s3,s4,s5,&k);
             else
             if (r3<=r1 && r3<=r2 && r3<=r4 && r3<=r5)
               bond5(numa,ac,s3,s3,s1,s2,s4,s5,&k);
             else
             if (r4<=r1 && r4<=r2 && r4<=r3 && r4<=r5)
               bond5(numa,ac,s4,s4,s1,s2,s3,s5,&k);
             else
             if (r5<=r1 && r5<=r2 && r5<=r3 && r5<=r4) 
               bond5(numa,ac,s5,s5,s1,s2,s3,s4,&k);
             if (!indep5(num,ac,ac0,s1,s2,s3,s4,s5,&k)) goto rep;
	     nrel=num-k; drel=num; rel=nrel/drel;
	     if (rel < 0.5-range || rel > 0.5+range) goto rep;
             fl11=0;fl12=0;fl13=0;fl14=0;fl15=0;
             fl21=0;fl22=0;fl23=0;fl24=0;fl25=0;
             for (h=0;h<num;h++)
                {
	         if (ac[h].v1 =='я')
	           {
	            if (ac0[h]==s1 || ac[h].v2==s1) fl11=1;
	            if (ac0[h]==s2 || ac[h].v2==s2) fl12=1;
	            if (ac0[h]==s3 || ac[h].v2==s3) fl13=1;
	            if (ac0[h]==s4 || ac[h].v2==s4) fl14=1;
	            if (ac0[h]==s5 || ac[h].v2==s5) fl15=1;
	           }
	         else
	           {
	            if (ac0[h]==s1 || ac[h].v2==s1) fl21=1;
	            if (ac0[h]==s2 || ac[h].v2==s2) fl22=1;
	            if (ac0[h]==s3 || ac[h].v2==s3) fl23=1;
	            if (ac0[h]==s4 || ac[h].v2==s4) fl24=1;
	            if (ac0[h]==s5 || ac[h].v2==s5) fl25=1;
	           }
                 if (   fl11 && fl12 && fl13 && fl14 && fl15
                     && fl21 && fl22 && fl23 && fl24 && fl25) goto net;
                }
             goto rep;
         net:
             if (pr1 > fabs(rel-0.5))
               {
                s1o=s1; s2o=s2; s3o=s3; s4o=s4;
		pr1=fabs(rel-0.5);
                for (h=0;h<numa;h++) ace[h]=ac[h];
                *flag=1;    goto fl;
               }
	 rep:
	     for (h=0;h<numa;h++) ac[h].v1=ac0[h];
         cnt: ;
	    }
	    }
        }
    }
fl:
 if (*flag)
   {
    printf(" ^ ");
    ext[0]=s1o; ext[1]=s2o; ext[2]=s3o; ext[3]=s4o; ext[4]='0';//s5; 
    next=4;
    for (k=0;k<num;k++) ac0[k]=ace[k].v1;
	form(num,ac0,next,ext,n,matr,act_n,act);
   }
ret:
 free(str); free(ac);
 free(ace); free(ac0);
}

