
#include "cirsym.h"

extern FILE *out,*inpa;

extern char *b,*c;
extern unsigned long leng;

/***************************************************************/
void fortra(int num,char *ac,SOURCE *act,
            SOURCE *act1,int key)
{
 int i,k=0,act1_n=0;
 for (i=0;i<num;i+=2)
	{
    	if (!key && ac[i] == 'п' || key && ac[i] != 'п')
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
void forte(int num,char *ac,int *act1_n,int key)
{
 int i;

 *act1_n=0;
 for (i=0;i<num;i+=2)
    if (!key && ac[i] == 'п' || key && ac[i] != 'п') (*act1_n)++;
}
/*************************************************************/
void transf(int act_n,SOURCE *act,int *num,GRAPH *ac)
{
 int i;

 *num=-1;
 for (i=0;i<act_n;i++)
    {
     ++(*num);
     ac[*num].v1=act[i].v1;
     ac[*num].v2=act[i].v2;
     ac[*num].ess=i+1;
     ++(*num);
     ac[*num].v1=act[i].v3;
     ac[*num].v2=act[i].v4;
     ac[*num].ess=-(i+1);
    }
 (*num)++;
}
/***************************************************************/
void formm(int num,char *ac,int next,char *ext,
           int act_n,SOURCE *act)
{
 int i,k,l,m,s,act1_n,n2,act2_n,na[2],na0,na1;
 unsigned long dl,dl1,dl2;
 char nus[16][5],sn[16],on[16],o,o1,o2;
 SOURCE *act1;
// 01
// 01
//    00  - 0      Десятичное представление ДВ размерности n=2 (2n=4)
//    01  - 1
//    10  - 2
//    11  - 3
//    011223      Десятичный вектор первой половины ДВ
//    012123      Десятичный вектор второй половины ДВ
//    00 11 12 21 22 33   Поэлементное объединение десятичных векторов подсхемы I
//    33 22 21 12 11 00   Дополняющее отображение по формуле: u-1-i (u-1-j) для подсхемы II
 o=ext[next];
 binvec(next,sn,on,nus,&s);
 forte(num,ac,&act1_n,0);
 forte(num,ac,&act2_n,1);
 dl2=leng+strlen(c);
 na0=act1_n; na1=act2_n;
 for (k=0;k<=s;k++)
      for (l=0;l<=s;l++)
		if(on[k]==on[l]) 
		{
		na[0]=na0+on[k]; na[1]=na1+next-on[k];		
        m=0;
rep:
	act1=( SOURCE *) malloc (na[m]*sizeof(SOURCE));
    fortra(num,ac,act,act1,m);	
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
	    if (act1_n > 1) strcat(b++,"(");
	    dl1=leng+strlen(c);
 	    gggfm(na[m],act1);
	    free(act1);
	    if (dl1>=leng+strlen(c)) {ster(dl); continue;}
	    else if (act1_n > 1) strcat(b++,")");
        if (   *(b-4) == '*' && *(b-3) == '(' && *(b-2) == '1'
            && *(b-1) == ')' )
          {b-=4; b[0]='\0';}			
		if (!m) {m++; goto rep;} 
		kontrol();
	   }
 freact(act_n,act);
}
/*  ********* ЏЋ„ЏђЋѓђЂЊЊЂ “—…’Ђ ЂђЌ€ђЋ‚ **********
    ************************************************ */
void bisec1m(int act_n,SOURCE *act,int *flag)
{
 int i,k,p,num,next=0;
 char *str,*ac0,ext[1];
 GRAPH *ac=( GRAPH *) malloc ((2*act_n)*sizeof(GRAPH));

 str=(char*) malloc(4*act_n);
 ac0=(char*) malloc(2*act_n);
 nodes(act_n,act,&p,str);
 transf(act_n,act,&num,ac);
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
     formm(num,ac0,next,ext,act_n,act);
     break;
rep:
     for (k=0;k<num;k++) ac[k].v1=ac0[k];
    }
 free(ac0); free(ac); free(str);
}
/**************************************************************/
 void bisec2m(int act_n,SOURCE *act,int *flag,float range)
{
 int i,j,i3,j3,h=0,m,k,p,num,numa,p1,p2,next=1;
 int fl11,fl12,fl21,fl22;
 char s1,s2,s1o,s2o,*str,*ac0,ext[2];
 float pr1,rel,nrel,drel;

 GRAPH *ac=( GRAPH *) malloc ((2*act_n+100)*sizeof(GRAPH));
 GRAPH *ace=( GRAPH *) malloc ((2*act_n+100)*sizeof(GRAPH));
 ac0=(char*) malloc(2*act_n+100);
 str=(char*) malloc(2*act_n+1);
 nodes(act_n,act,&p,str);
 pr1 = 1.;
 transf(act_n,act,&num,ac);
     for (i=0;i<p;i++)
       {
	k=0;
	for (j=0;j<num;j++)
	   {
	    if (ac[j].v1 == str[i] || ac[j].v2 == str[i]) k++;
	   }
	if (k > h) {h=k; s1=str[i];}
       }
 i=0;
 for (k=num;k<num+p-1;k++)
    {
     ac[k].ess=0;
     ac[k].v1=str[i];
     ac[k].v2=str[i+1];
     i++;
    }
 numa=num+p-1;
 for (k=0;k<numa;k++) ac0[k]=ac[k].v1;
 *flag=0;
     for (j=0;j<p;j++)
	{
	 if (str[j]==s1) continue;
	 s2=str[j];
	 k=0;
	 bond2(numa,ac,s1,s1,s2,&k);
	 if (   k==numa || k==1
	     || !indep2(numa,ac,ac0,s1,s2,&k)) goto rep;
	 nrel=num-k; drel=num; rel=nrel/drel;
	 if (k==numa-1) goto rep;
	 if (rel <= 0.5-range || rel >= 0.5+range) goto rep;
	 fl11=0;fl12=0;fl21=0;fl22=0;
	 for (m=0;m<num;m++)
	    {
	     if (ac[m].v1 =='п')
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
	    for (k=0;k<numa;k++) ace[k]=ac[k];
	    s1o=s1; s2o=s2;
	    pr1=fabs(rel-0.5);
	   }
    rep:
	 for (k=0;k<numa;k++) ac[k].v1=ac0[k];
	}
con:
 free(str);
 if (*flag)
   {
//    printf(" . ");
	for (k=0;k<num;k++) ac0[k]=ace[k].v1;
	   if (s2o=='0') {ext[0]=s1o; ext[1]=s2o;}
       else {ext[0]=s2o; ext[1]=s1o;}
	   formm(num,ac0,next,ext,act_n,act);
   }
 free(ac); free(ace); free(ac0);
}
/************************************************************/
 void bisec3m(int act_n,SOURCE *act,int *flag,float range)
{
 int i,j,h=0,k,k1,k2,k3,p,num,numa,next,r1,r2,r3,
     fl11,fl12,fl13,fl21,fl22,fl23;
 char s1,s2,s3,s1o,s2o,s3o,*str,*ac0,ext[3];
 float pr1=1.,rel,nrel,nrel1,nrel2,nrel3,drel;
 GRAPH *ac=( GRAPH *) malloc ((2*act_n+100)*sizeof(GRAPH));
 GRAPH *ace=( GRAPH *) malloc ((2*act_n+100)*sizeof(GRAPH));
// int kopt;
 ac0=(char*) malloc(2*act_n+100);
 str=(char*) malloc(2*act_n);
 *flag=0;
 nodes(act_n,act,&p,str);
 transf(act_n,act,&num,ac);
    for (i=0;i<p;i++)
       {
	k=0;
	for (j=0;j<num;j++)
	   {
	    if (ac[j].v1 == str[i] || ac[j].v2 == str[i]) k++;
	   }
	if (k > h) {h=k; s1=str[i];}
       }
 for (i=0;i<p;i++)
    {
     if (str[i]==s1) {for (j=i;j<p-1;j++) str[j]=str[j+1]; break;}
    }
 i=0;
 for (k=num;k<num+p-2;k++)
    {
     ac[k].ess=0;
     ac[k].v1=str[i];
     ac[k].v2=str[i+1];
     i++;
    }
 numa=num+p-2;
 for (k=0;k<numa;k++) ac0[k]=ac[k].v1;
     for (i=0;i<p-1;i++)
        {
         s2=str[i];
	 for (j=i+1;j<p-1;j++)
	    {
             s3=str[j];
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
             nrel=num-k;
	         drel=num;
             rel=nrel/drel;
             if (rel < 0.5-range || rel > 0.5+range) goto rep;
             fl11=0;fl12=0;fl13=0;fl21=0;fl22=0;fl23=0;
             for (h=0;h<num;h++)
                {
	         if (ac[h].v1 =='п')
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
                  s1o=s1; s2o=s2; s3o=s3;
		          pr1=fabs(rel-0.5);
                 // kopt=k;
                  for (h=0;h<numa;h++) ace[h]=ac[h];
                  *flag=1;
                 }
	  rep:
	      for (h=0;h<numa;h++) ac[h].v1=ac0[h];
	    }
	}
 free(str); free(ac); 
 if (*flag)
   {
    ext[0]=s1o; ext[1]=s2o; ext[2]=s3o; next=2;
 //   printf(" ' ");
    for (k=0;k<num;k++) ac0[k]=ace[k].v1;
    formm(num,ac0,next,ext,act_n,act);
   }
 else *flag=0;
 free(ace); free(ac0);
}
/************************************************************/
void bisec4m(int act_n,SOURCE *act,int *flag,float range)
{
 int i,j,l,h=0,k,k1,k2,k3,k4,p,num,numa,next,r1,r2,r3,r4;
 int fl11,fl12,fl13,fl14,fl21,fl22,fl23,fl24;
 char s1,s2,s3,s4,s1o,s2o,s3o,*str,*ac0,ext[4];
 float pr1=1.,rel,nrel,drel;
 GRAPH *ac=( GRAPH *) malloc ((2*act_n+100)*sizeof(GRAPH));
 GRAPH *ace=( GRAPH *) malloc ((2*act_n+100)*sizeof(GRAPH));

 ac0=(char*) malloc(2*act_n+100);
 str=(char*) malloc(2*act_n);
 *flag=0;
 nodes(act_n,act,&p,str);
 transf(act_n,act,&num,ac);
 for (i=0;i<p;i++)
    {
     k=0;
     for (j=0;j<num;j++)
	{
	 if (ac[j].v1 == str[i] || ac[j].v2 == str[i]) k++;
	}
     if (k > h) {h=k; s4=str[i];}
    }
 for (i=0;i<p;i++)
    {
     if (str[i]==s4) {for (j=i;j<p-1;j++) str[j]=str[j+1]; break;}
    }
 i=0;
 for (k=num;k<num+p-2;k++)
    {
     ac[k].ess=0;
     ac[k].v1=str[i];
     ac[k].v2=str[i+1];
     i++;
    }
 numa=num+p-2;
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
             k4=0;
             bond4(numa,ac,s4,s4,s1,s2,s3,&k4);
             r1=abs(num-2*k1);
             r2=abs(num-2*k2);
             r3=abs(num-2*k3);
             r4=abs(num-2*k4);
             k=0;
             for (h=0;h<numa;h++) ac[h].v1=ac0[h];
             if (r1<=r2 && r1<=r3 && r1<=r4) bond4(numa,ac,s1,s1,s2,s3,s4,&k);
             else
             if (r2<=r1 && r2<=r3 && r2<=r4) bond4(numa,ac,s2,s2,s1,s3,s4,&k);
             else
             if (r3<=r1 && r3<=r2 && r3<=r4) bond4(numa,ac,s3,s3,s1,s2,s4,&k);
             else
             if (r4<=r1 && r4<=r2 && r4<=r3) bond4(numa,ac,s4,s4,s1,s2,s3,&k);
             if (!indep4(num,ac,ac0,s1,s2,s3,s4,&k)) goto rep;
             nrel=num-k;
             drel=num;
             rel=nrel/drel;
             if (rel < 0.5-range || rel > 0.5+range) goto rep;
             fl11=0;fl12=0;fl13=0;fl14=0;fl21=0;fl22=0;fl23=0;fl24=0;
             for (h=0;h<num;h++)
                {
	         if (ac[h].v1 =='п')
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
             goto rep;
         net:
             if (pr1 > fabs(rel-0.5))
               {
                s1o=s1; s2o=s2; s3o=s3;
		pr1=fabs(rel-0.5);
                for (h=0;h<numa;h++) ace[h]=ac[h];
                *flag=1;
               }
	 rep:
	     for (h=0;h<numa;h++) ac[h].v1=ac0[h];
	    }
	}
    }
 free(str);  free(ac); 
 if (*flag)
   {
    printf(" : ");
	ext[0]=s1o; ext[1]=s2o; ext[2]=s3o; ext[3]=s4; next=3;
	for (k=0;k<num;k++) ac0[k]=ace[k].v1;
    formm(num,ac0,next,ext,act_n,act);
   }
 else *flag=0;
 free(ace); free(ac0);
}
/***********************************************************************/
void bisec5m(int act_n,SOURCE *act,int *flag,float range)
{
 int i,j,l,h=0,k,k1,k2,k3,k4,k5,m,p,num,numa,next,r1,r2,r3,r4,r5;
 int fl11,fl12,fl13,fl14,fl15,fl21,fl22,fl23,fl24,fl25;
 char s1,s2,s3,s4,s5,s1o,s2o,s3o,s4o,*str,*ac0,ext[5];
 float pr1=1.,rel,nrel,drel;
 GRAPH *ac=( GRAPH *) malloc ((2*act_n+100)*sizeof(GRAPH));
 GRAPH *ace=( GRAPH *) malloc ((2*act_n+100)*sizeof(GRAPH));

 ac0=(char*) malloc(2*act_n+100);
 str=(char*) malloc(2*act_n);
 *flag=0;
 nodes(act_n,act,&p,str);
 transf(act_n,act,&num,ac);
 for (i=0;i<p;i++)
    {
     k=0;
     for (j=0;j<num;j++)
	{
	 if (ac[j].v1 == str[i] || ac[j].v2 == str[i]) k++;
	}
     if (k > h) {h=k; s5=str[i];}
    }
 for (i=0;i<p;i++)
    {
     if (str[i]==s5) {for (j=i;j<p-1;j++) str[j]=str[j+1]; break;}
    }
 i=0;
 for (k=num;k<num+p-2;k++)
    {
     ac[k].ess=0;
     ac[k].v1=str[i];
     ac[k].v2=str[i+1];
     i++;
    }
 numa=num+p-2;
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
	     k1=0;
	     bond5(numa,ac,s1,s1,s2,s3,s4,s5,&k1);
             for (h=0;h<numa;h++) ac[h].v1=ac0[h];
             k2=0;
             bond5(numa,ac,s2,s2,s1,s3,s4,s5,&k2);
             for (h=0;h<numa;h++) ac[h].v1=ac0[h];
             k3=0;
             bond5(numa,ac,s3,s3,s1,s2,s4,s5,&k3);
             k4=0;
             bond5(numa,ac,s4,s4,s1,s2,s3,s5,&k4);
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
	     nrel=num-k;
	     drel=num;
	     rel=nrel/drel;
	     if (rel < 0.5-range || rel > 0.5+range) goto rep;
             fl11=0;fl12=0;fl13=0;fl14=0;fl15=0;
             fl21=0;fl22=0;fl23=0;fl24=0;fl25=0;
             for (h=0;h<num;h++)
                {
	         if (ac[h].v1 =='п')
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
                *flag=1;
               }
	 rep:
	     for (h=0;h<numa;h++) ac[h].v1=ac0[h];
	    }
	    }
        }
    }
con:
free(str); free(ac); 
 if (*flag)
   {
    printf(" ^ ");
    ext[0]=s1o; ext[1]=s2o; ext[2]=s3o; ext[3]=s4o; ext[4]=s5; 
	next=4;
    for (k=0;k<num;k++) ac0[k]=ace[k].v1;
	formm(num,ac0,next,ext,act_n,act);
   }
 else *flag=0;
 free(ace); free(ac0);
}