              REAL*8 FUNCTION TCRS(ZZ,n1,l1,m1)
              IMPLICIT NONE
              REAL*8 Z,ZZ,crs,DGAUSS,xmin,xmax,eps,alpha,pi
              REAL*8 mem,rb,cf
              INTEGER n,l,m,n1,l1,m1
              EXTERNAL FUN
              PARAMETER(xmin=0.0D0,xmax=1.0D0,eps=1.0D-8)
              PARAMETER(alpha=0.00729735D0,pi=3.14159265D0)
c             electron-muon mass ratio and Bohr's radius (in cm)
              PARAMETER(mem=0.00483633D0,rb=0.529 177 211D-8)

c             electron-charged pion mass ratio and Bohr's radius (in cm)
c             PARAMETER(mem=0.00366123298D0,rb=0.529 177 211D-8)

              COMMON /Z/Z
              COMMON /nl/n,l,m
              
              Z=ZZ
              n=n1
              l=l1
              m=m1
 
c             conversion factor
              cf=4*rb**2*mem**2 
              crs=cf*DGAUSS(FUN,xmin,xmax,eps)*(alpha/pi)

              TCRS=crs

              END

              REAL*8 FUNCTION TCRSM(Z,n,l)
              IMPLICIT NONE
              REAL*8 Z,Z1,crs,DGAUSS,xmin,xmax,eps,alpha,pi
              REAL*8 mem,rb,cf
              INTEGER n,l,n1,l1,m
              EXTERNAL FUNM
              PARAMETER(xmin=0.0D0,xmax=1.0D0,eps=1.0D-8)
              PARAMETER(alpha=0.00729735D0,pi=3.14159265D0)
c             electron-muon mass ratio and Bohr's radius (in cm)
              PARAMETER(mem=0.00483633D0,rb=0.529 177 211D-8)
              COMMON /nl/n1,l1,m
              COMMON /Z/Z1

              Z1=Z
              n1=n
              l1=l 
                  
c             conversion factor
              cf=4*rb**2*mem**2 
              crs=cf*DGAUSS(FUNM,xmin,xmax,eps)*(alpha/pi)
             
              TCRSM=crs

              END

              REAL*8 FUNCTION FUN(x)
              IMPLICIT NONE
              REAL*8 x,eps,G
              PARAMETER(eps=1.0D-8)
                        
              IF(x.GT.eps) THEN   
               FUN=G(x)+G(1.0D0/x)/x**2
              ELSE
               FUN=G(x)
              END IF

              END 

              REAL*8 FUNCTION FUNM(x)
              IMPLICIT NONE
              REAL*8 x,eps,GM
              PARAMETER(eps=1.0D-8)
                        
              IF(x.GT.eps) THEN   
               FUNM=GM(x)+GM(1.0D0/x)/x**2
              ELSE
               FUNM=GM(x)
              END IF
 
              END 


              REAL*8 FUNCTION G(x)
              IMPLICIT NONE
              REAL*8 x,U,F             

              G=U(x)**2*(1-F(x))*x
 
              END 

              REAL*8 FUNCTION GM(x)
              IMPLICIT NONE
              REAL*8 x,U,FM             

              GM=U(x)**2*(1-FM(x))*x
               
              END 

              REAL*8 FUNCTION F(x)
              IMPLICIT NONE
              REAL*8 x,FTRANS
              INTEGER n,l,m,iq(2)
              COMMON /nl/n,l,m
              COMMON /iq/iq
                   
              F=(-1)**(l+m)*FTRANS(iq,n,l,m,n,l,m,x)
                                             
              END

              REAL*8 FUNCTION LFACT(n)
              IMPLICIT NONE
              INTEGER i,n

              LFACT=0

              IF(n.GT.0) THEN
              DO i=1,n
                 LFACT=LFACT+DLOG(DBLE(i))
              END DO
              END IF
                 
              END

              REAL*8 FUNCTION Fhyp(n,mm,k,ph,alp)
              IMPLICIT NONE
              REAL*8 ph,JACOB,alp,bt
              INTEGER n,k,m,mm,nn 

                m=mm/2 

              IF(MOD(mm,2).EQ.0) THEN
                nn=n+m-k
                bt=0.5D0
          Fhyp=DCOS(ph)**(2*m)*JACOB(nn,alp,bt,DCOS(2*ph))
                Fhyp=Fhyp/JACOB(nn,alp,bt,1.0D0)
                
              ELSE
               nn=n+m+1-k 
               bt=-0.5D0
          Fhyp=DCOS(ph)**(2*m)*JACOB(nn,alp,bt,DCOS(2*ph))
                Fhyp=Fhyp/JACOB(nn,alp,bt,1.0D0)

              END IF

              END
 
              REAL*8 FUNCTION FM(x)
              IMPLICIT NONE
              REAL*8 x,ph,JACOB
              INTEGER n,l,m
              COMMON /nl/n,l,m

              IF(n.EQ.1.AND.l.EQ.0) THEN 
                FM=16.0D0/(4+x**2)**2

              ELSE IF(n.EQ.2.AND.l.EQ.1) THEN
                FM=(1-x**2)/(1+x**2)**4

              ELSE IF(n.EQ.2.AND.l.EQ.0) THEN
                FM=(1-3*x**2+2*x**4)/(1+x**2)**4

              ELSE IF(n.EQ.3.AND.l.EQ.0) THEN
           FM=(27*x**2-4)*(3*x**2-4)*(243*x**4-216*x**2+16)
                FM=16*FM/(4+9*x**2)**6

              ELSE IF(n.EQ.4.AND.l.EQ.0) THEN
        FM=(4*x**2-1)*(16*x**4-24*x**2+1)*(256*x**6-288*x**4+48*x**2-1)
                 FM=FM/(1+4*x**2)**8

              ELSE IF(n.EQ.3.AND.l.EQ.1) THEN
                FM=(27*x**2-4)*(3*x**2-4)*(1-9*x**2)
                FM=256*FM/(4+9*x**2)**6

              ELSE IF(n.EQ.3.AND.l.EQ.2) THEN
                 FM=(27*x**2-4)*(3*x**2-4)
                 FM=256*FM/(4+9*x**2)**6

              ELSE IF(n.EQ.4.AND.l.EQ.1) THEN
       FM=(1-4*x**2)*(16*x**4-24*x**2+1)*(160*x**4-40*x**2+1)
              FM=FM/(1+4*x**2)**8

              ELSE IF(n.EQ.4.AND.l.EQ.2) THEN
                 FM=(4*x**2-1)*(16*x**4-24*x**2+1)*(24*x**2-1)
                 FM=FM/(1+4*x**2)**8

               ELSE IF(n.EQ.4.AND.l.EQ.3) THEN
                 FM=(1-4*x**2)*(16*x**4-24*x**2+1)
                 FM=FM/(1+4*x**2)**8
             
               ELSE
                ph=DATAN(n*x/2.0D0)
                FM=JACOB(n-l-1,0.0D0,DBLE(2*l+1),DCOS(2*ph)) 
                FM=FM*DSIN(2*n*ph)*(DCOS(ph))**(2*l+4)
                FM=FM/(n*DSIN(2*ph))

              END IF
                           
              END

              REAL*8 FUNCTION JACOB(n,alp,bt,x)
              IMPLICIT NONE
              integer i,n
              REAL*8 x,alp,bt,P(0:n),a,b,c
               
              P(0)=1.0D0
              IF(n.EQ.0) GO TO 1
              P(1)=0.5*(alp-bt+(2+alp+bt)*x)

              IF(n.GE.2) THEN
              DO i=2,n
                a=2*i*(i+alp+bt)*(2*i+alp+bt-2)
                b=(2*i+alp+bt-2)*(2*i+alp+bt)*x+alp**2-bt**2
                b=(2*i+alp+bt-1)*b
                c=2*(i-1+alp)*(i-1+bt)*(2*i+alp+bt)
                P(i)=(b/a)*P(i-1)-(c/a)*P(i-2)
              END DO
              END IF

1             JACOB=P(n)
              
              END

              REAL*8 FUNCTION U(x)
              IMPLICIT NONE
              REAL*8 x,Z,a1,a2,a3,b1,b2,b3,pi,alpha,bt1,bt2,bt3,me,mem
              PARAMETER(a1=0.10D0,a2=0.55D0,a3=0.35D0)
              PARAMETER(b1=6.0D0,b2=1.2D0,b3=0.3D0)              
              PARAMETER(alpha=0.00729735D0,pi=3.14159265D0)
c             electron-muon mass ratio
              PARAMETER(mem=0.00483633D0)

c             electron-charged pion mass ratio              
c             PARAMETER(mem=0.00366123298D0)       

              COMMON /Z/Z
     
              me=(2/alpha)*mem
              bt1=me*(b1/121.0D0)*Z**(1.0D0/3.0D0)
              bt2=me*(b2/121.0D0)*Z**(1.0D0/3.0D0)              
              bt3=me*(b3/121.0D0)*Z**(1.0D0/3.0D0)
            
              U=a1/(x**2+bt1**2)+a2/(x**2+bt2**2)+a3/(x**2+bt3**2)
              U=4*pi*Z*DSQRT(alpha)*U
  
              END

c           i**(l1+l2)*(-1)**m1 phase is not included      
              REAL*8 FUNCTION FTRANS(iq,n1,l1,m1,n2,l2,m2,x)
              IMPLICIT NONE
              INTEGER n1,l1,m1,n2,l2,m2,l,m,s,sm,iq(2)
              REAL*8 N,Al,Il,a,b,sg,LFACT,x,DWIG3J,Blm,IlJN,IlW

              m=m1-m2

              a=DBLE(n2)/DBLE(n1+n2)
              b=DBLE(n1)/DBLE(n1+n2)
              sg=x*DBLE(n1*n2)/DBLE(n1+n2)

              N=LFACT(n1-l1-1)+LFACT(n2-l2-1)
              N=N-LFACT(n1+l1)-LFACT(n2+l2)
              N=DEXP(0.5D0*N)*DSQRT(DBLE((2*l1+1)*(2*l2+1)))
              N=(2*a)**(l1+1)*(2*b)**(l2+1)*N
              N=N/DBLE(n1+n2)
    
              sm=MIN0(l1,l2)
              FTRANS=0.0D0
              DO s=0,sm
              l=l1+l2-2*s
              Al=(-1)**s*(2*l+1)
              Al=Al*DWIG3J(DBLE(l1),DBLE(l2),DBLE(l),0.0D0,0.0D0,0.0D0)
              Al=Al*DWIG3J(DBLE(l1),DBLE(l2),DBLE(l),DBLE(m1),-DBLE(m2)
     &                     ,-DBLE(m))
              Al=Al*Blm(l,m,iq(1))
              IF(iq(2).EQ.1) THEN
              FTRANS=FTRANS+Al*Il(n1,l1,n2,l2,l,a,b,sg)
              ELSE IF(iq(2).EQ.2) THEN
              FTRANS=FTRANS+Al*IlJN(n1,l1,n2,l2,s,a,b,sg) 
              ELSE IF(iq(2).EQ.3) THEN
              FTRANS=FTRANS+Al*IlW(n1,l1,n2,l2,s,a,b,sg)
              END IF
              END DO
              FTRANS=N*FTRANS              

              END

c            Acta Phys. Pol. B 52, 1209 (2021)
             REAL*8 FUNCTION Il(n1,l1,n2,l2,l,a,b,sg)
             IMPLICIT NONE
             INTEGER n1,l1,n2,l2,l,m1,m2
             REAL*8 a,b,sg,N,Nm,Ir,LFACT,ph,alp,Fhyp,cs

             Il=0.0D0
             ph=DATAN(sg)
             alp=l+0.5D0
             cs=DCOS(ph)**(2*(l+2))

             DO m1=0,n1-l1-1
              DO m2=0,n2-l2-1
               Nm=LFACT(n1+l1)+LFACT(n2+l2)+LFACT(l+l1+l2+m1+m2+2) 
               Nm=Nm-LFACT(m1)-LFACT(m2)-LFACT(n1-l1-1-m1)
               Nm=Nm-LFACT(n2-l2-1-m2)-LFACT(2*l1+1+m1)
               Nm=Nm-LFACT(2*l2+1+m2)
               Nm=Dexp(Nm)               
               Ir=(-1)**(m1+m2)*(2*a)**m1*(2*b)**m2*Nm               
               Ir=Ir*Fhyp(l,m1+m2+l1+l2-l,l,ph,alp)
               Il=Il+Ir*cs*sg**l
              END DO
             END DO

             N=LFACT(l)-LFACT(2*l+1)
             N=2.0D0**l*DEXP(N)
             
             Il=N*Il
     
             END

c            EXP(i*m*phi) is not included
             REAL*8 FUNCTION Blm(l,m,iq)
             IMPLICIT NONE
             INTEGER l,m,iq,kp,km
             REAL*8 LFACT

             IF(iq.EQ.0) THEN
              IF(m.EQ.0) THEN
               Blm=1.0D0
              ELSE
               Blm=0.0D0
              END IF

             ELSE IF(iq.EQ.1) THEN
              IF(MOD(l+m,2).EQ.1) THEN
               Blm=0.0D0
              ELSE
               kp=(l+m)/2
               km=(l-m)/2
               Blm=0.5D0*(LFACT(l+m)+LFACT(l-m))
               Blm=Blm-LFACT(kp)-LFACT(km)
               Blm=DEXP(Blm)
               Blm=(-1)**(m+kp)*Blm/2.0D0**l
              END IF

             ELSE

             WRITE(*,1) 'Undefined quantization axis!'
             STOP

             END IF  

 1           FORMAT(A)

             END

c            Phys. Rep.511, 1 (2012)
             REAL*8 FUNCTION IlW(n1,l1,n2,l2,s,a,b,sg)
             IMPLICIT NONE
             INTEGER n1,l1,n2,l2,s,k
             REAL*8 a,b,sg,Ck,Jtk

             IlW=0.0D0
            
             DO k=0,n1+n2-l1-l2-2
              IlW=IlW+Ck(n1,l1,n2,l2,k,a,b)*Jtk(s,k,l1+l2,sg)
             END DO

             END

             REAL*8 FUNCTION Ck(n1,l1,n2,l2,k,a,b)
             IMPLICIT NONE
             INTEGER n1,l1,n2,l2,k,j
             REAL*8 a,b,LFACT,Ir

             Ck=0.0D0

             DO j=0,k
              IF(n1-l1-1-j.GE.0.AND.n2-l2-1+j-k.GE.0) THEN
               Ir=LFACT(n1+l1)+LFACT(n2+l2)-LFACT(j)-LFACT(k-j)
               Ir=Ir-LFACT(n1-l1-1-j)-LFACT(n2-l2-1+j-k)
               Ir=Ir-LFACT(2*l1+1+j)-LFACT(2*l2+1+k-j)
               Ir=DEXP(Ir)
              ELSE
               Ir=0.0D0
              END IF  
              Ck=Ck+Ir*(2*a)**j*(2*b)**(k-j)
             END DO

             Ck=(-1)**k*Ck

             END

             REAL*8 FUNCTION Jtk(s,k,l,sg)
             IMPLICIT NONE
             INTEGER s,k,l,p,pm
             REAL*8 sg, Ir,LFACT 

             pm=s+(k+1)/2

             Jtk=0.0D0

             DO p=0,pm
              Ir=LFACT(2*s+k+1)+LFACT(l-p+k+1)
              Ir=Ir+(2*(s-p)+k+1)*DLOG(2.0D0)
              Ir=Ir-LFACT(p)-LFACT(2*(s-p)+k+1)
              Ir=DEXP(Ir)
              Jtk=Jtk+(-1)**p*Ir/(1+sg**2)**(l-p+k+2)
             END DO

              Jtk=Jtk*(2*sg)**(l-2*s)  

             END

c            Phys. Atom. Nucl.59, 2130 (1996)
             REAL*8 FUNCTION IlJN(n1,l1,n2,l2,s,a,b,sg)
             IMPLICIT NONE
             INTEGER n1,l1,n2,l2,s,k,ll,km,p,nr1,nr2,alp,bt
             REAL*8 a,b,sg,LFACT,Bps,Hk,w,z,GGNBG,Ik

             IlJN=0.0D0
             ll=l1+l2
             km=n1+n2-ll-2
             nr1=n1-l1-1
             nr2=n2-l2-1
             alp=2*l1+1
             bt=2*l2+1
             w=1/(1+sg**2)
             z=1-2*w
                
             DO p=0,s
              DO k=0,km
               Bps=LFACT(s)+LFACT(ll-s-p+1)+LFACT(2*(ll-s+1))
               Bps=Bps-LFACT(2*(ll-s-p+1))-LFACT(ll-s+1)
               Bps=Bps-LFACT(s-p)-LFACT(p)
               Bps=(-1)**(s-p)*DEXP(Bps)
               Ik=LFACT(ll-p+1)+(ll-p+2)*DLOG(w)
               Ik=DEXP(Ik)
               Ik=2*(2*sg)**(ll-2*p)*Ik*GGNBG(k,DBLE(ll+2),DBLE(p),z)
               IlJN=IlJN+Bps*Ik*Hk(nr1,alp,nr2,bt,a,b,k)
              END DO
             END DO

          END

              REAL*8 FUNCTION Hk(n,alp,m,bt,a,b,k)
              IMPLICIT NONE
              INTEGER n,alp,m,bt,k
              REAL*8 a,b,Nr,LFACT,JACOBR

              Nr=LFACT(k)+LFACT(n+m-k)-LFACT(n)-LFACT(m)

c            Niukkanen form
              Nr=DEXP(Nr)*a**(k-m)*b**(k-n)
              Hk=Nr*JACOBR(n+m-k,k-m,k-n,b-a)
              Hk=Hk*JACOBR(n+m-k,alp+k-m,bt+k-n,b-a)

c             Afanasyev and Tarasov form
c              Nr=(-1)**(k+m)*DEXP(Nr) 
c              Hk=Nr*JACOBR(n+m-k,alp+k-m,bt+k-n,b-a)
c              Hk=Hk*JACOBR(k,m-k,n-k,b-a)
 
              END

              REAL*8 FUNCTION GGNBG(n,lm,pp,x)
              IMPLICIT NONE
              integer i,n
              REAL*8 x,lm,pp,P(-1:n),a,b,c

              P(-1)=0.0D0
              P(0)=1.0D0
              IF(n.EQ.0) GO TO 1
              P(1)=2*pp+2*(lm-pp)*x
              IF(n.EQ.1) GO TO 1
              P(2)=2*(pp-lm)*(pp-lm-1)*x**2-4*pp*(pp-lm)*x
              P(2)=P(2)+pp-lm+(1+2*pp)*pp  

              IF(n.GE.3) THEN
              DO i=3,n
                a=i-1+2*pp+2*x*(i-1+lm-pp)
                b=i-2+2*lm-2*pp+2*x*(i-2+lm+pp)
                c=i-3+2*lm
                P(i)=(a*P(i-1)-b*P(i-2)+c*P(i-3))/DBLE(i)
              END DO
              END IF

1             GGNBG=P(n)+P(n-1)

              END

              REAL*8 FUNCTION JACOBR(n,alp,bt,x)
              IMPLICIT NONE
              integer n,alp,bt
              REAL*8 x,JACOB,LFACT,Nr
              
              IF(alp.GE.0.AND.bt.GE.0) THEN
               JACOBR=JACOB(n,DBLE(alp),DBLE(bt),x)
              ELSE IF(alp.LT.0.AND.bt.LT.0) THEN
               JACOBR=((x-1)/2.0D0)**(-alp)*((x+1)/2.0D0)**(-bt)
               JACOBR=JACOBR*JACOB(n+alp+bt,-DBLE(alp),-DBLE(bt),x)
              ELSE IF(alp.LT.0.AND.bt.GE.0) THEN
               Nr=LFACT(n+bt)+LFACT(n+alp)-LFACT(n)-LFACT(n+alp+bt)
               Nr=DEXP(Nr)
               JACOBR=Nr*((x-1)/2.0D0)**(-alp)
               JACOBR=JACOBR*JACOB(n+alp,-DBLE(alp),-DBLE(bt),x)
              ELSE IF(alp.GE.0.AND.bt.LT.0) THEN
               Nr=LFACT(n+bt)+LFACT(n+alp)-LFACT(n)-LFACT(n+alp+bt)
               Nr=DEXP(Nr) 
               JACOBR=Nr*((x+1)/2.0D0)**(-bt)
               JACOBR=JACOBR*JACOB(n+bt,DBLE(alp),-DBLE(bt),x)
              END IF

              END


              REAL*8 FUNCTION wTRCRS(iq1,iq2,ZZ,nn1,ll1,mm1,nn2,ll2,mm2)
              IMPLICIT NONE
              REAL*8 ZZ, TRCRS
              INTEGER nn1,ll1,mm1,nn2,ll2,mm2
              INTEGER iq1,iq2,iq(2)
              COMMON /iq/iq
              iq(1)=iq1
              iq(2)=iq2
              wTRCRS=TRCRS(ZZ,nn1,ll1,mm1,nn2,ll2,mm2)
              END

              REAL*8 FUNCTION wTCRS(iq1,iq2,ZZ,nn,ll,mm)
              IMPLICIT NONE
              REAL*8 ZZ,TCRS
              INTEGER nn,ll,mm
              INTEGER iq1,iq2,iq(2)
              COMMON /iq/iq
              iq(1)=iq1
              iq(2)=iq2
              wTCRS=TCRS(ZZ,nn,ll,mm)
              END

c             transition cross section
              REAL*8 FUNCTION TRCRS(ZZ,nn1,ll1,mm1,nn2,ll2,mm2)
              IMPLICIT NONE
              REAL*8 Z,ZZ,crs,DGAUSS,xmin,xmax,eps,alpha,pi
              REAL*8 mem,rb,cf
              INTEGER n1,m1,l1,n2,m2,l2,nn1,ll1,mm1,nn2,ll2,mm2
              EXTERNAL FUNT
              PARAMETER(xmin=0.0D0,xmax=1.0D0,eps=1.0D-8)
              PARAMETER(alpha=0.00729735D0,pi=3.14159265D0)
c             electron-muon mass ratio and Bohr's radius (in cm)
              PARAMETER(mem=0.00483633D0,rb=0.529 177 211D-8)
              COMMON /Z/Z
              COMMON /nlm/n1,l1,m1,n2,l2,m2

              Z=ZZ
              n1=nn1
              l1=ll1 
              m1=mm1
              n2=nn2
              l2=ll2 
              m2=mm2
             
c             conversion factor
              cf=4*rb**2*mem**2

              IF((1-(-1)**(l2-l1)).EQ.0) THEN
                TRCRS=0.0D0
              ELSE
                crs=cf*DGAUSS(FUNT,xmin,xmax,eps)*(alpha/pi)
                crs=(1-(-1)**(l2-l1))*crs
                TRCRS=crs
              END IF

              END

              REAL*8 FUNCTION FUNT(x)
              IMPLICIT NONE
              REAL*8 x,eps,GT
              PARAMETER(eps=1.0D-8)
                        
              IF(x.GT.eps) THEN   
               FUNT=GT(x)+GT(1.0D0/x)/x**2
              ELSE
               FUNT=GT(x)
              END IF

              END 
 
              REAL*8 FUNCTION GT(x)
              IMPLICIT NONE
              REAL*8 x,U,FT

              GT=U(x)**2*FT(0.5D0*x)**2*x


              END
                                                        
              REAL*8 FUNCTION FT(x)
              IMPLICIT NONE
              REAL*8 x,FTRANS
              INTEGER n1,l1,m1,n2,l2,m2,iq(2)  
              COMMON /nlm/n1,l1,m1,n2,l2,m2
              COMMON /iq/iq

               FT=FTRANS(iq,n1,l1,m1,n2,l2,m2,x)
               
              END 



              REAL*8 FUNCTION DGAUSS(F,A,B,EPS)
              IMPLICIT NONE
              REAL*8 F,A,B,EPS,AA,BB,H,CONST,C1,C2,S8,S16,U
              REAL*8 W(12),X(12),Z1,HF,CST
              INTEGER I

              PARAMETER (Z1 = 1.0D0, HF = 0.5D0, CST = 5.0D-3)

              DATA X( 1) /9.6028985649753623D-1/, 
     &        W( 1) /1.0122853629037626D-1/
              DATA X( 2) /7.9666647741362674D-1/, 
     &        W( 2) /2.2238103445337447D-1/
              DATA X( 3) /5.2553240991632899D-1/, 
     &        W( 3) /3.1370664587788729D-1/
              DATA X( 4) /1.8343464249564980D-1/, 
     &        W( 4) /3.6268378337836198D-1/
              DATA X( 5) /9.8940093499164993D-1/, 
     &        W( 5) /2.7152459411754095D-2/
              DATA X( 6) /9.4457502307323258D-1/, 
     &        W( 6) /6.2253523938647893D-2/
              DATA X( 7) /8.6563120238783174D-1/, 
     &        W( 7) /9.5158511682492785D-2/
              DATA X( 8) /7.5540440835500303D-1/, 
     &        W( 8) /1.2462897125553387D-1/
              DATA X( 9) /6.1787624440264375D-1/, 
     &        W( 9) /1.4959598881657673D-1/
              DATA X(10) /4.5801677765722739D-1/, 
     &        W(10) /1.6915651939500254D-1/
              DATA X(11) /2.8160355077925891D-1/, 
     &        W(11) /1.8260341504492359D-1/
              DATA X(12) /9.5012509837637440D-2/, 
     &        W(12) /1.8945061045506850D-1/
              

              H=0
              IF(B .EQ. A) GO TO 99
              CONST=CST/DABS(B-A)
              BB=A
 1            AA=BB
              BB=B
 2            C1=HF*(BB+AA)
              C2=HF*(BB-AA)
              S8=0
              DO 3 I = 1,4
              U=C2*X(I)
 3            S8=S8+W(I)*(F(C1+U)+F(C1-U))
              S16=0
              DO 4 I = 5,12
              U=C2*X(I)
 4            S16=S16+W(I)*(F(C1+U)+F(C1-U))
              S16=C2*S16
              IF(DABS(S16-C2*S8) .LE. EPS*(1+DABS(S16))) THEN
                 H=H+S16
                 IF(BB .NE. B) GO TO 1
              ELSE
                 BB=C1
                 IF(1+CONST*DABS(C2) .NE. 1) GO TO 2
                 H=0
                 WRITE(*,*) 'DGAUSS ERROR: TOO HIGH ACCURACY REQUIRED'
                 GO TO 99
              END IF
 99           DGAUSS=H
              RETURN
              END
              


               
              REAL*8 FUNCTION DWIG3J(A1,B1,C1,X1,Y1,Z1)
              IMPLICIT NONE
              REAL*8  A1,B1,C1,X1,Y1,Z1,R1,HF,F,W,H,S,Q
              REAL*8 U(0:202)
              INTEGER IA,IB,IC,IX,IY,IZ,JX,JY,JZ,J0,KA,KZ,K,N
              INTEGER K0,K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11
              PARAMETER (R1 = 1.0D0, HF = 0.5D0)

              DATA U(0),U(2),(U(2*N-1),N=1,101) /103*0/
              DATA U(  4),U(  6) /6.931471805599453D-01, 
     &        1.791759469228055D+00/
              DATA U(  8),U( 10) /3.178053830347946D+00,
     &        4.787491742782046D+00/
              DATA U( 12),U( 14) /6.579251212010101D+00,
     &        8.525161361065414D+00/
              DATA U( 16),U( 18) /1.060460290274525D+01, 
     &        1.280182748008147D+01/
              DATA U( 20),U( 22) /1.510441257307552D+01, 
     &        1.750230784587389D+01/
              DATA U( 24),U( 26) /1.998721449566189D+01, 
     &        2.255216385312342D+01/
              DATA U( 28),U( 30) /2.519122118273868D+01,
     &        2.789927138384089D+01/
              DATA U( 32),U( 34) /3.067186010608067D+01,
     &        3.350507345013689D+01/
              DATA U( 36),U( 38) /3.639544520803305D+01,
     &        3.933988418719949D+01/
              DATA U( 40),U( 42) /4.233561646075349D+01,
     &        4.538013889847691D+01/
              DATA U( 44),U( 46) /4.847118135183522D+01,
     &        5.160667556776437D+01/
              DATA U( 48),U( 50) /5.478472939811232D+01,
     &        5.800360522298052D+01/
              DATA U( 52),U( 54) /6.126170176100200D+01,
     &        6.455753862700633D+01/
              DATA U( 56),U( 58) /6.788974313718153D+01,
     &        7.125703896716801D+01/
              DATA U( 60),U( 62) /7.465823634883016D+01,
     &        7.809222355331531D+01/
              DATA U( 64),U( 66) /8.155795945611504D+01,
     &        8.505446701758152D+01/
              DATA U( 68),U( 70) /8.858082754219768D+01,
     &        9.213617560368709D+01/
              DATA U( 72),U( 74) /9.571969454214320D+01,
     &        9.933061245478743D+01/
              DATA U( 76),U( 78) /1.029681986145138D+02,
     &        1.066317602606435D+02/
              DATA U( 80),U( 82) /1.103206397147574D+02,
     &        1.140342117814617D+02/
              DATA U( 84),U( 86) /1.177718813997451D+02,
     &        1.215330815154386D+02/
              DATA U( 88),U( 90) /1.253172711493569D+02,
     &        1.291239336391272D+02/
              DATA U( 92),U( 94) /1.329525750356163D+02,
     &        1.368027226373264D+02/
              DATA U( 96),U( 98) /1.406739236482343D+02,
     &        1.445657439463449D+02/
              DATA U(100),U(102) /1.484777669517730D+02,
     &        1.524095925844974D+02/
              DATA U(104),U(106) /1.563608363030788D+02,
     &        1.603311282166309D+02/
              DATA U(108),U(110) /1.643201122631952D+02,
     &        1.683274454484277D+02/
              DATA U(112),U(114) /1.723527971391628D+02,
     &        1.763958484069974D+02/
              DATA U(116),U(118) /1.804562914175438D+02,
     &        1.845338288614495D+02/
              DATA U(120),U(122) /1.886281734236716D+02,
     &        1.927390472878449D+02/
              DATA U(124),U(126) /1.968661816728900D+02,
     &        2.010093163992815D+02/
              DATA U(128),U(130) /2.051681994826412D+02,
     &        2.093425867525368D+02/
              DATA U(132),U(134) /2.135322414945633D+02,
     &        2.177369341139542D+02/
              DATA U(136),U(138) /2.219564418191303D+02,
     &        2.261905483237276D+02/
              DATA U(140),U(142) /2.304390435657770D+02,
     &        2.347017234428183D+02/
              DATA U(144),U(146) /2.389783895618343D+02,
     &        2.432688490029827D+02/
              DATA U(148),U(150) /2.475729140961869D+02,
     &        2.518904022097232D+02/
              DATA U(152),U(154) /2.562211355500095D+02,
     &        2.605649409718632D+02/
              DATA U(156),U(158) /2.649216497985528D+02,
     &        2.692910976510198D+02/
              DATA U(160),U(162) /2.736731242856937D+02,
     &        2.780675734403661D+02/
              DATA U(164),U(166) /2.824742926876304D+02,
     &        2.868931332954270D+02/
              DATA U(168),U(170) /2.913239500942703D+02,
     &        2.957666013507606D+02/
              DATA U(172),U(174) /3.002209486470141D+02,
     &        3.046868567656687D+02/
              DATA U(176),U(178) /3.091641935801469D+02,
     &        3.136528299498791D+02/
              DATA U(180),U(182) /3.181526396202093D+02,
     &        3.226634991267262D+02/
              DATA U(184),U(186) /3.271852877037752D+02,
     &        3.317178871969285D+02/
              DATA U(188),U(190) /3.362611819791985D+02,
     &        3.408150588707990D+02/
              DATA U(192),U(194) /3.453794070622669D+02,
     &        3.499541180407702D+02/
              DATA U(196),U(198) /3.545390855194408D+02,
     &        3.591342053695754D+02/
              DATA U(200),U(202) /3.637393755555635D+02,
     &        3.683544960724047D+02/


              H=0.0D0
              IA=IDNINT(2*A1)
              IB=IDNINT(2*B1)
              IC=IDNINT(2*C1)
              IX=IDNINT(2*X1)
              IY=IDNINT(2*Y1)
              IZ=IDNINT(2*Z1)
              IF(IA .LT. 0 .OR. IB .LT. 0 .OR. IC .LT. 0) GO TO 99
              IF(MOD(IA+IB+IC,2) .NE. 0) GO TO 99
              JX=IABS(IX)
              JY=IABS(IY)
              JZ=IABS(IZ)
              IF(IA .LT. JX .OR. IB .LT. JY .OR. IC .LT. JZ) GO TO 99
              IF(MOD(IA+JX,2) .NE. 0 .OR. MOD(IB+JY,2) .NE. 0) GOTO 99
              IF(MOD(IC+JZ,2) .NE. 0) GO TO 99
              J0=IA-IB-IZ
              F=1.0D0
              IF(IX+IY+IZ .NE. 0 .OR. MOD(J0,2) .NE. 0) GO TO 99
              K0=IA+IB+IC+2
              K1=IA+IB-IC
              K2=IA-IB+IC
              K3=IB+IC-IA
              IF(K1 .LT. 0 .OR. K2 .LT. 0 .OR. K3 .LT. 0) GO TO 99
              K4=IA+IX
              K5=IB+IY
              K6=IC+IZ
              K7=IA-IX
              K8=IB-IY
              K9=IC-IZ
              K10=IB-IC-IX
              K11=IA-IC+IY
              KA=MAX0(0,K10,K11)
              KZ=MIN0(K1,K5,K7)
              W=HF*(U(K1)+U(K2)+U(K3)+U(K4)+U(K5)+U(K6)+
     &        U(K7)+U(K8)+U(K9)-U(K0))
              S=0.0D0
              Q=(-1)**((KA+J0)/2)
              DO 1 K = KA,KZ,2
              S=S+Q*DEXP(W-(U(K)+U(K1-K)+U(K5-K)+U(K7-K)+
     &        U(K-K10)+U(K-K11)))
 1            Q=-Q
              H=F*S
              GO TO 99

 99           DWIG3J=H
              RETURN
              END                 
