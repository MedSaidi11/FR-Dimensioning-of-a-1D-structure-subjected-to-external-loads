% Les constantes
global l=input("Si votre structure est soumise a une flexion tapez 1 sinon si ell est soumise a une compresion ou traction tapez 0 \n")
if (l==1)
  global E=input("entrer E")
  global I=input("entrer I")
  global H=input("entrer H")
  global L=input("entrer L")
  global F0=input("entrer F0")
  global q0=input("entrer q0")
  global nh=input("entrer le nombre de noeuds")
else
  global Ea=input("entrer Ea")
  global F0=input("entrer F0")
  global q0=input("entrer q0")
  global A0=input("entrer A0")
  global p=input("tapez 1 si votre section est variable et 0 si elle est fixe \n")
  if (p==1)
    global a=input("entrer a")
    global b=input("entrer b")
  elseif(p==0) 
    global S=input("entrer la section")
  endif
  global fv=input("entrer la force volumique")
  global nh=input("entrer le nombre de noeuds")
endif
function y=psi(i,x,h)
  if (i==1) 
    y=1-3*(x/h)^2+2*(x/h)^3;
  elseif (i==2)
    y=x-2*(x^2)/h+(x^3/h^2);
  elseif (i==3)
    y=3*(x/h)^2-2*(x/h);
  else 
    y=-x^2/h+x^3/h^2;
  endif
endfunction
function y=dpsi(i,x,h)
  if (i==1) 
    y=-6/h**2 + 12*x/h**3;
  elseif (i==2)
    y=-4/h + 6*x/h**2;
  elseif (i==3)
    y=6/h**2 -12*x/h**3;
  else 
  y=-2/h + 6*x/h**2;
  endif
endfunction
function y=psii(i,x,h)
  if (i==1) 
    y=1-x/h ;
  else
    y=x/h;
  endif
endfunction 
function y=dpsii(i,x,h)
  if (i==1) 
     y=-1/h;
  else
     y=1/h;
  endif
endfunction
function y=EI(x)
  global E;
  global I;
  y=E*I;  
endfunction
function y=EA(x)
 global S;
 global Ea;
 global a;
 global b;
 if (s==1) 
  y=Ea*(a*x+b)
 else
  y=Ea*S
 endif
endfunction
function Ke=Kele(x1,x2,EI,EA)
 global l;
 if (l==1)
    for i=1:4
      for j=1:4
        x=0:0.005:x2-x1;
        for k=1:length(x)
          y(k)=EI(x(k))*dpsi(i,x(k),x2-x1)*dpsi(j,x(k),x2-x1);
        endfor
       Ke(i,j)=trapz(x,y);
      endfor
    endfor
 elseif (l==0)
    for i=1:2
     for j=1:2
       x=0:0.005:x2-x1;
       for k=1:length(x)
         y(k)=EA(x(k))*dpsii(i,x(k),x2-x1)*dpsii(j,x(k),x2-x1);
       endfor
       Ke(i,j)=trapz(x,y);
     endfor
    endfor
 endif 
endfunction        
function y=q(x)
  global q0;
  global L;
  y=(q0/L)*x;
endfunction
function[fe]=fele(x1,x2,q,S)
  global l;
  global fv;
  if (l==1)
    for i=1:4
      x=0:0.005:x2-x1;
      for k=1:length(x)
        y(k)=psi(i,x(k),x2-x1)*q(x(k)+x1);
      endfor
      fe(i)=trapz(x,y);
    endfor
 elseif (l==0)
    for i=1:2
      x=0:0.005:x2-x1;
      for k=1:length(x)
        y(k)=psi(i,x(k),x2-x1)*fv*S(x(k));
      endfor
      fe(i)=trapz(x,y);
    endfor
 endif
endfunction
function[x]=mailleregulier(nh)
  global L;
  for i=1:nh+1
    x(i)=(i-1)*(L/nh);
  endfor
endfunction
function[x]=maillevariable(nh)
  global L;
  for i=1:nh+1
    x(i)=(L/2)*(1 - @cos((i-1)*pi/nh));
  endfor
endfunction
function K=rigid_regulier(nh)
    global l;
    if (l==1)
      K=zeros(2*nh+2,2*nh+2);
      x=mailleregulier(nh);
      j=1;
      for i=1:2:2*nh-1
        K([i i+1 i+2 i+3],[i i+1 i+2 i+3])+=Kele(x(j),x(j+1),@EI,@EA);
        j=j+1;
      endfor
    else
      K=zeros(nh+1,nh+1);
      x=mailleregulier(nh);
      for i=1:1:nh
        K([i i+1],[i i+1])+=Kele(x(i),x(i+1),@EI,@EA);
      endfor
    endif 
endfunction
function K=rigid_variable(nh)
    global l;
    if (l==1)
      K=zeros(2*nh+2,2*nh+2);
      x=maillevariable(nh);
      j=1;
      for i=1:2:2*nh-1
        K([i i+1 i+2 i+3],[i i+1 i+2 i+3])+=Kele(x(j),x(j+1),@EI,@EA);
        j=j+1;
      endfor
    else 
      K=zeros(nh+1,nh+1);
      x=maillevariable(nh);
      for i=1:1:nh
        K([i i+1],[i i+1])+=Kele(x(i),x(i+1),@EI,@EA);
      endfor
    endif 
endfunction   
function f=force_regulier(nh)
     global l;
     if(l==1)
        f=zeros(2*nh+2,1);
        x=mailleregulier(nh);
        j=1;
        for i=1:2:2*nh-1
          f([i i+1 i+2 i+3])+=fele(x(j),x(j+1),@q)';
          j=j+1;
        endfor
     else 
        f=zeros(nh+1,1);
        x=mailleregulier(nh);
        for i=1:1:nh
          f([i i+1])+=fele(x(i),x(i+1),@q)';
        endfor
     endif 
endfunction
function f=force_variable(nh)
  global l;
  if(l==1)
     x=maillevariable(nh);
     f=zeros(2*nh+2,1);
     j=1;
     for i=1:2:2*nh-1
       f([i i+1 i+2 i+3])+=fele(x(j),x(j+1),@q,@S)';
       j=j+1;
     endfor 
  else
     x=maillevariable(nh);
     f=zeros(nh+1,1);
     for i=1:1:nh
         f([i i+1])+=fele(x(i),x(i+1),@q,@S)';
     endfor
  endif
endfunction
function U=Umaillereg(nh)
  global F0;
  global l;
  if (l==1)
    K=rigid_regulier(nh);
    f=force_regulier(nh);
    U=zeros(2*nh+2,1);
    Q=zeros(2*nh+2,1);
    Q(2*nh+1)=F0;
    F=f+Q;
    F([1 2])=[];
    K(1,:)=[];
    K(:,1)=[];
    K(2,:)=[];
    K(:,2)=[];
    U(3:2*nh+2)=K\F;
  else 
    f=force_variable(nh);
    K=rigid_variable(nh);
    U=zeros(nh+1,1);
    N=zeros(nh+1,1);
    N(nh+1)=F0;
    F=f+N
    F([1])=[];
    K(1,:)=[];
    K(:,1)=[];
    U(2:nh+1)=K\F;
  endif
endfunction 
function  U=Umaillevar(nh)
  global F0;
  global l;
  if (l==1)
    f=force_variable(nh);
    K=rigid_variable(nh);
    U=zeros(2*nh+2,1);
    Q=zeros(2*nh+2,1);
    Q(2*nh+1)=F0;
    F=f+Q
    F([1 2])=[];
    K(1,:)=[];
    K(:,1)=[];
    K(2,:)=[];
    K(:,2)=[];
    U(3:2*nh+2)=K\F;
  else 
    f=force_variable(nh);
    K=rigid_variable(nh);
    U=zeros(nh+1,1);
    N=zeros(nh+1,1);
    N(nh+1)=F0;
    F=f+N
    F([1])=[];
    K(1,:)=[];
    K(:,1)=[];
    U(2:nh+1)=K\F;
  endif
endfunction
function y=v(x)
  global E;
  global F0;
  global I;
  global L;
  global q0;
  y= (F0/(6*E*I))*(x**2)*(3*L -x) + ((q0*L**4)/(120*E*I))*( (20/L**2)*x**2 - (10/L**3)*x**3 + x**5/L**5);
endfunction
function y=dv(x)
  global I;
  global E;
  global F0;
  global q0;
  global L;
  y= (F0/(6*E*I))*((2*x)*(3*L -x)-x**2) + ((q0*L**4)/(120*E*I))*( (40/L**2)*x - (30/L**3)*x**2 + 5*x**4/L**5);
endfunction
function y=Ty(x)
  global F0;
  global q0;
  global L;
  y=F0+(q0*L/2)*(1-(x/L)**2);
endfunction
function y=Mfz(x)
  global F0;
  global q0;
  global L;
  y=F0*(L-x) +((q0*L**2)/6)*(2 -3*x/L +(x/L)**3);
endfunction
function Q=Qregulier(nh)
  U=Umaillereg(nh);
  K=rigid_regulier(nh);
  f=force_regulier(nh); 
  Q=K*U -f;
endfunction
function N=Nregulier(nh)
  U=Umaillereg(nh);
  K=rigid_regulier(nh);
  f=force_regulier(nh); 
  N=K*U -f;
endfunction
function N=Nvariable(nh)
  U=Umaillevar(nh);
  K=rigid_variable(nh);
  f=force_variable(nh); 
  N=K*U -f;
endfunction
function Q=Qereg(e,nh)
  x=mailleregulier(nh);
  U=Umaillereg(nh);
  Ke=Kele(x(e),x(e+1),@EI,@EA);
  fe=fele(x(e),x(e+1),@q,@S)';
  Q=Ke*U([2*(e-1)+1 2*(e-1)+2 2*(e-1)+3 2*(e-1)+4]) - fe;
endfunction
function Q=Qevar(e,nh)
  x=maillevariable(nh);
  U=Umaillevar(nh);
  Ke=Kele(x(e),x(e+1),@EI,@EA);
  fe=fele(x(e),x(e+1),@q,@S)';
  Q=Ke*U([2*(e-1)+1 2*(e-1)+2 2*(e-1)+3 2*(e-1)+4]) - fe;
endfunction
function Tyreg(nh)
  x=mailleregulier(nh);
  for i=1:length(x)-1
    Q=Qereg(i,length(x));
    T(i)=-Q(1);
  endfor
  T(length(x))=Q(3);
  plot(x(3:length(x)),T(3:length(x)))
endfunction
function Tyvar(nh)
  x=maillevariable(nh);
  for i=1:length(x)-1
    Q=Qevar(i,length(x));
    T(i)=-Q(1);
  endfor
  T(length(x))=Q(3);
  plot(x(3:length(x)),T(3:length(x)))
endfunction
function Mfzreg(nh)
  x=mailleregulier(nh);
  for i=1:length(x)-1
    Q=Qereg(i,length(x));
    M(i)=-Q(2);
  endfor
  M(length(x))=Q(4);
  plot(x(3:length(x)),M(3:length(x)))
endfunction
function Mfzvar(nh)
  x=maillevariable(nh);
  for i=1:length(x)-1
    Q=Qevar(i,length(x));
    M(i)=-Q(2);
  endfor
  M(length(x))=Q(4);
  plot(x(4:length(x)),M(4:length(x)))
endfunction
if (l==1)
  X=[0:0.01:1];
  figure('name','La déformée pour le deux types de maillages');
  subplot(3,2,1);
  for i=1:length(X)
    V(i)=1000*v(X(i));
  endfor
  plot(X,V,'b')
  hold on;
  U10=Umaillereg(10);
  X10=mailleregulier(10);
  plot(X10,1000*U10(1:2:22),'g');
  hold on;
  U20=Ureg(20);
  X20=maillereg(20);
  plot(X20,1000*U20(1:2:42),'g')
  U50=Ureg(50);
  X50=maillereg(50);
  plot(X50,1000*U50(1:2:102),'m')
  legend({'exacte','nh=10','nh=20','nh=50'});
  title('La déformée de la poutre pour le maillage régulier');
  subplot(3,2,2);
  for i=1:length(X)
    V(i)=1000*v(X(i));
  endfor
  plot(X,V,'b')
  hold on;
  U5=Uvar(5);
  X5=maillevar(5);
  plot(X5,1000*U5(1:2:12),'r');
  hold on;
  U10=Umaillevar(10);
  X10=maillevariable(10);
  plot(X10,1000*U10(1:2:22),'g')
  U20=Uvar(20);
  X20=maillevar(20);
  plot(X20,1000*U20(1:2:42),'m')
  legend({'exacte','nh=5','nh=10','nh=20'});
  title('La déformée de la poutre en pour le maillage variable');
  subplot(3,2,5)
  x=0:0.01:1;
  for i=1:length(x)
    T(i)=Ty(x(i));
  endfor
  plot(x,T)
  hold on;
  Tyreg(10)
  hold on;
  Tyvar(20)
  legend({'la solution exacte','nh=10','nh=20'});
  title('Effort tranchant pour un maillage régulier');
  subplot(3,2,6)
  x=0:0.01:1;
  for i=1:length(x)
    T(i)=Ty(x(i));
  endfor
  plot(x,T)
  hold on;
  Tyvar(10)
  hold on;
  Tyvar(20)
  legend({'la solution exacte','nh=10','nh=20'});
  title('Effort tranchant pour un maillage variable');
  subplot(3,2,3)
  x=0:0.001:1;
  for i=1:length(x)
    Mf(i)=Mfz(x(i));
  endfor
  plot(x,Mf)
  hold on;
  Mfzreg(10)
  hold on;
  Mfzreg(20)
  hold on;
  legend({'la solution exacte','nh=10','nh=20'});
  title('Moment fléchissant pour un maillage régulier');
  subplot(3,2,4)
  x=0:0.001:1;
  for i=1:length(x)
    Mf(i)=Mfz(x(i));
  endfor
  plot(x,Mf)
  hold on;
  Mfzvar(10)
  hold on;
  Mfzvar(20)
  legend({'exacte','nh=10','nh=20'});
  title('Moment flechissant pour un maillage variable'); 
  rigid_variable(nh)
  rigid_regulier(nh)
  force_variable(nh)
  force_regulier(nh)
  Umaillevar(nh)
  Umaillereg(nh)
else 
  rigid_variable(nh)
  rigid_regulier(nh)
  force_variable(nh)
  force_regulier(nh)
  Umaillevar(nh)
  Umaillereg(nh)
endif





