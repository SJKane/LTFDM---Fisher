%Works1FisherTalbotP2MittalJainPoints.m
%Thomas to run Fisher equation u_t=u_xx+Cu(u-1) with 
%u_0=(1 =exp(x))^-2 and1./(1+exp((sqrt(C/6))*x - (5*C/6)*t)).^2; .
%Inversion is via talbot.Thisworks

clear all
clc
n=755;%Talbot parameter always odd.
t =.5;
TT=20;
tt=t*TT
%Setting up the Talbot "Weights".


%for t =0.1:.1:.5
j =(1:n-1);
z = (sqrt(-1)*(2*pi*(2*j/n-1)));
S =( z./(1-exp(-z)));
Sprime =( (1-(1+z).*exp(-z))./(1-exp(-z)).^2);
%TT=10;
C=1;%Constant in Fisher equation above.

%Set up FD Scheme and solve in Laplace Transform Space.

L=60;% Length of bar.
x_start= -L;%Begining of bar.
x_end = L;%end of bar.
h=.25;%nodal distance.
N=L/h+1;%no of points along the bar.
%xx=(N-5)/2:1:(N+5)/2% middle point.
%xx=(N+1)/2
x=((-L/2:h:L/2));%points on the x axis.

%Diagonals for the tridiagonla matrix for Thomas,See TDMA functon file for
%further info.
a=ones(N-2,1);%(N-2:1) i.e. N-2 rows by 1 column
c=ones(N-2,1);
d=zeros(N-2,1);

%nn = 1:N;%Interger points.
U0 =(1./(1 +exp((sqrt(C/6).*x))).^2);%Initail condition
U_old = U0;%First approximation.~
hold on
for first = 1:TT
 %U0=Uold;
  for p=1:10%Iteration loop for nonlinear form
 %indirect iteration    
for k=j %LT loop
             
                 d=- h^2.*(U0((2:N-1)));%The rhs of the tridiag matrix equation.Matrix of
               

%U0 from U(0)(2) to U(0)(N-1) matrix demension 1,N-1 semi direct.
 b=(-2 -h^2.*(S(k) +U_old(2:N-1) -C));%This is the alpha in the tridiag matrix semi direct
               d(1)=( d(1)-1/S(k));%L.H node.
               d(N-2)= (d(N-2)-0);%R.H node
               F(k,:) = (TDMAsolver(a,b,c,d));%TDMA solver.
     % end %k loop
    w=(exp(S*t).*Sprime);%Talbot" Weights".

        %for k=j %Inversion loop
           G(k,:)=(w(k)*F(k,:));
        end %k loop
       tot =( 2/n*real(sum(G)));%final part of talbot inversion process.
     U_old=([1 tot 0]);
  end%Non linear p loop

U0=U_old;T=t*first;
Exact  = 1./(1+exp((sqrt(C/6))*x - (5*C/6)*T)).^2;
count=1;
for xx=(N-5)/2:1:(N+5)/2;
numerical=U_old(xx);
exact=Exact(xx);
percenterror(count)=(abs(((exact-numerical)/exact)*100));
count=count+1;
max(percenterror)
end%xx loop
plot(x,U_old,'r',x,Exact,'g')
%first loop or t loop
end
hold off
%for xx=(N-5)/2:1:(N+5)/2





%R=abs(exact-numerical)
%Percenterror = vpa(P,5)


%e
 %plot(x,Uold,'or',x,exact,'*g')
 
