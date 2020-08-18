%Works1FisherTalbotP2MittalJainPoints.m
%Thomas to run Fisher equation u_t=u_xx+Cu(u-1) with 
%u_0=(1 =exp(x))^-2 and1./(1+exp((sqrt(C/6))*x - (5*C/6)*t)).^2; .
%Inversion is via talbot.Thisworks
%Increasing the Talbot weights irons out the pertubations.
%The equation is taken from a
%the paper"Numerical solutions of nonlinear Fisher's reaction Diffusion
%Equations with modified cubic B spline collocation method.Example2.Here we
%use different parameters and a diffrent space domain for x results are
%comapred with the exact solution.
clear all
clc
format long
n=555;%Talbot parameter always odd.
t =.2;
TT=5;
%T=t*TT;
%Setting up the Talbot "Weights".


%for t =0.1:.1:.5
j =(1:n-1);
z = (sqrt(-1)*(2*pi*(2*j/n-1)));
S =( z./(1-exp(-z)));
Sprime =( (1-(1+z).*exp(-z))./(1-exp(-z)).^2);
%TT=10;
C=1;%Constant in Fisher equation above.

%Set up FD Scheme and solve in Laplace Transform Space.

L=80;% Length of bar.
x_start= -L;%Begining of bar.
x_end = L;%end of bar.
h=1;%nodal distance.
N=L/h+1;%no of points along the bar.
xx=(N-5)/2:1:(N+5)/2;%; middle point.
%xx=(N+1)/2
x=((-L/2:h:L/2));%points on the x axis.
%x=0:h:L/2;
%Diagonals for the tridiagonla matrix for Thomas,See TDMA functon file for
%further info.
a=ones(N-2,1);%(N-2:1) i.e. N-2 rows by 1 column
c=ones(N-2,1);
d=zeros(N-2,1);

%nn = 1:N;%Interger points.
U0 =(1./(1 +exp((sqrt(C/6).*x))).^2);%Initail condition


U_old = U0;%First approximation.~
%hold on
for first =1:TT
  for p=1:10%Iteration loop for nonlinear form the number of iterative loops enough for convergence
      
for k=j %LT loop
             %direct iteration.
% %                
%                  d = h^2*(C*U_old(2:N-1).^2-U0(2:N-1));
%                  b=-(2+h^2*S(k)+h^2*C)*a;
                  
 %semi direct.
%U0 from U(0)(2) to U(0)(N-1) matrix demension 1,N-1
   d=- h^2.*(U0((2:N-1)));%The rhs of the tridiag matrix equation.Matrix of
                 b=(-2 -h^2.*(S(k) +C*U_old(2:N-1) -C));%This is the alpha in the tridiag matrix semi direct
               d(1)=( d(1)-1/S(k));%L.H node.
               d(N-2)= (d(N-2)-0);%R.H node
               F(k,:) = (TDMAsolver(a,b,c,d));%TDMA solver.
     % end %k loop
    w=(exp(S*t).*Sprime);%Talbot " Weights".

        %for k=j %Inversion loop
           G(k,:)=(w(k)*F(k,:));
        end %k loop
       ff =( 2/n*real(sum(G)));%final part of talbot inversion process.
     
       u=([1 ff 0]);
  
       U_old =( u);


  end%Non linear loop
%    [aaa bbb] =  size(U_old);
 %U0=U_old; 
%   for i=1:1:bbb
%       if U_old(1,i) <=5*10^(-5);
%           U_old(1,i) = 0;
%       end
%    end

        %[aaa,bbb ] =  size(U_old);
   
%     for i=1:1:bbb
%         if U_old(1,i) >=1+5*10^(-5);
%            U_old(1,i) = 1;
%         end
%    end
%                   
%for i=1:1:
U0=U_old;
%
%hold on
size(U_old);
%T=t*first;
v= 5/sqrt(6);
T=t*first;
Y=-L:.5:L;
%Y=0:.5:L;
Exact  = 1./(1+exp((sqrt(1/6))*Y - (5*1/6)*T)).^2;

%plot(x,U_old,'r',x,Exact)

hold on
if T==1 
   plot(x,U_old,'r',Y,Exact,'b*')
end
 
% if T==2
%  plot(x,U_old,'r',Y,Exact,'b*')
% 
% end
% if T==3
%  plot(x,U_old,'r',Y,Exact,'b*')
%  
% end
% if T==4
%  plot(x,U_old,'r',Y,Exact,'b*')
%  
% end
% if T==5
%  plot(x,U_old,'r',Y,Exact,'b*')
% 
% end
hold off
end%first loop or t loop
T
legend('Talbot','Exact')
 
% Exact  = 1./(1+exp((sqrt(C/6))*x - (5*C/6)*T)).^2;
% for xx=10:20:100;
% x=xx*h
% Num=U_old(xx)
% exact=Exact(xx)
%   relerror=  (abs((U_old(xx)-Exact(xx))./Exact(xx)))
% end
%hold off
% max(abs( U_old(xx)- Exact(xx)));
%R=abs(exact-numerical)
%Percenterror = vpa(P,5)
 
%t*first
%e
% plot(x,U_old,'or',x,exact,'*g')
 
