%Works1FisherTalbot32MittalJain.m
%Thomas to run Fisher equation%U_t = U_xx -bu^2 + au with boundary conditions lim x->-inf = 0.5 , lim x->
%lim x->-inf = 0.5 , lim x->+ inf = 0.Withoutmiltipreciaion. With initial conditions as below.
%The values of the initial guess are being updated as t increases with time
%steps of 0.5 this is the tt loop below. Otherwise the algoritgm fails for
%large t. The equation together with all the parameters are taken from a
%the paper"Numerical solutions of nonlinear Fisher's reaction Diffusion
%Equations with modified cubic B spline collocation method.Example3.
% 
%

clear all
clc
n=55;%Talbot parameter always odd.
%mp.Digits(2*n);
t =.2;
TT=5;
T=t*TT;
aa=.5;
bb=1;
cc =1;

j =(1:n-1);
z = (sqrt(-1)*(2*pi*(2*j/n-1)));
S =( z./(1-exp(-z)));
Sprime =( (1-(1+z).*exp(-z))./(1-exp(-z)).^2);
%Setting up the Talbot "Weights".
L=300;% Length of bar.
x_start= -L;%Begining of bar.
x_end = L;%end of bar.
h=1;%nodal distance.

N=L/h+1;%no of points along the bar.
xx=(N+1)/2;% middle point.

x=((-L/2:h:L/2));%points on the x axis.


%Diagonals for the tridiagonla matrix for Thomas,See TDMA functon file for.
%further info.
a=ones(N-2,1);

c=ones(N-2,1);

d=zeros(N-2,1);

% U0 =(1./(1 +exp((sqrt(C/6)*x))).^2);%Initail condition.
%  P = ((aa/(24*cc)));
% 
%   A= (-1/4*(aa/bb));
%   C = ((sech(-sqrt(P)*x)).^2);
%   D = ((-2*tanh(-sqrt(P)*x))-2);
% 
%  
% 
%   U0=(A*(C+D));
for C= 4:4:20
%U0 = 1./(2*(exp(x./2./C)+1));
U0 = 1./(2*(exp(x./2./C)+1))+1/C^2.*exp(x./2./C)./(1+ exp(x./2./C)).^2.*log(sqrt(2)*exp(x./8./C)./...
    (1+exp(x./2./C)).^2);
  U_old = U0;


% hold on
% 
% 
%hold on
 for first = 1:TT
% 
  for pp=1:10%Iteration loop for nonlinear form.
%  %Semi Direct iteration . 
%   
     for k=j %LT loop.
%              
                  d=(-h^2.*(U0((2:N-1))));%The rhs of the tridiag matrix equation.
                   b=(-2-h^2.*(S(k) +bb*U_old(2:N-1) -aa));%This is the alpha in the tridiag matrix.
%             
                d(1)=( d(1)-0.5/S(k));%L.H node.
                d(N-2)= (d(N-2)-0);%R.H node
               F(k,:) = (TDMAsolver(a,b,c,d));%TDMA solver.

     w=(exp(S*t).*Sprime);%Talbot" Weights".
 
          %Inversion loop.
            G(k,:)=(w(k)*F(k,:));
     end % second k loop
        ff =( 2/n*real(sum(G)));%final part of talbot inversion process.
%      
        u=([0.5 ff 0]);
%   
        U_old = u;
% 
  end%pp loop
 U0=u;
% size(F)
% T=first*t
 T=first*t; 
% %numerical=Uold(xx)
 y=-L/2:10:L/2;
  %Exact = ( -1/4*(aa/bb)*(sech((-sqrt(P)*y)+5*aa*T/12).^2-2*tanh((-sqrt(P)*y)+5*aa*T/12)-2));
%  xx=285:1:310;
 %Maxerror=max(abs((U_old(xx)-Exact(xx))./Exact(xx))*100)




 

% exact = 1./(1+exp((sqrt(C/6))*x - (5*C/6)*T)).^2;
% % Exact=exact(xx)
% % P=abs(((Exact-numerical)/Exact)*100);
% % Percenterror = vpa(P,5)
%   plot(x,U_old,'-r',x,Exact,'g')
%  
% plot(x,U_old)
 
 end%first loop
z=y-C*T;
Exact= 1./(2*(exp(z./2./C)+1))+1/C^2.*exp(z./2./C)./(1+ exp(z./2./C)).^2.*log(sqrt(2)*exp(z./8./C)./...
    (1+exp(z./2./C)).^2); 
hold on
 if C==4
 plot(x,U_old,'r',y,Exact,'*')
end
if C==8
 plot(x,U_old,'r',y,Exact,'*')
 
end
if C==12
 plot(x,U_old,'r',y,Exact,'*')
 
end
if C==16
 plot(x,U_old,'r',y,Exact,'*')
 
end
hold off
end %C loop
 

%plot(x,U_old,'r', y,Exact,'*')




% plot(x,U_old,'-r',x,Exact,'g')
%title('Fisher Talbot p=10,t=5,dt=.2 n=55')
%shg
  %plot(x,Exact)
% size(x)
% size(u)
% hold off
%  xx=150:1:175;
%  Maxerror=max(abs((U_old(xx)-Exact(xx))./Exact(xx))*100)
 % plot(x,U_old,'-r',x,Exact,'g') 
%abs((U_old(44)-Exact(44)))
% x = linspace(30,5,30);
% y = exp(x/10).*sin(4*x);
% plot(x,y,'-o')