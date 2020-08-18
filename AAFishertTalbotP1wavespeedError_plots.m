%Works1FisherTalbot32MittalJain.m Error plot
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
n=555;%Talbot parameter always odd.
Time=1;
t =.1;%Time step
TT=Time/.2;
T=t*TT;%= Time
%Set constants.
aa=.5;
bb=1;
cc =1;
%Setting up the Talbot "Weights".
j =(1:n-1);
z = (sqrt(-1)*(2*pi*(2*j/n-1)));
S =( z./(1-exp(-z)));
Sprime =( (1-(1+z).*exp(-z))./(1-exp(-z)).^2);

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
C=20%Wave speed

U0 = 1./(2*(exp(x./2./C)+1))+1/C^2.*exp(x./2./C)./(1+ exp(x./2./C)).^2.*log(sqrt(2)*exp(x./8./C)./...
    (1+exp(x./2./C)).^2);
  U_old = U0;



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

 
 end%first loop
z=x-C*T;
Exact= 1./(2*(exp(z./2./C)+1))+1/C^2.*exp(z./2./C)./(1+ exp(z./2./C)).^2.*log(sqrt(2)*exp(z./8./C)./...
    (1+exp(z./2./C)).^2); 
% hold on
%  if C==10
%      y=-L/2:3:L/2;
%      z=y-C*T;
% Exact= 1./(2*(exp(z./2./C)+1))+1/C^2.*exp(z./2./C)./(1+ exp(z./2./C)).^2.*log(sqrt(2)*exp(z./8./C)./...
%     (1+exp(z./2./C)).^2); 
% 
%  plot(x,U_old,'r',y,Exact,'*')
% axis([-L/2 L/2 -0.03 0.503])
% end
% if C==8
%  plot(x,U_old,'r',y,Exact,'*')
%  
% end
% if C==10
%  plot(x,U_old,'r',y,Exact,'*')
%  
% end
% if C==16
%  plot(x,U_old,'r',y,Exact,'*')
%  
% end
% if C==20
%  plot(x,U_old,'r',y,Exact,'*')
%  
% end
% legend({'Numerical','Exact'})
%  hold off
% end %C loop
%  
%Tables of results.  
% z=x-C*T;
% Exact= 1./(2*(exp(z./2./C)+1))+1/C^2.*exp(z./2./C)./(1+ exp(z./2./C)).^2.*log(sqrt(2)*exp(z./8./C)./...
%     (1+exp(z./2./C)).^2); 
% for xx=(N+1)/2:1:N%these values are for +ve x values from zero  to L 
%  %for xx=2:1:(N+1)/2 %these values are for -ve x values from zero  to L/2 
% %       
% %  
%    x=-L/2+h*(xx-1)
% Num=U_old(xx);
% exact=(Exact(xx));
% Absolute=abs(exact-Num);
% 
% relerror=  (abs((Num-exact)./exact));
% location1 = [x];
% location2 = [Num];
% location3 = [exact];
% location4 = [Absolute];
% location5 = [relerror]
% T = table;
% T.x=location1
% T.Numerical = location2';
% T.Exact= location3';
% T.AbsoluteError = location4';
% T.RelError = location5'
% end
 z=x-C*T;
 Exact= 1./(2*(exp(z./2./C)+1))+1/C^2.*exp(z./2./C)./(1+ exp(z./2./C)).^2.*log(sqrt(2)*exp(z./8./C)./...
    (1+exp(z./2./C)).^2); 

 xx=(N+11)/2:1:N-11;
  x=-L/2+h*(xx-1);
Num=U_old(xx);
exact=(Exact(xx));
Absolute=abs(exact-Num);
percenterror=  (abs((Num-exact)./exact)*100);
plot(x,percenterror)














