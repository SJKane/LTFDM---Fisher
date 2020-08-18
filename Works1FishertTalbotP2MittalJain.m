%Works1FisherTalbotP2MittalJain.m
%Thomas to run Fisher equation u_t=u_xx+Cu(u-1) with 
%u_0=(1 =exp(x))^-2 and1./(1+exp((sqrt(C/6))*x - (5*C/6)*t)).^2; .
%Inversion is via talbot.Thisworks.This is not working.

clear all
clc
n=55;%Talbot parameter always odd.
t =.2;
TT=25;
%C=1;
r=1;
C=1;
T=t*TT
%Setting up the Talbot "Weights".
L=60;% Length of bar.
x_start= -L;%Begining of bar.
x_end = L;%end of bar.
h=.25;%nodal distance.
N=L/h+1;%no of points along the bar.
xx=(N+1)/2;% middle point.

x=((-L/2:h:L/2));%points on the x axis.

%Diagonals for the tridiagonla matrix for Thomas,See TDMA functon file for.
%further info.
a=ones(N-2,1);
c=ones(N-2,1);
d=zeros(N-2,1);
U0 =(1./(1 +exp((sqrt(C/6)*x))).^2);%Initail condition.
%U0=0.5*(1 + tanh(-x/4));
U_old = U0;
j =(1:n-1);
z = (sqrt(-1)*(2*pi*(2*j/n-1)));
S =( z./(1-exp(-z)));
Sprime =( (1-(1+z).*exp(-z))./(1-exp(-z)).^2);




for first = 1:TT

  for pp=1:10%Iteration loop for nonlinear form.
 %Semi Direct iteration . 
  
for k=j %LT loop.
             
                 d=- h^2.*(U0((2:N-1)));%The rhs of the tridiag matrix equation.
                 b=-2 -h^2.*(S(k)+ r*U_old(2:N-1) -r);%This is the alpha in the tridiag matrix.
               d(1)=( d(1)-1/S(k));%L.H node.
               d(N-2)= (d(N-2)-0);%R.H node
               F(k,:) = (TDMAsolver(a,b,c,d));%TDMA solver.
      %end %k loop
    w=(exp(S*t).*Sprime);%Talbot" Weights".

      % for k=j %Inversion loop.
           G(k,:)=(w(k)*F(k,:));
       end % second k loop
       ff =( 2/n*real(sum(G)));%final part of talbot inversion process.
     
       u=([1 ff 0]);
  
       U_old = u;

  end%pp loop
U0=u;

 T=first*t;
 hold on
if T==1 
   plot(x,U_old,'r')
   axes('box','on')
end

 if T==2
  plot(x,U_old,'r')
 end
 
if T==3
 plot(x,U_old,'r')
end 
if T==4
 plot(x,U_old,'r')
end 
if T==5
 plot(x,U_old,'r')
end 
%  plot(x,u,'r',x,exact,'g')

hold off

% Exact = 1./(1+exp((sqrt(C/6))*x - (5*C/6)*T)).^2;

end%first loop

% plot(x,u,'r',x,Exact,'g')
%plot(x,u)
