%Works1FisherTalbotP2MittalJain.m
%Non dimensionalised so that t=t'/s and
% x=x'/sqrt(s/v) Thomas to run Fisher equation u_t=v*u_xx+*su(u-1) With 
%A perturb solution to test for various wave speeds C in the code .


clear all
clc
n=555;%Talbot parameter always odd.
TT=4;
v=1;%diffusion constant.
s=1;%Intrinisic growth rate
t =.2/s;%Time non dimensionalised
j =(1:n-1);
z = (sqrt(-1)*(2*pi*(2*j/n-1)));
S =( z./(1-exp(-z)));
Sprime =( (1-(1+z).*exp(-z))./(1-exp(-z)).^2);

L=220;% Length of bar.
h=1;%nodal distance.
N=L/h+1;%no of points along the bar.
%xx=(N+11)/2% middle point.
x=((-L/2:h:L/2))./sqrt(s/v);%points on the x axis non dimensionlised
size(x);
L/2
%xx=700/2;
%Diagonals for the tridiagonla matrix for Thomas,See TDMA functon file for.
%further info.
a=ones(N-2,1);
c=ones(N-2,1);
d=zeros(N-2,1);
%Pertubation solution.
%C=6;
 %z=(x-C*t);
for C=4:4:20
% U0= 1/(1+exp((1/3).*x)) +...
% 1/c^2.*exp(x./c)*(1+exp(x./c))^(-2)*ln(4*exp(x./c))./(1+exp(x./c))^(-2);
%U0 =(1./(1 +exp((sqrt(1/6)*x))).^2);
%plot(x,U0)
U0= 1./(1 + exp(x./C))+1/C^2*exp(x/C)./(1 + exp(x./C)).^2.*log(4*exp(x/C)./(1 + exp(x./C)).^2) ;
%   U0 =1./(1 +exp(((1/C)*x)))+1/C^2*exp(((1/C)*x)).*(1+exp(((1/C)*x))).^(-2).*log...
%       ((4*exp(((1/C)*x)))./(1+exp(((1/C)*x))).^2);
U_old = U0;


  % U0= 1./(1 + exp(x./C))+1/C^2*exp(x/C)./(1 + exp(x./C)).^2.*log(4*exp(x/C)./(1 + exp(x./C)).^2) ; 
  
    for first =1:TT

  for p=1:5%Iteration loop for nonlinear form.
 %Semi Direct iteration . 
  
for k=j %LT loop.
             
                 d=- h^2.*(U0((2:N-1)));%The rhs of the tridiag matrix equation.
                 b=(-2 -h^2.*(S(k) +U_old(2:N-1) -1));%This is the alpha in the tridiag matrix.
               d(1)=( d(1)-1./S(k));%L.H node.
               d(N-2)= (d(N-2)-0);%R.H node
               F(k,:) = (TDMAsolver(a,b,c,d));%TDMA solver.
      
    w=(exp(S*t).*Sprime);%Talbot" Weights".

      % for k=j %Inversion loop.
           G(k,:)=(w(k)*F(k,:));
       end %  k loop
       ff =( 2/n*real(sum(G)));%final part of talbot inversion process.
     
       u=([1 ff 0]);
  
       U_old = u;

  end%p loop
U0=U_old;

T=t*first;
%T=t; 
%numerical=U_old(xx);
Y=-L/2:4:L/2;
Z=Y-C*T;

Pertub = 1./(1 + exp(Z./C))+1/C^2*exp(Z/C)./(1 + exp(Z./C)).^2.*log(4*exp(Z/C)./(1 + exp(Z./C)).^2) ;
P=Pertub;
%  Exact = 1./(1+exp((sqrt(1/6))*x - (5*1/6)*T)).^2;
%  exact = 1./(1+exp((sqrt(1/6))*x - (5*1/6)*T)).^2;

%Pertub=1./(1 +exp(((1/C)*x)-T))+1/C^2*exp(((1/C)*x)-T)./(1 +exp(((1/C)*x)-T)).^2.*log((4*exp(((1/C)*x)-T))./(1+exp(((1/C)*x)-T)).^2);

% Exact= 1/(1+exp(z./C)) +...
% 1/C^2.*exp(z./C)*(1+exp(z./C))^(-2)*ln(4*exp(z./C))/(1+exp(z./C))^(-2);
%exact = 0.5*(1 + tanh(-x/4 +(5/8)*t));
% Exact=exact(xx)
%  P=abs(((Pertub-numerical)/Pertub)*100);
%  Percenterror = vpa(P,5)
% plot(x,u,'or',x,exact,'*g')
T;

%Z=x-C*T;
% 
% hold on  
% if C==4
%     plot(x,P)
% end
% if C==8
%     plot(x,P)
% end
% % 
% if C==12
%     plot(x,P)
% end
% if C==16
%     plot(x,P)
% end

end%first loop
hold on  
if C==4
    plot(x,U_old,Y,P,'*')
    
    %axes ('box','on')
 axis([-L/2 L/2 -0.03 1.03])
end
%axes ('box','on')
if C==4
    plot(x,U_old,Y,P,'*')
    
end
if C==8
    plot(x,U_old,Y,P,'*') 
end

if C==12
  plot(x,U_old,Y,P,'*')   
end
if C==16
     plot(x,U_old,Y,P,'*') 
end
if C==20
     plot(x,U_old,Y,P,'*') 
end
hold off  
legend({'Numericval','Exact'})
 
%hold on
%Exact=1./(1 +exp(((1/C)*x)-T))+1/C^2*exp(((1/C)*x)-T)./(1 +exp(((1/C)*x)-T)).^2.*log((4*exp(((1/C)*x)-T))./(1+exp(((1/C)*x)-T)).^2);

%Pertubation=Pertub(xx)
 % plot(x,u,'r',x,Pertub,'g')
  %shg
%Pertubation = 1./(1 + exp(Z./C))+1/C^2*exp(Z/C)./(1 + exp(Z./C)).^2.*log(4*exp(Z/C)./(1 + exp(Z./C)).^2) ;
end%C loop

%     +1/C^2*exp(((1/C)*x)-T)./(1 +...
%     exp(((1/C)*x)-T)).^2.*log((4*exp(((1/C)*x)-T))./(1+exp(((1/C)*x)-T)).^2); 
 % plot(x,u)%shg
% numerical=U_old(xx)
% Pertubation(xx);
%Exact=Exact(xx);
%P=abs(((Pertubation-numerical)/Pertubation)*100);
%P=abs(((Exact-numerical)/Exact)*100);
%Rel=abs((Pertubation(xx)-numerical)/Pertubation(xx))
%Reletiveerror=vpa(Rel,5)
 %Percenterror = vpa(P,5)

%  for xx=10:20:600;
% x=xx*h;
% Num=U_old(xx);
% exact=(Pertubation(xx));
% %exact=vpa(exact,5);
%  % relerror=  (abs((U_old(xx)-Exact(xx))./Exact(xx)));
%  relerror=  (abs((Num-exact)./exact));
% location1 = [x];
% location2 = [Num];
% location3 = [exact];
% location4 = [relerror];
% T = table;
% T.x=location1
% T.Numerical = location2';
% T.Exact= location3';
% T.Error = location4'
% end
%plot(x,u,'or',x,exact,'*g')

%plot(x,u,x,Pertubation)
shg
