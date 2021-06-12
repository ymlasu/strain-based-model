function[sig,tau,sx,tx]=nonproportional(K, n, r, sig_y, sa, ta, phy, E, Epoisson)
sig=[0];
tau=[0];
sig_ai=[0 sig_y:r:800];
eps_ai=sig_ai/E+(sig_ai/K).^(1/n);
% phy=phy*pi/180;
for i=2:length(sig_ai)
    R(i-1)=sig_ai(i);
    H(i-1)=((eps_ai(i)-eps_ai(i-1))/(sig_ai(i)-sig_ai(i-1))-1/E)^(-1);
end
center=zeros(length(H),2);

for tt=2:1440
    xt=tt*2*pi/360;
    sx(1)=0;
    tx(1)=0;
    sx(tt)=sa*sin(xt);
    tx(tt)=ta*sin(xt-phy);
    if xt-phy<0
        tx(tt)=0;
    end
    i=1;
    while (sig(tt-1)-center(i,1))^2+(sqrt(3)*tau(tt-1)-center(i,2))^2>R(i)^2
        i=i+1;
    end
    dsx=sx(tt)-sx(tt-1);
    dtx=tx(tt)-tx(tt-1);
    dseq=(dsx^2+3/4*dtx^2/(1+Epoisson)^2)^0.5;
    dseqp=dseq-dseq*H(i)/E;
    dsigeq=sqrt(sig(tt-1)^2+3*tau(tt-1)^2);
    dlamda=dseqp/dsigeq;
    if i==1
        dlamda=0;
    end
    dsig=E*(dsx-dlamda*sig(tt-1));
    dtau=(dtx-3*dlamda*tau(tt-1))*E/2/(1+Epoisson);
    
%     global s t deps dgama h
%     s=sig(tt-1);
%     t=tau(tt-1);
%     deps=dsx;
%     dgama=dtx;
%     h=H(i);
%     x=fsolve(@(x) deltasig(x,E,Epoisson),[0;0]);
%     dsig=x(1);
%     dtau=x(2);
    sig(tt)=sig(tt-1)+dsig;
    tau(tt)=tau(tt-1)+dtau;
    if i>1
        if (sig(tt)-center(i,1))^2+(sqrt(3)*tau(tt)-center(i,2))^2>R(i)^2
            [x0,y0]=point(sig(tt-1),sig(tt),sqrt(3)*tau(tt-1),sqrt(3)*tau(tt),center(i,1),center(i,2),R(i));
            for k=1:i-1
                center(k,1)=(R(k+1)-R(k))/R(k+1)*(x0-center(k+1,1))+center(k+1,1);
                center(k,2)=(R(k+1)-R(k))/R(k+1)*(y0-center(k+1,2))+center(k+1,2);
            end
        end
    end
end
end


% function eq=deltasig(x,E,Epoisson)
% global s t deps dgama h
%     eq(1)=deps-x(1)/E-(sqrt((s+x(1))^2+3*(t+x(2))^2)-sqrt(s^2+3*t^2))/h/sqrt(s^2+3*t^2)*s;
%     eq(2)=dgama-(1+Epoisson)/E*x(2)+3/2*(sqrt((s+x(1))^2+3*(t+x(2))^2)-sqrt(s^2+3*t^2))/h/sqrt(s^2+3*t^2)*t;
% end


function[x,y]=point(x1,x2,y1,y2,c1,c2,r)
syms x y
[Sx,Sy]=solve((x1-x2)*(y-y1)-(y1-y2)*(x-x1)==0,(x-c1)^2+(y-c2)^2-r^2==0);
for k=1:length(Sx)
    if x1<x2
        if Sx(k)<x2 && Sx(k)>x1
           x=Sx(k);
        end
    else if x1==x2
            x=Sx(k);
        else if Sx(k)<x1 && Sx(k)>x2
                x=Sx(k);
            end
        end

    end
    
    if y1<y2
        if Sy(k)<y2 && Sy(k)>y1
           y=Sy(k);
        end
    else if y1==y2
        y=Sy(k);
    else if Sy(k)<y1 && Sy(k)>y2
            y=Sy(k);
        end
        end
    end
end
end