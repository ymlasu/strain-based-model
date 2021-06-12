%this program is to calculate the critical plane under 
%multiaxial fatigue and the amplitude of the maximum shear strain amplitude
%and normal strain amplitude on the critical plane

%%%for arbitry strain path
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set the initial data for calculation:
%sa:: normal strain amplitude of the x direction
%ta:: shear strain amplitude of the xy direction
%phy:: phase diffrence between the normal strain and shear strain
%the stress history is assumed to be the sinesodial function,say,
%sx=sa*sin(x); tx=ta*sin(x+phy)
clear all
%input cyclic material properties
matpro=xlsread('mat-pro-304.xls',1);
global Faxial Baxial Eaxial Caxial Eyoung Seq
global Kaxial Naxial sapoisson
Eyoung=matpro(1,1);
Gshear=matpro(1,2);
Epoisson=matpro(1,3);

Kaxial=matpro(2,1);
Naxial=matpro(2,2);
Faxial=matpro(3,1);
Baxial=matpro(3,2);
Eaxial=matpro(3,3);
Caxial=matpro(3,4);
Ktor=matpro(4,1);
Ntor=matpro(4,2);
Ftor=matpro(5,1);
Btor=matpro(5,2);
Etor=matpro(5,3);
Ctor=matpro(5,4);
%%%%%%%%%%%%%%%%%%%%%
%input strain amplitude data
Totaldata=xlsread('strain-304-step.xls');
fid = fopen('fatipre.dat','w');
totalnumber=size(Totaldata);
foptions=optimset('TolFun',1e-8);
for j=1:totalnumber(1,1)/2
        sapeak=Totaldata(2*j-1,1);
    tapeak=Totaldata(2*j,1);
    sm=Totaldata(2*j-1,2);
    tm=Totaldata(2*j,2);
    for pp=1:sapeak
    sa(pp)=Totaldata(2*j-1,2*(pp+1));
    satime(pp)=Totaldata(2*j-1,2*(pp+1)-1);
    end
    for pp=1:tapeak
    ta(pp)=Totaldata(2*j,2*(pp+1));
    tatime(pp)=Totaldata(2*j,2*(pp+1)-1);
    end 
    sapoisson=(max(sa)-min(sa))/2;
  Sigma=fsolve(@osgood1,400,foptions);
    if sapoisson>0.000001
    Epr=(Epoisson*Sigma/Eyoung+0.5*(sapoisson-Sigma/Eyoung))/sapoisson;
    %Epr=0.35;
    else
    Epr=Epoisson;
    end
%phy=0.0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %strain components on critical plane is calculated, both angle
               itsum=1.0;
               itcond=1.0;
               Flife(itsum)=6;           
               Tenlimit(itsum)=Faxial*(Flife(itsum))^(Baxial);
               Torlimit(itsum)=Ftor*(Flife(itsum))^(Btor);
                    while itcond==1.0
      %find the critical plane angle--maximum shear strain
xa=0.0;
xt=0.0;
anglestep=360;
tstep=360;
tampmem=(max(ta)+min(ta))/2;
tmeanmem=tm;
sampmem=(max(sa)+min(sa))/2;
smeanmem=sm;
smaxmem=max(sa)+sm;
angmem=0.0;
for i=1:anglestep
   xa=i*2*pi/anglestep;
   timecountsa=1;
   timecountta=1;
   for k=1:tstep
      % xt=k*2*pi/tstep;
      xt=k/tstep*360;
       %judge the time is within wich section
       if xt>tatime(timecountta+1)
           timecountta=timecountta+1;
           tx=ta(timecountta)+(ta(timecountta+1)-ta(timecountta))/(tatime(timecountta+1)-tatime(timecountta))*(xt-tatime(timecountta));
       else
           tx=ta(timecountta)+(ta(timecountta+1)-ta(timecountta))/(tatime(timecountta+1)-tatime(timecountta))*(xt-tatime(timecountta));
       end
       if xt>satime(timecountsa+1)
           timecountsa=timecountsa+1;
           sx=sa(timecountsa)+(sa(timecountsa+1)-sa(timecountsa))/(satime(timecountsa+1)-satime(timecountsa))*(xt-satime(timecountsa));
           sy=-Epr*sx;
       else
           sx=sa(timecountsa)+(sa(timecountsa+1)-sa(timecountsa))/(satime(timecountsa+1)-satime(timecountsa))*(xt-satime(timecountsa));
           sy=-Epr*sx;
       end
       sxangle(k)=(sx+sy)/2+(sx-sy)/2*cos(2*xa)+tx/2*sin(2*xa);
       txangle(k)=-(sy-sx)*sin(2*xa)+tx*cos(2*xa);  
       sahydro(k)=sx*(1-2*Epr)/3;
   end
   smaxamp=(max(sxangle)-min(sxangle))/2;
   smean=(max(sxangle)+min(sxangle))/2;
   tmaxamp=(max(txangle)-min(txangle))/2;
   tmean=(max(txangle)+min(txangle))/2;
   angt=xa;
   smax=smaxamp+smean;
   sahydromem=(max(sahydro)-min(sahydro))/2;
   if abs(smaxamp)>sampmem
   %if abs(smax)>smaxmem
       sampmem=smaxamp;
       smeanmem=smean;
       tampmem=tmaxamp;
       tmeanmem=tmean;
       angmem=angt;
       smaxmem=smax;
       %damagemem=damagemax;
   end
end
   s1=Torlimit(itsum)/Tenlimit(itsum);
   if s1>2
       alf=0.0;
       beta=s1/2;
       K=(0.25*s1^2-1)*9/(1-2*Epr)^2;
   elseif s1<2
       A=(1+Epr)^2+4-s1^2-4*(1+Epr)^2/s1^2;
       B=2*(1-Epr^2);
       C=(1-Epr)^2+4*(1+Epr)^2/s1^2-4;
       alf=acos((-B+(B^2-4*A*C)^0.5)/2/A)/2;
       %if alf>45/180*pi
           %alf=45/180*pi;
           %end
       beta=(0.25*s1^2*cos(2*alf)^2+sin(2*alf)^2)^0.5;
       %alf=45/180*pi;
       %beta=1.0;
   end
   angcri1=angmem-alf;
   angcri2=angmem+alf;
   timecountsa=1;
   timecountta=1;
   for k=1:tstep
      % xt=k*2*pi/tstep;
      xt=k/tstep*360;
       %judge the time is within wich section
       if xt>tatime(timecountta+1)
           timecountta=timecountta+1;
           tx=ta(timecountta)+(ta(timecountta+1)-ta(timecountta))/(tatime(timecountta+1)-tatime(timecountta))*(xt-tatime(timecountta));
       else
           tx=ta(timecountta)+(ta(timecountta+1)-ta(timecountta))/(tatime(timecountta+1)-tatime(timecountta))*(xt-tatime(timecountta));
       end
       if xt>satime(timecountsa+1)
           timecountsa=timecountsa+1;
           sx=sa(timecountsa)+(sa(timecountsa+1)-sa(timecountsa))/(satime(timecountsa+1)-satime(timecountsa))*(xt-satime(timecountsa));
           sy=-Epr*sx;
       else
           sx=sa(timecountsa)+(sa(timecountsa+1)-sa(timecountsa))/(satime(timecountsa+1)-satime(timecountsa))*(xt-satime(timecountsa));
           sy=-Epr*sx;
       end
       sahydro(k)=sx*(1-2*Epr)/3;
       sxangle1(k)=(sx+sy)/2+(sx-sy)/2*cos(2*angcri1)+tx/2*sin(2*angcri1);
       txangle1(k)=-(sy-sx)*sin(2*angcri1)+tx*cos(2*angcri1);   
       sxangle2(k)=(sx+sy)/2+(sx-sy)/2*cos(2*angcri2)+tx/2*sin(2*angcri2);
       txangle2(k)=-(sy-sx)*sin(2*angcri2)+tx*cos(2*angcri2); 
      end
   smaxamp1=(max(sxangle1)-min(sxangle1))/2;
   smean1=(max(sxangle1)+min(sxangle1))/2;
   tmaxamp1=(max(txangle1)-min(txangle1))/2;
   tmean1=(max(txangle1)+min(txangle1))/2;
   sampmem1=smaxamp1;
   smeanmem1=sm/2+sm/2*cos(2*angcri1);
   tampmem1=tmaxamp1;
   tmeanmem1=tmean1;
   smaxamp2=(max(sxangle2)-min(sxangle2))/2;
   smean2=(max(sxangle2)+min(sxangle2))/2;
   tmaxamp2=(max(txangle2)-min(txangle2))/2;
   tmean2=(max(txangle2)+min(txangle2))/2;
   sampmem2=smaxamp2;
   smeanmem2=sm/2+sm/2*cos(2*angcri2);
   tampmem2=tmaxamp2;
   tmeanmem2=tmean2;
   if s1<2
   Seq1=(smaxamp1^2+tmaxamp1^2/s1^2)^0.5/beta/(1-smeanmem1/900);
   Seq2=(smaxamp2^2+tmaxamp2^2/s1^2)^0.5/beta/(1-smeanmem2/900);
   else
   Seq1=(smaxamp1^2+tmaxamp1^2/s1^2+K*sahydromem^2)^0.5/beta/(1-smeanmem1/900);
   Seq2=(smaxamp2^2+tmaxamp2^2/s1^2+K*sahydromem^2)^0.5/beta/(1-smeanmem2/900);  
   end
   if Seq1>Seq2
          sampmem=sampmem1;
          smeanmem=smeanmem1;
          tampmem=tampmem1;
          tmeanmem=tmeanmem1;
          angcri=angcri1;
          Seq=Seq1;
   else
          sampmem=sampmem2;
          smeanmem=smeanmem2;
          tampmem=tampmem2;
          tmeanmem=tmeanmem2;  
          angcri=angcri2;
          Seq=Seq2;
   end
   %include the extra hadrdening
   phyeq=0.5;
   if abs(phyeq)>0.0
   if Seq>0.0015
       Seq=(1+(0.55+0.45*cos(pi*(s1-1)))*phyeq)*Seq;
       %Seq=(1+0.15*phyeq)*Seq;
   end
   end
   %calculate the fatigue life
   itsum=itsum+1;
   Flife(itsum)=10^((log10(Seq)-log10(Faxial))/(Baxial));
   Tenlimit(itsum)=Faxial*(Flife(itsum))^(Baxial);
   Torlimit(itsum)=Ftor*(Flife(itsum))^(Btor);
   %Epr=matpro(6,1)*Flife(itsum)+matpro(6,2);
   Lifediff=abs((Flife(itsum)-Flife(itsum-1))/Flife(itsum));
                   if itsum>10 || abs(Lifediff)<0.0001
                       itcond=0.0;
                   end
          end  
Fatipre(j,1)=sampmem;
Fatipre(j,2)=tampmem;
Fatipre(j,3)=smeanmem;
Fatipre(j,4)=Seq;
Fatipre(j,5)=angmem*180/pi;
Fatipre(j,6)=alf*180/pi;
Fatipre(j,7)=beta;
Fatipre(j,8)=Epr;
Fatipre(j,9)=Flife(itsum);
Fatipre(j,10)=s1;
Fatipre(j,11)=itsum;
%amonplane(j,5)=abs(angmem-2*pi)*180/pi;
%amonplane(j,6)=abs(angcri-2*pi)*180/pi;
end
  fprintf(fid,'%12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f %12.5f\n',Fatipre');
  fclose(fid);










