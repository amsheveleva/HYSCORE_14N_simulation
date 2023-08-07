% HYSCORE 
clear, clf, clc
xx=12;
yy=12;
% Experimental
subplot(2,1,1);
[~,~,par3] = eprload('15_HS_VO_ZrTCBPPpd_3417G_10K_9.7127GHz_16-128.DSC');
[par1,par2] = eprload('15_HS_VO_ZrTCBPPpd_3417G_10K_9.7127GHz_16-128_FTT.DSC'); 
[ff,dt] = meshgrid(1000*par1{2},1000*par1{1});

minLevel = 0.1;
maxLevel = 5;

nLevels = 30;            % 14; no. of contours
LevelType = 'log';       % "log" oder "lin"

if LevelType == 'log'; 
   ll=logspace(log10(minLevel),log10(maxLevel),nLevels);
else
   ll=linspace(minLevel,maxLevel,nLevels);
end;

contour(-dt,-ff,real(par2),real(ll*max(max(par2))));
grid;

set(gca,'DataAspectRatio',[1 1 1]); 
set(gca,'XTick',[-24:2:24],'YTick',[0:2:24]);
set(gca,'Linewidth',1,'FontSize',13);

xlabel('\nu_2 [MHz]','FontSize',13);
ylabel('\nu_1 [MHz]','FontSize',13);

xlim([-xx xx]);
ylim([0 yy]);

% xlim([8 18]); 
% ylim([8 18]);

line([0 0],[-24 24],'Color','k','LineStyle','-');
line([-24 24],[0 0],'Color','k','LineStyle','-');
line([-24 24],[-24 24],'Color','k','LineStyle','--');
line([24 -24],[-24 24],'Color','k','LineStyle','--');

hold all;



subplot(2,1,2);
[~,~,par6] = eprload('24_HS_VO_ZrTCBPPpd_3070G_10K_9.7145GHz_16-128.DSC');
[par4,par5] = eprload('24_HS_VO_ZrTCBPPpd_3070G_10K_9.7145GHz_16-128_FFT.DSC'); 
[ff1,dt1] = meshgrid(1000*par4{2},1000*par4{1});


minLevel = 0.02;
maxLevel = 5;

nLevels = 30;            % 14; no. of contours
LevelType = 'log';       % "log" oder "lin"

if LevelType == 'log'; 
   ll1=logspace(log10(minLevel),log10(maxLevel),nLevels);
else
   ll1=linspace(minLevel,maxLevel,nLevels);
end;

contour(-dt1,-ff1,real(par5),real(ll1*max(max(par5))));
grid;

set(gca,'DataAspectRatio',[1 1 1]); 
set(gca,'XTick',[-24:2:24],'YTick',[0:2:24]);
set(gca,'Linewidth',1,'FontSize',13);

xlabel('\nu_2 [MHz]','FontSize',13);
ylabel('\nu_1 [MHz]','FontSize',13);

xlim([-xx xx]);
ylim([0 yy]);

% xlim([8 18]); 
% ylim([8 18]);

line([0 0],[-24 24],'Color','k','LineStyle','-');
line([-24 24],[0 0],'Color','k','LineStyle','-');
line([-24 24],[-24 24],'Color','k','LineStyle','--');
line([24 -24],[-24 24],'Color','k','LineStyle','--');


hold all;




% -------------------------------------------------------------------------
% Dipolar contributions in MHz

% Experimental parameters for simulation

Exp.Sequence = 'HYSCORE';
Exp.Field = (par3.B0VL*1e3)-0.0; 
Exp.mwFreq = par3.MWFQ*(1e-9);
Exp.ExciteWidth = 1/16*1000;
Exp.tau = 0.128;
Exp.t1 = 0.116; 
Exp.t2 = 0.116; 
Exp.dt = 0.016; 
Exp.nPoints = 256;
Exp.Flip = [1 1 2 1];

Exp1.Sequence = 'HYSCORE';
Exp1.Field = par6.B0VL*1e3 -0; % 305 mT
Exp1.mwFreq = par6.MWFQ*1e-9; % X-band
Exp1.ExciteWidth = 1/16*1000;
Exp1.tau = 0.128;
Exp1.t1 = 0.100; 
Exp1.t2 = 0.100; 
Exp1.dt = 0.016; 
Exp1.nPoints = 256;
Exp1.Flip = [1 1 2 1];


% Opt.Method='perturb1';

[Sys1.g, Sys2.g, Sys3.g, Sys3b.g, Sys4.g, Sys5.g]  = deal([1.972 1.972 1.935]);
[Sys1.HStrain, Sys2.HStrain, Sys3.HStrain, Sys3b.HStrain, Sys4.HStrain, Sys5.HStrain] = deal([130 130 230]);

% % Nitrogen
eeQqh1 = 1.9;
eta1 = 0.6;
T=0.33*[1 1 -2];
aiso=-7.5*[1 1 1];
%Q2=[-0.2 -0.75 0.95];
An=T+aiso;

I = 1;         % nuclear spin must be known!
Q1 = eeQqh1/(4*I*(2*I-1)) * [-1+eta1, -1-eta1, 2];
Sys1.Nucs = '51V,14N';
Sys1.Q = [0 0 0;Q1];
Sys1.QFrame = [0 0 0;0 25 0]*pi/180;
Sys1.A = [193 193 530;An]; %total hyperfine matrix
% Sys1.A = N64 + 17*eye(3); %NO2 nitrogen
Sys1.AFrame = [0 0 0;0 0 0]*pi/180;



% I = 1;         % nuclear spin must be known!
% Q2 = eeQqh2/(4*I*(2*I-1)) * [-1+eta2, -1-eta2, 2]
% Sys2.Q = [Q2];
% Sys2.QFrame = [0 15 0]*pi/180;
% Sys2.Nucs = '14N';
% Sys2.weight = 1;
% 
% Sys2.A = -1*[0.0 0.0 18.7]; %total hyperfine matrix
% Sys2.AFrame = [0 0 0]*pi/180;
% 
% 
% % I = 1;         % nuclear spin must be known!
% % Q2 = eeQqh2/(4*I*(2*I-1)) * [-1+eta2, -1-eta2, 2]
% % Sys2.Q = [Q2];
% % Sys2.QFrame = [0 15 0]*pi/180;
% Sys3.Nucs = '14N';
% Sys3.weight = 10;
% 
% Sys3.A = 1*[0.1 0.1 0.1]; %total hyperfine matrix
% Sys3.AFrame = [0 0 0]*pi/180;


%%
% 
subplot(2,1,1);
[~,~,~,p] = saffron({Sys1},Exp);
%[~,~,~,d] = saffron({Sys2},Exp);
%[~,~,~,j] = saffron({Sys3},Exp);
% [x,y,p] = saffron({Sys2, Sys3, Sys3b, Sys4, Sys5},Exp);
[F1,F2] = meshgrid(p.f1,p.f2);
%[L1,L2] = meshgrid(d.f1,d.f2);
%[J1,J2] = meshgrid(j.f1,d.f2);
z1 = p.fd;
%zz1 = d.fd;
%jj = j.fd;
nLevelsOut = 15;            % no. of contours
LevelTypeOut = 'log';       % "log" or "lin"
 
% minLevelOut = 0.04; % for H and F
maxLevelOut = 5; 

minLevelOut = 0.03; % for N

if LevelTypeOut == 'log' 
    ss=logspace(log10(minLevelOut),log10(maxLevelOut),nLevelsOut);
else
   ss=linspace(minLevelOut,maxLevelOut,nLevelsOut);
end;

contour(F1,F2,abs(z1),ss*max(max(abs(z1))),'red','Linewidth',0.6);
%contour(L1,L2,abs(zz1),ss*max(max(abs(zz1))),'green','Linewidth',0.6);
%contour(J1,J2,abs(jj),ss*max(max(abs(jj))),'magenta','Linewidth',0.6);
%contour(F1,F2,abs(z1+zz1+jj),ss*max(max(abs(z1+zz1+jj))),'magenta','Linewidth',0.6);
title (string(Exp.Field) + "mT \tau = " + string(Exp.tau*1000) + "ns");
%title ('NiRad 1.3mM Tol/DCM 5K ' + string(Exp.Field) + "mT \tau = " + string(Exp.tau*1000) + "ns");  
nu_1H = larmorfrq('1H', Exp.Field);
line([0 2*nu_1H],[2*nu_1H 0],'Color','blue','LineStyle','--', 'Linewidth', 2);
  


nu_19F = larmorfrq('14N', Exp.Field);
line([0 2*nu_19F],[2*nu_19F 0],'Color','yellow','LineStyle','--', 'Linewidth', 2);


subplot(2,1,2);
[~,x,y,p2] = saffron(Sys1,Exp1);
%[~,~,~,d2] = saffron(Sys2,Exp1);
%[~,~,~,j2] = saffron(Sys3,Exp1);
[F3,F4] = meshgrid(p2.f1,p2.f2);
%[L3,L4] = meshgrid(d2.f1,p2.f2);
%[J3,J4] = meshgrid(j2.f1,p2.f2);
z2 = p2.fd;
% zz2 = d2.fd;
% jj2 = j2.fd;
nLevelsOut = 10;            % no. of contours
LevelTypeOut = 'log';       % "log" or "lin"
 
% minLevelOut = 0.08; % for H and F
maxLevelOut = 1; 

minLevelOut = 0.1; % for N

if LevelTypeOut == 'log' 
    ss1=logspace(log10(minLevelOut),log10(maxLevelOut),nLevelsOut);
else
   ss1=linspace(minLevelOut,maxLevelOut,nLevelsOut);
end;

%%
contour(F3,F4,abs(z2),ss1*max(max(abs(z2))),'red','Linewidth',0.6);
%contour(L3,L4,abs(zz2),ss1*max(max(abs(zz2))),'green','Linewidth',0.6);
% contour(J3,J4,abs(jj2),ss1*max(max(abs(jj2))),'magenta','Linewidth',0.6);
% contour(F3,F4,abs(z2+zz2+jj2),ss1*max(max(abs(z2+zz2+jj2))),'magenta','Linewidth',0.6);
title (string(Exp1.Field) + "mT \tau = " + string(Exp1.tau*1000) + "ns");
%title ('NiRad 1.3mM Tol/DCM 5K ' + string(Exp1.Field) + "mT \tau = " + string(Exp1.tau*1000) + "ns");

nu_19F = larmorfrq('14N', Exp1.Field);
line([0 2*nu_19F],[2*nu_19F 0],'Color','yellow','LineStyle','--', 'Linewidth', 2);

% subplot(2,2,3);
% [~,~,~,p4] = saffron({Sys1},Exp2);
% [F5,F6] = meshgrid(p4.f1,p4.f2);
% z3 = p4.fd;
% 
% [~,~,~,d4] = saffron({Sys2},Exp2);
% [L5,L6] = meshgrid(d4.f1,d4.f2);
% zz3 = d4.fd;
% 
% [~,~,~,j4] = saffron(Sys3,Exp2);
% [J5,J6] = meshgrid(j4.f1,j4.f2);
% jj4 = j4.fd;
% 
% nLevelsOut = 15;            % no. of contours
% LevelTypeOut = 'log';       % "log" or "lin"
%  
% % minLevelOut = 0.03; % for H and F
% maxLevelOut = 5; 
% 
% minLevelOut = 0.15; % for N
% 
% if LevelTypeOut == 'log' 
%     ss4=logspace(log10(minLevelOut),log10(maxLevelOut),nLevelsOut);
% else
%    ss4=linspace(minLevelOut,maxLevelOut,nLevelsOut);
% end;
% 
% %contour(F5,F6,abs(z3),ss4*max(max(abs(z3))),'red','Linewidth',0.6);
% %contour(L5,L6,abs(zz3),ss4*max(max(abs(zz3))),'green','Linewidth',0.6);
% contour(J5,J6,abs(jj4),ss4*max(max(abs(jj4))),'magenta','Linewidth',0.6);
% contour(F5,F6,abs(z3+zz3+jj4),ss4*max(max(abs(z3+zz3+jj4))),'magenta','Linewidth',0.6);
% title (string(Exp2.Field) + "mT \tau = " + string(Exp2.tau*1000) + "ns");
% %title ('NiRad 1.3mM Tol/DCM 5K ' + string(Exp2.Field) + "mT \tau = " + string(Exp2.tau*1000) + "ns");  
% nu_1H = larmorfrq('1H', Exp2.Field);
% line([0 2*nu_1H],[2*nu_1H 0],'Color','blue','LineStyle','--', 'Linewidth', 2);
  