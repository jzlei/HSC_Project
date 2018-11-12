function hopf_br = Bif_DDEibtool(k,min_bd,max_bd)
addpath(genpath('./dde_biftool'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Define functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


SNPmodel_rhs=@(xx,p)[(-0.1E1).*p(5).*xx(1,1,:)+0.1E-3.*0.271828E1.^((-0.1E1).*p(2).*p( ...
  4).*xx(1,2,:).*(p(6)+xx(1,2,:)).^(-1)).*p(1).*p(7).*((0.1E1+ ...
  0.1166E2.*xx(1,2,:).^0.129E1).^(-1).*xx(2,2,:)+(-0.1E1).* ...
  0.271828E1.^((-0.1E1).*p(3).*p(5)).*(0.1E1+0.1166E2.*xx(1,3,:) ...
  .^0.129E1).^(-1).*xx(2,3,:));0.375E0.*0.271828E1.^((-0.1E0).*p(8)) ...
  .*xx(2,4,:).*(0.625E-1+xx(2,4,:).^4).^(-1)+xx(2,1,:).*((-0.469E-2) ...
  +(-0.1E1).*p(1).*(0.1E1+0.1166E2.*xx(1,1,:).^0.129E1).^(-1)+( ...
  -0.1875E0).*(0.625E-1+xx(2,1,:).^4).^(-1)+(-0.144E0).*(0.36E0+xx( ...
  3,1,:)).^(-1));(-0.24E1).*xx(3,1,:)+0.101714E3.*0.271828E1.^(( ...
  -0.27E0).*p(9).*xx(3,5,:).*(0.797E2+xx(3,5,:)).^(-1)).*xx(2,5,:).* ...
  (0.36E0+xx(3,5,:)).^(-1)];

SPmodel_rhs=@(xx,p)[(-0.1E1).*p(5).*xx(1,1,:)+0.1E-3.*0.271828E1.^((-0.1E1).*p(2).*p( ...
  4).*xx(1,2,:).*(p(6)+xx(1,2,:)).^(-1)).*p(1).*p(7).*((0.1E1+ ...
  0.1166E2.*xx(1,2,:).^0.129E1).^(-1).*xx(2,2,:)+(-0.1E1).* ...
  0.271828E1.^((-0.1E1).*p(3).*p(5)).*(0.1E1+0.1166E2.*xx(1,3,:) ...
  .^0.129E1).^(-1).*xx(2,3,:));xx(2,1,:).*((-0.245247E-1)+(-0.1E1).* ...
  p(1).*(0.1E1+0.1166E2.*xx(1,1,:).^0.129E1).^(-1)+(-0.1875E0).*( ...
  0.625E-1+xx(2,1,:).^4).^(-1))+0.375E0.*0.271828E1.^((-0.1E0).*p(8) ...
  ).*xx(2,4,:).*(0.625E-1+xx(2,4,:).^4).^(-1)];

Pmodel_rhs=@(xx,p)[(-0.1E1).*p(5).*xx(1,1,:)+0.271828E1.^((-0.1E1).*p(2).*p(4).*xx( ...
  1,2,:).*(p(6)+xx(1,2,:)).^(-1)).*p(1).*p(7).*(0.11E-3.*(0.1E1+ ...
  0.1166E2.*xx(1,2,:).^0.129E1).^(-1)+(-0.11E-3).*0.271828E1.^(( ...
  -0.1E1).*p(3).*p(5)).*(0.1E1+0.1166E2.*xx(1,3,:).^0.129E1).^(-1)) ...
  ];

Smodel_rhs=@(xx,p)[xx(1,1,:).*((-0.245247E-1)+(-0.194814E-1).*p(1)+(-0.1875E0).*( ...
  0.625E-1+xx(1,1,:).^4).^(-1))+0.375E0.*0.271828E1.^((-0.1E0).*p(8) ...
  ).*xx(1,2,:).*(0.625E-1+xx(1,2,:).^4).^(-1)];

par0=[0.372E0,0.424E0,0.95E1,7,0.5E-1,0.359E2,487837,0.28E1,0.35E1];

x0=[0.31071E1;0.11E1;0.69E1];







model_functions={Pmodel_rhs,SPmodel_rhs,SNPmodel_rhs,Smodel_rhs};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameter Definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sys_ntau=@(k)[(k<4)*(k+1)+(k==4)*2];
sys_tau=@(ind,xx,par)[(ind==1)*par(4)+(ind==2)*(par(4)+par(3))+...
    (ind==3)*par(8)+(ind==4)*par(9)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% System Structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

funcs=set_funcs(...
    'sys_rhs',model_functions{k},...
    'sys_ntau',sys_ntau(k),...
    'sys_tau',sys_tau,...
    'x_vectorized',1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Steady State
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
par = par0;
par(4) = 0.01;
par(1) = min_bd;
ss.kind='stst';
ss.parameter = par;
ss.x = x0(1:k);
ssx = ss.x;


method=df_mthod(funcs,'stst',1);
method1=method;
[ss,success]=p_correc(funcs,ss,[],[],method.point);
ss.x;
ss.stability = p_stabil(funcs,ss,method.stability);
ss2 = ss;
ss2.parameter(4) = ss.parameter(4) + 0.01;
[ss2, ss_succ] = p_correc(funcs,ss2,[],[],method.point);
ss2.stability = p_stabil(funcs,ss2,method.stability);
ss2.x;
ss_br = df_brnch(funcs,[4],'stst',1);
ss_br.point = ss;
ss_br.point(2) = ss2;
ss_br.parameter.min_bound=[4 0.01];
ss_br.parameter.max_bound=[4 50];
ss_br.parameter.max_step=[4 0.1];
[ss_br, succ,fail,rjct] = br_contn(funcs,ss_br,500);
hopf_br=br_rvers(ss_br);
disp('Rvers 1')
[ss_br, succ,fail,rjct] = br_contn(funcs,ss_br,500);
ss_br=br_stabl(funcs,ss_br,0,0);
[xm,ym]=df_measr(1,ss_br);
figure;
br_plot(ss_br,xm,ym,'b');
ym.subfield='l0';
br_plot(ss_br,xm,ym,'c');
plot([0 5],[0 0],'-.');
axis([0 5 -2 1.5]);
xlabel('a21');ylabel('\Re\lambda');
figure;
br_plot(ss_br,[],ym,'b');
br_plot(ss_br,[],ym,'b.');
plot([0 30],[0 0],'-.');
xlabel('point number along branch');ylabel('\Re(\lambda)');


close all;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Hopf bifurcation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('############## Hopf Bifurcation ##############')
ind1 = find(arrayfun(@(x)real(x.stability.l0(1))>0,ss_br.point));
ind2 = find(arrayfun(@(x)real(x.stability.l0(1))<=0,ss_br.point));
if length(ind1)>0 & length(ind2)>0
    [ind0 indx] = find(arrayfun(@(x)min(abs(x-ind2)),ind1)==1);
    ind = ind1(indx(1));
elseif length(ind2) == 0
    ind = [];
    length(ind1)
    length(ind2)
    hopf_br=df_brnch(funcs,[1 4],'hopf');
    par1 = par;
    par1(4) = 0;
    hopf_br.point.parameter = par1;
    return;
elseif length(ind1) == 0
    ind = [];
    length(ind1)
    length(ind2)
    hopf_br=df_brnch(funcs,[1 4],'hopf');
    par1 = par;
    par1(4) = 1e5;
    hopf_br.point.parameter = par1;
    return;
else
    ind = []
    hopf_br=df_brnch(funcs,[1 4],'hopf');
    return;
end;
hf=p_tohopf(funcs,ss_br.point(ind));
method = df_mthod(funcs,'hopf',1);
method.stability.minimal_real_part=-1;
method.stability.minimal_time_step=1e-9;
[hopf,success1]=p_correc(funcs,hf,[4],[],method.point);




hopf_br=df_brnch(funcs,[1 4],'hopf');
hopf_br.parameter.min_bound(1:2,:)=[[1 min_bd]' [4 0]']'; 
hopf_br.parameter.max_bound(1:2,:)=[[1 max_bd]' [4 100]']'; 
hopf_br.parameter.max_step(1:2,:)=[[1 min(0.01, min_bd/5)]' [4 0.01]']';   
hopf_br.point=hopf;           

hf.parameter(1)=hf.parameter(1)+0.001;               
[hopf,success2]=p_correc(funcs,hf,4,[],method.point); 
hopf_br.point(2)=hopf;                               
if success1==0 | success2 == 0
    success1,success2
    hopf_br.point = [];
    return;
end;

[hopf_br,s,f,r]=br_contn(funcs,hopf_br,8500);
hopf_br = br_rvers(hopf_br);
[hopf_br,s,f,r]=br_contn(funcs,hopf_br,1500);
close all;

