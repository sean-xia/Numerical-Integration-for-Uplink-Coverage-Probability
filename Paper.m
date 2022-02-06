 %% I1(x)图像
 clear;
 l = 0.25; % 用户分布（泊松）的密度
 t = 10;   % SINR门限
 T = 10^(t/10); % db2mag
 a = 4; % 路径损耗指数
 e = 1;  % FPC 功率控制因子 e\in [0,1]

 D1 =  @(x,u) 2*pi^2*l^2.*exp(-l*pi.*u)./(1+T^(-1).*u.^(-a*e/2).*x.^a);
 I1 = @(x) integral(@(u) D1(x,u), 0 ,x.^2,'ArrayValued',true);
 fplot(I1,[0,20]);

%%  h(r)图像
 h = @(r) r.*exp(-pi*l.*r.^2);
 fplot(h,[0,10],'k','LineWidth',1.5);
 ylim([0,0.5]);
 xlabel('$r$');
 ylabel('$h(r)$');

 
 %% 比较数值积分方法和解析法
 t = linspace(-10,10,7);
 pcN = arrayfun(@(t) UplinkPC_N(t,1,4),t);
 pcA = arrayfun(@(t) UplinkPC_A(t),t);
 pcS = simul(t,1,4);
 save pc1.mat pcN pcA pcS;
 figure;
 h1 = plot(t, pcN,'o',t,pcA,'-',t,pcS,'-.');
 set(h1,'linewidth',2,'markersize',6);
 legend(h1,'数值积分','解析','仿真');
 xlabel('SINR');
 ylabel('覆盖率');
 
 %% 计算不同路径损耗指数下的覆盖率
 t = linspace(-10,10,7);
 pcN1 = arrayfun(@(t) UplinkPC_N(t,1,4),t);
 pcN2 = arrayfun(@(t) UplinkPC_N(t,1,6),t);
 save pcn.mat pcN1 pcN2
 figure;
 plot(t, pcN1,'o-'...
     ,t,pcN2,'d-');
 
 %% 比较不同sir下的计算时间l=1/4
 t = linspace(-10,10,7);
 nn = length(t);
 e = 0.5;
 a = 4;
 time1 = zeros(1,nn);
 time2 = zeros(1,nn);
 for i = 1:nn
     f1 = @() UplinkPC_N(t(i),e,a);
     time1(i) = timeit(f1);
     f2 = @() UplinkPC_N2(t(i),e,a);
     time2(i) = timeit(f2);
 end
save timel1.mat time1 time2
figure;
bar(t,[time1; time2]);
xlabel('SINR(dB)');
ylabel('计算时间（s)');
legend('有限数值积分','无穷数值积分');

 %% 比较不同sir下的计算时间 l=1e-5
 t = linspace(-10,10,7);
 nn = length(t);
 e = 0.5;
 a = 4;
 time3 = zeros(1,nn);
 time4 = zeros(1,nn);
 for i = 1:nn
     f1 = @() UplinkPC_N(t(i),e,a,1e-6);
     time3(i) = timeit(f1);
     f2 = @() UplinkPC_N2(t(i),e,a,1e-6);
     time4(i) = timeit(f2);
 end
save timel2.mat time3 time4
figure;bar(t,[time3; time4]);
xlabel('SINR(dB)');
ylabel('计算时间（s)');
legend('有限数值积分','无穷数值积分');