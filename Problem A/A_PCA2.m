clear;clc;tic;
%%%%%%%%%%%%%%%%%%%%%%%%%% 对训练样本数据预处理 %%%%%%%%%%%%%%%%%%%%%%%%%%%
[t1,t2,t3]=xlsread('附件4：门店交易训练数据');
t=zeros(10000,7);t0=1000;%t0为成交周期惩罚因子
t(:,1)=t1(:,1);t(:,3)=t1(:,3);
%将日期转换成日期数
t(:,2)=datenum(char(t3(:,2)),'yyyy/mm/dd');
t(:,5)=datenum(char(t3(:,5)),'yyyy/mm/dd');
for i=1:10000
    %读取改价信息
    if length(t3{i,4})~=2
       td=double(t3{i,4});k=0;l=0;kl=zeros(2,1);
       for j=1:length(t3{i,4})
        if td(j)==58
           k=k+1;
        end
        if td(length(td)-j+1)==34&&l<2
           l=l+1;kl(l)=length(td)-j;
        end
        if l==2
            %若已售出，则计算售出价格-原定价，作为预期利润增量
            if ~isnan(t3{i,6})
                t(i,7)=str2double(char(td(kl(2)+2:kl(1))))-t(i,3);
            end
            l=l+1;
       end
       t(i,4)=k;%储存改价次数
       end
    end
     if isnan(t3{i,6})
       t(i,6)=t(i,5)+t0;
    else
       t(i,6)=datenum(char(t3(i,6)),'yyyy/mm/dd');
    end
end
%读取处理过的有效训练数据
tt=xlsread('et_data');
[rtt,ctt]=size(tt);

ttt=zeros(10000,31);

for i=1:rtt
for j=1:10000
    if tt(i,1)==tt(j,1)
       ttt(i,:)=[t(i,:),tt(i,2:ctt-1)];
    end
end
end

%下面开始处理数据

x0=zeros(10000,27);
x0(:,1:2)=ttt(:,3:4);
x0(:,3)=ttt(:,7);
x0(:,4:end)=ttt(:,8:end);

y0=ttt(:,6)-ttt(:,2);

[x,y,ljj,lyy]=standardize(x0,y0);
[n,p]=size(x);
R=x'*x;
[v,e]=eig(R);
lambda=diag(e);[lambda,ii]=sort(lambda,'descend');
con=zeros(p,1);sumcon=zeros(p,1);
vv=v;%将v中特征向量重排
for i=1:length(ii)
    v(:,i)=vv(:,ii(i));
end

%计算主成分贡献率
for i=1:p
    con(i)=lambda(i)/p;
    sumcon(i)=sum(lambda(1:i))/p;
end
%前17个主成分贡献率到87.72%>85%
nz=17;
zv=v(:,1:nz);

Z=zeros(27,nz);
for i=1:nz
    Z(:,i)=Z(:,i)+zv(:,i)*lambda(i);%计算得分矩阵
end
%计算各因子的综合得分
Zx=zeros(27,1);
for i=1:27
for j=1:nz
    Zx(i)=Zx(i)+Z(i,j)*con(i);
end
end
%将各变量得分按绝对值排序
Zxrk=zeros(27,2);
[~,Zxrk(:,1)]=sort(abs(Zx),'descend');
Zxrk(:,2)=Zx(Zxrk(:,1));
%%%%%%%%%%%%%%%%%%%%%%%%%% 可视化解释变量 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
subplot(1,2,1);
h11=histogram2(x0(:,1),y0,'FaceColor','flat',...
'XBinLimits',[0 100],'YBinLimits',[0 100],'NumBins',[35 35]);
xlabel('Price');
ylabel('Transaction cycle');
title('Frequency distribution of price and transaction cycle');

subplot(1,2,2);
h12=histogram2(x0(:,2),y0,'FaceColor','flat','YBinLimits',[0 100],'NumBins',[35 35]);
xlabel('Price change times');
ylabel('Transaction cycle');
title('Frequency distribution of price change times and transaction cycle');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(2);
subplot(1,2,1);
h21=histogram2(x0(:,8),y0,'FaceColor','flat','YBinLimits',[0 100],'NumBins',[35 35]);
xlabel('mileage');
ylabel('Transaction cycle');
title('Frequency distribution of mileage and transaction cycle');

subplot(1,2,2);
h22=histogram2(x0(:,5),y0,'FaceColor','flat','YBinLimits',[0 100],'NumBins',[35 35]);
xlabel('brand');
ylabel('Transaction cycle');
title('Frequency distribution of brand and transaction cycle');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3);
subplot(1,2,1);
h31=histogram2(x0(:,14),y0,'FaceColor','flat','YBinLimits',[0 100],'NumBins',[35 35]);
xlabel('registerDate');
ylabel('Transaction cycle');
title('Frequency distribution of registerDate and transaction cycle');

subplot(1,2,2);
h32=histogram2(x0(:,15),y0,'FaceColor','flat','YBinLimits',[0 100],'NumBins',[35 35]);
xlabel('licenseDate');
ylabel('Transaction cycle');
title('Frequency distribution of licenseDate and transaction cycle');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4);
subplot(1,2,1);
h41=histogram2(x0(:,16),y0,'FaceColor','flat','YBinLimits',[0 100],'NumBins',[35 35]);
xlabel('modelyear');
ylabel('Transaction cycle');
title('Frequency distribution of modelyear and transaction cycle');

subplot(1,2,2);
h42=histogram2(x0(:,6),y0,'FaceColor','flat','YBinLimits',[0 100],'NumBins',[35 35]);
xlabel('serial');
ylabel('Transaction cycle');
title('Frequency distribution of serial and transaction cycle');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(5);
subplot(1,2,1);
h51=histogram2(x0(:,7),y0,'FaceColor','flat','YBinLimits',[0 100],'NumBins',[35 35]);
xlabel('model');
ylabel('Transaction cycle');
title('Frequency distribution of model and transaction cycle');

subplot(1,2,2);
h52=histogram2(x0(:,12),y0,'FaceColor','flat','YBinLimits',[0 100],'NumBins',[35 35]);
xlabel('transferCount');
ylabel('Transaction cycle');
title('Frequency distribution of transferCount and transaction cycle');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 部分所用函数 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%标准化函数
function[x_cs,y_cs,ljj,lyy]=standardize(x,y)
[n,q]=size(x);
x_cs=x-ones(n,1)*mean(x);
y_cs=y-ones(n,1)*mean(y);
ljj=zeros(q,1);lyy=y_cs'*y_cs;
y_cs=y_cs/sqrt(lyy);
for i=1:q
    ljj(i)=x_cs(:,i)'*x_cs(:,i);
    x_cs(:,i)=x_cs(:,i)/sqrt(ljj(i));
end
end
