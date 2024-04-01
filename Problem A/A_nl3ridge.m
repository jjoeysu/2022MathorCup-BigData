clear;clc;tic;
%%%%%%%%%%%%%%%%%%%%%%%%%% 对训练样本数据预处理 %%%%%%%%%%%%%%%%%%%%%%%%%%%
[t1,t2,t3]=xlsread('附件1：估价训练数据');
t1(isnan(t1))=0;
tb=zeros(36,1);ta=zeros(30000,35);
for i=1:30000
    for j=1:36
    ta(i,j)=sum(isnan(t3{i,j}));
    end
end

%删除数据缺失过多的列
b=[];
for i=1:36
    tb(i)=sum(ta(:,i));
    if tb(i)>1000
        b=[b,i];
    end
end
 t3(:,b)=[]; t1(:,b)=[];
%将日期转换成日期数
t33=zeros(30000,26);
t33(:,2)=datenum(t3(:,2),'yyyy/mm/dd');
t33(:,12)=datenum(t3(:,12),'yyyy/mm/dd');
t33(:,13)=datenum(t3(:,13),'yyyy/mm/dd');
for i=1:30000
     k1=0;
     k2=0;
     k3=0;
     for j=1:4
     k1=k1+(char(t3{i,24}(j))-48)*10^(4-j);
     k2=k2+(char(t3{i,24}(5+j))-48)*10^(4-j);
     k3=k3+(char(t3{i,24}(10+j))-48)*10^(4-j);
     end
     t33(i,24)=k1*k2*k3;
     
     if ta(i,31)==0
     if length(t3{i,23})==1
         t33(i,23)=0;
     else
     if length(t3{i,23})==3
     if t2{i,30}==1+2
        t33(i,23)=2;
     else
         t33(i,23)=3;
     end
     else
         t33(i,23)=4;
     end
     end
     end
end

t=t1+t33;

%观察部分有缺失项的指标分布
figure(1);
subplot(1,3,1);
h1=histogram(t(:,9),'BinWidth',0.4);
title('carCode');
subplot(1,3,2);
h2=histogram(t(:,14),'BinWidth',0.4,'BinLimits', [2000,2022]);
title('modelyear');
subplot(1,3,3);
h3=histogram(t(:,23),'BinWidth',0.4);
title('a_{11}');

c=[];d=[];
for i=1:30000
    %索引离群样本
    if t(i,26)>=50
        c=[c,i];
    else
    if t(i,26)<0.8
        d=[d,i];
    end
    end
    %填补缺失项
    if t(i,9)==0
       t(i,9)=1;%国标码用众数1填补
    end
    if t(i,14)==0
        t(i,14)=2017;%年款用众数2017填补
    end
    if t(i,18)==0
       t(i,18)=(t(i-1,18)+t(i+1,18))/2;%新车价格只缺失一个值，用前后项均值填补
    end
    if t(i,23)==0
       t(i,23)=3;%类型变量，用众数3填补
    end
end
cd=[c,d];

%把填补后的有效训练数据导出备用
xlswrite('et_data.xlsx',t);

%删除离群样本
t(cd,:)=[];
%剔除无用变量
t(:,1)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%% 对验证样本数据预处理 %%%%%%%%%%%%%%%%%%%%%%%%%%%
[s1,s2,s3]=xlsread('附件2：估价验证数据');
s1(isnan(s1))=0;s1=[s1,zeros(5000,1)];

%删除训练集已删除的列
 s3(:,b)=[]; s1(:,b)=[];
 
 sb=zeros(length(s3(1,:)),1);sa=zeros(5000,length(s3(1,:)));
 for i=1:5000
    for j=1:length(s3(1,:))
    sa(i,j)=sum(isnan(s3{i,j}));
    end
 end

 %删除数据缺失过多的列
bb=[];
for i=1:length(s3(1,:))
    sb(i)=sum(sa(:,i));
    if sb(i)>500
        bb=[bb,i];
    end
end
 s3(:,bb)=[]; s1(:,bb)=[];%发现没有缺失过多的列
 
%将日期转换成日期数，并把部分变量量化
s33=zeros(5000,25);
s33(:,2)=datenum(s3(:,2),'yyyy/mm/dd');
s33(:,12)=datenum(s3(:,12),'yyyy/mm/dd');
s33(:,13)=datenum(s3(:,13),'yyyy/mm/dd');
for i=1:5000
     if ~isnan(s3{i,24})
     k1=0;
     k2=0;
     k3=0;
     for j=1:4
     k1=k1+(char(s3{i,24}(j))-48)*10^(4-j);
     k2=k2+(char(s3{i,24}(5+j))-48)*10^(4-j);
     k3=k3+(char(s3{i,24}(10+j))-48)*10^(4-j);
     end
     s33(i,24)=k1*k2*k3;
     end
     if sa(i,23)==0
     if length(s3{i,23})==1
         s33(i,23)=0;
     else
     if length(s3{i,23})==3
     if s2{i,30}==1+2
        s33(i,23)=2;
     else
         s33(i,23)=3;
     end
     else
         s33(i,23)=4;
     end
     end
     end
end

s=s1+s33;

%观察部分有缺失项的指标分布

figure(2);
subplot(1,2,1);
h4=histogram(s(:,14),'BinWidth',0.4,'BinLimits', [2000,2022]);
title('modelyear(test)');
subplot(1,2,2);
h5=histogram(s(:,23),'BinWidth',0.4);
title('a_{11}(test)');


for i=1:5000
    %填补缺失项
    if s(i,14)==0
        s(i,14)=2017;%年款用众数2017填补
    end
    if s(i,23)==0
       s(i,23)=3;%类型变量，用众数3填补
    end
    if s(i,24)==0
       s(i,24)=(s(i-1,24)+s(i+1,24))/2;%此变量只缺失一个值，用前后项均值填补
    end
end

%剔除无用变量
s(:,1)=[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[n,p]=size(t);
m=ceil(0.8*n);%m个训练样本
x0=t(:,1:p-1);y0=t(:,p);%初始数据集

%%%%%%%%%%%%%%%%%%%%%%%% 多元非线性回归（三次） %%%%%%%%%%%%%%%%%%%%%%%%%%%

%根据样本生成设计矩阵
[temp,x1]=makex(x0,n,p);

%标准化样本数据
[x,y,ljj,lyy]=standardize(x1(1:m,2:temp),y0(1:m));

%使用岭估计来估计回归系数，经过调试选定了岭参数k
k=2e-08;
beta_e=(x'*x+k*eye(temp-1))^(-1)*x'*y;
beta_e=beta_e./sqrt(ljj)*sqrt(lyy);
beta_e=[mean(y0(1:m))-mean(x1(1:m,2:temp))*beta_e;beta_e];
y_hat=x1(1:m,:)*beta_e;

%计算训练集拟合指标
sse1=(y_hat-y0(1:m))'*(y_hat-y0(1:m));rxm=rank(x1(1:m,:));
sst1=(y0(1:m)-mean(y0(1:m)))'*(y0(1:m)-mean(y0(1:m)));
r2_1=1-sse1/sst1;
mse1=sse1/(m-rxm);

y_e=x1(m+1:n,:)*beta_e;
%计算测试集拟合指标
sse2=(y_e-y0(m+1:n))'*(y_e-y0(m+1:n));rxn=rank(x1(m+1:n,:));
sst2=(y0(m+1:n)-mean(y0(m+1:n)))'*(y0(m+1:n)-mean(y0(m+1:n)));
r2_2=1-sse2/sst2;
mse2=sse2/(n-m-rxn);

l=0;%计算预测负值数
for i=1:length(y_e)
    if y_e(i)<0
        l=l+1;
    end
end

%检测预测结果得分
[ape,mape,accuracy5,scores]=fscores(y_e,y0(m+1:n));
toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%% 预测验证数据集 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[sn,sp]=size(s);
[~,xs1]=makex(s,5000,sp+1);
y_test=xs1*beta_e;

sl=0;%计算验证集预测负值数
for i=1:length(y_test)
    if y_test(i)<0
        sl=sl+1;
    end
end

%将预测结果导入txt文档
pv(:,1)=s1(:,1);pv(:,2)=y_test;

fileid = fopen('predictive_value.txt','w');
fprintf(fileid,'%u\t%.2f\n',pv');
fclose(fileid);

toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 部分所用函数 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%生成设计矩阵
function[temp,x1]=makex(x0,n,p)
temp=0;
for i=1:p-1
    temp=temp+i*(i+1)/2;
end
temp=temp+p+(p-1)*p/2;

x1=zeros(n,temp);x1(:,1)=ones(n,1);
x1(:,2:p)=x0;
k=p;
for i=1:p-1
for j=i:p-1
    k=k+1;
    x1(:,k)=x0(:,i).*x0(:,j);
end
end
for i=1:p-1
for j=i:p-1
for l=j:p-1
    k=k+1;
    x1(:,k)=x0(:,i).*x0(:,j).*x0(:,l);
end
end
end
end

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

%评价函数
function[ape,mape,accuracy5,scores]=fscores(y_e,y)
n=length(y);
ape=abs(y_e-y)./y;
mape=sum(ape)/n;
k=0;
for i=1:n
    if ape(i)<=0.05
        k=k+1;
    end
accuracy5=k/n;
end
scores=0.2*(1-mape)+0.8*accuracy5;
end

