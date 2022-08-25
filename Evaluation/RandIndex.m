function [AR,RI,MI,HI]=RandIndex(c1,c2)
%RANDINDEX - calculates Rand Indices to compare two partitions
% ARI=RANDINDEX(c1,c2), where c1,c2 are vectors listing the 
% class membership, returns the "Hubert & Arabie adjusted Rand index".
% [AR,RI,MI,HI]=RANDINDEX(c1,c2) returns the adjusted Rand index, 
% the unadjusted Rand index, "Mirkin's" index and "Hubert's" index.
%
% See L. Hubert and P. Arabie (1985) "Comparing Partitions" Journal of 
% Classification 2:193-218

%(C) David Corney (2000)   		D.Corney@cs.ucl.ac.uk

if nargin < 2 | min(size(c1)) > 1 | min(size(c2)) > 1
   error('RandIndex: Requires two vector arguments')
   return
end

C=Contingency(c1,c2);	%form contingency matrix

n=sum(sum(C));
nis=sum(sum(C,2).^2);		%sum of squares of sums of rows
njs=sum(sum(C,1).^2);		%sum of squares of sums of columns

t1=nchoosek(n,2);		%total number of pairs of entities
t2=sum(sum(C.^2));	%sum over rows & columnns of nij^2
t3=.5*(nis+njs);

%Expected index (for adjustment)
nc=(n*(n^2+1)-(n+1)*nis-(n+1)*njs+2*(nis*njs)/n)/(2*(n-1));

A=t1+t2-t3;		%no. agreements
D=  -t2+t3;		%no. disagreements

if t1==nc
   AR=0;			%avoid division by zero; if k=1, define Rand = 0
else
   AR=(A-nc)/(t1-nc);		%adjusted Rand - Hubert & Arabie 1985
end

RI=A/t1;			%Rand 1971		%Probability of agreement
MI=D/t1;			%Mirkin 1970	%p(disagreement)
HI=(A-D)/t1;	%Hubert 1977	%p(agree)-p(disagree)

function Cont=Contingency(Mem1,Mem2)

if nargin < 2 | min(size(Mem1)) > 1 | min(size(Mem2)) > 1
   error('Contingency: Requires two vector arguments')
   return
end

Cont=zeros(max(Mem1),max(Mem2));

for i = 1:length(Mem1);
   Cont(Mem1(i),Mem2(i))=Cont(Mem1(i),Mem2(i))+1;
end



%---------------------------------------------------------------

% function adjrand=RandIndex(u,v)
% 
% %function adjrand=adjrand(u,v)
% %
% % Computes the adjusted Rand index to assess the quality of a clustering.
% % Perfectly random clustering returns the minimum score of 0, perfect
% % clustering returns the maximum score of 1.
% %
% %INPUTS
% % u = the labeling as predicted by a clustering algorithm
% % v = the true labeling
% %
% %OUTPUTS
% % adjrand = the adjusted Rand index
% %
% %
% %Author: Tijl De Bie, february 2003.
% 
% 
% 
% n=length(u);
% ku=max(u);
% kv=max(v);
% m=zeros(ku,kv);
% for i=1:n
%     m(u(i),v(i))=m(u(i),v(i))+1;
% end
% mu=sum(m,2);
% mv=sum(m,1);
% 
% a=0;
% for i=1:ku
%     for j=1:kv
%         if m(i,j)>1
%             a=a+nchoosek(m(i,j),2);
%         end
%     end
% end
% 
% b1=0;
% b2=0;
% for i=1:ku
%     if mu(i)>1
%         b1=b1+nchoosek(mu(i),2);
%     end
% end
% for i=1:kv
%     if mv(i)>1
%         b2=b2+nchoosek(mv(i),2);
%     end
% end
% 
% c=nchoosek(n,2);
% 
% adjrand=(a-b1*b2/c)/(0.5*(b1+b2)-b1*b2/c);