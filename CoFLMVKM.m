clear all;close all;clc

load prokaryotic.mat
points = {gene_repert proteome_comp text};
label=truth;
        
points_view=size(points,2);


cluster_n=length(unique(label));
points_n =size(gene_repert,1);


m     = 2;

%% generate initialization of cluster center
rng(232) %7 8 169
for h=1:points_view
    [n(h),d(h)]=size(points{h});
    clust_cen{h} = rand(cluster_n,d(h));
end

%% MEMBERSHIP INITIALIZATION

for h=1:points_view
    [n(h),d(h)]=size(points{h});
    u{h}=rand(cluster_n,n(h));
    
    u_sum{h}=sum(u{h});
    u{h}=u{h}./repmat(u_sum{h},cluster_n,1);
end 
u0=u;



%% VIEW WEIGHT INITIALIZATION
rng(13)
x=rand([1 points_view]);
wv=x/sum(x,2);


%% FEATURE WEIGHT INITIALIZATION
rng(56)
for h=1:points_view
    [n(h),d(h)]=size(points{h});
    wf{h} = rand(1,d(h));
    
    wf_sum{h}=sum(wf{h});
    wf{h}=wf{h}./repmat(wf_sum{h},1);

end


t_max=70;
index_wf_red_all{h}=[];
SSIGMA = 0;


for itr =1:t_max
    
    alpha  = itr/points_n;
    beta  = itr/points_view;


    % update feature weights 
    
    for h=1:points_view
        temp9=0;
        
        for k=1:cluster_n
            W1=bsxfun(@minus,points{h},clust_cen{h}(k,:)).^2;%bsxfun(@minus,points_temp,clust_cen_temp(k,:));
            for i=1:points_n
            W2=bsxfun(@times,u{h}(k,i).^m,W1);
            
            for hh=1:points_view
                if hh~=h 
                    temp9 = temp9+alpha* (u{h}(k,i) - u{hh}(k,i))^m*W1;
                end
            end
            
            W3=sum(W2,1);
            temp10= sum(temp9,1);
            end
        
        end
        
        W4=(W3+temp10);
        W5=sum(W4,2);
        new_wf{h}=bsxfun(@rdivide,W4,W5); 
    end  
    
    
    %% DISCARD WEIGHT

    for h=1:points_view
        
        [n(h),d(h)]=size(points{h});
        thres_reduce=1/sqrt(n(h)*d(h));   
        index_W_red{h}=find(new_wf{h}<thres_reduce);
        
        %adjust feature weight for H-th view
        new_wf{h}(index_W_red{h})=[];
        new_wf{h}=new_wf{h}/sum(new_wf{h});
        
        %adjust points
        new_points{h}=points{h};
        new_points{h}(:,index_W_red{h})=[];
        
        %adjust cluster center
        new_clust_cen{h}=clust_cen{h};
        new_clust_cen{h}(:,index_W_red{h})=[];
        index_wf_red_all{h}=[index_wf_red_all{h} index_W_red{h}];
        
    end
    
    points=new_points;
    wf=new_wf;
    clust_cen=new_clust_cen;


    %% updating view weights
           
    V5=[];
    for h=1:points_view
        temp12=0;
        temp11=0;
        for k=1:cluster_n
                
            for i=1:points_n
                V1=new_wf{h}.^2.*GetDistance(points{h}(i,:),clust_cen{h}(k,:)).^2; %new_wf{h}.^2*GetDistance(points{h}(i,:),clust_cen{h}(k,:)).^2;%bsxfun(@minus,points_temp,clust_cen_temp(k,:));
                temp12=temp12+(u{h}(k,i).^m.*V1);
                for hh=1:points_view
                    if hh~=h
                       temp11 = temp11+alpha* (u{h}(k,i) - u{hh}(k,i))^m*V1;
                    end
                end
            end
            V3=temp12+temp11;
            
        end
        V4=sum(V3);
        V5=[V5 V4];
%         new_wv=V4./V5;
    end
    V6=V5.^(1/(beta-1));
    V7=sum(V6,2);
    new_wv=bsxfun(@rdivide,V6,V7);
    
    
    %% Update cluster centers
    
    
    for h = 1:points_view
        for k = 1:cluster_n 
            for s = 1:size(points{h},2)
            temp1 = 0; 
            temp2 = 0; 
            temp3 = 0; 
            temp4 = 0; 
            for i = 1:points_n 
                temp1 = temp1 + new_wf{h}.^2*u{h}(k,i)^2*points{h}(i,s); 
                temp2 = temp2 + new_wf{h}.^2*u{h}(k,i)^2; 
                for hh = 1:points_view 
                    if hh ~= h 
                       temp3 = temp3 + alpha*new_wf{h}.^2*(u{h}(k,i) - u{hh}(k,i))^2*points{h}(i,s); 
                       temp4 = temp4 + alpha*new_wf{h}.^2*(u{h}(k,i) - u{hh}(k,i))^2; 
                    end 
                end 
            end 
            new_clust_cen{h}(k,s) = (temp1 + temp3)/(temp2 + temp4); 
            end
                    
             
        end
            
    end
    
    
    %% Update memberships
    
    for h=1:points_view
        for i=1:points_n
            for k=1:cluster_n
                temp1=0;
                temp2=0;
                
                for hh=1:points_view
                    if hh ~= h 
                        temp1 = temp1 + alpha*u{hh}(k,i); 
                        temp2 = temp2 + alpha; 
                    end 
                    
                end
                
                temp3 = 0; 
                for l = 1:cluster_n 
                    temp3 = temp3 + (new_wf{h}.^2*(GetDistance(points{h}(i,:),new_clust_cen{h}(k,:)).^2))/(new_wf{h}.^2*(GetDistance(points{h}(i,:),new_clust_cen{h}(l,:)).^2)); 
                end 
                
                temp4 = 0; 
                for l = 1:cluster_n 
                    temp = 0; 
                    for hh = 1:points_view 
                        if hh ~= h 
                            temp = temp + alpha*u{hh}(l,i); 
                        end 
                    end 
                     
                    temp4 = temp4 + temp/(1 + temp2); 
                end 
                new_u{h}(k,i) = temp1/(1 + temp2) + (1 - temp4)/temp3;
                
            end
        end

    end    
    
    
    u=new_u;
    clust_cen=new_clust_cen;
    wf=new_wf;
    wv=new_wv;
%     points=new_points;
    
    
        for h=1:points_view
            cost_temp(h)=sum(sum(((u0{h}-u{h}).^2)))/(points_n*cluster_n);   %.^0.15;   
        end
        cost(itr)=sum(cost_temp);
        fprintf('rate = %d, fitness = %f\n', itr, cost(itr))
        if itr>1,
            if abs(cost(itr)-cost(itr-1))<0.01 %0.01
                break;
            end
        end
    
    
end

    
uu=(u{1}.*wv(1))+(u{2}.*wv(2))+(u{3}.*wv(3));
    
    
clust=[];
for i=1:points_n
    [num idx]=max(uu(:,i));
    clust=[clust;idx];
end    

AR=1-ErrorRate(label,clust,cluster_n)/points_n;
[valid_external(label,clust) nmi(label,clust)]



