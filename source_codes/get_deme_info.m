function [deme, deme_t, deme_m, deme_r]=get_deme_info(y_max,x_max,ind,deme_count_base,method)

% get deme info
% deme is the number of individuals per deme
% deme_t is the average threshold per deme
% deme_m is the number of migrated individuals per deme
% deme_r is the number of offspring per deme

[~, sorti]=sortrows(ind(:,2:3)); y=ind(sorti,2:3); ind=ind(sorti,:);
index_c = find([true;sum(diff(y),2)~=0]);
values = y(index_c,:); instances = diff(index_c); % instances = nr of individuals with the same xy combination.
instances=[instances; size(y,1)-sum(instances)];
deme_count=[values instances]; clear values instances y

deme_count_base_str=1000*deme_count_base(:,1)+deme_count_base(:,2);
deme_count_str=1000*deme_count(:,1)+deme_count(:,2);
index=ismember(deme_count_base_str,deme_count_str); clear deme_count_base_str deme_count_str 

deme=deme_count_base; deme(index,3)=deme_count(:,3); clear deme_count
deme=reshape(deme(:,3),y_max, x_max);

switch method
    case 'complete'
        deme_t=zeros(y_max,x_max); deme_m=zeros(y_max,x_max); deme_r=zeros(y_max,x_max);
        for dm=1:length(index_c)-1
            row=mean(ind(index_c(dm):index_c(dm+1)-1,3));
            column=mean(ind(index_c(dm):index_c(dm+1)-1,2));
            deme_t(row,column)=mean(ind(index_c(dm):index_c(dm+1)-1,1));    %the individuals sorted by population nr, so averaging over the individuals per population
            deme_m(row,column)=sum(ind(index_c(dm):index_c(dm+1)-1,4));
            deme_r(row,column)=sum(ind(index_c(dm):index_c(dm+1)-1,5));
        end; 
        deme_t(end)=mean(ind(index_c(dm+1):end,1));
        deme_m(end)=sum(ind(index_c(dm+1):end,4)); 
        deme_r(end)=sum(ind(index_c(dm+1):end,5)); clear dm index_c row column
    
    case 'light'        
        deme_t=zeros(y_max,x_max); deme_m=zeros(y_max,x_max); deme_r=zeros(y_max,x_max);
    
    case 'deme_m_only'
        deme_m=zeros(y_max,x_max); deme_t=zeros(y_max,x_max);deme_r=zeros(y_max,x_max);
        for dm=1:length(index_c)-1
            row=mean(ind(index_c(dm):index_c(dm+1)-1,3));
            column=mean(ind(index_c(dm):index_c(dm+1)-1,2));
            deme_m(row,column)=sum(ind(index_c(dm):index_c(dm+1)-1,4));  
        end; 
        deme_m(end)=sum(ind(index_c(dm+1):end,4)); clear dm index_c row column
        
end
        