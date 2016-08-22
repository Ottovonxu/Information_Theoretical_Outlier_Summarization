function [ out_OF ] = OF_cal( input_attribute )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
attri_len=length(input_attribute);
out_OF=zeros(attri_len,1);
tab_attribute=tabulate(input_attribute);
[tablen,tabwid]=size(tab_attribute);
tab_label=zeros(tablen,1);
for i=1:1:tablen
    if tab_attribute(i,2)==1
        tab_label(i,1)=0;
    else tab_label(i,1)=(tab_attribute(i,2)-1)*log2((tab_attribute(i,2)-1))-tab_attribute(i,2)*log2(tab_attribute(i,2));
    end
end

for j=1:1:attri_len
    for k=1:1:tablen
        if input_attribute(j,1)==tab_attribute(k,1);
            out_OF(j,1)=tab_label(k,1);
        end
    end
end
max_label=max( out_OF);
for i=1:1:attri_len
    if out_OF(i,1)==max_label
        out_OF(i,1)=1;
    else out_OF(i,1)=0;
    end
end
end

