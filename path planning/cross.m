function [pop]=cross(pop)
[s,t]=size(pop);

pop1=pop;

for i=1:2:s
   m=randperm(t-3)+1;
   crosspoint(1)=min(m(1),m(2));
   crosspoint(2)=max(m(1),m(2));

%   middle=pop(i,crosspoint(1)+1:crosspoint(2));
%   pop(i,crosspoint(1)+1:crosspoint(2))=pop(i+1,crosspoint(1)+1:crosspoint(2));
%   pop(i+1,crosspoint(1)+1:crosspoint(2))=middle;
   
   for j=1:crosspoint(1)
      while find(pop(i,crosspoint(1)+1:crosspoint(2))==pop(i,j))
         zhi=find(pop(i,crosspoint(1)+1:crosspoint(2))==pop(i,j));
         y=pop(i+1,crosspoint(1)+zhi);
         pop(i,j)=y;
      end
   end
   
   for j=crosspoint(2)+1:t-1
      while find(pop(i,crosspoint(1)+1:crosspoint(2))==pop(i,j))
         zhi=find(pop(i,crosspoint(1)+1:crosspoint(2))==pop(i,j));
         y=pop(i+1,crosspoint(1)+zhi);
         pop(i,j)=y;
   end
end
end
[pop]=chap10_1dis(pop);

for i=1:s
    if pop1(i,t)<pop(i,t)
       pop(i,:)=pop1(i,:);
    end
end