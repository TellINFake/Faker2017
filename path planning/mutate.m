function [pop]=mutate(pop)

[s,t]=size(pop);

pop1=pop;

for i=1:2:s
   m=randperm(t-3)+1;
   mutatepoint(1)=min(m(1),m(2));
   mutatepoint(2)=max(m(1),m(2));

   mutate=round((mutatepoint(2)-mutatepoint(1))/2-0.5);
   for j=1:mutate
      zhong=pop(i,mutatepoint(1)+j);
      pop(i,mutatepoint(1)+j)=pop(i,mutatepoint(2)-j);
      pop(i,mutatepoint(2)-j)=zhong;
   end
end

[pop]=chap10_1dis(pop);
for i=1:s
   if pop1(i,t)<pop(i,t)
      pop(i,:)=pop1(i,:);
   end
end