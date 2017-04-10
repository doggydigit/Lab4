a1 = -20;
a2 = 0;
f1 = 1;
f2 = 4;
Na = 41;
Nf = 81;

a = linspace(a1,a2,Na);
f = linspace(f1,f2,Nf);

s = zeros(Na,Nf);

close all
%figure(1)
%plot([a1,a2],[0,0],[0,0],[f1,f2]);
%xlabel('coupling strength');
%ylabel('ration between intrinsic frequencies');
%hold on
for i=1:Na
    for j=1:Nf
          s(i,j)=ex7bsimul(a(i),f(j));
%          if (ex7bsimul(a(i),f(j)) == 1)
%             plot(a(i),f(j),'o');
%          else
%             plot(a(i),f(j),'x');
%          end
    end
end

imagesc([f1 f2],[a1 a2],s)
%set(gca,'YDir','normal')
grid on
ylabel('Coupling strength');
xlabel('Time constant ratio')