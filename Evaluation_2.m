format long

C = 360000; % C = c1.x1 + c2.x2
c1 = 50; %c1 = 10,50,100
c2 = 1;
Nv = 25000;
Nu = 500;
kv = 77; %kv = 40
ku = 3; % ku = 40
k  = 80;
d = 29;
M =2.375*10^9; % Total number of segments in the system 
maxX1= ceil(C/c1);
S2_r = zeros(1,maxX1);
PrF_r = zeros(1,maxX1);

for X1 = 7189:maxX1
    X2 = C-(c1.*X1)/(c2);

S2=0; %define S2
    for i= k-d+1:k 
        S1 = 0;
        for j = max(i-ku,0):min(i,kv) 
            h1 = hygepdf(j,Nv,X1,kv); 
            h2 = hygepdf(i-j,Nu,X2,ku);     
            P = h1.*h2; 
            S1 = S1+P;
        end
        S2 = S2+S1;
        S2(isnan(S2))=0;
        S2_r(X1)=S2;
        PrF = 1-(1-S2)^M;
        PrF_r(X1) = PrF;
    end
disp([X1, X2, X1+X2, S2,PrF]);
end

plot(1:maxX1,PrF_r)
xlim([7189 maxX1])
ylim([0 1])