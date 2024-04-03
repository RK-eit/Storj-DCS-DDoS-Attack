format long

C = 360000; % C = c1.x1 + c2.x2
c1 = 50; %c1 = 10,50,100
c2 = 1;
Nv = 25300;
Nu = 50000;
kv = 74; %kv = 40
ku = 6; % ku = 40
kv1 = 75;
ku1 = 5;
kv2 = 76;
ku2 = 4;
kv3 = 77;
ku3 = 3;
kv4 = 78;
ku4 = 2;
kv5 = 79;
ku5 = 1;
kv6 = 80;
ku6 = 0;
k = 80;
d = 29;  % d = 40
M = 2.375*10^9; % Total number of segments in the system 
maxX1= ceil(C/c1);
S2_r = zeros(1,maxX1);
S4_r = zeros(1,maxX1);
S6_r = zeros(1,maxX1);
S8_r = zeros(1,maxX1);
S10_r = zeros(1,maxX1);
S12_r = zeros(1,maxX1);
S14_r = zeros(1,maxX1);
PrF_r = zeros(1,maxX1);
PrF1_r = zeros(1,maxX1);
PrF2_r = zeros(1,maxX1);
PrF3_r = zeros(1,maxX1);
PrF4_r = zeros(1,maxX1);
PrF5_r = zeros(1,maxX1);
PrF6_r = zeros(1,maxX1);


% defining x1 from 1 to ceil C/c1

for X1 = 6200:maxX1
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

for X1 = 6200:maxX1
    X2 = C-(c1.*X1)/(c2);

S4=0; %define S2
    for i= k-d+1:k 
        S3 = 0;
        for j = max(i-ku1,0):min(i,kv1) 
            h1 = hygepdf(j,Nv,X1,kv1); 
            h2 = hygepdf(i-j,Nu,X2,ku1);     
            P = h1.*h2; 
            S3 = S3+P;
        end
        S4 = S4+S3;
        S4(isnan(S4))=0;
        S4_r(X1)=S4;
        PrF1 = 1-(1-S4)^M;
        PrF1_r(X1) = PrF1;
    end
disp([X1, X2, X1+X2, S4,PrF1]);
end

for X1 = 6200:maxX1
    X2 = C-(c1.*X1)/(c2);

S6=0; %define S2
    for i= k-d+1:k 
        S5 = 0;
        for j = max(i-ku2,0):min(i,kv2) 
            h1 = hygepdf(j,Nv,X1,kv2); 
            h2 = hygepdf(i-j,Nu,X2,ku2);     
            P = h1.*h2; 
            S5 = S5+P;
        end
        S6 = S6+S5;
        S6(isnan(S6))=0;
        S6_r(X1)=S6;
        PrF2 = 1-(1-S6)^M;
        PrF2_r(X1) = PrF2;
    end
disp([X1, X2, X1+X2, S6,PrF2]);
end

for X1 = 6200:maxX1
    X2 = C-(c1.*X1)/(c2);

S8=0; %define S2
    for i= k-d+1:k 
        S7 = 0;
        for j = max(i-ku3,0):min(i,kv3) 
            h1 = hygepdf(j,Nv,X1,kv3); 
            h2 = hygepdf(i-j,Nu,X2,ku3);     
            P = h1.*h2; 
            S7 = S7+P;
        end
        S8 = S8+S7;
        S8(isnan(S8))=0;
        S8_r(X1)=S8;
        PrF3 = 1-(1-S8)^M;
        PrF3_r(X1) = PrF3;
    end
disp([X1, X2, X1+X2, S8,PrF3]);
end

for X1 = 6200:maxX1
    X2 = C-(c1.*X1)/(c2);

S10=0; %define S2
    for i= k-d+1:k 
        S9 = 0;
        for j = max(i-ku4,0):min(i,kv4) 
            h1 = hygepdf(j,Nv,X1,kv4); 
            h2 = hygepdf(i-j,Nu,X2,ku4);     
            P = h1.*h2; 
            S9 = S9+P;
        end
        S10 = S10+S9;
        S10(isnan(S10))=0;
        S10_r(X1)=S10;
        PrF4 = 1-(1-S10)^M;
        PrF4_r(X1) = PrF4;
    end
disp([X1, X2, X1+X2, S10,PrF4]);
end

for X1 = 6200:maxX1
    X2 = C-(c1.*X1)/(c2);

S12=0; %define S2
    for i= k-d+1:k 
        S11 = 0;
        for j = max(i-ku5,0):min(i,kv5) 
            h1 = hygepdf(j,Nv,X1,kv5); 
            h2 = hygepdf(i-j,Nu,X2,ku5);     
            P = h1.*h2; 
            S11 = S11+P;
        end
        S12 = S12+S11;
        S12(isnan(S12))=0;
        S12_r(X1)=S12;
        PrF5 = 1-(1-S12)^M;
        PrF5_r(X1) = PrF5;
    end
disp([X1, X2, X1+X2, S12,PrF5]);
end

for X1 = 6200:maxX1
    X2 = C-(c1.*X1)/(c2);

S14=0; %define S2
    for i= k-d+1:k 
        S13 = 0;
        for j = max(i-ku6,0):min(i,kv6) 
            h1 = hygepdf(j,Nv,X1,kv6); 
            h2 = hygepdf(i-j,Nu,X2,ku6);     
            P = h1.*h2; 
            S13 = S13+P;
        end
        S14 = S14+S13;
        S14(isnan(S14))=0;
        S14_r(X1)=S14;
        PrF6 = 1-(1-S14)^M;
        PrF6_r(X1) = PrF6;
    end
disp([X1, X2, X1+X2, S14,PrF6]);
end


plot(1:maxX1,PrF_r,'b', 1:maxX1,PrF1_r,'r',1:maxX1,PrF2_r,'g',1:maxX1,PrF3_r,'c',1:maxX1,PrF4_r,'y',1:maxX1,PrF5_r,'m',1:maxX1,PrF6_r,'k');
xlim([6200 maxX1])
ylim([0 0.028])