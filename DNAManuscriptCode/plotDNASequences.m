clear
%Loading in parameters
load Elasticparams.mat
load Crystalparams.mat
load Molecularparams.mat
rotate3d on
%Put your sequence of choice here
seq = 'TTTAAAGGGGGTTTCAAGGAGGTTTAAAGGGGGTTTGAAGGGGTTTTCAAAAGGGGGTTTGAAAGGGGTTTGAAAGGGGTTTGAAAGGGGTTTGAAAGGGGTTTGAAAGGGGTTTGAAAGGGGTTCGAAGGGGGTTTAAAGGGGTTTAAAGGGGGGTTTGAAAGGGGTTTGAAAGGGGGTTTGAAGGGGGTTTGAAGGGGGTTTGAAGGGGTTTGAAAGGGGTTTGAAAGGGGTTTGAAGGGGGTTTAAAGGGGTTTAAAGGGGGTTTGAAAGGGGTTTGAAAGGGGGTTTGAAGGGGGTTTGAAAGGGGTTTGAAAGGGGTTTGAAAGGGGTTTGAAGGGGGTTTGAAAGGGGTTTGAAGGGGGTTTGAAAGGGGTTTGAAAGGGGTTTGAAAGGGGTTTGAAAGGGGTTTCAGAAGGGTTTAAAGAGGGTTCCAAAGAGGATTTGAAAGGGGTTTCAAGGGGTTCAAAAAAGATTTCAGAAAAGTTTAAAGGGGGTTTCAAGGAGGTTTAAAGGGGGTTTGAAGGGGTTT';
stran=size(seq);
N = stran(1,2)-1;
seqnum = ones(N,1)*1;
reversechecker=zeros(N,1);
grid on;
%Need to get map base-pair steps to appropriate indicies for identification
for k=1:N
    id = double(seq(k))*100+double(seq(k+1));
    if id==6771
        seqnum(k,1)=1;
    elseif (id==6765)||(id==8471)
        seqnum(k,1)=2;
        if id==8471
            reversechecker(k,1)=1;
        end
    elseif id==8465
        seqnum(k,1)=3;
    elseif (id==6571)||(id==6784)
        seqnum(k,1)=4;
        if id==6784
            reversechecker(k,1)=1;
        end
    elseif (id==7171)||(id==6767)
        seqnum(k,1)=5;
        if id==6767
            reversechecker(k,1)=1;
        end
    elseif (id==6565)||(id==8484)
        seqnum(k,1)=6;
        if id==8484
            reversechecker(k,1)=1;
        end
    elseif (id==7165)||(id==8467)
        seqnum(k,1)=7;
         if id==8467
            reversechecker(k,1)=1;
        end
    elseif id==6584
        seqnum(k,1)=8;
    elseif (id==6567)||(id==7184)
        seqnum(k,1)=9;
         if id==7184
            reversechecker(k,1)=1;
        end
    elseif id==7167
        seqnum(k,1)=10;
    end
end
Twist = repmat(34,N, 1);
Tilt = repmat(0.1,N, 1);
Rise = repmat(3.4,N, 1);
Slide = zeros(N, 1);
Shift = zeros(N, 1);
Roll = repmat(0.1,N,1);
A=theta0rho0new(1:10,:);
for i=1:N 
    FGH = MasterMatrixCopy(:,:,seqnum(i));
    FGH(1:3,4:6) = FGH(1:3,4:6)';
    thrho = A(seqnum(i),:)';
   
    if abs(imag(thrho)) > 1e-6
        keyboard
    end
    
    Tilt(i) = thrho(1);
    Roll(i) = thrho(2);
    Twist(i) = thrho(3);
    Shift(i) = thrho(4);
    Slide(i) = thrho(5);
    Rise(i) = thrho(6);
    if reversechecker(i)==1
        Tilt(i)=-thrho(1);
        Shift(i)=-thrho(4);
    end
end

angles = atand(Tilt./Roll);
rolltilt= sqrt(Tilt.^2+Roll.^2);
orient(:,:,1) = eye(3);
for x = 1:N
    orient(:,:,x+1) = orient(:,:,x)*(rotz(Twist(x)/2-angles(x))*roty(rolltilt(x))*rotz(Twist(x)/2+angles(x)));
end
for y = 1:N
    midorient(:,:,y) = orient(:,:,y)*(rotz(Twist(y)/2-angles(y))*roty(rolltilt(y)/2)*rotz(angles(y)));
end
transl=zeros(3,N);
for r = 1:N
    w=midorient(:,:,r);
    transl(:,r+1)=transl(:,r)+Shift(r,1)*w(:,1)+Slide(r,1)*w(:,2)+Rise(r,1)*w(:,3);
end

t=tiledlayout(1,2);
ax1=nexttile;
grid on;
plot3(ax1, transl(1,:),transl(2,:),transl(3,:))%,'-s','LineWidth', 3)
w=9;
l=18;
thickness=2;

for j = 1:N+1
    xv(:,1) = transl(:,j)+w/2*orient(:,1,j)+thickness/2*orient(:,3,j);
    xv(:,2) = transl(:,j)+w/2*orient(:,1,j)-thickness/2*orient(:,3,j);
    xv(:,3) = transl(:,j)+w/2*orient(:,1,j)+l/2*orient(:,2,j)+thickness/2*orient(:,3,j);
    xv(:,4) = transl(:,j)+w/2*orient(:,1,j)+l/2*orient(:,2,j)-thickness/2*orient(:,3,j);
    xv(:,5) = transl(:,j)-w/2*orient(:,1,j)+l/2*orient(:,2,j)+thickness/2*orient(:,3,j);
    xv(:,6) = transl(:,j)-w/2*orient(:,1,j)+l/2*orient(:,2,j)-thickness/2*orient(:,3,j);
    xv(:,7) = transl(:,j)-w/2*orient(:,1,j)+thickness/2*orient(:,3,j);
    xv(:,8) = transl(:,j)-w/2*orient(:,1,j)-thickness/2*orient(:,3,j);
    
    if seq(j)=='A'
        patch(xv(1,[1 3 5 7]),xv(2,[1 3 5 7]),xv(3,[1 3 5 7]),'b')
        patch(xv(1,[2 4 6 8]),xv(2,[2 4 6 8]),xv(3,[2 4 6 8]),'b')
        patch(xv(1,[1 2 4 3]),xv(2,[1 2 4 3]),xv(3,[1 2 4 3]),'b')
        patch(xv(1,[5 6 8 7]),xv(2,[5 6 8 7]),xv(3,[5 6 8 7]),'b')
        patch(xv(1,[5 6 4 3]),xv(2,[5 6 4 3]),xv(3,[5 6 4 3]),'b')
        patch(xv(1,[1 2 8 7]),xv(2,[1 2 8 7]),xv(3,[1 2 8 7]),'b')
    elseif seq(j)=='T'
        patch(xv(1,[1 3 5 7]),xv(2,[1 3 5 7]),xv(3,[1 3 5 7]),'y')
        patch(xv(1,[2 4 6 8]),xv(2,[2 4 6 8]),xv(3,[2 4 6 8]),'y')
        patch(xv(1,[1 2 4 3]),xv(2,[1 2 4 3]),xv(3,[1 2 4 3]),'y')
        patch(xv(1,[5 6 8 7]),xv(2,[5 6 8 7]),xv(3,[5 6 8 7]),'y')
        patch(xv(1,[5 6 4 3]),xv(2,[5 6 4 3]),xv(3,[5 6 4 3]),'y')
        patch(xv(1,[1 2 8 7]),xv(2,[1 2 8 7]),xv(3,[1 2 8 7]),'y')
    elseif seq(j)=='G'
        patch(xv(1,[1 3 5 7]),xv(2,[1 3 5 7]),xv(3,[1 3 5 7]),'g')
        patch(xv(1,[2 4 6 8]),xv(2,[2 4 6 8]),xv(3,[2 4 6 8]),'g')
        patch(xv(1,[1 2 4 3]),xv(2,[1 2 4 3]),xv(3,[1 2 4 3]),'g')
        patch(xv(1,[5 6 8 7]),xv(2,[5 6 8 7]),xv(3,[5 6 8 7]),'g')
        patch(xv(1,[5 6 4 3]),xv(2,[5 6 4 3]),xv(3,[5 6 4 3]),'g')
        patch(xv(1,[1 2 8 7]),xv(2,[1 2 8 7]),xv(3,[1 2 8 7]),'g')
    elseif seq(j)=='C'
        patch(xv(1,[1 3 5 7]),xv(2,[1 3 5 7]),xv(3,[1 3 5 7]),'r')
        patch(xv(1,[2 4 6 8]),xv(2,[2 4 6 8]),xv(3,[2 4 6 8]),'r')
        patch(xv(1,[1 2 4 3]),xv(2,[1 2 4 3]),xv(3,[1 2 4 3]),'r')
        patch(xv(1,[5 6 8 7]),xv(2,[5 6 8 7]),xv(3,[5 6 8 7]),'r')
        patch(xv(1,[5 6 4 3]),xv(2,[5 6 4 3]),xv(3,[5 6 4 3]),'r')
        patch(xv(1,[1 2 8 7]),xv(2,[1 2 8 7]),xv(3,[1 2 8 7]),'r')
    end
end
for j = 1:N+1
    xv(:,1) = transl(:,j)+w/2*orient(:,1,j)-l/2*orient(:,2,j)+thickness/2*orient(:,3,j);
    xv(:,2) = transl(:,j)+w/2*orient(:,1,j)-l/2*orient(:,2,j)-thickness/2*orient(:,3,j);
    xv(:,3) = transl(:,j)+w/2*orient(:,1,j)+thickness/2*orient(:,3,j);
    xv(:,4) = transl(:,j)+w/2*orient(:,1,j)-thickness/2*orient(:,3,j);
    xv(:,5) = transl(:,j)-w/2*orient(:,1,j)+thickness/2*orient(:,3,j);
    xv(:,6) = transl(:,j)-w/2*orient(:,1,j)-thickness/2*orient(:,3,j);
    xv(:,7) = transl(:,j)-w/2*orient(:,1,j)-l/2*orient(:,2,j)+thickness/2*orient(:,3,j);
    xv(:,8) = transl(:,j)-w/2*orient(:,1,j)-l/2*orient(:,2,j)-thickness/2*orient(:,3,j);
    
    if seq(j)=='T'
        patch(xv(1,[1 3 5 7]),xv(2,[1 3 5 7]),xv(3,[1 3 5 7]),'b')
        patch(xv(1,[2 4 6 8]),xv(2,[2 4 6 8]),xv(3,[2 4 6 8]),'b')
        patch(xv(1,[1 2 4 3]),xv(2,[1 2 4 3]),xv(3,[1 2 4 3]),'b')
        patch(xv(1,[5 6 8 7]),xv(2,[5 6 8 7]),xv(3,[5 6 8 7]),'b')
        patch(xv(1,[5 6 4 3]),xv(2,[5 6 4 3]),xv(3,[5 6 4 3]),'b')
        patch(xv(1,[1 2 8 7]),xv(2,[1 2 8 7]),xv(3,[1 2 8 7]),'b')
    elseif seq(j)=='A'
        patch(xv(1,[1 3 5 7]),xv(2,[1 3 5 7]),xv(3,[1 3 5 7]),'y')
        patch(xv(1,[2 4 6 8]),xv(2,[2 4 6 8]),xv(3,[2 4 6 8]),'y')
        patch(xv(1,[1 2 4 3]),xv(2,[1 2 4 3]),xv(3,[1 2 4 3]),'y')
        patch(xv(1,[5 6 8 7]),xv(2,[5 6 8 7]),xv(3,[5 6 8 7]),'y')
        patch(xv(1,[5 6 4 3]),xv(2,[5 6 4 3]),xv(3,[5 6 4 3]),'y')
        patch(xv(1,[1 2 8 7]),xv(2,[1 2 8 7]),xv(3,[1 2 8 7]),'y')
    elseif seq(j)=='C'
        patch(xv(1,[1 3 5 7]),xv(2,[1 3 5 7]),xv(3,[1 3 5 7]),'g')
        patch(xv(1,[2 4 6 8]),xv(2,[2 4 6 8]),xv(3,[2 4 6 8]),'g')
        patch(xv(1,[1 2 4 3]),xv(2,[1 2 4 3]),xv(3,[1 2 4 3]),'g')
        patch(xv(1,[5 6 8 7]),xv(2,[5 6 8 7]),xv(3,[5 6 8 7]),'g')
        patch(xv(1,[5 6 4 3]),xv(2,[5 6 4 3]),xv(3,[5 6 4 3]),'g')
        patch(xv(1,[1 2 8 7]),xv(2,[1 2 8 7]),xv(3,[1 2 8 7]),'g')
    elseif seq(j)=='G'
        patch(xv(1,[1 3 5 7]),xv(2,[1 3 5 7]),xv(3,[1 3 5 7]),'r')
        patch(xv(1,[2 4 6 8]),xv(2,[2 4 6 8]),xv(3,[2 4 6 8]),'r')
        patch(xv(1,[1 2 4 3]),xv(2,[1 2 4 3]),xv(3,[1 2 4 3]),'r')
        patch(xv(1,[5 6 8 7]),xv(2,[5 6 8 7]),xv(3,[5 6 8 7]),'r')
        patch(xv(1,[5 6 4 3]),xv(2,[5 6 4 3]),xv(3,[5 6 4 3]),'r')
        patch(xv(1,[1 2 8 7]),xv(2,[1 2 8 7]),xv(3,[1 2 8 7]),'r')
    end
end
axis equal
grid on
hold on

title('Crystal', 'FontSize',16)
N = stran(1,2)-1;
reversechecker=zeros(N,1);
for k=1:N
    id = double(seq(k))*100+double(seq(k+1));
    if id==6771
        seqnum(k,1)=1;
    elseif (id==6765)||(id==8471)
        seqnum(k,1)=2;
        if id==8471
            reversechecker(k,1)=1;
        end
    elseif id==8465
        seqnum(k,1)=3;
    elseif (id==6571)||(id==6784)
        seqnum(k,1)=4;
        if id==6784
            reversechecker(k,1)=1;
        end
    elseif (id==7171)||(id==6767)
        seqnum(k,1)=5;
        if id==6767
            reversechecker(k,1)=1;
        end
    elseif (id==6565)||(id==8484)
        seqnum(k,1)=6;
        if id==8484
            reversechecker(k,1)=1;
        end
    elseif (id==7165)||(id==8467)
        seqnum(k,1)=7;
         if id==8467
            reversechecker(k,1)=1;
        end
    elseif id==6584
        seqnum(k,1)=8;
    elseif (id==6567)||(id==7184)
        seqnum(k,1)=9;
         if id==7184
            reversechecker(k,1)=1;
        end
    elseif id==7167
        seqnum(k,1)=10;
    end
end
Twist = repmat(34,N, 1);
Tilt = repmat(0.1,N, 1);
Rise = repmat(3.4,N, 1);
Slide = zeros(N, 1);
Shift = zeros(N, 1);
Roll = repmat(0.1,N,1);
A = eye(10);

addpath('./Utilities');
[Buckle_Propeller_Opening, ...
Shear_Stretch_Stagger, ...
pho_W_rot , pho_W_tr, ...
Tilt_Roll_Twist, ...
Shift_Slide_Rise, ...
pho_C_rot , pho_C_tr ] = vector2shapes(nondim2cur(cgDNAp(seq).groundstate));
Tilt=Tilt_Roll_Twist(:,1);
Roll=Tilt_Roll_Twist(:,2);
Twist=Tilt_Roll_Twist(:,3);
Shift=Shift_Slide_Rise(:,1);
Slide=Shift_Slide_Rise(:,2);
Rise=Shift_Slide_Rise(:,3);
angles = atand(Tilt./Roll);
rolltilt= sqrt(Tilt.^2+Roll.^2);
orient2(:,:,1) = eye(3);
for x = 1:N
    orient2(:,:,x+1) = orient2(:,:,x)*(rotz(Twist(x)/2-angles(x))*roty(rolltilt(x))*rotz(Twist(x)/2+angles(x)));
end

for y = 1:N
    midorient(:,:,y) = orient2(:,:,y)*(rotz(Twist(y)/2-angles(y))*roty(rolltilt(y)/2)*rotz(angles(y)));
end
transl2=zeros(3,N);
for r = 1:N
    w=midorient(:,:,r);
    transl2(:,r+1)=transl2(:,r)+Shift(r,1)*w(:,1)+Slide(r,1)*w(:,2)+Rise(r,1)*w(:,3);
end


Rvec=transl;
Rvec(:,:,2)=transl2;
Nconfig = length(Rvec(1,1,:));
Npts = length(Rvec(1,:,1));
ax2 = nexttile;
distvec = [];
RvecT = zeros(size(Rvec));
RvecT(:,:,1) = Rvec(:,:,1);
clear xv
w=9;
l=18;
thickness=2;
for i = 1:Nconfig
    for k = i+1:Nconfig
        [dist, X, transform] = procrustes(Rvec(:,:,i)',Rvec(:,:,k)','scaling',false);
        transform.T
        if i == 1
                RvecT(:,:,k) = X';
                plot3(ax2, RvecT(1,:,k),RvecT(2,:,k),RvecT(3,:,k))
                for p=1:N+1
    neworient(:,:,p)=transform.T'*orient2(:,:,p);
end
for j = 1:N+1
    xv(:,1) = RvecT(:,j,k)+w/2*neworient(:,1,j)+thickness/2*neworient(:,3,j);
    xv(:,2) = RvecT(:,j,k)+w/2*neworient(:,1,j)-thickness/2*neworient(:,3,j);
    xv(:,3) = RvecT(:,j,k)+w/2*neworient(:,1,j)+l/2*neworient(:,2,j)+thickness/2*neworient(:,3,j);
    xv(:,4) = RvecT(:,j,k)+w/2*neworient(:,1,j)+l/2*neworient(:,2,j)-thickness/2*neworient(:,3,j);
    xv(:,5) = RvecT(:,j,k)-w/2*neworient(:,1,j)+l/2*neworient(:,2,j)+thickness/2*neworient(:,3,j);
    xv(:,6) = RvecT(:,j,k)-w/2*neworient(:,1,j)+l/2*neworient(:,2,j)-thickness/2*neworient(:,3,j);
    xv(:,7) = RvecT(:,j,k)-w/2*neworient(:,1,j)+thickness/2*neworient(:,3,j);
    xv(:,8) = RvecT(:,j,k)-w/2*neworient(:,1,j)-thickness/2*neworient(:,3,j);
    
    if seq(j)=='A'
        patch(xv(1,[1 3 5 7]),xv(2,[1 3 5 7]),xv(3,[1 3 5 7]),'b')
        patch(xv(1,[2 4 6 8]),xv(2,[2 4 6 8]),xv(3,[2 4 6 8]),'b')
        patch(xv(1,[1 2 4 3]),xv(2,[1 2 4 3]),xv(3,[1 2 4 3]),'b')
        patch(xv(1,[5 6 8 7]),xv(2,[5 6 8 7]),xv(3,[5 6 8 7]),'b')
        patch(xv(1,[5 6 4 3]),xv(2,[5 6 4 3]),xv(3,[5 6 4 3]),'b')
        patch(xv(1,[1 2 8 7]),xv(2,[1 2 8 7]),xv(3,[1 2 8 7]),'b')
    elseif seq(j)=='T'
        patch(xv(1,[1 3 5 7]),xv(2,[1 3 5 7]),xv(3,[1 3 5 7]),'y')
        patch(xv(1,[2 4 6 8]),xv(2,[2 4 6 8]),xv(3,[2 4 6 8]),'y')
        patch(xv(1,[1 2 4 3]),xv(2,[1 2 4 3]),xv(3,[1 2 4 3]),'y')
        patch(xv(1,[5 6 8 7]),xv(2,[5 6 8 7]),xv(3,[5 6 8 7]),'y')
        patch(xv(1,[5 6 4 3]),xv(2,[5 6 4 3]),xv(3,[5 6 4 3]),'y')
        patch(xv(1,[1 2 8 7]),xv(2,[1 2 8 7]),xv(3,[1 2 8 7]),'y')
    elseif seq(j)=='G'
        patch(xv(1,[1 3 5 7]),xv(2,[1 3 5 7]),xv(3,[1 3 5 7]),'g')
        patch(xv(1,[2 4 6 8]),xv(2,[2 4 6 8]),xv(3,[2 4 6 8]),'g')
        patch(xv(1,[1 2 4 3]),xv(2,[1 2 4 3]),xv(3,[1 2 4 3]),'g')
        patch(xv(1,[5 6 8 7]),xv(2,[5 6 8 7]),xv(3,[5 6 8 7]),'g')
        patch(xv(1,[5 6 4 3]),xv(2,[5 6 4 3]),xv(3,[5 6 4 3]),'g')
        patch(xv(1,[1 2 8 7]),xv(2,[1 2 8 7]),xv(3,[1 2 8 7]),'g')
    elseif seq(j)=='C'
        patch(xv(1,[1 3 5 7]),xv(2,[1 3 5 7]),xv(3,[1 3 5 7]),'r')
        patch(xv(1,[2 4 6 8]),xv(2,[2 4 6 8]),xv(3,[2 4 6 8]),'r')
        patch(xv(1,[1 2 4 3]),xv(2,[1 2 4 3]),xv(3,[1 2 4 3]),'r')
        patch(xv(1,[5 6 8 7]),xv(2,[5 6 8 7]),xv(3,[5 6 8 7]),'r')
        patch(xv(1,[5 6 4 3]),xv(2,[5 6 4 3]),xv(3,[5 6 4 3]),'r')
        patch(xv(1,[1 2 8 7]),xv(2,[1 2 8 7]),xv(3,[1 2 8 7]),'r')
    end
end
for j = 1:N+1
    xv(:,1) = RvecT(:,j,k)+w/2*neworient(:,1,j)-l/2*neworient(:,2,j)+thickness/2*neworient(:,3,j);
    xv(:,2) = RvecT(:,j,k)+w/2*neworient(:,1,j)-l/2*neworient(:,2,j)-thickness/2*neworient(:,3,j);
    xv(:,3) = RvecT(:,j,k)+w/2*neworient(:,1,j)+thickness/2*neworient(:,3,j);
    xv(:,4) = RvecT(:,j,k)+w/2*neworient(:,1,j)-thickness/2*neworient(:,3,j);
    xv(:,5) = RvecT(:,j,k)-w/2*neworient(:,1,j)+thickness/2*neworient(:,3,j);
    xv(:,6) = RvecT(:,j,k)-w/2*neworient(:,1,j)-thickness/2*neworient(:,3,j);
    xv(:,7) = RvecT(:,j,k)-w/2*neworient(:,1,j)-l/2*neworient(:,2,j)+thickness/2*neworient(:,3,j);
    xv(:,8) = RvecT(:,j,k)-w/2*neworient(:,1,j)-l/2*neworient(:,2,j)-thickness/2*neworient(:,3,j);
    
    if seq(j)=='T'
        patch(xv(1,[1 3 5 7]),xv(2,[1 3 5 7]),xv(3,[1 3 5 7]),'b')
        patch(xv(1,[2 4 6 8]),xv(2,[2 4 6 8]),xv(3,[2 4 6 8]),'b')
        patch(xv(1,[1 2 4 3]),xv(2,[1 2 4 3]),xv(3,[1 2 4 3]),'b')
        patch(xv(1,[5 6 8 7]),xv(2,[5 6 8 7]),xv(3,[5 6 8 7]),'b')
        patch(xv(1,[5 6 4 3]),xv(2,[5 6 4 3]),xv(3,[5 6 4 3]),'b')
        patch(xv(1,[1 2 8 7]),xv(2,[1 2 8 7]),xv(3,[1 2 8 7]),'b')
    elseif seq(j)=='A'
        patch(xv(1,[1 3 5 7]),xv(2,[1 3 5 7]),xv(3,[1 3 5 7]),'y')
        patch(xv(1,[2 4 6 8]),xv(2,[2 4 6 8]),xv(3,[2 4 6 8]),'y')
        patch(xv(1,[1 2 4 3]),xv(2,[1 2 4 3]),xv(3,[1 2 4 3]),'y')
        patch(xv(1,[5 6 8 7]),xv(2,[5 6 8 7]),xv(3,[5 6 8 7]),'y')
        patch(xv(1,[5 6 4 3]),xv(2,[5 6 4 3]),xv(3,[5 6 4 3]),'y')
        patch(xv(1,[1 2 8 7]),xv(2,[1 2 8 7]),xv(3,[1 2 8 7]),'y')
    elseif seq(j)=='C'
        patch(xv(1,[1 3 5 7]),xv(2,[1 3 5 7]),xv(3,[1 3 5 7]),'g')
        patch(xv(1,[2 4 6 8]),xv(2,[2 4 6 8]),xv(3,[2 4 6 8]),'g')
        patch(xv(1,[1 2 4 3]),xv(2,[1 2 4 3]),xv(3,[1 2 4 3]),'g')
        patch(xv(1,[5 6 8 7]),xv(2,[5 6 8 7]),xv(3,[5 6 8 7]),'g')
        patch(xv(1,[5 6 4 3]),xv(2,[5 6 4 3]),xv(3,[5 6 4 3]),'g')
        patch(xv(1,[1 2 8 7]),xv(2,[1 2 8 7]),xv(3,[1 2 8 7]),'g')
    elseif seq(j)=='G'
        patch(xv(1,[1 3 5 7]),xv(2,[1 3 5 7]),xv(3,[1 3 5 7]),'r')
        patch(xv(1,[2 4 6 8]),xv(2,[2 4 6 8]),xv(3,[2 4 6 8]),'r')
        patch(xv(1,[1 2 4 3]),xv(2,[1 2 4 3]),xv(3,[1 2 4 3]),'r')
        patch(xv(1,[5 6 8 7]),xv(2,[5 6 8 7]),xv(3,[5 6 8 7]),'r')
        patch(xv(1,[5 6 4 3]),xv(2,[5 6 4 3]),xv(3,[5 6 4 3]),'r')
        patch(xv(1,[1 2 8 7]),xv(2,[1 2 8 7]),xv(3,[1 2 8 7]),'r')
    end
end
        end
        distvec = [distvec dist];
    end
end

Z = linkage(distvec);
xsmall=min([ax1.XLim(1) ax2.XLim(1)]);
xlarge=max([ax1.XLim(2) ax2.XLim(2)]);
ysmall=min([ax1.YLim(1) ax2.YLim(1)]);
ylarge=max([ax1.YLim(2) ax2.YLim(2)]);
zsmall=min([ax1.ZLim(1) ax2.ZLim(1)]);
zlarge=max([ax1.ZLim(2) ax2.ZLim(2)]);
ax1.XAxis.FontSize = 14;
ax1.YAxis.FontSize = 14;
ax1.ZAxis.FontSize = 14;
ax2.XAxis.FontSize = 14;
ax2.YAxis.FontSize = 14;
ax2.ZAxis.FontSize = 14;
axis equal
title('cgDNA+', 'FontSize', 16)
grid on;
