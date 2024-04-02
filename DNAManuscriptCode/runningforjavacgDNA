function final_return = runningforjavacgDNA(seq,desiredindex)
load Elasticparams.mat
load Crystalparams.mat
load Molecularparams.mat
%Set up some initial parameters that will be summed
stran=size(seq);
N = stran(1,2)-1;
seqnum = ones(N,1);
beam_projection=0;
contour_length=0;
overall_compliance=0;
bending_compliance=0;
reversechecker=zeros(N,1);
%Map base pair steps to relevant indicies in parametrization matrix
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
A=A*molecdata;
if N>1
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
end
if N<2
    for i=1:N 
    FGH = MasterMatrixCopy(:,:,seqnum(i)); %Flexibility matrix
    FGH(1:3,4:6) = FGH(1:3,4:6)'; 

    thrho = A(seqnum(i),:)';
    %1 Angstrom is 16.85 degrees
    %Calculating associated flexibilities
    FGH(1:3,1:3)=FGH(1:3,1:3)*283.923;
    FGH(1:3,4:6)=FGH(1:3,4:6)*16.85;
    FGH(4:6,1:3)=FGH(4:6,1:3)*16.85;
    SigF=inv(FGH(1:2,1:2));
    SigFGH = inv(FGH); 
    [U, D] =  eig(SigFGH);
    [UF, DF] =  eig(SigF);

    
    overall_compliance=overall_compliance+trace(D);
    bending_compliance=overall_compliance+trace(DF);
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
%Performing parameter calculations
for met=2:stran(2)
    beam_projection=beam_projection+dot((transl(:,met)-transl(:,met-1))/norm(transl(:,met)-transl(:,met-1)),[0;0;1]);
    contour_length = contour_length + norm(transl(:,met)-transl(:,met-1));
end
rotsum=0;
for met=1:floor(N/2)
    rotsum=rotsum+dot(orient(:,3,met),orient(:,3,met+floor(N/2)));
end
rotsum=rotsum/(floor(N/2));
angoffset = acosd((transl(3,N+1)/(norm(transl(:,N+1)))));
if N>2
    persistence_length=real(-(contour_length/2)/log(rotsum));
else
    persistence_length=0;
end
angdev = acosd(orient(3,3,N+1)/(norm(orient(:,3,N+1))));
bend=norm(cross(transl(:,N+1),transl(:,ceil((N+1)/2))))/norm(transl(:,N+1));
extension=norm(transl(:,N+1));
curvature=diff(transpose(transl),2);
curvsum=0;
for a=1:length(curvature)-2
    curvsum=curvsum+norm(curvature(a,:));
end
possiblereturns=[angoffset persistence_length angdev bend extension contour_length beam_projection curvsum overall_compliance bending_compliance];
final_return=possiblereturns(desiredindex);
