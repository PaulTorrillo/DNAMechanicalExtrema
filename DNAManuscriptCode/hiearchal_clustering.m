%Loading in DNA parameter files

load Elasticparams.mat
load Crystalparams.mat
load Molecularparams.mat

%Need to compare all branches in the dendrogram eventually

comparisonmatrix=zeros(16,16);
%Nonrandom sequences given here
provided_seqs={'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA',
    'TTTAAAGGGGTTTAAAGGGGGTTTAAAGGGGTTTAAAGGGGGTTTAAAGGGGTTTAAAGGGGGTTTAAAGGGGTTTAAAGGGGGTTTAAAGGGGTTTAAA',
    'CATGGGAAAAACATGGGAAAAACATGGGAAAAACATGGGAAAAACATGGGAAAAACATGGGAAAAACATGGGAAAAACATGGGAAAAACATGGGAAAAAC',
    'GCCATAGCAAAGGTGATTAATAAATACACAATCTTCTTAGATGACAGGTAAGCTCCATCTTTATTTGTAGGACTACTATAAAATAAATACTCCAACATAA',
    'CAAGCAGCTGCCGGGCTGCGTCCAACTCACGTCGTAAGCAGAACGATTGGGCGCGGAGTCTGCGGAACCCGGGCGGTGGCCGAGCGGCTGCTGGGGGAGC',
    'ATCTCATATGGAGGAAACCGTAAAGAATGCAACCTATGGAAAAAAGCGTTCAACCTACTAATGTATACACTACTACTCTGAATTAGTTACAAACGTTGCA',
    'TCCGCTGCGCGAGGTTTGCCGCCTCGGCTGGTCCGGACTACATTGGTCCGACCCCACTCCAGCCAGCCATTGGAATAGTCATGCCTGGGTCAGAACTGCG'};

%Parametrize sequence
for branch=1:16
    A1 = eye(10);
    if branch<9
        A1=A1*theta0rho0new(1:10,:);
    else
        A1=A1*molecdata;
    end
    if (branch==8)||(branch==16)
            updateseq1=1;
            required_samples=100;
    else
            updateseq1=0;
            seq1=provided_seqs{mod(branch,8)};
            required_samples=1;
    end
    for branch_to_compare=1:(branch-1)
        A2 = eye(10);
        if branch_to_compare<9
            A2=A2*theta0rho0new(1:10,:);
        else
            A2=A2*molecdata;
        end
        if (branch_to_compare==8)||(branch_to_compare==16)
            updateseq2=1;
            required_samples=100;
        else
            updateseq2=0;
            seq2=provided_seqs{mod(branch_to_compare,8)};
            required_samples=1;
        end
samples=zeros(required_samples,1);
%If needed, make random sequence
for current_sample=1:required_samples
operators = {'A', 'T', 'G', 'C'};
indexes = randi(length(operators), 1, 100);
if updateseq1==1
    seq1=cell2mat(operators(indexes));
end

stran=size(seq1);
N = stran(1,2)-1;
seqnum = ones(N,1)*1;
reversechecker=zeros(N,1);
%This is checking to ensure some base-steps are parametrized correctly to
%be consistent in the reverse direction
for k=1:N
    id = double(seq1(k))*100+double(seq1(k+1));
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

%Adding in the parameters and using the Cambridge model
Twist = repmat(34,N, 1);
Tilt = repmat(0.1,N, 1);
Rise = repmat(3.4,N, 1);
Slide = zeros(N, 1);
Shift = zeros(N, 1);
Roll = repmat(0.1,N,1);

for i=1:N 
    FGH = MasterMatrixCopy(:,:,seqnum(i));
    FGH(1:3,4:6) = FGH(1:3,4:6)';
    
    thrho = A1(seqnum(i),:)';
    
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

%Do it again for the other sequence
operators = {'A', 'T', 'G', 'C'};
indexes = randi(length(operators), 1, 100);
if updateseq2==1   
    seq2 = cell2mat(operators(indexes));
end
stran=size(seq2);
N = stran(1,2)-1;
seqnum = ones(N,1)*1;
reversechecker=zeros(N,1);
for k=1:N
    id = double(seq2(k))*100+double(seq2(k+1));
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

for i=1:N 
    FGH = MasterMatrixCopy(:,:,seqnum(i));
    FGH(1:3,4:6) = FGH(1:3,4:6)';
    
    thrho = A2(seqnum(i),:)';

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

%Perform Procustes Analysis

Rvec=transl2;
Rvec(:,:,2)=transl;
Nconfig = length(Rvec(1,1,:));
Npts = length(Rvec(1,:,1));

distvec = [];
RvecT = zeros(size(Rvec));
RvecT(:,:,1) = Rvec(:,:,1);
for i = 1:Nconfig
    for j = i+1:Nconfig
        [dist, X] = procrustes(Rvec(:,:,i)',Rvec(:,:,j)','Scaling',false);
        if i == 1
                RvecT(:,:,j) = X';
        end
        distvec = [distvec dist];
    end
   
end
samples(current_sample)=dist;
end
Z = linkage(distvec);
comparisonmatrix(branch,branch_to_compare)=mean(samples);
    end
end
%Plot
A=6+log10(comparisonmatrix);
v = A(tril(true(size(A)), -1))';

tree=linkage(v,'ward');
labels=["Crystal (A) Repeat","Crystal (T_3A_3G_4T_3A_3G_5) Repeat","Crystal (CATG_3A_5) Repeat","Crystal AT Biased (Crystal)","Crystal GC Biased (Crystal)", "Crystal AT Biased (Molecular)", "Crystal GC Biased (Molecular)","Crystal Random","Molecular (A) Repeat","Molecular (T_3A_3G_4T_3A_3G_5) Repeat","Molecular (CATG_3A_5) Repeat","Molecular AT Biased (Crystal)","Molecular GC Biased (Crystal)", "Molecular AT Biased (Molecular)", "Molecular GC Biased (Molecular)","Molecular Random"];
h=dendrogram(tree,'Labels',labels,'Orientation','left');

set(h,'LineWidth',2)
set(gca, 'FontSize', 12); 


xlabel('RMSD (Angstroms^2)', 'FontSize', 14);  


    
    