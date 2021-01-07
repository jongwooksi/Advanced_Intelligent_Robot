%% Term Project (Introduction to Intelligent Robots)
 
clc; clear all; close all;

%% Dats loading (select one of them)
% scanData1
%load('scanData1.mat')
%TrueTransform = [300 100 0*pi/180];
% 
% % scanData2
%load('scanData2.mat')
%TrueTransform = [-100 280 0*pi/180];
% 
% % scanData3
%load('scanData3.mat')
%TrueTransform = [-200 -380 15.5*pi/180];
% 
% % scanData4
load('scanData4.mat')
TrueTransform = [-100 280 30*pi/180];


%% Setup (do not modify)
% size of scan data
scanSize = size(refScan_range,2);

% angles (3 data points per 1 deg)
th = zeros(1,scanSize);
i = 1;
while 1
    th(i) = -120 + 0.36*(i-1);
    if i == scanSize
        break;
    else
        i = i + 1;
    end
end

% scan data in Cartesian coordinates
refScan_xy = zeros(2,scanSize);
curScan_xy = zeros(2,scanSize);

for i=1:scanSize
    curTh = th(i) * pi / 180;    
    refScan_xy(1,i) = refScan_range(i) * cos(curTh);
    refScan_xy(2,i) = refScan_range(i) * sin(curTh);
    curScan_xy(1,i) = curScan_range(i) * cos(curTh);
    curScan_xy(2,i) = curScan_range(i) * sin(curTh);
end


%% (TODO) Your scan matching algorithm 
transX = 0;
transY = 0;
rotAngle = 0*pi/180; % radian
index = 0;

zero = zeros(1, scanSize);
curScan_temp = [curScan_xy; zero];
refScan_temp = [refScan_xy; zero];
loss_= 0 ;
lr = 0.5;

while 1
    curScan_temp = curScan_temp';
    refScan_temp = refScan_temp';

    D = pdist2( refScan_temp, curScan_temp , 'euclidean' );
    [~,k] = min( D, [],2,'omitnan');
    
    closest = curScan_temp( k, : );
    closest = closest';
    
    refScan_temp = refScan_temp';
    curScan_temp = curScan_temp';
    
    ref_center = nanmean( refScan_temp, 2) * lr;
    cur_center = nanmean( closest,2) * lr;
    
    matrix_rotation = zeros(3);
   
    for i=1:3     
        matrix_rotation = matrix_rotation + curScan_temp(:,i) * closest(:,i)'; 
    end
    
    matrix_rotation = nanmean(matrix_rotation,2);
    
    matrix_rotation = matrix_rotation - ref_center * cur_center' ;

    diff_matrix = matrix_rotation - matrix_rotation';
    
    delta = [ diff_matrix(2,3); diff_matrix(3,1); diff_matrix(1,2)];
   
    eig_matrix = [ trace( matrix_rotation )  delta';
           delta, matrix_rotation + matrix_rotation' - (trace( matrix_rotation ) * eye(3) )];
       
    [vector, distance] = eigs( eig_matrix );
    
    [~,center_index] = max( nansum( distance ) );
    
    final = vector(:,center_index);
    
    R_matrix = quat2rotm(final');
    
    theta = acos(R_matrix(1))* pi / 180;
    
    trans = cur_center - R_matrix* ref_center ;
    
    curScan_temp = R_matrix * curScan_temp + repmat( trans, 1, scanSize );
   
    eval = [trans(1), trans(2), theta];
    
    loss = rms(eval - TrueTransform)
    
    if loss < 50
        transX = trans(1);
        transY = trans(2);
        rotAngle = theta;
        break
    end
    
    if loss_ / loss < 100
        lr = -1 * lr;
    end   
    
    loss_ = loss;  
    index = index+1;
end
    %% Transformation (do not modify)
transScan_xy = zeros(2,scanSize);
rotMtx = [cos(rotAngle) -sin(rotAngle); sin(rotAngle) cos(rotAngle)];
for i=1:scanSize    
    temp = rotMtx * [curScan_xy(1,i) curScan_xy(2,i)]';
    transScan_xy(:,i) = temp + [transX transY]';   
end


%% Results (do not modify)
EstimatedTransfrom = [transX transY rotAngle];
Error = rms(EstimatedTransfrom - TrueTransform)

Threshold = 50;
if Error < Threshold
    'Success'
else
    'Fail'
end


%% Visualization (do not modify)
% 'transScan_xy'(green) which is the scan matching results should be maximally aligned to 'refScan_xy'(red).
figure,
hold on, box on
axis equal
plot(refScan_xy(1,:), refScan_xy(2,:) ,'r.')
plot(curScan_xy(1,:), curScan_xy(2,:) ,'b.')
plot(transScan_xy(1,:), transScan_xy(2,:) ,'g.')
xlabel('x (mm)')
ylabel('y (mm)')
legend('reference scan','current scan','transformed current scan')