% Main Function to Complete Phase Retrieval of Mid-frequency
%% Parameter Initialization
clear
tic

lambda = 632.8e-6;
k=2*pi/lambda;
pixNum=512;
Dim=10;
%% Phase Groundtruth
%Effective Dimension(non-zero)
EffDim=10;

%Amplitude
pupil=makepupil2017(pixNum,Dim/2,EffDim/2);

%Phase
zerVec=lambda*[0, 0, 0, 0,...
               0, 0, 0, 1];%The parameter is corresponding to the parameter of Zemax (Fringe Zernike Phase)
HalfDim = EffDim/2;
NormRadius = HalfDim;
SagTruth=generateZernike(zerVec, pixNum, HalfDim, NormRadius);
PhaseTruth = pupil.*exp(-1j*k*SagTruth);
%% Propagation
f=150;
focalSag = generateFocal(f,Dim,pixNum);
focalPhase = exp(1j*k*focalSag);

ComplexToProp = PhaseTruth.*focalPhase;
z_prop=[130,140,160];

imgNum = numel(z_prop);
Iout = cell(1,imgNum);
for i = 1:imgNum
    Iout{i}=abs(ASMDiffgpu(ComplexToProp,z_prop(i),lambda,Dim)).^2;
end


%% Preprocession of Intensity Data
ImgMeasured = Iout;
% ImgMeasured = Iin;
% for i=1:numel(IoutExtended)
%     IoutExtended{i}=padarray(Iout{i},[(pixNum-subImgSizeRow)/2,(pixNum-subImgSizeRow)/2],min(min(Iout{i})),'both');
% end
% ImgMeasured = IoutExtended;
%% Cells Initialization for Phase Retrieval
E_Obj = cell(1,imgNum);
for ii=1:imgNum;E_Obj{ii}=zeros(pixNum,pixNum);end
E_Diff = E_Obj;
E_invDiff = E_Obj;
E_Diff_sum = E_Obj;
E_Diff_temp = E_Obj;
E_zero = E_Obj;
%% Image Energy
% EM=pixNum;
% Energy_ImgMeasured = zeros(1,imgNum);
% for imgCount=1:imgNum
%     Energy_ImgMeasured(imgCount) = sum(sum(Iout{imgCount}(:,:)));
% end
Energy_Diff = zeros(1,imgNum);
CostFuncTemp = zeros(1,imgNum);
%% Initial Figure
figure(1)   % create figure for displaying the result through out the iteration
colormap(gray)
subplot(2,2,1)
imagesc(angle(phaseReal))
title('Phase Groundtruth');
axis off
axis equal
%% Initial Spatial Values
Lx=Dim;
Ly=Dim;
Mx=EM;
My=EM;
dx1=Lx/Mx;
dy1=Ly/My;%sample interval
X1=-Lx/2:dx1:Lx/2-dx1;
Y1=-Ly/2:dy1:Ly/2-dy1;%coordinate vector of 2D
[x1,y1]=meshgrid(X1,Y1);
Fx1=-1/(2*dx1):1/Lx:1/(2*dx1)-(1/Lx);
Fy1=-1/(2*dy1):1/Ly:1/(2*dy1)-(1/Ly);
[fx1,fy1]=meshgrid(Fx1,Fy1);
%% Initial Solve
% Radii=300;
% Amp=lambda/8;
% nxa=30;%ripple number across the aperture
% p=makepupil2017(pixNum,Dim/2,Dim/2);
% [sphereMidReal,rippleReal] = mixSurfaceFunc(Radii,Dim,pixNum,Amp,nxa);
% sphereMidReal=Radii-sqrt(Radii^2-(Dim/2)^2)- sphereMidReal;
% f=Radii/2;

f=150;
E_Obj_Result_0 = generateFocal(f,Dim,pixNum);
% E_Obj_Result_0 = exp(1j*k*(sphereMidReal*2));
[pupil,~,~,x,y]=makepupil2017(pixNum,Dim/2,EffDim/2);
E_Obj_Result_0 = E_Obj_Result_0.*pupil;
E_Obj_Result = E_Obj_Result_0;
%% Iteration Settings
iterNum=200;	% The number of iterations to be calculated
CostFunc = ones(1,iterNum);

p=1; % The counter of the loop
while p < iterNum&&(CostFunc(p)>1e-6)	% The number of iterations that will be
                                    % calculated (This method converge within
                                    % approximately 10-50 iterations).
    p = p+1;      % Increase the counter        

% 2006 Nonlinear
%% Procedure 1
    E_Diff_f = ASMDiffgpu(E_Obj_Result,f,lambda,Dim);
    E_Diff_f(E_Diff_f==0)=eps;
%% Procedure 2 3
    E_Diff{1} = ASMDiffgpu(E_Diff_f,z_prop(1)-f,lambda,Dim);
    E_Diff_sum=E_zero;
    for imgCount = 1:imgNum        
        % Master plane
        E_Diff{imgCount} = sqrt(ImgMeasured{imgCount}).*(E_Diff{imgCount}./abs(E_Diff{imgCount}));
        % Slaves planes
        for imgCountTemp = 1:imgNum
            if imgCountTemp~=imgCount
                E_Diff_temp{imgCountTemp}=ASMDiffgpu(E_Diff{imgCount},z_prop(imgCountTemp)-z_prop(imgCount),lambda,Dim);
                E_Diff_temp{imgCountTemp} = sqrt(ImgMeasured{imgCountTemp}).*(E_Diff_temp{imgCountTemp}./abs(E_Diff_temp{imgCountTemp}));
                E_Diff_sum{imgCount}=E_Diff_sum{imgCount}+ASMDiffgpu(E_Diff_temp{imgCountTemp},-(z_prop(imgCountTemp)-z_prop(imgCount)),lambda,Dim);
            end
        end
        E_Diff{imgCount}=E_Diff_sum{imgCount}/(imgNum-1);
        E_Diff{imgCount}=sqrt(ImgMeasured{imgCount}).*(E_Diff{imgCount}./abs(E_Diff{imgCount}));
%         E_Diff{imgCount}=sqrt(ImgMeasured{imgCount}).*((E_Diff_sum{imgCount}/(imgNum-1))./abs(E_Diff_sum{imgCount}/(imgNum-1)));
        E_theta_Gradient=2*imag(E_Diff{imgCount}.*E_Diff_sum{imgCount});
%         E_Diff{imgCount}=E_Diff{imgCount}-0.02*(E_Diff{imgCount}-E_theta_Gradient);
%         E_Diff{imgCount}=real(E_Diff{imgCount})+1j*(imag(E_Diff{imgCount})+0.01*(imag(E_Diff{imgCount})-E_theta_Gradient));
        if imgCount<imgNum
            E_Diff{imgCount+1}=ASMDiffgpu(E_Diff{imgCount},z_prop(imgCount+1)-z_prop(imgCount),lambda,Dim);
        else   
    % Inverse Diffraction
            E_Obj_Result = ASMDiffgpu(E_Diff{imgCount},-z_prop(imgCount),lambda,Dim);
        end
    % Statistics
%     CostFuncTemp(imgCount)=norm((abs(E_Diff{imgCount})-sqrt(ImgMeasured{imgCount})),2);
    end 

    %% Synthesize New Object Wavefront
    E_Obj_Result = E_Obj_Result.*pupil;
    
    Phase_Obj_Result = angle(E_Obj_Result);
    Amp_Obj_Result = abs(E_Obj_Result);
    E_Diff_test = ASMDiffgpu(E_Obj_Result,z_prop(1),lambda,Dim);
%     for ii=1:imgNum;CostFuncTemp(ii)=norm((abs(E_Diff{ii})-sqrt(ImgMeasured{ii})),2);end;
    CostFuncTemp(1)=norm((abs(E_Diff_test)-sqrt(ImgMeasured{1})),2);
    CostFunc(p)=sum(CostFuncTemp);
    disp([num2str(CostFunc(p)),'(',num2str(p),')']);
    
%% Plot Current Results
    subplot(2,2,2)								% Plot the phase levels
    imagesc(Phase_Obj_Result);
    title('Recovered Phase');
    axis equal
    axis off
    
    subplot(2,2,3)								% Plot the phase levels
%     imagesc(abs(Phase_Obj_Result-angle(E_Obj_Result_0)));
    imagesc(angle(Phase_Obj_Result./(E_Obj_Result_0/k)));
    title('MSF Error');
    axis equal
    axis off

    subplot(2,2,4)					% Calculate and plot the intensity quotient as a 
    plot((2:p),CostFunc(2:p),'b.-')	    % function of the number of iterations
    axis([1 iterNum -0.01 max(CostFunc)])
    xlabel('Iteration Times')
    ylabel('Least Square Error')    
    title('Lease Square Error Curve');
    drawnow
    
    
    
end
toc
phaseDiff=unwrap_phase(Phase_Obj_Result-angle(E_Obj_Result_0))/k;
figure;imagesc(phaseDiff);
% phaseDiff2=unwrap(Phase_Obj_Result-angle(E_Obj_Result_0))/k;
% figure;imagesc(phaseDiff2);
RMS=sqrt(sum(sum((phaseDiff-rippleReal).^2)));