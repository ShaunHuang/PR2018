% Main Function to Complete Phase Retrieval of Mid-frequency
%% Parameter Initialization
clear
tic

lambda = 632.8e-6;
k=2*pi/lambda;
pixNum=512;
Dim=5;
%% Phase Groundtruth
%Effective Dimension(non-zero)
EffDim=4;

%Amplitude
pupil=makepupil2017(pixNum,Dim/2,EffDim/2);

%Phase
zerVec=lambda*[0, 0, 0, 0.9,...
               0, 0, 0.6, 0.5];%The parameter is corresponding to the parameter of Zemax (Fringe Zernike Phase)
HalfDim = EffDim/2;
NormRadius = HalfDim;
SagTruth=generateZernike(zerVec, pixNum, HalfDim, NormRadius);
PhaseTruth = pupil.*exp(-1j*k*SagTruth);
%% Propagation
f=150;
FocalSag = generateFocal(f,Dim,pixNum);
FocalPhase = exp(1j*k*FocalSag);

ComplexToProp = PhaseTruth.*FocalPhase;
z_prop=[130,140,160];

imgNum = numel(z_prop);
Iout = cell(1,imgNum);
for i = 1:imgNum
    Iout{i}=abs(ASMDiffgpu(ComplexToProp,z_prop(i),lambda,Dim)).^2;
end
PSF=abs(ASMDiffgpu(ComplexToProp,f,lambda,Dim)).^2;

%% Preprocession of Intensity Data
ImgMeasured = Iout;
% ImgMeasured = Iin;
% for i=1:numel(IoutExtended)
%     IoutExtended{i}=padarray(Iout{i},[(pixNum-subImgSizeRow)/2,(pixNum-subImgSizeRow)/2],min(min(Iout{i})),'both');
% end
% ImgMeasured = IoutExtended;
%% Cells Initialization for Phase Retrieval
E_zeros = cell(1,imgNum);
E_ones = cell(1,imgNum);
for ii=1:imgNum
    E_zeros{ii}=zeros(pixNum,pixNum);
    E_ones{ii}=ones(pixNum,pixNum);
end
E_Diff = E_zeros;
E_Diff_sum = E_zeros;
E_Diff_temp = E_zeros;
E_Update = E_zeros;

CostFuncTemp = zeros(1,imgNum);
%% Initial Figure
figNum=2;
figure(figNum)   % create figure for displaying the result through out the iteration
% colormap(gray)
subplot(2,2,1)
imagesc(unwrap(angle(PhaseTruth)))
title('Phase Groundtruth');
axis off
axis equal
%% Initial Solve

    
% Radii=300;
% Amp=lambda/8;
% nxa=30;%ripple number across the aperture
% p=makepupil2017(pixNum,Dim/2,Dim/2);
% [sphereMidReal,rippleReal] = mixSurfaceFunc(Radii,Dim,pixNum,Amp,nxa);
% sphereMidReal=Radii-sqrt(Radii^2-(Dim/2)^2)- sphereMidReal;
% f=Radii/2;

%     %% Initializing to Ideal Lens
%     f_Init=Inf;
%     FocalSag_Init = generateFocal(f_Init,Dim,pixNum);
%     FocalPhase_Init = exp(1j*k*FocalSag_Init);
    %% Initializing to Zernike Surface
    zerVec_Init=lambda*[0, 0, 0, 0.9,...
                   0, 0, 0.6, 0.5];%The parameter is corresponding to the parameter of Zemax (Fringe Zernike Phase)

    Sag_Init=generateZernike(zerVec_Init, pixNum, HalfDim, NormRadius);
    Phase_Init = pupil.*exp(-1j*k*Sag_Init);


% Pupil_Init=makepupil2017(pixNum,Dim/2,EffDim/2);
Pupil_Init=ones(pixNum,pixNum);
E_Obj_Init = Phase_Init.*Pupil_Init;
E_Obj_Result = E_Obj_Init.*FocalPhase;
%% Iteration Settings
iterNum=60;	% The number of iterations to be calculated
CostFunc = ones(1,iterNum);
E_theta_Gradient_Last_Iter=cell(1,imgNum);
for ii=1:imgNum
    E_theta_Gradient_Last_Iter{ii}=ones(pixNum,pixNum)/pixNum/pixNum;%the sum equals to one
end
E_theta_Gradient = E_zeros;
Dq_Last_Iter=E_zeros;
Dq=E_zeros;

p=1; % The counter of the loop
while p < iterNum&&(CostFunc(p)>1e-4)	% The number of iterations that will be
                                    % calculated (This method converge within
                                    % approximately 10-50 iterations).
    p = p+1;      % Increase the counter        

% 2006 Nonlinear
%% Procedure 1
    E_Diff_f = ASMDiffgpu(E_Obj_Result,f,lambda,Dim);
    E_Diff_f(E_Diff_f==0)=eps;
%% Procedure 2 3
    E_Diff{1} = ASMDiffgpu(E_Diff_f,z_prop(1)-f,lambda,Dim);
    E_Diff_sum=E_zeros;
    for imgCount = 1:imgNum        
        % Master plane
        E_Diff{imgCount} = sqrt(ImgMeasured{imgCount}).*(E_Diff{imgCount}./abs(E_Diff{imgCount}));
        % Slaves planes
        for imgCountTemp = 1:imgNum
            if imgCountTemp~=imgCount
                E_Diff_temp{imgCountTemp}=ASMDiffgpu(E_Diff{imgCount},z_prop(imgCountTemp)-z_prop(imgCount),lambda,Dim);
                E_Diff_temp{imgCountTemp} = sqrt(ImgMeasured{imgCountTemp}).*(E_Diff_temp{imgCountTemp}./abs(E_Diff_temp{imgCountTemp}))-E_Diff_temp{imgCountTemp};
                E_Diff_sum{imgCount}=E_Diff_sum{imgCount}+conj(ASMDiffgpu(E_Diff_temp{imgCountTemp},-(z_prop(imgCountTemp)-z_prop(imgCount)),lambda,Dim));
            end
        end
%         E_Diff{imgCount}=E_Diff_sum{imgCount}/(imgNum-1);
        E_Diff{imgCount}=sqrt(ImgMeasured{imgCount}).*(E_Diff{imgCount}./abs(E_Diff{imgCount}));
        
        E_theta_Gradient{imgCount}=2*imag(E_Diff{imgCount}.*E_Diff_sum{imgCount});
        
        Dq{imgCount}=-1/2*E_theta_Gradient{imgCount}+(sum(sum(E_theta_Gradient{imgCount}.^2))/sum(sum(E_theta_Gradient_Last_Iter{imgCount}.^2)))*Dq_Last_Iter{imgCount};
        E_theta_Gradient_Last_Iter{imgCount}=E_theta_Gradient{imgCount};
        Dq_Last_Iter{imgCount}=Dq{imgCount};
%         E_Obj_Result=real(E_Obj_Result)+1j*(imag(E_Obj_Result)-0.00001*Dq);
        

        
        if imgCount<imgNum
%             E_Update{imgCount}=real(E_Diff{imgCount})+1j*(imag(E_Diff{imgCount})-0.001*Dq{imgCount});
            E_Update{imgCount}=abs(E_Diff{imgCount}).*exp(1j*(angle(E_Diff{imgCount})+0.0003*Dq{imgCount}));
%             E_Update{imgCount}=abs(E_Diff{imgCount}).*exp(1j*(angle(E_Diff{imgCount})));
            E_Diff{imgCount+1}=ASMDiffgpu(E_Update{imgCount},z_prop(imgCount+1)-z_prop(imgCount),lambda,Dim);
        else   
            
    % Inverse Diffraction
            E_Obj_Result = ASMDiffgpu(E_Diff{imgCount},-z_prop(imgCount),lambda,Dim);
        
            
        end
    % Statistics
%     CostFuncTemp(imgCount)=norm((abs(E_Diff{imgCount})-sqrt(ImgMeasured{imgCount})),2);
    end 

    %% Synthesize New Object Wavefront
    
    Phase_Obj_Result = angle(E_Obj_Result);
    Amp_Obj_Result = abs(E_Obj_Result);
    Amp_Obj_Result(Amp_Obj_Result<1e-6)=0;
    E_Obj_Result = Amp_Obj_Result.*exp(1j*Phase_Obj_Result);
%     Phase_Obj_Result = angle(E_Obj_Result).*pupil;
%     Amp_Obj_Result = abs(E_Obj_Result).*pupil;
    
    E_Diff_test = ASMDiffgpu(E_Obj_Result,z_prop(1),lambda,Dim);
%     for ii=1:imgNum;CostFuncTemp(ii)=norm((abs(E_Diff{ii})-sqrt(ImgMeasured{ii})),2);end;
    
    CostFuncTemp(1)=norm((abs(E_Diff_test)-sqrt(ImgMeasured{1})),2);
    CostFunc(p)=sum(CostFuncTemp);
    disp([num2str(CostFunc(p)),'(',num2str(p),')']);
    
%% Plot Current Results
    figure(figNum)    
    subplot(2,2,2)								% Plot the phase levels
    imagesc(unwrap((Phase_Obj_Result-angle(FocalPhase)).*Amp_Obj_Result));
    title('Recovered Phase');
    axis equal
    axis off
    
    figure(figNum)
    subplot(2,2,3)								% Plot the phase levels
    imagesc(angle(Phase_Obj_Result./(E_Obj_Result/k)));
    title('MSF Error');
    axis equal
    axis off

    figure(figNum)
    subplot(2,2,4)					% Calculate and plot the intensity quotient as a 
    plot((2:p),CostFunc(2:p),'b.-')	    % function of the number of iterations
    axis([1 iterNum -0.01 max(CostFunc)])
    xlabel('Iteration Times')
    ylabel('Least Square Error')    
    title('Lease Square Error Curve');
    drawnow
    
    
    
end
toc
final=CostFunc(end)

% RMS=sqrt(sum(sum((unwrap((angle(E_Obj_Result)-angle(FocalPhase)).*(abs(E_Obj_Result)>1e-6)))-unwrap(angle(PhaseTruth))).^2))/pixNum^2;
