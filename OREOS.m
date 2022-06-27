%% Optical REspOnse Simulator (OREOS) 

% Welcome to OREOS, the package that allows simulation of the optical response and 
% vibronic spectra of various molecular systems! If this is your first time using OREOS, 
% please refer to the README file for instructions, samples, and more information. 

% This hub should be navigated section-by-section 
% (use "ctrl+enter" instead of MATLAB's "Run" button)

% Please refer to <J. Phys. Chem. C 2022, 126, 1, 120–131> for more description of the model 

clear
close all

cd(fileparts(matlab.desktop.editor.getActiveFilename))
addpath(genpath(pwd));

%% Choose and calculate molecular Hamiltonian

% remember to set your model parameters in the "modelparamters.m" file!
model = 'monomer'; % As of this distribution, model options are: 'monomer', 'FE dimer'

p = modelparameters(model);
Hmat = RSgenHam(model,p);

%% Set phenomenological bath parameters:

bfluct = 3000;           % [cm-1]      
bfluct = bfluct*p.c;     % [fs-1] 
btime = 40;              % [fs]
bt2fluct = 125;          % [cm-1]
bt2fluct = bt2fluct*p.c; % [fs-1]  
bt2time = 300;           % [fs]

%% Prepare time and frequency vectors

% Res will be the primary structure for all simulated data

Res.dt = 3;           % timestep along coherence and rephasing time delays 
Res.tmax = Res.dt*32; % maximum time for coherence and rephasing time delays
Res.dt2 = 6;          % timestep along waiting time delay 
Res.t2max = 300;        %  waiting time delay 
Res.nsteps = Res.tmax/Res.dt;  
Res.numres = 4;

Res.t1 = (0:Res.dt:Res.tmax);        % coherence time
Res.t2 = (0:Res.dt2:Res.t2max);      % waiting time
Res.t3 = (0:Res.dt:Res.tmax);        % rephasing time
npoints = (length(Res.t2))*(Res.nsteps+1)^2;

% Define frequency w3 axis for Fourier Transform
Res.padnum = 4;                 % zero padding factor
Res.nw1 = (Res.nsteps+1)*Res.padnum; 
Res.nw3 = (Res.nsteps+1)*Res.padnum;  
Res.dw3 = 1 / (p.c*Res.dt*Res.nw3);   % defining the frequency step [cm-1]
Res.dw1 = 1 / (p.c*Res.dt*Res.nw1);   % defining the frequency step [cm-1]

% frequency axes [cm-1%]
Res.w3 = ((-Res.nw3*Res.dw3)/2:Res.dw3:((Res.nw3*Res.dw3)/2-Res.dw3));  
Res.w3=Res.w3+p.e1;     % correct for rotating frame

Res.w1 = ((-Res.nw1*Res.dw1)/2:Res.dw1:((Res.nw1*Res.dw1)/2-Res.dw1));   
Res.w1=Res.w1+p.e1;     % correct for rotating frame

%% Calculate linear response using full Hamiltonian

Res.Rt1 = zeros(size(Res.t1));
ket = Hmat.fock(:,1);               % assumption that all population starts in the global ground state
bra = ket';
U = @(t) expm(-1i*Hmat.H*t*2*pi);   % time propagator for numerical integration

for i = 1:length(Res.t1)
      Res.Rt1(i) = bra*(Hmat.mu)*U(Res.t1(i))*(Hmat.mu)*ket;
end

% phenomenological dephasing
t1ls = exp(-(bfluct^2)*(btime^2)*(exp(-(Res.t1)/btime)+((Res.t1)/btime)-1)); 
Res.Rt1d = Res.Rt1.*t1ls;   % multiply by dephasing

Res.Rw1 = real(fftshift(fft(Res.Rt1d,Res.nw1)));
Res.Rw1norm = (Res.Rw1-min(Res.Rw1))./(max(Res.Rw1)-min(Res.Rw1));
Res.w1linres =fliplr(Res.w1-p.e1)+p.e1;      % correct for rotating frame

f1 = fig_open;
plot(Res.w1linres./1000,Res.Rw1norm,'k','LineWidth',2)
xlim([p.e1/1000-5 p.e1/1000+5])
ylim([0 1.1])
xlabel('\omega/2\pic (10^{3} cm^{-1})');
ylabel('absorption cross section (norm.)');

%% Evaluate phenomenological dephasing

T1mesh = ndgrid(Res.t1,Res.t1)';
T3mesh = ndgrid(Res.t3,Res.t3);
t1t3ls = exp(-(bfluct^2)*(btime^2)*(exp(-(T1mesh+T3mesh)/btime)+((T1mesh+T3mesh)/btime)-1)); 
t1t3ls3D = repmat(t1t3ls,[1,1,length(Res.t2)]);
t2ls = exp(-(bt2fluct^2)*(bt2time^2)*(exp(-(Res.t2)/bt2time)+((Res.t2)/bt2time)-1)); 
t2ls3D = ones(length(Res.t3),length(Res.t1),length(Res.t2)).*permute(t2ls,[1, 3, 2]);

%% Set up nonlinear response functions

for i = 1:Res.numres
    Res.Rt1t3{i} = zeros(Res.nsteps+1,Res.nsteps+1,length(Res.t2));   %pre-allocate response function matrix
end
%R{1} = Non-rephasing ground state bleach
%R{2} = Rephasing stimulated emission
%R{3} = Rephasing ground state bleach
%R{4} = Non-rephasing stimulated emission

Ugg = @(t) expm(-1i*Hmat.Hgg*t*2*pi);  %incremental propagator for ground-state Hamiltonian block
Use = @(t) expm(-1i*Hmat.Hee*t*2*pi);  %incremental propagator for excited-state Hamiltonian block

%assuming no thermal distribution
ket = Hmat.fockgg(:,1);    
bra = ket';

%% Simulate nonlinear molecular response site basis

% this is the bottleneck section - I recommend running a minimal model
% (i.e., less vibrational quanta, less timepoints, etc.) to ensure that
% things are behaving prior to committing to the full simulation

count = 0;
progressbar2 % Initialize progress bar
for k = 1:length(Res.t2)
    for j = 1:length(Res.t1)    % column index is coherence time
        for i= 1:length(Res.t3) % row index is rephasing time
            
            Res.Rt1t3{1}(i,j,k) = bra*Hmat.MUeg'*Use(Res.t3(i))*Hmat.MUge'...
                *Ugg(Res.t2(k))*Hmat.MUeg'*Use(Res.t1(j))*Hmat.MUge'*ket;

            Res.Rt1t3{2}(i,j,k) = bra*Hmat.MUge*conj(Use(Res.t1(j)))...
                *conj(Use(Res.t2(k)))*Hmat.MUeg*conj(Ugg(Res.t3(i)))*(Hmat.MUeg')*Use(Res.t3(i))...
                *Use(Res.t2(k))*Hmat.MUge'*ket; 

            Res.Rt1t3{3}(i,j,k) = bra*Hmat.MUge*conj(Use(Res.t1(j)))*Hmat.MUeg...
                *conj(Ugg(Res.t2(k)))*Hmat.MUeg'*Use(Res.t3(i))*Hmat.MUge'*ket;

            Res.Rt1t3{4}(i,j,k) = bra*Hmat.MUge*conj(Use(Res.t2(k)))*Hmat.MUeg...
                *conj(Ugg(Res.t3(i)))*Hmat.MUeg'*Use(Res.t3(i))*Use(Res.t2(k))...
                *Use(Res.t1(j))*Hmat.MUeg*ket;  
            
            count=count+1;
        end
    end
    progressbar2(k/length(Res.t2)) % Update progress bar
end
progressbar2(1) % close progress bar

%% Apply artificial lineshapes

for i = 1:Res.numres
    Res.Rt1t3d{i} = zeros(size(Res.Rt1t3{i}));   %pre-allocate
end

for i = 1:Res.numres
    Res.Rt1t3d{i} = -Res.Rt1t3{i}.*t1t3ls3D.*t2ls3D;
end

%% Process time-domain signals
% Apply the trapezoidal rule correction (see Hamm and Zanni ch. 9)

Res.Rt1t3c = Res.Rt1t3d;

for i = 1:Res.numres
        Res.Rt1t3c{i}(:,1,:) = Res.Rt1t3c{i}(:,1,:)./2;   
        Res.Rt1t3c{i}(1,2:end,:) = Res.Rt1t3c{i}(1,2:end,:)./2;
end

%Convert to frequency domain

for i = 1:Res.numres
    Res.Rw1w3{i} = zeros(Res.nw1,Res.nw3,length(Res.t2));
end

for i = 1:Res.numres
    Res.Rw1w3{i} = rescalc(Res.Rt1t3c{i},Res,i);
end

%% Option to trim and sub-sample data matrix to a particular (w1,w3) window

maxmat = Res.Rw1w3{1}(:,:,1);
maxmat = abs(maxmat);

% center of window is defined as the point of maximum signal + a shift
% equal to the highest-frequency vibration
[row, col] = find(ismember(maxmat, max(maxmat(:))));

width = 2500;                   % width of window in w1 domain [cm-1]
height = 3500;                  % width of window in w3 domain [cm-1]
w1offset = 0.9*max(p.wfreq);    % offset based on highest frequency vibration [cm-1]
w3offset = -0*max(p.wfreq);     % offset based on highest frequency vibration [cm-1]
indw1 = Res.w1>Res.w1(col)-width+w1offset & Res.w1<Res.w1(col)+width+w1offset;
indw3 = Res.w3>Res.w3(row)-height+w3offset & Res.w3<Res.w3(row)+height+w3offset;

w1space = 1;  % resolution in w1 domain (sub sampling = any integer greater than 1)
w3space = 1;  % resolution in w3 domain (sub sampling = any integer greater than 1)
Res.w1cut = Res.w1(indw1);
Res.w1cut = Res.w1cut(1:w1space:end);
Res.w3cut = Res.w3(indw3);
Res.w3cut = Res.w3cut(1:w3space:end);

for i = 1:Res.numres
    Res.Rw1w3cut{i} = Res.Rw1w3{i}(indw3,indw1,:); 
    Res.Rw1w3cut{i} = Res.Rw1w3cut{i}(1:w3space:end,1:w1space:end,:);
end

% Calculate reponses more relevant to experimental measurements
Res.Rw1w3Rephasing = real(Res.Rw1w3cut{2}+Res.Rw1w3cut{3});
Res.Rw1w3NonRephasing = real(Res.Rw1w3cut{1}+Res.Rw1w3cut{4});
Res.Rw1w3Absorptive = real(Res.Rw1w3cut{1}+Res.Rw1w3cut{2}+Res.Rw1w3cut{3}+Res.Rw1w3cut{4});

%% Plot 2D spectra (visualization purposes)

plotting.tsel = 1;
plotting.sigs = [{'R4'} {'RABS'}]; 
plotting.opt = {};

resplot(Res,plotting)

%% Start of quantum beating analysis

% Checking to make sure beats are present
% click on 2D plot to choose (w1,w3) points
 
resnum = 2; % the number of the response function to analyze here

subplot(1,2,1);
contourf(Res.w1cut,Res.w3cut,normdim(real(Res.Rw1w3cut{1,resnum}(:,:,1))),'LineWidth',0.25,'LevelStep',0.1,'ButtonDownFcn',{@ButtonDownFcnTest,Res,resnum});
colormap(cmap2d(20))
ylabel('\omega_{3} (10^{3} cm^{-1})');
xlabel('\omega_{1} (10^{3} cm^{-1})');
set(gca,'BoxStyle','full','DataAspectRatio',[1 1 1],...
    'FontSize',18,'Layer','top');

%% Population subtraction

beatproc = cell(1,Res.numres);
npad = length(Res.t2)*Res.padnum;

for i = 1:Res.numres
        
    %real data
    beatproc{1,i}.rdatadimred = MDSDimRed(real(Res.Rw1w3cut{i}),Res.w1cut,Res.w3cut,Res.t2);
    [rexpfits,rOscFits] = OAS_sim(beatproc{1,i}.rdatadimred);
    beatproc{1,i}.rexpfits = unpackOAS(rexpfits,length(Res.w3cut),length(Res.w1cut),length(Res.t2));
    beatproc{1,i}.roscfits = unpackOAS(rOscFits,length(Res.w3cut),length(Res.w1cut),length(Res.t2));
    beatproc{1,i}.rfits = beatproc{1,i}.rexpfits; 
    beatproc{1,i}.riso = real(Res.Rw1w3cut{1,i})-beatproc{1,i}.rfits;
    
    %imaginary data
    beatproc{1,i}.idatadimred = MDSDimRed(imag(Res.Rw1w3cut{i}),Res.w1cut,Res.w3cut,Res.t2);
    [iexpfits,iOscFits] = OAS_sim(beatproc{1,i}.idatadimred);
    beatproc{1,i}.iexpfits = unpackOAS(iexpfits,length(Res.w3cut),length(Res.w1cut),length(Res.t2));
    beatproc{1,i}.ioscfits = unpackOAS(iOscFits,length(Res.w3cut),length(Res.w1cut),length(Res.t2));
    beatproc{1,i}.ifits = beatproc{1,i}.iexpfits;
    beatproc{1,i}.iiso = imag(Res.Rw1w3cut{1,i})-beatproc{1,i}.ifits;
    
    output = complexFFT_v2(beatproc{1,i}.riso,beatproc{1,i}.iiso,Res.t2,Res.w1cut,Res.w3cut,npad,0);
    beatproc{1,i}.dataw1w2w3 = output.data_w1w2w3;
    clear expfits OscFits output
end

%% Define beat frequency axis

Res.dw2 = 1 / (p.c*Res.dt2*npad);                                   %defining the frequency step [cm-1]
Res.w2 = ((-npad*Res.dw2)/2:Res.dw2:((npad*Res.dw2)/2-Res.dw2));    %frequency axis [cm-1%]
Res.w2 = -Res.w2;

%% Generate additional signals

%total rephasing (denoted as R7)
realRiso = beatproc{1,2}.riso+beatproc{1,3}.riso; %+beatproc{1,5}.riso;
imagRiso = beatproc{1,2}.iiso+beatproc{1,3}.iiso; %+beatproc{1,5}.iiso;
output = complexFFT_v2(realRiso,imagRiso,Res.t2,Res.w1cut,Res.w3cut,npad,0);
beatproc{1,5}.dataw1w2w3 = output.data_w1w2w3;
clear output

%total non-rephasing (denoted as R8)
realNRiso = beatproc{1,1}.riso+beatproc{1,4}.riso; %+beatproc{1,6}.riso;
imagNRiso = beatproc{1,1}.iiso+beatproc{1,4}.iiso; %+beatproc{1,6}.iiso;
output = complexFFT_v2(realNRiso,imagNRiso,Res.t2,Res.w1cut,Res.w3cut,npad,0);
beatproc{1,6}.dataw1w2w3 = output.data_w1w2w3;
clear output

%absorptive (denoted as R9)
Absiso = realRiso+realNRiso;
beatproc{1,7}.dataw1w2w3 = real_FFT_v1(Absiso,Res.t2,Res.w1cut,Res.w3cut,npad,0);
clear output

%% Examine quantum beats
%click on 2D plot to choose (w1,w3) points
 
resnum = 2;

h = figure;
subplot1 = subplot(2,2,1);
contourf(Res.w1cut,Res.w3cut,normdim(real(Res.Rw1w3cut{1,resnum}(:,:,1))),'LineWidth',0.25,'LevelStep',0.1,'ButtonDownFcn',{@ButtonDownFcnRes,Res,beatproc,resnum});
colormap(cmap2d(10))
hold on
line(min(Res.w1cut),max(Res.w1cut),'LineWidth',1.5);
ylabel('\omega_{3} (cm^{-1})');
xlabel('\omega_{1} (cm^{-1})');
title(Res.t2(1));
box on
set(gca,'BoxStyle','full','CLim',[-1 1],'DataAspectRatio',[1 1 1],...
    'FontSize',18,'Layer','top');
set (gcf, 'WindowButtonMotionFcn', @mouseMove);

%% Normalize to absolute maximum quantum beating

for i = 1:Res.numres+3
    beatproc{1,i}.dataw1w2w3norm = abs(beatproc{1,i}.dataw1w2w3)./max(max(max(abs(beatproc{1,i}.dataw1w2w3)))); %(:,:,148:188)
end

%% Power spectral norm over entire dataset

for i = 1:Res.numres+3
    beatproc{1,i}.norm = permute(vecnorm(vecnorm(beatproc{1,i}.dataw1w2w3norm,1),2),[1 3 2]);
end

leg = cell(1,Res.numres);

figure
for i = 1:Res.numres
   plot(Res.w2,beatproc{1,i}.norm,'LineWidth',1)
   leg{i} = (strcat("R",num2str(i)));
   hold on
end
set(gca, 'FontSize', 20);
xlim([-2000 2000])
legend(string(leg))
hold off

%% Plot beatmaps

opts.mapfreqs = [1350 1450];
opts.signal = 7;
opts.limcor = 0;
opts.misc = {};

plotbeatmaps(Res, beatproc, opts, p);