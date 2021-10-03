%=========================================================================
%                                                                     
%	TITLE: 
%       NUC - EXERCISE 1
%								
%	DESCRIPTION:						
%       GENERATE AND RECONSTRUCT PET/CT DATA
%
%	INPUT:								
%       NONE	
%
%	OUTPUT:							
%       DISPLAY
%			
%	VERSION HISTORY:						
%	    191101K INITIAL VERSION
%	    201081K UPDATE
%
%=========================================================================

%=========================================================================
%	M A I N  F U N C T I O N
%=========================================================================
function [] = NUC_EXERCISE1()

    clear all; close all; 
    
    % --------------------------------------------------------------------
    % Display title
    % -------------------------------------------------------------------- 
    fprintf ( '-----------------------------------------\n' );  
    fprintf ( ' BIOMEDICAL IMAGING - NUC-EXERCISE #1\n' );  
    fprintf ( '-----------------------------------------\n' );  
     
    % --------------------------------------------------------------------
    % Set tissue density
    % (see: Table 2 in www.nist.gov/pml/data/xraycoef/)
    % --------------------------------------------------------------------
    rho_blood       = 1.060;                % density blood     [g/cm3]
    rho_bone        = 1.450;                % density bone      [g/cm3]
    rho_lung        = 0.001;                % density lung/air  [g/cm3]
    rho_muscle      = 1.050;                % density muscle    [g/cm3]
    rho_water       = 1.000;                % density water     [g/cm3]
     
    % --------------------------------------------------------------------
    % Set mass attenuation coefficients for 100 keV (CT) and 511 keV (PET)
    % (see: Table 4 in www.nist.gov/pml/data/xraycoef/)
    % --------------------------------------------------------------------
    mac_blood(1)    = 0.169;                % blood  @ 100 keV  [cm2/g]
    mac_blood(2)    = 9.598E-02;            % blood  @ 511 keV  [cm2/g]
    mac_bone(1)     = 0.186;                % bone   @ 100 keV  [cm2/g]
    mac_bone(2)     = 9.022E-02;            % bone   @ 511 keV  [cm2/g]
    mac_lung(1)     = 0.154;                % lung   @ 100 keV  [cm2/g]
    mac_lung(2)     = 8.712E-02;            % lung   @ 511 keV  [cm2/g]
    mac_muscle(1)   = 0.169;                % muscle @ 100 keV  [cm2/g]
    mac_muscle(2)   = 9.598E-02;            % muscle @ 511 keV  [cm2/g]
    mac_water(1)    = 0.171;                % water  @ 100 keV  [cm2/g]
    mac_water(2)    = 9.687E-02;            % water  @ 511 keV  [cm2/g]
    
    % --------------------------------------------------------------------
    % Task 1.3. - for demonstration only; comment out otherwise
    % --------------------------------------------------------------------
    %mac_bone        = 10*mac_bone;               
    
    % --------------------------------------------------------------------
    % Calculate linear attenuation coefficients
    % --------------------------------------------------------------------
    mue_blood       = rho_blood.*mac_blood(:);    
    mue_bone        = rho_bone.*mac_bone(:); 
    mue_lung        = rho_lung.*mac_lung(:); 
    mue_muscle      = rho_muscle.*mac_muscle(:);
    mue_water       = rho_water.*mac_water(:);
      
    % --------------------------------------------------------------------
    % Set imaging parameters
    % --------------------------------------------------------------------
    matrix = 256;                           % image matrix [1pix=1mm]
    
    
    % --------------------------------------------------------------------
    % PET SIMULATION (begin)
    % --------------------------------------------------------------------
    % --------------------------------------------------------------------
    % Define analytical phantom using ellipses with [x0 y0 a b phi mue act]
    %
    %       x0,y0   - center point [cm] (+x -> left-right, +y -> bottom-up)
    %       a,b     - half axes [cm]
    %       theta   - rotation angle relative to x-axis [deg]
    %       mue     - linear attenuation coefficient [1/cm]
    %       act     - radioactivity [MBq]
    %
    % --------------------------------------------------------------------
    phantom.ellipse = [   0   0   90  80  0   mue_bone(2)               0;      % thorax
                          0   0   70  60  0   mue_lung(2)-mue_bone(2)   0;      % lung
                       +110   0   15  15  0   mue_muscle(2)             0;      % left arm muscle
                       +110   0    5   5  0   mue_bone(2)-mue_muscle(2) 0;      % left arm bone
                       -110   0   15  15  0   mue_muscle(2)             0;      % right arm muscle
                       -110   0    5   5  0   mue_bone(2)-mue_muscle(2) 0;      % right arm bone
                          0   0   10  10  0   mue_blood(2)-mue_lung(2)  0;      % aorta 
                        +30 +25   25  20 35   mue_muscle(2)-mue_lung(2) 0;      % heart
                        -30 +25   10  10  0   mue_muscle(2)-mue_lung(2) 300];   % lung tumor

    
    % --------------------------------------------------------------------
    % Compute and display PET phantom 
    % --------------------------------------------------------------------
    phantom.discrete = CalcPetPhantom(matrix,phantom);
    DisplayData(phantom.discrete,[2,4,1]); title('PET Activity'); colormap(jet);
    
    % --------------------------------------------------------------------
    % Display true PET signal value of lung tumor (=9) on command line
    % --------------------------------------------------------------------
    fprintf('PET tumor activity ground truth: %f\n',CalcROISignal(phantom,9,phantom.discrete));
    
    % --------------------------------------------------------------------
    % Compute PET sinogram
    % --------------------------------------------------------------------
    projection_angles = 0:179;       
    phantom.petsino   = CalcPetSinogram(projection_angles,matrix,phantom); 
   
    % --------------------------------------------------------------------
    % Add PET noise
    % --------------------------------------------------------------------    
    phantom.petsino = poissrnd(phantom.petsino,size(phantom.petsino));
    
    % --------------------------------------------------------------------
    % Display PET sinogram
    % --------------------------------------------------------------------
    DisplayData(phantom.petsino,[2,4,2]); title('PET Sinogram'); xlabel('phi'); ylabel('r'); 

    % --------------------------------------------------------------------
    % Compute and display PET reconstruction (FBP)
    % --------------------------------------------------------------------
    filter            = CalcFilter(matrix+1);
    phantom.petfbp    = CalcFBPRecon(projection_angles,matrix,phantom.petsino,filter);
    DisplayData(phantom.petfbp,[2,4,3]); title('PET Reconstruction'); xlabel('x'); ylabel('y');  

    % --------------------------------------------------------------------
    % Display reconstructed PET signal value of lung tumor (=9) 
    % --------------------------------------------------------------------
    fprintf('PET tumor activity w/o correction: %f\n',CalcROISignal(phantom,9,phantom.petfbp));

    
    % --------------------------------------------------------------------
    % CT SIMULATION (begin)
    % --------------------------------------------------------------------
    % --------------------------------------------------------------------
    % Define analytical phantom using ellipses with [x0 y0 a b phi mue act]
    %
    %       x0,y0   - center point [cm] (+x -> left-right, +y -> bottom-up)
    %       a,b     - half axes [cm]
    %       theta   - rotation angle relative to x-axis [deg]
    %       mue     - linear attenuation coefficient [1/cm]
    %       act     - radioactivity [MBq]
    %
    % --------------------------------------------------------------------
    % --------------------------------------------------------------------
    % Update [mue act] of analytical phantom [mue act]
    %
    %       mue     - linear attenuation coefficient [1/cm]
    %       act     - radioactivity [MBq]
    %
    % --------------------------------------------------------------------
    phantom.ellipse(:,6:7) = [  mue_bone(1)               0;      % thorax
                                mue_lung(1)-mue_bone(1)   0;      % lung
                                mue_muscle(1)             0;      % left arm muscle
                                mue_bone(1)-mue_muscle(1) 0;      % left arm bone
                                mue_muscle(1)             0;      % right arm muscle
                                mue_bone(1)-mue_muscle(1) 0;      % right arm bone
                                mue_blood(1)-mue_lung(1)  0;      % aorta 
                                mue_muscle(1)-mue_lung(1) 0;      % heart
                                mue_muscle(1)-mue_lung(1) 0];     % lung tumor
                    
    % --------------------------------------------------------------------
    % Compute and display CT phantom 
    % --------------------------------------------------------------------
    phantom.discrete = CalcCtPhantom(matrix,phantom,1:length(phantom.ellipse(:,1)));
    DisplayData(phantom.discrete,[2,4,5]); title('CT Object'); 
    
    % --------------------------------------------------------------------
    % Compute CT sinogram
    % --------------------------------------------------------------------
    projection_angles = 0:179;       
    phantom.ctsino   = CalcCtSinogram(projection_angles,matrix,phantom);  
    
    % --------------------------------------------------------------------
    % Add CT noise
    % --------------------------------------------------------------------    
    N0 = 1e10;
    N  = N0*exp(-phantom.ctsino);
    N  = poissrnd(N,size(N));
    phantom.ctsino = -log(N/N0);
    
    % --------------------------------------------------------------------
    % Display CT sinogram
    % --------------------------------------------------------------------
    DisplayData(phantom.ctsino,[2,4,6]); title('CT Sinogram'); xlabel('phi'); ylabel('r'); 
        
    % --------------------------------------------------------------------
    % Compute and display PET reconstruction (FBP)
    % --------------------------------------------------------------------
    filter           = CalcFilter(matrix+1);
    phantom.ctfbp    = CalcFBPRecon(projection_angles,matrix,phantom.ctsino,filter);
    DisplayData(phantom.ctfbp,[2,4,7]); title('CT Reconstruction'); xlabel('x'); ylabel('y'); 
    
    
    % --------------------------------------------------------------------
    % PET ATTENUATION CORRECTION (begin)
    % --------------------------------------------------------------------
    % --------------------------------------------------------------------
    % TASK 1.3 Map linear attenuation coefficients from 100 keV to 511 keV
    phantom.corr = (5.7*(mue_bone(1)-mue_water(1))).*phantom.ctfbp;
        
    % --------------------------------------------------------------------
    % Recompute CT sinogram from remapped discrete CT image 
    % --------------------------------------------------------------------
    phantom.ctsino = radon(phantom.corr,projection_angles);
    
    % --------------------------------------------------------------------
    % Normalize for number of projection angles in CT sinogram
    % --------------------------------------------------------------------
    phantom.ctsino = phantom.ctsino/180/sqrt(2);
    
    % --------------------------------------------------------------------
    % Recompute PET sinogram from discrete PET image
    % --------------------------------------------------------------------
    phantom.petsino = radon(phantom.petfbp,projection_angles);
      
    % --------------------------------------------------------------------
    % TASK 1.3 Apply attenuation correction to PET sinogram
    phantom.petsino = phantom.petsino.*exp(phantom.ctsino);
    DisplayData(phantom.petsino,[2,4,8]); title('Attenuation Corrected PET Sinogram'); xlabel('phi'); ylabel('r');
    
    % --------------------------------------------------------------------
    % Recompute PET image from attenuation-corrected sinogram
    phantom.petcorr = iradon(phantom.petsino,projection_angles);
 
    % --------------------------------------------------------------------
    % Display PET image with attenuation correction
    % --------------------------------------------------------------------
    DisplayData(phantom.petcorr,[2,4,4]); title('Recon with Attenuation Correction'); xlabel('x'); ylabel('y');  
 
    % --------------------------------------------------------------------
    % Display corrected PET signal value of lung tumor (=9)
    % --------------------------------------------------------------------
    fprintf('PET tumor activity with correction: %f\n',CalcROISignal(phantom,9,phantom.petcorr(1:matrix+1,1:matrix+1)));
 
end

%=========================================================================
function [signal] = CalcROISignal(phantom,k,image)

    % --------------------------------------------------------------------
    % Find region-of-interest (ROI)
    % --------------------------------------------------------------------
    [x,y]   = meshgrid(-fix(size(image,1)/2):+fix(size(image,1)/2)); 
    theta   = phantom.ellipse(k,5)*pi/180;

    
    X0      = [x(:)'-phantom.ellipse(k,1);y(:)'-phantom.ellipse(k,2)];
    D       = [1/(phantom.ellipse(k,3)-3) 0;0 1/(phantom.ellipse(k,4)-3)];
    Q       = [cos(theta) sin(theta); -sin(theta) cos(theta)];

    equ     = sum((D*Q*X0).^2);
    i       = find(equ<=1);

    % --------------------------------------------------------------------
    % Calculate mean signal
    % --------------------------------------------------------------------    
    signal = mean(image(i)); 
        
end


%=========================================================================
%=========================================================================