%=========================================================================
%                                                                     
%	TITLE: 
%       XCT - EXERCISE 2
%								
%	DESCRIPTION:						
%       COMPUTE SINOGRAM AND CT RECONSTRUCTIONS (FBP,FT)
%
%	INPUT:								
%       NONE	
%
%	OUTPUT:							
%       DISPLAY
%			
%	VERSION HISTORY:						
%	    150818SK INITIAL VERSION
%       170926SK UPDATE
%	    201020SK UPDATE
%
%=========================================================================

%=========================================================================
%	M A I N  F U N C T I O N
%=========================================================================
function [] = XCT_EXERCISE2()

    clear all; close all; 
    
   
    % --------------------------------------------------------------------
    % Display title
    % -------------------------------------------------------------------- 
    fprintf ( '-----------------------------------------\n' );  
    fprintf ( ' BIOMEDICAL IMAGING - XCT-EXERCISE #2\n' );  
    fprintf ( '-----------------------------------------\n' );  
    
    
    % --------------------------------------------------------------------
    % Set imaging parameters
    % --------------------------------------------------------------------
    matrix      = 256;                      % image matrix  [1pixel = 1mm]
    ua          = 50;                       % anode voltage [keV]
      
    
    % --------------------------------------------------------------------
    % Set density
    % (see: Table 2 in www.nist.gov/pml/data/xraycoef/)
    % --------------------------------------------------------------------
    rho_blood       = 1.060;                % density blood    [g/cm3]
    rho_bone        = 1.920;                % density bone     [g/cm3]
    rho_lung        = 0.001;                % density lung/air [g/cm3]
    rho_muscle      = 1.050;                % density muscle   [g/cm3]
    
   
    % --------------------------------------------------------------------
    % Set X-ray mass attenuation coefficients for 50 and 150 keV
    % (see: Table 4 in www.nist.gov/pml/data/xraycoef/)
    % --------------------------------------------------------------------
    mac_blood(1)    = 0.228;                % blood @  50 keV   [cm2/g]
    mac_blood(2)    = 0.149;                % blood @ 150 keV   [cm2/g]
    
    mac_bone(1)     = 0.424;                % bone @  50 keV    [cm2/g]
    mac_bone(2)     = 0.148;                % bone @ 150 keV    [cm2/g]
    
    mac_lung(1)     = 0.208;                % lung @  50 keV    [cm2/g]
    mac_lung(2)     = 0.136;                % lung @ 150 keV    [cm2/g]
    
    mac_muscle(1)   = 0.226;                % muscle @  50 keV  [cm2/g]
    mac_muscle(2)   = 0.149;                % muscle @ 150 keV  [cm2/g]
    
    
    % --------------------------------------------------------------------
    % Calculate linear attenuation coefficients
    % --------------------------------------------------------------------
    mue_blood       = rho_blood.*mac_blood(:);    
    mue_bone        = rho_bone.*mac_bone(:); 
    mue_lung        = rho_lung.*mac_lung(:); 
    mue_muscle      = rho_muscle.*mac_muscle(:); 
    
    idx = 1; if ua==150 idx = 2; end        % set index
   
    
    % --------------------------------------------------------------------
    % Define analytical phantom using ellipses with [x0 y0 a b phi mue]
    %
    %       x0,y0   - center point [cm] (+x -> left-right, +y -> bottom-up)
    %       a,b     - half axes [cm]
    %       theta   - rotation angle relative to x-axis [deg]
    %       mue     - linear attenuation coefficient [1/cm]
    % --------------------------------------------------------------------
    phantom.ellipse = [   0   0   90  80  0   mue_muscle(idx);                  % thorax
                          0   0   70  60  0   mue_lung(idx)-mue_muscle(idx);    % lung
                       +110   0   15  15  0   mue_muscle(idx);                  % left arm muscle
                       +110   0    5   5  0   mue_bone(idx)-mue_muscle(idx);    % left arm bone
                       -110   0   15  15  0   mue_muscle(idx);                  % right arm muscle
                       -110   0    5   5  0   mue_bone(idx)-mue_muscle(idx);    % right arm bone
                          0   0   10  10  0   mue_blood(idx)-mue_lung(idx);     % aorta 
                        +30 +25   25  20 35   mue_muscle(idx)-mue_lung(idx)];   % heart
    
                    
    % --------------------------------------------------------------------
    % Compute phantom and display
    % --------------------------------------------------------------------
    [x,y]            = meshgrid(-fix(matrix/2):+fix(matrix/2));
    phantom.discrete = CalcDiscretePhantom(x,y,phantom,ua);

    DisplayData(phantom.discrete,[2,3,1]); title('Phantom'); 
    
    
    % --------------------------------------------------------------------
    % Compute and display projections and sinogram
    % --------------------------------------------------------------------
    projection_angles   = 0:179; idx = 1;                                  
    phantom.sino        = zeros(matrix+1,length(projection_angles));

    for phi=projection_angles

        [r,phi_] = meshgrid(-fix(matrix/2):1:+fix(matrix/2),phi);  
        phantom.sino(:,idx) = CalcLineIntegrals(r,phi_,phantom,ua)/(matrix+1);
         
        DisplayData(phantom.sino(:,idx),[2,3,2]); 
        title(sprintf('Projection (phi=%d)',phi)); xlabel('r'); ylabel('mue(r)'); 
        
        DisplayData(phantom.sino,[2,3,3]); 
        title('Sinogram'); xlabel('phi'); ylabel('r'); drawnow;
        
        idx = idx+1;
    end
    
    
    % --------------------------------------------------------------------
    % Reconstruct image using simple back-projection (SBP)
    % --------------------------------------------------------------------
    [x,y]        = meshgrid(-fix(matrix/2):1:+fix(matrix/2)); idx = 1; 
    phantom.sbp  = zeros(matrix+1,matrix+1);
    
    for phi=projection_angles
        
        % TASK 2.2 FILL IN HERE
        
        r = x*cosd(phi)+y*sind(phi);
        rs = round(r)+129;
        ix = find((0<rs) & (rs<=257));
        
        % TASK 2.2 UNCOMMENT HERE
        projection = phantom.sino(:,idx);
        phantom.sbp(ix) = phantom.sbp(ix)+projection(rs(ix));
        
        DisplayData(phantom.sbp,[2,3,4]); 
        title('Backprojection'); xlabel('x'); ylabel('y'); drawnow;
        
        idx = idx+1;
    end
    
 
    % --------------------------------------------------------------------
    % Reconstruct image using Filtered Back-Projection (FBP)
    % --------------------------------------------------------------------
    [x,y]       = meshgrid(-fix(matrix/2):1:+fix(matrix/2)); idx = 1; 
    phantom.fbp = zeros(matrix+1,matrix+1);
    
    filter      = CalcFilter(matrix+1);   % TASK 2.3 EDIT 
    conv_freq_domain = zeros(matrix+1,length(projection_angles));
    
    for phi=projection_angles
        
        % TASK 2.3 FILL IN HERE
        projection  = phantom.sino(:,idx);
        projection_fourier  = R2U(projection');
        
        r = x*cosd(phi)+y*sind(phi);
        rs = round(r)+129;
        ix = find((0<rs) & (rs<=257));
        
        conv_freq_domain(rs(ix),idx) = projection_fourier(1, rs(ix)).*filter(1, rs(ix));
        filteredprojection = U2R(conv_freq_domain(:,idx));
        
        % TASK 2.3 UNCOMMENT HERE
        phantom.fbp(ix) = phantom.fbp(ix)+filteredprojection(rs(ix));
           
        DisplayData(phantom.fbp,[2,3,5]); 
        title('Filtered Backprojection'); xlabel('x'); ylabel('y'); drawnow;
        
        idx = idx+1;
    end  
 
    
    % --------------------------------------------------------------------
    % Reconstruct image using Fast Fourier Transform (FFT)
    % --------------------------------------------------------------------
    [r]         = meshgrid(-fix(matrix/2):1:+fix(matrix/2),1); idx = 1; 
    phantom.fft = zeros(matrix+1,matrix+1);
    
    for phi=projection_angles
    
        % TASK 2.4 FILL IN HERE
        
        r = x*cosd(phi)+y*sind(phi);
        rs = round(r)+129;
        ix = find((0<rs) & (rs<=257));
        
        projection = phantom.sino(:,idx);
        phantom.sbp(ix) = phantom.sbp(ix)+projection(rs(ix));
        
        phantom.fft(ix) = R2U(phantom.sbp(ix));
        
        idx = idx+1;
        
    end
    
    phantom.fft = U2R(phantom.fft);

    DisplayData(phantom.fft*length(projection_angles),[2,3,6]); 
    title('2D Fourier transform'); xlabel('x'); ylabel('y'); drawnow;
    
end


%=========================================================================
function [image] = CalcDiscretePhantom(x,y,phantom,ua)
    
    image = zeros(size(x));
    
    for k = 1:length(phantom.ellipse(:,1))
        
        theta   = phantom.ellipse(k,5)*pi/180;
        
        X0      = [x(:)'-phantom.ellipse(k,1);y(:)'-phantom.ellipse(k,2)];
        D       = [1/phantom.ellipse(k,3) 0;0 1/phantom.ellipse(k,4)];
        Q       = [cos(theta) sin(theta); -sin(theta) cos(theta)];
         
        % ----------------------------------------------------------------
        % Find inside of ellipse given by X0,D,Q
        % ----------------------------------------------------------------
        equ = sum((D*Q*X0).^2);
        i = find(equ <= 1);
         
        % ----------------------------------------------------------------
        % Assign linear attenuation coefficients
        % ----------------------------------------------------------------
        image(i) = image(i)+phantom.ellipse(k,6);
        
    end
end


%=========================================================================
function [projection] = CalcLineIntegrals(r,phi,phantom,ua)
    
    projection  = zeros(size(r));
    phi         = phi/180*pi;
    
    sinphi  = sin(phi(:)); 
    cosphi  = cos(phi(:));
    
    rx      = r(:).*cosphi; 
    ry      = r(:).*sinphi;
    
    for k=1:length(phantom.ellipse(:,1))
        
        x0      = phantom.ellipse(k,1); y0 = phantom.ellipse(k,2);
        a       = phantom.ellipse(k,3); b  = phantom.ellipse(k,4);
        
        theta   = phantom.ellipse(k,5)*pi/180; 
        mue     = phantom.ellipse(k,6);
        
        r0      = [rx-x0,ry-y0]';
        
        % ----------------------------------------------------------------
        % Find entry and exit points of each ellipse
        % ----------------------------------------------------------------
        DQ      = [cos(theta)/a sin(theta)/a; -sin(theta)/b cos(theta)/b];
        DQphi   = DQ*[sinphi,-cosphi]'; 
        DQr0    = DQ*r0;
        
        A       = sum(DQphi.^2); 
        B       = 2*sum(DQphi.*DQr0);
        C       = sum(DQr0.^2)-1; 
        equ     = B.^2-4*A.*C;
        
        i       = find(equ>0);
        s1      = 0.5*(-B(i)+sqrt(equ(i)))./A(i);
        s2      = 0.5*(-B(i)-sqrt(equ(i)))./A(i);
        
        % ----------------------------------------------------------------
        % Compute projection i.e. mue * length inside object
        % ----------------------------------------------------------------
        proj    = mue*abs(s1-s2);
        
        % ----------------------------------------------------------------
        % Sum projections
        % ----------------------------------------------------------------
        projection(i) = projection(i)+proj;
    end
    
    projection = reshape(projection,size(r));
    
end


%=========================================================================
function [filter] = CalcFilter(matrix)

    % --------------------------------------------------------------------
    % Define high-pass filter according to |u|
    % --------------------------------------------------------------------
    hpf = abs(-fix(matrix/2):+fix(matrix/2));
      
    % TASK 2.3 FILL IN HERE
    hpf = hpf/128; % High-Pass Filter Frequency Domain
    
    lpf = zeros(1,matrix); % Low-Pass Filter Frequency Domain
    for n=0:1:256
        lpf(1,n+1) = 0.54 - (1-0.54)*cos((2*pi*n)/256); % Hamming Window Function
    end
    
    filter = hpf.*lpf; % Combine HPF and LPF to create a BPF
    
%     % Plot Fourier Pairs of the Ideal HPF and the Modified BPF
%     figure;
%     
%     subplot(2,2,1);plot(U2R(hpf));title('High-Pass Filter - Image Domain');xlabel('r');ylabel('Filter(r)');axis tight;
%     subplot(2,2,2);plot(hpf);title('High-Pass Filter - Image Frequency Domain');xlabel('u');ylabel('FT(Filter)');axis tight;
%     subplot(2,2,3);plot(U2R(filter));title('Band-Pass Filter - Image Domain');xlabel('r');ylabel('Filter(r)');axis tight;
%     subplot(2,2,4);plot(filter);title('Band-Pass Filter - Image Frequency Domain');xlabel('u');ylabel('FT(Filter)');axis tight;
    
        
end


%=========================================================================
%=========================================================================