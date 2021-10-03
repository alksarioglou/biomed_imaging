%=========================================================================
%                                                                     
%	TITLE: 
%       XCT - EXERCISE 1
%								
%	DESCRIPTION:						
%       COMPUTE PROJECTIONS OF ANALYTICAL PHANTOM 
%
%	INPUT:								
%       NONE	
%
%	OUTPUT:							
%       DISPLAY
%			
%	VERSION HISTORY:						
%	    150816SK INITIAL VERSION
%	    191020SK UPDATE
%       201018SK UPDATE
%
%=========================================================================

%=========================================================================
%	M A I N  F U N C T I O N
%=========================================================================
function [] = XCT_EXERCISE1()

    clear all; close all; 
    
    
    % --------------------------------------------------------------------
    % Display title
    % -------------------------------------------------------------------- 
    fprintf ( '-----------------------------------------\n' );  
    fprintf ( ' BIOMEDICAL IMAGING - XCT-EXERCISE #1\n' );  
    fprintf ( '-----------------------------------------\n' );  
    
  
    % --------------------------------------------------------------------
    % Set imaging parameters
    % --------------------------------------------------------------------
    matrix      = 256;                      % image matrix  [1pixel = 1mm]
    ua          = 50;                       % anode voltage [keV]
    
    
    % --------------------------------------------------------------------
    % TASK 1.1 (begin)
    % --------------------------------------------------------------------
    % Set density 
    % (see: Table 2 in www.nist.gov/pml/data/xraycoef/)
    % --------------------------------------------------------------------
    
    % TASK 1.1 EDIT HERE
    rho_blood       = 1.060E+00;                % density blood    [g/cm3]
    rho_bone        = 1.920E+00;                % density bone     [g/cm3]
    rho_lung        = 0.001;                    % density lung/air [g/cm3]
    rho_muscle      = 1.050E+00;                % density muscle   [g/cm3]
    
     
    % --------------------------------------------------------------------
    % Set X-ray mass attenuation coefficients for 50 and 150 keV
    % (see: Table 4 in www.nist.gov/pml/data/xraycoef/)
    % --------------------------------------------------------------------
    
    % TASK 1.1 EDIT HERE
    mac_blood(1)    = 2.278E-01;                % blood @  50 keV   [cm2/g]
    mac_blood(2)    = 1.492E-01;                % blood @ 150 keV   [cm2/g]
    
    mac_bone(1)     = 4.242E-01;                % bone @  50 keV    [cm2/g]
    mac_bone(2)     = 1.480E-01;                % bone @ 150 keV    [cm2/g]
    
    mac_lung(1)     = 0.208;                    % lung @  50 keV    [cm2/g]
    mac_lung(2)     = 0.136;                    % lung @ 150 keV    [cm2/g]
    
    mac_muscle(1)   = 2.262E-01;                % muscle @  50 keV  [cm2/g]
    mac_muscle(2)   = 1.492E-01;                % muscle @ 150 keV  [cm2/g]
    
    
    % --------------------------------------------------------------------
    % Calculate linear attenuation coefficients
    % --------------------------------------------------------------------
    
    % TASK 1.1 EDIT HERE
    mue_blood       = mac_blood.*rho_blood;  
    mue_bone        = mac_bone.*rho_bone;  
    mue_lung        = mac_lung.*rho_lung;  
    mue_muscle      = mac_muscle.*rho_muscle;  
    
    idx = 1; if ua==150 idx = 2; end
   
    
    % --------------------------------------------------------------------
    % Define test phantom using ellipses with [x0 y0 a b phi mue]
    %
    %       x0,y0   - center point [cm] (+x -> left-right, +y -> bottom-up)
    %       a,b     - half axes [cm]
    %       theta   - rotation angle relative to x-axis [deg]
    %       mue     - linear attenuation coefficient [1/cm]
    % --------------------------------------------------------------------
    
    % TASK 1.1 EDIT HERE
    phantom.ellipse = [50 60 20 20 0 mue_muscle(idx);
                       0 0 20 20 0 mue_bone(idx);
                       -50 -60 20 20 0 mue_blood(idx);
    
    ];     
    
    
    % --------------------------------------------------------------------
    % Compute and display discrete phantom 
    % --------------------------------------------------------------------
   [x,y] = meshgrid(-fix(matrix/2):+fix(matrix/2));       
   
    phantom.discrete = CalcDiscretePhantom(x,y,phantom,ua); 
  
    DisplayData(phantom.discrete,[1,4,1]); title('Discrete phantom'); 
    
    
    % --------------------------------------------------------------------
    % Display profile of mue(x)
    % --------------------------------------------------------------------
    DisplayData(phantom.discrete(fix(matrix/2),:)',[1,4,2]); 
    title('Phantom (profile)'); xlabel('x'); ylabel('mue(x)');
    
    
    % --------------------------------------------------------------------
    % Display column-wise sum of mue(y)
    % --------------------------------------------------------------------
    
    % TASK 1.1 FILL IN HERE
    DisplayData(sum(phantom.discrete(:,:),1)',[1,4,3]);
    title('Phantom (projection)'); xlabel('x'); ylabel('mue(x)');
    
  
    % --------------------------------------------------------------------
    % Compute and display projection
    % --------------------------------------------------------------------
    [r,phi] = meshgrid(-fix(matrix/2):+fix(matrix/2),[0]);  
        
    phantom.projection = CalcLineIntegrals(r,phi,phantom,ua);
    
    DisplayData(phantom.projection',[1,4,4]); 
    title('Projection'); xlabel('r'); ylabel('P(r))'); 
    
    waitforbuttonpress;                                     % pause here
   
   
    % --------------------------------------------------------------------
    % TASK 1.3 (begin)
    % --------------------------------------------------------------------
    % Define thorax phantom using ellipses with [x0 y0 a b phi mue]
    %
    %       x0,y0   - center point [cm] (+x -> left-right, +y -> bottom-up)
    %       a,b     - half axes [cm]
    %       theta   - rotation angle relative to x-axis [deg]
    %       mue     - linear attenuation coefficient [1/cm]
    % --------------------------------------------------------------------
    mue_lung_new = mue_lung(idx)-mue_muscle(idx);
    mue_armbone_new = mue_bone(idx)-mue_muscle(idx);
    mue_aorta_new = mue_blood(idx)-mue_lung(idx);
    mue_heart_new = mue_muscle(idx)-mue_lung(idx);
    
    % TASK 1.3 EDIT HERE
    phantom.ellipse = [   0   0   90  80  0   mue_muscle(idx);                                      % thorax
                          0   0   70  60  0   mue_lung_new;                                         % lung
                         -110 0   15  10  0   mue_muscle(idx);                                      % left arm muscle
                         -110 0   5   5   0   mue_armbone_new;                                      % left arm bone
                          110 0   15  10  0   mue_muscle(idx);                                      % right arm muscle
                          110 0   5   5   0   mue_armbone_new;                                      % right arm bone
                          0   0   15  10  0   mue_aorta_new;                                        % aorta 
                          30  30  20  10  20  mue_heart_new];                                       % heart
    
     
    % --------------------------------------------------------------------
    % Compute phantom, projection and display
    % --------------------------------------------------------------------
    [x,y]   = meshgrid(-fix(matrix/2):+fix(matrix/2));
    [r,phi] = meshgrid(-fix(matrix/2):+fix(matrix/2),[0]);  
   
    phantom.discrete    = CalcDiscretePhantom(x,y,phantom,ua);
    phantom.projection  = CalcLineIntegrals(r,phi,phantom,ua);

    phantom_no_contrast = phantom.discrete;
    project_no_contrast = phantom.projection;
    
    % TASK 1.3 UNCOMMENT HERE
    DisplayData(phantom.discrete,[1,4,1]); 
    title('Phantom (w/o contrast agent)'); 
    
    % --------------------------------------------------------------------
    % Display profile of mue(x)
    % --------------------------------------------------------------------
    DisplayData(phantom.discrete(fix(matrix/2),:)',[1,4,2]); 
    title('Phantom (profile)'); xlabel('x'); ylabel('mue(x)');
    
    
    % --------------------------------------------------------------------
    % Display column-wise sum of mue(y)
    % --------------------------------------------------------------------
    
    DisplayData(phantom.discrete(:,fix(matrix/2)),[1,4,3]);
    title('Phantom (projection)'); xlabel('y'); ylabel('mue(y)');
    
  
    % --------------------------------------------------------------------
    % Compute and display projection
    % --------------------------------------------------------------------
    
    DisplayData(phantom.projection',[1,4,4]); 
    title('Projection (w/o contrast agent)'); xlabel('r'); ylabel('P(r))'); 
    
    waitforbuttonpress;                                     % pause here
   
    
    % --------------------------------------------------------------------
    % Recompute phantom and projections (blood with contrast agent)
    % --------------------------------------------------------------------
    
    % TASK 1.3 EDIT HERE
    phantom.ellipse(7,6) = 2*mue_blood(idx);   % mue of iodine contrast agent in blood              
    
    phantom.discrete    = CalcDiscretePhantom(x,y,phantom,ua);
    phantom.projection  = CalcLineIntegrals(r,phi,phantom,ua);
    
   
    % --------------------------------------------------------------------
    % Subtract images, projections and display (with contrast agent)
    % --------------------------------------------------------------------
    
    % TASK 1.3 FILL IN HERE
    phantom_dsa = phantom.discrete - phantom_no_contrast;
    project_dsa = phantom.projection - project_no_contrast;
   
    % TASK 1.3 UNCOMMENT HERE
    DisplayData(phantom.projection',[1,4,2]); 
    title('Projection (w/o contrast agent)'); xlabel('r'); ylabel('P(r))'); 
    
    DisplayData(phantom_dsa,[1,4,3]); 
    title('Phantom (DSA)'); 
                    
    DisplayData(project_dsa',[1,4,4]); 
    title('Projection (DSA)'); xlabel('r'); ylabel('P(r))');
          
end


%=========================================================================
function [image] = CalcDiscretePhantom(x,y,phantom,ua)
    
    % --------------------------------------------------------------------
    % See Figure 1 in task sheet for definitions
    % --------------------------------------------------------------------
    image = zeros(size(x));
    
    % --------------------------------------------------------------------
    % Loop over all ellipses defined in phantom.ellipse
    % --------------------------------------------------------------------
    for k = 1:length(phantom.ellipse(:,1))      
        
        % ----------------------------------------------------------------
        % Angle theta as defined in Figure 1 in task sheet
        % ----------------------------------------------------------------
        theta   = phantom.ellipse(k,5)*pi/180;
        
        % ----------------------------------------------------------------
        % Center point (x0,y0) of ellipse as defined in Figure 1 in task
        % ----------------------------------------------------------------
        X0      = [x(:)'-phantom.ellipse(k,1);y(:)'-phantom.ellipse(k,2)];
        
        % ----------------------------------------------------------------
        % Half axis (a,b) of ellipse
        % ----------------------------------------------------------------
        D       = [1/phantom.ellipse(k,3) 0;0 1/phantom.ellipse(k,4)];
        
        % ----------------------------------------------------------------
        % Rotation of ellipse with angle theta 
        % ----------------------------------------------------------------
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
    
    % --------------------------------------------------------------------
    % See Figure 1 in task sheet for definitions
    % --------------------------------------------------------------------
    projection  = zeros(size(r));
    phi         = phi/180*pi;
    
    sinphi  = sin(phi(:)); 
    cosphi  = cos(phi(:));
    
    rx      = r(:).*cosphi; 
    ry      = r(:).*sinphi;
     
    % --------------------------------------------------------------------
    % Loop over all ellipses defined in phantom.ellipse
    % --------------------------------------------------------------------
    for k=1:length(phantom.ellipse(:,1))
        
        % ----------------------------------------------------------------
        % Center point (x0,y0) of ellipse as defined in Figure 1 in task
        % ----------------------------------------------------------------
        x0      = phantom.ellipse(k,1); y0 = phantom.ellipse(k,2);
        
        % ----------------------------------------------------------------
        % Half axis (a,b) of ellipse as defined in Figure 1 in task
        % ----------------------------------------------------------------
        a       = phantom.ellipse(k,3); b  = phantom.ellipse(k,4);
        
        % ----------------------------------------------------------------
        % Angle theta as defined in Figure 1 in task
        % ----------------------------------------------------------------
        theta   = phantom.ellipse(k,5)*pi/180; 
        
        % ----------------------------------------------------------------
        % Linear attenuation coefficient
        % ----------------------------------------------------------------
        mue     = phantom.ellipse(k,6);
        
        % ----------------------------------------------------------------
        % Matrix indices 
        % ----------------------------------------------------------------
        r0      = [rx-x0,ry-y0]';
        
        % ----------------------------------------------------------------
        % Rotation including ellipse half axes
        % ----------------------------------------------------------------
        DQ      = [cos(theta)/a sin(theta)/a; -sin(theta)/b cos(theta)/b];
        
        % ----------------------------------------------------------------
        % Include projection angle phi
        % ----------------------------------------------------------------
        DQphi   = DQ*[sinphi,-cosphi]'; 

        % ----------------------------------------------------------------
        % Include ellipse position
        % ----------------------------------------------------------------
        DQr0    = DQ*r0;
        
        % ----------------------------------------------------------------
        % Solve quadratic equation to find intersection of ellipse with
        % line of projection (X-ray source - detector)
        % ----------------------------------------------------------------
        A       = sum(DQphi.^2); 
        B       = 2*sum(DQphi.*DQr0);
        C       = sum(DQr0.^2)-1; 
        equ     = B.^2-4*A.*C;
        i       = find(equ>0);
        
        sp      = 0.5*(-B(i)+sqrt(equ(i)))./A(i);
        sq      = 0.5*(-B(i)-sqrt(equ(i)))./A(i);
        
        
        % ----------------------------------------------------------------
        % TASK 1.2 (begin)
        % ----------------------------------------------------------------
        
        % TASK 1.2 FILL IN/EDIT HERE
  
% My solution        
%         lsf     = @(y) mue+0*y;
%         
%         for m=1:length(i)
%             proj    = ua*exp(-integral(lsf,sq(m),sp(m)));
%         end
        
        proj    = mue*abs(sp-sq); 

        % ----------------------------------------------------------------
        % TASK 1.2 (end)
        % ----------------------------------------------------------------
         
        projection(i) = projection(i)+proj;
    end
    
    projection = reshape(projection,size(r));
end


%=========================================================================
%=========================================================================