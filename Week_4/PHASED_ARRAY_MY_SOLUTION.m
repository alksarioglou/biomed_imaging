%=========================================================================
%                                                                     
%       BIOMEDICAL IMAGING
%       ULTRASOUND 2
%
%=========================================================================

%=========================================================================
%	PHASED ARRAY
%=========================================================================

function [] = PHASED_ARRAY_DEPTH_OF_5_CM()

    %clear all; close all; 
    
    fprintf ( '-----------------------------------------\n' );  
    fprintf ( ' PHASED ARRAY                            \n' );  
    fprintf ( '-----------------------------------------\n' );  
    
    
    % Ultrasound specs
    
    f = 1.5E6;                      % frequency [Hz]
    c = 1480;                       % speed of sound in water [m/s]
    k = 2*pi*f/c;                   % wave number [rad/m]
    lambda = c/f;                   % wavelength [m]
    
    
    % Transducer array
    
    nt    = 40;                         % number of transducers
    pitch = 0.5*lambda;                 % spacing of transducers
    
    % Definition of 2D grid for wave calculation
    
    n = 1024;                       % grid size in each dimension
    x_range = [0,0.4];              % covered range in x [m]
    y_range = [-0.2,0.2];           % covered range in y [m]
    
    x_values = x_range(1)+(x_range(2)-x_range(1))*[0:(n-1)]./(n-1);    % vector of x values on the grid
    y_values = y_range(1)+(y_range(2)-y_range(1))*[0:(n-1)]./(n-1);    % vector of y values on the grid
    
    [x,y] = meshgrid(x_values,y_values);    % the grid, given by 2D arrays x,y 
    
    % Gaussian bell function
      
    a = 1;
    b = 20;
    c = 9.32;
    
    % Loop through transducers and sum up their wave contributions
    
    wave = zeros(n,n);                      % container for pressure wave
    
    for t = 1:nt  
        
        x0 = -0.001;                        % x coordinate of transducer, 1 mm outside the container 
        y0 = -pitch*(t-(nt+1)/2);            % y coordinate of transducer
        
        r = sqrt((x-x0).^2+(y-y0).^2);      % distance from transducer, on the grid
        
        gaussian_amp = a*exp(-((t-b)^2)/(2*c^2)); % gaussian amplitude calculation
            
        phase = 0;                          % phase of transducer input [rad]
        amplitude = gaussian_amp;           % relative amplitude of transducer input          
        wave = wave + exp(1i*(k*r))./sqrt(r) .* amplitude .* exp(1i*phase);   % add partial wave to net wave
        
    end
    

    figure('position',[100 100 800 500])    
    imagesc(real(wave)), title('real part');  
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    
    figure('position',[100 100 800 500])
    imagesc(abs(wave)), title('magnitude')
    set(gca, 'XTick', []);
    set(gca, 'YTick', []);
    
    
end
    