%=========================================================================
%                                                                     
%	TITLE: 
%       DisplayData.m				
%								
%	DESCRIPTION:						
%	    Display data using subplot
%
%	INPUT:								
%       1D or 2D or 3D data array 		
%
%	OUTPUT:							
%       Display
%			
%	VERSION HISTORY:						
%	    120216SK INITIAL VERSION 
%	    191020SK UPDATE
%
%=========================================================================

%=========================================================================
%	M A I N  F U N C T I O N
%=========================================================================f
function DisplayData(data,subplotno,scale)

    numberimages = numel(data)/(size(data,1)*size(data,2));
    data = reshape(data,size(data,1),size(data,2),numberimages);
    imageposno = 1;
    
    if nargin>1 && numel(subplotno)>=2,
        gridsizex = subplotno(1);
        gridsizey = subplotno(2);
        if numel(subplotno)>=3,
            imageposno = subplotno(3);
        end
    else
        gridsizex = ceil(sqrt(numberimages));
        gridsizey = ceil(numberimages/gridsizex);
    end
    
    if nargin<3, 
        defaultScale = 1; 
    else
        defaultScale = 0;
        minScale = scale(1);
        maxScale = scale(2);
    end

    if isreal(data),
        RealData = 1;       % real data
        maxint = max(data(:));
        minint = min(data(:));
    elseif isreal(data*i),
        RealData = -1;      % imaginary data
        maxint = max(imag(data(:)));
        minint = min(imag(data(:))); 
    else
        RealData = 0;       % complex data
        maxint = max(abs(data(:)));
        minint = min(abs(data(:))); 
    end
    
    if ~defaultScale
        minint = minScale;
        maxint = maxScale;
    end

    if (maxint==minint),
        maxint=maxint+1;
        minint=minint-1;
    end

    for imageno=1:numberimages
        subplot(gridsizex,gridsizey,imageposno);
        
        if imageposno==1
            h = gcf;
            set(h,'Position',[500,500,size(data,1)*gridsizey,size(data,2)*gridsizex]);
        end
        
        imageposno = imageposno+1;
        
        if size(data,2)>10
            if RealData==1,         % real data
                imagesc(      data(:,:,imageno), [minint, maxint]);
            elseif RealData==-1,    % imaginary data
                imagesc(imag(data(:,:,imageno)), [minint, maxint]);
            else                    % complex data
                imagesc( abs(data(:,:,imageno)), [minint, maxint]);
            end
            graphtitle = sprintf('%d / %d',imageno,numberimages);
            title(graphtitle);
            axis xy;
            xlabel('X');
            ylabel('Y');
            colormap(gray);
        else
            plot(abs(data));
            axis tight; 
        end
    end
end


%=========================================================================
%=========================================================================     