function im = runCamera()
%function [xaxreal, yaxreal, imzoomNormalized] = runCamera()
global fittype
global toFWHM

% add .mex functions for camera to path
if ~isdeployed
    addpath('matcam');
end

% initialize camera, check if it's closed first
try
    CloseCamera(uint32(1));
end
[id, wi, hi] = LoadCamera();
ccdWi = double(wi);
ccdHi = double(hi);
pix = 5.2;      % pixel size in µm
magn = 66;      % (estimated) magnification
dx = pix/magn;
dy = dx;
airyInt = 0.302593;         %Airy disk of NA=0.8 objective integrated over full surface

% get camera parameters
maxgain = GetMaxGain(id);
gain = GetGain(id);
expset = GetExposureInfo(id);
gboost = SetGainBoost(id,-1);

% set exposure to some reasonable value
myExposure = 0.7;
[~] = SetExposure(id, myExposure);
expset(1) = myExposure;

% get display parameters
dum = get(0, 'ScreenSize');
scrWi = dum(3);
scrHi = dum(4);

% initialize figure
fiWi = 0.8*scrWi;
fiHi = 0.8*scrHi;
fh = figure('Position', [(scrWi-fiWi)/2, (scrHi-fiHi)/2, fiWi, fiHi]);
im = GetImage(id, wi, hi);
im = rot90(reshape(im, wi, hi));
ih = imagesc(im);
ah = gca;
ah.Position(1) = 0.07;
ah.Position(3) = ah.Position(4) * scrHi/scrWi * ccdWi/ccdHi;
status = title('Free running');
colormap(gray);
ah.YDir = 'normal';

xlim([0, 1279]);
ylim([0, 1027]);

% ui controls
% exposure time
expth = uicontrol('Style','text','Position',[80 30 370 15],'String',['Exposure Time: ' num2str(expset(1)) ' ms']);
set(expth,'BackgroundColor',get(fh,'Color'));
exph = uicontrol('Style','slider','Position',[80 50 370 15]);
minexpstep = expset(4)/(expset(3) - expset(2));
maxexpstep = 10*minexpstep;
exph.Callback = {@ThorCam_ExpCallback,expth,expset(4),id};
exph.Min = expset(2);
exph.Max = expset(3);
exph.Value = expset(1);
exph.SliderStep = [minexpstep maxexpstep];

% Zoom axes
zoomwin = 0.55;
posaxzoom = [ah.Position(1)+ah.Position(3)+0.04, ah.Position(2), ...
    zoomwin*scrHi/scrWi, zoomwin];
axzoom = axes('Position', posaxzoom);
box on;
% Zoom parameters
region = 30;

% Fitted Gaussian axes
posaxave = [posaxzoom(1), posaxzoom(2)+posaxzoom(4)+0.05, 0.215*scrHi/scrWi, 0.215];
fitax = axes('Position', posaxave);
box on;
annotsize = [0.8 0.74 0.2 0.15];
ant = annotation('textbox', annotsize, 'string', '', ...
    'EdgeColor', 'none', 'Color', 'k', 'FontSize', 40, ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle');

% toggle for PSF measurement
pospsfgo = [1000 25 100 50];
fitgo = uicontrol('Style', 'togglebutton', ...
    'Position', pospsfgo, 'String', 'FIT GO');

% toggle for zoom
poszoombut = [800 25 100 50];
zoombut = uicontrol('Style', 'togglebutton', ...
    'Position', poszoombut, 'String', 'Zoom Enhance!');

% radio dials for fit type
buttongroupPosition = [0.32 0.03 0.1 0.05];
bg = uibuttongroup('Visible', 'off', ...
    'Position', buttongroupPosition, ...
    'SelectionChangeFcn', @buttonSelection);
r1 = uicontrol(bg, 'Style', 'radiobutton', ...
    'String', 'Gaussian', ...
    'Position', [0 0 100 50]);
r2 = uicontrol(bg, 'Style', 'radiobutton', ...
    'String', 'Airy', ...
    'Position', [100 0 100 50]);
bg.Visible = 'on';
% defaults here:
fittype = 'a * (2*besselj(1, sqrt((x-b1)^2/c^2 + (y-b2)^2/e^2)) / (sqrt((x-b1)^2/c^2 + (y-b2)^2/e^2)))^2 + d';
toFWHM = 2.355;

% exit button
posexit = [1200 25 100 50];
exitbut = uicontrol('Style', 'pushbutton', 'Position', posexit, ...
    'String', 'EXIT', 'Callback', {@close_down, id, fh});

% start loop
while ishandle(fh)
    try
        im = double(GetImage(id, wi, hi));
        im = rot90(reshape(im, wi, hi));
    catch
        % Sometimes the camera crashes; restart it like this.
        oldExposure = exph.Value;
        [id, wi, hi] = LoadCamera();
        [~] = SetExposure(id, oldExposure);
        continue
    end
    
    % Cut out the part around the brightest part of the image, if need be.
    if zoombut.Value == 1
        [~, b] = max(max(im));
        [maxval, a] = max(im(:, b));
        xaxpix = max(b-region, 1):min(b+region, wi);
        yaxpix = max(a-region, 1):min(a+region, hi);
        xaxreal = double((xaxpix - b)) * pix / magn;
        yaxreal = double((yaxpix - a)) * pix / magn;
        imzoom = im(yaxpix, xaxpix);
        axes(axzoom);
        imagesc(xaxreal, yaxreal, imzoom(:, :, 1), [0 255]);
        axzoom.YDir = 'normal';
        colormap(axzoom, 'jet');
        
        % Calculate Strehl ratio
        background = mean(mean(imzoom(1:10, 1:10, 1)));
        disp(maxval);
        imMinusBackground = imzoom(:, :, 1) - background;
        imzoomNormalized = airyInt * imMinusBackground ...
            / trapz(yaxreal, trapz(xaxreal, imMinusBackground, 2));
        strehlRat = max(max(imzoomNormalized));
        if fitgo.Value == 0
            ant.String = sprintf('SR = %.2f', strehlRat);
        end
    end
    
    ih.CData = im;
    drawnow;
    if ishandle(fh)
        if fitgo.Value == 1 && zoombut.Value == 1
            try
                status.String = 'Fitting';
                initconds2D = [maxval, 0.001, 0.001, 0.2, 10, 0.2];
                [Xo, Yo, Zo] = prepareSurfaceData(xaxreal, yaxreal, imzoom);
                fit2D = fit([Xo, Yo], Zo, fittype, 'StartPoint', initconds2D);
                axes(fitax);
                [X, Y] = meshgrid(xaxreal, yaxreal);
                imagesc(xaxreal, yaxreal, fit2D(X, Y));
                colormap(fitax, 'jet');
                
                %Alternatively calculate Strehl ratio from Airy fit:
                %strehlRat = fit2D.a*0.3026/(4*pi*fit2D.a*fit2D.c*fit2D.e);
                
                % Update fitted FWHMs in figure
                ant.String = [sprintf('FWHM = %.2f µm x %.2f µm\n', ...
                    [toFWHM*fit2D.c toFWHM*fit2D.e]), ...
                    sprintf('SR = %.2f', strehlRat)];
            catch
                status.String = 'Something went wrong, try again!';
                fitgo.Value = 0;
            end
        elseif fitgo.Value == 1 && zoombut.Value == 0
            zoombut.Value = 1;
        end
    end
end
end

function ThorCam_ExpCallback(src,~,expth,inc,id)
exptime = get(src,'Value');
factor = floor(exptime/inc);
exptime = factor*inc;
set(src,'Value',exptime);
set(expth,'String',['Exposure Time: ' num2str(exptime) ' ms']);
[~] = SetExposure(id,exptime);
end

function buttonSelection(~, event)
global fittype
global toFWHM
gaussEqn2D = 'a*exp( -(x - b1)^2 / (2*c^2) - (y - b2)^2 / (2*e^2)) + d';
airyEqn2D = 'a * (2*besselj(1, sqrt((x-b1)^2/c^2 + (y-b2)^2/e^2)) / (sqrt((x-b1)^2/c^2 + (y-b2)^2/e^2)))^2 + d';
if strcmp(event.NewValue.String, 'Gaussian')
    fittype = gaussEqn2D;
    toFWHM = 2.355;
elseif strcmp(event.NewValue.String, 'Airy')
    fittype = airyEqn2D;
    toFWHM = 2*1.6;
end
end

function close_down(~, ~, id, fh)
CloseCamera(id);
close(fh);
end