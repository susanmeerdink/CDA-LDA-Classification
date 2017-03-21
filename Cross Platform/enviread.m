function I=envireadnostrip(varargin)
%Input:
%varargin{1} - input file name to be read in (this file can be with no
%extension or with a .img extension; the function will use this filename to
%look for the corresponding header file

%varargin{2 ...} - function options, the option type followed by any needed
%parameters for the option
%E.G. - I=enviread(inputFilename, 'lines', [1000:1500], 'samples', [2000:2500]);

%OPTIONS:
%lines - this option requires a parameter consisting of a range of lines
%(e.g. [1000:1500]; function will read in requested lines

%samples - this option requires a parameter consisting of a range of samples
%(e.g. [1000:1500]; function will read in requested samples


%Original version by Ian Howat, Ohio State Universtiy, ihowat@gmail.com
%Thanks to Yushin Ahn and Ray Jung

%Version change history by CSI Lab, University of Florida

%12/22/2010 - changed format for all input interleave types to
%strcat('*',format), this change insures that Matlab will store the data in
%the same filetype as it is read from...  RRC (rrclose@cise.ufl.edu)

%1/4/2011 - added ability for function to strip off .img extension to given
%file names for use in building the header file name... RRC

%1/4/2011 - added ability for function to read in only a specific set of
%lines and samples (in BIP format only)... RRC

file=varargin{1};
%[st1] = regexp(varargin{1}, '.img', 'split');
[st1] = {varargin{1}};

hdrfile=[deblank(st1{1}),'.hdr'];
info=read_envihdr(hdrfile);
info.hdrname=hdrfile;

LinesToProcess = 1:info.lines;
SamplesToProcess = 1:info.samples;
samplesOptionFlag = 0;

if nargin > 1
    processedArgs = 1;
    while processedArgs < nargin
        processedArgs = processedArgs + 1;
        switch lower(varargin{processedArgs})
            case {'lines'}
                processedArgs = processedArgs + 1;
                LinesToProcess = varargin{processedArgs};
                if strcmpi(info.interleave, 'bil') || strcmpi(info.interleave,  'bsq')
                    fprintf('BIL and BSQ format not supported for selecting lines!\n');
                    return;
                end
            case {'samples'}
                processedArgs = processedArgs + 1;
                SamplesToProcess = varargin{processedArgs};
                samplesOptionFlag = 1;
                if strcmpi(info.interleave, 'bil') || strcmpi(info.interleave,  'bsq')
                    fprintf('BIL and BSQ format not supported for selecting samples!\n');
                    return;
                end
            otherwise
                fprintf('Unknown option: %s\n', varargin{processedArgs});
        end
        
        
    end
    
end




% Make geo-location vectors
if isfield(info, 'map_info')
    if isfield(info.map_info,'mapx') && isfield(info.map_info,'mapy')
        xi = info.map_info.image_coords(1);
        yi = info.map_info.image_coords(2);
        xm = info.map_info.mapx;
        ym = info.map_info.mapy;
        %adjust points to corner (1.5,1.5)
        if yi > 1.5
            ym =  ym + ((yi*info.map_info.dy)-info.map_info.dy);
        end
        if xi > 1.5
            xm = xm - ((xi*info.map_info.dy)-info.map_info.dx);
        end
        
        I.x= xm + ((0:info.samples-1).*info.map_info.dx);
        I.y = ym - ((0:info.lines-1).*info.map_info.dy);
    end
end

% Set binary format parameters
switch info.byte_order
    case {0}
        machine = 'ieee-le';
    case {1}
        machine = 'ieee-be';
    otherwise
        machine = 'n';
end
switch info.data_type
    case {1}
        format = 'uint8';
        numBytes = 1;
    case {2}
        format= 'int16';
        numBytes = 2;
    case{3}
        format= 'int32';
        numBytes = 4;
    case {4}
        format= 'single';
        numBytes = 4;
    case {5}
        format= 'double';
        numBytes = 8;
    case {6}
        disp('>> Sorry, Complex (2x32 bits)data currently not supported');
        disp('>> Importing as double-precision instead');
        format= 'double';
        numBytes = 8;
    case {9}
        error('Sorry, double-precision complex (2x64 bits) data currently not supported');
    case {12}
        format= 'uint16';
        numBytes = 2;
    case {13}
        format= 'uint32';
        numBytes = 4;
    case {14}
        format= 'int64';
        numBytes = 8;
    case {15}
        format= 'uint64';
        numBytes = 8;
    otherwise
        error(['File type number: ',num2str(dtype),' not supported']);
end

% file read
% Version 2 code by Yushin Ahn - replaces resize calls with loops (except
% for BIP formats) to work on big arrays.
tmp=zeros(numel(LinesToProcess), numel(SamplesToProcess),info.bands,format);
fid=fopen(file,'r');

switch lower(info.interleave)
    
    case {'bsq'}
        % Format:
        % [Band 1]
        % R1: C1, C2, C3, ...
        % R2: C1, C2, C3, ...
        %  ...
        % RN: C1, C2, C3, ...
        %
        % [Band 2]
        %  ...
        % [Band N]
        
        for b=1:info.bands
            for i=1:info.lines
                t=fread(fid,info.samples,strcat('*',format));
                tmp(i,:,b)=t;
            end
        end
        
    case {'bil'}
        % Format:
        % [Row 1]
        % B1: C1, C2, C3, ...
        % B2: C1, C2, C3, ...
        %
        %  ...
        % [Row N]
        
        for i=1:info.lines
            for b=1:info.bands
                t=fread(fid,info.samples,strcat('*',format));
                tmp(i,:,b)=t;
            end
        end
        
    case {'bip'}
        
        % Row 1
        % C1: B1 B2 B3, ...
        % C2: B1 B2 B3, ...
        % ...
        % Row N
        %This section authored by Ray Jung, APL-Johns Hopkins
        %         Z = fread(fid,info.samples*info.lines*info.bands,format,0,machine);
        
        STATUS = fseek(fid, info.samples*(min(LinesToProcess) - 1)*info.bands*numBytes, 0);
        if STATUS == -1
            fprintf( 'Erorr with fseek\n');
        end
        
        if samplesOptionFlag
            
            STATUS = fseek(fid, (min(SamplesToProcess) - 1)*info.bands*numBytes, 0);
            if STATUS == -1
                fprintf( 'Erorr with fseek\n');
            end
            Z = zeros(numel(SamplesToProcess)*numel(LinesToProcess)*info.bands,1, format);
            counter = 1;
            for i=1:numel(LinesToProcess)
                startIndex = ((counter-1)*numel(SamplesToProcess)*1*info.bands + 1);
                endIndex = startIndex + numel(SamplesToProcess)*1*info.bands - 1;
                Z(startIndex:endIndex ,1) = fread(fid,numel(SamplesToProcess)*1*info.bands,strcat('*',format),0,machine);
                counter = counter + 1;
                if i ~= max(LinesToProcess)
                    STATUS = fseek(fid, (info.samples - max(SamplesToProcess) + min(SamplesToProcess) - 1)*info.bands*numBytes, 0);
                    if STATUS == -1
                        fprintf( 'Erorr with fseek\n');
                    end
                end
            end
        else
            
            Z = fread(fid,info.samples*numel(LinesToProcess)*info.bands,strcat('*',format),0,machine);
        end
       
        Z = reshape(Z, [info.bands, numel(SamplesToProcess), numel(LinesToProcess)]);
        
        for k=1:info.bands
            tmp(:,:,k) = squeeze(Z(k,:,:))';
        end
        
        if isfield(info, 'map_info')
            I.x = I.x(SamplesToProcess);
            I.y = I.y(LinesToProcess);
        end
end
fclose(fid);
I.z=tmp;
I.info =info;



% sub function
function info = read_envihdr(hdrfile)
% READ_ENVIHDR read and return ENVI image file header information.
%   INFO = READ_ENVIHDR('HDR_FILE') reads the ASCII ENVI-generated image
%   header file and returns all the information in a structure of
%   parameters.
%
%   Example:
%   >> info = read_envihdr('my_envi_image.hdr')
%   info =
%          description: [1x101 char]
%              samples: 658
%                lines: 749
%                bands: 3
%        header_offset: 0
%            file_type: 'ENVI Standard'
%            data_type: 4
%           interleave: 'bsq'
%          sensor_type: 'Unknown'
%           byte_order: 0
%             map_info: [1x1 struct]
%      projection_info: [1x102 char]
%     wavelength_units: 'Unknown'
%           pixel_size: [1x1 struct]
%           band_names: [1x154 char]
%
%   NOTE: This function is used by ENVIREAD to import data.
% Ian M. Howat, Applied Physics Lab, University of Washington
% ihowat@apl.washington.edu
% Version 1: 19-Jul-2007 00:50:57
fid = fopen(hdrfile);
while fid;
    line = fgetl(fid);
    if line == -1
        break
    else
        eqsn = findstr(line,'=');
        if ~isempty(eqsn)
            param = strtrim(line(1:eqsn-1));
            param(findstr(param,' ')) = '_';
            value = strtrim(line(eqsn+1:end));
            if isempty(str2num(value))
                if ~isempty(findstr(value,'{')) && isempty(findstr(value,'}'))
                    while isempty(findstr(value,'}'))
                        line = fgetl(fid);
                        value = [value,strtrim(line)];
                    end
                end
                eval(['info.',param,' = ''',value,''';'])
            else
                eval(['info.',param,' = ',value,';'])
            end
        end
    end
end
fclose(fid);

if isfield(info,'wavelength')
    [st1] = regexp(info.wavelength, ',', 'split');
    info.wavelength = zeros(info.bands,1);
    
    for i = 1:info.bands
        [st2] = regexp(st1{i}, '{', 'split');
        [st2] = regexp(st2{end}, '}', 'split');
        info.wavelength(i,1) = str2num(st2{1});
    end
end

if isfield(info,'map_info')
    line = info.map_info;
    line(line == '{' | line == '}') = [];
    line = strtrim(split(line,','));
    info.map_info = [];
    info.map_info.projection = line{1};
    info.map_info.image_coords = [str2num(line{2}),str2num(line{3})];
    info.map_info.mapx = str2num(line{4});
    info.map_info.mapy = str2num(line{5});
    info.map_info.dx  = str2num(line{6});
    info.map_info.dy  = str2num(line{7});
    if length(line) == 9
        info.map_info.datum  = line{8};
        info.map_info.units  = line{9}(7:end);
    elseif length(line) == 11
        info.map_info.zone  = str2num(line{8});
        info.map_info.hemi  = line{9};
        info.map_info.datum  = line{10};
        info.map_info.units  = line{11}(7:end);
    end
end

if isfield(info,'pixel_size')
    line = info.pixel_size;
    line(line == '{' | line == '}') = [];
    line = strtrim(split(line,','));
    info.pixel_size = [];
    info.pixel_size.x = str2num(line{1});
    info.pixel_size.y = str2num(line{2});
    info.pixel_size.units = line{3}(7:end);
end

%
function A = split(s,d)
%This function by Gerald Dalley (dalleyg@mit.edu), 2004
A = {};
while (length(s) > 0)
    [t,s] = strtok(s,d);
    A = {A{:}, t};
end


