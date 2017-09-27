%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mod by H Urakubo                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
%%%
%%%
function  main_OneSpine8
	clear;
	Addpaths();
	dd  = char(datetime('now','Format','yyMMdd_HHmm'));
	dd  = '170706_2220';
	fprintf('%s\n',dd);
%	for i = 1:27;
%		if ismember(i,[9,10,11]);
%			continue;
%		end;
%		OneSpine(i, dd)
%	end;

OneSpine(18, dd)

%%%
%%%
%%%
function  OneSpine(NUM, dd)
%%%
%%%
%%%

	HeadName    = sprintf('d6sp%02ghead',NUM);
	NeckName    = sprintf('d6sp%02gneck',NUM);
	DendName    = 'Dendrite6';
	
	
	TargetDomain  = {HeadName, NeckName};
	OutputDir   = sprintf('./DATA/%02g',NUM);{HeadName, NeckName}

	fprintf('HeadName:  %s\n', HeadName);
	fprintf('NeckName:  %s\n', NeckName);
	fprintf('OutputDir: %s\n', OutputDir);

	MorphDir  = 'E:\170727Laxmi\CA1_HU_mod\Dendrite6';
	ddd = dd;
%%%
%%%

%%%
%%% Parameters
%%%
	az = -24;
	el = 66;
	maxvol = 100;


%% Coarse level 1.5
	opt       = struct( ...
		'radbound' , 4, ...    % Max tetrahedra size
		'distbound', 1);        % Distance from surface structure to volume mesh
	Radius    = [1 2];  %% Two domains
	RadiusPSD = 2;  %% Two domains


%% Coarse level 1
%	opt       = struct( ...
%		'radbound' , 3, ...    % Max tetrahedra size
%		'distbound', 1);        % Distance from surface structure to volume mesh
%	Radius    = [1 4];  %% Two domains
%	RadiusPSD = 2;  %% Two domains

	%% zpitch/xypitch must be an integer;
	zpitch  = 0.04; % 0.04
	xypitch = 0.013333333; 
	zmult   = round(zpitch/xypitch);


%%%
%%% Obtain 'dxf' filenames in a target directory.
%%%

	SortedFileNames  = ObtainDxfFilenames(MorphDir);
	SortedTargetFileNames = SpecifyTargetFiles(MorphDir, SortedFileNames,TargetDomain);

%%%
%%% Create canpus
%%%

%  	[minX, maxX, minY, maxY] = MinMax(MorphDir, SortedFileNames);
%	fprintf('Min X: %f\n', minX);
%	fprintf('Max X: %f\n', maxX);
%	fprintf('Min Y: %f\n', minY);
%	fprintf('Max Y: %f\n', maxY);
%
%	return;

	minX    = 5.0;
	maxX    = 5.2+9.1;

	minY    = 4.3;
	maxY    = 4.5+8.8;
	

	XX = 0:xypitch:(maxX-minX);
	YY = 0:xypitch:(maxY-minY);

	tmp = zeros(numel(XX),numel(YY));
	Canps = mat2gray(tmp);
	% size(tmp)
%%%
%%% Volume data aquition and Smoothing
%%%

	[SpineHeadI, ZZ] = VolumeDataAcquisition(MorphDir, SortedTargetFileNames, ...
		Canps, minX, minY, xypitch, zmult, TargetDomain{1});
	[SpineNeckI, ZZ] = VolumeDataAcquisition(MorphDir, SortedTargetFileNames, ...
		Canps, minX, minY, xypitch, zmult, TargetDomain{2});
	[DendriteI, ZZ] = VolumeDataAcquisition(MorphDir, SortedTargetFileNames, ...
		Canps, minX, minY, xypitch, zmult, DendName);
	[PSDI, ZZ] = VolumeDataAcquisition(MorphDir, SortedTargetFileNames, ...
		Canps, minX, minY, xypitch, zmult, 'PSD');
	[ERI, ZZ]  = VolumeDataAcquisition(MorphDir, SortedTargetFileNames, ...
		Canps, minX, minY, xypitch, zmult, 'ER');

	
	%%%
	%%% Extract connection between dend and neck,
	%%% and dend area extraction
	%%%

	SpineHeadI_dilate = imdilate(SpineHeadI, strel('sphere', 1));
	SpineNeckI_dilate = imdilate(SpineNeckI, strel('sphere', 1));
	DendriteI_dilate  = imdilate( DendriteI, strel('sphere', 1));
	tmp               = (SpineHeadI_dilate + SpineNeckI_dilate) .* DendriteI_dilate;
	tmp 	          = sum(tmp,1);
	tmp 	          = sum(tmp,2);
	Neck_Dend_Section = squeeze(tmp);
	Neck_Dend_Connect_ID = find(Neck_Dend_Section)
	Dend_ID_min = max(min(Neck_Dend_Connect_ID)-20,1)
	Dend_ID_max = min(max(Neck_Dend_Connect_ID)+20,numel(ZZ))


	%%%
	%%%
	DendriteI(:,:,Dend_ID_max:numel(ZZ)) = 0;
	DendriteI(:,:,1:Dend_ID_min)         = 0;
	%%%
	%%%


	%%%
	%%% Segment dendrite
	%%%


	StackedI = SpineHeadI + SpineNeckI + DendriteI;
	StackedI = (StackedI > 0.5);
	ERI      = (ERI > 0.5);
	StackedI = StackedSmoothing(StackedI, ERI, Radius); %% Radius;
	ExceptI  = SimpleSmoothing(PSDI,     RadiusPSD);


	SpineI   = SpineHeadI + SpineNeckI;
	SpineI   = (SpineI > 0.5);
	fname    = sprintf('%s/%sStep1.mat',OutputDir,dd);
	save(fname);
	fname    = sprintf('%s/%sStep1.mat',OutputDir,ddd);
	load(fname);

	PlotStacked(YY,XX,ZZ,StackedI, ExceptI, az,el);
%	PlotStacked(YY,XX,ZZ,ERI, ExceptI, az,el);
%	PlotStacked(YY,XX,ZZ,SpineHeadI, ExceptI, az,el);
%	return;
	
%%%
%%% Mesh the volume data
%%%

	[node,  elem,  face] = cgalv2m(StackedI, opt, maxvol);
	face  = unique(face(:,1:3),'rows');

%%%
%%% Obtain a center line
%%%

	bw   = logical(StackedI);
	skel = skeleton(bw);

	fname    = sprintf('%s/%sStep2.mat',OutputDir,dd);
	save(fname);
	fname    = sprintf('%s/%sStep2.mat',OutputDir,ddd);
	load(fname);

%%%
%%% Select PSD Surface
%%%
	nodeEL  = node; % x,y,z,?
	nodeEL  = ceil(nodeEL);
	numEL   = numel(nodeEL(:,1));
	PSDnode = zeros(numEL,1);
	for i = 1:numEL;
		PSDnode(i,1) = ExceptI(nodeEL(i,1),nodeEL(i,2),nodeEL(i,3));
	end;

	PSDface   = [];
	PSDfaceID = [];
	for i = 1:numel(face(:,1));
		tmp = PSDnode(face(i,1),1) + PSDnode(face(i,2),1)+ PSDnode(face(i,3),1);
		if (tmp > 1);
			PSDface = [PSDface; face(i,:)];
			PSDfaceID = [PSDfaceID, i];
		end;
	end;

%%%
%%% Save abaqus datas
%%%

	fname    = sprintf('./%s/%sFinal.mat',OutputDir,dd);
	save(fname);
	fname    = sprintf('./%s/%sFinal.mat',OutputDir,ddd);
	load(fname);

	FILENAME = sprintf('./%s/%s.inp',OutputDir,dd);
    saveabaqus(node, PSDface, elem, FILENAME);

%%%
%%% Plot and save
%%%

  	%% X and Y should be replaced ??
	tmpnode = [node(:,2), node(:,1),node(:,3)];
	figure;
%	plotmesh(node,elem,'x<400');
	plotmesh(tmpnode,face); 
	view(az,el);	
	hold on;
	plotmesh(tmpnode,PSDface,'facecolor','r');
	FILENAME = sprintf('./%s/%s.fig',OutputDir,dd);
	saveas(gcf,FILENAME);
	


%%%
%%% Plot a target spine and centerline
%%%

	figure;
	fv1   = isosurface(YY,XX,ZZ,SpineI,0.5);
	p1    = patch(fv1,'FaceColor','b','EdgeColor','none','FaceAlpha',.5);
	hold on;
	for i=1:length(skel);
		L=skel{i};
		plot3(L(:,2)*xypitch,L(:,1)*xypitch,L(:,3)*xypitch,...
			'-','Color',rand(1,3),'LineWidth',2);
	end;
	view(az,el);

function Addpaths()

	addpath('C:\Users\urakubo\Desktop\Prog_Okabe_Laximi_Mizutani\iso2mesh');
 	addpath('C:\Users\urakubo\Desktop\Prog_Okabe_Laximi_Mizutani\STLRead');
 	addpath('C:\Users\urakubo\Desktop\Prog_Okabe_Laximi_Mizutani\FastMarching_version3b');
	addpath('C:\Users\urakubo\Desktop\Prog_Okabe_Laximi_Mizutani\FastMarching_version3b\functions');
	addpath('C:\Users\urakubo\Desktop\Prog_Okabe_Laximi_Mizutani\FastMarching_version3b\shortestpath');
	addpath('C:\Users\urakubo\Desktop\Prog_Okabe_Laximi_Mizutani\distancePointLine');

%%%
%%% DataLoader
%%%

function [namepoly, polylines, polynum] = DataLoader(FileDir, FileName)

	%% Data load
	fname = sprintf('%s\\%s',FileDir,FileName);
	fid   = fopen(fname);
	C     = textscan(fid,'%s','Delimiter','\n');
	fclose(fid);
	C=C{1,1};

	%% ID & NUM of Polyline
	idpoly     = strcmp('POLYLINE',C);
	polynum    = sum(idpoly);
	idpoly     = find(idpoly == 1);
	namepoly   = C(idpoly+8);

	%% ID of end of Polyline
	idployend  = strcmp('SEQEND',C);
	idployend  = find(idployend == 1);

	%% ID of Vertices
	idvert     = strcmp('VERTEX',C);
	idvert     = find(idvert == 1);

	%% preallocate variable
	polylines  = cell(polynum,1);

	%%
	for i=1:polynum;
	%%
		clear id targid xpoly ypoly polyline;
    	id     = (idvert > idpoly(i)) .* (idvert < idployend(i));
		targid = idvert(find(id));
		xpoly  = str2double(C(targid+2));
		ypoly  = str2double(C(targid+4));
    	polyline     = [xpoly ypoly]; % create polylinesubset    
    	polylines{i} = polyline;
	%%
	end;
	%%


%%%
%%% Obtain Max Min
%%%

function [minX, maxX, minY, maxY] = MinMax(FileDir, SortedFileNames)

	minX = [];
	maxX = [];
	minY = [];
	maxY = [];

	for i = 1:numel(SortedFileNames);
		% fprintf('%s\n',char(SortedFileNames(i)));
		[namepoly, polylines, polynum] = ...
			DataLoader(FileDir, char(SortedFileNames(i)));
		for j = 1:polynum;
			tmpX = polylines{j}(:,1);
			tmpY = polylines{j}(:,2);
 			minX = min([minX; tmpX]);
 			maxX = max([maxX; tmpX]);
 			minY = min([minY; tmpY]);
 			maxY = max([maxY; tmpY]);
		end;
	end;


%%%
%%% Obtain 'dxf' filenames of a target directory.
%%%

function SortedFileNames = ObtainDxfFilenames(FileDir)

	tmp         = dir(sprintf('%s\\*.dxf',FileDir));
	FileNames   = {tmp.name};
	FileNo      = [];
	for i = 1:numel(FileNames);
		DecompFile = textscan(FileNames{i},'%s %d %s', 'Delimiter', '.');
		FileNo = [FileNo, DecompFile{end-1}];
	end;
	[FileNo, I] = sort(FileNo);
	SortedFileNames = FileNames(I);


%%%
%%% Specify target FileNames for Z stacks
%%%


function SortedTargetFileNames = ...
	SpecifyTargetFiles(FileDir, SortedFileNames,TargetDomain)

	%%%
	TargetFileNameID = [];
	for i = 1:numel(SortedFileNames);
	%%%
		%%% Load data
		[namepoly, polyl, polynum] = ...
			DataLoader(FileDir, char(SortedFileNames(i)));
		t = 0;
		for j = 1:polynum;
			for k = 1:numel(TargetDomain);
				if strcmp(char(namepoly{j}), char(TargetDomain{k}))
					t = t + 1;
				end;
			end;
		end;
		if (t > 0)
			TargetFileNameID = [TargetFileNameID,i];
		end;
	%%%
	end;
	%%%
	mmax = max(TargetFileNameID)+20;
	mmin = min(TargetFileNameID)-20;
	mmax = min(mmax, numel(SortedFileNames));
	mmin = max(mmin, 1);
	SortedTargetFileNames = SortedFileNames(mmin:mmax);
	%%%


%%%
%%% Volume data acquisition
%%%

function [TargetDomainI, ZZ] = ...
		VolumeDataAcquisition(FileDir, SortedTargetFileNames,...
			Canps, minX, minY, xypitch, zmult,TargetDomain);
	%%%
	Addpaths()
	%%%
	TargetDomainI = [];
	for j = 1:zmult*4;
		TargetDomainI = cat(3, TargetDomainI, Canps);
	end;
	%%%
	for i = 1:numel(SortedTargetFileNames);
	%%%
		%%% Load data
		[namepoly, polyl, polynum] = ...
			DataLoader(FileDir, char(SortedTargetFileNames(i)));
		%%% Fill holes and binarization
		Cyt   = Canps;
		for j = 1:polynum;
			if strcmp(char(namepoly{j}), TargetDomain) %% PSD
				tmp = roipoly(Cyt,(polyl{j}(:,1)-minX)/xypitch, (polyl{j}(:,2)-minY)/xypitch );
				Cyt = Cyt + tmp;
			end;
		end;
		%%% Resample
		Cyt = (Cyt > 0.5);
		for j = 1:zmult;
			TargetDomainI = cat(3, TargetDomainI, Cyt);
		end;
	%%%
	end;
	%%%
	for j = 1:zmult*4;
		TargetDomainI = cat(3, TargetDomainI, Canps);
	end;
	%%%
	ZZ = (0: numel(TargetDomainI(1,1,:))-1)*xypitch;
	%%%

%%%
%%% PlotStacked
%%%

function PlotStacked(YY,XX,ZZ,StackedI, ExceptI, az,el)
	figure;
	fv1   = isosurface(YY,XX,ZZ,StackedI,0.5);
	fv2   = isosurface(YY,XX,ZZ,StackedI,1.5);
	fvPSD = isosurface(YY,XX,ZZ,ExceptI,0.5);
	p1    = patch(fv1,'FaceColor','b','EdgeColor','none','FaceAlpha',.5);
	p2    = patch(fv2,'FaceColor','b','EdgeColor','none','FaceAlpha',.5);
	pPSD  = patch(fvPSD,'FaceColor','r','EdgeColor','none','FaceAlpha',.3);
	view(az,el);


%%%
%%% Smoothing
%%%

function StackedI = SimpleSmoothing(StackedI, Radius);

	Stacked1st = logical(StackedI > 0.5);
	Stacked1st = imdilate(Stacked1st, strel('sphere', Radius));
	Stacked1st = imerode(Stacked1st, strel('sphere', Radius));
	StackedI   = uint8(Stacked1st > 0.5);

function StackedI = StackedSmoothing(StackedI, ER,  Radius);


	%%%
	%%% Separate cytosol and SER
	%%%

	Stacked1st = logical(StackedI > 0.5);
	Stacked2nd = logical(ER > 0.5);
	
	%%%%
	%%%% Erosion of the 1st area
	%%%%
	Stacked1st =  imerode(Stacked1st, strel('sphere', Radius(1)));		

	%%%%
	%%%% Deletion of unconnected volumes
	%%%%
	CC = bwconncomp(Stacked1st);
	numPixels = cellfun(@numel,CC.PixelIdxList);
	[biggest,idx] = max(numPixels);
	Stacked1st = zeros(size(StackedI),'logical');
	Stacked1st(CC.PixelIdxList{idx}) = 1;
	
	%%%
	%%% Dilution of the 1st area
	%%%
	Stacked1st = imdilate(Stacked1st, strel('sphere', Radius(1)));

	
	%%%
	%%% Dilution of the 1st area
	%%%
	Stacked1st = imdilate(Stacked1st, strel('sphere', Radius(1)));

	%%%%
	%%%% Selection of Stacked2nd within Stacked1st
	%%%%
	Stacked2nd = Stacked2nd .* Stacked1st;

	%% The dilated 2nd stack (SER/Golgi) enlarges the area of 1st stack (cytosol).
	%% This is beacuse the cytosol must completely wrap such intracellular structures
	%% for keeping topology of surface mesh (each triangle must not have >3 neighbors).

	Stacked2nd_dil = imdilate(Stacked2nd, strel('sphere', Radius(2)));

	%%%%
	%%%% Final Rearrangement
	%%%%

	StackedI   = uint8((Stacked1st + Stacked2nd_dil) > 0.5) + uint8(Stacked2nd > 0.5);


