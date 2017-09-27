%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mod by H Urakubo                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%
%%%
function  main_distance4
	clear;
	Addpaths();
%%
%%
%%
	dd  = '170706_2220';
	fprintf('%s\n',dd);
%	for i = 1:27;
%		if ismember(i,[9,10,11]);
%			continue;
%		end;
%		distance4(i, dd)
%	end;

distance4(18, dd)

%%%
function  distance4(NUM, dd)
%%%
%%%
%%%

	OutputDir   = sprintf('./DATA/%02g',NUM);
	ddd  = dd;
	grey = [0.6,0.6,0.6];


%%% -27, 67
%%% -16, 80


%%%
%%% Load data and save output files
%%%
	fprintf('./%s/%sFinal.mat',OutputDir,ddd);
	fname    = sprintf('./%s/%sFinal.mat',OutputDir,ddd);
	load(fname);
	FILENAME = sprintf('./%s/%s.inp',OutputDir,dd);
    [inode, ihedra, iface] = saveabaqus_mod2(node, PSDface, elem, FILENAME);

%%%
%%% Dilution of areas of SpineHead and SpineNeck
%%%
	Radius      = 3;
	SpineHeadI  = logical(SpineHeadI > 0.5);
	SpineHeadI2 = imdilate(SpineHeadI, strel('sphere', Radius));
	SpineNeckI  = logical(SpineNeckI > 0.5);
	SpineNeckI2 = imdilate(SpineNeckI, strel('sphere', Radius));


%%%
%%% Detect lines crossing a target spine
%%%

% Select branches across the spine.
% Obtain a single center line
% Obtain target segment

	skel = SelectBranchCrossingDomain(skel, SpineI);
	fprintf('Number of branches for a centerline: %g\n', numel(skel));
	CenterLine = [];
	for i=1:length(skel);
		CenterLine = [CenterLine; flipud(skel{i})];
	end;
	if (length(skel) == 1);
		CenterLine = skel{1};
	elseif (length(skel) == 2)
		CenterLine =  [skel{1}(1:end-1,:); flipud(skel{2})] % 
	elseif (length(skel) == 3) % for No. 18
		CenterLine =  [skel{1}(1:end-1,:); skel{2}] % 
		% fprintf('\n\n');
		% CenterLine =  skel{1}
		% fprintf('Please manually set the centerline. \n');
		
		figure;
		plot3(	skel{1}(:,2)*xypitch, skel{1}(:,1)*xypitch, skel{1}(:,3)*xypitch, ...
			'-','Color','b','LineWidth',4);
		hold on;
		plot3(	skel{2}(:,2)*xypitch, skel{2}(:,1)*xypitch, skel{2}(:,3)*xypitch, ...
			'-','Color','g','LineWidth',4);
		plot3(	CenterLine(:,2)*xypitch, CenterLine(:,1)*xypitch, CenterLine(:,3)*xypitch, ...
			'-','Color','r','LineWidth',2);

	else		
		fprintf('Please manually set the centerline. \n');
		return;
	end;

	CenterLine = SelectSegmentCrossingDomain(CenterLine, SpineI);


%%%
%%% "PSD-near end" will be set as the first index of center line.
%%%

	BC_PSD = SurfTriBaryCenter(node, PSDface);

	SS     = norm(CenterLine(1,:)   - BC_PSD);
	SG     = norm(CenterLine(end,:) - BC_PSD);
	if (SG < SS)
		CenterLine = flipud(CenterLine);
	end;


%%%
%%% Graph plot for confirmation
%%%

	GraphPlot1(YY,XX,ZZ,SpineHeadI,SpineNeckI,DendriteI,grey,az,el, CenterLine, xypitch);

	tmpnode = [node(:,2)* xypitch, node(:,1)* xypitch,node(:,3)* xypitch];



%%%
%%% Obtain barycenters of tetrahedra
%%%	

	inode     = inode(:,2:4);
	barycenter_hedra = TetBaryCenter(inode, ihedra{1});

%%%
%%%  Tetrahedra are assigned to Domains 'SpineHead' or 'SpineNeck'
%%%

	ihedra_id = ihedra{1}(:,1);

	Xid = ceil(barycenter_hedra(1,:));
	Yid = ceil(barycenter_hedra(2,:));
	Zid = ceil(barycenter_hedra(3,:));

	tetra_domain = ones(numel(ihedra_id), 1);
	for i = 1:numel(ihedra_id);
		if( SpineHeadI2(Xid(i),Yid(i), Zid(i)) > 0.5 );...
				tetra_domain(i,1) = 0.2;
		elseif( SpineNeckI2(Xid(i),Yid(i), Zid(i)) > 0.5 );
				tetra_domain(i,1) = 0.5;
		end;
	end;

	Head_ID     = find(tetra_domain(:,1) == 0.2);
	Head_Bc     = barycenter_hedra(:,Head_ID);
	Neck_ID     = find(tetra_domain(:,1) == 0.5);
	Neck_Bc     = barycenter_hedra(:,Neck_ID);
	Dendrite_ID = find(tetra_domain(:,1) == 1);
	Dendrite_Bc = barycenter_hedra(:,Dendrite_ID);

%%%
%%% 170616 CHECK
%%%
%	figure;
%	plotmesh(tmpnode,elem(Neck_ID,:));
%%%
%%%

%%%
%%% Location along the centerline
%%%


	LengthCenterLine = TotalDistance(CenterLine);
	
	HeadLoc = Derive_1D_Loc(CenterLine, Head_ID, barycenter_hedra);
	NeckLoc = Derive_1D_Loc(CenterLine, Neck_ID, barycenter_hedra);

%	figure;
%	plot(NeckLoc);

	HeadLoc = HeadLoc * xypitch;
	NeckLoc = NeckLoc * xypitch;


%%%
%%% PSD location along the centerline
%%%

	barycenter_PSDs = TriBaryCenter(inode, iface{1});
	PSDLoc          = Derive_1D_Loc(CenterLine, [1:numel(barycenter_PSDs(1,:))], barycenter_PSDs);
	Face_ID         = iface{1}(:,1);
	Face_Trg        = iface{1}(:,2:4);

	PSD_Area = [];
	for i = 1:numel(Face_Trg(:,1))
		tmp = TriArea(inode, Face_Trg(i,:));
		PSD_Area = [PSD_Area; tmp];
	end;

	PSDLoc   = PSDLoc   .* xypitch;
	PSD_Area = PSD_Area .* xypitch .* xypitch;


%%%
%%%  Obtain section area along the center line.
%%%

	dist_area = ObtainSectionAreaAlongLine(CenterLine, SpineI, xypitch);

%	figure;
%	plot(dist_area(:,1), dist_area(:,2), 'o-');
%	xlabel('Distance from tip (um)');
%	ylabel('area (um2)');
	

%%%
%%% Plot location
%%%
% 	GraphPlot2(YY,XX,ZZ, Head_Bc,Neck_Bc,Dendrite_Bc,...
% 		Head_ID, Neck_ID, HeadLoc, NeckLoc, grey, az, el, CenterLine, xypitch, cnum);

%%%
%%%
%%% Save files
%%%
%%%

	%%
	FILENAME = sprintf('./%s/%sHead_ID.dat',OutputDir,dd);
    dlmwrite(FILENAME, Head_ID,'precision','%g','delimiter',' ');
	FILENAME = sprintf('./%s/%sNeck_ID.dat',OutputDir,dd);
    dlmwrite(FILENAME, Neck_ID,'precision','%g','delimiter',' ');
    FILENAME = sprintf('./%s/%sDendrite_ID.dat',OutputDir,dd);
    dlmwrite(FILENAME, Dendrite_ID,'precision','%g','delimiter',' ');
	FILENAME = sprintf('./%s/%sHead_Loc.dat',OutputDir,dd);
    dlmwrite(FILENAME, HeadLoc,'precision','%g','delimiter',' ');
	FILENAME = sprintf('./%s/%sNeck_Loc.dat',OutputDir,dd);
    dlmwrite(FILENAME, NeckLoc,'precision','%g','delimiter',' ');
	FILENAME = sprintf('./%s/%sDist_Area.dat',OutputDir,dd);
	dlmwrite(FILENAME, dist_area,'precision','%g','delimiter',' ');
	%%
	FILENAME = sprintf('./%s/%sFace_ID.dat',OutputDir,dd);
    dlmwrite(FILENAME, Face_ID,'precision','%g','delimiter',' ');
	FILENAME = sprintf('./%s/%sFace_Loc_Area.dat',OutputDir,dd);
    dlmwrite(FILENAME, [PSDLoc, PSD_Area],'precision','%e','delimiter',' ');
	%%

%%%
%%%
%%% Functions
%%%
%%%


function Addpaths()

	addpath('C:\Users\urakubo\Desktop\Prog_Okabe_Laximi_Mizutani\iso2mesh');
 	addpath('C:\Users\urakubo\Desktop\Prog_Okabe_Laximi_Mizutani\STLRead');
 	addpath('C:\Users\urakubo\Desktop\Prog_Okabe_Laximi_Mizutani\FastMarching_version3b');
	addpath('C:\Users\urakubo\Desktop\Prog_Okabe_Laximi_Mizutani\FastMarching_version3b\functions');
	addpath('C:\Users\urakubo\Desktop\Prog_Okabe_Laximi_Mizutani\FastMarching_version3b\shortestpath');
	addpath('C:\Users\urakubo\Desktop\Prog_Okabe_Laximi_Mizutani\distancePointLine');
	addpath('C:\Users\urakubo\Desktop\Prog_Okabe_Laximi_Mizutani\affine3d');

	set(groot,'defaultAxesFontName',    'Arial');
	set(groot,'defaultTextFontName',    'Arial');
	set(groot,'defaultLegendFontName',  'Arial');
	set(groot,'defaultColorbarFontName','Arial');

function TotalDist = TotalDistance(skeleton)

	tmpd1 = skeleton;
	tmpd2 = skeleton;
	tmpd1(1,:)    = [];
	tmpd2(end,:)  = [];
	TotD          = tmpd1 - tmpd2;
	TotalDist     = 0;
	for i = 1:numel(TotD(:,1));
		tmp = sqrt(TotD(i,:)*TotD(i,:)');
		TotalDist = TotalDist + tmp;
	end;


function Target_BC_Loc = Derive_1D_Loc(skeleton, Target_BC_IDs, barycenter_hedra);
	
	Target_BC_Loc = zeros(numel(Target_BC_IDs),1);
	for i = 1:numel(Target_BC_IDs);
		TargetBC_XYZ = barycenter_hedra(:,Target_BC_IDs(i))';
		Dist = [];
		CloseP = []; 
		for j = 2:numel(skeleton(:,1)); %%
			[D, C, t0] = distancePoint2Line(skeleton(j-1,:), skeleton(j,:), TargetBC_XYZ, 'segment');
			Dist = [Dist;D];
			CloseP = [CloseP;C];
		end;
		[M,j] = min(Dist);
		% fprintf('ID j: %g; ', j);
		PartSkeleton = [skeleton(1:j,:); CloseP(j,:)];
		Target_BC_Loc(i,1) = TotalDistance(PartSkeleton);
	end;


function MM = Translation_matrix(CenterAxis);
	CenterAxis = CenterAxis'./norm(CenterAxis);
	P1_CenterAxis = [CenterAxis(2,1); -CenterAxis(1,1); CenterAxis(3,1)];
	M = makehgtform('axisrotate',CenterAxis,3.14159/2);
	M = M(1:3,1:3);
	P2_CenterAxis = M * P1_CenterAxis;
	MM = [P1_CenterAxis, P2_CenterAxis, CenterAxis];


function PSDcenter = SurfTriBaryCenter(node,PSDface);
	PSDcenter = [];
	for i = 1:numel(PSDface(:,1));
		tmp = node(PSDface(i,1),:) + node(PSDface(i,2),:) + node(PSDface(i,3), :);
		PSDcenter = [PSDcenter; tmp/3];
	end;
	PSDcenter = sum(PSDcenter)/length(PSDcenter(:,1));
	PSDcenter = PSDcenter(:,1:3);


function CenterL = SelectSegmentCrossingDomain(CenterLine, SpineI);
	L = CenterLine;
	len = numel(L(:,1));
	flag = zeros(len,1);
	for j=2:len;
		SS = SpineI(floor(L(j-1,1)), floor(L(j-1,2)), floor(L(j-1,3)));
		SG = SpineI(floor(L(j,1)), floor(L(j,2)), floor(L(j,3)));
		if (SS > 0.5)|(SG > 0.5);
			flag(j,1) = 1;
		end;
	end;
	CenterL = CenterLine(find(flag),:); %% This may cause trouble


function skel = SelectBranchCrossingDomain(skel, SpineI);
	len = numel(skel);
	flag = zeros(len,1);
	for i=1:len;
		L=skel{i};
		for j = 1:numel(L(:,1))
			SS = SpineI(floor(L(j,1)), floor(L(j,2)), floor(L(j,3)));
			if (SS > 0);
				flag(i,1) = 1;
				break;
			end;
		end;
	end;
	flag'
	skel = skel(find(flag));


function dist_area = ObtainSectionAreaAlongLine(CenterLine, SpineI, xypitch)

	snum = numel(CenterLine(:,1));

	[x1,y1,z1] = size(SpineI);
	xh = floor(x1/2);
	yh = floor(y1/2);
	zh = floor(z1/4);
	xst = (1-xh):(x1-xh);
	yst = (1-yh):(y1-yh);

	Stacked_Area   = [];
	dist_area = [];

	for j = 2:snum;
		CenterAxis = CenterLine(j,:) - CenterLine(j-1,:);
		MM = Translation_matrix(CenterAxis);
		bc_CA = (CenterLine(j,:) + CenterLine(j-1,:))./ 2;
		tmp = [MM,bc_CA';0,0,0,1];
		tmp_img = affine3(SpineI, tmp,xst,yst,1);
		tmp_img = logical(tmp_img > 0.5);
		tmpArea = bwarea(tmp_img) * xypitch * xypitch;
		Stacked_Area = cat(3, Stacked_Area, tmp_img);
		
		PartLine = [CenterLine(1:j-1,:); bc_CA];
		tmpDist  = TotalDistance(PartLine) * xypitch;
		dist_area = [dist_area; tmpDist, tmpArea];
	end;

	[x1,y1,z1] = size(Stacked_Area);
	zst = 1:z1;
	% figure
	% fv2   = isosurface(yst,xst,zst,Stacked_Area,0.5);
	% p2    = patch(fv2,'FaceColor',[0.3,0.3,1],'EdgeColor','none','FaceAlpha',.5);


function barycenter_hedra = TetBaryCenter(inode, ihedra);

	ihedra_id = ihedra(:,1);
	barycenter_hedra = zeros(3, numel(ihedra_id));
	for i = 1:numel(ihedra_id);
		tmp = inode(ihedra(i,2),:) + inode(ihedra(i,3),:) ...
			+ inode(ihedra(i,4),:) + inode(ihedra(i,5),:);
		barycenter_hedra(:,i) = tmp/4 ;
	end;


function barycenter_triangles = TriBaryCenter(inode, iface);

	iface_id = iface(:,1);
	barycenter_triangles = zeros(3, numel(iface_id));
	for i = 1:numel(iface_id);
		tmp = inode(iface(i,2),:) + inode(iface(i,3),:) + inode(iface(i,4),:);
		barycenter_triangles(:,i) = tmp/3 ;
	end;

function Area = TriArea(inode, Face_Trg)
	A = inode(Face_Trg(1),:);
	B = inode(Face_Trg(2),:);
	C = inode(Face_Trg(3),:);
	x = [A(1) B(1) C(1)];
	y = [A(2) B(2) C(2)];
	z = [A(3) B(3) C(3)];
	ons = [1 1 1];
	Area = 0.5*sqrt(det([x;y;ons])^2 + det([y;z;ons])^2 + det([z;x;ons])^2);


function GraphPlot1(YY,XX,ZZ,SpineHeadI,SpineNeckI,DendriteI,grey,az,el, CenterLine, xypitch);

	fig = figure;
	view(az,el);
	fv1   = isosurface(YY,XX,ZZ,DendriteI,0.5);
	p1    = patch(fv1,'FaceColor',grey,'EdgeColor','none','FaceAlpha',.5);
	hold on;
	fv2   = isosurface(YY,XX,ZZ,SpineNeckI,0.5);
	p2    = patch(fv2,'FaceColor',[0.3,0.3,1],'EdgeColor','none','FaceAlpha',.5);
	fv3   = isosurface(YY,XX,ZZ,SpineHeadI,0.5);
	p3    = patch(fv3,'FaceColor','b','EdgeColor','none','FaceAlpha',.5);
%	fv4   = isosurface(YY,XX,ZZ,PSDI,0.5);
%	p3    = patch(fv4,'FaceColor','r','EdgeColor','none','FaceAlpha',.5);
	hold on;
	plot3(	CenterLine(:,2)*xypitch, ...
			CenterLine(:,1)*xypitch, ...
			CenterLine(:,3)*xypitch, ...
			'-','Color','r','LineWidth',4);
	view(az,el);
	axis([4.5, 8, 4, 6, 0, 1.5]);
	axis equal;



function GraphPlot2(YY,XX,ZZ,Head_Bc,Neck_Bc,Dendrite_Bc,...
 	Head_ID, Neck_ID, HeadLoc, NeckLoc, ...
 	grey,az,el, CenterLine, xypitch, cnum);

	MarkerSize = 4;
	Opacity = 0.3;
	figure;
	view(az,el);
	cmap = colormap(jet(cnum));
	plot3(CenterLine(:,2)* xypitch,CenterLine(:,1)* xypitch,CenterLine(:,3)* xypitch,...
			'-','Color','r','LineWidth',4);
	hold on;
	size(Dendrite_Bc)
	s1 = scatter3(	Dendrite_Bc(2,:)* xypitch,...
				Dendrite_Bc(1,:)* xypitch,...
				Dendrite_Bc(3,:)* xypitch,...
		3,'MarkerEdgeColor','none','MarkerFaceColor',grey);
	%alpha(s1,Opacity)
	for i = 1:numel(Head_ID);
		s2 = scatter3(Head_Bc(2,i)* xypitch,Head_Bc(1,i)* xypitch,Head_Bc(3,i)* xypitch, MarkerSize,...
			'MarkerEdgeColor','none','MarkerFaceColor',cmap(1+ceil(abs(HeadLoc(i))),:) );
		%alpha(s2,Opacity);
	end;
	
	size(Neck_ID)
	size(Neck_Bc)
	
	for i = 1:numel(Neck_ID);
		s3 = scatter3(Neck_Bc(2,i)* xypitch,Neck_Bc(1,i)* xypitch,Neck_Bc(3,i)* xypitch, MarkerSize,...
			'MarkerEdgeColor','none','MarkerFaceColor',cmap(1+ceil(NeckLoc(i)),:) );
		%alpha(s3,Opacity);
	end;
	view(az,el);

