
	% ---------------------
	% optional contours

	with_contours = true;

	% projection = 'miller';
	projection = 'robinson';

	% ---------------
	% ---------------

	% ---------------

	% cm1 = turbo();
	cm1 = turbo;
	cm2 = flipud(plasma);
	% cm1_min = 0;
	% cm1_max = 2 * threshold; 	% threshold 1.18

	titleFontSize = 12;
	colorbarFontSize = 14;
	colorbarTickFontSize = 12;
	% contourFontSize = 9;
	contourFontSize = 7;
	% lonlatFontSize = 6;
	lonlatFontSize = 7;
	lab1FontSize=18;	% 24
	contours = [0:25:200];

	lontick = [-360:90:360];

	figure(1)
	hFig=gcf;
	clf(hFig);
	%

	tiledlayout(1,2,'TileSpacing','tight','Padding','tight');

	% [ (inner vert) (inner horiz) ]. [ (Bottom margin) (Top margin) ]. [ (L margin) (R margin) ] 
	% ha = tight_subplot(1,2,[.02 .03],[.12 .03],[.02 .02]);

	% - - -
	% a. etamax with dtalk = 1 umol/kg
	% - - -

	% ha(1) = nexttile(T);
	ha(1) = nexttile;
	% axes(ha(1));
	X = lon-360;
	Y = lat;
	% Z = dpct2020_htot_ssp126;
	Z = co2sys.etamax_2010_dtalk_1;
	%
	hold on

	if strcmp(projection,'robinson')
		m_proj('robinson','lon',[-340 20]);
	else
		m_proj('miller','lon',[-340 20],'lat',[-80 89.5]);
	end

	m_pcolor(X,Y,Z);
	m_coast('patch',[.7 .7 .7],'edgecolor','none');
	m_grid('tickdir','out','linewi',2, ...
		'xtick',lontick, ...
		'ytick',[-90 -60 -30 0 30 60 90], ...
		'linewidth',.5, ...
		'XaxisLocation','bottom','YaxisLocation','left','fontsize',lonlatFontSize);

	colormap(ha(1),cm1);
	% hcb1 = colorbar
	% hcb1.Label.String=(['NaOH treatment \DeltaTA \mumol kg^{-1}']);
	% caxis(cminmax);

	% xtext = 0.02;		% .13
	% ytext = 0.95;		% .84
	% text(xtext,ytext,'a','Units','Normalized','fontsize',lab1FontSize,'fontweight','bold')

	if with_contours
	[C,h] = m_contour(X,Y,Z,[10 25 50 75 100 150 200 250 300 400 500],'LineColor','k','Fill','off','ShowText','on','LineWidth',1);		% N Pacific
	clabel(C,h,'color','k','FontName','Arial Narrow','FontSize',contourFontSize);
	end

	% title({'a. OAE using NaOH'; '(unequilibrated eqn 17)'},'FontSize',titleFontSize)
	title('a._ \etamax in 2010 with \DeltaTA=1 \mumol kg^{-1}','FontSize',titleFontSize)
	c1 = colorbar(ha(1),'southoutside')		
	c1.Label.String = '\etamax';
	c1.Label.FontSize = colorbarFontSize;
	c1.FontSize = colorbarTickFontSize;
	hold off

	% - - -
	% b. delta_etamax using dtalk=100 umol/kg vs dtalk=1 umol/kg
	% - - -

	% ha(2) = nexttile(T);
	ha(2) = nexttile;
	% axes(ha(2));
	X = lon-360;
	Y = lat;
	% Z = dpct2020_htot_ssp126;
	Z = co2sys.etamax_2010_dtalk_100 - co2sys.etamax_2010_dtalk_1;
	%
	hold on

	if strcmp(projection,'robinson')
		m_proj('robinson','lon',[-340 20]);
	else
		m_proj('miller','lon',[-340 20],'lat',[-80 89.5]);
	end

	m_pcolor(X,Y,Z);
	m_coast('patch',[.7 .7 .7],'edgecolor','none');
	m_grid('tickdir','out','linewi',2, ...
		'xtick',lontick, ...
		'ytick',[-90 -60 -30 0 30 60 90], ...
		'linewidth',.5, ...
		'XaxisLocation','bottom','YaxisLocation','left','fontsize',lonlatFontSize);

	colormap(ha(2),cm2);
	% hcb1 = colorbar
	% hcb1.Label.String=(['Na2CO3 treatment \DeltaTA \mumol kg^{-1}']);
	% caxis(cminmax);

	% xtext = 0.02;		% .13
	% ytext = 0.95;		% .84
	% text(xtext,ytext,'a','Units','Normalized','fontsize',lab1FontSize,'fontweight','bold')

	if with_contours
	[C,h] = m_contour(X,Y,Z,[10 25 50 75 100 150 200 250 300 400 500],'LineColor','k','Fill','off','ShowText','on','LineWidth',1);		% N Pacific
	clabel(C,h,'color','k','FontName','Arial Narrow','FontSize',contourFontSize);
	end

	% title({'b. OAE using Na2CO3'; '(unequilibrated eqn 19)'},'FontSize',titleFontSize)
	title('b._ \Delta\etamax with \DeltaTA=100 vs \DeltaTA=1 \mumol kg^{-1}','FontSize',titleFontSize)
	c2 = colorbar(ha(2),'southoutside')		
	c2.Label.String = '\Delta\etamax';
	c2.Label.FontSize = colorbarFontSize;
	c2.FontSize = colorbarTickFontSize;
	hold off

	% % % - - -
	% % % Set colormap and color limits for all subplots
	% % % set(ha, 'Colormap', turbo, 'CLim', cminmax)
	% % set(ha, 'Colormap', plasma)
	% % % assign color bar to one tile 
	% % cbh = colorbar(ha(end)); 
	% % % To position the colorbar as a global colorbar representing
	% % % all tiles, 
	% % cbh.Layout.Tile = 'south'; 
	% % cbh.Label.String=(['OAE treatment \DeltaTA \mumol kg^{-1}']);
	% % cbh.FontSize=colorbarFontSize;

	% % - - -
	% % common color bar outside of subplots
	% h = axes(hFig,'visible','off'); 
	% h.Title.Visible = 'on';
	% h.XLabel.Visible = 'on';
	% h.YLabel.Visible = 'on';

	% % [left, bottom, width, height]
	% % c = colorbar(h,'southoutside','Position',[0.1 0.1 0.8 0.02]);  % attach colorbar to h
	% c = colorbar(h,'southoutside','Position',[0.1 0.09 0.8 0.015]);  % attach colorbar to h

	% colormap(c,cm1)

	% % caxis(h,[cm1_min cm1_max]);             % set colorbar limits
	% caxis(h,cminmax);             % set colorbar limits
	% % c.Label.String = 'Aragonite saturation state (\Omega_a_r)';
	% c.Label.String=(['OAE treatment required to restore TA-DIC in 2010 to pre-industrial (\DeltaTA \mumol kg^{-1})']);
	% c.Label.FontSize = colorbarFontSize;
	% % c.Ticks = round(linspace(cm1_min,cm1_max,9),2);
	% c.FontSize = colorbarTickFontSize;


	k=1.5;
	set(gcf, 'PaperPosition', [0 0 k*9 k*3])   
	print(gcf, [pwd '/png/global_map_2010_etamax_dTA1_and_delta_etamax_dTA100vs1_v21.png'], '-dpng', '-r600' );  


% - - -
% - - -
% - - -

% 5-step method

% - - -
% - - -
% - - -


	% ---------------------
	% optional contours

	with_contours = true;

	% projection = 'miller';
	projection = 'robinson';

	% ---------------
	% ---------------

	% ---------------

	% cm1 = turbo();
	cm1 = turbo;
	% cm1_min = 0;
	% cm1_max = 2 * threshold; 	% threshold 1.18

	titleFontSize = 12;
	colorbarFontSize = 14;
	colorbarTickFontSize = 14;
	% contourFontSize = 9;
	contourFontSize = 7;
	% lonlatFontSize = 6;
	lonlatFontSize = 7;
	lab1FontSize=18;	% 24
	contours = [0:25:200];

	Z1 = dTAtrt.unequil_2010_NaOH;
	% Z1(~in_CCE) = nan;
	Z2 = dTAtrt.cdreff80_2010_Na2CO3_step5;
	% Z2(~in_CCE) = nan;
	% cminmax = [min(nanmin2(Z1(:)),nanmin2(Z2(:))) max(nanmax2(Z1(:)),nanmax2(Z2(:)))];
	cminmax = [min(nanmin2(Z1(:)),nanmin2(Z2(:))) max(nanmax2(Z1(:)),prctile(Z2(:),99))];

	lontick = [-360:90:360];

	figure(1)
	hFig=gcf;
	clf(hFig);
	%

	% % T = tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
	% T = tiledlayout(2,2,'TileSpacing','Compact','Padding','Compact');
	% t_ab = tiledlayout(T,2,2,'TileSpacing','Compact','Padding','Compact');		% inner layout for a and b
	% t_ab.Layout.TileSpan = [2 2];

	% [ (inner vert) (inner horiz) ]. [ (Bottom margin) (Top margin) ]. [ (L margin) (R margin) ] 
	% ha = tight_subplot(2,2,[.03 .03],[.12 .03],[.02 .02]);
	% ha = tight_subplot(2,2,[.005 .03],[.12 .03],[.02 .02]);
	% ha = tight_subplot(2,2,[.015 .03],[.12 .03],[.02 .02]);
	ha = tight_subplot(2,2,[.02 .03],[.12 .03],[.02 .02]);

	% - - -
	% a. Un-equilibrated NaOH
	% - - -

	% ha(1) = nexttile(T);
	% ha(1) = nexttile(t_ab);
	axes(ha(1));
	X = lon-360;
	Y = lat;
	% Z = dpct2020_htot_ssp126;
	Z = dTAtrt.unequil_2010_NaOH;
	%
	hold on

	if strcmp(projection,'robinson')
		m_proj('robinson','lon',[-340 20]);
	else
		m_proj('miller','lon',[-340 20],'lat',[-80 89.5]);
	end

	m_pcolor(X,Y,Z);
	m_coast('patch',[.7 .7 .7],'edgecolor','none');
	m_grid('tickdir','out','linewi',2, ...
		'xtick',lontick, ...
		'ytick',[-90 -60 -30 0 30 60 90], ...
		'linewidth',.5, ...
		'XaxisLocation','bottom','YaxisLocation','left','fontsize',lonlatFontSize);

	colormap(ha(1),cm1);
	% hcb1 = colorbar
	% hcb1.Label.String=(['NaOH treatment \DeltaTA \mumol kg^{-1}']);
	caxis(cminmax);

	% xtext = 0.02;		% .13
	% ytext = 0.95;		% .84
	% text(xtext,ytext,'a','Units','Normalized','fontsize',lab1FontSize,'fontweight','bold')

	if with_contours
	[C,h] = m_contour(X,Y,Z,[10 25 50 75 100 150 200 250 300 400 500],'LineColor','k','Fill','off','ShowText','on','LineWidth',1);		% N Pacific
	clabel(C,h,'color','k','FontName','Arial Narrow','FontSize',contourFontSize);
	end

	% title({'a. OAE using NaOH'; '(unequilibrated eqn 17)'},'FontSize',titleFontSize)
	title('a._ Unequilibrated NaOH (0% CDR_{eff})','FontSize',titleFontSize)
	hold off

	% - - -
	% b. Un-equilibrated Na2CO3
	% - - -

	% ha(2) = nexttile(T);
	% ha(2) = nexttile(t_ab);
	axes(ha(2));
	X = lon-360;
	Y = lat;
	% Z = dpct2020_htot_ssp126;
	Z = dTAtrt.unequil_2010_Na2CO3;
	%
	hold on

	if strcmp(projection,'robinson')
		m_proj('robinson','lon',[-340 20]);
	else
		m_proj('miller','lon',[-340 20],'lat',[-80 89.5]);
	end

	m_pcolor(X,Y,Z);
	m_coast('patch',[.7 .7 .7],'edgecolor','none');
	m_grid('tickdir','out','linewi',2, ...
		'xtick',lontick, ...
		'ytick',[-90 -60 -30 0 30 60 90], ...
		'linewidth',.5, ...
		'XaxisLocation','bottom','YaxisLocation','left','fontsize',lonlatFontSize);

	colormap(ha(2),cm1);
	% hcb1 = colorbar
	% hcb1.Label.String=(['Na2CO3 treatment \DeltaTA \mumol kg^{-1}']);
	caxis(cminmax);

	% xtext = 0.02;		% .13
	% ytext = 0.95;		% .84
	% text(xtext,ytext,'a','Units','Normalized','fontsize',lab1FontSize,'fontweight','bold')

	if with_contours
	[C,h] = m_contour(X,Y,Z,[10 25 50 75 100 150 200 250 300 400 500],'LineColor','k','Fill','off','ShowText','on','LineWidth',1);		% N Pacific
	clabel(C,h,'color','k','FontName','Arial Narrow','FontSize',contourFontSize);
	end

	% title({'b. OAE using Na2CO3'; '(unequilibrated eqn 19)'},'FontSize',titleFontSize)
	title('b._ Unequilibrated Na_2CO_3 (0% CDR_{eff})','FontSize',titleFontSize)
	hold off


	% - - -
	% c. Equilibrated NaOH
	% - - -

	% ha(3) = nexttile(t_ab);
	axes(ha(3));
	X = lon-360;
	Y = lat;
	% Z = dpct2020_htot_ssp126;
	Z = dTAtrt.cdreff80_2010_NaOH_step5;
	%
	hold on

	if strcmp(projection,'robinson')
		m_proj('robinson','lon',[-340 20]);
	else
		m_proj('miller','lon',[-340 20],'lat',[-80 89.5]);
	end

	m_pcolor(X,Y,Z);
	m_coast('patch',[.7 .7 .7],'edgecolor','none');
	m_grid('tickdir','out','linewi',2, ...
		'xtick',lontick, ...
		'ytick',[-90 -60 -30 0 30 60 90], ...
		'linewidth',.5, ...
		'XaxisLocation','bottom','YaxisLocation','left','fontsize',lonlatFontSize);

	colormap(ha(3),cm1);
	% hcb3 = colorbar
	% hcb3.Label.String=(['NaOH treatment \DeltaTA \mumol kg^{-1}']);
	caxis(cminmax);

	% xtext = 0.02;		% .13
	% ytext = 0.95;		% .84
	% text(xtext,ytext,'a','Units','Normalized','fontsize',lab1FontSize,'fontweight','bold')

	if with_contours
	[C,h] = m_contour(X,Y,Z,[10 25 50 75 100 150 200 250 300 400 500],'LineColor','k','Fill','off','ShowText','on','LineWidth',1);		% N Pacific
	clabel(C,h,'color','k','FontName','Arial Narrow','FontSize',contourFontSize);
	end

	% title({'c. OAE using NaOH'; '(equilibrated eqn 11 and 22)'},'FontSize',titleFontSize)
	title('c._ Equilibrated NaOH (80% CDR_{eff})','FontSize',titleFontSize)
	hold off


	% - - -
	% d. Equilibrated Na2CO3
	% - - -

	% ha(4) = nexttile(t_ab);
	axes(ha(4));
	X = lon-360;
	Y = lat;
	% Z = dpct2020_htot_ssp126;
	Z = dTAtrt.cdreff80_2010_Na2CO3_step5;
	%
	hold on

	if strcmp(projection,'robinson')
		m_proj('robinson','lon',[-340 20]);
	else
		m_proj('miller','lon',[-340 20],'lat',[-80 89.5]);
	end

	m_pcolor(X,Y,Z);
	m_coast('patch',[.7 .7 .7],'edgecolor','none');
	m_grid('tickdir','out','linewi',2, ...
		'xtick',lontick, ...
		'ytick',[-90 -60 -30 0 30 60 90], ...
		'linewidth',.5, ...
		'XaxisLocation','bottom','YaxisLocation','left','fontsize',lonlatFontSize);

	colormap(ha(4),cm1);
	% hcb4 = colorbar
	% hcb4.Label.String=(['NaOH treatment \DeltaTA \mumol kg^{-1}']);
	caxis(cminmax);

	% xtext = 0.02;		% .13
	% ytext = 0.95;		% .84
	% text(xtext,ytext,'a','Units','Normalized','fontsize',lab1FontSize,'fontweight','bold')

	if with_contours
	[C,h] = m_contour(X,Y,Z,[10 25 50 75 100 150 200 250 300 400 500],'LineColor','k','Fill','off','ShowText','on','LineWidth',1);		% N Pacific
	clabel(C,h,'color','k','FontName','Arial Narrow','FontSize',contourFontSize);
	end

	% title({'d. OAE using Na2CO3'; '(equilibrated eqn 15 and 28)'},'FontSize',titleFontSize)
	title('d._ Equilibrated Na_2CO_3 (80% CDR_{eff})','FontSize',titleFontSize)
	hold off



	% % - - -
	% % Set colormap and color limits for all subplots
	% % set(ha, 'Colormap', turbo, 'CLim', cminmax)
	% set(ha, 'Colormap', plasma)
	% % assign color bar to one tile 
	% cbh = colorbar(ha(end)); 
	% % To position the colorbar as a global colorbar representing
	% % all tiles, 
	% cbh.Layout.Tile = 'south'; 
	% cbh.Label.String=(['OAE treatment \DeltaTA \mumol kg^{-1}']);
	% cbh.FontSize=colorbarFontSize;

	% - - -
	% common color bar outside of subplots
	h = axes(hFig,'visible','off'); 
	h.Title.Visible = 'on';
	h.XLabel.Visible = 'on';
	h.YLabel.Visible = 'on';

	% [left, bottom, width, height]
	% c = colorbar(h,'southoutside','Position',[0.1 0.1 0.8 0.02]);  % attach colorbar to h
	c = colorbar(h,'southoutside','Position',[0.1 0.09 0.8 0.015]);  % attach colorbar to h

	colormap(c,cm1)

	% caxis(h,[cm1_min cm1_max]);             % set colorbar limits
	caxis(h,cminmax);             % set colorbar limits
	% c.Label.String = 'Aragonite saturation state (\Omega_a_r)';
	c.Label.String=(['OAE treatment required to restore TA-DIC in 2010 to pre-industrial (\DeltaTA \mumol kg^{-1})']);
	c.Label.FontSize = colorbarFontSize;
	% c.Ticks = round(linspace(cm1_min,cm1_max,9),2);
	c.FontSize = colorbarTickFontSize;


	k=1.5;
	set(gcf, 'PaperPosition', [0 0 k*8 k*5])   
	print(gcf, [pwd '/png/global_map_dTA_to_restore_2010_alkstar_with_NaOH_and_Na2CO3_v21.png'], '-dpng', '-r600' );  



