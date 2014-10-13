function plot_video_frame_V1(y,Ys,g_gam,DS,fCutLP,ordLP,rho,fc,fs,fr,DSy,DSA,saveDir,saveName,printOutput)


T= length(y);
D = length(g_gam);

[ATar,YHWTar,YTar] = aud_mod_V1(y,g_gam,DS,fCutLP,ordLP,rho);
[A,YHW,Y] = aud_mod_V1(Ys(:,fr),g_gam,DS,fCutLP,ordLP,rho);

delta = 1e-6;
ATar(ATar<delta) = delta;
A(A<delta) = delta;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fonts
FontName = 'Arial';
FSsm = 10;
FSmed = 10;
FSlg = 12;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
LWy = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

left = 0.05;
vspace1 = 0.2;
top = 0.1;
bottom = 0.1;
width1 = 0.5;
right = 0.01;

height1 = (1-top-bottom-vspace1)/2;

hspace1 = 0.1;
width2 = 1-width1-hspace1-left-right;

height2 = height1;

pos11 = [left,1-top-height1,width1,height1];
pos21 = [left,bottom,width1,height1];

pos12 = [left+width1+hspace1,1-top-height2,width2,height1];
pos22 = [left+width1+hspace1,bottom,width2,height2];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tim = [1:T]/fs;
tlim = [min(tim),max(tim)];
ylim = max(abs([Ys(:)]))*[-1,1];
flim = [min(fc),max(fc)]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure1 = figure;

PP = [0,0,23.00,15.40]; %*** paper position in centimeters
PS = PP(end-1:end); % paper size in centimeters

set(figure1,'paperpositionmode','manual','paperposition', ...
        PP,'papersize',PS, 'paperunits','centimeters');

% So the figure is the same size on the screen as when it is printed:
pu = get(gcf,'PaperUnits');
pp = get(gcf,'PaperPosition');
set(gcf,'Units',pu,'Position',pp)

PR = PS(1)/PS(2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax11 = axes('position',pos11);
hold on
plot(tim(1:DSy:end),y(1:DSy:end),'-k','linewidth',LWy)
set(gca,'xlim',tlim,'ylim',ylim)
set(gca,'TickDir','out')
xlabel('time /s','fontname',FontName,'FontSize',FSlg)
title('original','fontname',FontName,'FontSize',FSlg)              
set(gca,'fontname',FontName,'FontSize',FSlg)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax21 = axes('position',pos21);
hold on
plot(tim(1:DSy:end),Ys(1:DSy:end,fr),'-k','linewidth',LWy)
set(gca,'xlim',tlim,'ylim',ylim)
xlabel('time /s','fontname',FontName,'FontSize',FSlg)
title('new signal','fontname',FontName,'FontSize',FSlg) 
set(gca,'TickDir','out')
set(gca,'fontname',FontName,'FontSize',FSlg)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax12 = axes('position',pos12);
hold on
tDS = [1:DSA:T];
surf(tim(tDS),[1:D],log(ATar(tDS,:)'),'edgecolor','none')
set(gca,'xlim',tlim,'ylim',[1,D]);
view(0,90)
xlabel('time /s','fontname',FontName,'FontSize',FSlg)
ylabel('freq /kHz','fontname',FontName,'FontSize',FSlg)              
set(gca,'ytick',[1:2:D])
set(gca,'yticklabel',num2str(round(fc([1:2:D])'/100)/10))
set(gca,'fontname',FontName,'FontSize',FSlg)
clim = get(gca,'clim');

%set(gca,'zlim',[0,])

title('target envelopes','fontname',FontName,'FontSize',FSlg)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ax22 = axes('position',pos22);
hold on
%imagesc(log(A'))
surf(tim(tDS),[1:D],log(A(tDS,:)'),'edgecolor','none')
set(gca,'xlim',tlim,'ylim',[1,D]);
view(0,90)
xlabel('time /s','fontname',FontName,'FontSize',FSlg)
ylabel('freq /kHz','fontname',FontName,'FontSize',FSlg)              
set(gca,'fontname',FontName,'FontSize',FSlg)
title('new envelopes','fontname',FontName,'FontSize',FSlg)
set(gca,'ytick',[1:2:D])
set(gca,'yticklabel',num2str(round(fc([1:2:D])'/100)/10))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
set(gca,'clim',clim)

saveName = [saveName,'_',num2str(fr)];

if printOutput==1
  print(figure1,'-depsc','-painters',[saveDir,saveName,'.eps'],'-loose');
  %plot2svg([filename,'ColBar','.svg'],figure2,'png')
  %str = ['! gv ',filename,'ColBar','.eps &'];
  %eval(str)
  plot2svg([saveDir,saveName,'.svg'],figure1);
end
