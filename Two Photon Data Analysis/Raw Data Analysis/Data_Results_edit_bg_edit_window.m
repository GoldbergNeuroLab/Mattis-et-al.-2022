if exist(fullfile([folder_name_wr,'average_Image_ROIdata.mat']))
    load (fullfile([folder_name_wr,'average_Image_ROIdata.mat']));
else
    load (fullfile([folder_name_wr,'average_selected_Image_ROIdata.mat']));
end

% Frequency=30;
Interval=2;
% th_pk_amp=0.5;%if using change to 0.3
% mov_data = handles.mov_data;
% ParametersOutput = handles.parasaved;
Pixels=ParametersOutput.Pixels;
nFrames=size(M_rg,3);
Conditions=nFrames/500;
Cellnum=length(ParametersOutput.xypos);
F=zeros([nFrames Cellnum]);
vidHeight = size(M_rg,2);
vidWidth  = size(M_rg,1);
[bb,aa]=butter(3,[0.15]);
F=[];
for fi=1:nFrames
    M_rgc=M_rg(:,:,fi);
  for ci=1:Cellnum
      xi = Pixels{ci}(1,:);
      yi = Pixels{ci}(2,:);
      ind=sub2ind([vidHeight vidWidth],yi,xi);
      F(fi,ci)=mean(M_rgc(ind));
  end
 
end
F_filt=[];
for ci=1:Cellnum-1
     F_filt(:,ci)=filtfilt(bb,aa,F(:,ci));%-F(:,Cellnum)+mean(F(:,Cellnum))
end
F_reshaped=reshape(F_filt,[500,Conditions,Cellnum-1]);
F0=squeeze(mean(F_reshaped(5:50,:,:)));% changed from 5:15
% F0_std=squeeze(std(F_reshaped(5:50,:,:)));

y_labe=num2cell(0:5:Cellnum-1);
  for cci=1:length(y_labe)
     y_labe{cci}=num2str(y_labe{cci});
  end
  dF_data=[];
  dF_std=[];
  dF_mean=[];
  dF_data=[];
  count=[];
  val=[];
  loc=[];
  wid=[];
for cond=1:size(F_reshaped,2)
figure;
title(['Calcium traces for Condition' num2str(cond)])
hold on;
ymax=0;
ymin=0;
ind=find(cond_list==cond);%%%%%%included list if condition list is not in order
% count(cond)=0;
for ii=1:Cellnum-1
    F_temp=(F_reshaped(:,ind,ii)-F0(ind,ii))/F0(ind,ii);
    if isempty(strfind(folder_name_wr,'DW'))
    dF_std(cond,ii)=std(F_temp(:,:));%5:50 ,considered all first 400:end for DW
    dF_mean(cond,ii)=mean(F_temp(:,:));%5:50,considered all first 400:end for DW
    else
        dF_std(cond,ii)=std(F_temp(400:end,:));%5:50 ,considered all first 400:end for DW
        dF_mean(cond,ii)=mean(F_temp(400:end,:));%5:50,considered all first 400:end for DW
    end
%     plot((1:size(F_reshaped,1))/Frequency,F_temp+(ii)*Interval);
    F_temp(F_temp>100)=0;%%%% just to avoid noise at the ends
    plot((1:size(F_reshaped,1)),F_temp+(ii)*Interval);
    dF_data(:,cond,ii)=F_temp;
    ymin = min(min(F_temp+(ii)*Interval),ymin);
    ymax = max(max(F_temp+(ii)*Interval),ymax);
    
    [v,l,w,p]=findpeaks(dF_data(:,cond,ii),'WidthReference','halfheight');%,'MinPeakProminence',0.25);%,'MinPeakProminence',0.3);%0.3
%     if isempty(loc_flag)
    if isempty(pk_loc_mat)
        Pk_loc=200;
    else
        Pk_loc=pk_loc_mat(ind);
    end
    loc_t=l(find(l>=Pk_loc,1,'first'));
%     [v,l]=findpeaks(dF_data(:,cond,ii),'MinPeakWidth',5);%,'MinPeakHeight',th_pk_amp
%     if max(v)>=th_pk_amp & l(v==max(v))<300
    if (~isempty(loc_t) && loc_t>=Pk_loc && loc_t<=Pk_loc+15) && ((v(l==loc_t)-dF_mean(cond,ii))>3*dF_std(cond,ii))%50 300
        count(cond,ii)=1;
        val(cond,ii)=v(l==loc_t);
        loc(cond,ii)=loc_t;
        wid(cond,ii)=w(l==loc_t);
        hold on;plot(loc(cond,ii)-5,val(cond,ii)+(ii)*Interval+0.1,'k*','MarkerSize',4)
%         [val(cond,ii),loc(cond,ii)]=findpeaks(dF_data(:,cond,ii),'MinPeakWidth',15);
    else
        count(cond,ii)=0;
        val(cond,ii)=NaN;
        loc(cond,ii)=NaN;
        wid(cond,ii)=NaN;
    end
%     else
%         if l(v==max(v))>100 && l(v==max(v))<215 && ((max(v)-dF_mean(cond,ii))>3*dF_std(cond,ii))%50 300
%             count(cond,ii)=1;
%             val(cond,ii)=max(v);
%             loc(cond,ii)=l(v==max(v));
%             wid(cond,ii)=w(v==max(v));
%             hold on;plot(loc(cond,ii)-5,val(cond,ii)+(ii)*Interval+0.1,'k*','MarkerSize',4)
%             %         [val(cond,ii),loc(cond,ii)]=findpeaks(dF_data(:,cond,ii),'MinPeakWidth',15);
%         else
%             count(cond,ii)=0;
%             val(cond,ii)=NaN;
%             loc(cond,ii)=NaN;
%             wid(cond,ii)=NaN;
%         end
%     end
%     if l(p==max(p))>100 && l(p==max(p))<300 && v(p==max(p))>(dF_mean(cond,ii)+4*dF_std(cond,ii))
%         count(cond,ii)=1;
%         val(cond,ii)=v(p==max(p));
%         loc(cond,ii)=l(p==max(p));
%         wid(cond,ii)=w(p==max(p));
%         hold on;plot(loc(cond,ii)-5,val(cond,ii)+(ii)*Interval+0.1,'k*','MarkerSize',4)
% %         [val(cond,ii),loc(cond,ii)]=findpeaks(dF_data(:,cond,ii),'MinPeakWidth',15);
%     else
%         count(cond,ii)=0;
%         val(cond,ii)=NaN;
%         loc(cond,ii)=NaN;
%         wid(cond,ii)=NaN;
%     end
end

 xlabel('Frame Number');
 ylabel('Cell Number');
 xmin = 0;
 xmax = size(F_reshaped,1);
%  xmax = size(F_reshaped,1)/Frequency;
 ymin=0;
 axis([xmin xmax ymin ymax+Interval]);

 set(gca, 'YTick', [Interval Interval*5:Interval*5:Cellnum*Interval]);
 set(gca,'YTickLabel',y_labe);
 set(gca,'FontName','Times New Roman','FontSize',14);
 saveas(gca,fullfile([folder_name_wr,'_condition_' num2str(cond) '.tif']));
end
fname_xl=fullfile([folder_name_wr,'_results.xls']);
xlswrite(fname_xl,{'Condition';'Activated Cells';'Total Cells';'Peak Amplitude';'Peak Width'},'Summary','A1');
xlswrite(fname_xl,[1:Conditions],'Summary','B1');
xlswrite(fname_xl,[sum(count,2)'],'Summary','B2');
xlswrite(fname_xl,[repmat(Cellnum-1,1,Conditions)],'Summary','B3');
xlswrite(fname_xl,[nanmean(val,2)'],'Summary','B4');
xlswrite(fname_xl,nanmean(wid,2)','Summary','B5');

xlswrite(fname_xl,{'Condition'},'Peak Amplitude','B1');
xlswrite(fname_xl,{'Cell'},'Peak Amplitude','A2');
xlswrite(fname_xl,[1:Conditions],'Peak Amplitude','B2');
xlswrite(fname_xl,[0:Cellnum-2]','Peak Amplitude','A3');
xlswrite(fname_xl,[val]','Peak Amplitude','B3');

xlswrite(fname_xl,{'Condition'},'Activation','B1');
xlswrite(fname_xl,{'Cell'},'Activation','A2');
xlswrite(fname_xl,[1:Conditions],'Activation','B2');
xlswrite(fname_xl,[0:Cellnum-2]','Activation','A3');
xlswrite(fname_xl,[count]','Activation','B3');

    
% F_ave=mean(F);
% y_labe=num2cell(0:5:Cellnum);
%  
%  for cci=1:length(y_labe)
%      y_labe{cci}=num2str(y_labe{cci});
%  end
%  
% figure,
% title('Calcium traces')
% hold on;
% ymax=0;
% ymin=0;
% dF_data=zeros(size(F));
% for ii=1:Cellnum
%     F_temp=(F(:,ii)-F_ave(ii))/F_ave(ii);
%     plot((1:nFrames)/Frequency,F_temp+(ii)*Interval);
%     dF_data(:,ii)=F_temp;
%     ymin = min(min(F_temp+(ii)*Interval),ymin);
%     ymax = max(max(F_temp+(ii)*Interval),ymax);
% end
% 
%  xlabel('t(s)');
%  ylabel('Cell NO');
%  xmin = 0;
%  xmax = nFrames/Frequency;
%  ymin=0;
%  axis([xmin xmax ymin ymax+Interval]);
% 
%  set(gca, 'YTick', [Interval Interval*5:Interval*5:Cellnum*Interval]);
%  set(gca,'YTickLabel',y_labe);
%  set(gca,'FontName','Times New Roman','FontSize',14);

 
 
 
