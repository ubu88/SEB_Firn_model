function [data_site] = CompactionAnalysis(DataCompaction, MetadataCompaction, ID)
   if ischar(ID)
       site = ID;
        site_num = find(strcmp(MetadataCompaction.sitename,site));

        disp('Instruments ID at that site:')
        instrument_list = MetadataCompaction.FC_instrument_IDs(site_num,MetadataCompaction.FC_instrument_IDs(site_num,:)~=0);
        disp(instrument_list)
        ID = instrument_list;
        data_site = CompactionAnalysis(DataCompaction,MetadataCompaction,ID);
   else
       for i = 1:length(ID)
           data_site{i} = DataCompaction{ID(i)};
            i_site=[];
            for ii =1:length(MetadataCompaction.sitename)
                temp = strfind(MetadataCompaction.sitename(ii),unique(data_site{i}.sitename));
                if ~isempty(temp{1})
                    i_site = ii;
                    break
                end
            end
            time_obs = datetime(datestr(data_site{i}.time));
            save{i}=time_obs;
            [ii, jj] = find(MetadataCompaction.FC_instrument_IDs==ID(i));
            
            ind_nan = isnan(data_site{i}.Compaction_Distance_m);

            % The "Compaction_Distance_m" given in Mike's files is the
            % length of a 2 m cable linked to the rewinding system. This
            % cable i further extended to a secondary cable down to the
            % bottom of the borehole. Therefor the borehole length reads as:            
            
            % Borehole Distance = Wire Measurement + (*Original* Borehole
            % Distance [metadata] - Original Wire Measurement at Day-0)
        
            data_site{i}.Compaction_Distance_m(~ind_nan) = ...
                data_site{i}.Compaction_Distance_m(~ind_nan) ...
                + (abs(MetadataCompaction.FC_borehole_initial_length_m(ii,jj)) ...
                - data_site{i}.Compaction_Distance_m(find(~ind_nan,1,'first')));
            data_site{i}.Compaction_Distance_m_1 = smooth(data_site{i}.Compaction_Distance_m,2*7,'lowess');
            data_site{i}.Compaction_Distance_m_2 = smooth(data_site{i}.Compaction_Distance_m,4*7,'lowess');

            data_site{i}.Compaction_Distance_m_1(ind_nan) = NaN;
            data_site{i}.Compaction_Distance_m_2(ind_nan) = NaN;
            a = data_site{i}.time(1:end-1)-data_site{i}.time(2:end);
            if sum(a>=0)>0
                fprintf('Duplicate time detected for instrument %i',ID(i))
                ind_err = find(a==0);
                format longg
                tab_err = table(data_site{i}.daynumber_YYYYMMDD(ind_err), ...
                    data_site{i}.daynumber_YYYYMMDD(ind_err +1), ...
                    data_site{i}.daynumber_YYYYMMDD(ind_err -1), ...
                    data_site{i}.Compaction_Distance_m(ind_err), ...
                    data_site{i}.Compaction_Distance_m(ind_err+1),...
                    'VariableName',{'Time_1','Time_2','Time_prev','Val_1','Val_2'})
            end
            
            data_site{i}.Compaction_Rate_md = ...
                [data_site{i}.Compaction_Distance_m(1:end-1)-data_site{i}.Compaction_Distance_m(2:end); 0];
            data_site{i}.Compaction_Rate_md_1 = ...
                [data_site{i}.Compaction_Distance_m_1(1:end-1)-data_site{i}.Compaction_Distance_m_1(2:end); 0];
            data_site{i}.Compaction_Rate_md_2 = ...
                [data_site{i}.Compaction_Distance_m_2(1:end-1)-data_site{i}.Compaction_Distance_m_2(2:end); 0];
            ind_err = ([data_site{i}.time(1:end-1)-data_site{i}.time(2:end); 1] ~= -1);
            data_site{i}.Compaction_Rate_md_2(ind_err) = NaN;
            data_site{i}.Compaction_Rate_md_1(ind_err) = NaN;
            data_site{i}.Compaction_Rate_md_0(ind_err) = NaN;
            
       end
       
    f= figure;
    num_subplot = 4;
    ha = tight_subplot(num_subplot,1,0.035, [.07 .05], 0.06);
    count = 0;
    for i = 1:length(ID)
        count = count+1;
        if count >num_subplot
            axes(ax(1))
            xlabel(gca,'Date')
            datetick('x','mm-yyyy','keepticks')
            xlim(gca,[min(data_site{1}.time),max(data_site{1}.time)])

            axes(ax(2))
            xlabel(gca,'Date')
            datetick('x','mm-yyyy','keepticks')
            xlim(gca,[min(data_site{1}.time),max(data_site{1}.time)])

            i_file = 1;
            NameFile = sprintf('./Output/%s_%i',...
                MetadataCompaction.sitename{i_site},i_file);
            while exist(strcat(NameFile,'.tif'), 'file') == 2
                i_file = i_file + 1;
                NameFile = sprintf('%s_%i',NameFile(1:end-2),i_file);
            end
            print(f,NameFile,'-dtiff'); 

            f= figure;
            ha = tight_subplot(num_subplot,1,0.03, [.07 .05], 0.06);
            count = 1;
        end
        
        time_obs = datetime(datestr(data_site{i}.time));
        [ii, jj] = find(MetadataCompaction.FC_instrument_IDs==ID(i));

        text_title = sprintf('%s - Instr. #%i - Length: %0.2f m - Depth of top: %0.2f m',...
            MetadataCompaction.sitename{ii} , ID(i),...
            -MetadataCompaction.FC_borehole_initial_length_m(ii,jj),...
            -MetadataCompaction.FC_borehole_top_from_surface_m(ii,jj));
            
        axes(ha(count));
        if sum(~isnan(data_site{i}.Compaction_Distance_m_2))<10
            text_title = sprintf('%s - Instr. #%i - No data available',...
                MetadataCompaction.sitename{ii} , ID(i));  
            title(text_title)
        else
            hold on
            [ax,h1, h2] = plotyy(data_site{i}.time,data_site{i}.Compaction_Distance_m_2,...
                data_site{i}.time,data_site{i}.Compaction_Rate_md_2.*1000);
            ylim(ax(2),[-0.5, 2])
            h1.LineWidth = 2;
            h2.LineWidth = 2;

            ind_nan = isnan(data_site{i}.Compaction_Distance_m);
            ii=find(~ind_nan,1,'first');

            axes(ax(1))
            hold on
            plot(data_site{i}.time(ii:(ii+60)), ...
                data_site{i}.Compaction_Distance_m_2(ii:(ii+60)), ...
                'Color',RGB('light light blue'),'LineWidth',2);
            set(gca,'XTickLabel',[],...
                'YLim',[min(data_site{i}.Compaction_Distance_m_2), max(data_site{i}.Compaction_Distance_m_2)],...
                'XLim',[min(data_site{1}.time),max(data_site{1}.time)])
            set_monthly_tick(time_obs,gca);

            axes(ax(2))
            hold on
            plot(data_site{i}.time, ...
                data_site{i}.Compaction_Rate_md.*1000, ...
                'Color',[0.8 0.8 0.8]);
    %         plot(data_site{i}.time, ...
    %             data_site{i}.Compaction_Rate_md_1.*1000, ...
    %             'Color',[0.5 0.5 0.5]);
            plot(data_site{i}.time, ...
                data_site{i}.Compaction_Rate_md_2.*1000, ...
                'r','LineWidth',2);

            plot(data_site{i}.time(ii:(ii+60)), ...
                data_site{i}.Compaction_Rate_md_2(ii:(ii+60)).*1000, ...
                'Color',RGB('light light red'),'LineWidth',2);

            set(gca,'XTickLabel',[],...
                'YLim',[-0.5 2],...
                'YTick',-0.5:0.5:2,...
                'YTickLabel',-0.5:0.5:2,...
                'XLim',[min(data_site{1}.time),max(data_site{1}.time)])
            set_monthly_tick(time_obs,gca);

            title(text_title)
            ylabel(ax(1),'Comp. distance (m)')
            ylabel(ax(2),'Comp. rate (mm/day)')
        end
    end
    
    if count >num_subplot
        axes(ax(1))
        xlabel(gca,'Date')
        datetick('x','mm-yyyy','keepticks')

        axes(ax(2))
        xlabel(gca,'Date')
        datetick('x','mm-yyyy','keepticks')

        i_file = 1;
        NameFile = sprintf('./Output/%s_%i',...
            MetadataCompaction.sitename{i_site},i_file);
        while exist(strcat(NameFile,'.tif'), 'file') == 2
            i_file = i_file + 1;
            NameFile = sprintf('%s_%i',NameFile(1:end-2),i_file);
        end
        print(f,NameFile,'-dtiff'); 

        f= figure;
        ha = tight_subplot(num_subplot,1,0.03, [.07 .05], 0.06);
        count = 1;
    end
        
    axes(ha(count+1));
    plot(data_site{1}.time,data_site{1}.AirTemp_median_C,'k','LineWidth',2)
    axis tight
    ylabel('Air temperature (^oC)','Interpreter','tex')
    set(gca,'XTickLabel',[],'XLim',[min(data_site{1}.time),max(data_site{1}.time)])
    set_monthly_tick(time_obs,gca);

    axes(ha(count+2));
    ind_nonan = ~isnan(data_site{1}.SonicRangeDist_Corrected_m);
    ii=find(ind_nonan,1,'first');
    plot(data_site{1}.time,...
        hampel(data_site{1}.SonicRangeDist_Corrected_m(ii)-data_site{1}.SonicRangeDist_Corrected_m,14),...
        'k','LineWidth',2)
    axis tight
    ylabel('Surface Height (m)','Interpreter','tex')
    set(gca,'YAxisLocation','right','XLim',[min(data_site{1}.time),max(data_site{1}.time)])
    set_monthly_tick(time_obs,gca);
    datetick('x','mm-yyyy','keeplimits','keepticks')            
    xlabel('Date')
    
    if count+2<num_subplot
        for ii = (count+3) :num_subplot
            set(ha(ii),'Visible','off')
        end
    end
    
    i_file = 1;
    NameFile = sprintf('./Output/%s_%i',...
        MetadataCompaction.sitename{i_site},i_file);
    while exist(strcat(NameFile,'.tif'), 'file') == 2
        i_file = i_file + 1;
        NameFile = sprintf('%s_%i',NameFile(1:end-2),i_file);
    end
    print(f,NameFile,'-dtiff'); 
   end
end