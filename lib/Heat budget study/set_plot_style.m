function [] = set_plot_style(ii, station, time_mod, ylab)

    axis tight
    box on 
    time_ext = datenum(1995,1,1):1/24:datenum(2016,12,1);
    set_monthly_tick(time_ext);
    temp = get(gca,'XTickLabel');
    for k = 1:length(temp)
        if k/2==floor(k/2)
            temp(k,:)=' ';
        end
    end
    h_text = text(time_ext(24*60), ...
       -18,...
        sprintf('%s) %s',char(ii+96),station{ii}));
    h_text.FontSize = 15;
    h_text.FontWeight = 'bold';
    h_text.Units = 'Normalized';
    h_text.Position(1:2) = [0.05 0.9];
    xlim([datenum(1996,1,1) datenum(2015,6,1)])

    if ismember(ii,1:6)
       set(gca,'XTickLabel','')
    else
       xlabel('Year')
       set_monthly_tick(time_mod{ii})
    end
   xlim([datenum(1996,1,1) datenum(2015,1,1)])
   if ismember(ii,[2 3 5 6 8 9])
%        set(gca,'YTickLabel','')
   end
   if ii==4
        ylabel(ylab,'Interpreter','tex')
   end
end