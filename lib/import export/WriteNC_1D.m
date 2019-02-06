function [] = WriteNC_1D(namefile, time, data, varname, unit, long_varname)

if time(1)<3000
    % then it is decimal year
    time = datenum(time,1,1) - datenum(1900,1,1,0,0,0);
elseif time(1)>401767
    % then it is in Matlab format
    time = time - datenum(1900,1,1,0,0,0);
else
    % then it is already in the good format
    warning('Time vector already in days since 1900-1-1 0:0:0')
end

     WriteNetCDF2(namefile,time, ...
         'VariableName','time', ...
         'DimensionNames',{'time'},...
         'Unit', 'days since 1900-1-1 0:0:0', ...
         'LongVariableName','Time (UTC)');


         for i=1:length(data)
                 WriteNetCDF2(namefile,data{i}, ...
                     'VariableName',varname{i}, ...
                     'DimensionNames',{'time'},...
                     'Unit', unit{i}, ...
                     'LongVariableName',long_varname{i});
         end
end