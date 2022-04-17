function [t,y_interp] = util_interp(y)
    y_interp = y;
    nandata = isnan(y_interp);
    t = 1:numel(y_interp);
    if ~isempty(y_interp(~nandata))
        y_interp(nandata) = interp1(t(~nandata), y_interp(~nandata), t(nandata));
    end
    