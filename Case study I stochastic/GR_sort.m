% This is a function for sorting the parameter according to the GR-50 from
% small to large (This function only works for 2 sub-group hl-dyn-sto model)


function ret = GR_sort(par,d_max)
    ret = par;
    if length(par) == 12
        indi = get_indi(par,d_max);
        GR = indi(4:5);
        GR_sort = sort(GR);
        if GR ~= GR_sort
            ret(1) = 1 - ret(1);
            temp = ret(2:6);
            ret(2:6) = ret(7:11);
            ret(7:11) = temp;
        end
    elseif length(par) == 13
        indi = get_indi(par,d_max);
        GR = indi(4:5);
        GR_sort = sort(GR);
        if GR ~= GR_sort
            ret(1) = ret(7);
            ret(7) = 1 - ret(1);
            temp = ret(2:6);
            ret(2:6) = ret(8:12);
            ret(8:12) = temp;
        end
    elseif length(par) == 11
        indi = get_indi(par,d_max);
        GR = indi(4:5);
        GR_sort = sort(GR);
        if GR~= GR_sort
            ret(5) = 1 - ret(5);
            temp = ret(1:4);
            ret(1:4) = ret(6:9);
            ret(6:9) = temp;
        end
    end
end