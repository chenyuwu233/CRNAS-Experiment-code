function ret = get_indi(theta,d_max)
    if length(theta) ==12 
        GR50_1 = GR_50(theta(2) - theta(3), theta(4), theta(5), theta(6),d_max);
        GR50_2 = GR_50(theta(7) - theta(8), theta(9), theta(10), theta(11),d_max);
        ret = [theta(1),theta(2) - theta(3), theta(7)- theta(8), GR50_1, GR50_2];
    elseif length(theta) == 13
        GR50_1 = GR_50(theta(2) - theta(3), theta(4), theta(5), theta(6),d_max);
        GR50_2 = GR_50(theta(8) - theta(9), theta(10), theta(11), theta(12),d_max);
        ret = [theta(1),theta(2) - theta(3), theta(8)- theta(9), GR50_1, GR50_2];
    else
        GR50_1 = GR_50(theta(1) , theta(2), theta(3), theta(4),d_max);
        GR50_2 = GR_50(theta(6) , theta(7), theta(8), theta(9),d_max);
        ret = [theta(5),theta(1), theta(6), GR50_1, GR50_2];
    end
end