function selectAUCdata = calculateAUC(select_TIC,retention_range) % calculate area between left and right index under TIC curve
if length(select_TIC(:,1)) > 1
    leftind = find(select_TIC(:,1) == retention_range(1,1));
    rightind = find(select_TIC(:,1) == retention_range(1,2));
    selectAUCdata = trapz(select_TIC(leftind:rightind,1),select_TIC(leftind:rightind,2)); % calculate area using trapezoidal method
else
    selectAUCdata = -1;
end
end