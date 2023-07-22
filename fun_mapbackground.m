function [lonconv, latconv] = fun_mapbackground(lolim, lalim, lon, lat); 

m_proj('miller', 'long',lolim + [-2 2],'lat',lalim + [-2 2]);
m_coast('patch',[1 1 1]); % m_coast('patch',[1 .85 .7]);

[latbord, lonbord] = borders('states'); % add states map
for iplace = 1:length(lonbord); 
    m_line(lonbord{iplace}, latbord{iplace}, 'LineWidth',1,'color',0*[1 1 1])
end

m_grid('box','fancy','linestyle','-','gridcolor','w','backcolor',.8.*[1, 1, 1]);

[lonconv, latconv] = m_ll2xy(lon, lat);  

end