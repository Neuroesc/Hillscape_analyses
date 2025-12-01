%% >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> Sub function for map masks
function [msk] = get_map_mask(m1,type)
    switch type
        case {'overlapping','overlap'}
            % correlations limited to overlapping bins
            ps = linspace(0,size(m1,2),7);
            i = knnsearch(ps(:),(1:size(m1,2))');
            t_idx = 1:2:length(ps);
            msk = ismember(repmat(i',size(m1,1),1),t_idx); % matrix same size as m1, 1s for points closes to troughs, 0s otherwise

        case {'edges','edge'}
            % correlations limited to edges
            e_idx = zeros(size(m1));
            e_idx(:,1) = 1;
            e_idx(:,end) = 1;
            e_idx(1,:) = 1;
            e_idx(end,:) = 1;
            e_dist = bwdist(e_idx);
            msk = e_dist<5;

        case {'center','centre'}
            % correlations limited to center
            e_idx = zeros(size(m1));
            e_idx(:,1) = 1;
            e_idx(:,end) = 1;
            e_idx(1,:) = 1;
            e_idx(end,:) = 1;
            e_dist = bwdist(e_idx);
            msk = e_dist>10;

        case {'corners','corner'}
            % correlations limited to corners
            c_idx = zeros(size(m1));
            c_idx(1,1) = 1;
            c_idx(1,end) = 1;
            c_idx(end,1) = 1;
            c_idx(end,end) = 1;
            c_dist = bwdist(c_idx);
            msk = c_dist<15;

        otherwise
            keyboard
    end

end