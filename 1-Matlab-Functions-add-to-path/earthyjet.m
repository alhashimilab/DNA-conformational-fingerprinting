function cmap_OS = earthyjet(m)
    if nargin < 1
        m = size(get(gcf,'colormap'),1);
    end

    % Define key earthy colors (dark blue → teal → olive → ochre → brick red)
   % base_colors = [ ...
   %      45 85 170;     % vivid blue
   %      80 160 170;    % teal
   %      110 150 120;    % olive green
   %      240 240 100;   % soft warm yellow
   %      255 170 60;    % golden orange
   %      235 135 90;    % coral / salmon
   %      204 102 119];  % earthy rose

   base_colors = [ ...
    45 85 170;      % vivid blue
    80 160 170;     % teal
    110 150 120;    % soft olive green
    240 255 50;    % soft yellow
    235 170 60;     % golden orange
    215 110 95;     % coral (new bridge)
    204 102 119];   % earthy rose

    % Normalize to [0,1]
    base_colors = base_colors ./ 255;

    % Interpolate to m colors
    cmap_OS = interp1(linspace(0,1,size(base_colors,1)), base_colors, linspace(0,1,m), 'linear');
end
