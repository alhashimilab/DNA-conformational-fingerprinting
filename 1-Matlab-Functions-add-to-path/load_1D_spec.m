function [axis_nmr,spec_nmr] = load_1D_spec(parent_path,filename,F2_left,F2_right)

filepath = sprintf('%s/%s',parent_path,filename);

if isfile(filepath)

    specdata = importdata(filepath);

    spec_nmr     = specdata(:,2);

    F2 = length(specdata);

    % finding axis tick values from edges of spectral window (defined by EXT in
    % nmrPipe, or by the actual SW limits in Topspin)

    sw2 = F2_left - F2_right;

    axis_nmr = (F2_left:-sw2/(F2-1):F2_right)';

else

    axis_nmr = [];
    spec_nmr = [];

end

end