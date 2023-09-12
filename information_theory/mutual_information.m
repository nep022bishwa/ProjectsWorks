function mi =mutual_information(x,y,L)
%  returns CMI(x,y|z) and MI(mi_option)
%  Inputs are three time series or images
%  mi_option = 1 = (x,y)
%  mi_option = 2 = (x,z)
%  mi_option = 3 = (y,z)
%  mi_option = 4 = H(y)+H(x,z)-H(x,y,z)

counts = length(x(:));
eps = 1e-8;
xn = normalize(x(:),'range')*(1-eps);
yn = normalize(y(:),'range')*(1-eps);
%zn = normalize(z(:),'range')*(1-eps);

%quantize the data

xL = floor(xn*L);
yL = floor(yn*L);
%zL = floor(zn*L);

%convert to a 1d label

xy_bin = xL + L*yL;
%yz_bin = yL + L*zL;
%xz_bin = xL + L*zL;
%xyz_bin = xL + L*yL + L^2*zL;

% now compute probabilities

[s c_x] = mmrepeat(sort(xL));
[s c_y] = mmrepeat(sort(yL));
%[s c_z] = mmrepeat(sort(zL));
[s c_xy] = mmrepeat(sort(xy_bin));
%[s c_yz] = mmrepeat(sort(yz_bin));
%[s c_xz] = mmrepeat(sort(xz_bin));
%[s c_xyz] = mmrepeat(sort(xyz_bin));

p_x = c_x/counts;
p_y = c_y/counts;
%p_z = c_z/counts;
p_xy = c_xy/counts;
%p_yz = c_yz/counts;
%p_xz = c_xz/counts;
%p_xyz = c_xyz/counts;

entropy = @(p) - sum(p.*log2(p));

%if mi_option == 1
    mi = entropy(p_x) + entropy(p_y) - entropy(p_xy);
    %{
elseif mi_option == 2
    mi = entropy(p_x) + entropy(p_z) - entropy(p_xz);
elseif mi_option == 3
    mi = entropy(p_y) + entropy(p_z) - entropy(p_yz);
elseif mi_option == 4
    mi = entropy(p_y) + entropy(p_xz) - entropy(p_xyz);
end

cmi = entropy(p_xz) + entropy(p_yz) - entropy(p_xyz) - entropy(p_z);
    %}
end

