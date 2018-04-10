function afterzone = mayaafterzones(timearray,window)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

maxtime = size(timearray,1);
z = [];

for i = 1:maxtime
            if fix(timearray(i)) + window < 1500
                z(i,1) = fix(timearray(i)) + 1;
                z(i,2) = fix(timearray(i)) + window;
            else
    end
afterzone = z;

end

