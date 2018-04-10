function beforezone = mayabeforezones(timearray, window)
%March 14, 2018 Maya Erler
%Takes as an input an array with times of significant events and a window
%which determines the size of the zone
%As 

maxtime = size(timearray,1);
z = [];
counter = 1;
for i = 1:maxtime
            if fix(timearray(i))-window > 0
                z(counter,2) = fix(timearray(i)) - 1;
                z(counter,1) = fix(timearray(i)) - window;
                counter = counter + 1;
            else
                
    end
beforezone = z;
end

