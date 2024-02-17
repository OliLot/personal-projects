% flooring number to sf significant digits
function result = floor_sf(num, sf)
if (num < 0)
    num = -num;
    correction = floor(log10(num));
    scale_factor = 10^(sf - correction);
    fullstring = string(- num * scale_factor);
    charArray = char(fullstring);
    result = string(charArray(1:5));
else
    correction = floor(log10(num));
    scale_factor = 10^(sf - correction);
    fullstring = string(num * scale_factor);
    charArray = char(fullstring);
    result = string(charArray(1:4));
end

