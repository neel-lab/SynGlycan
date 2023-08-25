function output = base64decode(input)
% BASE64DECODE: Decode Base64 string to a byte array.
% 
% Syntax:
% output = base64decode(input)
% 
% Input:
% input: string. Data encoded in Base64
% 
% Output:
% output: unit8 array. Decoded data.
% 
% Note:
% JAVA must be available for function to work.
% 
% Example:
% N/A. Set breakpoints while running PREPROCESSGUI
% 
% Children function: 
% n/a
% 
% See Also:
% BASE64ENCODE  GETMSDATA  PREPROCESSGUI

% GlycoPAT 2 authors: Kai Cheng, Gang Liu, Gabrielle Pawlowski, Yusen Zhou, Sriram Neelamegham
% (c) 2020, Research Foundation for State University of New York. All rights reserved
% Date Lastly Updated: 11/10/2020

error(nargchk(1, 1, nargin));
error(javachk('jvm'));
if ischar(input)
    input = uint8(input);
end
output = typecast(org.apache.commons.codec.binary.Base64.decodeBase64(input), 'uint8')';

end

