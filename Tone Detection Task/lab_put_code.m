%function lab_put_code(L,ecode)
%Currently only works for ecodes 1-8 
%Requires lab_init first to create the L structure
%
 function [e] = lab_put_code(L,ecode)

try

L.ljudObj.ePut(L.ljhandle, LabJack.LabJackUD.IO.PUT_DIGITAL_BIT, L.chan{ecode}, 1, 0);
WaitSecs(.002);
L.ljudObj.ePut(L.ljhandle, LabJack.LabJackUD.IO.PUT_DIGITAL_BIT, L.chan{ecode}, 0, 0);

e = 0;
catch e
end