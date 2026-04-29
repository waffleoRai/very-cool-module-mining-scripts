%
%%
function val = alignDown(val, snap)
    a = floor(val / snap);
    val = a * snap;    
end