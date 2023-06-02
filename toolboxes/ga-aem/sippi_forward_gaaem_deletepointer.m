% sippi_forward_gaaem_deletepointer: Delete pointer in forward
function forward = sippi_forward_gaaem_deletepointer(forward)

if isfield(forward.LM,'hS')
    sippi_verbose(sprintf('%s: deleting pointer for LM',mfilename),1)
    forward.LM=rmfield(forward.LM,'hS');
end
if isfield(forward.HM,'hS')
    sippi_verbose(sprintf('%s: deleting pointer for HM',mfilename),1)
    forward.HM=rmfield(forward.HM,'hS');
end


