function T_out = loadTable(outputMatFile)

if exist(outputMatFile, 'file') == 2
    S = load(outputMatFile);
else
    error('Give correct output file path');    
end
%rename what we loaded
fn = fieldnames(S);
T_out = S.(fn{1});



end