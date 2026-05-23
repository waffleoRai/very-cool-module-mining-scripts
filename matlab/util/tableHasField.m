%
%%
function boolRes = tableHasField(inputTable, fieldName)
    vnames = inputTable.Properties.VariableNames;
    boolRes = ismember(fieldName, vnames);
end