function idx = strIdx(CellWithStrings,string)
% use this function to find a specific string in a cell containing only strings
% function returns error if the cell does not contain the searched string or if the cell does contain the string several times.
% Example: CellWithString={'apples','peaches','strawberries','apples'}; string='strawberries'; idx = strIdx(CellWithStrings,string) returns 3.
% Example 2: CellWithString={'apples','peaches','strawberries','apples'}; string='apples'; idx = strIdx(CellWithStrings,string) returns
% error.
% Example 3: CellWithString={'apples','peaches','strawberries','apples'}; string='raspberries'; idx = strIdx(CellWithStrings,string) returns
% error.

if sum(ismember(CellWithStrings,string))==0
    error('string is not in CellWithStings');
elseif sum(ismember(CellWithStrings,string))==1
    idx = find(ismember(CellWithStrings,string));
else
    error('string is more than once in CellWithStrings');
end