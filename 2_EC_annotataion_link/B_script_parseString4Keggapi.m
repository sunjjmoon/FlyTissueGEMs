%  Goal:
%   - To parse the KEGGid to run the KEGG API.
% 
%  After obatining the query_str, writhe the following in the URL:
%       rest.kegg.jp/link/ec/ ~~~
%           ~~~ is the query_str (copy and paste)
%       This will generate the kegg-id linked with EC. Non-linked EC will
%       not be included.
%  Save it as txt file, and  use 'C_script_run_linkEC2rxn.m' to associate
%  linked EC to reactions.
%   expression: 
%       -> https://www.kegg.jp/entry/~~
%           ~~ is the KeggID. 
%   Input:
%    param: range of the uniqKeggID that will be input in the KEGG API
%    e.g.) param = 1:80;
%    e.g.) param = 81:end;
%       since the KEGG API cannot handle ~ more than 100 queries, splite
%       the list and iterate.

%Procedure
%1.
param = 1:80;
query_str = parseString4Keggapi(uniqKEGGid,param);
%2. 
% Go to "rest.kegg.jp/link/ec/(copy and paste of the query_str)
% 3. Save it as a text file and in the folder B
% Repeat with the rest of the param
param = 81:length(uniqKEGGid);
query_str = parseString4Keggapi(uniqKEGGid,param);

function query_str = parseString4Keggapi(uniqKEGGid,param)
str = uniqKEGGid(param);
query_str= strjoin(str,'+');
end