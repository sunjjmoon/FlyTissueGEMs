%% KEGG API Query String Generator
% --------------------------------------------------------
% Purpose:
%   Generate query strings from a list of KEGG reaction IDs (e.g., RXXXX)
%   to query the KEGG API and retrieve associated EC numbers.
%
% Usage:
%   1. Divide `uniqKEGGid` into manageable chunks (≤100 IDs).
%   2. Call the function with a specific index range (e.g., 1:80).
%   3. Paste the output query string to the KEGG API URL:
%      → https://rest.kegg.jp/link/ec/[query_str]
%   4. Save the response as a `.txt` file.
%   5. Use a script (e.g., `C_script_run_linkEC2rxn.m`) to map ECs back to reactions.
%
% Example:
%   param = 1:80;
%   query_str = parseString4Keggapi(uniqKEGGid, param);
%   → paste into URL: https://rest.kegg.jp/link/ec/[query_str]
%
%   Repeat for: param = 81:length(uniqKEGGid)
% --------------------------------------------------------

%% Define Parameters
uniqKEGGid = readtable(strcat(pwd,'\1_rxn_woEC\uniqKEGGid_list.csv'));
uniqKEGGid = table2cell(uniqKEGGid);
param = 1:80;
query_str_1 = parseString4Keggapi(uniqKEGGid, param);

param = 81:length(uniqKEGGid);
query_str_2 = parseString4Keggapi(uniqKEGGid, param);

%% Function Definition
function query_str = parseString4Keggapi(uniqKEGGid, param)
    % Extract a subset of KEGG IDs and join with '+' for KEGG API
    subset_ids = uniqKEGGid(param);
    query_str = strjoin(subset_ids, '+');
end
