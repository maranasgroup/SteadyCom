function indCom = infoCom2indCom(modelCom,infoCom,revFlag,spAbbr,spName)
%Transform between community reaction IDs and reaction names
%
%indCom = infoCom2indCom(modelCom)
%   Get the IDs from names. modelCom must have .infoCom
%indCom = infoCom2indCom(modelCom, infoCom)
%   Supply infoCom if it is not a field of modelCom
%infoCom = infoCom2indCom(modelCom,indCom,true,spAbbr,spName)
%   Get the name structure infoCom from ID structure indCom
%   spAbbr (must be supplied if revFlag = true) : the abbreviation of each organism
%   spName (optional): the name of each organism
if nargin < 2
    if ~isfield(modelCom,'infoCom')
        error('infoCom must be provided.\n');
    end
    infoCom = modelCom.infoCom;
end
if nargin < 3
    revFlag = false;
end
indCom = struct();
if ~revFlag
    %from infoCom to indCom
    indCom.spBm = findRxnIDs(modelCom,infoCom.spBm);
    if isfield(infoCom,'spATPM')
        indCom.spATPM = findRxnIDs(modelCom,infoCom.spATPM);
    end
    if isfield(infoCom,'rxnSD')
        indCom.rxnSD = findRxnIDs(modelCom,infoCom.rxnSD);
    end
    indCom.EXcom = findRxnIDs(modelCom,infoCom.EXcom);
    indCom.EXsp = zeros(size(infoCom.EXsp));
    SpCom = ~cellfun(@isempty,infoCom.EXsp);
    indCom.EXsp(SpCom) = findRxnIDs(modelCom,infoCom.EXsp(SpCom));
    indCom.Mcom = findMetIDs(modelCom,infoCom.Mcom);
    indCom.Msp = zeros(size(infoCom.Msp));
    SpCom = ~cellfun(@isempty,infoCom.Msp);
    indCom.Msp(SpCom) = findMetIDs(modelCom,infoCom.Msp(SpCom));
    [~,indCom.rxnSps] = ismember(infoCom.rxnSps,infoCom.spAbbr);
    [~,indCom.metSps] = ismember(infoCom.metSps,infoCom.spAbbr);
else
    %from indCom to infoCom
    if nargin < 4
        error('spAbbr must be provided to get the organisms'' abbreviations');
    end
    if nargin < 5
        spName = spAbbr;
    end
    indCom.spBm = modelCom.rxns(infoCom.spBm);
    if isfield(infoCom,'spATPM')
        indCom.spATPM = modelCom.rxns(infoCom.spATPM);
    end
    if isfield(infoCom,'rxnSD')
        indCom.rxnSD = modelCom.rxns(infoCom.rxnSD);
    end
    indCom.EXcom = repmat({''},size(infoCom.EXcom,1),2);
    indCom.EXcom(infoCom.EXcom~=0) = modelCom.rxns(infoCom.EXcom(infoCom.EXcom~=0));
    indCom.EXsp = repmat({''},size(infoCom.EXsp,1),size(infoCom.EXsp,2));
    SpCom = infoCom.EXsp ~= 0;
    indCom.EXsp(SpCom) = modelCom.rxns(infoCom.EXsp(SpCom));
    indCom.Mcom = modelCom.mets(infoCom.Mcom);
    indCom.Msp = repmat({''},size(infoCom.Msp,1),size(infoCom.Msp,2));
    SpCom = infoCom.Msp ~= 0;
    indCom.Msp(SpCom) = modelCom.mets(infoCom.Msp(SpCom));
    indCom.spAbbr = spAbbr;
    indCom.spName = spName;
    indCom.rxnSps = repmat({'com'},numel(modelCom.rxns),1);
    indCom.rxnSps(infoCom.rxnSps > 0) = spAbbr(infoCom.rxnSps(infoCom.rxnSps > 0));
    indCom.metSps = repmat({'com'},numel(modelCom.mets),1);
    indCom.metSps(infoCom.metSps > 0) = spAbbr(infoCom.metSps(infoCom.metSps > 0));
end
end