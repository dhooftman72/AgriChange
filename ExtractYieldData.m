function [LTETesterSelect,AverageYield,MaxYieldTreat,MinimumYieldTreat] = ExtractYieldData(File)
Filename = [File,'.mat'];
load(Filename)
Array = eval(File);
FileShort = File(1:3);
ListExclude = find(strcmp(Array.Exclude,'Y')==1);
Array(ListExclude,:) = [];
isFieldResult = myIsField (Array, 'Season');
if isFieldResult ~=1
    Array.Season(:,1) = {'LR'};
end

LTETester = dataset({FileShort},'Varnames', 'LTE');
if length(File)>3
    SubLTE = File(5:end);
else
    SubLTE = 'One Crop';
end
LTETester.SubLTE = {SubLTE};
LTETester.Treatment = unique(Array.TreatType);
for i = 1: length(LTETester.Treatment)
    Tester = find(strcmp(Array.TreatType,LTETester.Treatment(i))==1);
    LTETester.SubLTE(i,1) = {SubLTE};
    LTETester.FullLTE(i,1) = {File};
    LTETester.MeanYields(i,1) = nanmean(Array.Yield(Tester)); % for checking purposes
    LTETester.Years(i,1) = length(unique(Array.Year(Tester)));
    LTETester.StartYear (i,1) = min(unique(Array.Year(Tester)));
    LTETester.EndYear (i,1) = max(unique(Array.Year(Tester)));
    LTETester.Nutrients(i,1) = (Array.Nutrients(Tester(1,1)));
    LTETester.Residue(i,1) = (Array.Residue(Tester(1,1)));
    LTETester.CType(i,1) = (Array.CType(Tester(1,1)));
    LTETester.LTE(i,1) = {FileShort};
end

MaximumYears = max(LTETester.Years);
LTETester.PercYears =  LTETester.Years./MaximumYears;
NinetyList = find(LTETester.PercYears >= 0.9);
LTETesterSelect =  LTETester (NinetyList,:);

clear LTETester NinetyList MaximumYears i Tester ListExclude
SelectedTreatments = LTETesterSelect.Treatment;
for i = 1:length(SelectedTreatments)
    List = find(strcmp(Array.TreatType,SelectedTreatments(i))==1);
    if i == 1
        SelectedTreatmentsperLTE = Array(List,:);
    else
        SelectedTreatmentsperLTE = [SelectedTreatmentsperLTE;Array(List,:)];  %#ok<*AGROW>
    end
end
clear i List


% Stats
Yield = SelectedTreatmentsperLTE.Yield;
Year = SelectedTreatmentsperLTE.Year;
TreatType = SelectedTreatmentsperLTE.TreatType;
Season = SelectedTreatmentsperLTE.Season;
Replica1 = SelectedTreatmentsperLTE.Replica1;
Replica2 = SelectedTreatmentsperLTE.Replica2;

% Selection for coefficients per Treatment type
[~,StatsOutput,stats]= anovan((Yield),{Year, TreatType, Season, Replica1, Replica2},'sstype',1,...
    'model',[1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0 ; 0 0 0 1 0; 0 0 0 0 1],'nested',[0 0 0 0 0; 0 0 0 0 0; 0 1 0 0 0; 0 1 0 0 0; 0 0 0 1 0 ],'display', 'off',...
    'continuous',1, 'random', [2 3 4 5],'varnames', {'Year', 'TreatType' , 'Season', 'Replica1', 'Replica2'}); 
Endvalue = (3 + size(LTETesterSelect,1)) -1;
LTETesterSelect.NameCheck = stats.coeffnames(3:Endvalue);
LTETesterSelect.Coeffs = stats.coeffs(3:Endvalue);

% Per year statistics 
for i = 1:length(SelectedTreatments)
    List = find(strcmp(SelectedTreatmentsperLTE.TreatType,SelectedTreatments(i))==1);
    PerTreat = SelectedTreatmentsperLTE(List,:);
    Yield = PerTreat.Yield;
    Year = PerTreat.Year-1970;
    Season = PerTreat.Season;
    Replica1 = PerTreat.Replica1;
    Replica2 = PerTreat.Replica2;
    [~,StatsOutput2,stats2]= anovan((Yield),{Year, Season, Replica1, Replica2},'sstype',1,...
        'model',[1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1],'nested',[0 0 0 0; 0 0 0 0; 0 1 0 0; 0 0 1 0 ],'display', 'off',...
        'continuous',1, 'random', [2 3 4],'varnames', {'Year', 'Season', 'Replica1', 'Replica2'});
    LTETesterSelect.YearConstant(i,1) =  stats2.coeffs(1);
    LTETesterSelect.YearCoef(i,1) =  stats2.coeffs(2);
    LTETesterSelect.Year_FValue(i,1) =  StatsOutput2(2,6);
    LTETesterSelect.Year_PValue(i,1) =  StatsOutput2(2,7);
end
LTETesterSelect(i+1,[1,2,4,11]) = LTETesterSelect(i,[1,2,4,11]);
LTETesterSelect.Years(i+1,1) = round(nanmedian(LTETesterSelect.Years));
LTETesterSelect.StartYear(i+1,1) = round(nanmedian(LTETesterSelect.StartYear));
LTETesterSelect.EndYear(i+1,1) = round(nanmedian(LTETesterSelect.EndYear));
LTETesterSelect.Treatment(i+1,1) = {'Combined'};
LTETesterSelect.NameCheck(i+1,1) = {'Combined'};
LTETesterSelect.YearConstant(i+1,1) =  stats.coeffs(1);
LTETesterSelect.YearCoef(i+1,1) =  stats.coeffs(2);
LTETesterSelect.Year_FValue(i+1,1) =  StatsOutput(2,6);
LTETesterSelect.Year_PValue(i+1,1) =  StatsOutput(2,7);


%% MEDIAN YIELD
% Refernce Yield Median
[AverageYield,FullYears] = CalculateYearMeans(SelectedTreatmentsperLTE,-999,File);
%% MAXIMUM YIELD TREATMENT
NoCombi = (length(LTETesterSelect.Coeffs)-1);
MaxTreat = find(LTETesterSelect.Coeffs(1:NoCombi) == max(LTETesterSelect.Coeffs(1:NoCombi)));
ListMax = find(strcmp(SelectedTreatmentsperLTE.TreatType,SelectedTreatments(MaxTreat))==1); %#ok<*FNDSB>
SelectedforMax = SelectedTreatmentsperLTE(ListMax,:);
[MaxYieldTreat,~] = CalculateYearMeans(SelectedforMax,FullYears,File);


%% NO NUTRIENT OR MINIMUM YIELD]
NoCombi = (length(LTETesterSelect.Coeffs)-1);
NoNutList = find(strcmp(LTETesterSelect.Nutrients(1:NoCombi),'Y')~=1 & strcmp(LTETesterSelect.Residue(1:NoCombi),'Y')~=1);
NoNutListMin = find(LTETesterSelect.Coeffs(1:NoCombi) == min(LTETesterSelect.Coeffs(1:NoCombi)));
if isempty(NoNutList) == 1
    NoNutList = NoNutListMin;
end
if length(NoNutList) > 1
    NoNutList = find(LTETesterSelect.Coeffs(1:NoCombi) == min(LTETesterSelect.Coeffs(NoNutList)));
end
ListMin = find(strcmp(SelectedTreatmentsperLTE.TreatType,LTETesterSelect.Treatment(NoNutList))==1);
SelectedforMin= SelectedTreatmentsperLTE(ListMin,:);
[MinimumYieldTreat,~] = CalculateYearMeans(SelectedforMin,FullYears,File);

%%
save(Filename,File,'LTETesterSelect','SelectedTreatmentsperLTE','AverageYield','MaxYieldTreat','MinimumYieldTreat','StatsOutput')
end

function [Yields,ListerYears] = CalculateYearMeans(InputTreatments,ListerYears,File)
% Regression when 2010 doesn't appear
mdl = LinearModel.fit(InputTreatments.Year,InputTreatments.Yield);
c = table2array(mdl.Coefficients(1,1));
a = table2array(mdl.Coefficients(2,1));
RegValue = c + (a*2010);
clear mdl a c

% Reference Value
ValuesPerYear = find(InputTreatments.Year ==2010);
if isempty(ValuesPerYear) ~= 1
    PerYearPerLTE = InputTreatments(ValuesPerYear,:);
    Tester = unique(PerYearPerLTE.TreatType);
    PerTreatNorm = [];
    for x = 1:length(Tester)
        ListAvg = find(strcmp(PerYearPerLTE.TreatType,Tester(x))==1);
        PerTreatNorm(x) = nanmedian(PerYearPerLTE.Yield(ListAvg));
        % Normalised
    end
    ReferenceYield = nanmedian(PerTreatNorm);
else
    disp('Calcuated Reference Yield');
    ReferenceYield = RegValue;
end
clear ValuesPerYear PerYearPerLTE ListAvg PerTreatNorm Tester

if ListerYears == -999
    ListerYearMin = min(InputTreatments.Year);
    if ListerYearMin < 1970
        ListerYearMin = 1970;
    end
    ListerYearMax = max(InputTreatments.Year);
    if ListerYearMax > 2020
        ListerYearMax = 2020;
    end
    ListerYears = ListerYearMin:1:ListerYearMax;
end

TreatTyper = unique(InputTreatments.TreatType);
if length(TreatTyper)>1
    TreatTyper = {'Combined'};
end
FileShort = File(1:3);
Yields = dataset({FileShort},'Varnames', 'LTE'); %#ok<*DTSET>
if length(File)>3
    SubLTE = File(5:end);
else
    SubLTE = 'One Crop';
end
Yields.SubLTE = {SubLTE};
warning off 
for i = 1: length(ListerYears)
    ValuesPerYear = find(InputTreatments.Year ==ListerYears(i));
    Yields.LTE(i,1) = {FileShort};
    Yields.SubLTE(i,1) = {SubLTE};
    Yields.FullLTE(i,1) = {File};
    if isempty(ValuesPerYear) ~= 1
        PerYearPerLTE = InputTreatments(ValuesPerYear,:);
        Tester = unique(PerYearPerLTE.TreatType);
        PerTreat = [];
        PerTreatNorm = [];
        for x = 1:length(Tester)
            ListAvg = find(strcmp(PerYearPerLTE.TreatType,Tester(x))==1);
            PerTreat(x) = nanmedian(PerYearPerLTE.Yield(ListAvg));
            % Normalised
            PerTreatNorm(x) = nanmedian(PerYearPerLTE.Yield(ListAvg)./ReferenceYield);
        end
        Yields.TreatType(i,1) = TreatTyper;
        Yields.Year(i,1) = ListerYears(i);
        Variety = unique(PerYearPerLTE.Variety);
        if length(Variety) >1
            disp('Multiple varieties in a year')
            disp(ListerYears(i))
        end
        if (sum(strcmp(Variety,'Unknown')) ~= 0) && (i ~= 1)
            Variety = Yields.Variety((i-1),1);
        end
        Yields.Variety(i,1) = Variety(1);
        Yields.Ctype(i,1) = PerYearPerLTE.CType(1);
        Yields.NrValuesYear(i,1) = length(PerTreat);
        Yields.MedianYield(i,1) = nanmedian(PerTreat); %#ok<*NANMEDIAN>
        Yields.MedianYieldCV(i,1) = nanstd(PerTreat)./nanmean(PerTreat); %#ok<*NANMEAN,*NANSTD>
        Yields.MedianYieldNorm(i,1) = nanmedian(PerTreatNorm);
        clear PerTreat PerTreatNorm
    else
        Yields.TreatType(i,1) = {'No values'};
        Yields.Year(i,1) = ListerYears(i);
        Yields.Ctype(i,1) = InputTreatments.CType(1);
        Yields.Variety(i,1) = {'No crop'};
        Yields.NrValuesYear(i,1) = NaN;
        Yields.MedianYield(i,1) = NaN;
        Yields.MedianYieldCV(i,1) = NaN;
        Yields.MedianYieldNorm(i,1) = NaN;
    end
end
end

% Code from Internet
% from https://nl.mathworks.com/matlabcentral/answers/95923-is-there-a-matlab-function-that-can-check-if-a-field-exists-in-a-matlab-structure
function isFieldResult = myIsField (inStruct, fieldName)
% inStruct is the name of the structure or an array of structures to search
% fieldName is the name of the field for which the function searches
isFieldResult = 0;
f = fieldnames(inStruct);
for i=1:length(f)
    if(strcmp(f{i},strtrim(fieldName)))
        isFieldResult = 1;
        return;
    elseif isstruct(inStruct.(f{i}))
        isFieldResult = myIsField(inStruct.(f{i}), fieldName);
        if isFieldResult
            return;
        end
    end
end
end


