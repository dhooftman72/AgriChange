function [CombinedYields,PerYearFunctions] = PerLTEcombine
load('Continents')
List = (Continents(:,1));

% Calculate yield arrays 
disp('Running per LTE arrays')
for i = 1: size(List,1) 
    disp(List(i))
    [~,~,~,~] = ExtractYieldData(char(List(i)));
end

warning off
for Var = 1:3
    Types = [{'MaxYieldTreat'},'AverageYield','MinimumYieldTreat'];
    disp('Running per Type combination')
    disp(Types(Var))
    % Combine in one file
    Type = Types(Var);
    for i = 1: length(List)
        %disp(List(i))
        File = (char(List(i)));
        Continent = Continents((find(strcmp(Continents(:,1),File)==1)),2);
        [Array,LTETesterSelect] = OpenArray(File,Type);
        if i == 1
            Combined = Array;
            PerLTEYearFunction = dataset({File},'varnames', {'LTE'}); %#ok<*DTSET>    
        else
            Combined = [Combined;Array]; %#ok<*AGROW>
        end
        PerLTEYearFunction.LTE(i,1) = {File};
        ToTest = FirstTraitTypeMention(Array);
        Tester = find(strcmp(LTETesterSelect.Treatment,ToTest.Treatment)==1);
        PerLTEYearFunction.TreatType(i,1) =ToTest.Treatment;
        PerLTEYearFunction.Crop(i,1) =ToTest.CType;
        PerLTEYearFunction.Continent(i,1) =Continent;
        PerLTEYearFunction.YearConstant(i,1) = LTETesterSelect.YearConstant(Tester);
        PerLTEYearFunction.YearCoef(i,1) = LTETesterSelect.YearCoef(Tester);
        PerLTEYearFunction.Year_FValue(i,1) = LTETesterSelect.Year_FValue(Tester);
        PerLTEYearFunction.Year_PValue(i,1) = LTETesterSelect.Year_PValue(Tester);
        PerLTEYearFunction.Years(i,1) = LTETesterSelect.Years(Tester);
        PerLTEYearFunction.StartYear(i,1) = LTETesterSelect.StartYear(Tester);
        PerLTEYearFunction.EndYear(i,1) = LTETesterSelect.EndYear(Tester);
        clear Array LTETesterSelect
    end
    disp('Adding Environmental info')
    load('Environment.mat');
    load('GeoLocationsStart');
    for i = 1:size(Combined,1)
        Year =  Combined.Year(i);
        YearChar = {['Y_', mat2str(Year)]};
        LTETest = Combined.LTE(i);
        if strcmp(LTETest,'ND2') == 1
            LTETest = 'D2';
        end
        if strcmp(LTETest,'M29') == 1
            LTETest = 'M209';
        end
        if strcmp(LTETest,'M27') == 1
            LTETest = 'M217';
        end
        save('all')
        Combined.Longitude_X(i,1) = double(GeoLocationsStart.Longitude_X(LTETest,1));
        Combined.Lattitude_Y(i,1) = double(GeoLocationsStart.Lattitude_Y(LTETest,1));
        Combined.Elevation(i,1) = double(GeoLocationsStart.Elevation(LTETest,1));
        Combined.Continent(i,1) = GeoLocationsStart.Continent(LTETest,1);
        Combined.SeasonMinTemp(i,1) = double(SeasonMinTemperatureYear(LTETest,YearChar));
        Combined.SeasonMaxTemp(i,1) = double(SeasonMaxTemperatureYear(LTETest,YearChar));
        Combined.SeasonDailyTemp(i,1) = double(SeasonDailyMeanTemperatureYear(LTETest,YearChar));
        Combined.YearMinTemp(i,1) = double(MinTemperatureYear(LTETest,YearChar));
        Combined.YearMaxTemp(i,1) = double(MaxTemperatureYear(LTETest,YearChar));
        Combined.YearDailyTemp(i,1) = double(DailyMeanTemperatureYear(LTETest,YearChar));
        Combined.DailyTempRange(i,1) = double(DailyTempRange(LTETest,YearChar));
        Combined.MovingDailyTempRange(i,1) = double(MovingDailyTempRange(LTETest,YearChar));
        Combined.AnnualPET(i,1) = double(AnnualPET(LTETest,YearChar));
        Combined.AnnualRainFall(i,1) = double(AnnualRainfall(LTETest,YearChar));
        Combined.MovingAnnualPET(i,1) = double(MovingAnnualPET(LTETest,YearChar));
        Combined.MovingAnnualRainFall(i,1) = double(MovingAnnualRainFall(LTETest,YearChar));
        Combined.CO2ManuaLoa(i,1) = double(CO2ManuaLoa(LTETest,YearChar)); 
        Combined.CWDperCallenderYear(i,1) = double(CWDperCallenderYear(LTETest,YearChar)); 
        Combined.PopPressure(i,1) = double(PopPressure(LTETest,YearChar)); 
        Combined.Radiation(i,1) = double(Radiation(LTETest,YearChar)); 
        Combined.Ozone(i,1) = double(Ozone(LTETest,YearChar)); 

    end
    clearvars -except Combined* List Var Types Per* Continents
    CombinedYields.(genvarname(char(Types(Var)))) = Combined;
    PerYearFunctions.(genvarname(char(Types(Var)))) = PerLTEYearFunction;
    save('CombinedYields','CombinedYields','PerYearFunctions')
end
end

function [Array,LTETesterSelect] = OpenArray(File,Type)
Filename = [File,'.mat'];
load(Filename)
Array = eval(char(Type));
end

function ToTest = FirstTraitTypeMention(Array)
count = 1;
ToTest = [];
while isempty(ToTest)==1
    if strcmp(Array.TreatType(count,1),'No values')~=1
        ToTest.Treatment = Array.TreatType(count,1);
        ToTest.CType = Array.Ctype(count,1);
    end
    count = count + 1;
end
end
