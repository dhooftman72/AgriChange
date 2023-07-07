function AgrichangeStatistcalAnalyses_July_2023
warning off
load('CombinedYields.mat') %#ok<*LOAD>
for Var = 1:1:3

    Types = [{'MaxYieldTreat'},'AverageYield','MinimumYieldTreat'];
    OutputFileName = ['Outputs_',char(Types(Var))];
    disp(Types(Var))

    WhichYield = 'MedianYieldNorm' ;

    ArrayToUse = CombinedYields.(genvarname(char(Types(Var))));
    PerYearFunctionToUse = PerYearFunctions.(genvarname(char(Types(Var)))); 
    Tester = unique(ArrayToUse.FullLTE);
    MedianperFullLTE = [];
    for i = 1:length(Tester)
        Values = find(strcmp(ArrayToUse.FullLTE,Tester(i))==1);
        MedianperFullLTE.Yield(i,1) = nanmedian(ArrayToUse.MedianYield(Values)); %#ok<*NANMEDIAN>
        MedianperFullLTE.Xcoor(i,1) = ArrayToUse.Longitude_X(Values(1));
        MedianperFullLTE.Observedoor(i,1) = ArrayToUse.Lattitude_Y(Values(1));
    end
    AutoCorrelationCoef = SpatialAutoCorrelation2023(MedianperFullLTE.Xcoor,MedianperFullLTE.Observedoor,MedianperFullLTE.Yield);
    for i = 1:length(Tester)
        Values = find(strcmp(ArrayToUse.FullLTE,Tester(i))==1);
        ArrayToUse.AutoCor(Values,1) = AutoCorrelationCoef(i);
        WhichOne = find(strcmp(PerYearFunctionToUse.LTE,Tester(i))==1);
        ArrayToUse.YearsRunning(Values,1) = PerYearFunctionToUse.Years(WhichOne);
    end
    IsnanList = find(isnan(ArrayToUse.(genvarname(WhichYield)))==1); %#ok<COMPNOP>
    ArrayToUse(IsnanList,:) = []; %#ok<FNDSB>
    IsnanList = find(ArrayToUse.(genvarname(WhichYield))<=0);
    ArrayToUse(IsnanList,:) = []; %#ok<FNDSB>

    Parameters = [{'TempratureRange'},'SeasonMinTemp','SeasonMaxTemp','Radiation','CO2ManuaLoa',...
        'Changes in PET','Changes in Rainfall','CWDperCallenderYear','PopPressure','Ozone','Years'];

    [NormalisationArray,Tmp] = NormalisationArrayFunc(ArrayToUse);
    DailyTemp = Tmp.Daily./max(Tmp.Daily);
    MinTemp = Tmp.Min./max(Tmp.Min);
    MaxTemp = Tmp.Max./max(Tmp.Max);
    Radiation = Tmp.Radiation./max(Tmp.Radiation);
    CO2 = Tmp.CO2./max(Tmp.CO2);
    Pet = Tmp.Pet./max(Tmp.Pet);
    Rain = Tmp.Rain./max(Tmp.Rain);
    CWD = Tmp.CWD./max(Tmp.CWD);
    PopulationPres = Tmp.Pop./max(Tmp.Pop);
    Ozone = Tmp.Ozone./max(Tmp.Ozone); % Don't TRUST THIS ONE, Locations west and south, minus sign didn't seem to work and probably high atmosphere
    Years = Tmp.Year ./max(Tmp.Year);

    EnvironmentalArray = [DailyTemp,MinTemp,MaxTemp,Radiation,CO2,Pet,Rain,CWD,PopulationPres, Ozone];
    EnvironmentalArrayPlus = [DailyTemp,MinTemp,MaxTemp,Radiation,CO2,Pet,Rain,CWD,PopulationPres,Ozone, Years];

    AutoCorrelation = ArrayToUse.AutoCor;
    LTE = cellstr(ArrayToUse.LTE);
    SubLTE = cellstr(ArrayToUse.SubLTE);
    Variety = cellstr(ArrayToUse.Variety);

    Yield = ArrayToUse.(genvarname(WhichYield));
    CropType = cellstr(ArrayToUse.Ctype);

    Weights = ArrayToUse.YearsRunning;

    % Model for Year series
    [~,PerYearRegresStats.Output,PerYearRegresStats.Stats]= anovan((Yield),{AutoCorrelation, Years,CropType, LTE, SubLTE, Variety},'sstype',1,...
        'model',[1 0 0 0 0 0; 0 1 0 0 0 0; 0 0 1 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1],'nested',[0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 0 0 0; 0 0 0 1 0 0; 0 0 0 0 1 0],'display', 'off',...
        'continuous',[1 2], 'random', [3 4 5 6],'varnames', {'AutoCorrelation','Year', 'CType', 'LTE' , 'SubLTE', 'Variety'});

    % REMOVING THE NOISE: creating marginla values without variables of
    % interest
    [~,MarginalStats.Output,MarginalStats.Stats]= anovan((Yield),{AutoCorrelation, CropType, LTE, SubLTE, Variety},'sstype',1,...
        'model',[1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1],'nested',[0 0 0 0 0; 0 0 0 0 0; 0 0 0 0 0; 0 0 1 0 0; 0 0 0 1 0],'display', 'off',...
        'continuous',[1], 'random', [2 3 4 5],'varnames', {'AutoCorrelation', 'CType', 'LTE' , 'SubLTE', 'Variety'});


    %SecondTierVar = statsFirst.resid';
    MarginalValues = Yield-(MarginalStats.Stats.resid);
    ArrayToUse.MarginalValues = MarginalValues;
    %AV_PerEnvironment = dataset({'dummy'}, 'varnames',{'Parameters'});
    NLM_PerEnvironment = dataset({'dummy'}, 'varnames',{'Parameters'});
    Year_PerEnvironmentCorrelation = dataset({'dummy'}, 'varnames',{'Parameters'});
    for Env = 1:1:length(Parameters)
        if Env > size(EnvironmentalArray,2)
            EnvironVar = Years;
        else
            EnvironVar =EnvironmentalArray(:,Env);
        end
        % [~,StatsOutputSecond,statsSecond]= anovan((MarginalValues),{Years,EnvironVar},'sstype',3,...
        %     'model',[1 0; 0 1; 1 1],'display', 'off','continuous',[1 2], 'varnames', {'Years', 'Environment'});
        % AV_PerEnvironment.Parameters(Env,1) = Parameters(Env);
        % AV_PerEnvironment.YearF(Env,1) = StatsOutputSecond(2,6);
        % AV_PerEnvironment.YearP(Env,1) = StatsOutputSecond(2,7);
        % AV_PerEnvironment.EnvironmentF(Env,1) = StatsOutputSecond(3,6);
        % AV_PerEnvironment.EnvironmentP(Env,1) = StatsOutputSecond(3,7);
        % AV_PerEnvironment.InteractionF(Env,1) = StatsOutputSecond(4,6);
        % AV_PerEnvironment.InteractionP(Env,1) = StatsOutputSecond(4,7);
        % AV_PerEnvironment.YearCoeff(Env,1) = statsSecond.coeffs(2,1);
        % AV_PerEnvironment.EnvCoeff(Env,1) = statsSecond.coeffs(3,1);
        % AV_PerEnvironment.InteractionCoeff(Env,1) = statsSecond.coeffs(4,1);
        % AV_PerEnvironment.Constant(Env,1) = statsSecond.coeffs(1,1);
        % clear StatsOutputSecond statsSecond

        wnlm = nlmRegression(EnvironVar,MarginalValues,Weights);
        NLM_PerEnvironment.Parameters(Env,1) = Parameters(Env);
        NLM_PerEnvironment.MinEnvValue(Env,1) = min(EnvironVar);
        NLM_PerEnvironment.MaxEnvValue(Env,1) = max(EnvironVar);
        NLM_PerEnvironment.Tstat(Env,1) = table2array(wnlm.Coefficients(2,3));
        NLM_PerEnvironment.PValue(Env,1) = table2array(wnlm.Coefficients(2,4));
        NLM_PerEnvironment.Rsquared(Env,1) = wnlm.Rsquared.Ordinary(1,1);
        NLM_PerEnvironment.LogLikelihood(Env,1) = wnlm.LogLikelihood(1,1);
        NLM_PerEnvironment.Coefficient(Env,1) = table2array(wnlm.Coefficients(2,1));
        NLM_PerEnvironment.CoefficientSE(Env,1) = table2array(wnlm.Coefficients(2,2));
        NLM_PerEnvironment.Constant(Env,1) = table2array(wnlm.Coefficients(1,1));
        clear wnlm
        for cropType = 1:1:2
            if cropType == 1
                OnetoFind = 'C3';
            else
                OnetoFind = 'C4';
            end
            List = find(strcmp(CropType,OnetoFind)==1);
            EnvironType = EnvironVar(List);
            MarginalType = MarginalValues(List);
            WeightType = Weights(List);
            wnlm = nlmRegression(EnvironType,MarginalType,WeightType);
            NLM_PerEnvironment.(genvarname(char({[OnetoFind,'_Coefficient',]})))(Env,1) = table2array(wnlm.Coefficients(2,1));
            NLM_PerEnvironment.(genvarname(char({[OnetoFind,'_Constant',]})))(Env,1) = table2array(wnlm.Coefficients(1,1));
            NLM_PerEnvironment.(genvarname(char({[OnetoFind,'_PValue',]})))(Env,1) = table2array(wnlm.Coefficients(2,4));
            NLM_PerEnvironment.(genvarname(char({[OnetoFind,'_MinValue',]})))(Env,1)= min(EnvironType);
            NLM_PerEnvironment.(genvarname(char({[OnetoFind,'_MaxValue',]})))(Env,1)= max(EnvironType);
        end
        clear list
        wnlm = nlmRegression(Years,EnvironVar,Weights);
        Year_PerEnvironmentCorrelation.Parameters(Env,1) = Parameters(Env);
        Year_PerEnvironmentCorrelation.Tstat(Env,1) = table2array(wnlm.Coefficients(2,3));
        Year_PerEnvironmentCorrelation.PValue(Env,1) = table2array(wnlm.Coefficients(2,4));
        Year_PerEnvironmentCorrelation.Rsquared(Env,1) = wnlm.Rsquared.Ordinary(1,1);
        Year_PerEnvironmentCorrelation.LogLikelihood(Env,1) = wnlm.LogLikelihood(1,1);
        Year_PerEnvironmentCorrelation.Coefficient(Env,1) = table2array(wnlm.Coefficients(2,1));
        Year_PerEnvironmentCorrelation.CoefficientSE(Env,1) = table2array(wnlm.Coefficients(2,2));
        clear wnlm
    end
    % per croptype



    %%
    disp('PCA')

    % PLS regression Full data

    [NLMOuts,PLS] = PLSandNLMRegres(EnvironmentalArray,MarginalValues,Weights);
    [PCAOverview,PCA_perEnv.FullAxis1,EstimatedValues] = WhatToOutputPCA(NLMOuts,PLS,Parameters(1:(length(Parameters)-1)),{'Full_Axis'},1);
    [~,PCA_perEnv.FullAxis2,~] = WhatToOutputPCA(NLMOuts,PLS,Parameters(1:(length(Parameters)-1)),{'Full_Axis'},2);
    ArrayToUse.PCA_XFull = EstimatedValues.X;
    ArrayToUse.PCA_YFull = EstimatedValues.X;
    clear NLMOuts PLS EstimatedValues

    for AxNum = 1:2
        % Make loop of take one out each, than two etc. Record highest R2 and
        % store.
        ShortendEnvArray = EnvironmentalArray;
        LogLikeStore = PCAOverview.LogLikelihood(1,1);
        ParametersTmp = Parameters(1:(length(Parameters)-1));
        StopLike = 0;

        loop = 0;
        while StopLike == 0
            loop = loop + 1;
            clear TmpShort
            for Env = 1:size(ShortendEnvArray,2)
        	    List = 1:size(ShortendEnvArray,2);
                List(Env) = [];
                Tmp_Array = ShortendEnvArray(:,List);
                clear NLMOutsLoop 
                [NLMOutsLoop,~] = PLSandNLMRegres(Tmp_Array,MarginalValues,Weights);
                TmpShort.Parameters(Env,1) = ParametersTmp(Env);
                TmpShort.Rsquared(Env,1) = NLMOutsLoop.Rsquared.Ordinary(1,1);
                TmpShort.LogLikelihood(Env,1) = NLMOutsLoop.LogLikelihood(1,1);
            end
            MaxLogLike = find(TmpShort.LogLikelihood==max(TmpShort.LogLikelihood));
            if TmpShort.LogLikelihood(MaxLogLike) > LogLikeStore
                ShortendEnvArray(:, MaxLogLike) = [];
                ParametersTmp( MaxLogLike) = [];
                Message = ['Removed is ',char(TmpShort.Parameters( MaxLogLike))];
                disp(Message)
                Message = ['Increase from ' num2str(LogLikeStore) ' to ' num2str(TmpShort.LogLikelihood(MaxLogLike))];
                disp(Message)
                LogLikeStore = TmpShort.LogLikelihood( MaxLogLike);       
                save ('all')
            else
                StopLike = 1;
            end
            [NLMOutsLoop,PLSLoop] = PLSandNLMRegres(ShortendEnvArray,MarginalValues,Weights);

        end % while
        [NLM_Clean,PCA_perEnv.(genvarname(char({['CleanedAxis',num2str(AxNum)]}))),EstimatedValues] = WhatToOutputPCA(NLMOutsLoop,PLSLoop,ParametersTmp,{['CleanedAxis',num2str(AxNum)]},AxNum);
        if AxNum == 1
            ArrayToUse.PCA_XSelected = EstimatedValues.X;
            ArrayToUse.PCA_YSelected = EstimatedValues.X;
        end
        PCAOverview(2,:) = NLM_Clean(1,:);
        PCAOverview.Parameters(2,1) = {'CleanedAxis'};
    end
    PerLTETimesSeries = PerYearFunctions;
    Rest.Year_PerEnvironmentCorrelation = Year_PerEnvironmentCorrelation;
    Rest.EnvironmentalArrayPlus = EnvironmentalArrayPlus;
    Rest.PerYearRegresStats = PerYearRegresStats;
    Rest.NormalisationArray = NormalisationArray;
    Rest.MarginalStats = MarginalStats;
    
    save(OutputFileName, 'Rest', 'PerLTETimesSeries', 'NLM_PerEnvironment','PCAOverview','ArrayToUse','PCA_perEnv');
end % for crop yield types
end % End of main function


function [VarOut,PLS] = PLSandNLMRegres(X_array,Y_values,Weights)
[~,~,PLS.Axis,PLS.YAxis,PLS.BETA,PLS.PCTVAR,~,stats] = plsregress(X_array,Y_values,3);
PLS.Weights = stats.W;
VarOut = nlmRegression(PLS.Axis(:,1),Y_values,Weights);
end

function VarOut = nlmRegression(XVar,YVar,Weights)
modelFun = @(b,x) b(1)+ b(2).*x;
start = [1; 0];
VarOut = fitnlm(XVar,YVar,modelFun,start,'Weight',Weights);
end

function [OutOverall,OutParameter,EstimatedValues] = WhatToOutputPCA(NLMOuts,PLSOuts,Parameters,FullAxName,AxNum)
EstimatedValues.X = PLSOuts.Axis(:,1);
EstimatedValues.Y = PLSOuts.YAxis(:,1);
OutOverall = dataset(FullAxName, 'varnames',{'Parameters'});
OutOverall.Tstat(1,1) = table2array(NLMOuts.Coefficients(2,3));
OutOverall.PValue(1,1) = table2array(NLMOuts.Coefficients(2,4));
OutOverall.ExplainedAxis1X(1,1) = PLSOuts.PCTVAR(1,1);
OutOverall.ExplainedAxis1Y(1,1) = PLSOuts.PCTVAR(2,1);
OutOverall.ExplainedAxis2X(1,1) = PLSOuts.PCTVAR(1,2);
OutOverall.ExplainedAxis2Y(1,1) =  PLSOuts.PCTVAR(2,2);
OutOverall.Rsquared(1,1) = NLMOuts.Rsquared.Ordinary(1,1);
OutOverall.LogLikelihood(1,1) = NLMOuts.LogLikelihood(1,1);
OutOverall.Coefficient(1,1) = table2array(NLMOuts.Coefficients(2,1));
OutOverall.CoefficientSE(1,1) = table2array(NLMOuts.Coefficients(2,2));
OutOverall.Constant(1,1) = table2array(NLMOuts.Coefficients(1,1));
OutOverall.ConstantSE(1,1) = table2array(NLMOuts.Coefficients(1,2));

% [~,StatsOutputSecond,statsSecond]= anovan((MarginalValues),{Years,PLSOuts.Axis(:,AxNum)},'sstype',3,...
%     'model',[1 0; 0 1; 1 1],'display', 'off','continuous',[1 2], 'varnames', {'Years', 'Axis'});
% OutOverall.YearF(1,1) = StatsOutputSecond(2,6);
% OutOverall.YearP(1,1) = StatsOutputSecond(2,7);
% OutOverall.EnvironmentF(1,1) = StatsOutputSecond(3,6);
% OutOverall.EnvironmentP(1,1) = StatsOutputSecond(3,7);
% OutOverall.InteractionF(1,1) = StatsOutputSecond(4,6);
% OutOverall.InteractionP(1,1) = StatsOutputSecond(4,7);
% OutOverall.YearCoeff(1,1) = statsSecond.coeffs(2,1);
% OutOverall.EnvCoeff(1,1) = statsSecond.coeffs(3,1);

OutParameter = dataset((Parameters(1:length(Parameters)))', 'varnames',{'Parameters'});
OutParameter.Weightfactors = PLSOuts.Weights(:,AxNum);
absWeights = abs(OutParameter.Weightfactors);
SumWeights = sum(absWeights);
OutParameter.RelativeWeights = absWeights./SumWeights;
OutParameter.Betas = PLSOuts.BETA(2:end,1);
end

function  [NormalisationArray,Tmp] = NormalisationArrayFunc(Array)
    NormalisationArray.Max = dataset(max(Array.MovingDailyTempRange),'varnames',{'DailyTemp'});
    NormalisationArray.Min = dataset(min(Array.MovingDailyTempRange),'varnames',{'DailyTemp'});

    NormalisationArray.Max.MinTemp = max(Array.SeasonMinTemp);
    NormalisationArray.Max.MaxTemp = max(Array.SeasonMaxTemp);
    NormalisationArray.Max.Radiation = max(Array.Radiation);
    NormalisationArray.Max.CO2ManuaLoa = max(Array.CO2ManuaLoa);
    NormalisationArray.Max.PET = max(Array.MovingAnnualPET);
    NormalisationArray.Max.Rain = max(Array.MovingAnnualRainFall);
    NormalisationArray.Max.CWD = max(Array.CWDperCallenderYear);
    NormalisationArray.Max.PopPressure = max(Array.PopPressure);
    NormalisationArray.Max.Ozone = max(Array.Ozone);
    NormalisationArray.Max.Years =max(Array.Year);
    NormalisationArray.Min.MinTemp = min(Array.SeasonMinTemp);
    NormalisationArray.Min.MaxTemp = min(Array.SeasonMaxTemp);
    NormalisationArray.Min.Radiation = min(Array.Radiation);
    NormalisationArray.Min.CO2ManuaLoa = min(Array.CO2ManuaLoa);
    NormalisationArray.Min.PET = min(Array.MovingAnnualPET);
    NormalisationArray.Min.Rain = min(Array.MovingAnnualRainFall);
    NormalisationArray.Min.CWD = min(Array.CWDperCallenderYear);
    NormalisationArray.Min.PopPressure = min(Array.PopPressure);
    NormalisationArray.Min.Ozone = min(Array.Ozone);
    NormalisationArray.Min.Years =min(Array.Year);

    Tmp.Daily = Array.MovingDailyTempRange - NormalisationArray.Min.DailyTemp;
    Tmp.Min = Array.SeasonMinTemp - NormalisationArray.Min.MinTemp;
    Tmp.Max = Array.SeasonMaxTemp - NormalisationArray.Min.MaxTemp;
    Tmp.Radiation = Array.Radiation - NormalisationArray.Min.Radiation;
    Tmp.CO2 = Array.CO2ManuaLoa - NormalisationArray.Min.CO2ManuaLoa;
    Tmp.Pet = Array.MovingAnnualPET - NormalisationArray.Min.PET;
    Tmp.Rain = Array.MovingAnnualRainFall - NormalisationArray.Min.Rain;
    Tmp.CWD = Array.CWDperCallenderYear - NormalisationArray.Min.CWD;
    Tmp.Pop = Array.PopPressure - NormalisationArray.Min.PopPressure;
    Tmp.Ozone = Array.Ozone - NormalisationArray.Min.Ozone;
    Tmp.Year = (Array.Year)- NormalisationArray.Min.Years;
end

function AutoCorrelationCoef = SpatialAutoCorrelation2023(XCoor,YCoor,VarIn)
Size = length(XCoor);
btsfact = 5; %5 degrees maximum interaction, else autocorrelation coefficient = 0
[wij] = Morans(XCoor,YCoor,5,(Size.*btsfact));
save('wij.mat','wij');

AutoCorrelationCoef = AutoCorFunc(VarIn,Size);
end

function AutoCorrelationCoef = AutoCorFunc(InVariation,Size)
wij = load('wij.mat');
parfor i = 1:Size
    Autot = 0;
    for j = 1:Size
        if i ~= j
            if isnan(InVariation(j))~= 1
                Autot = Autot + (wij.wij(i,j).*InVariation(j));  %#ok<*NODEF,*PFIIN>
            end
        end
    end
    if Autot > 0
        AutoCorrelationCoef(i,1) = Autot./nansum(wij.wij(i,:)); %#ok<*PFBNS,*AGROW>
    else
        AutoCorrelationCoef(i,1) = 0;
    end
end
end

function [wij] = Morans(X,Y,MaxD,bts) %#ok<INUSD>
if length(X) ~= length(Y)
    display('Fatal Error, unequal grid size')
    return
end
Size = length(X);
if Size > 1000
    poolobj = gcp('nocreate');
    if isempty(poolobj) == 1
        parpool 'Processes'
    end
end
parfor i = 1:Size
    for j = 1:Size
        if i == j
            wij(i,j) = 0;  %#ok<*AGROW>
        else
            Dist =  log10((sqrt ((((X(i)-X(j))^2) + (Y(i)-Y(j))^2)))+1);
            wijtmp = (log10(MaxD)-Dist)./log10(MaxD);
            wijtmp(wijtmp<0) = 0;
            wij(i,j) = wijtmp;
        end
    end
end
poolobj = gcp('nocreate');
if isempty(poolobj)~= 0
    delete(poolobj);
end
clear poolobj
end
