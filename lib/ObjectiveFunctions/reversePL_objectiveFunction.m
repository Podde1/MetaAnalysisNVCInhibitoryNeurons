function[ v ] = reversePL_objectiveFunction(theta, modelName, thetaConstants, OptParamBool, circuitConstants, DataStructure,ExperimentIdx,variableNames,modelObservablesIndex, parameterMappingTable, modelPropertyNames, options, parameterIdx, polarity, threshold)
fval = objectiveFunction(theta, modelName, thetaConstants, OptParamBool, circuitConstants, DataStructure, ExperimentIdx, variableNames ,modelObservablesIndex, parameterMappingTable ,modelPropertyNames,options, 0,[]);

v = polarity*theta(parameterIdx); % return the parameter value at index parameterIdx. multiply by polarity = {-1,1} to swap between finding maximum and minimum parameter value. 

if fval>threshold % check if solution is under limit 
    v = abs(v) + (fval-threshold).^2; % add penelty if the solution is above the limit. penelty grows the more over the limit the solution is.  
end

end 