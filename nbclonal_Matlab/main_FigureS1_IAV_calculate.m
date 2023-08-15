function void = main_FigureS1_IAV_calculate(void)

clear all; close all; clc;

for thresh = 1:3
    switch thresh
        case 1
            infile = 'figure_IAV_T005'; load(infile);
        case 2
            infile = 'figure_IAV_T03'; load(infile);
        case 3
            infile = 'figure_IAV_T07'; load(infile);
    end

    emp_clonal_list = 0:max(results_emp.nclonal);
    for n_clonal = emp_clonal_list
        emp_n(n_clonal+1) = length(find(results_emp.nclonal == n_clonal));
    end

    for R0val = 1:3

        switch R0val
            case 1
                params.R0 = 4.4;
            case 2
                params.R0 = 11.1;
            case 3
                params.R0 = 37.7;
        end
        params.mulist = 0.01:0.01:4;
        params.lambdalist = 0.01:0.02:4;

        % since can't go up to Inf, use these max values:
        params.maxFinalSize = 50;
        params.maxGeneratedMutantLineages = 50;
        params.maxEstablishedMutantLineages = 50;
        params.maxClonal = 10;

        outfileS1_IAV = strcat('figureS1_IAV_results_thres', int2str(thresh), '_R0', int2str(R0val));

        cntrmu = 1;
        for mu = params.mulist
            mu
            params.mu = mu;
            N0_matrix = [];
            cntrN = 1; params.N0list = 1:20;
            for N0 = params.N0list
                params.N = N0;
                [n_clonals, rho] = calculateClonalPMF(params, 0);
                pmf_matrix(cntrN,:) = rho;
                cntrN = cntrN + 1;
            end
            cntrlambda = 1;
            for lambda = params.lambdalist
                logL_val = GetLogL(lambda, params.N0list, n_clonals, pmf_matrix, results_emp.nclonal, params.R0);
                logL_val
                results.logLMatrix(cntrmu, cntrlambda) = logL_val;
                cntrlambda = cntrlambda + 1;
            end
            cntrmu = cntrmu + 1;
        end

        save(outfileS1_IAV, 'params', 'results');

        [params.meanNlist, params.meanNblist] = ConvertLambdaListToMeanNandMeanNbList(params.lambdalist, params.R0);

        save(outfileS1_IAV, 'params', 'results');
    end
end
