function [prob1, post1, post2] = posteriorProb2(model1_lambda, model2_lambda, testDat)
    post1 = cell2mat(arrayfun(@(a,b) poisspdf(a,b), round(testDat), model1_lambda, 'un', 0));
    post2 = cell2mat(arrayfun(@(a,b) poisspdf(a,b), round(testDat), model2_lambda, 'un', 0));
    prob1 = post1./(post1+post2);
end