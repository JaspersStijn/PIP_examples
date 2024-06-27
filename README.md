This repository contains R code related to the probability of improved prediction (PIP), which is a newly introduced statistical concept for model comparison. In general, the PIP is a probabilistic measure for directly comparing two competing models. Based on a user-defined loss function, this is achieved by quantifying the frequency of instances 
that one model gives better predictions than the other model. 

Two papers were prepared, for which this repository contains the R code related to their respective data examples.

--> Paper 1: 'The Probability of Improved Prediction: a new concept in statistical inference'  (available on arXiv: https://arxiv.org/abs/2405.17064)
            
            -- R file "Functions_paper1.R" contains the required methodology for determinining the PIP using (repeated) k-fold cross-validation. 
            -- R file "Reproducibility papers.R" shows an application and outputs Table 4 of the paper


--> Paper 2: 'The Probability of Improved Prediction as a New Concept for Model Selection in the Presence of Outliers' (submitted)
            -- R file "Functions_paper2.R" contains the required methodology for performing model comparison using the PIP, calculated using m out of n bootstrapping.
            -- R file "Model Comparison Data Examples.R" shows two data applications:
                    * Stackloss data (model building using robust regression models with three covariates) 
                    * Ozone data (model building with gradient boosting machine)
            -- R file "EMS presentation code.R" provides graphs and summaries that were shown during the European Meeting of Statisticians (2023). Presentation abailable in presentations folder.

            
                            
