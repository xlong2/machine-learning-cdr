Project Title <br />
This tool implements a training scheme for a classifier predicting the structural class of antibody complementary determining regions. However the training scheme can be used for training a Gradient Boosted Machine classifier for any already clustered sequences data. The implementation allows the comparison of the performance by the GBM method to that of the prediction method using BLAST to find a sequence of the highest bitscore and using its class as the predicted class of query. The tool generate figures enabling visualization of the accuracies and results.. 

File list:<br />
Major scripts:<br />
	trainingGBM.R<br />  
	trainingBLAST.R <br />
	making_plots.R<br />
Functions files:<br />
	functions.R<br />
	utility_function.R<br />


Getting Started <br />
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.<br />
1. Dowload the package into your local directory. <br />
2. Enter machine_learning_cdr<br />
3. Run # Rscript trainingBLAST.R --data training_table.tsv , which will install all dependent R libraries. The script can be terminated once all libraries are installed without finishing the model building.  <br />


Prerequisites <br />
R-3.4 or newer version installed. <br />
Dependent R packages:<br />
 "Matrix", "grid","caret","RWebLogo","MLmetrics","parallel","pryr","protr", "gbm","ggplot2","ggrepel","reshape2","gridExtra","doMC","RColorBrewer","scales","e1071"  <br />


Help <br />
## By invoking any of the three script "trainingGBM.R", "trainingBLAST.R" and "making_plots.R", it will install all dependent packages and take a tab separated table as input data for training. More options for passing arguments to any of the scripts specific to training and data directory can be specified. Help can be viewed by involking "# Rscript trainingGBM.R  --help"  <br />


Test examlple <br />

Example 1: Automated test <br />
In the test case data training_table.tsv, CDR sequences are their corresponding cluster information are recorde in the tsv file and used for training the prediction models. By default, a minimal 5-trees GBM is generated for each modeling building. The accuracy of the GBM method and BLAST method are each obtained by 3 repeats 10 folds cross validation.  
# Rscript trainingGBM.R  --data training_table.tsv
# Rscript trainingBLAST.R  --data training_table.tsv
# Rscript making_plots.R  
"training_table.tsv" should have a "sequence", "class", "sequence_type", "sequence_length". A single class should have at least 8 cases. 
Results will be in default Data_dir and Plot_dir. <br /> 
For visualizing the accuracies of the two prediction method BLAST and GBM, please refer to the resultant figure in Plot_dir. <br />

Example 2:<br />
This example pass arguments to scripts to specify the number of repeats and folds in the x repeats n folds cross validation scheme, and also the results data directory and the plots directory. The GBM modeling building will encompass grid searching for the number of trees of 5, 10 and 15, also for the complexity of each tree of 3 and 6 branches. The arguments passed to different scripts should not contradict in a single test. For example, the number of repeats and the data saving directory should be the same for all three scripts.  <br />   
# Rscript trainingGBM.R  --data training_table.tsv  --n_repeats 3 --n_folds 10  --output_dir New_dir --n_trees 5:10:15 --complexity 3:6<br />
# Rscript trainingBLAST.R  --data training_table.tsv  --n_repeats 3 --n_folds 10  --output_dir New_dir <br />
# Rscript making_plots.R --n_repeats 3  --plot_dir  New_plot <br />
#data.table should have a "sequence", "class", "sequence_type", "sequence_length". A single class should have at least 8 cases. <br />
Results will be in specified New_dir and New_plot. <br />
For visualizing the accuracies of the two prediction method BLAST and GBM and error cases comparison, please refer to the resultant figure in New_plot. <br />


Running the method <br />
The script trainingGBM.R and trainingBLAST.R can be run independently. The script making_plots.R rely on the outs from the other two scripts therefore should only be run after running the trainingGBM.R and trainingBLAST.R scripts.  <br />
If the data size is large and grid searching is extensive. The table recording the data to be trained can be splitted so long as all the data needed for generating a single model is in the same table. The grid searching parameters can be further divided into smaller taskes for each running script. For instance, number of trees searching from 100 to 2000 in 200 increments can be separated to 100 to 1000, then 1000 to 2000 for runnning the scripts in different cores.   <br/ >
Generally, the prediction accuracy will increase as the n_trees and n_complexity increases up to a point. It is generally suggested to keep the value of n_complexity between 3 and 10. <br/ >

Specific project: <br/ >
The tool is a simplified and beautified version of the method used in my paper in the process of submission: "Non-H3 CDR template selection in antibody modeling through machine learning". The original implementation, data and plots can be found in cdr_classifier_result, cdr_classifier_result/Data_processed and cdr_classifier_result/Plots. <br/ >


Author <br/ >
Xiyao Long  <xlong2@jhu.edu>; <bibilong1111@gmail.com>  <br/ >
See also the list of contributors who participated in this project. <br/ >

License <br/ >
This project is licensed under the MIT License - see the LICENSE.md file for details. <br/ >



Inspiration <br/ >
The tool is inspired from a review(1) comparing different machine learning methods building sequence based classifiers and concluding the Gradient Boosted Machine gives the best model accuracy. The implementation heavily depend on the GBM implementation of "gbm" and "caret" packages but circumvent some problems when I try to perform extensive grid searching in a clustering environment. This implementation allows parallelization of different datasets by initiating the scripts with multiple data sets and enable specified multiple cores for a single script running instance. And it save the x-repeats-n-folds model estimation for each indiviual set of model parameter which allows results curation from multiple running scripts.<br/ > 

The model training and testing scheme is employed for testing the data from PyIgClassify(2), which is a antibody CDR database with each CDR sequence assigned to a specific structural class within a specific loop type and length. The dataset is chosen for the study in the hope of improving the accuracy of selecting a good structure template for a CDR sequence during the antibody structure modeling method implemented in RosettaAntibody(3). <br/ >


1. Jain P., Garibaldi JM., Hirst JD. 2009. Supervised machine learning algorithms for protein structure classification. Computational Biology and Chemistry 33:216-223. DOI: 10.1016/j.compbiolchem.2009.04.004.
2. Adolf-Bryfogle J., Xu Q., North B. et al. 2015. PyIgClassify: a database of antibody CDR structural classifications. Nucleic Acids Research 43:D432-D438. DOI: 10.1093/nar/gku1106.
3. Weitzner BD., Kuroda D., Marze N. et al. 2014. Blind prediction performance of RosettaAntibody 3.0: Grafting, relaxation, kinematic loop modeling, and full CDR optimization. Proteins: Structure, Function and Bioinformatics 82:1611-1623. DOI: 10.1002/prot.24534.

