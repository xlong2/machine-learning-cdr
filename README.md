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
If the data size is large and grid searching is extensive. The table recording the data to be trained can be splitted so long as all the data needed for generating a single model is in the same table. The grid searching parameters can be further divided into smaller taskes for each running script. For instance, number of trees searching from 100 to 2000 in 200 increments can be separated to 100 to 1000, then 1000 to 2000 for runnning the scripts in different cores. 
Generally, the prediction accuracy will increase as the n_trees and n_complexity increases up to a point. It is generally suggested to keep the value of n_complexity between 3 and 10. 

Specific project:
The tool is a simplified and beautified version of the method used in my paper in the process of submission: "Non-H3 CDR template selection in antibody modeling through machine learning". The original implementation, data and plots can be found in cdr_classifier_result, cdr_classifier_result/Data_processed and cdr_classifier_result/Plots.


Author
Xiyao Long  xlong2@jhu.edu; bibilong1111@gmail.com
See also the list of contributors who participated in this project.

License
This project is licensed under the MIT License - see the LICENSE.md file for details

Acknowledgments


Inspiration




