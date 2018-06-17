Project Title
This tool implements a training scheme for a classifier predicting the structural class of antibody complementary determining regions. However the training scheme can be used for training a Gradient Boosted Machine classifier for any already clustered sequences data. The implementation allows the comparison of the performance by the GBM method to that of the prediction method using BLAST to find a sequence of the highest bitscore and using its class as the predicted class of query. The tool generate figures enabling visualization of the accuracies and results.. 

Getting Started
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. See deployment for notes on how to deploy the project on a live system.

Prerequisites
R(figure out the least version required) installed. 
What things you need to install the software and how to install them

Give examples
Installing
# run "Rscript trainingGBM.R  --trainingT training_table.tbl   

#data.table should have a "sequence", "class", "sequence_type", "sequence_length". 

Rscript GBM_model_building.R -input_sequence_table  data1.table
Rscript BLAST_result.R -input_sequence_table  data2.table
Rscript demonstrating_results.R 
cd ./Plots 


Give the example

until finished
End with an example of getting some data out of the system or using it for a little demo

Running the tests
Explain how to run the automated tests for this system

Break down into end to end tests
Explain what these tests test and why

Give an example
And coding style tests
Explain what these tests test and why

Give an example
Deployment
Add additional notes about how to deploy this on a live system

Built With
Dropwizard - The web framework used
Maven - Dependency Management
ROME - Used to generate RSS Feeds
Contributing
Please read CONTRIBUTING.md for details on our code of conduct, and the process for submitting pull requests to us.

Versioning
We use SemVer for versioning. For the versions available, see the tags on this repository.

Authors
Billie Thompson - Initial work - PurpleBooth
See also the list of contributors who participated in this project.

License
This project is licensed under the MIT License - see the LICENSE.md file for details

Acknowledgments
Hat tip to anyone whose code was used
Inspiration
etc
