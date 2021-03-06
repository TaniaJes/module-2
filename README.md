# Description taken from http://portals.broadinstitute.org/cgi-bin/cancer/publications/pub_paper.cgi?mode=view&paper_id=43"

Description of sample tables, datasets and prediction results for 

Golub et al "Molecular Classification of Cancer: Class Discovery and Class
Prediction by Gene Expression Monitoring"

  There are two datasets containing the initial (training, 38 samples)
and independent (test, 34 samples) datasets used in the paper.

A description of the samples can be found in table_ALL_AML_samples.rtf
(MS Word version) and table_ALL_AML_samples.txt (ASCII version)

The prediction results can be found in table_ALL_AML_predic.rft (MS
Word version) and table_ALL_AML_predic.txt (ASCII version).

These datasets contain measurements corresponding to ALL and AML
samples from Bone Marrow and Peripheral Blood. Details about the
experimental method and protocol can be found in
Experimental_protocol.html. Intensity values have been re-scaled such
that overall intensities for each chip are equivalent. This is done by
fitting a linear regression model using the intensities of all genes
with "P" (present) calls in both the first sample (baseline) and each
of the other samples. The inverse of the "slope" of the linear
regression line becomes the (multiplicative) re-scaling factor for the
current sample. This is done for every chip (sample) in the dataset
except the baseline which gets a re-scaling factor of one. File
table_ALL_AML_rfactors.txt contains a table with the chip re-scaling
factors.

  Initial (train) Dataset (38 samples)
 
  The dataset can be viewed as a text (tab delimited) file:
data_set_ALL_AML_train.txt, or loaded into your local disk or Excel
spreadsheet: data_set_ALL_AML_train.tsv

  Independent (test) Dataset (34 samples)

 
  The dataset can be viewed as a text (tab delimited) file:
data_set_ALL_AML_independent.txt, or loaded into your local disk or
Excel spreadsheet: data_set_ALL_AML_independent.tsv

