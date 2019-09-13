The scripts in this folder were used to generate the figures from the data that was analyzed in matlab ("data.csv").

When using the scripts with a different datatable, either keep the format of 11 columns as metadata, or change the code 
at the following command: gather(feature, value, 11:ncol(data), na.rm = TRUE) %>%
In some cases, the name of the last column of the metadata ("Section") was used: gather(feature, value, (1+which(colnames(data) == "Section" )):ncol(data), na.rm=TRUE)
Change it according to the name of the last column. 