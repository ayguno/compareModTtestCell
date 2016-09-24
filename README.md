# compareModTtestCell
Multiple experiments can be compared, experiments from different cell lines can be integrated

This function is designed to help us to interpret similarities between various treatments based on the proteins being co-effected by the individual treatments
Multiple TMT10 experiments can be compared, experiments from different cell lines can be integrated
"ModT test table results" printed by the Shiny server will be provided into the same working directory, along with the class vector files specific for each experiment (both in csv format)
The function will take a customly prepared csv file as input, which matches the data tables with class vector files
The function will generate heatmaps of top 10 up/down regulated proteins for each experiment
The function will also attempt to plot the entire set of experiments into a single scatterplot
ModT_ClassV_match (chr) = csv table that matches ModT table file names (first column), with the respective Class vector files (second column),experiment code (third column) and cell line the experiment was carried out (4th column) 
 keratins= T/F? (logical) False by default. If filtering for human Keratins, choose TRUE and also include "human_keratins_as_downloaded_from_HGNC_06242016.csv" in the working directory

 Functional as of 07/25/2016 (last update)
