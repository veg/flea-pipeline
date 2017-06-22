LoadFunctionLibrary("chooseGeneticCode", {"0" : "Universal"});

/* get filenames and read data */
fprintf(stdout, "MSA file:");
fscanf(stdin, "String", _msaFile);

fprintf(stdout, "Tree file:");
fscanf(stdin, "String", _treeFile);

fprintf(stdout, "Result file:");
fscanf(stdin, "String", _ancestorsFile);

DataSet msa = ReadDataFile(_msaFile);
DataSetFilter filteredData = CreateFilter(msa, 3, "", "", GeneticCodeExclusions);
/* DataSetFilter filteredData = CreateFilter(msa, 1, "", "", ""); */

/* set up model */
SelectTemplateModel(filteredData);

/* read tree */
ACCEPT_ROOTED_TREES = 1;
ACCEPT_BRANCH_LENGTHS = 1;

fscanf(_treeFile, "Tree", tree);

/* set up likelihood */
LikelihoodFunction lf = (filteredData, tree);

/* ancestor reconstruction */
DataSet ancestors = ReconstructAncestors(lf);
DataSetFilter ancestral_filter = CreateFilter(ancestors, 1);
DATA_FILE_PRINT_FORMAT = 9;
fprintf(_ancestorsFile, CLEAR_FILE, ancestral_filter);
