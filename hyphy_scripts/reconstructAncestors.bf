LoadFunctionLibrary ("ReadDelimitedFiles");
LoadFunctionLibrary ("tools.bf");

/* UNCOMMENT FOR CODON MODEL */
/* LoadFunctionLibrary("chooseGeneticCode", {"0" : "Universal"}); */

fprintf(stdout, "\nMSA file:");
fscanf(stdin, "String", _msaFile);

fprintf(stdout, "\nTree file:");
fscanf(stdin, "String", _treeFile);

fprintf(stdout, "\nResult file:");
fscanf(stdin, "String", _ancestorsFile);

fprintf(stdout, "\nRunning\n");

DataSet msa = ReadDataFile(_msaFile);

/* UNCOMMENT FOR CODON MODEL */
/* DataSetFilter filteredData = CreateFilter(msa, 3, "", "", GeneticCodeExclusions); */

/* COMMENT FOR CODON MODEL */
DataSetFilter filteredData = CreateFilter(msa, 1, "", "", "");

SelectTemplateModel(filteredData);

ACCEPT_ROOTED_TREES = 1;
ACCEPT_BRANCH_LENGTHS = 1;
fscanf(_treeFile, "Tree", tree);

LikelihoodFunction lf = (filteredData, tree);

DataSet ancestors = ReconstructAncestors(lf);
DataSetFilter ancestral_filter = CreateFilter(ancestors, 1);
DATA_FILE_PRINT_FORMAT = 9;
fprintf(_ancestorsFile, CLEAR_FILE, ancestral_filter);
