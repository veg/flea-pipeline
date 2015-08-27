LoadFunctionLibrary ("HXB2Mapper", {"0": "Universal"});
LoadFunctionLibrary ("lib/dates.bf");
LoadFunctionLibrary ("DescriptiveStatistics");
LoadFunctionLibrary ("CodonTools");
LoadFunctionLibrary ("NJ");
LoadFunctionLibrary ("ReadDelimitedFiles");
LoadFunctionLibrary ("lib/tools.bf");

/* important -- specify as absolute paths, otherwise FUBAR breaks */
fprintf(stdout,"\Enter the alignment file:");
fscanf  (stdin, "String", _nucSequences);

fprintf(stdout,"\Enter the dates file:");
fscanf  (stdin, "String", _dates);

fprintf(stdout,"\Enter the mrca file:");
fscanf  (stdin, "String", _mrca_from);

fprintf(stdout,"\Enter the working directory:");
fscanf  (stdin, "String", _fubar_directory);

fprintf(stdout,"\Enter the rates output file:");
fscanf  (stdin, "String", _rates_to);

_c2p_mapping = defineCodonToAA ();
COUNT_GAPS_IN_FREQUENCIES = 0;

sequences_by_date = {};


DataSet         allData         = ReadDataFile (_nucSequences);
DataSetFilter   codonData       = CreateFilter (allData,3,"","",GeneticCodeExclusions);
DataSetFilter   nucData         = CreateFilter (allData,1);

DataSet         mrca             = ReadDataFile  (_mrca_from);
DataSetFilter   mrcaFilter       = CreateFilter (allData,3,"","",mrca);

GetDataInfo     (mrcaSequence, mrcaFilter, 0);
sequences_by_date ["MRCA"] =   translateCodonToAA (mrcaSequence, _c2p_mapping, 0);

GetDataInfo (codonMap, codonData);

dateInfo     = _loadDateInfo (_dates);
uniqueDates  = _getAvailableDates (dateInfo);

call_count = 0;
GetString          (inputSequenceOrder, allData, -1);
COUNT_GAPS_IN_FREQUENCIES     = 0;

/*
GetDataInfo (cons, nucData, "CONSENSUS");
cons = translateCodonToAA (cons, defineCodonToAA(), 0);
hxb2coord = mapSequenceToHXB2Aux (cons, _HXB2_AA_ENV_, 1);
*/


byDate = {}; // the JSON output with rates



for (date_index = -1; date_index < Abs (uniqueDates); date_index += 1) {
    if (date_index < 0) {
        thisDate = "Combined";
        sequences = "";
    } else {
        thisDate                    = uniqueDates["INDEXORDER"][date_index];
        sequences                   = _selectSequencesByDate (thisDate, dateInfo, inputSequenceOrder);
    }

    DataSetFilter filteredData  = CreateFilter (allData,1,"",Join(",",sequences));

    if (date_index < 0) {
        aa_seqs      = {};

        for (s_id = 0; s_id < filteredData.species; s_id += 1) {
            GetString (seq_name, filteredData, s_id);
            GetDataInfo (codon_seq, filteredData, s_id);
            aa_seqs [seq_name] = translateCodonToAA (codon_seq, _c2p_mapping, 0);
        }

        buffer = ""; buffer * 128;
        aa_seqs ["addSeqToBuffer"][""];
        buffer * 0;

        DataSet aa_data = ReadFromString (buffer);
        aa_sequences = "";
        GetString          (inputSequenceOrder_aa, aa_data, -1);

    } else {
        aa_sequences    = _selectSequencesByDate (thisDate, dateInfo, inputSequenceOrder_aa);
    }

    DataSetFilter aa_filter  = CreateFilter (aa_data,1,"",Join(",",aa_sequences));
    if (date_index >= 0) {
         sequences_by_date [thisDate] = {};
         for (s = 0; s < aa_filter.species; s += 1) {
            GetString   (s_name, aa_filter, s);
            GetDataInfo (s_data, aa_filter, s);
            (sequences_by_date [thisDate])[s_name] = s_data;
         }
    }
    entropies = {aa_filter.sites, 1};
    for (s = 0; s < aa_filter.sites; s += 1) {
        DataSetFilter site_filter = CreateFilter (aa_filter, 1, "" + s);
        HarvestFrequencies (site_freqs, site_filter, 1, 1, 1);
        entropies[s] = entropy (site_freqs);
    }

    _date_slice                 = _fubar_directory + thisDate + ".nex";
    _date_FUBAR                 = _date_slice + ".fubar.csv";

    if (!_date_FUBAR) {

    } else {
        IS_TREE_PRESENT_IN_DATA     = 1;
        DATAFILE_TREE               = InferTreeTopology (1);
        fprintf                       (_date_slice, CLEAR_FILE, filteredData);
        ExecuteAFile                  (HYPHY_LIB_DIRECTORY + "TemplateBatchFiles" +
                                       DIRECTORY_SEPARATOR + "FUBAR.bf",
                                      {
                                       "00" : "Universal",
                                       "01" : "1",
                                       "02" : _date_slice,
                                       "03" : "12",
                                       "04" : "2",
                                       "05" : "510000",
                                       "06" : "255000",
                                       "07" : "100",
                                       "08" : "0.5"});
    }



    fubar_rates = (ReadCSVTable (_date_FUBAR,1))[1];
    n_sites   = Rows(fubar_rates);
    fubar_rates_n = {n_sites,5};
    mean_alpha = (+(fubar_rates[-1][1]))/n_sites;

    for (k = 0; k < n_sites; k += 1) {
        fubar_rates_n [k][0] = fubar_rates[k][1] / mean_alpha;
        fubar_rates_n [k][1] = fubar_rates[k][2] / mean_alpha;
        fubar_rates_n [k][2] = fubar_rates[k][4];
        fubar_rates_n [k][3] = fubar_rates[k][5];
        fubar_rates_n [k][4] = entropies  [k];
    }

    byDate [thisDate] = tools.matrix_to_JSON (fubar_rates_n);
}

//SetDialogPrompt ("Save to file (rates.json):");
fprintf (_rates_to, CLEAR_FILE, byDate);
//SetDialogPrompt ("Save to file (sequences.json):");
//fprintf (PROMPT_FOR_FILE, CLEAR_FILE, sequences_by_date);

/*----------------------------------------------------------------------------------------------------------*/

function addSeqToBuffer (key, value) {
    buffer * (">" + key + "\n" + value + "\n");
    return 0;
}


/*----------------------------------------------------------------------------------------------------------*/

function entropy (mx) {
    non_zero = mx[mx["_MATRIX_ELEMENT_VALUE_ > 0"]];
    return +(non_zero ["-_MATRIX_ELEMENT_VALUE_*Log(_MATRIX_ELEMENT_VALUE_)/Log(2)"]);
}
