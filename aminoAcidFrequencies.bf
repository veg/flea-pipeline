LoadFunctionLibrary ("HXB2Mapper", {"0": "Universal"});
LoadFunctionLibrary ("lib/dates.bf");
LoadFunctionLibrary ("DescriptiveStatistics");
LoadFunctionLibrary ("CodonTools");
LoadFunctionLibrary ("ReadDelimitedFiles");
LoadFunctionLibrary ("lib/dates.bf");

_use_copy_numbers = "_([0-9]+)$";

fprintf(stdout,"\nPlease enter the base path:");
fscanf  (stdin, "String", basePath);

_nucSequences = "`basePath`/input/merged.fas";
_dates        = "`basePath`/input/merged.dates";
_mrca_from    = "`basePath`/input/mrca.seq";
_frequencies_to = "`basePath`/results/frequencies.json";



DataSet         allData         = ReadDataFile (_nucSequences);
DataSetFilter   codonData       = CreateFilter (allData,3,"","",GeneticCodeExclusions);
DataSetFilter   nucData         = CreateFilter (allData,1);
GetDataInfo (codonMap, codonData);

dateInfo     = _loadDateInfo (_dates);
uniqueDates  = _getAvailableDates (dateInfo);

call_count = 0;
GetString          (inputSequenceOrder, allData, -1);
COUNT_GAPS_IN_FREQUENCIES     = 0;

GetDataInfo (cons, nucData, "CONSENSUS");
cons = translateCodonToAA (cons, defineCodonToAA(), 0);
_hxb_alignOptions_prot ["SEQ_ALIGN_NO_TP"] = 0;
hxb2coord = mapSequenceToHXB2Aux (cons, _HXB2_AA_ENV_, 1);

senseMap = defineSenseCodonToAA();

bySite = {};

sequence_weights = ({codonData.species,1})["1"];


if (Abs(_use_copy_numbers)) {
    for (s = 0; s < codonData.species; s += 1) {
        GetString (sName, codonData, s);
        m = sName $ _use_copy_numbers;
        if (m[0] > 0) {
            sequence_weights [s] =  0 + sName[m[2]][m[3]];
        }
    }
}

//fprintf (stdout, sequence_weights, "\n");

for (date_index = 0; date_index < Abs (uniqueDates); date_index += 1) {
    thisDate = uniqueDates["INDEXORDER"][date_index];
    
    sequences = _selectSequencesByDate (thisDate, dateInfo, inputSequenceOrder);

    for (s = 0; s < codonData.sites; s += 1) {
        if (date_index == 0) {
            bySite [s+1] = {"HXB2": hxb2coord__[s__]+1};
        }
        aa_counts = {};
        for (q = 0; q < Columns (sequences); q+=1) {
        
            GetDataInfo (char, codonData, sequences[q], codonMap[s]);
            total = +char;
            if (total > 0 && total < Rows (char)) {
                char = char["_MATRIX_ELEMENT_VALUE_*(1+_MATRIX_ELEMENT_ROW_)"];
                char = char[char["_MATRIX_ELEMENT_VALUE_"]];
                for (c = 0; c < Columns (char); c+=1) {
                   aa_counts [senseMap[char[c]-1]] += sequence_weights[sequences[q]]/total;
                }
            }
        }
        (bySite[s+1])[thisDate] = aa_counts;
    }
}

fprintf (_frequencies_to, CLEAR_FILE, bySite);

/*----------------------------------------------------------------------------------------------------------*/

function defineSenseCodonToAA ()
{	
	codonToAAMap = {};
	
	p2 = 0;
	for (p1=0; p1<64; p1 += 1)
	{
	    if (_Genetic_Code[p1] == 10) {
	        continue;
	    }
		codonToAAMap[p2] = _hyphyAAOrdering[_Genetic_Code[p1]];
		p2 += 1;
	}
	
	
	
	return codonToAAMap;
}