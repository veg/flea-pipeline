LoadFunctionLibrary ("HXB2Mapper", {"0": "Universal"});
LoadFunctionLibrary ("DescriptiveStatistics");
LoadFunctionLibrary ("CodonTools");
LoadFunctionLibrary ("ReadDelimitedFiles");

fprintf(stdout,"\nPlease enter the path of the final nucleotide alignment:");
fscanf  (stdin, "String", _nucSequences);

fprintf(stdout,"\nOutput file:");
fscanf  (stdin, "String", _outputFile);

DataSet         allData         = ReadDataFile (_nucSequences);
DataSetFilter   codonData       = CreateFilter (allData,3,"","",GeneticCodeExclusions);
DataSetFilter   nucData         = CreateFilter (allData,1);

GetDataInfo (codonMap, codonData);

call_count = 0;
GetString          (inputSequenceOrder, allData, -1);
COUNT_GAPS_IN_FREQUENCIES     = 0;
GetDataInfo (cons, nucData, "CONSENSUS");

cons = translateCodonToAA (cons, defineCodonToAA(), 0);
_hxb_alignOptions_prot ["SEQ_ALIGN_NO_TP"] = 0;
hxb2coord = mapSequenceToHXB2Aux (cons, _HXB2_AA_ENV_, 1);

dict={};
for (k = 0; k < Columns (_HXB_aa_offset_matrix); k += 1) {
  for(h = 0; h < Rows(hxb2coord); h=h+1){
    if(hxb2coord[h][0] == _HXB_aa_offset_matrix[0][k])
      {
        dict[""+k]=h;
        break;
      }
  }
}
dict[""+(Columns(_HXB_aa_offset_matrix)-1)]=Rows(hxb2coord);
fprintf(_outputFile,CLEAR_FILE,dict,"\n\n");

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
