/*----------------------------------------------------------------------------------------------------------*/

function tools.matrix_to_JSON (mx) {
    json_out = ""; json_out * 128;
    json_out * "[";
    for (k = 0; k < Rows (mx); k+=1) {
        if (k) {
            json_out * ",";
        }
        json_out * ("[" + Join (",", mx[k][-1]) + "]");
    }
    json_out * "]";
    json_out * 0;
    return json_out; 
}

function tools.addKeyIfMissing (array, key, value) {
	if (Abs (array[key]) == 0 && Type (array[key]) == "Number") {
		array[key] = value;
	}
}