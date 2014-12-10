LoadFunctionLibrary ("GrabBag");

function _loadDateInfo (file) {
    fscanf (file, "Raw", info);
    return Eval (info);
}

function _getAvailableDates (info) {
    _toReturn = {};
    info["_accessor1"][""];
    return _toReturn;
}

function _selectSequencesByDate (date, dates, inputOrder) {
    set = {};
    dates ["_accessor2"][""];
    return mapSets (Rows (set), inputOrder);
}


//-----------------------------------------------------------------

function _accessor1 (key, value) {
    if (_toReturn[value] == 0) {
        _toReturn[value] = 1; 
    }
    return 0;
}

function _accessor2 (key, value) {
    if (value == date) {
        set [key] = 1;
    }
    return 0;
}

