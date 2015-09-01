using ArgParse
using Dates
using JSON

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "Filename"
            help = "The JSON file"
            required = true
            arg_type = String
        "OutFile"
            help = "The JSON output file"
            required = true
            arg_type = String
        "Regularization"
            help = "A number between 0 and 1 - default: 0.01"
            required = true
            arg_type = Real
    end
    return parse_args(s)
end

parsed_args = parse_commandline()

reg=parsed_args["Regularization"];
filename=string(parsed_args["Filename"]);
fileout=string(parsed_args["OutFile"]);

AAs="FLIMVSPTAYXHQNKDECWRG-";

function getAAVec(dateDict)
         return map(aa->get(dateDict,string(aa),0), collect(AAs))
end

function norm(vec)
         return float(vec/sum(vec))
end

function getReadCount(dict,key)
        count=0;
        for k in sortedDates
                count=count+sum(getAAVec(dict[key][k]));
        end
        return count
end

function entropy(vec)
        gZero=filter(x->x>0.0,vec)
        return -sum(gZero.*log(2, gZero))
end

function rSqr(vec1,vec2,psFreq)
        freq1=norm(vec1+psFreq);
        freq2=norm(vec2+psFreq);
        avg=0.5*freq1+0.5*freq2;
        ent=entropy(avg)+entropy([0.5,0.5])-(entropy(freq1/2)+entropy(freq2/2));
        return(ent)
end

f=open(filename);
jdict=JSON.parse(f,ordered=true);
dateStrings=filter(x->x!="HXB2",collect(keys(jdict["1"])));
sortedDates=sort(dateStrings,by=y->DateTime(y,"yyyymmdd"));

keys(jdict)

n_positions = maximum(map(i->int(i), keys(jdict)))
if n_positions != length(keys(jdict))
    throw("error")
end

newDict = Dict()
for i in keys(jdict)
    tempDict=Dict();
    tempDict["readCount"]=getReadCount(jdict,i);
    for k in 2:length(sortedDates)
        tempDict[sortedDates[k]]=round(rSqr(getAAVec(jdict[i][sortedDates[k-1]]),
                                            getAAVec(jdict[i][sortedDates[k]]),
                                            reg), 4)
    end
    newDict[string(i)]=tempDict
end

finalDict = Dict()
for k in 2:length(sortedDates)
    date = sortedDates[k]
    newarr = Array(Float32, n_positions)
    for i in keys(jdict)
        newarr[int(i)] = newDict[string(i)][date]
     end
    finalDict[date] = newarr
end

fout=open(fileout,"w");
JSON.print(fout, finalDict)
