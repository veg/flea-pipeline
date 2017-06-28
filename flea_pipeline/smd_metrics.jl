using DataFrames
using Glob
using BeNGS

function read_files(fasta_file, copynumber_file="")
    templates = BeNGS.read_fasta_records(fasta_file)

    template_df = DataFrame()
    template_df[:timepoint] = [Symbol(split(s.name, "_")[1]) for s in templates]
    template_df[:name] = [s.name for s in templates]
    template_df[:sequence] = [String(s.seq) for s in templates]

    fulldf = if length(copynumber_file) > 0
        copynumber_df = readtable(copynumber_file, header=false,
                                  names=[:name, :copynumber])
        join(template_df, copynumber_df, on=:name, kind=:outer)
    else
        template_df[:copynumber] = fill(1, size(template_df[1]))
        template_df
    end

    if !all(completecases(fulldf))
        error("sequences do not match copynumbers")
    end

    nrows, _ = size(fulldf)
    if length(unique(fulldf[:name])) < nrows
        error("non-unique sequence names")
    end
    return fulldf
end

function run_metric(true_seqs, true_cn, inf_seqs, inf_cn)
    # run for inferred sequences
    distmat = BeNGS.default_dist_matrix(true_seqs, inf_seqs)
    emd1, _ = BeNGS.smd_mf(distmat; freq1=true_cn, freq2=inf_cn)
    emd2, _ = BeNGS.smd_mf(distmat; freq1=true_cn, freq2=inf_cn,
                           unbounded_first=true)
    emd3, _ = BeNGS.smd_mf(distmat; freq1=true_cn, freq2=inf_cn,
                           unbounded_second=true)
    return distmat, emd1, emd2, emd3
end

function run_all(timepoint, true_fasta_file, inf_fasta_file, inf_cn_file)
    true_df = read_files(true_fasta_file)
    inf_df = read_files(inf_fasta_file, inf_cn_file)

    true_df_sub = true_df[true_df[:timepoint] .== timepoint, :]
    true_seqs = true_df_sub[:sequence]
    true_cn = convert(Array{Float64}, true_df_sub[:copynumber])
    
    inf_df_sub = inf_df[inf_df[:timepoint] .== timepoint, :]
    inf_seqs = inf_df_sub[:sequence]
    inf_cn = convert(Array{Float64}, inf_df_sub[:copynumber])

    distmatrix, smd1, smd2, smd3 = run_metric(true_seqs, true_cn, inf_seqs, inf_cn)
    return [timepoint, smd1, smd2, smd3]
end

function main()
    timepoint = Symbol(ARGS[1])
    true_fasta_file = ARGS[2]
    inf_fasta_file = ARGS[3]
    inf_cn_file = ARGS[4]
    outfile = ARGS[5]

    result = run_all(timepoint, true_fasta_file, inf_fasta_file, inf_cn_file)
    final_df = DataFrame()
    final_df[:timepoint] = [result[1]]
    final_df[:smd] = [result[2]]
    final_df[:false_positive] = [result[3]]
    final_df[:false_negative] = [result[4]]

    writetable(outfile, final_df)
end

main()