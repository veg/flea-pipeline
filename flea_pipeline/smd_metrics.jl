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

    incomplete = ~completecases(fulldf)
    if any(incomplete)
        missing_seqnames = fulldf[incomplete, :name]
        error("different sequences in fasta and copynumber files: $(missing_seqnames)")
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

function run_all(true_fasta_file, inf_fasta_file, inf_cn_file)
    true_df = read_files(true_fasta_file)
    inf_df = read_files(inf_fasta_file, inf_cn_file)
    results = []
    for timepoint in unique(true_df[:timepoint])
        true_df_sub = true_df[true_df[:timepoint] .== timepoint, :]
        true_seqs = true_df_sub[:sequence]
        true_cn = convert(Array{Float64}, true_df_sub[:copynumber])

        inf_df_sub = inf_df[inf_df[:timepoint] .== timepoint, :]
        inf_seqs = inf_df_sub[:sequence]
        inf_cn = convert(Array{Float64}, inf_df_sub[:copynumber])

        distmatrix, smd1, smd2, smd3 = run_metric(true_seqs, true_cn, inf_seqs, inf_cn)
        push!(results, [timepoint, smd1, smd2, smd3])
    end
    return results
end

function main()
    true_fasta_file = ARGS[1]
    inf_fasta_file = ARGS[2]
    inf_cn_file = ARGS[3]
    outfile = ARGS[4]

    results = run_all(true_fasta_file, inf_fasta_file, inf_cn_file)
    final_df = DataFrame()
    final_df[:timepoint] = [r[1] for r in results]
    final_df[:smd] = [r[2] for r in results]
    final_df[:false_positive] = [r[3] for r in results]
    final_df[:false_negative] = [r[4] for r in results]

    writetable(outfile, final_df)
end

main()