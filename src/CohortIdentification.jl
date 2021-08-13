module CohortIdentification

import Pkg
using ArgParse
using CSV
using ArgParse
using DataFrames
using DataFramesMeta


function metadata_variable_type(metadata::DataFrame)::DataType

	map(typeof, metadata[!,2]) |> Set |> only

end


function count_elements(elements::Array{T, 1})::Dict{T, Int64} where T <: Union{String, Int64}

	Dict(x => count(elements .== x) for x in Set(elements))

end


function count_metadata(
	cluster_name::Int64,
	clusters::Array{Int64, 1}, isolates::Array{String, 1},
	metadata::DataFrame
	)

	mask::BitArray = clusters .== cluster_name

	selected_isolates = Set(isolates[mask])

	metadata_isolates = @select(metadata, :isolate)[!, 1]
	selected_metadata::Array{Bool, 1} = [
		x in selected_isolates
		for x in metadata_isolates
		]
	metadata_counts = begin
		metadata[selected_metadata, 2] |>
		count_elements
	end

	metadata_counts
end


function percent_positive(positive_value, metadata_counts::Dict)::Float64

	total::Int64 = sum(values(metadata_counts))

	positive::Int64 = try
		metadata_counts[positive_value]
	catch KeyError
		0
	end

	positive / total

end


# need better name & return type
function count_positives_for_all_clusters(
		clusters::DataFrame,
		metadata::DataFrame
		)::Dict{String, Dict{Int64, Dict{Any, Int64}}}


	drop_isolate_names = select(clusters, Not(:isolate))
	threshold_clusters = eachcol(drop_isolate_names)
	threshold_names = names(drop_isolate_names)
	isolates_array::Array{String, 1} = @select(clusters, :isolate)[!, 1]

	positive_counts = Dict(n => Dict() for n in threshold_names)

	# for each threshold
	for (n::String, c::Array{Int64,1}) in zip(threshold_names, threshold_clusters)


		# for each cluster in that threshold
		for cluster in Set(c)

			# count the number of strains that fall into
			# each value for the metadata variable
			metadata_counts = count_metadata(cluster,
											 c,
											 isolates_array,
											 metadata)

			positive_counts[n][cluster]  = metadata_counts
		end

	end
	positive_counts
end


function get_cluster_sizes(clusters::DataFrame)::Dict{String, Dict{Any, Int64}}

	drop_isolate_names = select(clusters, Not(:isolate))
	threshold_clusters = eachcol(drop_isolate_names)
	threshold_names = names(drop_isolate_names)


	Dict(n => count_elements(c)
		 for (n, c)
		 in zip(threshold_names, threshold_clusters))
end


# function to get positivity rates for each cluster
function get_cluster_positivity_rates(positive_value, results)

	out = Dict(threshold => Dict() for  threshold in keys(results))

	for threshold in keys(results)

		for cluster in keys(results[threshold])
			counts = results[threshold][cluster]

			out[threshold][cluster] = percent_positive(positive_value, counts)

		end
	end
	out
end


function format_cohort(
		cohort::Dict{String, Set{Any}},
		group::String,
		sizes::Dict{String, Dict{Any, Int64}},
		positive_counts::Dict{String, Dict{Int64, Int64}}
		)::Array{NamedTuple, 1}

	rows = Array{NamedTuple, 1}()

	thresholds = keys(cohort)

	for threshold_name in thresholds

		for cluster_name_ in cohort[threshold_name]

			row = (threshold = threshold_name,
				   cluster = cluster_name_,
				   cohort = group,
				   size = sizes[threshold_name][cluster_name_],
				   positive = positive_counts[threshold_name][cluster_name_]
			)

			push!(rows, row)

		end

	end

	rows

end


function format_output(
		positive_clusters::Dict{String,Set{Any}},
		negative_clusters::Dict{String,Set{Any}},
		sizes::Dict{String, Dict{Any, Int64}},
		counts::Dict{String, Dict{Int64, Dict{Any, Int64}}},
		positive_value,
		)::DataFrame


	positive_counts = Dict{String, Dict{Int64, Int64}}()

	for (threshold, clusters_) in counts

		positive_counts[threshold] = Dict()

		for cluster in keys(clusters_)
			pos = try
				counts[threshold][cluster][positive_value]
			catch KeyError
				0
			end
			positive_counts[threshold][cluster] = pos

		end
	end

	vcat(format_cohort(positive_clusters, "positive", sizes, positive_counts),
		 format_cohort(negative_clusters, "negative", sizes, positive_counts)) |>
	DataFrame

end


function write_output(out_table::DataFrame, out_path::String, delimiter)

	sort!(out_table, [:cohort, :threshold, :cluster])
	CSV.write(out_path, out_table; delim=delimiter)
end



function pos_adj(cluster_size::Int64)::Float64

	Δsize = cluster_size - size_1

	pos_1 - (Δsize * posΔ_sizeΔ_ratio)

end


function neg_adj(cluster_size::Int64)::Float64

	Δsize = cluster_size - size_1

	neg_1 + (Δsize * negΔ_sizeΔ_ratio)

end


function is_positive_cohort(
	cluster_size::Int64,
	cluster_positivity::Float64
	)::Bool

	# This isn't required, buit I want it clear
	# that small clusters are not considered
	if cluster_size  < size_1
		false

	elseif cluster_size == size_1 && cluster_positivity >= pos_1
		true

	elseif cluster_size >= size_2 && cluster_positivity >= pos_2
		true

	elseif size_1 < cluster_size < size_2 &&
		   (cluster_positivity >= pos_adj(cluster_size))
		true

	else
		false

	end
end


function is_negative_cohort(
	cluster_size::Int64,
	cluster_positivity::Float64
	)::Bool

	# This isn't required, buit I want it clear
	# that small clusters are not considered
	if cluster_size < size_1
		false

	elseif (cluster_size == size_1) && (cluster_positivity <= neg_1)
		true

	elseif (cluster_size >= size_2) && (cluster_positivity <= neg_2)
		true

		# NB I flipped the >= in Ed's description to <= here
	elseif (size_1 < cluster_size < size_2) &&
		(cluster_positivity <= neg_adj(cluster_size))
		true

	else
		false

	end
end


function populate_cohorts(sizes, cluster_positivities)

		thresholds = keys(sizes)

		positive_clusters = Dict(threshold => Set() for threshold in thresholds)
		negative_clusters = Dict(threshold => Set() for threshold in thresholds)

		for threshold in thresholds

			cluster_names = keys(sizes[threshold])

			for cluster_name in cluster_names

				cluster_size = sizes[threshold][cluster_name]
				cluster_positivity = cluster_positivities[threshold][cluster_name]

				if is_positive_cohort(cluster_size, cluster_positivity)
					push!(positive_clusters[threshold], cluster_name)
				end

				if is_negative_cohort(cluster_size, cluster_positivity)
					push!(negative_clusters[threshold], cluster_name)
				end

			end
		end

	positive_clusters, negative_clusters
end

function arguments(cli)

	s = ArgParseSettings()

	s.version = "0.2.2"
	s.add_version = true

	@add_arg_table! s begin

		"--size-1"
			help = "The lower bound size for consideration"
			arg_type = Int64
			required = true

		"--size-2"
			arg_type = Int64
			required = true

		"--pos-1"
			arg_type = Float64
			required = true

		"--pos-2"
			arg_type = Float64
			required = true

		"--neg-1"
			arg_type = Float64
			required = true

		"--neg-2"
			arg_type = Float64
			required = true

		"--clusters"
			arg_type = String
			required = true

		"--metadata"
			arg_type = String
			required = true

		"--variable"
			arg_type = String
			required = true

		"--positive-value"
			arg_type = String
			required = true

		"--output", "-o"
			arg_type = String
			required = true

		"--delimiter", "-d"
			arg_type = String
			default = "\t"
	end

	args = parse_args(cli, s)

	global size_1 = args["size-1"]
	global size_2 = args["size-2"]
	global pos_1 = args["pos-1"]
	global pos_2 = args["pos-2"]
	global neg_1 = args["neg-1"]
	global neg_2 = args["neg-2"]

	global posΔ = pos_1 - pos_2
	global negΔ = neg_2 - neg_1
	global sizeΔ = size_2 - size_1

	global posΔ_sizeΔ_ratio = (posΔ / sizeΔ)
	global negΔ_sizeΔ_ratio = (negΔ / sizeΔ)

	args
end

function main(cli)

	args = arguments(cli)

	# Load the cluster data
	# Strain names are in leftmost column, "isolate"

	clusters = begin
		cluster_path = args["clusters"]
		CSV.File(cluster_path, delim=args["delimiter"]) |>  DataFrame
	end


	# Load the metadata
	# return a dataframe with two variables: the isolate name,
	# and the cohort definition variable

	metadata = begin
		column_of_interest = args["variable"]
		metadata_path = args["metadata"]
		_metadata = CSV.File(metadata_path, delim=args["delimiter"]) |> DataFrame
		select(_metadata, "isolate", column_of_interest)
	end

	sizes = get_cluster_sizes(clusters)

	counts = count_positives_for_all_clusters(clusters, metadata)
	pos_val_type = metadata_variable_type(metadata)

	positive_value = parse(pos_val_type, args["positive-value"])

	cluster_positivities = get_cluster_positivity_rates(positive_value, counts)

	positive_clusters, negative_clusters = populate_cohorts(sizes, cluster_positivities)

	output_df = format_output(
		positive_clusters,
		negative_clusters,
		sizes,
		counts,
		positive_value
	)

	write_output(output_df, args["output"])
end

end
