using ArgParse

include("lib.jl")

function arguments()

	s = ArgParseSettings()

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
	end

	return parse_args(s)
end


args = arguments()

const size_1 = args["size-1"]
const size_2 = args["size-2"]
const pos_1 = args["pos-1"]
const pos_2 = args["pos-2"]
const neg_1 = args["neg-1"]
const neg_2 = args["neg-2"]

const posΔ = pos_1 - pos_2
const negΔ = neg_2 - neg_1
const sizeΔ = size_2 - size_1

const posΔ_sizeΔ_ratio = (posΔ / sizeΔ)
const negΔ_sizeΔ_ratio = (negΔ / sizeΔ)


# Load the cluster data
# Strain names are in leftmost column, "isolate"
const clusters = begin
	cluster_path = args["clusters"]
	CSV.read(cluster_path, DataFrame)
end


# Load the metadata
# return a dataframe with two variables: the isolate name,
# and the cohort definition variable

const metadata = begin
	column_of_interest = args["variable"]
	metadata_path = args["metadata"]
	_metadata = CSV.read(metadata_path, DataFrame)
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
