#!/usr/bin/env ruby


######################################################
######## INIT
######################################################

require "optparse" 
require 'fileutils'


CTRL = 0
TREAT = 1

######################################################
######## METHODS
######################################################

def load_list(input)
	if File.file?(input)
		list = File.readlines(input).map {|line| line.chomp!}
	else
		list = input.split(",")
	end
	return list
end



def load_table(input_file, blacklist = nil, whitelist = nil, filter = nil)
	header = nil 
	table = {}
	blacklist = [] if blacklist.nil?
	File.open(input_file).each do |line|
		line = line.chomp.split("\t")
		if header.nil?
			line.shift
			header = line
			next
		end
		next if blacklist.include?(line[0])
		if !whitelist.nil?
			next unless whitelist.include?(line[0])
		end 
		unless filter.nil?
			next unless filter[1].include?(line[header.index(filter[0]) + 1]) 
		end
		table[line[0].to_sym] = {}
		header.each_with_index do |feature_name, feat_index|
			table[line[0].to_sym][feature_name.to_sym] = line[feat_index + 1]
		end
	end
	return table
end


#target_name>
# 		feature_name:
#			feature_ctrl, feature_tr1/feature_tr2
		;
#		feature_nameN:
#			feature_ctrl, featue_trn 
def parse_target(target_string)
	parsed_target = {}
	target_string = target_string.split(">")
	target_name = target_string.shift
	all_features = target_string.shift.split(";")
	all_features.each do |factor|
		
		feature_name, features = factor.split(":")
		features = features.split(",")
		features = features.map{|feature| feature = feature.split("/")}
		parsed_target[feature_name.to_sym] = features
	end
	return target_name, parsed_target
end


def build_target(table, target)
		new_target = {"Ctrl" => [], "Treat" => []}
		ctrl_features = filter_features(target, CTRL)
		new_target["Ctrl"] = find_features(ctrl_features, table)
		treat_features = filter_features(target, TREAT)
		new_target["Treat"] = find_features(treat_features, table)
	return new_target

end

def filter_features(features, factor)
	filtered_ft = features.each_with_object({}) do |(ft_name, feature), result| 
		result[ft_name]= feature[factor] 
	end
	return filtered_ft
end

def find_features(features, table)
	samples_list = []
	table.each do |sample, all_features|
		include_sample = true
		features.each do |feature_name, ft_values|
			if !ft_values.include?(all_features[feature_name.to_sym])  
				include_sample = false
				break
			end
		end
		samples_list << sample if include_sample
	end
	return samples_list
end

def save_target(target_name, treats, output_path, experiment_design, additional_columns)
	File.open("#{output_path}/#{target_name.to_s}_target.txt",'w') do |out_file|
		header = "sample\ttreat"
		header = "sample\ttreat\t#{additional_columns.join("\t")}" if additional_columns.length > 0

		out_file.puts header

		treats.each do |treat, samples|
			samples.each do |sample|
				if additional_columns.length > 0
					features_by_sample = []
					additional_columns.each do |additional_feature|
						features_by_sample << experiment_design[sample][additional_feature]
					end
					out_file.puts "#{sample}\t#{treat}\t#{features_by_sample.join("\t")}" if !sample.nil?
				else
					out_file.puts "#{sample}\t#{treat}" if !sample.nil?
				end

			end
		end
	end
end

def save_aux_options(target_name, output_path, aux_options)
	FileUtils.mkdir_p output_path
	aux_file = "#{output_path}/#{target_name}_target.aux"
	File.open(aux_file, 'w') do |out|
		out.puts aux_options
	end
end


def parse_filter(string)
	feature_name, features = string.split("=")
	filter = [feature_name, features.split(",")]
	return filter
end
######################################################
######## OPTONS
######################################################
options = {}

OptionParser.new do |opts|

	options[:table] = nil
	opts.on("-e FILE", "--exp_file FILE", "Tabulated file which describes current experiment") do |file|
		options[:table] = file
	end

	options[:filter] = nil
	opts.on("-f STRING", "--filter STRING", "Set filter as string 'FEATURE_NAME=feature'") do |string|
		options[:filter] = string
	end

	options[:target] = nil
	opts.on("-t STRING", "--target STRING", "String which describes target. EXAMPLE: 'TARGET_A>COLUMN_A:FEAT_CTL1,FEAT_TRT1;TARGET_B>COLUMN_B:FEAT_CTL1/FEAT_CTL2,FEAT_TRT1/FEAT_TRT2'") do |string|
		options[:target] = string
	end

	options[:additional_features] = []
	opts.on("--additional_features STRING", "String with extra factors separated by commas to be added to target.") do |string|
		options[:additional_features] = string.split(",")
	end

	options[:aux_options] = []
	opts.on("--aux_options STRING", "String with extra options for DEGenesHunter.") do |string|
		options[:aux_options] = string
	end

	options[:blacklist] = nil
	opts.on("-b FILE/STRING", "--blacklist FILE/STRING", "List with samples name to exclude from target. File or comma separated string") do |file|
		options[:blacklist] = file
	end

	options[:whitelist] = nil
	opts.on("-w FILE/STRING", "--whitelist FILE/STRING", "List with samples name to acept from target. File or comma separated string") do |file|
		options[:whitelist] = file
	end

	options[:output_path] = "."
	opts.on("-o PATH", "--output_path PATH", "Set the output path") do |path|
		options[:output_path] = path
	end

end.parse!

######################################################
######## MAIN
######################################################
blacklist = load_list(options[:blacklist]) if !options[:blacklist].nil?
whitelist = load_list(options[:whitelist]) if !options[:whitelist].nil?

filter = parse_filter(options[:filter]) if !options[:filter].nil?
experiment_design = load_table(options[:table], blacklist, whitelist, filter)
target_name, target = parse_target(options[:target])
target = build_target(experiment_design, target) # meter aqui lo de guardar columnas y quitar la columna replicate

save_target(target_name, target, options[:output_path], experiment_design, options[:additional_features].map!{|feature| feature.to_sym})
save_aux_options(target_name, options[:output_path], options[:aux_options]) if !options[:aux_options].empty?