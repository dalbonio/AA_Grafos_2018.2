using DelimitedFiles

if lenght(ARGS) != 1
	println("no file passed by parameter")
	exit
end

input_file = open(ARGS[1])
dimacs_file = open(string(ARGS[1], ".dimacs"))

write(dimacs_file, string("c File:" ARGS[1], ".dimacs"))
write(dimacs_file, string("p edge ", readline(input_file)))

adjacency_list = readdlm(input_file)

new_dict = Dict{String, String}()

for (i, neibh) in enumerate(adjacency_list)