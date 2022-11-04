require './gene.rb'

genelist = Array.new
f = File.open("./gene_list.txt", "r")
f.readlines.each do |line|
    line.strip!  # remove all spaces and newline characters
    g = Gene.new(id: line)
    puts "created Gene #{g.id}"

end

Gene.get_genes.each do |gene|
    puts "I have retrieved gene #{gene.id}"
end
