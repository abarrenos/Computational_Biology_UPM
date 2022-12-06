
require_relative './Interactions_Adrian'

# Checking if the arguments required are specified 
'''
unless ARGV.length == 2
    abort("FATAL ERROR: Files pathways are required. \nHELP MESSAGE: Check README.md for more information.")
end

#Checking if the input file exists

unless File.file?(ARGV[0])
    abort("FATAL ERROR: File #{ARGV[0]} does not exist or the pathway specified is not correct")
end


#Checking if the output file already exists and asking if it should be overwrite

if File.file?(ARGV[1])
    puts "#{ARGV[1]} already exists, indicate if you want to overwrite [Y/N]" 
    stdin = ""
    until stdin == "n" || stdin == "N" || stdin == "y" || stdin == "Y"
        stdin = STDIN.gets.strip
        if stdin == "N" || stdin == "n"
            abort("Run cancelled")
        end
    end
end
'''
### ---------------- READ GENES FROM THE FILE ------------------- ###

gene_list = Array.new       # Creating an empty array for saving each gene of the document

# Source: https://www.rubyguides.com/2015/05/working-with-files-ruby/

File.foreach("./documents/ArabidopsisSubNetwork_GeneList.txt"){ |line|
    gene = line.gsub("\n",'')        # We eliminate metacharacter \n
    unless gene.match(/AT\dG\d{5}/i) # Check if genes belong to Arabidopsis and save each gene in the array.
        abort("ERROR: the gene list have some errors. #{gene} has not correct format") 
    end
    gene_list <<  gene.upcase!}

########### --------------- Main Cycle ---------------- ##########

InteractionNetwork.find_interactions(gene_list: gene_list)
    
gene_list.each {|gene|

    network = InteractionNetwork.new(gene: gene, max_depth: 25)
    if network.interactors_within(gene_list: gene_list).length > 2
        print network.query_gene
        print "\t", network.interactors.length
        print "\t", network.depth
        '''
        print "\n", network.interactors
        print "\t", @@interaction_dict[gene.to_sym].length unless @@interaction_dict[gene.to_sym].nil?
        '''
        puts 

        print "\n", network.interactors_within(gene_list: gene_list)
        puts
        puts
    end
}
#InteractionNetwork.find_interactions(gene_list[0..-1])

