require_relative './InteractionNetwork'
require_relative './goterms'
require_relative './keggterms'

# Checking if the arguments required are specified 

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


### ---------------- READ GENES FROM THE FILE ------------------- ###

gene_list = Array.new       # Creating an empty array for saving each gene of the document

# Source: https://www.rubyguides.com/2015/05/working-with-files-ruby/

File.foreach(ARGV[0].to_s){ |line|
    gene = line.gsub("\n",'')               # We eliminate metacharacter \n
    unless gene.match(/AT\dG\d{5}/i)        # Check if genes belong to Arabidopsis and save each gene in the array.
        abort("ERROR: the gene list have some errors. #{gene} has not correct format") 
    end
    gene_list <<  gene.upcase!}


### ---------- Find interactions and generate Interaction Network report ----------- ###

puts "\nFinding Interaction Networks..."

InteractionNetwork.find_interactions(gene_list: gene_list)      # Find full interaction network of our input genes.

output_interactions = Array.new     # Initialize Array to store the interactions between genes in the list

gene_list.each {|gene|

    network = InteractionNetwork.new(gene: gene, max_depth: 3)      # Generate interaction networks for each query gene

    if network.interactors.length > 1

        inlist_interactors = network.interactors_within(gene_list: gene_list)   # Find interactors within the input gene list

        if inlist_interactors.length > 1
    
            inlist_interactors.sort!                                    # To avoid the presence of interactions with the same genes but different order
            unless output_interactions.include?(inlist_interactors)     # Check of the within-list interaction has already been saved

                output_interactions << inlist_interactors       # Save interactions between genes of the list

            end
        end
    end
}


### --------- Write Output Interaction Report ---------- ###

puts "\nGenerating the report..."

outfile = File.new(ARGV[1], "w+")   # Create output file

outfile.write("#####################################\n")
outfile.write("Final Report:\n")
outfile.write("A total number of #{InteractionNetwork.get_all_networks.length.to_s} protein-protein Interaction Networks have been identified where input genes are involved, with a min. quality IntactMiscore > 0.3 and max. depth = 3.")
outfile.write("#{output_interactions.length.to_s} interactions have been identified between different some of the #{gene_list.length.to_s} genes of the list:")

output_interactions.each { |inlist_interactors|

    interactors_out = inlist_interactors.map { |x| x.to_s }

    outfile.write("\n\n/-----------------------------------------------/")
    outfile.write("\n\n\t#{interactors_out.join("-")}\n\n")
    
    go_terms = GoTerms.new(inlist_interactors)
    kegg_terms = KeggTerms.new(inlist_interactors)
    
    unless go_terms.go_terms.empty?     # Write GO IDs and terms if exist
        outfile.write("Gene Ontology:\n")
        go_terms.go_terms.uniq.each {|go| outfile.write(" - #{go.keys[0]} => #{go.values[0].split(":")[1]}\n")}
    end
    unless kegg_terms.pathways.empty?   # Write KEGG IDs and terms if exist
        outfile.write("\nKEGG Pathways:\n")
        kegg_terms.pathways.uniq.each {|kegg| outfile.write(" - #{kegg.keys[0]} => #{kegg.values[0]}\n")}
    end
}        

outfile.close()    # Close output file

puts "\nInteraction Report is completed!"
   
