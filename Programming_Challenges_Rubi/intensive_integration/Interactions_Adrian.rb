require 'rest-client'

class InteractionNetwork
    
    attr_accessor :query_gene
    attr_accessor :interactors
    attr_accessor :depth

    @@interaction_dict = {}
    @@network_depth = 0

    def initialize(gene:, interaction_dict: @@interaction_dict, max_depth: 3)
    
        @query_gene = gene.to_sym
        @interactors = InteractionNetwork.network(gene_id: gene, \
                                                from_dict: interaction_dict, \
                                                max_depth: max_depth)
        @depth = @@network_depth
        @@network_depth = 0
    
    end

    def interactors_within(gene_list:)
        
        int_network = Array.new
        gene_list = gene_list.map { |x| x.to_sym }
        self.interactors.each {|gene|
            int_network << gene if gene_list.include?(gene)
        }
        return int_network

    end

    def self.find_interactions(gene_list:)

        intact_hash = Hash.new
    
        gene_list.each do |gene_id|
            gene_id.upcase!       # Convert gene ID to upercase
            address = "http://www.ebi.ac.uk/Tools/webservices/psicquic/intact/webservices/current/search/interactor/#{gene_id}?format=tab25"
            response = RestClient::Request.execute(method: :get, url: address)
            intact_data = response.body.split(/\n/)

            record_list = Array.new

            intact_data.each do |record|
                intact_gene1 = record.split(/\t/)[2].split(/\|/).grep(/^ensemblplants/)  # We filter to extract ensembleplants id
                intact_gene2 = record.split(/\t/)[3].split(/\|/).grep(/^ensemblplants/)  # from columns 2 and 3 of each record.
                intact_genes = intact_gene1 + intact_gene2                               # Combine records from both columns. 

                score = record.split(/\t/).pop.split(/\:/).pop   # Extract the interaction score of each record.

                # For some genes, we obtain interaction reports for different splicing variants (.1, .2, .3)
                # For each gene record, we extract only the gene ID using regular expressions and convert gene IDs
                # into symbol format.

                intact_genes_filtered = Array.new
                intact_genes.each {|splicing| intact_genes_filtered << splicing.match(/AT\dG\d*/).to_s.to_sym}

                # Finally, we remove duplicates from our final interaction records.

                intact_genes_filtered.uniq!

            ''' With this process, we obtain a set of inteaction records for each gene, each record containing one
                or more interactors. Different records might share common interactors and even contain the query gene.
                We want to obtain a final list of unique interactors for each query gene, so we need to remove interactor
                redundancy and interaction of the query gene with itself.
            '''
                # Introduce a quality filter

                unless score.to_f < 0.3
                    record_list += intact_genes_filtered        # Combine the interactors from different records
                    record_list.uniq!                           # Remove redundant interactors
                    record_list.delete(gene_id.to_sym)          # Remove the query gene from its interactor list
                    record_list.delete("".to_sym)               # Remove empty gene IDs (observed within the results)
                    intact_hash[gene_id.to_sym] = record_list   # Append the results for every query gene to a Hash
                end
            end
    
        #print "\n\n", gene_id.to_s, " interactors: ", record_list.length
        #print "\n", Hash[gene_id.to_sym => record_list]
    
        end
        
        @@interaction_dict = intact_hash
        return intact_hash
        
        # Results are contained in intact_hash
        # We can omit genes without reported interactions:
        # clean_interactions = first_interactions.delete_if {|key,value| value.empty? }
    end

    private

    def InteractionNetwork.network(gene_id:, from_dict: @@interaction_dict, max_depth: 3)
      
        if @@interaction_dict.empty? then
            abort("Aborting: No interactions provided, find interactions before building the network.")
        end

        if @@network_depth >= max_depth
            return []
        else
            @@network_depth += 1
            puts "Finding interactions with depth = #{@@network_depth}"
            network = [gene_id.to_sym]
            interactors = @@interaction_dict[gene_id.to_sym]

            unless interactors.nil?
                network += interactors
                interactors.each {|int| network += network(gene_id: int, max_depth: max_depth)}
            end
            return network.uniq
        end
    end
    
    gene_list = Array.new
    File.foreach("./documents/ArabidopsisSubNetwork_GeneList.txt"){ |line|
        gene = line.gsub("\n",'')        # We eliminate metacharacter \n
        unless gene.match(/AT\dG\d{5}/i) # Check if genes belong to Arabidopsis and save each gene in the array.
            abort("ERROR: the gene list have some errors. #{gene} has not correct format") 
        end
        gene_list <<  gene.upcase!}
    
    
    InteractionNetwork.find_interactions(gene_list: gene_list)
    
    gene_list.each {|gene|

        network = InteractionNetwork.new(gene: gene, max_depth: 3)
        print network.query_gene
        print "\t", network.interactors.length
        print "\t", network.depth

        print "\n", network.interactors
        print "\t", @@interaction_dict[gene.to_sym].length unless @@interaction_dict[gene.to_sym].nil?
        puts 

        print "\n", network.interactors_within(gene_list: gene_list)
        puts
        puts
      }

    
end
