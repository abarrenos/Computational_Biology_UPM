require 'rest-client'

# == InteractionNetwork
#
# InteractionNetwork class allows us to identify full Interaction Networks in which input genes parcipate and 
# local interactions within a particular gene set.
#
# == Summary
#
# For each query gene, interactors are fetched from Psiquic IntAct database and appended to an interaction Hash. 
# This interaction information is used to elaborate complex Interaction Networks through a recursive function. 
# Interaction Networks can be compared to an input gene list to analyse gene interactions within a particular set of genes. 
#

class InteractionNetwork

    # @!attribute [rw]
    # @return [Symbol] The query gene
    attr_accessor :query_gene
    # @!attribute [rw]
    # @return [Array<Symbol>] The full network of interactors for the query gene in Symbol format
    attr_accessor :interactors
    # @!attribute [rw]
    # @return [Array<Integer>] The maximum depth of the network
    attr_accessor :max_depth

    # The Hash containing the whole set of query genes and their interactors.
    @@interaction_dict = Hash.new
    
    # An Array containing the whole set of Interaction Networks built.
    @@all_networks = Array.new

    # Create a new instance of the complete InteractionNetwork given a query gene, a Hash with gene interaction information (default = @@interaction_dict), and a maximum network depth.
    # @param gene [Symbol] the query gene id.
    # @param interaction_dict [Hash] set of query genes and their interactors.
    # @param max_depth [Integer] the maximum depth of the network. Consecutive neighboring nodes to be analysed.
    # @return [InteractionNetwork] the query gene Interaction Network instance.
    def initialize(gene:, interaction_dict: @@interaction_dict, max_depth: 3)
    
        @query_gene = gene.to_sym
        @interactors = InteractionNetwork.network(gene_id: gene, \
                                                from_hash: interaction_dict, \
                                                max_depth: max_depth)
        @max_depth = max_depth
        @@all_networks << self
    
    end
    
    # Return the Array containing the whole set of Interaction Networks built.
    # @return [Array<InteractionNetwork>] all Interaction Networks.
    def self.get_all_networks
        return @@all_networks
    end
    
    # Analyse which genes conforming the Interaction Network are found within an input gene list.
    # @return [Array<Symbol>] interactions within an input gene list.
    def interactors_within(gene_list:)

        int_network = Array.new
        gene_list = gene_list.map { |x| x.to_sym }
        self.interactors.each {|gene|
            int_network << gene if gene_list.include?(gene)
        }
        return int_network

    end

    # For a given gene list, interactors for each gene are fetched from Psiquic IntAct database and appended to an interaction Hash.
    # @param gene_list [Array<String>] the list with query gene IDs.
    # @return [Hash] the Hash with interaction information for each gene.
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
    
    # Private method only called from the Class initializer. It generated an Array containing the genes that conform the Interaction Network given an input gene id.
    # @param gene:id [Symbol] the query gene id.
    # @param from_hash [Hash] set of query genes and their interactors (default = @@interaction_dict).
    # @param max_depth [Integer] the maximum depth of the network. Consecutive neighboring nodes to be analysed.
    # @param analysed_genes [Array] an array to keep track of the genes that have already been analysed.
    # @return [Array] the genes conforming the Interaction Network.
    def InteractionNetwork.network(gene_id:, from_hash: @@interaction_dict, max_depth: 3, analysed_genes: Array.new)
      
        if from_hash.empty? then
            abort("Aborting: No interaction Hash provided, find interactions before building the network.")
        end
                
        gene_id = gene_id.to_sym        # Convert gene_id to symbol

        network = [gene_id]             # Create an interaction network array

        # If the current gene is not in the interaction hash or the list of interactors is empty return the network array
        return network if !(from_hash.include?(gene_id)) || from_hash[gene_id].empty?
 
        analysed_genes << gene_id

        # Get the list of genes that interact with the current gene
        interactors = from_hash[gene_id]

        # If the depth is greater than 0, continue recursively calling the function.
        if max_depth > 0        
            interactors.each { |interactor|

                # In the interactors have not been previously analysed, run the function recursively over them
                unless analysed_genes.include?(interactor)
                    network += network(gene_id: interactor, max_depth: max_depth-1,\
                                        analysed_genes: analysed_genes)
                end }                 
        end
        return network
    end
        
end
