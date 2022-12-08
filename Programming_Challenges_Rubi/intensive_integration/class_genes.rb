
require 'rest-client'
require 'json'

class GeneFeatures
    attr_accessor :gene_id 
    attr_accessor :go_terms
    attr_accessor :kegg_id
    attr_accessor :pathway

    def search_go(gene_id:)
        gene_id = gene_id.to_s
        address ="http://togows.dbcls.jp/entry/uniprot/#{@gene_id}/dr.json"
        response = RestClient::Request.execute(method: :get,  url: address) 
        go_terms=[]
        response.body.split(/\[/).grep(/GO/)[1..response.body.length].each do |go| # We eliminate the first row (it is like a header)
            go_terms << {go.gsub(/\"/,"").split(/\,/)[0].gsub(/\n/,"").strip => go.gsub(/\"/,"").split(/\,/)[1].gsub(/\n/,"").strip}
        end
        return go_terms
    end

    def search_kegg(gene_id:)
        gene_id = gene_id.to_s
        address ="http://togows.org/entry/uniprot/#{@gene_id}/dr.json"
        response = RestClient::Request.execute(method: :get, url: address)  
        JSON.parse(response.body)[0]['KEGG'].each do |kegg|
            kegg_id << kegg[0]
        end
        kegg_id.each do |id|
            address ="http://togows.org/entry/kegg-genes/#{id}.json"
            response = RestClient::Request.execute(
            method: :get, url: address)  
             JSON.parse(response.body)[0]['pathways'].each do |path|
               pathway << {path[0] => path[1]}
             end
        end
        return pathway
    end
end
