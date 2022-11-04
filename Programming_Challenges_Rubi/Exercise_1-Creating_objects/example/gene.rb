class Gene
    @@genelist = Array.new  # initialize a class variable to hold the genes
    attr_accessor :id  # create an "attribute accessor" 
  
    def initialize(id:) # get a name from the "new" call, or set a default
        @id = id  
        @@genelist.append(self)  # self is the current object
    end

    def get_genes
        return @@genelist
    end

    def self.get_genes
        return @@genelist
    end

end
