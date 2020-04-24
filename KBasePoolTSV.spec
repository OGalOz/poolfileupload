/*
This module is for storing Pool Files and PoolCount Files as used in 
Random Barcode Transposon Sequencing (RBTnSeq)
Apps: Design Random Pool, Reads to PoolCount, BarSeqR
*/
module KBasePoolTSV {

    /*
    A handle id from the Handle Service for a shock node.
    @id handle
    */
    typedef string handle_id;
    
    
    /*
    @id ws KBaseGenomes.Genome
    */
    typedef string genome_ref;
    
    
    /*
    A header for a  column
    */
    typedef string col_header;
    
    
    /*
    The list of column headers
    */
    typedef list<col_header> col_list;
    
    /*
    file_type - KBasePoolTSV.PoolFile, the name of the file type.
    poolfile - handle that allows to download file, and get info re. shock node, shock url,
    handle_type - the type of the handle. This should always be ‘shock’.
    shock_url - the url of the shock server
    shock_node_id - the id of the shock node in the server
    compression_type - the type of compression used
    file_name - the name of the file
    column_header_list - a list of the headers of the columns, the length of this 
        list should be the num of columns in the file. Currently: 
        <"barcode", "rcbarcode", "nTot", "n", "scaffold", "strand", 
        "pos", "n2", "scaffold2", "strand2", "pos2", "nPastEnd">
        making a total of 12 columns.
    run_method - keeps track of what poolfile is initially used for,
        for now it will always be "poolcount"
    related_genome_ref -  the genome which is related to the pool file.
    related_organism_scientific_name -  the related scientific_name from the genome_ref
    description - A description given by the uploader as to what the
        pool file means.
    
    @metadata ws handle_type as handle_type
    @metadata ws run_method as run_method
    @metadata ws shock_url as shock_url
    @metadata ws shock_node_id as shock_node_id
    @metadata ws related_genome_ref as related_genome_ref
    @metadata ws related_organism_scientific_name as related_organism_scientific_name
    @metadata ws description
    */
    typedef structure {
        string file_type;
        handle_id poolfile;
        string handle_type;
        string shock_url;
        string shock_node_id;
        string compression_type;
        string file_name;
        col_list column_header_list;
        string run_method;
        genome_ref related_genome_ref;
        string related_organism_scientific_name;
        string description;
    
    } PoolFile;
    
    
    /*
    file_type KBasePoolTSV.PoolCount
    handle_id will be poolcount file handle
    col_list column_header_list will be
        barcode, rcbarcode, scaffold, strand, pos, and an unknown number of columns
    
    @metadata ws handle_type as handle_type
    @metadata ws run_method as run_method
    @metadata ws shock_url as shock_url
    @metadata ws shock_node_id as shock_node_id
    @metadata ws related_genome_ref as related_genome_ref
    @metadata ws related_organism_scientific_name as related_organism_scientific_name
    @metadata ws description
    */
    typedef structure {
        string file_type;
        handle_id poolcount;
        string handle_type;
        string shock_url;
        string shock_node_id;
        string compression_type;
        string file_name;
        col_list column_header_list;
        string run_method;
        genome_ref related_genome_ref;
        string related_organism_scientific_name;
        string description;
    
    } PoolCount;

};

