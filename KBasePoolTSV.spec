/*
This module is for storing Pool Files, PoolCount Files, and Experiment Files 
as used in Random Barcode Transposon Sequencing (RBTnSeq)
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
    @id ws KBaseFile.SingleEndLibrary 
    */
    typedef string fastq_ref;
   
    /*
    A list of fastq_refs
    */
    typedef list<fastq_ref> fastqs;

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
    utc_created - the Coordinated Universal Time of creation
    column_header_list - a list of the headers of the columns, the length of this 
        list should be the num of columns in the file. Currently: 
        <"barcode", "rcbarcode", "nTot", "n", "scaffold", "strand", 
        "pos", "n2", "scaffold2", "strand2", "pos2", "nPastEnd">
        making a total of 12 columns.
    num_lines - the number of lines in the file - keeps track of the general size
    related_genome_ref -  the genome which is related to the pool file.
    related_organism_scientific_name -  the related scientific_name from the genome_ref
    fastqs_used - the fastqs which were used to create the poolfile
    description - A description given by the uploader as to what the
        pool file means.
    
    @optional fastqs_used
    @metadata ws utc_created as utc_created
    @metadata ws handle_type as handle_type
    @metadata ws shock_url as shock_url
    @metadata ws shock_node_id as shock_node_id
    @metadata ws num_lines 
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
        string utc_created;
        col_list column_header_list;
        string num_lines;
        genome_ref related_genome_ref;
        string related_organism_scientific_name;
        fastqs fastqs_used; 
        string description;
    
    } PoolFile;
    
    
    /*
    file_type KBasePoolTSV.PoolCount
    handle_id will be poolcount file handle
    handle_type - the type of the handle. This should always be ‘shock’.
    col_list column_header_list will be
        barcode, rcbarcode, scaffold, strand, pos, and an unknown number of columns
    shock_url - the url of the shock server
    shock_node_id - the id of the shock node in the server
    compression_type - the type of compression used
    file_name - the name of the file
    utc_created - the Coordinated Universal Time of creation
    set_name - the name of the set
    num_lines - the number of lines in the file - keeps track of the general size
    related_genome_ref -  the genome which is related to the pool file.
    related_organism_scientific_name -  the related scientific_name from the genome_ref
    fastqs_used - the fastqs used to create the poolcount file
    poolfile_ref - the ref for the poolfile used to create the poolcount file
    description - A description given by the uploader as to what the
        pool file means.
    
    @optional poolfile_ref fastqs_used
    @metadata ws utc_created as utc_created
    @metadata ws handle_type as handle_type
    @metadata ws shock_url as shock_url
    @metadata ws shock_node_id as shock_node_id
    @metadata ws set_name
    @metadata ws related_genome_ref as related_genome_ref
    @metadata ws related_organism_scientific_name as related_organism_scientific_name
    @metadata ws description
    @metadata ws num_lines
    */
    typedef structure {

        string file_type;
        handle_id poolcount;
        string handle_type;
        col_list column_header_list;
        string shock_url;
        string shock_node_id;
        string compression_type;
        string file_name;
        string utc_created;
        string set_name;
        string num_lines;
        genome_ref related_genome_ref;
        string related_organism_scientific_name;
        fastqs fastqs_used; 
        string poolfile_ref;
        string description;
    
    } PoolCount;


    /*
    file_type KBasePoolTSV.Experiments
    exps_file will be experiments file handle
    handle_type - the type of the handle. This should always be ‘shock’.
    col_list column_header_list will have required parts:
        SetName, Index, Description, Date_pool_expt_started
        Not required headers have the following sometimes:
            Drop, Person, Mutant Library,
            gDNA plate, gDNA well, Sequenced at, Media, Growth Method,
            Group, Temperature, pH, Liquid v. solid, Aerobic_v_Anaerobic, Shaking,
            Condition_1, Concentration_1, Units_1, Condition_2, Concentration_2,
            Units_2, Timecourse, Timecourse Sample, Growth Plate ID, 
            Growth Plate wells, StartOK, EndOD, Total Generations
    shock_url - the url of the shock server
    shock_node_id - the id of the shock node in the server
    compression_type - the type of compression used
    file_name - the name of the file
    utc_created - the Coordinated Universal Time of creation
    num_lines - the number of lines in the file - keeps track of the general size
    related_genome_ref -  the genome which is related to the pool file.
    related_organism_scientific_name -  the related scientific_name from the genome_ref
    description - A description given by the uploader as to what the
        pool file means.

    @optional poolfile_ref
    @metadata ws utc_created as utc_created
    @metadata ws handle_type as handle_type
    @metadata ws shock_url as shock_url
    @metadata ws shock_node_id as shock_node_id
    @metadata ws related_genome_ref as related_genome_ref
    @metadata ws related_organism_scientific_name as related_organism_scientific_name
    @metadata ws description
    @metadata ws num_lines
    */
    typedef structure {

        string file_type;
        handle_id expsfile;
        string handle_type;
        string shock_url;
        string shock_node_id;
        string compression_type;
        string file_name;
        string utc_created;
        col_list column_header_list;
        string num_lines;
        genome_ref related_genome_ref;
        string related_organism_scientific_name;
        string poolfile_ref;
        string description;

    } Experiments;

};

