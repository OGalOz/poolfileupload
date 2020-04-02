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
typedef string col_header


/*
The list of column headers
*/
typedef list<col_header> col_list;

/*
file_type - KBaseRBTnSeq.PoolTSV, the name of the file type.
handle_id - allows to get info re. shock node, shock url,
handle_type - the type of the handle. This should always be ‘shock’.
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
@metadata ws related_genome_ref as related_genome_ref
@metadata ws related_organism_scientific_name as related_organism_scientific_name
*/
typedef structure {
    string file_type;
    handle_id hid;
    string handle_type;
    string compression_type;
    string file_name;
    col_list column_header_list;
    string run_method;
    genome_ref related_genome_ref;
    string related_organism_scientific_name;
    string description;

} PoolFile;
