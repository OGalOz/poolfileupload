#!python3

import os
import sys
import logging

def check_output_name(output_name):
    #output_name is string

    if ' ' in output_name:
        output_name.replace(' ', '_')
        logging.info("Output Name contains spaces, changing to " + output_name)
    
    return output_name

def catch_NaN(val):
    # val is a string. If not a number returns NaN
    if val in ["NaN", "NA"]:
        val = "NaN"
    else:
        try:
            val = float(val)
        except ValueError:
            val = "NaN"
        if val - math.floor(val) == 0:
            val = int(val)
    return val

# We download Genome Files: gfu is Genome File Util
def DownloadGenomeToFNA(gfu, genome_ref, scratch_dir):
    """
    Inputs: GFU Object, str (A/B/C), str path
    Outputs: [fna_fp (str), gbk_fp (str)]
    """

    GenomeToGenbankResult = gfu.genome_to_genbank({'genome_ref': genome_ref})

    logging.info("Genome File Util download result:")
    logging.info(GenomeToGenbankResult)

    genome_fna_fp = get_fa_from_scratch(scratch_dir)

    if genome_fna_fp is None:
        raise Exception("GFU Genome To Genbank did not download Assembly file in expected Manner.")

    return genome_fna_fp

def get_fa_from_scratch(scratch_dir):
    """
    Careful... May not work in the Future
    Args:
        scratch_dir: (str) Path to work dir/ tmp etc..
    Returns:
        FNA fp: (str) Automatic download through GenbankToGenome
    """
    
    fna_fp = None
    scratch_files = os.listdir(scratch_dir)
    for f in scratch_files:
        if f[-2:] == "fa":
            fna_fp = os.path.join(scratch_dir, f)
            break
    
    if fna_fp is None:
        raise Exception("Could not find Assembly FNA file in scratch (work) dir")

    return fna_fp


def GetScaffoldLengths(genome_fna_fp):
    """ This function gets the lengths of the scaffolds, returns a dict

    Args:
        genome_fna_fp: (str) Path to genome fna file (FASTA)

    Returns:
        Scaffold_To_Length: (dict) 
            scaffold_name: (str) -> length (int)
    """

    Scaffold_To_Length = {}

    FNA_FH = open(genome_fna_fp)

    c_line = FNA_FH.readline().strip()
    c_scaffold_name = ""
    while c_line != "":
        if c_line[0] == ">":
            if c_scaffold_name != "":
                Scaffold_To_Length[c_scaffold_name] = cs_len
            if " " in c_line:
                logging.warning(f"A space found in scaffold name: '{c_line}'."
                                " This might cause an error.")
                c_scaffold_name = (c_line.split(' ')[0])[1:]
                logging.warning(f"Instead using scaffold name {c_scaffold_name}")
            else:
                c_scaffold_name = c_line[1:]
            # Current scaffold length is reset
            cs_len = 0
        else:
            cs_len += len(c_line)

        c_line = FNA_FH.readline().strip()

    FNA_FH.close()

    if c_scaffold_name != "":
        Scaffold_To_Length[c_scaffold_name] = cs_len

    if len(Scaffold_To_Length.keys()) == 0:
        logging.warning("No Scaffolds found in " + genome_fna_fp)

    return Scaffold_To_Length





def main():

    return None

if __name__ == "__main__":
    main()
