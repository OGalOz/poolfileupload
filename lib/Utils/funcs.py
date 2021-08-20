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



def main():

    return None

if __name__ == "__main__":
    main()
