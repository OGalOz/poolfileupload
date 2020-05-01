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


def main():

    return None

if __name__ == "__main__":
    main()
