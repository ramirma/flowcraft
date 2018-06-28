#!/usr/bin/env python3


"""
Purpose
-------

This module is intended to generate a json output for mapping results that
can be imported in pATLAS.

Expected input
--------------

The following variables are expected whether using NextFlow or the
:py:func:`main` executor.

- ``depth_file`` : String with the name of the mash screen output file.
    - e.g.: ``'samtoolsDepthOutput_sampleA.txt'``
- ``json_dict`` : the file that contains the dictionary with keys and values for
        accessions and their respective lengths.
    - e.g.: ``'reads_sample_result_length.json'``
- ``cutoff`` : The cutoff used to trim the unwanted matches for the minimum
        coverage results from mapping. This value may range between 0 and 1.
    - e.g.: ``0.6``


Code documentation
------------------

"""

__version__ = "1.0.1"
__build__ = "20022018"
__template__ = "mapping2json-nf"

import os
import json
import sys

from flowcraft_utils.flowcraft_base import get_logger, MainWrapper

logger = get_logger(__file__)

if __file__.endswith(".command.sh"):
    DEPTH_TXT = '$depthFile'
    JSON_LENGTH = '$lengthJson'
    SAMPLE_ID = '$sample_id'
    CUTOFF = '$params.cov_cutoff'
else:
    DEPTH_TXT = sys.argv[1]
    JSON_LENGTH = sys.argv[2]
    SAMPLE_ID = sys.argv[3]
    CUTOFF = sys.argv[4]

# check if all variables are assigned
if DEPTH_TXT and JSON_LENGTH and SAMPLE_ID and CUTOFF:
    logger.debug("Running {} with parameters:".format(
        os.path.basename(__file__)))
    logger.debug("DEPTH_TXT: {}".format(DEPTH_TXT))
    logger.debug("JSON_LENGHT: {}".format(JSON_LENGTH))
    logger.debug("CUTOFF: {}".format(CUTOFF))
else:
    logger.error("Args should be given to this template, either from sys.argv"
                 " or through nextflow variables")


def depth_file_reader(depth_file, plasmid_length, cutoff):
    """
    Function that parse samtools depth file and creates 3 dictionaries that
    will be useful to make the outputs of this script, both the tabular file
    and the json file that may be imported by pATLAS

    Parameters
    ----------
    depth_file: textIO
        the path to depth file for each sample
    plasmid_length: dict
        a dictionary that stores length of all plasmids in fasta given as input
    cutoff: float
        the cutoff used to trim the unwanted matches for the minimum coverage
        results from mapping. This is then converted into a float within this
        function in order to compare with the value returned from the
        perc_value_per_ref.

    Returns
    -------
    percentage_bases_covered: dict
            stores the percentage of the total sequence of a
            reference/accession (plasmid) in a dictionary
    """

    # dict to store the mean coverage for each reference
    depth_dic_coverage = {}
    # the number of points to generate the plot of coverage in flowcraft-webapp
    number_of_points = 10000
    # dict to store coverage results for a given interval of points
    dict_cov = {}
    # temporary array to store coverage values for a given interval
    array_of_cov = []

    counter = 0
    times = 1
    for line in depth_file:
        tab_split = line.split()  # split by any white space
        reference = "_".join(tab_split[0].strip().split("_")[0:3])  # store
        # only the gi for the reference
        position = tab_split[1]
        num_reads_align = float(tab_split[2].rstrip())

        if reference not in depth_dic_coverage:
            depth_dic_coverage[reference] = {}
            # get the number of intervals for this reference sequence
            interval = int(plasmid_length[reference]) / number_of_points
            # some plasmids can be smaller than 10000
            if interval < 1:
                interval = 1
            dict_cov[reference] = {
                "xticks": [],
                "values": [],
            }

        depth_dic_coverage[reference][position] = num_reads_align

        if interval and float(counter) < interval:
            counter += 1
            array_of_cov.append(num_reads_align)
        else:
            array_of_cov.append(num_reads_align)
            dict_cov[reference]["xticks"].append(interval * times)
            dict_cov[reference]["values"].append(
                round(sum(array_of_cov) / len(array_of_cov), 2)
            )
            counter = 0
            times += 1

    logger.info("Finished parsing depth file.")

    percentage_bases_covered = {}
    for ref in depth_dic_coverage:
        # calculates the percentage value per each reference
        perc_value_per_ref = float(len(depth_dic_coverage[ref])) / \
                             float(plasmid_length[ref])
        # checks if percentage value is higher or equal to the cutoff defined
        if perc_value_per_ref >= float(cutoff):
            percentage_bases_covered[ref] = perc_value_per_ref

    logger.info("Successfully generated dicts necessary for output json file "
                "and .report.json depth file.")

    return percentage_bases_covered, dict_cov


@MainWrapper
def main(depth_file, json_dict, cutoff, sample_id):
    """
    Function that handles the inputs required to parse depth files from bowtie
    and dumps a dict to a json file that can be imported into pATLAS.

    Parameters
    ----------
    depth_file: str
         the path to depth file for each sample
    json_dict: str
        the file that contains the dictionary with keys and values for
        accessions
        and their respective lengths
    cutoff: str
        the cutoff used to trim the unwanted matches for the minimum coverage
        results from mapping. This value may range between 0 and 1.
    sample_id: str
        the id of the sample being parsed

    """

    # check for the appropriate value for the cutoff value for coverage results
    try:
        cutoff_val = float(cutoff)
    except ValueError:
        logger.error("Cutoff value should be a string such as: '0.6'. "
                     "The outputted value: {}. Make sure to provide an "
                     "appropriate value for --cov_cutoff".format(cutoff))

    # loads dict from file, this file is provided in docker image

    plasmid_length = json.load(open(json_dict))
    if plasmid_length:
        logger.info("Loaded dictionary of plasmid lengths")
    else:
        logger.error("Something went wrong and plasmid lengths dictionary"
                     "could not be loaded. Check if process received this"
                     "param successfully.")
        sys.exit(1)

    # read depth file
    depth_file_in = open(depth_file)

    # first reads the depth file and generates dictionaries to handle the input
    # to a simpler format
    logger.info("Reading depth file and creating dictionary to dump.")
    percentage_bases_covered, dict_cov = depth_file_reader(depth_file_in,
                                                           plasmid_length,
                                                           cutoff_val)

    if percentage_bases_covered and dict_cov:
        logger.info("percentage_bases_covered length: {}"
                    "dict_cov length: {}".format(
                        str(len(percentage_bases_covered)),
                        str(len(dict_cov))
                    ))
    else:
        logger.error("Both dicts that dump to JSON file or .report.json are "
                     "empty.")

    # then dump do file
    output_json = open("{}_mapping.json".format(depth_file), "w")
    logger.info("Dumping to {}".format("{}_mapping.json".format(depth_file)))
    output_json.write(json.dumps(percentage_bases_covered))
    output_json.close()

    json_dic = {
        "sample": sample_id,
        "patlas_mapping": percentage_bases_covered,
        "plotData": {
            "header": "Coverage results for plasmids",
            "data": dict_cov,
            "plot": "mappingPlasmids"
        }
    }

    logger.info("Writting to .report.json")
    with open(".report.json", "w") as json_report:
        json_report.write(json.dumps(json_dic, separators=(",", ":")))


if __name__ == "__main__":
    main(DEPTH_TXT, JSON_LENGTH, CUTOFF, SAMPLE_ID)
