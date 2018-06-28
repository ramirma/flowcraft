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


def depth_file_reader(depth_file):
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
    depth_dic_coverage: dict
            dictionary with the coverage per position for each plasmid
    """

    # dict to store the mean coverage for each reference
    depth_dic_coverage = {}

    for line in depth_file:
        tab_split = line.split()  # split by any white space
        reference = "_".join(tab_split[0].strip().split("_")[0:3])  # store
        # only the gi for the reference
        position = tab_split[1]
        num_reads_align = float(tab_split[2].rstrip())

        if reference not in depth_dic_coverage:
            depth_dic_coverage[reference] = {}

        depth_dic_coverage[reference][position] = num_reads_align

    logger.info("Finished parsing depth file.")
    return depth_dic_coverage


def generate_jsons(depth_dic_coverage, plasmid_length, cutoff):
    """

    Parameters
    ----------
    depth_dic_coverage: dict
         dictionary with the coverage per position for each plasmid

    Returns
    -------
    percentage_bases_covered: dict
    dict_cov:  dict

    """

    # initializes the dictionary with the mean coverage results per plasmid
    percentage_bases_covered = {}
    # the number of points to generate the plot of coverage in flowcraft-webapp
    number_of_points = 10000
    # dict to store coverage results for a given interval of points
    dict_cov = {}
    for ref in depth_dic_coverage:
        # calculates the percentage value per each reference
        perc_value_per_ref = float(len(depth_dic_coverage[ref])) / \
                             float(plasmid_length[ref])
        # checks if percentage value is higher or equal to the cutoff defined
        if perc_value_per_ref >= float(cutoff):
            percentage_bases_covered[ref] = perc_value_per_ref

        # starts parser to get the array with the coverage for all the positions
        array_of_cov = []
        last_position = 0
        for pos in depth_dic_coverage[ref]:
            current_position = float(pos)
            # if the first element is being parsed
            if last_position == 0:
                array_of_cov.append(depth_dic_coverage[ref][pos])
                last_position = 1
            # if current_position is different from the last_position + 1
            # it means that there is a gap in the coverage information
            elif current_position != (last_position + 1):
                diff_position = int(current_position - last_position)
                array_of_cov.extend([0] * diff_position)
                last_position = current_position
            # otherwise if data is continuous just add it to array_of_cov
            else:
                array_of_cov.append(depth_dic_coverage[ref][pos])
                last_position = current_position

        # then finally if array_of_cov doesnt have elements to the end of the
        # full plasmid length, add 0 to the remaining positions
        diff_arrays = int(plasmid_length[ref] - len(array_of_cov))
        if diff_arrays:
            array_of_cov.extend([0] * diff_arrays)

        interval = round(int(plasmid_length[ref]) / number_of_points,
                         ndigits=0)
        # some plasmids can be smaller than 10000
        if interval < 1:
            interval = 1

        # starts dict cov for the reference
        dict_cov[ref] = {
            "xticks": [],
            "values": [],
        }

        counter = 1
        previous_counter = 0
        for x, entry in enumerate(array_of_cov):
            max_xtick = x + 1
            if counter == interval:
                dict_cov[ref]["xticks"].append(max_xtick)
                dict_cov[ref]["values"].append(round(
                    (sum(array_of_cov[previous_counter:max_xtick])/interval)))
                previous_counter = max_xtick
                counter = 1
            else:
                counter += 1

        dict_cov[ref]["xticks"].append(max_xtick)
        dict_cov[ref]["values"].append(round(
            (sum(array_of_cov[previous_counter:max_xtick]) / interval)))

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
    depth_dic_coverage = depth_file_reader(depth_file_in)
    percentage_bases_covered, dict_cov = generate_jsons(depth_dic_coverage,
                                                        plasmid_length, cutoff)

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
