#!/usr/bin/env python

"""
assembly_stat

Report assembly file (fasta format) statistics


Created by Tae-Hyuk (Ted) Ahn on 01/15/2014.
Copyright (c) 2013 Tae-Hyuk Ahn (ORNL). Allrights reserved.
"""


import sys, warnings, os, re
from datetime import datetime, date, time
from subprocess import Popen, PIPE, check_call, STDOUT
import getopt
from Bio import SeqIO
from bitarray import bitarray

## Version control
version = "1.0.1"

## Help message
help_message = '''

  [Usage]
    assembly_stat denovo -i <contigs fasta file> 
    assembly_stat mapped -i <contigs fasta file> -r <reference fasta file>

  [Requirements]
    Biopython (>= V1.6.0)
    bitarray python module

  [Options]
    1. If you work with de-novo assembly contigs:
    -m <int>    : minimum length of fragments
                    (Default: 0 --> it means that no filtering will be processed in default)

    2. If you work with BWA mapped assembly contigs:
    -m <int>    : minimum length of fragments
                    (Default: 0 --> it means that no filtering will be processed in default)
    -q <float>  : filter out fragments that have lower mapping rate using edit distance
                    Ex) if -q 0.95, then only consider 1 - (EditDistance/MappedLength) >= 0.95 mapped fragments
                    (Default: 0 --> it means that no filtering will be processed in default)

    -h/--help
    -v/--version

  [Output]
    BaseFilename.filtered.fasta (if -m or -q option is provided)
    BaseFilename.stat.txt
'''

## Class Usage
class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg 


## Exit system with error message
def die(msg=None):
    if msg is not None:
        print >> sys.stderr, msg
        sys.exit(1)


## string class
class OptionString:

     program_mode_str       = "program_mode"
     input_filename_str     = "input_filename"
     ref_filename_str       = "ref_filename"
     map_quality_cutoff_str = "map_quality_cutoff"
     min_length_cutoff_str  = "min_length_cutoff"


## Parse options
def parse_option(argv, program_mode, version, help_message):

    # option map dictionay
    option_map = {}

    if program_mode == "denovo":
        option_map[OptionString.program_mode_str] = "denovo"
    elif program_mode == "mapped":
        option_map[OptionString.program_mode_str] = "mapped"
    else:
        raise Usage(help_message)

    # get arguments
    try:
        opts, args = getopt.getopt(argv[2:], "hvVi:r:q:m:",
                                   ["help","version"])

    # Error handling of options
    except getopt.error as msg:
        raise Usage(msg)

    # get program name
    program_name =  sys.argv[0].split("/")[-1]

    # Basic options
    for option, value in opts:
        # help
        if option in ("-h", "--help"):
            raise Usage(help_message)
        # version
        if option in ("-v", "-V", "--version"):
            #print "\n%s V%s\n" % (program_name, version)
            sys.exit(0)
        # -i
        if option in ("-i"):
            option_map[OptionString.input_filename_str] = value
        # -r
        if option in ("-r"):
            option_map[OptionString.ref_filename_str] = value
        # -q: for map quality cutoff
        if option in ("-q"):
            option_map[OptionString.map_quality_cutoff_str] = value
        # -m: for minimum length cutoff
        if option in ("-m"):
            option_map[OptionString.min_length_cutoff_str] = value

    # check you got input file
    if not OptionString.input_filename_str in option_map:
        raise Usage(help_message)

    if program_mode == "mapped":
        if not OptionString.ref_filename_str in option_map:
            raise Usage(help_message)

    return (option_map)


def get_N_Length(teo_N, contigs_length):
    N_val = 0
    test_sum = 0
    for contig_length in contigs_length:
        test_sum += contig_length
        if teo_N < test_sum:
            N_val = contig_length
            break
    nN_val = 1
    for index, contig_length in enumerate(contigs_length):
        if contig_length == N_val:
            nN_val = index
            break

    return N_val,nN_val


def get_N80(sum_val, contigs_length):
    N80_val = 0
    teo_N80 = sum_val * 0.8
    test_sum = 0
    for contig_length in contigs_length:
        test_sum += contig_length
        if teo_N80 < test_sum:
            N80_val = contig_length
            break

    return N80_val


def get_N50(sum_val, contigs_length):
    N50_val = 0
    teo_N50 = sum_val * 0.5
    test_sum = 0
    for contig_length in contigs_length:
        test_sum += contig_length
        if teo_N50 < test_sum:
            N50_val = contig_length
            break

    return N50_val


def get_nN50(N50_val, contigs_length):
    nN50_val = 1
    for index, contig_length in enumerate(contigs_length):
        if contig_length == N50_val:
            nN50_val = index
            break

    return nN50_val


def get_N20(sum_val, contigs_length):
    N20_val = 0
    teo_N20 = sum_val * 0.2
    test_sum = 0
    for contig_length in contigs_length:
        test_sum += contig_length
        if teo_N20 < test_sum:
            N20_val = contig_length
            break

    return N20_val


def work_denovo(option_map):

    # input
    input_filename = option_map[OptionString.input_filename_str]

    # get base and output filename
    base_filename = os.path.basename(input_filename)
    #dir_filename = os.getcwd()
    dir_filename = os.path.dirname(input_filename)
    base_wo_ext = os.path.splitext(base_filename)[0]
    output_stat_filename = base_wo_ext + ".stat.txt"
    output_assembly_filename = base_wo_ext + ".filtered.fasta"
    if dir_filename != "":
        output_stat_filename = dir_filename + "/" + base_wo_ext + ".stat.txt"
        output_assembly_filename = dir_filename + "/" + base_wo_ext + ".filtered.fasta"

    # minimum length
    min_length_cutoff = 0
    if OptionString.min_length_cutoff_str in option_map:
        try:
            min_length_cutoff = int(option_map[OptionString.min_length_cutoff_str])
        except:
            die("Check -m option value whether it is integer value!")

    # open stat output 
    output_stat = open(output_stat_filename, 'wb')

    # open assembly output if map_quality_cutoff != 0.0
    if (min_length_cutoff != 0):
        output_assembly = open(output_assembly_filename, 'wb')

    # initialize 
    sum_val = 0
    contigs_length = []

    # loop fasta
    for seq_record in SeqIO.parse(open(input_filename), "fasta"):
        seq_description_list = seq_record.description.split(",")

        # check contigs length
        if len(seq_record.seq) >= min_length_cutoff:
            sum_val += len(seq_record.seq)
            contigs_length.append(len(seq_record.seq))
            if (min_length_cutoff != 0):
                SeqIO.write(seq_record, output_assembly, "fasta")

    # if map_quality_cutoff != 0.0, 
    if min_length_cutoff != 0:
        output_assembly.close()

    # sort list
    contigs_length.sort()
    contigs_length.reverse()

    # initialize
    n_val = 0
    min_val = 0
    max_val = 0
    N80_val = 0
    N50_val = 0
    nN50_val = 0
    N20_val = 0

    # get statistics
    if len(contigs_length) > 0:
        n_val = len(contigs_length)
        min_val = min(contigs_length)
        max_val = max(contigs_length)
        N80_val = get_N80(sum_val, contigs_length)
        N50_val = get_N50(sum_val, contigs_length)
        nN50_val = get_nN50(N50_val, contigs_length)
        N20_val = get_N20(sum_val, contigs_length)
        N50M_val, n50M_val = get_N_Length(50000000, contigs_length)
        N100M_val, n100M_val = get_N_Length(100000000, contigs_length)
        N200M_val, n200M_val = get_N_Length(200000000, contigs_length)
        N300M_val, n300M_val = get_N_Length(300000000, contigs_length)
        N500M_val, n500M_val = get_N_Length(500000000, contigs_length)
        N800M_val, n800M_val = get_N_Length(800000000, contigs_length)
        N1000M_val, n1000M_val = get_N_Length(1000000000, contigs_length)
    
    # print to termial and file
    #print "n\tn:N50\tmin\tN80\tN50\tN20\tN50M\tN100M\tN200M\tN300M\tN500M\tN800M\tN1000M\tmax\tsum\tfilename"
    #print str(n_val)+"\t"+str(nN50_val)+"\t"+str(min_val)+"\t"+str(N80_val)+"\t"+str(N50_val)+"\t"+str(N20_val)+"\t"+\
          #  str(N50M_val)+"\t"+str(N100M_val)+"\t"+str(N200M_val)+"\t"+str(N300M_val)+"\t"+str(N500M_val)+"\t"+str(N800M_val)+"\t"+str(N1000M_val)+"\t"\
          #  +str(max_val)+"\t"+str(sum_val)+"\t"+input_filename
    #output_stat.write("n\tn:N50\tmin\tN80\tN50\tN20\tmax\tsum\tfilename\n")
    #output_stat.write(str(n_val)+"\t"+str(nN50_val)+"\t"+str(min_val)+"\t")
    #output_stat.write(str(N80_val)+"\t"+str(N50_val)+"\t"+str(N20_val)+"\t")
    #output_stat.write(str(max_val)+"\t"+str(sum_val)+"\t"+input_filename+"\n")
    output_stat.write("n\tn:N50\tmin\tN80\tN50\tN20\tN50M\tN100M\tN200M\tN300M\tN500M\tN800M\tN1000M\tmax\tsum\tfilename\n")
    output_stat.write(str(n_val)+"\t"+str(nN50_val)+"\t"+str(min_val)+"\t"+str(N80_val)+"\t"+str(N50_val)+"\t"+str(N20_val)+"\t"+\
            str(N50M_val)+"\t"+str(N100M_val)+"\t"+str(N200M_val)+"\t"+str(N300M_val)+"\t"+str(N500M_val)+"\t"+str(N800M_val)+"\t"+str(N1000M_val)+"\t"\
            +str(max_val)+"\t"+str(sum_val)+"\t"+input_filename+"\n")

def work_mapped(option_map):

    # input
    input_filename = option_map[OptionString.input_filename_str]
    ref_filename = option_map[OptionString.ref_filename_str]

    # get base and output filename
    base_filename = os.path.basename(input_filename)
    dir_filename = os.path.dirname(input_filename)
    base_wo_ext = os.path.splitext(base_filename)[0]
    output_stat_filename = base_wo_ext + ".stat.txt"
    output_assembly_filename = base_wo_ext + ".filtered.fasta"
    if dir_filename != "":
        output_stat_filename = dir_filename + "/" + base_wo_ext + ".stat.txt"
        output_assembly_filename = dir_filename + "/" + base_wo_ext + ".filtered.fasta"

    # qual
    map_quality_cutoff = 0.0
    if OptionString.map_quality_cutoff_str in option_map:
        try:
            map_quality_cutoff = float(option_map[OptionString.map_quality_cutoff_str])
        except:
            die("Check -q option value whether it float value!")
            exit(1)

    # minimum length
    min_length_cutoff = 0
    if OptionString.min_length_cutoff_str in option_map:
        try:
            min_length_cutoff = int(option_map[OptionString.min_length_cutoff_str])
        except:
            die("Check -m option value whether it is integer value!\n")

    # open stat output 
    output_stat = open(output_stat_filename, 'wb')

    # open assembly output if map_quality_cutoff != 0.0
    if (map_quality_cutoff != 0.0) or (min_length_cutoff != 0):
        output_assembly = open(output_assembly_filename, 'wb')

    # initialize 
    sum_val = 0
    contigs_length = []
    ref_id_list = []
    ref_len_list = []

    # loop reference fasta
    for seq_record in SeqIO.parse(open(ref_filename), "fasta"):
        ref_id_list.append(seq_record.id)
        ref_len_list.append(len(seq_record.seq))

    # initialize 2D arr including bitarray vectors
    ref_arrays = [bitarray(ref_len_list[idx]) for idx, ref_id in enumerate(ref_id_list)]
    for idx, ref_array in enumerate(ref_arrays):
        ref_arrays[idx].setall(False)

    # loop contigs fasta
    for seq_record in SeqIO.parse(open(input_filename), "fasta"):
        seq_description_list = seq_record.description.split(",")

        # convertBWA format result
        if len(seq_description_list) == 10:
            # Col1: ContigName=contig_11,
            # Col2: ContigLength=525700,
            # Col3: MappedStartPositionInContig=524968,
            # Col4: ReferenceID=gi|283856168|ref|NC_006526.2|
            # Col5: MappedStartPositionInReference=2437351,
            # Col6: MappedLength=388,
            # Col7: Insertion=1,
            # Col8: Deletion=1,
            # Col9: Mismatch=49,
            # Col10: EditDistance=51
            mapped_length = int(seq_description_list[5].split("=")[1])
            if mapped_length < 1:
                mapped_length = 1
            edit_distance = int(seq_description_list[9].split("=")[1])

            reference_id  = seq_description_list[3].split("=")[1]
            mapped_start_ref = int(seq_description_list[4].split("=")[1])

            # if map_quality_cutoff = 0.0 and min_length_cutoff = 0, 
            # then consider 1 - (EditDistance/MappedLength) >= 0.95 mapped fragments
            # and consider (len(seq_record.seq) >= min_length_cutoff)
            if (map_quality_cutoff != 0.0) and (min_length_cutoff != 0):
                if (1 - (float(edit_distance)/float(mapped_length)) >= map_quality_cutoff):
                    if (len(seq_record.seq) >= min_length_cutoff):
                        # update sum and contigs length list
                        sum_val += len(seq_record.seq)
                        contigs_length.append(len(seq_record.seq))
                        # make "True" for mapped regions
                        target_ref_idx = ref_id_list.index(reference_id)
                        ref_arrays[target_ref_idx][mapped_start_ref-1:mapped_start_ref-1+mapped_length]=True
                        # write filtered fasta
                        SeqIO.write(seq_record, output_assembly, "fasta")

            # only consider 1 - (EditDistance/MappedLength) >= 0.95 mapped fragments
            elif (map_quality_cutoff != 0.0):
                if (1 - (float(edit_distance)/float(mapped_length)) >= map_quality_cutoff):
                    # update sum and contigs length list
                    sum_val += len(seq_record.seq)
                    contigs_length.append(len(seq_record.seq))
                    # make "True" for mapped regions
                    target_ref_idx = ref_id_list.index(reference_id)
                    ref_arrays[target_ref_idx][mapped_start_ref-1:mapped_start_ref-1+mapped_length]=True
                    # write filtered fasta
                    SeqIO.write(seq_record, output_assembly, "fasta")
            # only consider sequence length cutoff
            elif (min_length_cutoff != 0.0):
                if len(seq_record.seq) >= min_length_cutoff:
                    # update sum and contigs length list
                    sum_val += len(seq_record.seq)
                    contigs_length.append(len(seq_record.seq))
                    # make "True" for mapped regions
                    target_ref_idx = ref_id_list.index(reference_id)
                    ref_arrays[target_ref_idx][mapped_start_ref-1:mapped_start_ref-1+mapped_length]=True
                    # write filtered fasta
                    SeqIO.write(seq_record, output_assembly, "fasta")
            # no options
            else:
                # update sum and contigs length list
                sum_val += len(seq_record.seq)
                contigs_length.append(len(seq_record.seq))
                # make "True" for mapped regions
                target_ref_idx = ref_id_list.index(reference_id)
                ref_arrays[target_ref_idx][mapped_start_ref-1:mapped_start_ref-1+mapped_length]=True
            
        # de novo assembly fasta result
        else:
            die("Check contig fasta file.ID should have 10 elements!")

    # if map_quality_cutoff != 0.0, 
    if (map_quality_cutoff != 0.0) or (min_length_cutoff != 0):
        output_assembly.close()

    # sort list
    contigs_length.sort()
    contigs_length.reverse()

    # initialize
    n_val = 0
    min_val = 0
    max_val = 0
    N80_val = 0
    N50_val = 0
    nN50_val = 0
    N20_val = 0
    cov_rate_pct = 0

    # get statistics
    if len(contigs_length) > 0:
        n_val = len(contigs_length)
        min_val = min(contigs_length)
        max_val = max(contigs_length)
        N80_val = get_N80(sum_val, contigs_length)
        N50_val = get_N50(sum_val, contigs_length)
        nN50_val = get_nN50(N50_val, contigs_length)
        N20_val = get_N20(sum_val, contigs_length)
    
    # get coverage
    total_ref_seq_len = sum(ref_len_list)
    total_ref_cov_len = 0
    for ref_array in ref_arrays:
        total_ref_cov_len += sum(ref_array)
    if float(total_ref_seq_len) != 0:
        cov_rate_pct = 100*(float(total_ref_cov_len)/float(total_ref_seq_len))

    # print to termial and file
    #print "n\tn:N50\tmin\tN80\tN50\tN20\tN50M\tN100M\tN200M\tN300M\tmax\tsum\tcoverage(%)\tfilename"
    #print str(n_val)+"\t"+str(nN50_val)+"\t"+str(min_val)+"\t"+str(N80_val)+"\t"+str(N50_val)+"\t"+str(N20_val)+"\t"+\
     #       str(N50M_val)+"\t"+str(N100M_val)+"\t"+str(N200M)+"\t"+str(N300M)+"\t"\
      #      +str(max_val)+"\t"+str(sum_val)+"\t"+str("{0:.2f}".format(cov_rate_pct))+"\t"+input_filename
    output_stat.write("n\tn:N50\tmin\tN80\tN50\tN20\tmax\tsum\tcoverage(%)\tfilename\n")
    output_stat.write(str(n_val)+"\t"+str(nN50_val)+"\t"+str(min_val)+"\t")
    output_stat.write(str(N80_val)+"\t"+str(N50_val)+"\t"+str(N20_val)+"\t")
    output_stat.write(str(max_val)+"\t"+str(sum_val)+"\t"+str("{0:.2f}".format(cov_rate_pct))+"\t"+input_filename+"\n")


def main(argv=None):

    # try to get arguments and error handling
    try:
        if argv is None:
            argv = sys.argv
        try:
            if len(sys.argv) < 2:
                raise Usage(help_message)

            # get program name
            program_name =  os.path.basename(sys.argv[0])
            program_mode =  os.path.basename(sys.argv[1])

            # parse option
            (option_map) = parse_option(argv, program_mode, version, help_message)

            # display work start and time record
            start_time = datetime.now()
            sys.stdout.write("\n=================================================================================\n")
            sys.stdout.write("** Beginning %s run (V%s)\n" % ( program_name, version))

            # work
            sys.stdout.write("---------------------------------------------------------------------------------\n")
            if program_mode == "denovo":
                work_denovo(option_map)
            elif program_mode == "mapped":
                work_mapped(option_map)
            else:
                raise Usage(help_message)
            sys.stdout.write("---------------------------------------------------------------------------------\n")

            # time record, calculate elapsed time, and display work end
            finish_time = datetime.now()
            duration = finish_time - start_time
            sys.stdout.write("** Ending %s run\n" % (program_name))
            sys.stdout.write("** Total Elapsed Time =  %s [seconds]\n" % (duration))
            sys.stdout.write("=================================================================================\n\n")


        # Error handling
        except Usage as err:
            sys.stderr.write("%s: %s\n" %(os.path.basename(sys.argv[0]), str(err.msg)))
            return 2


    # Error handling
    except Usage as err:
        sys.stderr.write("%s: %s\n" %(os.path.basename(sys.argv[0]), str(err.msg)))
        sys.stderr.write("for help use -h/--help")
        return 2


## If this program runs as standalone, then exit.
if __name__ == "__main__":
    sys.exit(main())