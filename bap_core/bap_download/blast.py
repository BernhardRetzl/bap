import glob
import numpy as np
import subprocess
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML


def collect_sequences_for_blast_db(in_path, out_file,  with_limits=False, minimum_values=None, maximum_values=None):
    file_list = glob.glob(in_path)
    to_write = list()
    if with_limits:
        for file in file_list:
            handler = open(file)
            info = handler.readline().strip()
            sequence = handler.readline().strip()
            while info != '':
                spacing = info.split(' ')[-1].split('-')[1:-1]
                spacing = np.array([int(i) for i in spacing])
                spacing_excess_length = len(spacing) - 4
                for i in range(spacing_excess_length):
                    new_spacing = spacing[0 + i:5 + i]
                    if all(new_spacing >= minimum_values) and all(new_spacing <= maximum_values):
                        to_write.append(info)
                        to_write.append(sequence)
                        break
                info = handler.readline().strip()
                sequence = handler.readline().strip()

    for file in file_list:
        handler = open(file)
        info = handler.readline().strip()
        sequence = handler.readline().strip()
        while info != '':
            to_write.append(info)
            to_write.append(sequence)

    with open(out_file, 'wt') as out_file:
        for i in to_write:
            out_file.write(i+'\n')


def make_blast_db(path_to_blast_plus, in_path, out_file):
    subprocess.call([path_to_blast_plus,
                     '-in', in_path,
                     '-out', out_file,
                     '-dbtype', 'prot', ])



def make_blast(path_to_blast_plus, query_file, out_path, database):
    blastp_cline = NcbiblastpCommandline(cmd=path_to_blast_plus,
                                         query=query_file,
                                         db=database,
                                         evalue=10, outfmt=5,
                                         out=out_path,
                                         word_size=2,
                                         num_alignments=1000000000)
    stdout, stderr = blastp_cline()
    print(stdout, stderr)


def parse_blast_output(input_file, output_file):
    hit_dict = dict()
    result_handle = open(input_file)
    blast_records = list(NCBIXML.parse(result_handle))
    for blast_record in blast_records:
        print('next')
        alignments = list(blast_record.alignments)
        query_title = blast_record.query
        for alignment in alignments:
            alignment_haps = list(alignment.hsps)
            alignment_title = alignment.title
            for hsp in alignment_haps:
                if hsp.sbjct.count('C') != 6:
                    continue
                hsp.sbjct = ''.join(hsp.sbjct.split('-'))
                e_value = hsp.expect
                match = hsp.match
                if match.count('C') != 6:
                    continue
                if query_title in hit_dict:
                    hit_dict[query_title].append((alignment_title, e_value, hsp.sbjct, query_title), )
                else:
                    hit_dict[query_title] = [(alignment_title, e_value, hsp.sbjct, query_title)]
    for item in hit_dict:
        with open(output_file+item, 'wt') as out_file:
            for thing in hit_dict[item]:
                out_file.write(item +'\t'+ thing[0] +'\t'+ str(thing[1])+'\t'+ thing[2]+'\t'+ thing[3] +'\n')
