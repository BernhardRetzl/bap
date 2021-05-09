from Bio.Blast import NCBIXML


class Alignment:
    def __init__(self, seq_number, plant_name, contig_number, pattern, e_value, hsp_sbjct, hsp_counter):
        self.seq_number = seq_number
        self.plant_name = plant_name
        self.contig_number = contig_number
        self.pattern = pattern
        self.e_value = e_value
        self.hsp_sbjct = hsp_sbjct
        self.hsp_counter = hsp_counter


def parse_blast_output(in_path, out_path):
    hit_dict = dict()
    for blast_record in NCBIXML.parse(open(in_path)):
        for alignment in blast_record.alignments:
            hsp_counter = 0
            alignment_title = alignment.title
            for hsp in alignment.hsps:
                if hsp.sbjct.count('C') != 6:
                    continue
                hsp_counter += 1
                hsp.sbjct = ''.join(hsp.sbjct.split('-'))
                e_value = hsp.expect
                junk, seq_number, plant_name, contig_number, pattern = alignment_title.split(' ')
                if hsp_counter > 1:
                    alignment_title += '_hsp_' + str(hsp_counter)
                if alignment_title in hit_dict:
                    if hit_dict[alignment_title][0] > e_value:
                        hit_dict[alignment_title] = (e_value, hsp.sbjct)
                    else:
                        continue
                else:
                    hit_dict[alignment_title] = Alignment(seq_number=1)


    my_tuple = list(zip(hit_dict.keys(), hit_dict.values()))
    my_tuple.sort(key=lambda x: x[1][0])
    with open(out_path, 'wt') as file:
        for i in my_tuple:
            file.write(i[0]+'\t'+str(i[1][0])+'\t'+i[1][1]+'\n')


def main():
    parse_blast_output(in_path=r'C:\Users\b\PycharmProjects\bap_data\anaylsis_data\Cyclotides_1\4_blast_output\my_blast_output',
                       out_path=r'C:\Users\b\PycharmProjects\bap_data\anaylsis_data\Cyclotides_1\5_formated_blast_output\my_output.txt')


if __name__ == '__main__':
    main()
