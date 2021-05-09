import re
import glob
import os


class Hit:
    def __init__(self, re_obj, seq, e_value):
        self.re_obj = re_obj
        self.seq = seq
        self.e_value = e_value


def yield_next_file(file_path):
    file_list = glob.glob(file_path + os.sep + '*')
    for file in file_list:
        my_file = glob.glob(file + os.sep + '*')
        my_file = [i for i in my_file if '.log' not in i]
        yield my_file[0]


def get_complete_contigs(in_path, out_path):
    to_write = []
    re_dict = dict()
    with open(in_path) as file:
        for line in file:
            description = line.strip().split('\t')[0].split(' ')
            e_value = float(line.strip().split('\t')[1])
            info_description = '_'.join(description[3:])
            name = description[2].split(' ')[1]
            sequence = line.strip().split('\t')[-1]
            if name in re_dict:
                re_dict[name].append(Hit(re.compile(info_description), sequence, e_value))
            else:
                re_dict[name] = [Hit(re.compile(info_description), sequence, e_value), ]
    already_done = 0
    for my_file in yield_next_file(file_path=r'C:\Users\b\PycharmProjects\bap_data\test_data\genbank\plant'):
        for name in re_dict:
            if name in my_file:
                already_done += 1
                print(my_file)
                handler = open(my_file)
                info = handler.readline().strip()
                sequence = handler.readline().strip()
                while info and sequence != '':
                    poss_matches = [i for i in [j for j in re_dict[name]] if i.re_obj.search(info)]
                    if poss_matches:
                        match = poss_matches[0]
                        to_write.append((info, sequence, match.e_value, match.seq, ))
                    info = handler.readline().strip()
                    sequence = handler.readline().strip()

    to_write.sort(key= lambda x: x[2])
    with open(out_path, 'wt') as file:
        for j in to_write:
            file.write(j[0]+'\t'+j[1]+'\t'+str(j[2])+'\t'+j[3]+'\n')


def main():
    get_complete_contigs(in_path=r'C:\Users\b\PycharmProjects\bap_data\anaylsis_data\Cyclotides_1\5_formated_blast_output\my_output.txt',
                         out_path=r'C:\Users\b\PycharmProjects\bap_data\anaylsis_data\Cyclotides_1\6_get_complete_contigs\complete_contigs.txt')


if __name__ == '__main__':
    main()
