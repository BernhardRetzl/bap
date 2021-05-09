import os
import glob
import numpy as np


def yield_next_file(file_path):
    file_list = glob.glob(file_path + os.sep + '*')
    for file in file_list:
        my_file = glob.glob(file + os.sep + '*')
        my_file = [i for i in my_file if '.log' not in i]
        yield my_file[0]


def get_sequences_for_blast_db(file_path, minimum_values, maximum_values, out_path):
    to_write = list()
    handler = open(file_path)
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
    handler.close()

    with open(out_path, 'wt') as out_file:
        for thing in to_write:
            out_file.write(thing + '\n')


def get_recommended_sequences(path_to_queries):
    with open(path_to_queries) as in_file:
        final_dict = dict()
        for line in in_file:
            if line.startswith('>'):
                continue
            elif line.count('C') == 6:
                for number, part in enumerate(line.split('C')):
                    if number in final_dict:
                        final_dict[number].append(len(part))
                    else:
                        final_dict[number] = [len(part),]
        minimum_list = list()
        maximum_list = list()
        print(final_dict)
        for item in final_dict:
            if item == 0 or item == 6:
                continue
            else:
                print('Loop : '+str(item)+'\n'+'The minimum value is: '+str(min(final_dict[item]))+
                      '\n'+'The maximum value is: '+str(max(final_dict[item])))
                minimum_list.append(min(final_dict[item]))
                maximum_list.append(max(final_dict[item]))


def main():
    get_recommended_sequences(path_to_queries=r'C:\Users\b\PycharmProjects\bap_data\anaylsis_data\Cyclotides_1\1_queries\queries.txt')
    new_minimum_values = np.array([3, 4, 4, 1, 4])
    new_maximum_values = np.array([3, 5, 7, 1, 7])
    for file in yield_next_file(file_path=r'C:\Users\b\PycharmProjects\bap_data\test_data\genbank\plant'):
        print(file)
        get_sequences_for_blast_db(file_path=file, minimum_values=new_minimum_values,
                                   maximum_values = new_maximum_values,
                                   out_path=r'C:\Users\b\PycharmProjects\bap_data\anaylsis_data\Cyclotides_1\2_seq_for_blast_db\seq_for_blast_db.txt')


if __name__ == '__main__':
    main()
