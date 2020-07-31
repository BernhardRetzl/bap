import os
import glob
import gzip
import shutil
from Bio import SeqIO
from ftplib import FTP, error_perm


def ftp_downloader(database, temp_path, out_path, database_name):
    """Downloads genomes or transcriptomes from NCBI's FTP-server"""
    def download_item(item):
        try:
            my_ftp = FTP('ftp.ncbi.nlm.nih.gov')
            my_ftp.login()
            path = database + item
            my_ftp.cwd(path)
            available = my_ftp.nlst()
            if 'representative' in available:
                path = database+item+'/representative/'
            elif 'latest_assembly_versions' in available:
                path = database+item+'/latest_assembly_versions/'
            else:
                path = database+item+'/all_assembly_versions/'
            my_ftp.cwd(path)
            plant = my_ftp.nlst()[0]
            file_path = path + plant + '/' + plant + '_genomic.fna.gz'
            with open(temp_path + item + '#' + plant, 'wb') as file:
                my_ftp.retrbinary('RETR ' + file_path, file.write)
            return temp_path + item + '#' + plant
        except error_perm:
            plant = my_ftp.nlst()[0]
            print('A download error occurred')
            with open(item + '#' + plant + '_error', 'wt'):
                pass
            return

    def unzip_item(item):
        with gzip.open(item, 'rb') as in_file:
            with open(item+'.temp', 'wb') as out_file:
                shutil.copyfileobj(in_file, out_file)
        os.remove(item)
        os.rename(item+'.temp', item)
        return str(os.path.getsize(item))

    def translate(item):
        cys_rich_sequence_counter = 0
        not_cys_rich_sequence_counter = 0
        to_write = []
        my_seq = SeqIO.parse(item, 'fasta')
        own_identifier = item.split('#')[0].split('/')[-1]
        record_number = 0
        for record in my_seq:
            record_number += 1
            gen_bank_id = record.description.split(' ')[0]
            frame_number = 0
            for nuc in [record.seq, record.seq.reverse_complement()]:
                for frame in range(3):
                    frame_number += 1
                    length = 3 * ((len(record) - frame) // 3)
                    orf_number = 0
                    for pro in nuc[frame: frame + length].translate().split('*'):
                        orf_number += 1
                        pattern = [len(i) for i in pro.split('C')][1:-1]
                        if len(pattern) > 4:
                            numbers = '_'.join([str(record_number), str(frame_number), str(orf_number)])
                            pattern = '-' + '-'.join([str(i) for i in pattern]) + '-'
                            to_write.append(
                                '>'+gen_bank_id+' '+own_identifier+' '+numbers+' '+pattern+'\n'+str(pro)+'\n')
                            cys_rich_sequence_counter += 1
                        else:
                            not_cys_rich_sequence_counter += 1
        with open(out_path+item.split('/')[-1], 'wt') as out_file:
            for i in to_write:
                out_file.write(i)
        os.remove(item)
        return cys_rich_sequence_counter, not_cys_rich_sequence_counter

    def write_statistics(sample_name, size, cys_rich_sequences, rest):
        with open('/home/b/Blast_project/Data/_1/statistics/'+database_name, 'at') as out_file:
            out_file.write(date+'\t'+sample_name+'\t'+str(size)+'\t'+database_name +
                           '\t'+str(cys_rich_sequences)+'\t'+str(rest)+'\n')

    def next_item():
        my_ftp = FTP('ftp.ncbi.nlm.nih.gov')
        my_ftp.login()
        my_ftp.cwd(database)
        plant_list = my_ftp.nlst()
        plant_list = [i for i in plant_list if '.txt' not in i]
        already_done_list = glob.glob(out_path+'/' + '*') + glob.glob(temp_path + '*')
        already_done_list = [i.split('/')[-1].split('#')[0] for i in already_done_list]
        plant_list = [i for i in plant_list if i not in already_done_list]
        plant_list.sort()
        print(len(plant_list))
        if plant_list:
            return plant_list[0]
        else:
            return

    item = next_item()
    while item:
        print(item)
        my_ftp = FTP('ftp.ncbi.nlm.nih.gov')
        my_ftp.login()
        my_ftp.cwd(database)
        item = download_item(item)
        if item:
            size = unzip_item(item)
            print(size)
            cys_rich_sequences, rest = translate(item)
            write_statistics(item, size, cys_rich_sequences, rest)
        item = next_item()
