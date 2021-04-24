import os
import glob
import gzip
from Bio import SeqIO
from ftplib import FTP
import datetime
import psutil
import time

date = str(datetime.datetime.now()).split(' ')[0]


def check_for_memory():
    while True:
        memory = psutil.virtual_memory().percent
        print(memory)
        if memory < 50:
            return
        print('Sleeping')
        time.sleep(300)


def login_and_cwd(path):
    my_ftp = FTP('ftp.ncbi.nlm.nih.gov')
    my_ftp.login()
    my_ftp.cwd(path)
    return my_ftp


def next_ftp_path(database, out_path, temp_path):
    check_for_memory()
    plant_list = login_and_cwd(database).nlst()
    plant_list = [i for i in plant_list if '.txt' not in i]
    already_done_list = [i.split(os.sep)[-1] for i in glob.glob(out_path+'*')] +\
                        [i.split(os.sep)[-1] for i in glob.glob(temp_path + '*')]

    plant_list = sorted([i for i in plant_list if i not in already_done_list])
    print(len(plant_list))
    if plant_list:
        return database+plant_list[0]


def download_item(ftp_path, temp_path, plant_name):
    available = login_and_cwd(path=ftp_path).nlst()
    if 'representative' in available:
        path = ftp_path+'/representative/'
    elif 'latest_assembly_versions' in available:
        path = ftp_path+'/latest_assembly_versions/'
    else:
        path = ftp_path+'/all_assembly_versions/'
    my_ftp = login_and_cwd(path=path)
    file = my_ftp.nlst()[0]
    if file == 'suppressed':
        os.mkdir(temp_path + os.sep + plant_name)
        return
    file_path = path + file + '/' + file + '_genomic.fna.gz'
    os.mkdir(temp_path+os.sep+plant_name)
    local_path = temp_path+os.sep+plant_name+os.sep+file
    my_ftp.cwd('/')
    with open(local_path, 'wb') as out_file:
        my_ftp.retrbinary('RETR ' + file_path, out_file.write)
    return local_path


def unzip_item(local_path):
    my_input = gzip.GzipFile(local_path, 'rb')
    s = my_input.read()
    my_input.close()
    output = open(local_path, 'wb')
    output.write(s)
    output.close()
    print("done")
    return str(os.path.getsize(local_path))


def translate(local_path, plant_name, out_path):
    cys_rich_sequence_counter = 0
    not_cys_rich_sequence_counter = 0
    my_seq = SeqIO.parse(local_path, 'fasta')
    record_number = 0
    os.mkdir(out_path + os.sep + plant_name)
    handler = open(out_path+os.sep+plant_name+os.sep+local_path.split('/')[-1], 'wt')
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
                        handler.write(
                            '>'+gen_bank_id+' '+plant_name+' '+numbers+' '+pattern+'\n'+str(pro)+'\n')
                        cys_rich_sequence_counter += 1
                    else:
                        not_cys_rich_sequence_counter += 1

    handler.close()
    os.remove(local_path)
    os.rmdir(os.sep.join(local_path.split(os.sep)[:-1]))
    return cys_rich_sequence_counter, not_cys_rich_sequence_counter


def update_log_file(local_path, out_path, plant_name, size, cys_rich, rest):
    with open(out_path+os.sep+plant_name+os.sep+plant_name+'.log', 'wt') as out_file:
        out_file.write('\t'.join([plant_name, str(size), str(cys_rich), str(rest), local_path.split('/')[-1]]))


def ftp_downloader(database, temp_path, out_path):
    """Downloads genomes or transcriptomes from NCBI's FTP-server. At the beginning the temporary path (temp_path) and
    the out_path are checked if the genome/transcriptome was already downloaded or if an error file was created. If yes
    the genome/transcriptome will be skipped.

    Parameters:
    database (str): path to the folder where the genetic data is located (e.g. '/genomes/genbank/plant/') on NCBI´s
    FTP-server ('ftp.ncbi.nlm.nih.gov').
    temp_path (str): path (starting from the root directory) where temporary files can be stored (e.g.
    '/home/user/Data/temporary/').

    Returns:
    In case a new entry was found in the specified folder of NCBI´s FTP server the file will be downloaded to the
    output path (out_path). In case of an error an empty file will be created in the temp_path"
    """
    ftp_path = next_ftp_path(database=database, out_path=out_path, temp_path=temp_path)
    while ftp_path:
        plant_name = ftp_path.split('/')[-1]
        print(plant_name)
        local_path = download_item(ftp_path=ftp_path, temp_path=temp_path, plant_name=plant_name)
        if local_path:
            size = unzip_item(local_path=local_path)
            cys_rich, rest = translate(local_path=local_path, plant_name=plant_name, out_path=out_path)
            update_log_file(local_path=local_path, out_path=out_path,
                            plant_name=plant_name, size=size, cys_rich=cys_rich, rest=rest)
        else:
            print('error')
        ftp_path = next_ftp_path(database=database, out_path=out_path, temp_path=temp_path)


def main():
    # ftp_downloader(database='genomes/genbank/vertebrate_mammalian/',
    #                temp_path='/home/b/PycharmProjects/bap_data/temporary/vertebrate_mammalian/',
    #                out_path='/home/b/PycharmProjects/bap_data/genbank/vertebrate_mammalian/')

    ftp_downloader(database='genomes/genbank/plant/',
                   temp_path=os.sep.join(['G:', 'bap_data', 'temporary', 'plant', '']),
                   out_path=os.sep.join(['G:', 'bap_data', 'genbank', 'plant', '']))

    # ftp_downloader(database='genomes/genbank/invertebrate/',
    #                temp_path='/home/b/PycharmProjects/bap_data/temporary/invertebrate/',
    #                out_path='/home/b/PycharmProjects/bap_data/genbank/invertebrate/')


if __name__ == '__main__':
    main()
