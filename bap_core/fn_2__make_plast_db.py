import subprocess


def make_blast_db(path_to_blast_execuatable, input_path, output_path):
    subprocess.call([path_to_blast_execuatable,
                     '-in', input_path,
                     '-out', output_path,
                     '-title', 'Peptides',
                     '-dbtype', 'prot', ])


def main():
    make_blast_db(path_to_blast_execuatable=r'C:\Program Files\NCBI\blast-2.11.0+\bin\makeblastdb.exe',
                  input_path=r'C:\Users\b\PycharmProjects\bap_data\anaylsis_data\Cyclotides_1\2_seq_for_blast_db\seq_for_blast_db.txt',
                  output_path=r'C:\Users\b\PycharmProjects\bap_data\anaylsis_data\Cyclotides_1\3_blast_db\my_blast_db')


if __name__ == '__main__':
    main()
