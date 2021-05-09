from Bio.Blast.Applications import NcbiblastpCommandline


def make_blast():
    blastp_cline = NcbiblastpCommandline(cmd=r'C:\Program Files\NCBI\blast-2.11.0+\bin\blastp.exe',
                                         query=r'C:\Users\b\PycharmProjects\bap_data\anaylsis_data\Cyclotides_1\1_queries\queries.txt',
                                         db=r'C:\Users\b\PycharmProjects\bap_data\anaylsis_data\Cyclotides_1\3_blast_db\my_blast_db',
                                         evalue=10, outfmt=5,
                                         out=r'C:\Users\b\PycharmProjects\bap_data\anaylsis_data\Cyclotides_1\4_blast_output\my_blast_output',
                                         num_alignments=50)
    stdout, stderr = blastp_cline()


def main():
    make_blast()


if __name__ == '__main__':
    main()
