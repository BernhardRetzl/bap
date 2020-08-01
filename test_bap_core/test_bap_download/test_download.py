import bap_core.bap_download.download


def test_next_ftp_path_1():

    actual = bap_core.bap_download.download.next_ftp_path(database='genomes/genbank/plant/', out_path='')
    assert actual.startswith('genomes/genbank/plant/A')

def test_2():
    from ftplib import FTP
    def login_and_cwd(path):
        my_ftp = FTP('ftp.ncbi.nlm.nih.gov')
        my_ftp.login()
        my_ftp.cwd(path)
        return my_ftp

    plant_list = login_and_cwd('genomes/genbank/plant/').nlst()
    plant_list = [i for i in plant_list if '.txt' not in i]
    already_done_list = list()
    plant_list = sorted([i for i in plant_list if i not in already_done_list])
    assert len(plant_list) >= 706

