import logging
import os
from functools import partial
from multiprocessing import Pool
from pathlib import Path
from urllib.request import urlopen, Request

logger = logging.getLogger(__name__)


class ILBSMDatabaseDownloader:
    def __init__(self, search_url, download_urls, directory='out', create_separate_dirs=False, http_method='GET', sep=';'):
        self.directory = Path(directory)
        self.setup_download_dir()
        self.search_url = search_url
        self.http_method = http_method
        self.proteins = self.get_proteins(sep=sep)
        self.download_urls = download_urls
        self.create_separate_dirs = create_separate_dirs

    def get_proteins(self, sep=';'):
        req = Request(self.search_url, method=self.http_method)
        proteins = []
        with urlopen(req) as resp:
            data = resp.read().decode('utf-8')
            for line in data.splitlines():
                if len(line) > 4 and not line.startswith('#'):
                    cells = line.strip().split(sep)
                    proteins.append((cells[0], cells[1]))
        return proteins

    def download_link(self, download_url, protein):
        file_name = download_url.split('/')[-1]
        if not file_name.find('{0}') >= 0:
            file_name = '{0}_{1}_' + file_name
        directory_dir = os.path.join(self.directory, protein[0], protein[1]) if self.create_separate_dirs else self.directory
        if self.create_separate_dirs and not os.path.exists(directory_dir):
            os.makedirs(directory_dir, exist_ok=True)
        download_path = os.path.join(directory_dir, os.path.basename(file_name.format(protein[0], protein[1])))
        link = download_url.format(protein[0], protein[1])
        with urlopen(link) as image, open(download_path, 'wb') as f:
            f.write(image.read())
        logger.debug('Downloaded %s', link)

    def setup_download_dir(self):
        if not self.directory.exists():
            self.directory.mkdir()

    def get_all(self):
        for url in self.download_urls:
            download = partial(self.download_link, url)
            with Pool(8) as pl:
                pl.map(download, self.proteins)


def get_proteins_from_alphaknot(search_url):
    req = Request(search_url, method='GET')
    proteins = []
    with urlopen(req) as resp:
        data = resp.read().decode('utf-8')
        for line in data.splitlines():
            if len(line) > 4 and not line.startswith('#'):
                cells = line.strip().split(';')
                proteins.append((cells[0], cells[1]))
    return proteins

class AlphaKnotDatabaseDownloader(ILBSMDatabaseDownloader):
    def get_proteins(self, sep=None):
        cols = ['organism', 'uniprot', 'version', 'category']
        col_num_dict = {}
        req = Request(self.search_url, method=self.http_method)
        proteins = []
        with urlopen(req) as resp:
            data = resp.read().decode('utf-8')
            header = data.splitlines()[0]
            if header:
                header = header[1:] if header.startswith('#') else header
                if header.find('#') > 0:
                    header = header[:header.find('#')]
                cells = header.split()
                i=0
                for cell in cells:
                    if cell.strip().lower() in cols:
                        col_num_dict[cell.strip().lower()] = i
                    i += 1
            if 'uniprot' not in col_num_dict.keys():
                print(f'ERROR Missing "Uniprot" in the data returned from the Search URL {self.search_url}')
                return
            if 'category' not in col_num_dict.keys() and 'category' not in col_num_dict.keys():
                print(f'ERROR Missing "Category" or "Version" in the data returned from the Search URL {self.search_url}. At least one of those columns is necessary.')
                return
            if 'organism' not in col_num_dict.keys():
                print(f'WARNING Missing "Organism" in the data returned from the Search URL {self.search_url}. It may be necessary to retrieve data for AF1 models.')

            for line in data.splitlines():
                if len(line) > 4 and not line.startswith('#'):
                    cells = line.strip().split()
                    uniprot, entid = self.get_uniprot_and_model(cells[col_num_dict['uniprot']])
                    version = cells[col_num_dict['version']] if 'version' in col_num_dict.keys() else self.get_version_from_category(cells[col_num_dict['category']])
                    proteins.append((uniprot, entid, version, cells[col_num_dict['organism']]))
        return proteins

    def get_uniprot_and_model(self, uniid):
        uniprot, entid = uniid.split('-')
        return uniprot, entid.replace('F', '')

    def get_version_from_category(self, cat):
        if cat.startswith('AF'):
            return cat[2:3]
        elif cat.startswith('ESM'):
            return 'E' + cat[3:4]
        else:
            return None

    def download_link(self, download_url, protein):
        file_name = download_url.split('/')[-1]
        if not file_name.find('{0}') >= 0:
            file_name = '{0}_{1}_{2}' + file_name
        directory_dir = os.path.join(self.directory, protein[0], protein[1]) if self.create_separate_dirs else self.directory
        if self.create_separate_dirs and not os.path.exists(directory_dir):
            os.makedirs(directory_dir, exist_ok=True)
        download_path = os.path.join(directory_dir, os.path.basename(file_name.format(protein[0], protein[1], protein[2])))
        link = download_url.format(protein[0], protein[1], protein[2], protein[3])
        with urlopen(link) as image, open(download_path, 'wb') as f:
            f.write(image.read())
        logger.debug('Downloaded %s', link)
