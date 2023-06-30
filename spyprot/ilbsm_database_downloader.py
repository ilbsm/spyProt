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
    def get_proteins(self, sep=';'):
        cols = ['organism', 'uniprot', 'version', 'category']
        col_num_dict = {}
        req = Request(self.search_url, method=self.http_method)
        proteins = []
        with urlopen(req) as resp:
            data = resp.read().decode('utf-8')
            header = data.splitlines()[0]
            if header:
                cells = header.strip().split(sep)
                i=0
                for cell in cells:
                    if cell.lower() in cols:
                        col_num_dict[cell] = i
                    i += 1
            for line in data.splitlines():
                if len(line) > 4 and not line.startswith('#'):
                    cells = line.strip().split(sep)
                    proteins.append((cells[0], cells[1]))
        return proteins






