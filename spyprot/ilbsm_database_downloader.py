import logging
import os
from functools import partial
from multiprocessing import Pool
from urllib2 import Request, urlopen

logger = logging.getLogger(__name__)


def get_proteins(search_url):
    req = Request(search_url)
    proteins = []
    resp = urlopen(req)
    data = resp.read().decode('utf-8')
    for line in data.splitlines():
        if len(line) > 4 and not line.startswith('#'):
            cells = line.strip().split(';')
            proteins.append((cells[0], cells[1]))
    return proteins


class ILBSMDatabaseDownloader:
    def __init__(self, search_url, download_urls, directory='out', create_separate_dirs=False):
        self.directory = directory
        self.setup_download_dir()
        self.proteins = get_proteins(search_url)
        self.download_urls = download_urls
        self.create_separate_dirs = create_separate_dirs

    def download_link(self, download_url, protein):
        file_name = download_url.split('/')[-1]
        if not file_name.find('{0}') >= 0:
            file_name = '{0}_{1}_' + file_name
        directory_dir = os.path.join(self.directory, protein[0], protein[1]) if self.create_separate_dirs else self.directory
        if self.create_separate_dirs and not os.path.exists(directory_dir):
            os.makedirs(directory_dir)
        download_path = os.path.join(directory_dir, os.path.basename(file_name.format(protein[0], protein[1])))
        link = download_url.format(protein[0], protein[1])
        image = urlopen(link)
        with open(download_path, 'wb') as f:
            f.write(image.read())
        logger.debug('Downloaded %s', link)

    def setup_download_dir(self):
        if not os.path.exists(self.directory):
            os.makedirs(self.directory)

    def get_all(self):
        for url in self.download_urls:
            download = partial(self.download_link, url)
            for p in self.proteins:
                self.download_link(url, p)

