import hashlib
import urllib.request
import os
import tempfile
import zipfile


class DownloadData:
    """Utility class for downloading test data to the website"""
    FILES = {
        "errors.bin":
        "47a1ef9cb699113ff34f74ac761561ba19f2cb97bfab1896fd28a4c0adce3716",
        "gen_signal_1d.bin":
        "7d18e6dca7e11c7c6be79c5b4c5d6f6cb7e6d6177c43f2e70d4efb143682c037",
        "gen_signal_2d_rectangle.bin":
        "9f8c8578cc5ff40991d4ea0ec2e92858b84c7d52e921345d1428e515a7e5b03c",
    }
    URL = ("https://github.com/CNES/swot_simulator/releases/download"
           "/0.0.0/TestDatasets.zip")

    def __init__(self):
        self.prefix = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   "data")
        os.makedirs(self.prefix, exist_ok=True)
        while not self.check():
            archive = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                   "TestDatasets.zip")
            try:
                with open(archive, mode="wb") as stream:
                    self.download(stream)
                self.extract(archive)
            finally:
                os.unlink(archive)

    def extract(self, name):
        """Extraction of the downloaded archive."""
        zip_file = zipfile.ZipFile(name)
        for name in zip_file.namelist():
            info = zip_file.getinfo(name)
            deflated = os.path.join(self.prefix, name)
            with open(deflated, "wb") as stream:
                stream.write(zip_file.read(info))

    def check(self):
        """Verify data integrity"""
        for item in self.FILES:
            path = os.path.join(self.prefix, item)
            if not os.path.exists(path) or \
                    self.sha256sum(path) != self.FILES[item]:
                return False
        return True

    @classmethod
    def download(cls, stream):
        """Download data from bitbucket"""
        response = urllib.request.urlopen(cls.URL)

        while True:
            data = response.read(65536)
            if not data:
                break
            stream.write(data)

    @staticmethod
    def sha256sum(path):
        """
        Computes the SHA256 hash for a file
        """
        sha256 = hashlib.sha256()
        with open(path, 'rb') as stream:
            for block in iter(lambda: stream.read(65536), b''):
                sha256.update(block)
        return sha256.hexdigest()


DownloadData()