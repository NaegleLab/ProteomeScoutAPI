import tarfile
import requests
import os
import json

class FigshareInterface:
    """
    Class for interacting with FigShare API and downloading datasets.
    
    Parameters
    ----------
    article_id: int 
        The FigShare article ID. this can be found in the FigShare DOI link.
    """
    def __init__(self, article_id):
        self.article_id = article_id

    def get_available_versions(self):
        """
        Get version information for a FigShare dataset (article).

        Returns
        -------
            List of version info dicts containing api link and version numbers
        """
        url = f"https://api.figshare.com/v2/articles/{self.article_id}/versions"
        response = requests.get(url)
        if response.status_code != 200:
            raise Exception(f"Error fetching versions: {response.status_code}")
        
        return response.json()

    def get_latest_dataset_version(self):
        """
        Get the latest version number for a FigShare dataset (article).

        Returns
        -------
        int
            Latest version number of the dataset.
        """
        url = f"https://api.figshare.com/v2/articles/{self.article_id}/versions"
        response = requests.get(url)
        if response.status_code != 200:
            raise Exception(f"Error fetching versions: {response.status_code}")
        
        #iterate through versions to find the latest
        versions = response.json()
        latest_version = max(versions, key=lambda v: v['version'])
        return latest_version['version']

    def get_version_download_link(self, version_number):
        """
        Get download link for a specific version of a FigShare dataset.

        Parameters
        ----------
        str
            download link for the specified version.
        
        Raises
        ------
        Exception
            If there is an error fetching version files from FigShare
        """
        url = f"https://api.figshare.com/v2/articles/{self.article_id}/versions/{version_number}/files"
        response = requests.get(url)
        if response.status_code != 200:
            raise Exception(f"Error fetching version files: {response.status_code}")
        else:
            file_info = response.json()
            return file_info[0]['download_url']
        
    def get_version_meta(self, version_number):
        """
        Get metadata for a specific version of a FigShare dataset.
        
        Parameters
        ----------
        version_number: int
            The version number to get metadata for.

        Returns
        -------
        dict
            Metadata dictionary for the specified version, including date updated, file size, figshare URL, download link, and DOI.
        """
        url = f"https://api.figshare.com/v2/articles/{self.article_id}/versions/{version_number}"
        response = requests.get(url)
        meta = response.json()
        
        meta_output = {
            "version_number": meta['version'],
            "figshare_url": meta['figshare_url'],
            "download_link": meta['files'][0]['download_url'],
            "size_in_GB": meta['size'] / (1024 ** 3),
            "data_last_updated": meta['created_date'],
            'doi': meta['doi']
        }
        return meta_output  
    
    def download(self, version_number, output_path):
        """
        Download a specific version of a FigShare dataset and decompress it to the specified output path. Also saves metadata for version.
        
        Parameters
        ----------
        version_number: int 
            The version number to download.
        output_path: str
            The path to save the downloaded file.

        Postconditions
        --------------
        The dataset is downloaded, extracted to the output path in the "ProteomeScout_Dataset" directory, and metadata is saved to a metadata.json file in the same directory.
        """
        latest_version = self.get_latest_dataset_version()
        if isinstance(version_number, str):
            version_number= int(version_number)

        if version_number is None:
            print("Identifying latest version of ProteomeScout dataset on FigShare...")
            version_number = latest_version
        elif version_number <= 0:
            raise ValueError("Version number must be a positive integer")
        elif version_number > latest_version:
            raise ValueError(f"Requested version {version_number} exceeds latest available version {latest_version} and does not exist")
        else:
            print(f"Downloading ProteomeScout dataset version {version_number} from FigShare...")
            download_link = self.get_version_download_link(version_number)

            #download and write to file
            r = requests.get(download_link, allow_redirects=True)
            outputFile =  output_path + "/ProteomeScout_Dataset.tar.gz"
            open(outputFile, 'wb').write(r.content)

            #extract tar file
            t = tarfile.open(outputFile, 'r')
            print("Extracting %s" % outputFile)
            t.extractall(output_path)
            t.close()
            #remove tar file
            os.remove(outputFile)

            #extract dataset metadata and save to meta_file
            meta_file = os.path.join(output_path, "ProteomeScout_Dataset", "metadata.json")
            meta = self.get_version_meta(version_number)
            with open(meta_file, 'w') as f:
                json.dump(meta, f, indent=4)