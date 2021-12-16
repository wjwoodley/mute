import os


def download_file(url):

    """Download the data files from GitHub."""

    import requests
    from zipfile import ZipFile

    # Download the zip file

    response = requests.get(url, stream=True)
    file_out = "mute_data_files.zip"
    dir_out = os.path.dirname(__file__)

    with open(file_out, "wb") as f:

        for line in response.iter_content(chunk_size=1024 * 1024):

            f.write(line)

    # Unzip the file

    with ZipFile(file_out, "r") as f:

        f.extractall(dir_out)

    # Delete the zip file

    os.remove(file_out)


# Download the files from GitHub

if not os.path.isfile(os.path.join(os.path.dirname(__file__), "data", "data_20211216.txt")):
    pass
    
#     download_file("https://github.com/wjwoodley/mute/releases/download/1.0.0/data_20211216.zip")
