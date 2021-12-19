import os


def download_file(url, dir_out):

    """Download the data files from GitHub."""

    import requests
    from zipfile import ZipFile

    # Download the zip file

    response = requests.get(url, stream=True)
    file_out = "mute_data_files.zip"

    with open(file_out, "wb") as f:

        for line in response.iter_content(chunk_size=1024 * 1024):

            f.write(line)

    # Unzip the file

    with ZipFile(file_out, "r") as f:

        f.extractall(dir_out)

    os.chdir(dir_out)
    os.rename("data_20211219", "data")
    os.chdir("..")

    # Delete the zip file

    os.remove(file_out)


# Download the files from GitHub

if not os.path.isfile(
    os.path.join(os.path.dirname(__file__), "data", "data_20211219.txt")
):

    download_file(
        "https://github.com/wjwoodley/mute/releases/download/0.1.0/data_20211219.zip",
        os.path.dirname(__file__),
    )
