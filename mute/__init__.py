import os
from shutil import copytree, rmtree

GitHub_data_file = "data_20220715"


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

    # Rename the data directory "data"
    # If it already exists, move the unzipped files into it and replace existing ones with the same name
    # If it does not, rename the unzipped directory "data"

    os.chdir(dir_out)

    if os.path.isdir("data"):

        copytree(GitHub_data_file, "data", dirs_exist_ok=True)
        rmtree(GitHub_data_file)

    else:

        os.rename(GitHub_data_file, "data")

    os.chdir("..")

    # Delete the zip file

    if os.path.isfile(file_out):
        os.remove(file_out)


# Download the files from GitHub to the directory MUTE is being run from

if not os.path.isfile(
    os.path.join(os.path.dirname(__file__), "data", f"{GitHub_data_file}.txt")
):

    download_file(
        f"https://github.com/wjwoodley/mute/releases/download/0.1.0/{GitHub_data_file}.zip",
        os.path.dirname(__file__),
    )
