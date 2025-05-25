import os
from shutil import copytree, rmtree

mute_message = """*************************************************************************
*                                                                       *
*                ███████████████████████████████████████                *
*                ▓  ▓▓▓▓  ▓▓  ▓▓▓▓  ▓▓        ▓▓       ▓                *
*                ▓   ▓▓   ▓▓  ▓▓▓▓  ▓▓▓▓▓  ▓▓▓▓▓  ▓▓▓▓▓▓                *
*                ▒        ▒▒  ▒▒▒▒  ▒▒▒▒▒  ▒▒▒▒▒       ▒                *
*                ▒  ▒  ▒  ▒▒  ▒▒▒▒  ▒▒▒▒▒  ▒▒▒▒▒  ▒▒▒▒▒▒                *
*                ░  ░░░░  ░░░░    ░░░░░░░  ░░░░░       ░                *
*                ░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░                *
*                   https://github.com/wjwoodley/mute                   *
*                                                                       *
* Author:  William Woodley                                              *
* Version: 3.0.0                                                        *
*                                                                       *
* Please cite:                                                          *
*  - https://inspirehep.net/literature/1927720                          *
*  - https://inspirehep.net/literature/2799258                          *
*  - The models used for daemonflux, MCEq, PROPOSAL, and mountain maps  *
*                                                                       *
*************************************************************************"""

print(mute_message)

GitHub_data_file = "data_20250524"


def download_file(url, dir_out):
    """Download the data files from GitHub."""

    import requests
    from zipfile import ZipFile

    # Download the zip file

    response = requests.get(url, stream=True)
    file_out = os.path.join(dir_out, "mute_data_files.zip")

    with open(file_out, "wb") as f:

        for chunk in response.iter_content(chunk_size=1024 * 1024):

            f.write(chunk)

    # Unzip the file

    with ZipFile(file_out, "r") as f:

        f.extractall(dir_out)

    # If the "data" directory already exists, move the unzipped files into it and replace existing ones with the same name
    # If it does not, rename the unzipped directory "data"

    extracted_path = os.path.join(dir_out, GitHub_data_file)
    data_path = os.path.join(dir_out, "data")

    if os.path.isdir(data_path):

        copytree(extracted_path, data_path, dirs_exist_ok=True)
        rmtree(extracted_path)

    else:

        os.rename(extracted_path, data_path)

    # Delete the zip file

    if os.path.isfile(file_out):
        os.remove(file_out)


# Download the data files from GitHub to the directory MUTE is being run from if needed

_data_initialised = False


def initialise_data():

    global _data_initialised

    if _data_initialised:

        return

    else:

        if not os.path.isfile(
            os.path.join(os.path.dirname(__file__), "data", f"{GitHub_data_file}.txt")
        ):

            download_file(
                f"https://github.com/wjwoodley/mute/releases/download/0.1.0/{GitHub_data_file}.zip",
                os.path.dirname(__file__),
            )

        _data_initialised = True


initialise_data()
