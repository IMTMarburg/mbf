#!/usr/bin/python3
import sys
import os
import subprocess
import json

# This is the register with website code"""


def load_meta_data():
    try:
        with open("web/scb/metadata.json") as op:
            return json.load(op)
    except OSError:
        raise ValueError(
            "web/scb/metadata.json not found - make sure your run script calls mbf.publish.scb.prep_scb()"
        )


def _rsync_to_server(project_name):
    input_path = os.path.abspath(".") + "/"
    if "ANYSNAKE2_PROJECT_DIR" in os.environ:
        id_file = "/home/ffs/.ssh/id_rsa"
        input_path = os.environ["ANYSNAKE2_PROJECT_DIR"] + input_path[len("/project") :]
    else:
        id_file = "/.ffs_ssh/id_rsa"
    cmd = [
        "sudo",
        # '-u',
        # 'ffs',
        "rsync",
        input_path,
        "ffs@mbf.imt.uni-marburg.de:/mf/scb/%s/@@@chmod_after@@@chmod=o+rwX,g+rwX@@@chown=1000:2000"
        % (project_name),
        "--files-from=web/scb/rsync_list.txt",
        "-e",
        f"ssh -p 223 -i {id_file} -o StrictHostKeyChecking=no",
        "--rsync-path=rprsync",
        "-r",
        "-P",
        "-t",  # otherwise we retrigger sync because of the chmod after transfer...
        "-v",
    ]
    if "ANYSNAKE2_PROJECT_DIR" in os.environ:
        import sys
        import shlex

        print(
            "you have to execute this in your project path outside the container (sorry)"
        )
        print("")
        print(" ".join([shlex.quote(x) for x in cmd]))
        print("")
        print("press enter when you're done")
        sys.stdin.readline()
    else:
        subprocess.check_call(cmd)


def _register_with_server(accepted_server, path, revision):
    import requests

    auth = requests.auth.HTTPBasicAuth(
        os.environ["MBF_AUTH_USER"], os.environ["MBF_AUTH_PASSWORD"]
    )
    if "scb_server" in os.environ:
        top_level_url = os.environ["scb_server"]
    else:
        top_level_url = accepted_server
    url = top_level_url + "/register/%s?revision=%s" % (path, revision)
    print(url)
    req = requests.get(url, auth=auth)
    if req.status_code == 200:
        print("registered")
    else:
        print("error registring")
        print(req.text)


accepted_servers = {
    "scb": "http://mbf.imt.uni-marburg.de/scb",
    "scb_dev": "http://mbf.imt.uni-marburg.de/scb_dev",
    "localhost": "http://localhost:8080/scb",
}


def print_usage(msg=""):
    print("Usage:")
    print("scb_submit.py")
    if msg:
        print(msg)
    sys.exit()


def get_current_repo_revision():
    """Does not require auto commit"""
    x = subprocess.check_output(["hg", "log", "-r", "tip", "-q"]).decode("utf-8")
    return x[: x.find(":")]


def main():
    try:
        path = os.environ["ANYSNAKE_PROJECT_PATH"]
        project_name = path.split("/")[-1]
    except KeyError:
        try:
            path = os.environ["ANYSNAKE2_PROJECT_DIR"]
            project_name = path.split("/")[-1]
        except KeyError:
            print_usage("Must be run from inside anysnake")
    print("submitting %s to scb..." % project_name)
    print("now rsyncing")
    print("sudo password is test123 in anysnake container!")
    _rsync_to_server(project_name)

    server = os.environ.get("SCB_SERVER_URL", "http://mbf.imt.uni-marburg.de/scb")
    print("calling webserver")
    _register_with_server(server, project_name, get_current_repo_revision())


if __name__ == "__main__":
    main()
