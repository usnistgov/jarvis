"""Module to access or upload data in JARVIS-API."""
import requests
import pprint
import os


class Api(object):
    """Class for handling functions for JARVIS-API."""

    def __init__(
        self,
        base_url="https://jarvis.nist.gov",
        user_info_file_path="/users/knc6/myinfo",
    ):
        """
        Intialize the class.

        Args:
            base_url: For JARVIS

            Request an account prior to access.
            Store your username and password locally
            in your a file. Format:
            ---------------------------------------
            username
            password
        """
        f = open(user_info_file_path, "r")
        lines = f.read().splitlines()
        f.close()
        self.base_url = base_url
        self.username = lines[0]
        self.password = lines[1]

    def get_data_by_hash_id(self, id=""):
        """Gat data file info for a hash_id."""
        template_upload_url = "/rest/data/" + str(id) + str("/")
        turl = self.base_url + template_upload_url + str("/?format=json")
        response = requests.get(
            turl, verify=False, auth=(self.username, self.password)
        )
        if response.ok:
            return response.json()
        else:
            raise Exception("Problem occurred while getting the file.")

    def send_keyword_query(
        self, data={"query": "JVASP-1002.xml", "title": "JVASP-1002.xml"}
    ):
        """Post general query to API."""
        template_upload_url = "/rest/data/query/keyword/"
        turl = self.base_url + template_upload_url
        response = requests.post(
            turl, verify=False, data=data, auth=(self.username, self.password)
        )
        if response.ok:
            return response.json()
        else:
            raise Exception("Problem occurred while getting the file.")
            raise Exception("Check username,password, ... details")

    def upload_file(self, filename="", template_id=""):
        """Upload file using path of the file and schema/template id."""
        print("status: uploading data file...")
        xml_file = open(filename, "rb")
        xml_content = xml_file.read()

        xml_upload_url = "/rest/data/"
        turl = self.base_url + xml_upload_url
        base = os.path.basename(filename)
        print("Filename:" + base)
        print("Template_id:" + template_id)
        data = {
            "title": base,
            "template": template_id,
            "xml_content": xml_content,
        }
        response = requests.post(
            turl, data=data, verify=False, auth=(self.username, self.password)
        )
        response_code = response.status_code
        print("upload response_code", response_code)
        return response_code

    def delete_file_by_hash_id(self, id=""):
        """Delete files for a hash id."""
        turl = self.base_url + "/rest/data/" + str(id) + "/"
        response = requests.delete(
            turl, verify=False, auth=(self.username, self.password)
        )
        response_code = response.status_code
        print("delete response_code", response_code)
        return response_code

    def get_all_template_ids(self, title=None, only_user=False):
        """Get templates for a user or global templates/schema/xsd."""
        # title: exmaple:'testxy'
        if only_user:
            turl = self.base_url + "/rest/template-version-manager/user/"
        else:
            turl = self.base_url + "/rest/template-version-manager/global/"
        params = {}
        if title is not None:
            params["title"] = title
        response = requests.get(
            turl, params=params, auth=(self.username, self.password)
        )
        print("response", pprint.pprint(response.json()))

        response_code = response.status_code
        if response_code == requests.codes.ok:
            print("- status: successful.")
            return response.json()
        else:
            response.raise_for_status()
            raise Exception(
                "- error: a problem occurred when getting the template"
            )


# TODO:
# blob operations
# Adding certifiates
# workspace operations
# Add a url in views:JARVIS-DFT,FF etc.
# federated vs oai_pmh
# Make download function work

"""
if __name__ == "__main__":
    a = Api()
    id = "5df7e8c7eaf3b300328be81b"
    x_file = "/users/knc6/Software/Devs/jarvis/jarvis/db/testxy.xml"
    t_id = "5f395d85ece4b00031fb5cfe"
    t_id = "5f395de8ece4b00028fb5d92"
    # templates=a.get_all_template_ids()
    # templates=a.get_all_template_ids(only_user=True,title='testxy')
    id = "5f396b0dece4b00036fb5dab"
    id = "5f396ae9ece4b00035fb5d42"
    # document=a.get_data_by_hash_id(id=id)
    # print('document',document)
    # a.upload_file(filename=x_file,template_id=t_id)
    data={"query": "JVASP-1002.xml", "title": "JVASP-1002.xml"}
    # x=a.send_keyword_query(data=data)
    id = "5f396ae9ece4b00035fb5d42"
    id = "5f396a44ece4b00035fb5d3e"
    # a.delete_file_by_hash_id(id=id)
    data = {"query": "testxy.xml", "title": "testxy.xml"}
    data = {"query": "testxy", "template": "5f395de8ece4b00028fb5d92"}
    data = {"query": "testxy"}
    result = a.send_keyword_query(data=data)
    print("result", result)
"""
