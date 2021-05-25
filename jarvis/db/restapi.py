"""Module to access or upload data in JARVIS-API."""
import requests
import pprint
import os
import glob
from jarvis.db.jsonutils import dumpjson
import string
from collections import OrderedDict
import json


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

    def get_global_workspace_id(self):
        """Get workspace ID."""
        turl = self.base_url + "/rest/workspace/"
        response = requests.get(
            turl, verify=False, auth=(self.username, self.password)
        )
        workspace_list = json.loads(response.text)
        for workspace in workspace_list:
            if workspace["title"] == "Global Public Workspace":
                global_workspace_id = workspace["id"]
        return global_workspace_id

    def upload_xml_file(
        self,
        filename="",
        workspace_id="5df7b6defb7e53000c4652aa",
        template_id="",
    ):
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
            "workspace": workspace_id,
            "xml_content": xml_content,
        }
        response = requests.post(
            turl, data=data, verify=False, auth=(self.username, self.password)
        )
        out = response.json()
        # pprint.pprint(out)
        response_code = response.status_code
        print("code:", response_code, requests.codes.ok)
        if response_code == requests.codes.created:
            print("status: uploaded.")
            upload_id = out["id"]
            print("upload_id", upload_id)
        else:
            response.raise_for_status()
            pprint.pprint(out)
            raise Exception(
                "A problem occurred when uploading the schema (Error ",
                response_code,
                ")",
            )

        return upload_id

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

    def get_all_xmls(self):
        """Get all XML data. Used for deleting."""
        return self.send_query(query=".xml")

    def delete_chunks(self, array=[]):
        """Delete a list of records."""
        # print ('array',array[0]['id'])
        for i in array:
            try:
                id = i["id"]
                self.delete_file_by_hash_id(id)
            except Exception:
                print("Failed deleting", i)
                print()
                pass

    def get_data_by_template_id(
        self, id="5f6162b4ece4b00034e4fe19", query=".xml", max_record=5
    ):
        """
        Get data by template IDs.

        For template IDs see:
        https://jarvis.nist.gov/rest/template-version-manager/global/
        and is_disabled=False
        """
        xml_upload_url = "/rest/data/query/keyword/"
        turl = self.base_url + xml_upload_url
        print("turl", turl)
        input = {"query": query}
        response = requests.post(
            turl, data=input, auth=(self.username, self.password)
        )
        out = response.json()
        params = {"page": 2}
        mem = []
        for i in out["results"]:
            if i["template"] == id:
                mem.append(i)
        if max_record is not None and len(mem) < max_record:
            while out["next"] is not None:
                response = requests.post(
                    turl,
                    data=input,
                    params=params,
                    auth=(self.username, self.password),
                )
                out = response.json()
                # mem.extend(out["results"])
                params["page"] += 1
                for i in out["results"]:
                    if i["template"] == id:
                        mem.append(i)
                        if max_record is not None:
                            if len(mem) == max_record:
                                break

        return mem

    def delete_all_records(self, query=".xml"):
        """Delete all XML files."""
        print("Deleteing all XMLs.")
        print("Hope you backed up the data.")
        xml_upload_url = "/rest/data/query/keyword/"
        turl = self.base_url + xml_upload_url
        print("turl", turl)
        input = {"query": query}
        response = requests.post(
            turl, data=input, auth=(self.username, self.password)
        )
        out = response.json()
        f = open("tmpout", "w")
        f.write(json.dumps(out))
        f.close()
        data = out["results"]
        self.delete_chunks(data)
        params = {"page": 2}

        while out["next"] is not None:
            response = requests.post(
                turl,
                data=input,
                params=params,
                auth=(self.username, self.password),
            )
            out = response.json()
            # data.extend(out["results"])
            data = out["results"]
            params["page"] += 1
            print("params", params)
            try:
                self.delete_chunks(data)
            except Exception:
                print("Failed for", data)
                pass
        return data

    def send_query(self, query="JVASP-1002.xml"):
        """Send query to API."""
        xml_upload_url = "/rest/data/query/keyword/"
        turl = self.base_url + xml_upload_url
        print("turl", turl)
        # input={"query":"JLMP-1254.xml"}
        # input={"query":"JVASP-1002.xml"}
        # input={"query":".xml"}
        input = {"query": query}
        response = requests.post(
            turl, data=input, auth=(self.username, self.password)
        )
        out = response.json()
        data = out["results"]
        params = {"page": 2}
        while out["next"] is not None:
            response = requests.post(
                turl,
                data=input,
                params=params,
                auth=(self.username, self.password),
            )
            out = response.json()
            data.extend(out["results"])
            params["page"] += 1
        return out

    def upload_blob(self, filepath=""):
        """Upload binary blob files."""
        turl = self.base_url + "/rest/blob/"

        oldid = os.path.basename(filepath)
        if all(c in string.hexdigits for c in oldid):
            newfile = oldid[25:]
        else:
            newfile = oldid

        files = {"blob": open(filepath, "rb")}

        data = {"filename": newfile}

        response = requests.post(
            turl, files=files, data=data, auth=(self.username, self.password)
        )
        if response.status_code == 201:
            print("Uploaded blob.")
            return response.json()["id"]
        else:
            response.raise_for_status()
            raise Exception("- error: a problem occurred when uploading blob")

    def get_all_blobs(self):
        """Download binary blob files."""
        turl = self.base_url + "/rest/blob/"
        response = requests.get(
            turl, verify=False, auth=(self.username, self.password)
        )
        # print ('response.json()',response.json())
        print("All blob doownload response.status_code", response.status_code)
        return response.json()

    def download_blob(self, save_file=True, bid="5f6a4bb7ece4b0002de4fb60"):
        """Download a binary blob file."""
        turl = self.base_url + "/rest/blob/" + bid + "/"
        # turl = self.base_url + "/rest/blob/download/" + bid + "/"
        response = requests.get(
            turl, verify=False, auth=(self.username, self.password)
        )
        resp = OrderedDict(response.json())
        data = requests.get(
            resp["handle"], verify=False, auth=(self.username, self.password)
        ).content
        filename = resp["filename"]
        if save_file:
            open(filename, "wb").write(data)
        print(len(data), len(filename))
        return filename, data

    def delete_blob(self, bid=""):
        """Delete a blob by id."""
        turl = self.base_url + "/rest/blob/" + str(bid) + "/"
        response = requests.delete(
            turl, verify=False, auth=(self.username, self.password)
        )
        response_code = response.status_code
        print("delete response_code", response_code)
        return response_code

    def delete_all_blobs(self):
        """Delete binary blob files."""
        data = self.get_all_blobs()
        for i in data:
            print("Deleting", i["id"])
            self.delete_blob(i["id"])

    def upload_jarvisff_xmls(
        self,
        files="/rk2/knc6/DB/XMLs/LAMMPS/*.xml",
        template_id="5f6162b4ece4b00034e4fe19",
        json_name="jarvisff_xmls.json",
    ):
        """Upload JARVIS-FF XML files."""
        mem = []
        for i in glob.glob(files):
            jid = i.split("/")[-1].split(".xml")[0]
            try:
                upload_id = self.upload_xml_file(
                    filename=i, template_id=template_id
                )
                info = {}
                info["jid"] = jid
                info["api_id"] = upload_id
                mem.append(info)

            except Exception:
                info = {}
                info["jid"] = jid
                info["api_id"] = "Failed"

                print("Failed for", i)
                pass
        if json_name is not None:
            dumpjson(filename=json_name, data=mem)
        return mem

    def upload_jarvisdft_xmls(
        self,
        files="/rk2/knc6/DB/XMLs/VASPDASK/*.xml",
        template_id="5f626925ece4b00035e5277f",
        json_name="jarvisdft_xmls.json",
    ):
        """Upload JARVIS-DFT XML files."""
        mem = []
        for i in glob.glob(files):
            jid = i.split("/")[-1].split(".xml")[0]
            print(jid)
            try:
                upload_id = self.upload_xml_file(
                    filename=i, template_id=template_id
                )
                print(jid)
                info = {}
                info["jid"] = jid
                info["api_id"] = upload_id
                mem.append(info)
            except Exception:
                info = {}
                info["jid"] = jid
                info["api_id"] = "Failed"

                print("Failed for", i)
                pass
        if json_name is not None:
            dumpjson(filename=json_name, data=mem)
        return mem


# TODO:
# workspace operations
# federated vs oai_pmh


"""
if __name__ == "__main__":

    a = Api()
    a = Api(base_url="http://test-jarvis.nist.gov")
    tid="5f624bc547d1ac011a224672"
    tid='5f6909278c6d9e011c8cc8c9'
    a.upload_xml_file(filename='JVASP-1004.xml',template_id=tid)

    a = Api()
    a = Api(base_url="https://jarvis.nist.gov")
    x = a.get_global_workspace_id()
    print("x=", x)

    a = Api()
    a = Api(base_url="https://jarvis.nist.gov")

    # a.upload_jarvisff_xmls()
    # a.upload_jarvisff_xmls()
    # a.upload_jarvisdft_xmls()
    # tid="5f626925ece4b00035e5277f"

    # a.upload_xml_file(filename='JVASP-1067.xml',template_id=tid)
    # a.upload_xml_file(filename='JVASP-664.xml',template_id=tid)
    # a.upload_xml_file(filename='JVASP-1002.xml',template_id=tid)

    # a.delete_all_records()

    # filepath="/rk2/knc6/DB/RAW_FILES/JARVIS-DFT-DFPT/JVASP-10088.zip",
    # filepath="/rk2/knc6/JARVIS-DFT/Elements-bulkk/
    # mp-149_bulk_PBEBO/MAIN-RELAX-bulk@mp-149/CHGCAR"
    # x = a.upload_blob(filepath)
    print("id", x)

    #a.upload_xml_file(filename='JVASP-1067.xml',template_id=tid)
    #a.upload_xml_file(filename='JVASP-664.xml',template_id=tid)
    #a.upload_xml_file(filename='JVASP-1002.xml',template_id=tid)

    # a.delete_all_records()

    #filepath="/rk2/knc6/DB/RAW_FILES/JARVIS-DFT-DFPT/JVASP-10088.zip",
    #filepath="/rk2/knc6/JARVIS-DFT/Elements-bulkk/
    #mp-149_bulk_PBEBO/MAIN-RELAX-bulk@mp-149/CHGCAR"
    #x = a.upload_blob(filepath)
    print ('id',x)
    a.delete_blob(x)
    # a.get_all_blobs()

    # a.upload_jarvisdft_xmls()
    # a.delete_all_blobs()
    # print(x)  # 5f61ba30ece4b00030e500c6
    # a.delete_file_by_hash_id('5f616c20ece4b00035e5135d')
    # upload_jarvisff_xmls()
    # out=a.send_query(query='JLMP-1254.xml')
    # print ('out',len(out))
    # data = a.get_data_by_template_id("5f6162b4ece4b00034e4fe19")
    # print(len(data))
    # a.upload_file(filename=x_file,template_id=id)
    # data=a.get_data_by_hash_id('5f6130e9ece4b0002fe4f8f1')
    # print (data)

    # templates=a.get_all_template_ids(title='jarvisff')
    # print ('templates',templates)
    # doc=a.get_data_by_hash_id('5f5d1578ece4b00035e4f9f5')
    # print (doc)

    import sys

    sys.exit()

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
    data = {"query": "JVASP-1002.xml", "title": "JVASP-1002.xml"}
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
