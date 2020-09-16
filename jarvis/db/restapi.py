"""Module to access or upload data in JARVIS-API."""
import requests
import pprint
import os
import glob
from lxml import etree


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

    def upload_xml_file(self, filename="", template_id=""):
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
        out = response.json()
        upload_id = out["id"]
        # pprint.pprint(out)
        print("upload_id", upload_id)
        response_code = response.status_code
        # print ('code:',response_code, requests.codes.ok)
        if response_code == requests.codes.created:
            print("status: uploaded.")
        else:
            response.raise_for_status()
            pprint.pprint(out)
            raise Exception(
                "- error: a problem occurred when uploading the schema (Error ",
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

    def delet_chunks(self, array=[]):
        """Delete a list of records."""
        for i in array:
            id = i["id"]
            self.delete_file_by_hash_id(id)

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
                data.extend(out["results"])
                params["page"] += 1
                for i in out["results"]:
                    if i["template"] == id:
                        mem.append(i)
                        if max_record is not None:
                            if len(mem) == max_record:
                                break

        return mem

    def delete_all_records(self):
        """Delete all XML files."""
        print("Deleteing all XMLs.")
        print("Hope you know what you are doing.")
        print("At least store the function return to JSON for backup.")
        xml_upload_url = "/rest/data/query/keyword/"
        turl = self.base_url + xml_upload_url
        print("turl", turl)
        input = {"query": query}
        response = requests.post(
            turl, data=input, auth=(self.username, self.password)
        )
        out = response.json()
        data = out["results"]
        self.delete_chunks(out)
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
            self.delete_chunks(out)

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

    def upload_jarvisff_xmls(
        self,
        files="/rk2/knc6/DB/XMLs/LAMMPS/*.xml",
        template_id="5f6162b4ece4b00034e4fe19",
        json_name="jarvisff_xmls.json",
    ):
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


# TODO:
# blob operations
# Adding certifiates
# workspace operations
# Add a url in views:JARVIS-DFT,FF etc.
# federated vs oai_pmh
# Make download function work
# upload_dataa()



def is_xml_valid(xsd="jarvisdft.xsd", xml="JVASP-1002.xml"):
    """Check if XML is valid."""
    xml_file = etree.parse(xml)
    xml_validator = etree.XMLSchema(file=xsd)
    is_valid = xml_validator.validate(xml_file)
    return is_valid


"""
if __name__ == "__main__":
    a = Api()
    a = Api(base_url="https://jarvis.nist.gov")
    a.upload_jarvisff_xmls()
    # a.upload_xml_file(filename='out1.xml',template_id='5f6162b4ece4b00034e4fe19')

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
