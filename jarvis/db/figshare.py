"""
Downloads files from https://figshare.com/authors/Kamal_Choudhary/4445539
"""

import zipfile, os, requests
from jarvis.db.jsonutils import loadjson

def datasets(dataset=''):
    if dataset == 'dft_2d':
        url = "https://ndownloader.figshare.com/files/22471019"
        js_tag = 'jdft_2d-4-26-2020.json'
        print ('Downloading 2D dataset ...')
    elif dataset == 'dft_3d':
        url = "https://ndownloader.figshare.com/files/22471022"
        js_tag = 'jdft_3d-4-26-2020.json' 
        print ('Downloading 3D dataset ...')

    elif dataset == 'cfid_3d':
        url = "https://ndownloader.figshare.com/files/22470818"
        js_tag = 'jml_3d-4-26-2020.json' 
        print ('Downloading 3D CFID dataset ...')
    else:
        ValueError('Dataset doesnt exist',dataset)
    return url, js_tag

    

def data(dataset='dft_2d') :
    url, js_tag = datasets(dataset)
    path = str(os.path.join(os.path.dirname(__file__),js_tag ))
    if not  os.path.isfile(path):
        zfile = str(os.path.join(os.path.dirname(__file__), "tmp.zip"))
        r = requests.get(url)
        f = open(zfile, "wb")
        f.write(r.content)
        f.close()

        with zipfile.ZipFile(zfile, 'r') as zipObj:
            #zipObj.extract(path)
            zipObj.extractall(os.path.join(os.path.dirname(__file__)))
        os.remove(zfile)
    data = loadjson(path)
    return data




def get_ff_eneleast():
    jff1 = str(os.path.join(os.path.dirname(__file__), "jff1.json"))
    if not os.path.isfile(jff1):
        r = requests.get("https://ndownloader.figshare.com/files/10307139")
        f = open(jff1, "wb")
        f.write(r.content)
        f.close()
    data_ff1 = loadjson(jff1)
    return data_ff1



"""
if __name__ == "__main__":

    data_2d = data(dataset='dft_2d')
    print('2d',len(data_2d))
    data_3d = data(dataset='dft_3d')
    print('3d',len(data_3d))
    data_ml = data(dataset='cfid_3d')
    print('cfid3d',len(data_ml))
    data_ff = get_ff_eneleast()
    print ('ff',len(data_ff))
"""
