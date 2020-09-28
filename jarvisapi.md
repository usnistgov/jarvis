# API

## Overview

JARVIS-API is a Django-REST API framework for uploading and downloading data from [NIST-JARVIS](https://jarvis.nist.gov/).
A user can request an account through the web-interface to access a large and diverse amount of materials-data.
Click on the [Log In/Sign Up button](https://jarvis.nist.gov/login). The data is generated from several JARVIS resources 
such as JARVIS-DFT, JARVIS-FF, and JARVIS-ML etc. A user can also upload there own data using the API. 
Materials data can be uploaded in XML format or as binary blob files. The XML files are converted to HTML pages
using XSLT. XML upload requires a corresponding XSD file schema, but the blob uploads are possible without any fixed schema.
JARVIS-API can be used either using command line or using the web-interface. For beginners we recommend using the web-interface
for uploading and downloading data. 


## Methodology


- A user can send GET/POST requests to the JARVIS-API using the tools provided here: [REST-API](https://github.com/usnistgov/jarvis/blob/master/jarvis/db/restapi.py)

- A user should upload using the [Data Curation button](https://jarvis.nist.gov/curate/) on the web and select a particular template. 
If you need a custom made template, you can send an email to [Contact](https://jarvis.nist.gov/contact).

- There are several modules available right now to interact with the API such as: get XML data using ID, post text queries, upload files, delete your documents etc.

- After curating the data, you can find it using [Data Exploration](https://jarvis.nist.gov/explore/keyword/).




## References
1.      [JARVIS: An Integrated Infrastructure for Data-driven Materials Design, arXiv:2007.01831 (2020)](https://arxiv.org/abs/2007.01831).
