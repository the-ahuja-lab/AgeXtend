# Docker container

This is the repo of the official [Docker image](https://hub.docker.com/r/sanjayk741/agextend) for **AgeXtend**.

# Table of Contents  

1. [Introduction](#docker-container)  
2. [Installation Instructions](#installation-instructions)  
   - [Windows](#windows)  
   - [MacOS](#macos)  
   - [Linux](#linux)  
3. [Pulling and Running AgeXtend Docker Image](#pulling-and-running-agextend-docker-image)  
4. [Using AgeXtend](#using-agextend)  
   - [STEP-I: License Key Application](#step-i-license-key-application)  
   - [STEP-II: Import and Run AgeXtend](#step-ii-import-and-run-agextend)  
     - [Command Line Mode](#command-line-mode)  
     - [Jupyter Notebook Mode](#jupyter-notebook-mode)  
5. [Work Packages](#work-packages)  
   - [Predictor](#predictor)  
     - [Geroprediction and Explainability Modules](#geroprediction-and-explainability-modules)  
     - [Sub Explainability and Target Modules](#sub-explainability-and-target-modules)  
     - [Efficient Bulk Prediction / Custom AgeXtend Database](#efficient-bulk-prediction--custom-agextend-database)  
   - [Browser](#browser)  
6. [Pre-Complied AgeXtend Prediction Databases](#pre-complied-agextend-prediction-databases)  
7. [Input/Output Instructions](#inputoutput-instructions)  
   - [Input](#input)  
   - [Output](#output)  

**AgeXtend** is a multimodal, bioactivity-based, and mechanism-backed Explainable AI Predictor for geroprotectors. AgeXtend supports 2 distinct work packages: Predictor and Browser.
The first work package, **Predictor**, involves a bioactivity-based classifier (**Geroprediction** module) for predicting geroprotective compounds, coupled with an **Explainability** module providing mechanistic insights into biological processes behind the predictions, a **Toxicity** module to evaluate their toxicity in the biological context, and finally, a **Target** module that suggests protein targets of the putative geroprotective compounds.
The second work package, **Browser**, allows users to explore through the 20 pre-screened databases specifically curated for these functionalities.

Instructions for installing and running docker on any PC can be found [here](https://docs.docker.com/engine/install/) 
1. [Windows](https://docs.docker.com/desktop/install/windows-install/)
2. [MacOS](https://docs.docker.com/desktop/install/mac-install/)
3. [Linux](https://docs.docker.com/desktop/install/linux-install/)

**The Docker image for AgeXtend is currently only available for Linux systems. Please stay tuned for Windows and MacOS compatible versions.**

Pull the **AgeXtend** image from Docker Hub by running the following command in your terminal:
```
$ docker pull sanjayk741/agextend
```
Verify that the new image has been created using the **docker images** command.
```
$ docker images
```
To access the terminal of a Docker image, you can use the **docker run** command with the **-it** option.
```
$ docker run -it <image-name> /bin/bash
```
Replace **<image-name\>** with the name or ID of the Docker image of **AgeXtend**.

Find the ID of the currently running container for **input** and **output**.
```
$ docker ps -a
```
To start the container again and access its terminal. 
```
$ docker start <container-ID>
$ docker exec -it <container-ID> bash
```

## Using **AgeXtend**

**AgeXtend is free for academic institutions, however, for commercial utilization a commercial license key is required. Academic users may apply for a valid "License Key" [here](https://forms.gle/y1sCpSGEAML8XWGGA). Commerical users may request the license key from Dr. Gaurav Ahuja (gaurav.ahuja@iiitd.ac.in).**


### **STEP-I License Key Application**


Run the following command to get the hostname (container ID) and the IP address of your running container.

```
cat /etc/hostname
hostname -I
```

### **STEP-II Import and Run AgeXtend**
#### EITHER
Run AgeXtend by python **command line** mode directly in the terminal.
```
$ python
>>> from AgeXtend import Predictor
```
#### OR 

Run AgeXtend in **Jupyter Notebook**. 
Launch Jupyter Notebook by running the following command:
```
$ jupyter-notebook --ip=<ip-address> --allow-root
```
Here, replace **<ip-address\>** with the IP address that you got in **STEP-I**

This should start the Jupyter Notebook server and display a URL in the terminal that you can use to access the Jupyter Notebook interface.

## Work Packages
AgeXtend supports 2 distinct work packages:<br/>
1. Predictor
2. Browser

### Predictor
Activate AgeXtend license
```
>>> Predictor.license('license key') #Example: Predictor.license('KKKVFZ41111WF6RTQ')
```
#### Geroprediction and Explainability Modules ####

Predicts the anti-aging potential for the input SMILES:<br/>
```
>>> from AgeXtend import Predictor
```
Prepare a list of canonical SMILES (Open Babel generated) strings
```
>>> smiles = ['ClCC=C', 'C=CCOC(=O)CC(C)C'] 
```
Create an AgeXtend type object to featurize the query compounds (i.e. list of canonical SMILES)
```
>>> agex_obj = Predictor.featurize(smiles)
```
Use the AgeXtend object for predictions
```
>>> output = Predictor.predict(agex_obj)
```
Get the list of resulting data frames that are part of the rest of the modules
```
>>> output.keys()
dict_keys(['Anti_Aging_Prediction', 'Explainability_Status', 'Explainability_Probabilities'])
```
**Note:** Explainability_Probabilities will be empty by default unless selected otherwise by supplying additional arguments.

Get the result of a specific output data frame
```
>>> output['Explainability_Status']
```
#### Sub Explainability and Target Module ####

##### Additional arguments:
**AgeXtend** also supports the following modules along with the Prediction and Explainability modules, which can be chosen using the boolean option (True)

| Parameter Name | Description | Type | Default value | **Output (If True)** |
| -------- | -------- | -------- | -------- | -------- |
| probs | Probabilities of Explainability and/or Toxicity Module | boolean  | False | Explainability_Probabilities / Explainability_Toxicity_Probabilities |
| HC | Run Toxicity Module (Health Check) | boolean  | False | Explainability_Toxicity_Status |
| TS | Sub Explainability Level Tanimoto Similarity Test | boolean  | False | Explainability_response |
| BDL | BindingDB Target Information and Druggability (Lipinski Rule) | boolean  | False | Druggability_and_Potential_Targets |


**Example**
```
>>> output = Predictor.predict(agex_obj, TS=True)
>>> output.keys()
dict_keys(['Anti_Aging_Prediction', 'Explainability_Status', 'Explainability_Probabilities', 'Explainability_response'])
>>> output = Predictor.predict(agex_obj, HC=True, BDL=True)
>>> output.keys()
dict_keys(['Anti_Aging_Prediction', 'Explainability_Toxicity_Status', 'Explainability_Toxicity_Probabilities', 'Druggability_and_Potential_Targets'])
```

####  Efficient Bulk Prediction / Custom AgeXtend Database
Predictor module also provides functionalities for efficient bulk data prediction
```
>>> Predictor.bulk_predict(input=smiles_list) 
```
Or user can also use pre-calculated AgeXtend type object as input
```
>>> agex_obj = Predictor.featurize(smiles_list)
>>> Predictor.bulk_predict(input=agex_obj)
```

##### Additional arguments:

| Parameter Name | Description | Type | Default value |
| -------- | -------- | -------- | -------- |
| outfolder | Output Database type folder name | string | AgeXtendDB |
| chunksize | Number of predictions input per iteration for faster job completion | Integer | 10000 |

**Note**
The output folder of Bulk prediction function can be used as a Custom database for input in AgeXtend Browser module

### Browser

To explore the AgeXtend pre-complied predictions of various databases:<br/>
```
>>> from AgeXtend import Browser
```
Use Open Babel to Generate Canonical SMILES format of the query compound as input
```
>>> Browser.search(query='OC(=O)CCCc1c[nH]c2c1cccc2', output='/path/to/output/folder/')
```

Unzip the **AgeXtendBrowserOut.zip** file to visualize/print the generated report (HTML format)

**Note:** the **report file** (AgeXtend_BrowserOut.html) must be in the same folder with the **images** folder


#### Additional arguments:
**AgeXtend** also supports the use of locally complied Predictor module outputs (Folder containing CSV format outputs)

| Parameter Name | Description | Default Database |
| -------- | -------- | -------- |
| path | Path to the AgeXtend pre-complied prediction database | gutMGene |

**Example**
```
>>> Browser.search(path='/path/to/Database/Folder/', query='OC(=O)CCCc1c[nH]c2c1cccc2', output='/path/to/output/folder/')
```

### Pre-complied AgeXtend Prediction Databases
| DB Name | DB version | FTP size | FTP link |
| -------- | -------- | -------- | -------- |
| FooDB | Pre-release 1.0 | 63M | [http://agextend.ahujalab.iiitd.edu.in:8080/FOODB/](url) |
| HMDB | v5.0 | 231M | [http://agextend.ahujalab.iiitd.edu.in:8080/HMDB/](url) |
| IMPPAT | v2.0 | 4K | [http://agextend.ahujalab.iiitd.edu.in:8080/IMPPAT/](url) |
| AfroDb | - | 756K | [http://agextend.ahujalab.iiitd.edu.in:8080/AfroDB/](url) |
| AgingAtlas | v1.0 | 748K | [http://agextend.ahujalab.iiitd.edu.in:8080/AgingAtlas/](url) |
| ChemBridge | - | 788K | [http://agextend.ahujalab.iiitd.edu.in:8080/Chembridge/](url) |
| ChemDiv BBlocks | - | 39M | [http://agextend.ahujalab.iiitd.edu.in:8080/ChemdivBBlocks/](url) |
| CMNPD | v1.0 | 20M | [http://agextend.ahujalab.iiitd.edu.in:8080/CMNPD/](url) |
| DDPEDB | - | 4K | [http://agextend.ahujalab.iiitd.edu.in:8080/DDPD/](url) |
| ECMDB | v2.0 | 4K | [http://agextend.ahujalab.iiitd.edu.in:8080/ECMDB/](url) |
| RepoHub | release-3/24/2020 | 3.9M | [http://agextend.ahujalab.iiitd.edu.in:8080/RepoHub/](url) |



## Input/Output Instructions :-

### Input 
Find the ID of the currently running container, just executed using the **docker ps -a** command.
```
$ docker ps -a
```
To write a file to the container, use the **docker cp** command to copy the file from the host to the container.
```
$ docker cp file container_id:/root/CDir/
```
This command will copy the **file** file (pdb/tsv) from the host's current directory to the **AgeXtend** container with ID **container_id** at the CDir/ directory inside the container.

### Output
Find the ID of the currently running container, just executed using the **docker ps -a** command.
```
$ docker ps -a
```
To write a file from the container, use the **docker cp** command to copy the file from the container to the host.
```
$ docker cp container_id:/root/CDir/file-name .
```
This command will copy the **file** file (csv/pdf)from the **AgeXtend** container with ID **container_id** under the CDir/ directory inside the container to the host's current directory.
