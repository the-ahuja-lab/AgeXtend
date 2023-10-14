import os
import numpy as np
from PIL import Image
from rdkit import Chem
import matplotlib.pyplot as plt
import csv
import pandas as pd
import zipfile
import shutil
from importlib import resources
import io
import pkg_resources
import socket


print('---- Please read the documentation carefully for information about required input format (OpenBabel generated Canonical SMILES) and different functionalities of the AgeXtend Package. ----')

def license(key='',package_path=''):
        d={'A':'C','B':'D','C':'E','D':'F','E':'G','F':'H','G':'I','H':'J','I':'K','J':'L','K':'M','L':'N','M':'O','N':'P','O':'Q','P':'R','Q':'S','R':'T','S':'U','T':'V','U':'W','V':'X','W':'Y','X':'Z','Y':'A','Z':'B'}
        if len(key)==0:
                with resources.open_binary('AgeXtend','license.lic') as fp:
                        lcdata = fp.read()
                lcdata=io.BytesIO(lcdata)
                host_name = socket.gethostname()
                clean_nm=''.join([x for x in host_name if x.isalnum()])
                clean_key=[]
                for i in clean_nm:
                        if i.isalpha():
                                clean_key.append(d[i.upper()])
                        if i.isnumeric():
                                clean_key.append(str(int(i)+2))
                ob_key=''.join(clean_key)
                if lcdata.getvalue().decode() == 'MASTERBAIISHERE':
                        return 'GO'
                if lcdata.getvalue().decode() == ob_key:
                        return 'GO'
                else:
                        print('AgeXtend unlicensed: Please use license key for activation')
                        return 'STOP'
        else:
                host_name = socket.gethostname()
                clean_nm=''.join([x for x in host_name if x.isalnum()])
                clean_key=[]
                for i in clean_nm:
                        if i.isalpha():
                                clean_key.append(d[i.upper()])
                        if i.isnumeric():
                                clean_key.append(str(int(i)+2))
                ob_key=''.join(clean_key)
                if ob_key == key:
                        if len(package_path)==0:
                                package_path = pkg_resources.get_distribution('AgeXtend').location
                                write_byte = io.BytesIO(key.encode('ascii'))
                                with open(package_path+"/AgeXtend/license.lic", "wb") as f:
                                        f.write(write_byte.getbuffer())
                                print('License installed successfuly :)')
                        else:
                                if os.path.exists(package_path):
                                        write_byte = io.BytesIO(key.encode('ascii'))
                                        with open(package_path+"/license.lic", "wb") as f:
                                                f.write(write_byte.getbuffer())
                                        print('License installed successfuly :)')
                                else:
                                        print('Package path couild not be found !! please provide path along the licens :(')
                else:
                        print('Invalid license key provided !! AgeXtend activation failed :(')




'''
def predictor(path='',query='Oc1ccc(cc1)C1CC(=O)c2c(O1)cc(cc2O)O'):
  x=0
  if len(path)==0:
    with open('canage_predictions.csv')as fin:
      for line in fin:
        x+=1
        if x==1:
          continue
        else:
          dat=line.split(',')
          sm=dat[1]
          if sm == query:
            ft=dat[-1]
            return dat[2:-1],format(float(ft), ".2f")
  else:
    with open(path+'/canage_predictions.csv')as fin:
      for line in fin:
        x+=1
        if x==1:
          continue
        else:
          dat=line.split(',')
          sm=dat[1]
          if sm == query:
            ft=dat[-1]
            return dat[2:-1],format(float(ft), ".2f")
'''
def canpred(path):
    with open(path+'/file1.pkl', 'rb') as file:
        dt=pickle.load(file)
    for i in range(2,4):
        with open(path+'/file'+str(i)+'.pkl', 'rb') as file:
            dt1=pickle.load(file)
        dt=dt.merge(dt1, on='smiles', how='inner')
    dt.drop(columns=['Anti_Aging_Status'],inplace=True)
    dt.to_csv(path+'/canage_probabilities.csv') #.rename(columns={'Anti_Aging_Prob_y':'Anti_Aging_Prob'})

def probabilities(path='',query='Oc1ccc(cc1)C1CC(=O)c2c(O1)cc(cc2O)O'):
  x=0
  if len(path)==0:
    with resources.open_text('AgeXtend','canage_probabilities.csv') as fp:
      F_names = fp.readlines()
      for line in F_names:
        line=line.rstrip()
        x+=1
        if x==1:
          hd=['Anti Aging','     AIC','CS','DNS','EA','GI','LP','MD','SCE','TA','AMES','BBB','CYP1A2','CYP2C19','CYP2C9','CYP2D6','CYP3A4','DILI','Hepato','hERG','HLM','MMP','P-gp Inhibitor','P-gp Substrate']
          continue
        else:
          dat=line.split(',')
          sm=dat[1]
          if sm == query:
            ft=[]
            for i in dat[2:]:
              ft.append(float(format(float(i), ".4f")))
            return hd,ft,"Yes"
    return hd,[],"No"
  else:
    if os.path.exists(path+'/canage_probabilities.csv'):
      pass
    else:
      canpred(path)
    with open(path+'/canage_probabilities.csv')as fin:
      for line in fin:
        line=line.rstrip()
        x+=1
        if x==1:
          hd=['Anti Aging','     AIC','CS','DNS','EA','GI','LP','MD','SCE','TA','AMES','BBB','CYP1A2','CYP2C19','CYP2C9','CYP2D6','CYP3A4','DILI','Hepato','hERG','HLM','MMP','P-gp Inhibitor','P-gp Substrate']
          continue
        else:
          dat=line.split(',')
          sm=dat[1]
          if sm == query:
            ft=[]
            for i in dat[2:]:
              ft.append(float(format(float(i), ".4f")))
            return hd,ft,"Yes"
    return hd,[],"No"

def hallmark_tanimoto(path='',query='Oc1ccc(cc1)C1CC(=O)c2c(O1)cc(cc2O)O'):
  x=0
  if len(path)==0:
    with resources.open_text('AgeXtend','hallmark_tanimoto_test.csv') as fp:
      F_names = fp.readlines()
      for line in F_names:
        line=line.rstrip()
        x+=1
        if x==1:
          hd=['AIC_EA_SIRT_activator','AIC_NFkB_inhibitor','AIC_NLRP3_inhibitor','CS_HSP90_inhibitor','CS_Kinase_inhibitor','CS_Senolytic','DNS_Caloric_restriction_mimetic','DNS_mTOR_inhibitor','EA_Antagonist_for_HDAC_inhibitor','GI_CD38_inhibitor','GI_farnesyltranserase_inhibitor','LP_HSP70_activator','LP_HSR_activator','LP_Macroautophagy_inducer','MD_PINK1/Parkin_mitophagy_pathway_activator','MD_mitochondrial_dysfunction_protector/VDAC1_inhibitor','MD_mitochondrial_permeability_transition_pore_inhibitor','MD_mitophagy_inducer','SCE_JAK1_inhibitor','SCE_JAK3_inhibitor','SCE_Rho-associated_kinase_inhibitor_(anti-apoptosis)','SCE_Stem_cell_proliferators','SCE_prevents_stem_cell_differentiation','SCE_promote_stem_cell_reprogramming','SCE_self_renewal_of_stemcells','SCE_self_renewal_of_stemcells/Stem_cell_proliferation','TA_Activators_of_telomere_damage_signalling','TA_Telomerase_activator','TA_anti_telomere_shortening/preserving_telomere_length','TA_blocks_telomere_shortening','TA_derepresses_hTERT_expression']
          continue
        else:
          dat=line.split(',')
          sm=dat[1]
          if sm == query:
            ft_sm=[]
            ft_sc=[]
            for i in range(2,len(dat)):
              if i%2!=0:
                ft_sm.append(dat[i].rstrip())
              else:
                ft_sc.append(format(float(dat[i]), ".3f"))
            return hd,ft_sm,ft_sc
  else:
    with open(path+'/hallmark_tanimoto_test.csv')as fin:
      for line in fin:
        line=line.rstrip()
        x+=1
        if x==1:
          hd=['AIC_EA_SIRT_activator','AIC_NFkB_inhibitor','AIC_NLRP3_inhibitor','CS_HSP90_inhibitor','CS_Kinase_inhibitor','CS_Senolytic','DNS_Caloric_restriction_mimetic','DNS_mTOR_inhibitor','EA_Antagonist_for_HDAC_inhibitor','GI_CD38_inhibitor','GI_farnesyltranserase_inhibitor','LP_HSP70_activator','LP_HSR_activator','LP_Macroautophagy_inducer','MD_PINK1/Parkin_mitophagy_pathway_activator','MD_mitochondrial_dysfunction_protector/VDAC1_inhibitor','MD_mitochondrial_permeability_transition_pore_inhibitor','MD_mitophagy_inducer','SCE_JAK1_inhibitor','SCE_JAK3_inhibitor','SCE_Rho-associated_kinase_inhibitor_(anti-apoptosis)','SCE_Stem_cell_proliferators','SCE_prevents_stem_cell_differentiation','SCE_promote_stem_cell_reprogramming','SCE_self_renewal_of_stemcells','SCE_self_renewal_of_stemcells/Stem_cell_proliferation','TA_Activators_of_telomere_damage_signalling','TA_Telomerase_activator','TA_anti_telomere_shortening/preserving_telomere_length','TA_blocks_telomere_shortening','TA_derepresses_hTERT_expression']
          continue
        else:
          dat=line.split(',')
          sm=dat[1]
          if sm == query:
            ft_sm=[]
            ft_sc=[]
            for i in range(2,len(dat)):
              if i%2!=0:
                ft_sm.append(dat[i].rstrip())
              else:
                ft_sc.append(format(float(dat[i]), ".3f"))
            return hd,ft_sm,ft_sc


def bindingDB_ids(path='',query='Oc1ccc(cc1)C1CC(=O)c2c(O1)cc(cc2O)O'):
    with resources.open_binary('AgeXtend','BINDINGDB_db_MFPs_lb.pkl') as fp:
        db_file = fp.read()
    dbMFP = pd.read_pickle(io.BytesIO(db_file))
    x=0
    ft=[]
    for filn in os.listdir(path):
        if filn.startswith('bindingdb_info_lb'):
            with open(path+'/'+filn)as fin:
                for line in fin:
                    line=line.rstrip()
                    x+=1
                    if x==1:
                        continue
                    else:
                        reader=csv.reader(line.splitlines())
                        dat=list(reader)[0]
                        if dat[0]==query:
                            ct=dat[1].count(':')
                            sstng=dat[1]
                            up=1
                            end=len(sstng)
                            for t in range(ct):
                                sstng=sstng[up:end]
                                g=1
                                while g<sstng.index(':'):
                                    try:
                                        int(sstng[g:sstng.index(':')-1])
                                    except:
                                        g+=1
                                        continue
                                    else:
                                        mpid=int(sstng[g:sstng.index(':')-1])
                                        break
                                if t+1==ct:
                                    scor=format(float(sstng[sstng.index(':')+1:sstng.index(')')-1]),".3f")
                                else:
                                    scor=format(float(sstng[sstng.index(':')+1:sstng.index(',')-1]),".3f")
                                    up=sstng.index(',')+2
                                print(mpid)
                                ft.append(dbMFP['BindingDB Ligand Name'].to_list()[mpid])
                                ft.append(scor)
                                ft.append(dbMFP['Target Name'].to_list()[mpid])
                                ft.append(dbMFP['Target Source Organism'].to_list()[mpid])
                            return ft


def lipinsky(path='',query='Oc1ccc(cc1)C1CC(=O)c2c(O1)cc(cc2O)O'):
  x=0
  if len(path)==0:
    with resources.open_text('AgeXtend','lipinsky_info.csv') as fp:
      F_names = fp.readlines()
      for line in F_names:
        line=line.rstrip()
        x+=1
        if x==1:
          continue
        else:
          dat=line.split(',')
          sm=dat[1]
          if sm == query:
            ft=[]
            i=2
            while i < len(dat):
              if dat[i].startswith('\"Over 5'):
                ft.append(str(int(dat[i+1].split(' ')[-1].split("\"")[0])))
                i+=2
                continue
              if dat[i].startswith('\"Over 10'):
                ft.append(str(int(dat[i+1].split(' ')[-1].split("\"")[0])))
                i+=2
                continue
              if dat[i].startswith('Found'):
                ft.append(str(int(dat[i].split(' ')[1])))
              if ':' in dat[i]:
                ft.append(format(float(dat[i].split(' ')[-1]), ".8f"))
              if dat[i].endswith('over 500'):
                ft.append(format(float(dat[i+1].split(' ')[-1].split("\"")[0]), ".8f"))
                i+=2
                continue 
              if dat[i].endswith('over 5'):
                ft.append(format(float(dat[i+1].split(' ')[-1].split("\"")[0]), ".8f"))                       
                i+=2
                continue 
              if dat[i].startswith('True'):
                ft.append(str(1))
              if dat[i].startswith('False'):
                ft.append(str(0))
              i=i+1
            return ft
  else:
    for filn in os.listdir(path):
      if filn.startswith('lipinsky_info'):
        with open(path+'/'+filn)as fin:
          for line in fin:
            line=line.rstrip()
            x+=1
            if x==1:
              continue
            else:
              dat=line.split(',')
              sm=dat[1]
              if sm == query:
                ft=[]
                i=2
                while i < len(dat):
                  if dat[i].startswith('\"Over 5'):
                    ft.append(str(int(dat[i+1].split(' ')[-1].split("\"")[0])))
                    i+=2
                    continue
                  if dat[i].startswith('\"Over 10'):
                    ft.append(str(int(dat[i+1].split(' ')[-1].split("\"")[0])))
                    i+=2
                    continue
                  if dat[i].startswith('Found'):
                    ft.append(str(int(dat[i].split(' ')[1])))
                  if ':' in dat[i]:
                    ft.append(format(float(dat[i].split(' ')[-1]), ".8f"))
                  if dat[i].endswith('over 500'):
                    ft.append(format(float(dat[i+1].split(' ')[-1].split("\"")[0]), ".8f"))
                    i+=2
                    continue 
                  if dat[i].endswith('over 5'):
                    ft.append(format(float(dat[i+1].split(' ')[-1].split("\"")[0]), ".8f"))                       
                    i+=2
                    continue 
                  if dat[i].startswith('True'):
                    ft.append(str(1))
                  if dat[i].startswith('False'):
                    ft.append(str(0))
                  i=i+1
                return ft

def zip_folder(folder_path, output_path):
    with zipfile.ZipFile(output_path, 'w', zipfile.ZIP_DEFLATED) as zipf:
        for root, dirs, files in os.walk(folder_path):
            for file in files:
                file_path = os.path.join(root, file)
                zipf.write(file_path, os.path.relpath(file_path, folder_path))

    # Remove the original folder
    shutil.rmtree(folder_path)


def bindingDB(path='',query='Oc1ccc(cc1)C1CC(=O)c2c(O1)cc(cc2O)O'):
  x=0
  sm=''
  ft=[]
  if len(path)==0:
    with resources.open_text('AgeXtend','bindingdb_info.csv') as fp:
      F_names = fp.readlines()
      for line in F_names:
        line=line.rstrip()
        x+=1
        if x==1:
          continue
        else:
          reader = csv.reader(line.splitlines())
          dat = list(reader)[0]
          if dat[4]==query:
            if len(dat[1])>0:
              ft.append(dat[1])
            else:
              ft.append('NA')
            if len(dat[-1])>0:
              ft.append(format(float(dat[-1]), ".3f"))
            else:
              ft.append("0.0")
            if len(dat[2])>0:
              ft.append(dat[2])
            else:
              ft.append('NA')
            if len(dat[3])>0:
              ft.append(dat[3]) 
            else:
              ft.append('NA')
      return ft
  else:
    for filn in os.listdir(path):
      if filn.startswith('bindingdb_info_lb'):
        return bindingDB_ids(path,query)
        continue
      if filn.startswith('bindingdb_info'):
        with open(path+'/'+filn)as fin:
          for line in fin:
            line=line.rstrip()
            x+=1
            if x==1:
              continue
            else:
              reader = csv.reader(line.splitlines())
              dat = list(reader)[0]
              if dat[4]==query:
                if len(dat[1])>0:
                  ft.append(dat[1])
                else:
                  ft.append('NA')
                if len(dat[-1])>0:
                  ft.append(format(float(dat[-1]), ".3f"))
                else:
                  ft.append("0.0")
                if len(dat[2])>0:
                  ft.append(dat[2])
                else:
                  ft.append('NA')
                if len(dat[3])>0:
                  ft.append(dat[3]) 
                else:
                  ft.append('NA')
    return ft

def search(path='',query='Oc1ccc(cc1)C1CC(=O)c2c(O1)cc(cc2O)O',output=''):
  if len(path)>0:
    st=license()
    if st=='STOP':
      return 0
    else:
      for dirpath, dirnames, files in os.walk(path):
        for file_name in files:
          if file_name.endswith('pkl'):
            with open(dirpath+'/'+file_name, "rb") as file:
              data_frame = pd.read_pickle(file)
              data_frame.to_csv(dirpath+'/'+file_name.split('.')[0]+'.csv')
  if len(output)==0:
    if not os.path.exists('AgeXtendBrowserOut'):
      os.mkdir('AgeXtendBrowserOut')
      outpath='AgeXtendBrowserOut/'
    else:
      outpath='AgeXtendBrowserOut/'
  else:
    if not os.path.exists(output+'/AgeXtendBrowserOut'):
      os.mkdir(output+'/AgeXtendBrowserOut')
      outpath=output+'/AgeXtendBrowserOut/'
    else:
      outpath=output+'/AgeXtendBrowserOut/'
  if not os.path.exists(outpath+'images/'):
    os.mkdir(outpath+'images')
  
  opath=path
  if 'Chunk_0' in os.listdir(opath):
    for c in range(len(os.listdir(path))):
      if opath.endswith('/'):
        path=opath+'Chunk_'+str(c)
      else:
        path=opath+'/Chunk_'+str(c)
      lbs,val,fnd=probabilities(path=path,query=query)
      if fnd=="Yes":
        break
    if fnd=="No":
      return "Sorry, Query compound prediction data not found\nPlease proceed with Predictor module"
  else:
    lbs,val,fnd=probabilities(path=path,query=query)
    if fnd=="No":
      return "Sorry, Query compound prediction data not found\nPlease proceed with Predictor module"
  sizes = [val[0]*100,(1-val[0])*100]
  labels = ['Geroprotector', 'Neutral']
  colors = ['#A42c41', '#666666']
  fig, ax = plt.subplots()
  wedges= ax.pie(sizes, labels=labels, colors=colors, wedgeprops=dict(width=0.5), autopct='%1.2f%%', startangle=90)
  ax.set(aspect="equal", title='Prediction Module')
  plt.savefig(outpath+'images/image3.png', dpi=600)

  
  pb1=[num * 100 for num in val[1:-1]]
  pb0=[(1-num) * 100 for num in val[1:-1]]
  lb=lbs[1:]
  
  categories = lb
  categories.reverse()
  values = pb1
  values.reverse()
  colors = ['#FE9079','#F2C5B7','#FED766','#9CE09D','#18A68D','#3490DC','#AF81C9','#F66D9B','#A42C41','#C5577C','#EC8A7F','#E0B58A','#F0E1CE','#F5DF4E','#FDAB55','#B65930','#9B8A52','#A1DBA9','#82B27D','#00A070','#76C1BD','#9CB7D4','#0172B6']
  colors.reverse()
  bar_width = 0.8
  y = np.arange(len(categories[13:]))
  fig, ax = plt.subplots()
  rects = ax.barh(y, values[13:], bar_width, color=colors[13:], edgecolor='black')
  ax.set_xlabel('Percentage')
  ax.set_xlim(0,100)
  ax.set_title('Mechanistic Level Explainability')
  ax.set_yticks(y)
  ax.set_yticklabels(categories[13:])
  fig.set_size_inches(10, 8)
  plt.savefig(outpath+'images/image7.png', dpi=600)

  bar_width = 0.8
  y = np.arange(len(categories[:13]))
  fig, ax = plt.subplots()
  rects = ax.barh(y, values[:13], bar_width, color=colors[:13], edgecolor='black')
  ax.set_xlabel('Percentage')
  ax.set_xlim(0,100)
  ax.set_title('Toxicity Module')
  ax.set_yticks(y)
  ax.set_yticklabels(categories[:13])
  fig.set_size_inches(10, 8)
  plt.savefig(outpath+'images/image8.png', dpi=600)

  with resources.open_binary('AgeXtend',"image6.png") as fp:
    F_names = fp.read()
  image = Image.open(io.BytesIO(F_names))
  image.save(outpath+'images/image6.png')
  image.close()

  with resources.open_binary('AgeXtend',"image5.png") as fp:
    F_names = fp.read()
  image = Image.open(io.BytesIO(F_names))
  image.save(outpath+'images/image5.png')
  image.close()

  with resources.open_binary('AgeXtend',"image2.png") as fp:
    F_names = fp.read()
  image = Image.open(io.BytesIO(F_names))
  image.save(outpath+'images/image2.png')
  image.close()

  with resources.open_binary('AgeXtend',"image1.png") as fp:
    F_names = fp.read()
  image = Image.open(io.BytesIO(F_names))
  image.save(outpath+'images/image1.png')
  image.close()

  with resources.open_binary('AgeXtend',"image4.png") as fp:
    F_names = fp.read()
  image = Image.open(io.BytesIO(F_names))
  image.save(outpath+'images/image4.png')
  image.close()

  with resources.open_text('AgeXtend',"AgeXtend_BrowserOut.html") as fp:
    F_names = fp.readlines()
  with open(outpath+'AgeXtend_BrowserOut.html','w')as fout:
    with resources.open_text('AgeXtend',"AgeXtend_BrowserOut.html") as fp:
      F_names = fp.readlines()
      for line in F_names:
        line=line.replace('Qrys',query)
        try:
          mol = Chem.MolFromSmiles(query)
        except:
          pass
        else:
          line=line.replace('image6','image2')
        
        if len(path)>0:
          line=line.replace('gutMGene',path.split('/')[-1].rstrip())

        line=line.replace('CONCENTRATION',str(val[-1]))
        lipout=lipinsky(path=path,query=query)
        if int(lipout[0])>5:
          line=line.split('<span>LHD</span>')
          line[1]='<span style="overflow: hidden; display: inline-block; margin: 0.00px 0.00px; border: 0.00px solid #000000; transform: rotate(0.00rad) translateZ(0px); -webkit-transform: rotate(0.00rad) translateZ(0px); width: 20.00px; height: 25.43px;"><img alt="" src="images/image1.png" style="width: 24.12px; height: 25.43px; margin-left: -1.76px; margin-top: -1.51px; transform: rotate(0.00rad) translateZ(0px); -webkit-transform: rotate(0.00rad) translateZ(0px);" title=""></span>'
          line=''.join(line)
        else:
          line=line.split('<span>LHD</span>')
          line=''.join(line)
        if int(lipout[1])>10:
          line=line.split('<span>LHA</span>')
          line[1]='<span style="overflow: hidden; display: inline-block; margin: 0.00px 0.00px; border: 0.00px solid #000000; transform: rotate(0.00rad) translateZ(0px); -webkit-transform: rotate(0.00rad) translateZ(0px); width: 20.00px; height: 25.43px;"><img alt="" src="images/image1.png" style="width: 24.12px; height: 25.43px; margin-left: -1.76px; margin-top: -1.51px; transform: rotate(0.00rad) translateZ(0px); -webkit-transform: rotate(0.00rad) translateZ(0px);" title=""></span>'
          line=''.join(line)
        else:
          line=line.split('<span>LHA</span>')
          line=''.join(line)
        if float(lipout[2])<500:
          line=line.split('<span class="c17">LMW</span>')
          line[1]='<span style="overflow: hidden; display: inline-block; margin: 0.00px 0.00px; border: 0.00px solid #000000; transform: rotate(0.00rad) translateZ(0px); -webkit-transform: rotate(0.00rad) translateZ(0px); width: 20.00px; height: 20.00px;"><img alt="" src="images/image4.png" style="width: 20.00px; height: 20.00px; margin-left: 0.00px; margin-top: 0.00px; transform: rotate(0.00rad) translateZ(0px); -webkit-transform: rotate(0.00rad) translateZ(0px);" title=""></span>'
          line=''.join(line)
        else:
          line=line.split('<span class="c17">LMW</span>')
          line=''.join(line)
        if float(lipout[3])<5:
          line=line.split('<span class="c17">LLP</span>')
          line[1]='<span style="overflow: hidden; display: inline-block; margin: 0.00px 0.00px; border: 0.00px solid #000000; transform: rotate(0.00rad) translateZ(0px); -webkit-transform: rotate(0.00rad) translateZ(0px); width: 20.00px; height: 20.00px;"><img alt="" src="images/image4.png" style="width: 20.00px; height: 20.00px; margin-left: 0.00px; margin-top: 0.00px; transform: rotate(0.00rad) translateZ(0px); -webkit-transform: rotate(0.00rad) translateZ(0px);" title=""></span>'
          line=''.join(line)
        else:
          line=line.split('<span class="c17">LLP</span>')
          line=''.join(line)
        if lipout[4]=='1':
          line=line.replace('LSTPass','Pass')
        else:
          line=line.replace('LSTPass','Fail')
        
        bdout=bindingDB(path=path,query=query)
        if len(bdout)< 12:
          bdout.extend(['NA']*(12-len(bdout)))
        i=0
        j=1
        while i< len(bdout):
          line=line.replace('Rank'+str(j)+'_Lname',bdout[i])
          line=line.replace('Rank'+str(j)+'_TS',bdout[i+1])
          line=line.replace('Rank'+str(j)+'_Tname',bdout[i+2])
          line=line.replace('Rank'+str(j)+'_Sname',bdout[i+3])
          i+=4
          j+=1

        v=["A","B","C","D","E","F","G","H","6","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z","1","2","3","4","5"]
        hmout=hallmark_tanimoto(path=path,query=query)
        for t in range(len(hmout[1])):
          line=line.replace('SM'+v[t],hmout[1][t])
        for t in range(len(hmout[2])):
          line=line.replace('TM'+v[t],hmout[2][t])
        
        fout.write(line)
  
  zip_folder(outpath, outpath[:-1]+'.zip')
