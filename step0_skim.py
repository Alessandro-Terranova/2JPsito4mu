

#%%
import os
import ROOT

files_idx = {
    "Run1":
        ["data/Run2012B_MuOnia_NanoAODRun1.txt",
        "data/Run2012C_MuOnia_NanoAODRun1.txt"],
    "Run2":
        ["data/CMS_Run2016G_Charmonium_NANOAOD_UL2016_MiniAODv2_NanoAODv9-v1.txt",
        "data/CMS_Run2016H_Charmonium_NANOAOD_UL2016_MiniAODv2_NanoAODv9-v1.txt"]
}

def read_files_idx(files_idx):
    file_dict = {k : [] for k in files_idx.keys()}
    for era, era_list in files_idx.items():
        for txt in era_list:
            with open(os.path.join(os.path.dirname(__file__),txt), 'r') as f:
                file_dict[era] += [line.strip() for line in f if line.strip()]
    return file_dict

def skim(files_dict, save_path = "data/{eras}.root"):
    collections = ["Muon", "Dimu", "TrigObj"]
    for era, files_era in files_dict.items():
        print("Processing era:", era)
        df = ROOT.RDataFrame("Events", files_era)
        df = df.Filter("GoodLumisection").Filter("HLT_Dimuon0_Jpsi_Muon").Filter("nMuon>=4")
        df = df.Define("nMuon_p", "Sum(Muon_charge==1)")
        df = df.Define("nMuon_m", "Sum(Muon_charge==-1)")
        df = df.Filter("nMuon_p>=2 && nMuon_m>=2")
        list_of_columns = df.GetColumnNames()
        column_to_save=[]
        dimu_cols = []
        for branch in list_of_columns:
            for col in collections:
                if branch.startswith(f"{col}_") or branch==f"n{col}":
                    column_to_save.append(branch)
                if branch.startswith("Dimu_"):
                    dimu_cols.append(branch)
        for dimu_col in dimu_cols:
            if dimu_col != "Dimu_charge":
                df = df.Redefine(dimu_col, f"{dimu_col}[Dimu_charge==0]")
        df = df.Redefine(dimu_col, f"Dimu_charge[Dimu_charge==0]")
        df = df.Redefine("nDimu", "Dimu_charge.size()")
        df.Snapshot("Events", save_path.format(eras=era), column_to_save)

files_dict = read_files_idx(files_idx)
skim(files_dict)


