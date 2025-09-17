

#%%
import os
import ROOT

ROOT.EnableImplicitMT()  # Abilita il multithreading implicito in ROOT

# Definisco dizionari con i nomi dei file di input per le diverse ere
# Questi file contengono i nomi dei file NanoAOD da processare per ogni era
# I file sono organizzati in due categorie: Run1 e Run2
files_idx = {
    "Run1":
        ["data/Run2012B_MuOnia_NanoAODRun1.txt",
        "data/Run2012C_MuOnia_NanoAODRun1.txt"],
    "Run2":
        ["data/CMS_Run2016G_Charmonium_NANOAOD_UL2016_MiniAODv2_NanoAODv9-v1.txt",
        "data/CMS_Run2016H_Charmonium_NANOAOD_UL2016_MiniAODv2_NanoAODv9-v1.txt"]
}

path = os.path.dirname(__file__)

def read_files_idx(files_idx):
    """ Funzione che legge i file di input e restituisce un dizionario con le liste dei file per ogni era"""

    #creo un dizionario vuoto con le chiavi delle ere
    file_dict = {k : [] for k in files_idx.keys()}

    # Per ogni era, leggo i file di testo che contengono i nomi dei file
    for era, era_list in files_idx.items():
        # Ciclo su ogni file di testo nella lista dell'era
        for txt in era_list:
            # Apro il file e leggo le righe, aggiungendo i nomi dei file al dizionario
            with open(os.path.join(path,txt), 'r') as f:
                file_dict[era] += [line.strip() for line in f if line.strip()]
    return file_dict


def skim(files_dict, save_path = "data/{eras}.root"):
    """ Funzione che effettua lo skim dei file di input"""

    # Definisco le collezioni da usare
    collections = ["Muon", "Dimu", "TrigObj"]

    for era, files_era in files_dict.items():
        print("Processing era:", era)
        df = ROOT.RDataFrame("Events", files_era)
        df = df.Filter("GoodLumisection").Filter("HLT_Dimuon0_Jpsi_Muon").Filter("nMuon>=4")
        df = df.Define("nMuon_p", "Sum(Muon_charge==1)")
        df = df.Define("nMuon_m", "Sum(Muon_charge==-1)")
        df = df.Filter("nMuon_p>=2 && nMuon_m>=2")

        # Ottengo la lista di tutte le colonne nel DataFrame
        list_of_columns = df.GetColumnNames()
        # Creo liste per le colonne da salvare
        column_to_save=[]
        dimu_cols = []

        # Seleziono le colonne da salvare, includendo tutte le colonne delle collezioni di interesse
        for branch in list_of_columns:
            for col in collections:
                if branch.startswith(f"{col}_") or branch==f"n{col}":
                    column_to_save.append(branch)
            if branch.startswith("Dimu_"):
                    dimu_cols.append(branch)
        # Ridefinisco le colonne di Dimu per mantenere solo quelle con carica 0 (dimuoni neutri)
        for dimu_col in dimu_cols:
            if dimu_col != "Dimu_charge":
                df = df.Redefine(dimu_col, f"{dimu_col}[Dimu_charge==0]")
        df = df.Redefine("Dimu_charge", f"Dimu_charge[Dimu_charge==0]")
        df = df.Redefine("nDimu", "Dimu_charge.size()").Filter("nDimu>=2")


        print(f"Sto salvando: {column_to_save}")
        df.Snapshot("Events", save_path.format(eras=era), column_to_save)

if __name__ == '__main__':
    files_dict = read_files_idx(files_idx)

    skim(files_dict)


