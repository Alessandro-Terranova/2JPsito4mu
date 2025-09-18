import os
import ROOT

# ROOT.EnableImplicitMT()  # Abilita il multithreading implicito in ROOT

# Creo il percorso per la cartella contente i file
path = os.path.dirname(__file__)



def FilterCollection(rdf, collection, new_col=None, mask=None, indices=None):
    """Funzione che filtra una collezione (colonna specifica) in un RDataFrame.
        Si può filtrare usando una maschera booleana (mask) o usando degli indici (indices).
        Se si specifica new_col, la collezione filtrata viene salvata con un nuovo nome.
        Se non si specifica new_col, la collezione originale viene sovrascritta.
        Non si possono specificare entrambe le opzioni mask e indices. 
    """

    exclude = ["Muon_nNano"]  # Escludo questa variabile dal filtro

    # Controllo che almeno uno tra mask e indices sia None, ma non entrambi, sennò solleva un errore
    if mask is None:
        assert indices is not None
    else:
        assert indices is None

    # Ottengo la lista di tutte le colonne nel DataFrame
    columns = rdf.GetColumnNames()

    redefine = True
    # Se new_col è specificato, controllo che non esistano già colonne con quel prefisso
    if new_col is not None:
        redefine = False
        assert not any([name.startswith(new_col + "_") for name in columns])
    else:
        new_col = collection # Se non specificato, sovrascrivo con la stessa colonna

    # Funzione per definire o ridefinire una colonna in base alla sua esistenza nel DataFrame
    # name è il nome della colonna, expr è l'espressione (mask o indices)
    def define_or_redefine(name, expr, redefine):
        if redefine:
            return rdf.Redefine(name, expr)
        else:
            return rdf.Define(name, expr)

    # Controllo mask
    if mask is not None:
        mask_name = f"__mask__{collection}__{new_col}"
        rdf = define_or_redefine(mask_name, mask, mask_name in columns)
        mask = mask_name  # aggiorna il riferimento

    # Controllo indices
    if indices is not None:
        indices_name = f"__indices__{collection}__{new_col}"
        rdf = define_or_redefine(indices_name, indices, indices_name in columns)
        indices = indices_name  # aggiorna il riferimento

    # Creo la lista delle variabili da filtrare senza includere quelle in exclude
    vars = [
        name
        for name in columns
        if name.startswith(collection + "_") and name not in exclude
    ]

    # Applico il filtro a tutte le variabili della collezione
    for var in vars:
        target = var if redefine else var.replace(collection, new_col)
        if mask:
            rdf = define_or_redefine(target, f"{var}[{mask}]", redefine)
        if indices:
            rdf = define_or_redefine(target, f"Take({var}, {indices})", redefine)

    # Definisco o ridefinisco la variabile che conta il numero di elementi nella collezione filtrata
    if mask:
        size = f"Sum({mask})"
    if indices:
        size = f"{indices}.size()"

    size_name = f"n{collection}" if redefine else f"n{new_col}"
    rdf = define_or_redefine(size_name, size, redefine)

    return rdf

if __name__ == "__main__":

    # Definisco il ROOT file da usare
    input = "Run1.root"

    #-------------------------------------------#
    
    # Creo il percorso per il file di input
    input_file = os.path.join(path, f"data/{input}")

    # Verifico che il file di input esista
    if not os.path.isfile(input_file):
        raise FileNotFoundError(f"Il file specificato non esiste: {input_file}")

    ROOT.gInterpreter.Declare('#include "utils.cpp"') # Importo le funzioni definite in utils.cpp
    
    #-------------------------------------------#

    # Creo il DataFrame dal file di input
    rdf = ROOT.RDataFrame("Events", f"{input_file}")


    # Applico i filtri

    rdf = FilterCollection(rdf, "TrigObj", new_col="TrigMu", mask="abs(TrigObj_id) == 13") # filtro i TrigObj

    # Definisco le variabili che indicano se i muoni sono triggermatched
    rdf = rdf.Define(
        "Muon_isTrigMatched",
        "dRpTMatch(Muon_eta, Muon_phi, Muon_pt, TrigMu_eta, TrigMu_phi, TrigMu_pt)",
    ) 

    # Definisco le variabili che indicano se i muoni sono in accettanza tight o loose
    rdf = rdf.Define(
        "Muon_isInTightAccept", "tightMuAcceptance(Muon_isTrigMatched, Muon_pt, Muon_eta)"
    )
    rdf = rdf.Define("Muon_isInLooseAccept", "looseMuAcceptance(Muon_pt, Muon_eta)")

    # Filtro gli eventi che hanno almeno due muoni in accettanza tight o loose
    rdf = rdf.Define(
        "AreMuInAccept", "MuonsAcceptance(Muon_isInTightAccept, Muon_isInLooseAccept)"
    ).Filter("AreMuInAccept")

    # Creo la collezione dei muoni che compongono i dimuon
    rdf = FilterCollection(rdf, "Muon", new_col="Dimu_m1", indices="Dimu_t1muIdx")
    rdf = FilterCollection(rdf, "Muon", new_col="Dimu_m2", indices="Dimu_t2muIdx")

    # Definisco la probabilità del vertice del dimuon ricostruendo il vertice della J/psi
    rdf = rdf.Define("Dimu_vertexProb", "ROOT::VecOps::Map(Dimu_chi2, getProb)")

    rdf = rdf.Define("Dimu_rapidity", "rapidity(Dimu_pt, Dimu_eta)")  # Definisco la rapidità del dimuon

    # Filtro i dimuon in base a vari criteri di qualità
    rdf = FilterCollection(
        rdf,
        "Dimu",
        mask="""(Dimu_m1_isGlobal || Dimu_m1_isTracker) && (Dimu_m2_isGlobal || Dimu_m2_isTracker) &&
            Dimu_m1_isGood && Dimu_m2_isGood &&
            Dimu_m1_nPix >= 2 && Dimu_m2_nPix >= 2 &&
            Dimu_m1_nValid >= 11 && Dimu_m2_nValid >= 11 &&
            Dimu_m1_dxy < 3 && Dimu_m2_dxy < 3 &&
            Dimu_m1_dz < 30 && Dimu_m2_dz < 30 &&
            Dimu_m1_Chi2 < 1.8 && Dimu_m2_Chi2 < 1.8 &&
            (Dimu_m1_isInTightAccept || Dimu_m1_isInLooseAccept) && (Dimu_m2_isInTightAccept || Dimu_m2_isInLooseAccept) &&
            Dimu_mass > 2.8 && Dimu_mass < 3.35 &&
            Dimu_vertexProb > 0.005 &&
            JpsiAcceptance(Dimu_pt, Dimu_rapidity)
            """,
    ).Filter("nDimu >= 2")

    rdf = rdf.Define("Dimu_CandidateIdx", "JpsiCandidates(Dimu_vertexProb, Dimu_t1muIdx, Dimu_t2muIdx)")  # Indici dei dimuon candidati J/psi
    #Sommo i candidati J/psi + 1 così i non candidati sono 0, e i candidati sono 1 o 2
    rdf = rdf.Filter("Sum(Dimu_CandidateIdx + 1) == 3")  # Filtro gli eventi con esattamente due candidati J/psi

    rdf = FilterCollection(rdf, "Dimu", new_col="Jpsi", mask = "Dimu_CandidateIdx != -1")  # Creo collezione dei dimuon candidati J/psi
    
    rdf = FilterCollection(rdf, "Jpsi", indices = "ROOT::VecOps::Argsort(Jpsi_CandidateIdx)")  # Ordino i dimuon J/psi in base agli indici dei candidati
    rdf = rdf.Filter("""
                     Jpsi_dl[0] > -0.05 && Jpsi_dl[0] < 0.1
                     """)
 
 
    rdf.Display(["Jpsi_mass", "Jpsi_CandidateIdx", "Jpsi_t1muIdx", "Jpsi_t2muIdx", "Jpsi_vertexProb" ]).Print()  # Stampo le variabili del DataFrame

