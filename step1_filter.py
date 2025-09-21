"""
Script per filtrare gli eventi con due candidati J/psi che decadono in 4 muoni, secondo i criteri di selezioni richiesti.
I file di input sono file ROOT in formato NanoAOD, precedentemente skimmati.
"""
import os
import ROOT

ROOT.gInterpreter.Declare('#include "utils.cpp"') # Importo le funzioni definite in utils.cpp

ROOT.EnableImplicitMT()  # Abilita il multithreading implicito in ROOT

# Creo il percorso per la cartella contente i file
path = os.path.dirname(__file__)


def FilterCollection(rdf, collection, new_col=None, mask=None, indices=None):
    """

        Filtra una collezione di oggetti in un RDataFrame basato su una maschera booleana o su indici specifici.
        Può creare una nuova collezione o sovrascrivere quella esistente.

        Args:
            rdf (ROOT.RDataFrame): Il DataFrame da cui filtrare la collezione
            collection (str): Il prefisso delle colonne della collezione da filtrare (es. "Muon" per le colonne Muon_pt, Muon_eta, ecc.)
            new_col (str, optional): Il prefisso delle colonne della nuova collezione filtrata. Defaults to None.
            mask (str, optional): Espressione booleana per filtrare la collezione. Defaults to None.
            indices (str, optional): Espressione che restituisce gli indici per filtrare la collezione. Defaults to None.
        Returns:
            ROOT.RDataFrame: Il DataFrame con la collezione filtrata

    """

    exclude = ["Muon_nNano"]  # Escludo questa variabile dal filtro che ha lunghezza diversa dalle altre

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
    
    #-------------------------------------------#

    # Creo il DataFrame dal file di input
    rdf = ROOT.RDataFrame("Events", f"{input_file}")

#-----------------------------------------------------------#

    # Istogramma 2d del pt dei muoni con pt minore per evento vs la loro pseudo-rapidità prima della selezione
    #filtro la pseudorapidità in valore assoluto per avere un istogramma più leggibile
    rdf = rdf.Define("Muon_Index", "MuonPtOrdering(Muon_pt, Muon_trkId)")  # Creo un vettore con gli indici dei muoni
    rdf = FilterCollection(rdf, "Muon", new_col="MuonSorted", indices="ROOT::VecOps::Argsort(Muon_Index)")  # Creo collezione dei muoni ordinati per pt crescente
    rdf = FilterCollection(rdf, "MuonSorted", mask="MuonSorted_eta = abs(MuonSorted_eta)")  # Creo collezione dei muoni ordinati per pt crescente

    # Muon_lowestpt[0]
    rdf = rdf.Define("Muon0_lowestpt", "MuonSorted_pt[0]")
    rdf = rdf.Define("Muon0_lowestpteta", "MuonSorted_eta[0]")
    hist6 = rdf.Histo2D(("Muon0_lowestpt_vs_Muon0_lowestpteta", "Muon0_lowestpt_vs_Muon0_lowestpteta; Muon0_lowestpteta; Muon0_lowestpt", 30, 0, 3, 30, 0, 6), "Muon0_lowestpteta", "Muon0_lowestpt")  # Creo l'istogramma
    hist6.GetXaxis().SetTitle("|#eta|")
    hist6.GetYaxis().SetTitle("p_{T} [GeV/c]")
    hist6.SetTitle("Distribuzione muoni con p_{T} piu' basso per evento")
    hist6.SetStats(0)
    c6 = ROOT.TCanvas("c6", "c6", 800, 600)
    c6.SetRightMargin(0.15)
    hist6.Draw("COLZ")
    c6.SaveAs(os.path.join(path, "Muon0_lowestpt_vs_Muon0_lowestpteta.png"))

#-----------------------------------------------------------#

    ## Applico i filtri

    rdf = FilterCollection(rdf, "TrigObj", new_col="TrigMu", mask="abs(TrigObj_id) == 13") # filtro i TrigObj

    # Definisco i filtri che indicano se i muoni sono triggermatched
    rdf = rdf.Define(
        "Muon_isTrigMatched",
        "dRpTMatch(Muon_eta, Muon_phi, Muon_pt, TrigMu_eta, TrigMu_phi, TrigMu_pt)",
    ) 

    # Definisco i filtri che indicano se i muoni sono in accettanza tight o loose
    rdf = rdf.Define(
        "Muon_isInTightAccept", "tightMuAcceptance(Muon_isTrigMatched, Muon_pt, Muon_eta)"
    )
    rdf = rdf.Define("Muon_isInLooseAccept", "looseMuAcceptance(Muon_pt, Muon_eta)")

    # Filtro gli eventi che hanno almeno 4 muoni in accettanza tight o loose
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
 
    rdf = rdf.Define("DistanceBetweenJpsi", """sqrt( (Jpsi_z[0] - Jpsi_z[1])*(Jpsi_z[0] - Jpsi_z[1]) + 
                                                (Jpsi_y[0] - Jpsi_y[1])*(Jpsi_y[0] - Jpsi_y[1]) + 
                                                (Jpsi_x[0] - Jpsi_x[1])*(Jpsi_x[0] - Jpsi_x[1]) )
                     """)  # Definisco la distanza tra i due vertici dei J/psi


    # rdf.Display(["Jpsi_mass", "Jpsi_CandidateIdx", "Jpsi_t1muIdx", "Jpsi_t2muIdx", "Jpsi_vertexProb" ]).Print()  # Stampo le variabili del DataFrame

    hist = rdf.Histo1D(("DistanceBetweenJpsi", "Distanza tra i candidati J/#psi", 100, 0, 1), "DistanceBetweenJpsi")  # Creo l'istogramma della distanza tra i due vertici dei J/psi
    hist.GetXaxis().SetTitle("#Delta r [cm]")
    hist.GetYaxis().SetTitle("Numero di eventi candidati")
    hist.SetStats(0)
    c0 = ROOT.TCanvas("c0", "c0", 800, 600)
    hist.Draw()
    c0.SaveAs(os.path.join(path, "DistanceBetweenJpsi.png"))

    rdf = rdf.Filter("DistanceBetweenJpsi < 0.1")  # Filtro gli eventi con distanza tra i due vertici dei J/psi minore di 0.1 cm

    print(f"Numero di eventi dopo tutti i filtri: {rdf.Count().GetValue()}")  # Stampo il numero di eventi rimanenti dopo tutti i filtri

#-----------------------------------------------------------#

    # Plot: creo gli istogrammi delle variabili di interesse
    # Jpsi_mass[0]
    rdf = rdf.Define("Jpsi0_mass", "Jpsi_mass[0]")
    hist1 = rdf.Histo1D(("Jpsi_mass_1", "Massa J/#psi 1", 100, 2.8, 3.35), "Jpsi0_mass")  # Creo l'istogramma della massa del primo J/psi
    hist1.GetXaxis().SetTitle("M_{J/#psi} [GeV/c^{2}]")
    hist1.GetYaxis().SetTitle("Numero di candidati eventi ")
    c1 = ROOT.TCanvas("c1", "c1", 800, 600)
    hist1.Draw()
    c1.SaveAs(os.path.join(path, "Jpsi_mass_1.png"))

    # Jpsi_mass[1]
    rdf = rdf.Define("Jpsi1_mass", "Jpsi_mass[1]")
    hist2 = rdf.Histo1D(("Jpsi_mass_2", "Massa J/#psi 2", 100, 2.8, 3.35), "Jpsi1_mass")  # Creo l'istogramma della massa del secondo J/psi
    hist2.GetXaxis().SetTitle("M_{J/#psi} [GeV/c^{2}]")
    hist2.GetYaxis().SetTitle("Numero di candidati eventi ")
    c2 = ROOT.TCanvas("c2", "c2", 800, 600)
    hist2.Draw()
    c2.SaveAs(os.path.join(path, "Jpsi_mass_2.png"))

    # Jpsi_dl
    hist4 = rdf.Histo1D(("Jpsi_dl", "Lunghezza di decadimento Jpsi", 100, -0.05, 0.1), "Jpsi_dl")  # Creo l'istogramma della distanza di volo del primo J/psi
    hist4.GetXaxis().SetTitle("L [cm]")
    hist4.GetYaxis().SetTitle("Numero di candidati eventi ")
    c4 = ROOT.TCanvas("c4", "c4", 800, 600)
    hist4.Draw()
    c4.SaveAs(os.path.join(path, "Jpsi_dl.png"))

    # Jpsi_mass[0] vs Jpsi_mass[1]
    hist3 = rdf.Histo2D(("Jpsi_mass_1_vs_Jpsi_mass_2", "Jpsi_mass_1_vs_Jpsi_mass_2; Jpsi0_mass; Jpsi1_mass", 100, 2.8, 3.35, 100, 2.8, 3.35), "Jpsi0_mass", "Jpsi1_mass")  # Creo l'istogramma 2D della massa del primo J/psi vs la massa del secondo J/psi
    hist3.SetStats(0)
    hist3.GetXaxis().SetTitle("M_{J/#psi 1} [GeV/c^{2}]")
    hist3.GetYaxis().SetTitle("M_{J/#psi 2} [GeV/c^{2}]")
    hist3.SetTitle("Distribuzione delle masse delle due J/#psi per cardidato evento")
    c3 = ROOT.TCanvas("c3", "c3", 800, 600)
    hist3.Draw("COLZ")
    c3.SaveAs(os.path.join(path, "Jpsi_mass_1_vs_Jpsi_mass_2.png"))

    # Jpsi_pt vs Jpsi_rapidity
    #filtro la rapidità in valore assoluto per avere un istogramma più leggibile
    rdf = FilterCollection(rdf, "Jpsi", mask="Jpsi_rapidity = abs(Jpsi_rapidity)")  # Creo collezione dei dimuon candidati J/psi")
    hist5 = rdf.Histo2D(("Jpsi_pt_vs_Jpsi_rapidity", "Jpsi_pt_vs_Jpsi_rapidity; Jpsi_rapidity; Jpsi_pt", 30, 0, 2.2, 30, 4.5, 10), "Jpsi_rapidity", "Jpsi_pt")  # Creo l'istogramma 2D della rapidità del primo J/psi vs il pt del primo J/psi
    hist5.SetStats(0)
    hist5.GetXaxis().SetTitle("J/#psi |y|")
    hist5.GetYaxis().SetTitle("J/#psi p_{T} [GeV/c]")
    hist5.SetTitle("Distribuzione p_{T} vs |y| delle J/#psi")
    c5 = ROOT.TCanvas("c5", "c5", 800, 600)
    hist5.Draw("COLZ")
    c5.SaveAs(os.path.join(path, "Jpsi_pt_vs_Jpsi_rapidity.png"))
