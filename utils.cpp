/*
Script contente delle funzioni utili per l'analisi
*/

#ifndef UTILS_CPP
#define UTILS_CPP
#include "TMath.h"
#include "ROOT/RVec.hxx"

//-------------------------------------------------------------------------------------//

/*

Funzione che verifica se ogni oggetto in un vettore ha un match in un altro vettore in dR e dpT

Args:
    eta1, phi1, pt1: vettori di float rappresentanti il primo oggetto
    eta2, phi2, pt2: vettori di float rappresentanti il secondo oggetto
    dRmax: float, distanza massima in eta-phi per considerare un match (default 0.1)
    dPtMax: float, differenza massima in pT per considerare un match (default 10)
Returns:
    vettore di booleani che indica se ogni oggetto del primo vettore ha un match nel secondo vettore

*/
ROOT::VecOps::RVec<bool> dRpTMatch(const ROOT::VecOps::RVec<float> &eta1,
                                   const ROOT::VecOps::RVec<float> &phi1,
                                   const ROOT::VecOps::RVec<float> &pt1,
                                   const ROOT::VecOps::RVec<float> &eta2,
                                   const ROOT::VecOps::RVec<float> &phi2,
                                   const ROOT::VecOps::RVec<float> &pt2,
                                   float dRmax = 0.1, float dPtMax = 10) {
  ROOT::VecOps::RVec<bool> matches(eta1.size(), false);
  for (size_t i = 0; i < eta1.size(); ++i)
    for (size_t j = 0; j < eta2.size(); ++j)
      if (ROOT::VecOps::DeltaR(eta1[i], eta2[j], phi1[i], phi2[j]) < dRmax &&
          fabs(pt1[i] - pt2[j]) < dPtMax) {
        matches[i] = true;
        break;
      };
  return matches;
};


//-------------------------------------------------------------------------------------//

/*

Funzione che verifica se un muone è in accettanza tight

Args:
    isTrigMatched: vettore di booleani che indica se il muone è trigger matched
    pt: vettore di float rappresentante il pT dei muoni
    eta: vettore di float rappresentante l'eta dei muoni
Returns:
    vettore di booleani che indica se ogni muone è in accettanza tight

*/
ROOT::VecOps::RVec<bool>
tightMuAcceptance(const ROOT::VecOps::RVec<bool> &isTrigMatched,
                  const ROOT::VecOps::RVec<float> &pt,
                  const ROOT::VecOps::RVec<float> &eta) {
  ROOT::VecOps::RVec<bool> acceptances(pt.size(), false);
  for (size_t i = 0; i < pt.size(); ++i) {
    if (!isTrigMatched[i]) {
      continue;
    }

    if ((fabs(eta[i]) < 1.2 && pt[i] > 3.5) ||
        (fabs(eta[i]) > 1.2 && eta[i] < 1.6 &&
         pt[i] > -15 * fabs(eta[i]) / 4 + 8) || // retta che passa per i punti (1.2,3.5) e (1.6,2.5)
        (fabs(eta[i]) > 1.6 && fabs(eta[i]) < 2.4 && pt[i] > 3.5)) {
      acceptances[i] = true;
    }
  }
  return acceptances;
}

//-------------------------------------------------------------------------------------//

/*

Funzione che verifica se un muone è in accettanza loose, analoga alla tight

*/
ROOT::VecOps::RVec<bool>
looseMuAcceptance(const ROOT::VecOps::RVec<float> &pt,
                  const ROOT::VecOps::RVec<float> &eta) {
  ROOT::VecOps::RVec<bool> acceptances(pt.size(), false);
  for (size_t i = 0; i < pt.size(); ++i) {
    if (fabs(eta[i]) < 1.2 && pt[i] > 3 ||
        fabs(eta[i]) >= 1.2 && fabs(eta[i]) < 2.4 &&
            sqrt(pt[i] * pt[i] + pt[i] * sinh(fabs(eta[i]))) > 3) { 
      acceptances[i] = true;
    }
  }
  return acceptances;
}

//-------------------------------------------------------------------------------------//

/*

Funzione che verifica se un evento con 4 muoni soddisfa i criteri di accettanza:
- almeno 4 muoni in accettanza tight
- almeno 1 muone in accettanza loose (ma non tight)

Args:
    isInTightAccept: vettore di booleani che indica se ogni muone è in accettanza tight
    isInLooseAccept: vettore di booleani che indica se ogni muone è in accettanza loose
Returns:
    booleano che indica se l'evento soddisfa i criteri di accettanza

*/
bool MuonsAcceptance(const ROOT::VecOps::RVec<bool> &isInTightAccept,
                     const ROOT::VecOps::RVec<bool> &isInLooseAccept) {
  
  if (Sum(isInTightAccept) >= 4 and Sum(isInLooseAccept) >= 1) { // caso in cui tutti e quattro i muoni sono in accettanza tight
    return true;
  } else if (Sum(isInTightAccept) == 3 and
             Sum(!isInTightAccept && isInLooseAccept) >= 1) { // caso in cui tre muoni sono in accettanza tight e uno in accettanza loose, ma non tight
    return true;
    } else {
    return false;
      }
}

//-------------------------------------------------------------------------------------//

/* 

Lambda function per calcolare la probabilità di vertex fit associata a un valore di chi2 con 1 grado di libertà

Args:
    chi2: float, valore di chi2
Returns:
    float, probabilità associata al valore di chi2

*/

auto getProb = [](float chi2){return TMath::Prob(chi2,1);};

/*
Funzione che calcola la rapidità di un oggetto a partire da pt ed eta

Args:
    pt: vettore di float rappresentante il pT delle J/psi
    eta: vettore di float rappresentante l'eta delle J/psi
Returns:
    vettore di float che indica la rapidità delle J/psi

*/
ROOT::VecOps::RVec<float>
rapidity(const ROOT::VecOps::RVec<float> &pt, const ROOT::VecOps::RVec<float> &eta) {
  ROOT::VecOps::RVec<float> y(pt.size());
  for ( int i = 0; i < pt.size(); ++i) {
    float pz = pt[i] * sinh(eta[i]);
    float energy = sqrt(pt[i] * pt[i] + pz*pz + 0.1056580 * 0.1056580);
    y[i] = 0.5 * log((energy + pz) / (energy - pz));
  }
  return y;
}

//-------------------------------------------------------------------------------------//

/*

Funzione che verifica se una J/psi è in accettanza

Args:
    pt: vettore di float rappresentante il pT delle J/psi
    y: vettore di float rappresentante la rapidità delle J/psi
Returns:
    vettore di booleani che indica se ogni J/psi è in accettanza

*/
ROOT::VecOps::RVec<bool>
JpsiAcceptance(const ROOT::VecOps::RVec<float> &pt, const ROOT::VecOps::RVec<float> &y) {
  ROOT::VecOps::RVec<bool> acceptances(pt.size(), false); // Inizializzo il vettore delle accettanze a false
  for (size_t i = 0; i < pt.size(); ++i) {
    if ((fabs(y[i]) < 1.2 &&  pt[i] > 6.5) ||
        (fabs(y[i]) >= 1.2 && fabs(y[i]) < 1.43 &&
         pt[i] > -200 * fabs(y[i])/23  + 779/46) || // retta che passa per i punti (1.2,6.5) e (1.43,4.5)
        (fabs(y[i]) >= 1.43 && fabs(y[i]) < 2.2 && pt[i] > 4.5)) {
      acceptances[i] = true;
    }
  }
  return acceptances;
}


//-------------------------------------------------------------------------------------//

/*

Funzione che seleziona i due migliori candidati J/psi in un evento con 4 muoni in base alla probabilità del vertex fit

Args:
    Vtxprob: vettore di float rappresentante la probabilità del vertex fit dei dimuon
    firstindex: vettore di interi rappresentante l'indice del primo muone del dimuon
    secondindex: vettore di interi rappresentante l'indice del secondo muone del dimuon
Returns:
    vettore di interi che indica gli indici dei due candidati J/psi (0 e 1), -1 se non ci sono candidati

*/
ROOT::VecOps::RVec<int> JpsiCandidates(ROOT::VecOps::RVec<float> Vtxprob, ROOT::VecOps::RVec<int> firstindex, ROOT::VecOps::RVec<int> secondindex) {
  int nDimu = Vtxprob.size();
  ROOT::VecOps::RVec<int> CandidateIndexs(nDimu, -1); // Inizializzo il vettore degli indici dei candidati a -1
  ROOT::VecOps::RVec<int> SortedIndexs = ROOT::VecOps::Argsort(-Vtxprob); // Ordino gli indici in modo decrescente in base alla probabilità del vertice
  // ROOT::VecOps::RVec<float> VertexProbSorted = ROOT::VecOps::Take(Vtxprob, SortedIndexs); // Creo un vettore con le probabilità ordinate
  ROOT::VecOps::RVec<int> FirstIndexSorted = ROOT::VecOps::Take(firstindex, SortedIndexs); // Creo un vettore con i firstindex ordinati
  ROOT::VecOps::RVec<int> SecondIndexSorted = ROOT::VecOps::Take(secondindex, SortedIndexs); // Creo un vettore con i secondindex ordinati
  for(int i1 = 0; i1 < nDimu; ++i1) {
    for (int i2 = i1+1; i2 < nDimu; ++i2) { // Non confronto mai un dimuon con se stesso
      if (FirstIndexSorted[i1] != FirstIndexSorted[i2] && FirstIndexSorted[i1] != SecondIndexSorted[i2] && 
          SecondIndexSorted[i1] != FirstIndexSorted[i2] && SecondIndexSorted[i1] != SecondIndexSorted[i2]) {
        // Se i due dimuon condividono almeno un muone, esco dal ciclo interno
        CandidateIndexs[SortedIndexs[i1]] = 0; // Segno il dimuon come candidato
        CandidateIndexs[SortedIndexs[i2]] = 1; // Segno il dimuon come candidato
        return CandidateIndexs; // Esco dalla funzione restituendo gli indici dei candidati
      }
    }
  }  
  return CandidateIndexs;
}

//-------------------------------------------------------------------------------------//


/*

Funzione che ordina i muoni in base al loro pT e restituisce gli indici ordinati

Args:
    muon_pt: vettore di float rappresentante il pT dei muoni
    muon_index: vettore di interi rappresentante l'indice dei muoni
Returns:
    vettore di interi che indica gli indici dei muoni ordinati per pT

*/
ROOT::VecOps::RVec<int> MuonPtOrdering(const ROOT::VecOps::RVec<float> &muon_pt, const ROOT::VecOps::RVec<int> &muon_index) {
  int nMuons = muon_pt.size();
  ROOT::VecOps::RVec<int> MuonIndexSorted(nMuons, 0); // Inizializzo il vettore degli indici dei muoni ordinati a 0
  ROOT::VecOps::RVec<int> SortedIndexs = ROOT::VecOps::Argsort(muon_pt); // Ordino gli indici in modo crescente in base al pt
  ROOT::VecOps::RVec<int> IndexSortedMuon = ROOT::VecOps::Take(muon_index, SortedIndexs); // Creo un vettore con gli indici dei muoni ordinati per pt crescente
  for(int i = 0; i < MuonIndexSorted.size(); ++i) {
    MuonIndexSorted[i] = muon_index[SortedIndexs[i]];
  }
  return MuonIndexSorted;
}

#endif // UTILS_CPP