#ifndef UTILS_CPP
#define UTILS_CPP
#include "TMath.h"
#include "ROOT/RVec.hxx"

// Definisco la funzione dRpTmatch che verifica se due oggetti sono matched in dR e pT
// I due oggetti sono rappresentati da due vettori di eta, phi e pt
// La funzione restituisce un vettore di booleani che indica se ogni oggetto del primo vettore ha un match nel secondo vettore
// I parametri di matching sono dRmax e dPtMax
// La funzione utilizza la funzione DeltaR di ROOT per calcolare la distanza in eta-phi
// La funzione utilizza la funzione fabs di ROOT per calcolare la differenza in pT
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

// Funzione che verifica se un muone è in accettanza tight
// La funzione prende in input il vettore dei muoni trigger matched, il vettore dei pt e il vettore degli eta
// Restituisce un vettore di booleani che indica se ogni muone è in accettanza tight
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

// Funzione che verifica se un muone è in accettanza loose, analogamente alla tight
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

// Funzione che verifica se l'evento ha almeno 4 muoni in accettanza
// Almeno 3 in accettanza tight e almeno 1 in accettanza loose
// La funzione prende in input il vettore dei booleani che indicano se un muone è in accettanza tight o loose
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

// Lambda function per calcolare la probabilità associata a un valore di chi2 con 1 grado di libertà
auto getProb = [](float chi2){return TMath::Prob(chi2,1);};

// Funzione per la rapidità della J/psi
// La funzione prende in input il vettore dei pt e il vettore degli eta
// Restituisce un vettore di float con le rapidità calcolate
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

// Funzione per l'accettanza della J/psi
// La funzione prende in input il vettore dei pt e il vettore delle rapidità
// Restituisce un vettore di booleani che indica se ogni J/psi è in accettanza
ROOT::VecOps::RVec<bool>
JpsiAcceptance(const ROOT::VecOps::RVec<float> &pt, const ROOT::VecOps::RVec<float> &y) {
  ROOT::VecOps::RVec<bool> acceptances(pt.size(), false);
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

// Funzione che resituisce un Rvecinteger con gli indici degli elementi candidati a formare la coppia di J/psi
// La funzione prende in input il vettore delle probabilità dei vertici, il vettore dei firstindex e il vettore dei secondindex
// Restituisce un Rvecinteger con gli indici dei due dimuon candidati a formare la coppia di J/psi
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

#endif // UTILS_CPP