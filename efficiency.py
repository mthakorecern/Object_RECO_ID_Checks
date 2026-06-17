import os
import json
from collections import OrderedDict
import multiprocessing as mp
import math
import ROOT
ROOT.PyConfig.IgnoreCommandLineOptions = True

from PhysicsTools.NanoAODTools.postprocessing.framework.postprocessor import *
from PhysicsTools.NanoAODTools.postprocessing.framework.eventloop import Module
from PhysicsTools.NanoAODTools.postprocessing.framework.datamodel import Collection
from PhysicsTools.NanoAODTools.postprocessing.examples.exampleModule import *

class TruthMatchLeptonEfficiencyProducer(Module):

    def __init__(
        self,
        isData=False,
        debug=False,
        higgs_pdgid=25,
        resonance_pdgid=35,
        max_dr_ele=0.3,
        max_dr_muon=0.3,
        max_dr_hps_tau=0.1,
        max_dr_boosted_tau=0.3,
        json_path=None,
    ):
        self.isData = isData
        self.isMC = not isData
        self.debug = debug

        self.higgs_pdgid = abs(higgs_pdgid)
        self.resonance_pdgid = abs(resonance_pdgid) if resonance_pdgid is not None else None
        self.max_dr_ele = float(max_dr_ele)
        self.max_dr_muon = float(max_dr_muon)
        self.max_dr_hps_tau = float(max_dr_hps_tau)
        self.max_dr_boosted_tau = float(max_dr_boosted_tau)

        self.json_path = json_path

        self.file_summaries = {}
        self.current_file = None
        self.current_summary = None

    def beginJob(self):
        pass

    def endJob(self):
        if self.json_path is not None:
            outdir = os.path.dirname(os.path.abspath(self.json_path))
            if outdir and not os.path.exists(outdir):
                os.makedirs(outdir)

            payload = {
                "higgs_pdgid": self.higgs_pdgid,
                "resonance_pdgid": self.resonance_pdgid,
                "matching": {
                    "ele_deltaR_max": self.max_dr_ele,
                    "muon_deltaR_max": self.max_dr_muon,
                    "hps_tau_deltaR_max": self.max_dr_hps_tau,
                    "boosted_tau_deltaR_max": self.max_dr_boosted_tau,
                },
                "files": self.file_summaries,
            }

            with open(self.json_path, "w") as f:
                json.dump(payload, f, indent=4, sort_keys=True)

    def beginFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        self.out = wrappedOutputTree
        self.current_file = inputFile.GetName()

        self.current_summary = {
            "events_processed": 0,
            "raw_counts_no_fiducial_cuts": {
                "gen_particles_from_higgs_tau": {
                    "electrons": 0,
                    "muons": 0,
                    "visible_hadronic_taus": 0,
                },

                "gen_particles_with_matched_reco": {
                    "electrons_to_reco_electrons": 0,
                    "muons_to_reco_muons": 0,
                    "visible_taus_to_hps_taus": 0,
                    "visible_taus_to_boosted_taus": 0,
                    "electron_duplicate_match_attempts_before_unique_assignment": 0,
                    "events_with_electron_duplicate_match_attempts": 0,
                },
            },
        }

        ## GEN vis Taus kinematic variables
        self.out.branch("nGenVisTauFromHiggsTau", "I")
        self.out.branch("GenVisTauFromHiggsTau_pt", "F", lenVar="nGenVisTauFromHiggsTau")
        self.out.branch("GenVisTauFromHiggsTau_eta", "F", lenVar="nGenVisTauFromHiggsTau")
        self.out.branch("GenVisTauFromHiggsTau_phi", "F", lenVar="nGenVisTauFromHiggsTau")
        self.out.branch("GenVisTauFromHiggsTau_mass", "F", lenVar="nGenVisTauFromHiggsTau")
        self.out.branch("GenVisTauFromHiggsTau_status", "I", lenVar="nGenVisTauFromHiggsTau")
        self.out.branch("GenVisTauFromHiggsTau_genVisTauIdx", "I", lenVar="nGenVisTauFromHiggsTau")
        self.out.branch("GenVisTauFromHiggsTau_tauAncestorIdx", "I", lenVar="nGenVisTauFromHiggsTau")
        self.out.branch("GenVisTauFromHiggsTau_higgsAncestorIdx", "I", lenVar="nGenVisTauFromHiggsTau")
        self.out.branch("GenVisTauFromHiggsTau_resonanceAncestorIdx", "I", lenVar="nGenVisTauFromHiggsTau")

        # Matched HPS Tau kinematic variables
        self.out.branch("GenVisTauFromHiggsTau_matchedRecoTauIdx", "I", lenVar="nGenVisTauFromHiggsTau")
        self.out.branch("GenVisTauFromHiggsTau_hasMatchedRecoTau", "I", lenVar="nGenVisTauFromHiggsTau")
        self.out.branch("GenVisTauFromHiggsTau_matchedRecoTau_deltaR", "F", lenVar="nGenVisTauFromHiggsTau")
        self.out.branch("GenVisTauFromHiggsTau_matchedRecoTau_pt", "F", lenVar="nGenVisTauFromHiggsTau")
        self.out.branch("GenVisTauFromHiggsTau_matchedRecoTau_eta", "F", lenVar="nGenVisTauFromHiggsTau")
        self.out.branch("GenVisTauFromHiggsTau_matchedRecoTau_phi", "F", lenVar="nGenVisTauFromHiggsTau")
        self.out.branch("GenVisTauFromHiggsTau_matchedRecoTau_dz", "F", lenVar="nGenVisTauFromHiggsTau")
        self.out.branch("GenVisTauFromHiggsTau_matchedRecoTau_mass", "F", lenVar="nGenVisTauFromHiggsTau")
        self.out.branch("GenVisTauFromHiggsTau_matchedRecoTau_decayMode", "I", lenVar="nGenVisTauFromHiggsTau")
        self.out.branch("GenVisTauFromHiggsTau_matchedRecoTau_idDecayModeNewDMs", "I", lenVar="nGenVisTauFromHiggsTau")
        self.out.branch("GenVisTauFromHiggsTau_matchedRecoTau_idDeepTau2018v2p5VSjet", "I", lenVar="nGenVisTauFromHiggsTau")
        self.out.branch("GenVisTauFromHiggsTau_matchedRecoTau_idDeepTau2018v2p5VSe", "I", lenVar="nGenVisTauFromHiggsTau")
        self.out.branch("GenVisTauFromHiggsTau_matchedRecoTau_idDeepTau2018v2p5VSmu", "I", lenVar="nGenVisTauFromHiggsTau")
        self.out.branch("nGenMatchedRecoTauFromHiggsTau", "I")


        # Matched boosted Tau kinematic variables
        self.out.branch("GenVisTauFromHiggsTau_matchedRecoBoostedTauIdx", "I", lenVar="nGenVisTauFromHiggsTau")
        self.out.branch("GenVisTauFromHiggsTau_hasMatchedRecoBoostedTau", "I", lenVar="nGenVisTauFromHiggsTau")
        self.out.branch("GenVisTauFromHiggsTau_matchedRecoBoostedTau_deltaR", "F", lenVar="nGenVisTauFromHiggsTau")
        self.out.branch("GenVisTauFromHiggsTau_matchedRecoBoostedTau_pt", "F", lenVar="nGenVisTauFromHiggsTau")
        self.out.branch("GenVisTauFromHiggsTau_matchedRecoBoostedTau_eta", "F", lenVar="nGenVisTauFromHiggsTau")
        self.out.branch("GenVisTauFromHiggsTau_matchedRecoBoostedTau_phi", "F", lenVar="nGenVisTauFromHiggsTau")
        self.out.branch("GenVisTauFromHiggsTau_matchedRecoBoostedTau_mass", "F", lenVar="nGenVisTauFromHiggsTau")
        self.out.branch("GenVisTauFromHiggsTau_matchedRecoBoostedTau_decayMode", "I", lenVar="nGenVisTauFromHiggsTau")
        self.out.branch("GenVisTauFromHiggsTau_matchedRecoBoostedTau_idAntiEle2018", "I", lenVar="nGenVisTauFromHiggsTau")
        self.out.branch("GenVisTauFromHiggsTau_matchedRecoBoostedTau_idAntiMu", "I", lenVar="nGenVisTauFromHiggsTau")
        self.out.branch("GenVisTauFromHiggsTau_matchedRecoBoostedTau_idMVAnewDM2017v2", "I", lenVar="nGenVisTauFromHiggsTau")
        self.out.branch("GenVisTauFromHiggsTau_matchedRecoBoostedTau_rawBoostedDeepTauRunIIv2p0VSjet", "F", lenVar="nGenVisTauFromHiggsTau")
        self.out.branch("nGenMatchedRecoBoostedTauFromHiggsTau", "I")

        ## GEN electrons kinematic variables
        self.out.branch("nGenElectronFromHiggsTau", "I")
        self.out.branch("GenElectronFromHiggsTau_pt", "F", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_eta", "F", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_phi", "F", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_mass", "F", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_charge", "I", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_genPartIdx", "I", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_tauAncestorIdx", "I", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_higgsAncestorIdx", "I", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_resonanceAncestorIdx", "I", lenVar="nGenElectronFromHiggsTau")

        ## RECO electron's (matched to the GEN electron) kinematic variables
        self.out.branch("nGenMatchedRecoElectronFromHiggsTau", "I")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectronIdx", "I", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_hasMatchedRecoElectron", "I", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_pt", "F", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_eta", "F", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_phi", "F", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_superclusterEta", "F", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_vidNestedWPBitmap", "I", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_vidNestedWPBitmapHEEP", "I", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_cutBased", "I", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_cutBased_HEEP", "I", lenVar="nGenElectronFromHiggsTau")

        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_genPartFlav", "I", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_pfRelIso03_all", "F", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_pfRelIso03_chg", "F", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_pfRelIso04_all", "F", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_miniPFRelIso_all", "F", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_miniPFRelIso_chg", "F", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_convVeto", "I", lenVar="nGenElectronFromHiggsTau",)
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_sieie", "F",lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_deltaEtaSC", "F",lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_eInvMinusPInv", "F", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_hoe", "F", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_lostHits","I",lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_deltaR","F", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_rawEnergy","F", lenVar="nGenElectronFromHiggsTau")

        
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_nearestRecoTau_deltaR","F", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_nearestRecoTauIdx", "I", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_nearestRecoTau_pt", "F", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_nearestRecoTau_decayMode","I", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_nearestRecoBoostedTau_deltaR","F", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_nearestRecoBoostedTauIdx", "I", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_nearestRecoBoostedTau_pt", "F", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_nearestRecoBoostedTau_decayMode", "I", lenVar="nGenElectronFromHiggsTau")

        # Nearest GEN-matched RECO HPS Tau to matched RECO electron
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoTau_deltaR", "F", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoTauIdx", "I", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoTau_pt", "F", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoTau_eta", "F", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoTau_phi", "F", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoTau_decayMode", "I", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoTau_idDecayModeNewDMs", "I", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoTau_idDeepTau2018v2p5VSjet", "I", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoTau_idDeepTau2018v2p5VSe", "I", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoTau_idDeepTau2018v2p5VSmu", "I", lenVar="nGenElectronFromHiggsTau")


        # Nearest GEN-matched RECO boostedTau to matched RECO electron
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_deltaR", "F", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoBoostedTauIdx", "I", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_pt", "F", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_eta", "F", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_phi", "F", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_decayMode", "I", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_idAntiEle2018", "I", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_idAntiMu", "I", lenVar="nGenElectronFromHiggsTau")
        self.out.branch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_rawBoostedDeepTauRunIIv2p0VSjet", "F", lenVar="nGenElectronFromHiggsTau")


        ## GEN muons kinematic variables
        self.out.branch("nGenMuonFromHiggsTau", "I")
        self.out.branch("nGenMatchedRecoMuonFromHiggsTau", "I")
        self.out.branch("GenMuonFromHiggsTau_pt", "F", lenVar="nGenMuonFromHiggsTau")
        self.out.branch("GenMuonFromHiggsTau_eta", "F", lenVar="nGenMuonFromHiggsTau")
        self.out.branch("GenMuonFromHiggsTau_phi", "F", lenVar="nGenMuonFromHiggsTau")
        self.out.branch("GenMuonFromHiggsTau_mass", "F", lenVar="nGenMuonFromHiggsTau")
        self.out.branch("GenMuonFromHiggsTau_charge", "I", lenVar="nGenMuonFromHiggsTau")
        self.out.branch("GenMuonFromHiggsTau_genPartIdx", "I", lenVar="nGenMuonFromHiggsTau")
        self.out.branch("GenMuonFromHiggsTau_tauAncestorIdx", "I", lenVar="nGenMuonFromHiggsTau")
        self.out.branch("GenMuonFromHiggsTau_higgsAncestorIdx", "I", lenVar="nGenMuonFromHiggsTau")
        self.out.branch("GenMuonFromHiggsTau_resonanceAncestorIdx", "I", lenVar="nGenMuonFromHiggsTau")

        ## RECO muon's (matched to the GEN muon) kinematic variables
        self.out.branch("GenMuonFromHiggsTau_matchedRecoMuonIdx", "I", lenVar="nGenMuonFromHiggsTau")
        self.out.branch("GenMuonFromHiggsTau_hasMatchedRecoMuon", "I", lenVar="nGenMuonFromHiggsTau")
        self.out.branch("GenMuonFromHiggsTau_matchedRecoMuon_pt", "F", lenVar="nGenMuonFromHiggsTau")
        self.out.branch("GenMuonFromHiggsTau_matchedRecoMuon_eta", "F", lenVar="nGenMuonFromHiggsTau")
        self.out.branch("GenMuonFromHiggsTau_matchedRecoMuon_phi", "F", lenVar="nGenMuonFromHiggsTau")
        self.out.branch("GenMuonFromHiggsTau_matchedRecoMuon_looseId", "I", lenVar="nGenMuonFromHiggsTau")
        self.out.branch("GenMuonFromHiggsTau_matchedRecoMuon_mediumId", "I", lenVar="nGenMuonFromHiggsTau")
        self.out.branch("GenMuonFromHiggsTau_matchedRecoMuon_tightId", "I", lenVar="nGenMuonFromHiggsTau")
        self.out.branch("GenMuonFromHiggsTau_matchedRecoMuon_genPartFlav", "I", lenVar="nGenMuonFromHiggsTau")

       

    def endFile(self, inputFile, outputFile, inputTree, wrappedOutputTree):
        if self.current_file is not None and self.current_summary is not None:
            self.file_summaries[self.current_file] = dict(self.current_summary)

    def _valid_index(self, idx, coll):
        return 0 <= int(idx) < len(coll)

    def _mother_idx(self, genparts, idx):
        if not self._valid_index(idx, genparts):
            return -1

        mom = int(genparts[idx].genPartIdxMother)
        if self._valid_index(mom, genparts):
            return mom

        return -1

    def _walk_to_first_ancestor(self, genparts, start_idx, target_pdgid):
        target_pdgid = abs(int(target_pdgid))
        visited = set()
        idx = int(start_idx)

        while self._valid_index(idx, genparts) and idx not in visited:
            visited.add(idx)

            if abs(int(genparts[idx].pdgId)) == target_pdgid:
                return idx

            idx = self._mother_idx(genparts, idx)

        return -1

    def _find_tau_higgs_resonance_ancestors_from_lepton(self, genparts, start_idx):

        visited = set()
        idx = self._mother_idx(genparts, start_idx)

        tau_idx = -1
        higgs_idx = -1
        resonance_idx = -1

        while self._valid_index(idx, genparts) and idx not in visited:
            visited.add(idx)
            pdgid = abs(int(genparts[idx].pdgId))

            if pdgid == 15 and tau_idx < 0:
                tau_idx = idx

            if tau_idx >= 0 and pdgid == self.higgs_pdgid and higgs_idx < 0:
                higgs_idx = idx
                if self.resonance_pdgid is None:
                    break

            if higgs_idx >= 0 and self.resonance_pdgid is not None:
                if pdgid == self.resonance_pdgid:
                    resonance_idx = idx
                    break

            idx = self._mother_idx(genparts, idx)

        return tau_idx, higgs_idx, resonance_idx

    def _find_higgs_resonance_from_tau(self, genparts, tau_idx):

        if not self._valid_index(tau_idx, genparts):
            return -1, -1

        higgs_idx = self._walk_to_first_ancestor(
            genparts,
            tau_idx,
            self.higgs_pdgid,
        )

        resonance_idx = -1

        if higgs_idx >= 0 and self.resonance_pdgid is not None:
            resonance_idx = self._walk_to_first_ancestor(
                genparts,
                self._mother_idx(genparts, higgs_idx),
                self.resonance_pdgid,
            )

        if self.resonance_pdgid is not None and resonance_idx < 0:
            return -1, -1

        return higgs_idx, resonance_idx

    def _is_stable_gen_lepton_from_higgs_tau(self, genparts, idx, abs_lepton_pdgid):
        if not self._valid_index(idx, genparts):
            return False, -1, -1, -1

        gp = genparts[idx]

        if abs(int(gp.pdgId)) != abs_lepton_pdgid:
            return False, -1, -1, -1

        if int(getattr(gp, "status", -999)) != 1:
            return False, -1, -1, -1

        tau_idx, higgs_idx, resonance_idx = self._find_tau_higgs_resonance_ancestors_from_lepton(
            genparts,
            idx,
        )

        if tau_idx < 0 or higgs_idx < 0:
            return False, -1, -1, -1

        if self.resonance_pdgid is not None and resonance_idx < 0:
            return False, -1, -1, -1

        return True, tau_idx, higgs_idx, resonance_idx

    def _delta_phi(self, phi1, phi2):
        dphi = phi1 - phi2
        while dphi > math.pi:
            dphi -= 2.0 * math.pi
        while dphi <= -math.pi:
            dphi += 2.0 * math.pi
        return dphi

    def _delta_r(self, eta1, phi1, eta2, phi2):
        deta = eta1 - eta2
        dphi = self._delta_phi(phi1, phi2)
        return math.sqrt(deta * deta + dphi * dphi)

    def _find_nearest_reco_object(self, reco_coll, gen_eta, gen_phi, max_dr, used_indices=None):
        best_idx = -1
        best_dr = 999.0

        if used_indices is None:
            used_indices = set()

        for i, obj in enumerate(reco_coll):
            if i in used_indices:
                continue
            
            dr = self._delta_r(
                float(gen_eta),
                float(gen_phi),
                float(getattr(obj, "eta", 999.0)),
                float(getattr(obj, "phi", 999.0)),
            )

            if dr < best_dr:
                best_dr = dr
                best_idx = i

        if best_idx >= 0 and best_dr < max_dr:
            return best_idx, best_dr

        return -1, 999.0

    def _find_nearest_reco_object_no_dr_cut(self, reco_coll, ref_eta, ref_phi):
        best_idx = -1
        best_dr = 999.0

        for i, obj in enumerate(reco_coll):
            dr = self._delta_r(
                float(ref_eta),
                float(ref_phi),
                float(getattr(obj, "eta", 999.0)),
                float(getattr(obj, "phi", 999.0)),
            )

            if dr < best_dr:
                best_dr = dr
                best_idx = i

        return best_idx, best_dr

    def _find_nearest_reco_object_from_indices(self, reco_coll, allowed_indices, ref_eta, ref_phi):
        """
        Find nearest reco object to ref_eta/ref_phi, but only among allowed reco indices.

        This is useful for finding the nearest GEN-matched RECO tau to an electron.
        """
        best_idx = -1
        best_dr = 999.0

        for i in allowed_indices:
            if not self._valid_index(i, reco_coll):
                continue

            obj = reco_coll[i]

            dr = self._delta_r(
                float(ref_eta),
                float(ref_phi),
                float(getattr(obj, "eta", 999.0)),
                float(getattr(obj, "phi", 999.0)),
            )

            if dr < best_dr:
                best_dr = dr
                best_idx = i

        return best_idx, best_dr

    def _find_matched_reco_electron(self, electrons, gen_idx):
        for i, ele in enumerate(electrons):
            if int(getattr(ele, "genPartIdx", -1)) == int(gen_idx):
                return i
        return -1

    def _find_matched_reco_muon(self, muons, gen_idx):
        for i, mu in enumerate(muons):
            if int(getattr(mu, "genPartIdx", -1)) == int(gen_idx):
                return i
        return -1

    def _electron_tau_higgs_resonance_ancestors(self, ele, genparts):
        gen_idx = int(getattr(ele, "genPartIdx", -1))
        if gen_idx < 0:
            return -1, -1, -1

        if int(getattr(ele, "genPartFlav", 0)) != 15:
            return -1, -1, -1

        tau_idx, higgs_idx, resonance_idx = self._find_tau_higgs_resonance_ancestors_from_lepton(
            genparts,
            gen_idx,
        )

        if self.resonance_pdgid is not None and resonance_idx < 0:
            return tau_idx, -1, -1

        return tau_idx, higgs_idx, resonance_idx

    def _muon_tau_higgs_resonance_ancestors(self, mu, genparts):
        gen_idx = int(getattr(mu, "genPartIdx", -1))
        if gen_idx < 0:
            return -1, -1, -1

        if int(getattr(mu, "genPartFlav", 0)) != 15:
            return -1, -1, -1

        tau_idx, higgs_idx, resonance_idx = self._find_tau_higgs_resonance_ancestors_from_lepton(
            genparts,
            gen_idx,
        )

        if self.resonance_pdgid is not None and resonance_idx < 0:
            return tau_idx, -1, -1

        return tau_idx, higgs_idx, resonance_idx

    def _get_hps_tau_deeptau_vsjet(self, tau):
        for name in [
            "idDeepTau2018v2p5VSjet"
        ]:
            if hasattr(tau, name):
                return int(getattr(tau, name))
        return 0

    def _get_hps_tau_deeptau_vse(self, tau):
        for name in [
            "idDeepTau2018v2p5VSe"
        ]:
            if hasattr(tau, name):
                return int(getattr(tau, name))
        return 0

    def _get_hps_tau_deeptau_vsmu(self, tau):
        for name in [
            "idDeepTau2018v2p5VSmu"
        ]:
            if hasattr(tau, name):
                return int(getattr(tau, name))
        return 0

    def _get_boosted_tau_deeptau_vsjet(self, btau):
        for name in [
            "rawBoostedDeepTauRunIIv2p0VSjet",
        ]:
            if hasattr(btau, name):
                return float(getattr(btau, name))
        return 0

    # ----------------------------------------------------------------------
    # Main event loop
    # ----------------------------------------------------------------------
    def analyze(self, event):
        electrons = Collection(event, "Electron")
        muons = Collection(event, "Muon")
        taus = Collection(event, "Tau")
        boosted_taus = Collection(event, "boostedTau")

        genparts = Collection(event, "GenPart") if self.isMC else []
        genvistau = Collection(event, "GenVisTau") if self.isMC else []

        raw = self.current_summary["raw_counts_no_fiducial_cuts"]
        self.current_summary["events_processed"] += 1

        event_raw_reco_ele_to_gen_count = {}
        event_duplicate_ele_attempts = 0


        genElectron_pt = []
        genElectron_eta = []
        genElectron_phi = []
        genElectron_mass = []
        genElectron_charge = []
        genElectron_genPartIdx = []
        genElectron_tauAncestorIdx = []
        genElectron_higgsAncestorIdx = []
        genElectron_resonanceAncestorIdx = []

        genElectron_matchedRecoElectronIdx = []
        genElectron_hasMatchedRecoElectron = []
        genElectron_matchedRecoElectron_pt = []
        genElectron_matchedRecoElectron_eta = []
        genElectron_matchedRecoElectron_phi = []
        genElectron_matchedRecoElectron_superclusterEta = []
        genElectron_matchedRecoElectron_vidNestedWPBitmap = []
        genElectron_matchedRecoElectron_cutBased = []
        genElectron_matchedRecoElectron_vidNestedWPBitmapHEEP	 = []
        genElectron_matchedRecoElectron_cutBased_HEEP = []
        genElectron_matchedRecoElectron_genPartFlav = []
        genElectron_matchedRecoElectron_convVeto = []
        genElectron_matchedRecoElectron_pfRelIso03_all = []
        genElectron_matchedRecoElectron_pfRelIso03_chg = []
        genElectron_matchedRecoElectron_pfRelIso04_all = []
        genElectron_matchedRecoElectron_miniPFRelIso_all = []
        genElectron_matchedRecoElectron_miniPFRelIso_chg = []
        genElectron_matchedRecoElectron_sieie = []
        genElectron_matchedRecoElectron_deltaEtaSC = []
        genElectron_matchedRecoElectron_eInvMinusPInv = []
        genElectron_matchedRecoElectron_hoe = []
        genElectron_matchedRecoElectron_lostHits = []
        genElectron_matchedRecoElectron_rawEnergy = []



        genElectron_matchedRecoElectron_deltaR = []
        genElectron_matchedRecoElectron_nearestRecoTau_deltaR = []
        genElectron_matchedRecoElectron_nearestRecoTauIdx = []
        genElectron_matchedRecoElectron_nearestRecoTau_pt = []
        genElectron_matchedRecoElectron_nearestRecoTau_decayMode = []
        genElectron_matchedRecoElectron_nearestRecoBoostedTau_deltaR = []
        genElectron_matchedRecoElectron_nearestRecoBoostedTauIdx = []
        genElectron_matchedRecoElectron_nearestRecoBoostedTau_pt = []
        genElectron_matchedRecoElectron_nearestRecoBoostedTau_decayMode = []

        genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_deltaR = []
        genElectron_matchedRecoElectron_nearestGenMatchedRecoTauIdx = []
        genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_pt = []
        genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_eta = []
        genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_phi = []
        genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_decayMode = []
        genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_idDecayModeNewDMs = []
        genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_idDeepTau2018v2p5VSjet = []
        genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_idDeepTau2018v2p5VSe = []
        genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_idDeepTau2018v2p5VSmu = []

        genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_deltaR = []
        genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTauIdx = []
        genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_pt = []
        genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_eta = []
        genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_phi = []
        genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_decayMode = []
        genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_idAntiEle2018 = []
        genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_idAntiMu = []
        genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_rawBoostedDeepTauRunIIv2p0VSjet = []
        

        # GEN muon containers
        genMuon_pt = []
        genMuon_eta = []
        genMuon_phi = []
        genMuon_mass = []
        genMuon_charge = []
        genMuon_genPartIdx = []
        genMuon_tauAncestorIdx = []
        genMuon_higgsAncestorIdx = []
        genMuon_resonanceAncestorIdx = []

        genMuon_matchedRecoMuonIdx = []
        genMuon_hasMatchedRecoMuon = []
        genMuon_matchedRecoMuon_pt = []
        genMuon_matchedRecoMuon_eta = []
        genMuon_matchedRecoMuon_phi = []
        genMuon_matchedRecoMuon_looseId = []
        genMuon_matchedRecoMuon_mediumId = []
        genMuon_matchedRecoMuon_tightId = []
        genMuon_matchedRecoMuon_genPartFlav = []


        # GEN visible tau containers
        genVisTau_pt = []
        genVisTau_eta = []
        genVisTau_phi = []
        genVisTau_mass = []
        genVisTau_status = []
        genVisTau_genVisTauIdx = []
        genVisTau_tauAncestorIdx = []
        genVisTau_higgsAncestorIdx = []
        genVisTau_resonanceAncestorIdx = []

        genVisTau_matchedRecoTauIdx = []
        genVisTau_hasMatchedRecoTau = []
        genVisTau_matchedRecoTau_deltaR = []
        genVisTau_matchedRecoTau_pt = []
        genVisTau_matchedRecoTau_dz = []

        genVisTau_matchedRecoTau_eta = []
        genVisTau_matchedRecoTau_phi = []
        genVisTau_matchedRecoTau_mass = []
        genVisTau_matchedRecoTau_decayMode = []
        genVisTau_matchedRecoTau_idDecayModeNewDMs = []
        genVisTau_matchedRecoTau_idDeepTau2018v2p5VSjet = []
        genVisTau_matchedRecoTau_idDeepTau2018v2p5VSe = []
        genVisTau_matchedRecoTau_idDeepTau2018v2p5VSmu = []

        genVisTau_matchedRecoBoostedTauIdx = []
        genVisTau_hasMatchedRecoBoostedTau = []
        genVisTau_matchedRecoBoostedTau_deltaR = []
        genVisTau_matchedRecoBoostedTau_pt = []
        genVisTau_matchedRecoBoostedTau_eta = []
        genVisTau_matchedRecoBoostedTau_phi = []
        genVisTau_matchedRecoBoostedTau_mass = []
        genVisTau_matchedRecoBoostedTau_decayMode = []
        genVisTau_matchedRecoBoostedTau_idAntiEle2018 = []
        genVisTau_matchedRecoBoostedTau_idAntiMu = []
        genVisTau_matchedRecoBoostedTau_idMVAnewDM2017v2 = []
        genVisTau_matchedRecoBoostedTau_rawBoostedDeepTauRunIIv2p0VSjet = []

        if self.isMC:
            nGenMatchedRecoElectronFromHiggsTau = 0
            nGenMatchedRecoMuonFromHiggsTau = 0

            # GenVisTau from Higgs -> tau, then DeltaR match to Tau/boostedTau
            used_reco_tau_indices = set()
            used_reco_boosted_tau_indices = set()

            genmatched_reco_tau_indices = []
            genmatched_reco_boosted_tau_indices = []
            
            for ivis, vis in enumerate(genvistau):

                tau_idx = int(getattr(vis, "genPartIdxMother", -1))

                if not self._valid_index(tau_idx, genparts):
                    continue

                if abs(int(genparts[tau_idx].pdgId)) != 15:
                    continue

                higgs_idx, resonance_idx = self._find_higgs_resonance_from_tau(
                    genparts,
                    tau_idx,
                )

                if higgs_idx < 0:
                    continue

                raw["gen_particles_from_higgs_tau"]["visible_hadronic_taus"] += 1

                # HPS tau match
                reco_tau_idx, reco_tau_dr = self._find_nearest_reco_object(
                    taus,
                    float(vis.eta),
                    float(vis.phi),
                    self.max_dr_hps_tau,
                    used_indices=used_reco_tau_indices,

                )

                if reco_tau_idx >= 0:
                    used_reco_tau_indices.add(reco_tau_idx)
                    genmatched_reco_tau_indices.append(reco_tau_idx)

                has_reco_tau = int(reco_tau_idx >= 0)
                if has_reco_tau:
                    raw["gen_particles_with_matched_reco"]["visible_taus_to_hps_taus"] += 1
                    reco_tau = taus[reco_tau_idx]
                else:
                    reco_tau = None

                # boostedTau match
                reco_btau_idx, reco_btau_dr = self._find_nearest_reco_object(
                    boosted_taus,
                    float(vis.eta),
                    float(vis.phi),
                    self.max_dr_boosted_tau,
                    used_indices=used_reco_boosted_tau_indices,

                )

                
                if reco_btau_idx >= 0:
                    used_reco_boosted_tau_indices.add(reco_btau_idx)
                    genmatched_reco_boosted_tau_indices.append(reco_btau_idx)

                has_reco_btau = int(reco_btau_idx >= 0)
                if has_reco_btau:
                    raw["gen_particles_with_matched_reco"]["visible_taus_to_boosted_taus"] += 1
                    reco_btau = boosted_taus[reco_btau_idx]
                else:
                    reco_btau = None

                # Fill GenVisTau info
                genVisTau_pt.append(float(vis.pt))
                genVisTau_eta.append(float(vis.eta))
                genVisTau_phi.append(float(vis.phi))
                genVisTau_mass.append(float(vis.mass))
                genVisTau_status.append(int(getattr(vis, "status", -1)))

                genVisTau_genVisTauIdx.append(int(ivis))
                genVisTau_tauAncestorIdx.append(int(tau_idx))
                genVisTau_higgsAncestorIdx.append(int(higgs_idx))
                genVisTau_resonanceAncestorIdx.append(int(resonance_idx))

                # Fill matched HPS Tau info
                genVisTau_matchedRecoTauIdx.append(int(reco_tau_idx))
                genVisTau_hasMatchedRecoTau.append(has_reco_tau)
                genVisTau_matchedRecoTau_deltaR.append(float(reco_tau_dr))

                if reco_tau is not None:
                    genVisTau_matchedRecoTau_pt.append(float(getattr(reco_tau, "pt", -999.0)))
                    genVisTau_matchedRecoTau_eta.append(float(getattr(reco_tau, "eta", -999.0)))
                    genVisTau_matchedRecoTau_phi.append(float(getattr(reco_tau, "phi", -999.0)))
                    genVisTau_matchedRecoTau_dz.append(float(getattr(reco_tau, "dz", -999.0)))

                    genVisTau_matchedRecoTau_mass.append(float(getattr(reco_tau, "mass", -999.0)))
                    genVisTau_matchedRecoTau_decayMode.append(int(getattr(reco_tau, "decayMode", -1)))
                    genVisTau_matchedRecoTau_idDecayModeNewDMs.append(int(getattr(reco_tau, "idDecayModeNewDMs", 0)))
                    genVisTau_matchedRecoTau_idDeepTau2018v2p5VSjet.append(self._get_hps_tau_deeptau_vsjet(reco_tau))
                    genVisTau_matchedRecoTau_idDeepTau2018v2p5VSe.append(self._get_hps_tau_deeptau_vse(reco_tau))
                    genVisTau_matchedRecoTau_idDeepTau2018v2p5VSmu.append(self._get_hps_tau_deeptau_vsmu(reco_tau))
                else:
                    genVisTau_matchedRecoTau_pt.append(float(-999.0))
                    genVisTau_matchedRecoTau_eta.append(float(-999.0))
                    genVisTau_matchedRecoTau_phi.append(float(-999.0))
                    genVisTau_matchedRecoTau_dz.append(float(-999.0))
                    genVisTau_matchedRecoTau_mass.append(float(-999.0))
                    genVisTau_matchedRecoTau_decayMode.append(int(-1))
                    genVisTau_matchedRecoTau_idDecayModeNewDMs.append(int(0))
                    genVisTau_matchedRecoTau_idDeepTau2018v2p5VSjet.append(int(0))
                    genVisTau_matchedRecoTau_idDeepTau2018v2p5VSe.append(int(0))
                    genVisTau_matchedRecoTau_idDeepTau2018v2p5VSmu.append(int(0))

                # Fill matched boostedTau info
                genVisTau_matchedRecoBoostedTauIdx.append(int(reco_btau_idx))
                genVisTau_hasMatchedRecoBoostedTau.append(has_reco_btau)
                genVisTau_matchedRecoBoostedTau_deltaR.append(float(reco_btau_dr))

                if reco_btau is not None:
                    anti_ele = int(getattr(reco_btau, "idAntiEle2018", 0))
                    anti_mu = int(getattr(reco_btau, "idAntiMu", 0))


                    genVisTau_matchedRecoBoostedTau_pt.append(float(getattr(reco_btau, "pt", -999.0)))
                    genVisTau_matchedRecoBoostedTau_eta.append(float(getattr(reco_btau, "eta", -999.0)))
                    genVisTau_matchedRecoBoostedTau_phi.append(float(getattr(reco_btau, "phi", -999.0)))
                    genVisTau_matchedRecoBoostedTau_mass.append(float(getattr(reco_btau, "mass", -999.0)))
                    genVisTau_matchedRecoBoostedTau_decayMode.append(int(getattr(reco_btau, "decayMode", -1)))
                    genVisTau_matchedRecoBoostedTau_idAntiEle2018.append(anti_ele)
                    genVisTau_matchedRecoBoostedTau_idAntiMu.append(anti_mu)
                    genVisTau_matchedRecoBoostedTau_idMVAnewDM2017v2.append(int(getattr(reco_btau, "idMVAnewDM2017v2", 0)))
                    genVisTau_matchedRecoBoostedTau_rawBoostedDeepTauRunIIv2p0VSjet.append(self._get_boosted_tau_deeptau_vsjet(reco_btau))
                else:
                    genVisTau_matchedRecoBoostedTau_pt.append(float(-999.0))
                    genVisTau_matchedRecoBoostedTau_eta.append(float(-999.0))
                    genVisTau_matchedRecoBoostedTau_phi.append(float(-999.0))
                    genVisTau_matchedRecoBoostedTau_mass.append(float(-999.0))
                    genVisTau_matchedRecoBoostedTau_decayMode.append(int(-1))
                    genVisTau_matchedRecoBoostedTau_idAntiEle2018.append(int(0))
                    genVisTau_matchedRecoBoostedTau_idAntiMu.append(int(0))
                    genVisTau_matchedRecoBoostedTau_idMVAnewDM2017v2.append(int(0))
                    genVisTau_matchedRecoBoostedTau_rawBoostedDeepTauRunIIv2p0VSjet.append(float(-999.0))



            used_reco_electron_indices = set()
            used_reco_muon_indices = set()
            for igen, gp in enumerate(genparts):

                # GEN electron from Higgs -> tau -> electron
                is_ele_target, tau_idx, higgs_idx, resonance_idx = (
                    self._is_stable_gen_lepton_from_higgs_tau(
                        genparts,
                        igen,
                        abs_lepton_pdgid=11,
                    )
                )

                if is_ele_target:
                    raw["gen_particles_from_higgs_tau"]["electrons"] += 1

                   
                    genElectron_pt.append(float(gp.pt))
                    genElectron_eta.append(float(gp.eta))
                    genElectron_phi.append(float(gp.phi))
                    genElectron_mass.append(float(gp.mass))
                    genElectron_charge.append(int(-1 if int(gp.pdgId) == 11 else 1))

                    genElectron_genPartIdx.append(int(igen))
                    genElectron_tauAncestorIdx.append(int(tau_idx))
                    genElectron_higgsAncestorIdx.append(int(higgs_idx))
                    genElectron_resonanceAncestorIdx.append(int(resonance_idx))

                    
                    # Diagnostic only:
                    # What would this GEN electron match to if we did not enforce unique RECO usage?
                    raw_reco_idx, raw_reco_dr = self._find_nearest_reco_object(
                        electrons,
                        float(gp.eta),
                        float(gp.phi),
                        self.max_dr_ele,
                        used_indices=None,
                    )

                    if raw_reco_idx >= 0:
                        event_raw_reco_ele_to_gen_count[raw_reco_idx] = (event_raw_reco_ele_to_gen_count.get(raw_reco_idx, 0) + 1)

                        if event_raw_reco_ele_to_gen_count[raw_reco_idx] > 1:
                            event_duplicate_ele_attempts += 1

                    reco_idx, reco_ele_dr = self._find_nearest_reco_object(
                        electrons,
                        float(gp.eta),
                        float(gp.phi),
                        self.max_dr_ele,
                        used_indices=used_reco_electron_indices
                        )
                    genElectron_matchedRecoElectron_deltaR.append(float(reco_ele_dr))
                    
                    if reco_idx >= 0:
                        nGenMatchedRecoElectronFromHiggsTau += 1
                        used_reco_electron_indices.add(reco_idx)
                        ele = electrons[reco_idx]
                        cut_based = int(getattr(ele, "cutBased", 0))
                        cutBased_HEEP = int(getattr(ele, "cutBased_HEEP", 0))

                        raw["gen_particles_with_matched_reco"]["electrons_to_reco_electrons"] += 1
    
                        genElectron_matchedRecoElectronIdx.append(int(reco_idx))
                        genElectron_hasMatchedRecoElectron.append(int(reco_idx >= 0))
                        genElectron_matchedRecoElectron_pt.append(float(getattr(ele, "pt", -999.0)))
                        genElectron_matchedRecoElectron_eta.append(float(getattr(ele, "eta", -999.0)))
                        genElectron_matchedRecoElectron_phi.append(float(getattr(ele, "phi", -999.0)))
                        genElectron_matchedRecoElectron_superclusterEta.append(float(getattr(ele, "superclusterEta", -999.0)))
                        genElectron_matchedRecoElectron_vidNestedWPBitmap.append(int(getattr(ele, "vidNestedWPBitmap", 0)))
                        genElectron_matchedRecoElectron_cutBased.append(cut_based)
                        genElectron_matchedRecoElectron_vidNestedWPBitmapHEEP.append(int(getattr(ele, "vidNestedWPBitmapHEEP", 0)))
                        genElectron_matchedRecoElectron_cutBased_HEEP.append(cutBased_HEEP)
                        genElectron_matchedRecoElectron_genPartFlav.append(int(getattr(ele, "genPartFlav", 0)))


                        genElectron_matchedRecoElectron_convVeto.append(int(bool(getattr(ele, "convVeto", False))))
                        genElectron_matchedRecoElectron_pfRelIso03_all.append(float(getattr(ele, "pfRelIso03_all", -999.0)))
                        genElectron_matchedRecoElectron_pfRelIso03_chg.append(float(getattr(ele, "pfRelIso03_chg", -999.0)))
                        genElectron_matchedRecoElectron_pfRelIso04_all.append(float(getattr(ele, "pfRelIso04_all", -999.0)))
                        genElectron_matchedRecoElectron_miniPFRelIso_all.append(float(getattr(ele, "miniPFRelIso_all", -999.0)))
                        genElectron_matchedRecoElectron_miniPFRelIso_chg.append(float(getattr(ele, "miniPFRelIso_chg", -999.0)))
                        genElectron_matchedRecoElectron_sieie.append(float(getattr(ele, "sieie", -999.0)))
                        genElectron_matchedRecoElectron_deltaEtaSC.append(float(getattr(ele, "deltaEtaSC", -999.0)))
                        genElectron_matchedRecoElectron_eInvMinusPInv.append(float(getattr(ele, "eInvMinusPInv", -999.0)))
                        genElectron_matchedRecoElectron_hoe.append(float(getattr(ele, "hoe", -999.0)))
                        genElectron_matchedRecoElectron_lostHits.append(int(getattr(ele, "lostHits", -1)))
                        genElectron_matchedRecoElectron_rawEnergy.append(float(getattr(ele, "rawEnergy", -999)))

                        
                        ### Nearest Taus are not GEN matched

                        nearest_tau_idx, nearest_tau_dr = self._find_nearest_reco_object_no_dr_cut(taus, float(ele.eta), float(ele.phi))
                        genElectron_matchedRecoElectron_nearestRecoTauIdx.append(int(nearest_tau_idx))
                        genElectron_matchedRecoElectron_nearestRecoTau_deltaR.append(float(nearest_tau_dr))

                        if nearest_tau_idx >= 0:
                            nearest_tau = taus[nearest_tau_idx]
                            genElectron_matchedRecoElectron_nearestRecoTau_pt.append(float(getattr(nearest_tau, "pt", -999.0)))
                            genElectron_matchedRecoElectron_nearestRecoTau_decayMode.append(int(getattr(nearest_tau, "decayMode", -1)))
                        else:
                            genElectron_matchedRecoElectron_nearestRecoTau_pt.append(float(-999.0))
                            genElectron_matchedRecoElectron_nearestRecoTau_decayMode.append(int(-1))

                        nearest_btau_idx, nearest_btau_dr = self._find_nearest_reco_object_no_dr_cut(
                            boosted_taus,
                            float(ele.eta),
                            float(ele.phi),
                        )

                        genElectron_matchedRecoElectron_nearestRecoBoostedTauIdx.append(int(nearest_btau_idx))
                        genElectron_matchedRecoElectron_nearestRecoBoostedTau_deltaR.append(float(nearest_btau_dr))

                        if nearest_btau_idx >= 0:
                            nearest_btau = boosted_taus[nearest_btau_idx]
                            genElectron_matchedRecoElectron_nearestRecoBoostedTau_pt.append(
                                float(getattr(nearest_btau, "pt", -999.0))
                            )
                            genElectron_matchedRecoElectron_nearestRecoBoostedTau_decayMode.append(
                                int(getattr(nearest_btau, "decayMode", -1))
                            )
                        else:
                            genElectron_matchedRecoElectron_nearestRecoBoostedTau_pt.append(float(-999.0))
                            genElectron_matchedRecoElectron_nearestRecoBoostedTau_decayMode.append(int(-1))


                        ########### Nearest Taus are GEN matched

                        nearest_gm_tau_idx, nearest_gm_tau_dr = (self._find_nearest_reco_object_from_indices(taus, genmatched_reco_tau_indices, float(ele.eta), float(ele.phi)))
                        genElectron_matchedRecoElectron_nearestGenMatchedRecoTauIdx.append(int(nearest_gm_tau_idx))
                        genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_deltaR.append(float(nearest_gm_tau_dr))

                        if nearest_gm_tau_idx >= 0:
                            gm_tau = taus[nearest_gm_tau_idx]

                            genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_pt.append(float(getattr(gm_tau, "pt", -999.0)))
                            genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_eta.append(float(getattr(gm_tau, "eta", -999.0)))
                            genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_phi.append(float(getattr(gm_tau, "phi", -999.0)))
                            genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_decayMode.append(int(getattr(gm_tau, "decayMode", -1)))
                            genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_idDecayModeNewDMs.append(int(getattr(gm_tau, "idDecayModeNewDMs", 0)))
                            genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_idDeepTau2018v2p5VSjet.append(self._get_hps_tau_deeptau_vsjet(gm_tau))
                            genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_idDeepTau2018v2p5VSe.append(self._get_hps_tau_deeptau_vse(gm_tau))
                            genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_idDeepTau2018v2p5VSmu.append(self._get_hps_tau_deeptau_vsmu(gm_tau))
                        else:
                            genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_pt.append(float(-999.0))
                            genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_eta.append(float(-999.0))
                            genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_phi.append(float(-999.0))
                            genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_decayMode.append(int(-1))
                            genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_idDecayModeNewDMs.append(int(0))
                            genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_idDeepTau2018v2p5VSjet.append(int(0))
                            genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_idDeepTau2018v2p5VSe.append(int(0))
                            genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_idDeepTau2018v2p5VSmu.append(int(0))

                        nearest_gm_btau_idx, nearest_gm_btau_dr = (self._find_nearest_reco_object_from_indices(boosted_taus, genmatched_reco_boosted_tau_indices, float(ele.eta), float(ele.phi)))

                        genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTauIdx.append(int(nearest_gm_btau_idx))
                        genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_deltaR.append(float(nearest_gm_btau_dr))

                        if nearest_gm_btau_idx >= 0:
                            gm_btau = boosted_taus[nearest_gm_btau_idx]

                            genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_pt.append(float(getattr(gm_btau, "pt", -999.0)))
                            genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_eta.append(float(getattr(gm_btau, "eta", -999.0)))
                            genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_phi.append(float(getattr(gm_btau, "phi", -999.0)))
                            genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_decayMode.append(int(getattr(gm_btau, "decayMode", -1)))
                            genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_idAntiEle2018.append(int(getattr(gm_btau, "idAntiEle2018", 0)))
                            genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_idAntiMu.append(int(getattr(gm_btau, "idAntiMu", 0)))
                            genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_rawBoostedDeepTauRunIIv2p0VSjet.append(self._get_boosted_tau_deeptau_vsjet(gm_btau))
                        else:
                            genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_pt.append(float(-999.0))
                            genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_eta.append(float(-999.0))
                            genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_phi.append(float(-999.0))
                            genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_decayMode.append(int(-1))
                            genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_idAntiEle2018.append(int(0))
                            genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_idAntiMu.append(int(0))
                            genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_rawBoostedDeepTauRunIIv2p0VSjet.append(float(-999.0))


                    else:
                        genElectron_matchedRecoElectron_pt.append(float(-999.0))
                        genElectron_matchedRecoElectron_eta.append(float(-999.0))
                        genElectron_matchedRecoElectron_phi.append(float(-999.0))
                        genElectron_matchedRecoElectron_superclusterEta.append(float(-999.0))
                        genElectron_matchedRecoElectron_vidNestedWPBitmap.append(int(0))
                        genElectron_matchedRecoElectron_cutBased.append(int(0))
                        genElectron_matchedRecoElectron_vidNestedWPBitmapHEEP.append(int(0))
                        genElectron_matchedRecoElectron_cutBased_HEEP.append(int(0))
                        genElectron_matchedRecoElectron_genPartFlav.append(int(0))
                        genElectron_matchedRecoElectronIdx.append(int(-1))
                        genElectron_hasMatchedRecoElectron.append(int(0))
                        genElectron_matchedRecoElectron_convVeto.append(int(0))
                        genElectron_matchedRecoElectron_pfRelIso03_all.append(float(-999.0))
                        genElectron_matchedRecoElectron_pfRelIso03_chg.append(float(-999.0))
                        genElectron_matchedRecoElectron_pfRelIso04_all.append(float(-999.0))
                        genElectron_matchedRecoElectron_miniPFRelIso_all.append(float(-999.0))
                        genElectron_matchedRecoElectron_miniPFRelIso_chg.append(float(-999.0))
                        genElectron_matchedRecoElectron_sieie.append(float(-999.0))
                        genElectron_matchedRecoElectron_deltaEtaSC.append(float(-999.0))
                        genElectron_matchedRecoElectron_eInvMinusPInv.append(float(-999.0))
                        genElectron_matchedRecoElectron_hoe.append(float(-999.0))
                        genElectron_matchedRecoElectron_lostHits.append(int(-1))
                        genElectron_matchedRecoElectron_rawEnergy.append(float(-999.0))


                        genElectron_matchedRecoElectron_nearestRecoTauIdx.append(int(-1))
                        genElectron_matchedRecoElectron_nearestRecoTau_deltaR.append(float(-999.0))
                        genElectron_matchedRecoElectron_nearestRecoTau_pt.append(float(-999.0))
                        genElectron_matchedRecoElectron_nearestRecoTau_decayMode.append(int(-1))

                        genElectron_matchedRecoElectron_nearestRecoBoostedTauIdx.append(int(-1))
                        genElectron_matchedRecoElectron_nearestRecoBoostedTau_deltaR.append(float(-999.0))
                        genElectron_matchedRecoElectron_nearestRecoBoostedTau_pt.append(float(-999.0))
                        genElectron_matchedRecoElectron_nearestRecoBoostedTau_decayMode.append(int(-1))
                        
                        genElectron_matchedRecoElectron_nearestGenMatchedRecoTauIdx.append(int(-1))
                        genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_deltaR.append(float(999.0))
                        genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_pt.append(float(-999.0))
                        genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_eta.append(float(-999.0))
                        genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_phi.append(float(-999.0))
                        genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_decayMode.append(int(-1))
                        genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_idDecayModeNewDMs.append(int(0))
                        genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_idDeepTau2018v2p5VSjet.append(int(0))
                        genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_idDeepTau2018v2p5VSe.append(int(0))
                        genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_idDeepTau2018v2p5VSmu.append(int(0))
                        
                        genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTauIdx.append(int(-1))
                        genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_deltaR.append(float(999.0))
                        genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_pt.append(float(-999.0))
                        genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_eta.append(float(-999.0))
                        genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_phi.append(float(-999.0))
                        genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_decayMode.append(int(-1))
                        genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_idAntiEle2018.append(int(0))
                        genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_idAntiMu.append(int(0))
                        genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_rawBoostedDeepTauRunIIv2p0VSjet.append(float(-999.0))



                # GEN muon from Higgs -> tau -> muon
                is_mu_target, tau_idx, higgs_idx, resonance_idx = (
                    self._is_stable_gen_lepton_from_higgs_tau(
                        genparts,
                        igen,
                        abs_lepton_pdgid=13,
                    )
                )

                if is_mu_target:
                    raw["gen_particles_from_higgs_tau"]["muons"] += 1

                    reco_idx, reco_muon_dr = self._find_nearest_reco_object(
                        muons,
                        float(gp.eta),
                        float(gp.phi),
                        self.max_dr_muon,
                        used_indices=used_reco_muon_indices
                        )

                    genMuon_pt.append(float(gp.pt))
                    genMuon_eta.append(float(gp.eta))
                    genMuon_phi.append(float(gp.phi))
                    genMuon_mass.append(float(gp.mass))
                    genMuon_charge.append(int(-1 if int(gp.pdgId) == 13 else 1))

                    genMuon_genPartIdx.append(int(igen))
                    genMuon_tauAncestorIdx.append(int(tau_idx))
                    genMuon_higgsAncestorIdx.append(int(higgs_idx))
                    genMuon_resonanceAncestorIdx.append(int(resonance_idx))

                    genMuon_matchedRecoMuonIdx.append(int(reco_idx))
                    genMuon_hasMatchedRecoMuon.append(int(reco_idx >= 0))

                    if reco_idx >= 0:
                        used_reco_muon_indices.add(reco_idx)
                        nGenMatchedRecoMuonFromHiggsTau += 1
                        mu = muons[reco_idx]
                        loose_id = int(bool(getattr(mu, "looseId", False)))
                        medium_id = int(bool(getattr(mu, "mediumId", False)))
                        tight_id = int(bool(getattr(mu, "tightId", False)))
                        raw["gen_particles_with_matched_reco"]["muons_to_reco_muons"] += 1

                        genMuon_matchedRecoMuon_pt.append(float(getattr(mu, "pt", -999.0)))
                        genMuon_matchedRecoMuon_eta.append(float(getattr(mu, "eta", -999.0)))
                        genMuon_matchedRecoMuon_phi.append(float(getattr(mu, "phi", -999.0)))
                        genMuon_matchedRecoMuon_looseId.append(loose_id)
                        genMuon_matchedRecoMuon_mediumId.append(medium_id)
                        genMuon_matchedRecoMuon_tightId.append(tight_id)
                        genMuon_matchedRecoMuon_genPartFlav.append(int(getattr(mu, "genPartFlav", 0)))
                    else:
                        genMuon_matchedRecoMuon_pt.append(float(-999.0))
                        genMuon_matchedRecoMuon_eta.append(float(-999.0))
                        genMuon_matchedRecoMuon_phi.append(float(-999.0))
                        genMuon_matchedRecoMuon_looseId.append(int(0))
                        genMuon_matchedRecoMuon_mediumId.append(int(0))
                        genMuon_matchedRecoMuon_tightId.append(int(0))
                        genMuon_matchedRecoMuon_genPartFlav.append(int(0))

            if event_duplicate_ele_attempts > 0:
                raw["gen_particles_with_matched_reco"]["events_with_electron_duplicate_match_attempts"] += 1

                raw["gen_particles_with_matched_reco"]["electron_duplicate_match_attempts_before_unique_assignment"] += event_duplicate_ele_attempts


        
        self.out.fillBranch("nGenVisTauFromHiggsTau", len(genVisTau_pt))
        self.out.fillBranch("GenVisTauFromHiggsTau_pt", genVisTau_pt)
        self.out.fillBranch("GenVisTauFromHiggsTau_eta", genVisTau_eta)
        self.out.fillBranch("GenVisTauFromHiggsTau_phi", genVisTau_phi)
        self.out.fillBranch("GenVisTauFromHiggsTau_mass", genVisTau_mass)
        self.out.fillBranch("GenVisTauFromHiggsTau_status", genVisTau_status)
        self.out.fillBranch("GenVisTauFromHiggsTau_genVisTauIdx", genVisTau_genVisTauIdx)
        self.out.fillBranch("GenVisTauFromHiggsTau_tauAncestorIdx", genVisTau_tauAncestorIdx)
        self.out.fillBranch("GenVisTauFromHiggsTau_higgsAncestorIdx", genVisTau_higgsAncestorIdx)
        self.out.fillBranch("GenVisTauFromHiggsTau_resonanceAncestorIdx", genVisTau_resonanceAncestorIdx)
        self.out.fillBranch("GenVisTauFromHiggsTau_matchedRecoTauIdx", genVisTau_matchedRecoTauIdx)
        self.out.fillBranch("GenVisTauFromHiggsTau_hasMatchedRecoTau", genVisTau_hasMatchedRecoTau)
        self.out.fillBranch("GenVisTauFromHiggsTau_matchedRecoTau_deltaR", genVisTau_matchedRecoTau_deltaR)
        self.out.fillBranch("GenVisTauFromHiggsTau_matchedRecoTau_pt", genVisTau_matchedRecoTau_pt)
        self.out.fillBranch("GenVisTauFromHiggsTau_matchedRecoTau_eta", genVisTau_matchedRecoTau_eta)
        self.out.fillBranch("GenVisTauFromHiggsTau_matchedRecoTau_phi", genVisTau_matchedRecoTau_phi)
        self.out.fillBranch("GenVisTauFromHiggsTau_matchedRecoTau_dz", genVisTau_matchedRecoTau_dz)

        self.out.fillBranch("GenVisTauFromHiggsTau_matchedRecoTau_mass", genVisTau_matchedRecoTau_mass)
        self.out.fillBranch("GenVisTauFromHiggsTau_matchedRecoTau_decayMode", genVisTau_matchedRecoTau_decayMode)
        self.out.fillBranch("GenVisTauFromHiggsTau_matchedRecoTau_idDecayModeNewDMs", genVisTau_matchedRecoTau_idDecayModeNewDMs)
        self.out.fillBranch("GenVisTauFromHiggsTau_matchedRecoTau_idDeepTau2018v2p5VSjet", genVisTau_matchedRecoTau_idDeepTau2018v2p5VSjet)
        self.out.fillBranch("GenVisTauFromHiggsTau_matchedRecoTau_idDeepTau2018v2p5VSe", genVisTau_matchedRecoTau_idDeepTau2018v2p5VSe)
        self.out.fillBranch("GenVisTauFromHiggsTau_matchedRecoTau_idDeepTau2018v2p5VSmu", genVisTau_matchedRecoTau_idDeepTau2018v2p5VSmu)

        self.out.fillBranch("GenVisTauFromHiggsTau_matchedRecoBoostedTauIdx", genVisTau_matchedRecoBoostedTauIdx)
        self.out.fillBranch("GenVisTauFromHiggsTau_hasMatchedRecoBoostedTau", genVisTau_hasMatchedRecoBoostedTau)
        self.out.fillBranch("GenVisTauFromHiggsTau_matchedRecoBoostedTau_deltaR", genVisTau_matchedRecoBoostedTau_deltaR)
        self.out.fillBranch("GenVisTauFromHiggsTau_matchedRecoBoostedTau_pt", genVisTau_matchedRecoBoostedTau_pt)
        self.out.fillBranch("GenVisTauFromHiggsTau_matchedRecoBoostedTau_eta", genVisTau_matchedRecoBoostedTau_eta)
        self.out.fillBranch("GenVisTauFromHiggsTau_matchedRecoBoostedTau_phi", genVisTau_matchedRecoBoostedTau_phi)
        self.out.fillBranch("GenVisTauFromHiggsTau_matchedRecoBoostedTau_mass", genVisTau_matchedRecoBoostedTau_mass)
        self.out.fillBranch("GenVisTauFromHiggsTau_matchedRecoBoostedTau_decayMode", genVisTau_matchedRecoBoostedTau_decayMode)
        self.out.fillBranch("GenVisTauFromHiggsTau_matchedRecoBoostedTau_idAntiEle2018", genVisTau_matchedRecoBoostedTau_idAntiEle2018)
        self.out.fillBranch("GenVisTauFromHiggsTau_matchedRecoBoostedTau_idAntiMu", genVisTau_matchedRecoBoostedTau_idAntiMu)
        self.out.fillBranch("GenVisTauFromHiggsTau_matchedRecoBoostedTau_idMVAnewDM2017v2", genVisTau_matchedRecoBoostedTau_idMVAnewDM2017v2)
        self.out.fillBranch("GenVisTauFromHiggsTau_matchedRecoBoostedTau_rawBoostedDeepTauRunIIv2p0VSjet", genVisTau_matchedRecoBoostedTau_rawBoostedDeepTauRunIIv2p0VSjet)
        
        self.out.fillBranch("nGenElectronFromHiggsTau", len(genElectron_pt))
        self.out.fillBranch("GenElectronFromHiggsTau_pt", genElectron_pt)
        self.out.fillBranch("GenElectronFromHiggsTau_eta", genElectron_eta)
        self.out.fillBranch("GenElectronFromHiggsTau_phi", genElectron_phi)
        self.out.fillBranch("GenElectronFromHiggsTau_mass", genElectron_mass)
        self.out.fillBranch("GenElectronFromHiggsTau_charge", genElectron_charge)
        self.out.fillBranch("GenElectronFromHiggsTau_genPartIdx", genElectron_genPartIdx)
        self.out.fillBranch("GenElectronFromHiggsTau_tauAncestorIdx", genElectron_tauAncestorIdx)
        self.out.fillBranch("GenElectronFromHiggsTau_higgsAncestorIdx", genElectron_higgsAncestorIdx)
        self.out.fillBranch("GenElectronFromHiggsTau_resonanceAncestorIdx", genElectron_resonanceAncestorIdx)

        self.out.fillBranch("nGenMatchedRecoElectronFromHiggsTau",int(nGenMatchedRecoElectronFromHiggsTau))
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectronIdx", genElectron_matchedRecoElectronIdx)
        self.out.fillBranch("GenElectronFromHiggsTau_hasMatchedRecoElectron", genElectron_hasMatchedRecoElectron)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_pt", genElectron_matchedRecoElectron_pt)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_eta", genElectron_matchedRecoElectron_eta)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_phi", genElectron_matchedRecoElectron_phi)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_superclusterEta", genElectron_matchedRecoElectron_superclusterEta)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_vidNestedWPBitmap", genElectron_matchedRecoElectron_vidNestedWPBitmap)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_cutBased", genElectron_matchedRecoElectron_cutBased)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_vidNestedWPBitmapHEEP", genElectron_matchedRecoElectron_vidNestedWPBitmapHEEP)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_cutBased_HEEP", genElectron_matchedRecoElectron_cutBased_HEEP)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_genPartFlav", genElectron_matchedRecoElectron_genPartFlav)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_convVeto", genElectron_matchedRecoElectron_convVeto)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_pfRelIso03_all", genElectron_matchedRecoElectron_pfRelIso03_all)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_pfRelIso03_chg",genElectron_matchedRecoElectron_pfRelIso03_chg)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_pfRelIso04_all",genElectron_matchedRecoElectron_pfRelIso04_all)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_miniPFRelIso_all",genElectron_matchedRecoElectron_miniPFRelIso_all)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_miniPFRelIso_chg",genElectron_matchedRecoElectron_miniPFRelIso_chg)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_sieie",genElectron_matchedRecoElectron_sieie)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_deltaEtaSC",genElectron_matchedRecoElectron_deltaEtaSC)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_eInvMinusPInv",genElectron_matchedRecoElectron_eInvMinusPInv)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_hoe",genElectron_matchedRecoElectron_hoe)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_lostHits",genElectron_matchedRecoElectron_lostHits)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_rawEnergy",genElectron_matchedRecoElectron_rawEnergy)


        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_deltaR",genElectron_matchedRecoElectron_deltaR)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_nearestRecoTau_deltaR",genElectron_matchedRecoElectron_nearestRecoTau_deltaR)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_nearestRecoTauIdx",genElectron_matchedRecoElectron_nearestRecoTauIdx)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_nearestRecoTau_pt",genElectron_matchedRecoElectron_nearestRecoTau_pt)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_nearestRecoTau_decayMode",genElectron_matchedRecoElectron_nearestRecoTau_decayMode)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_nearestRecoBoostedTau_deltaR",genElectron_matchedRecoElectron_nearestRecoBoostedTau_deltaR)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_nearestRecoBoostedTauIdx",genElectron_matchedRecoElectron_nearestRecoBoostedTauIdx)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_nearestRecoBoostedTau_pt",genElectron_matchedRecoElectron_nearestRecoBoostedTau_pt)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_nearestRecoBoostedTau_decayMode",genElectron_matchedRecoElectron_nearestRecoBoostedTau_decayMode)

        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoTau_deltaR", genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_deltaR)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoTauIdx", genElectron_matchedRecoElectron_nearestGenMatchedRecoTauIdx)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoTau_pt", genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_pt)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoTau_eta", genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_eta)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoTau_phi", genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_phi)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoTau_decayMode", genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_decayMode)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoTau_idDecayModeNewDMs", genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_idDecayModeNewDMs)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoTau_idDeepTau2018v2p5VSjet", genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_idDeepTau2018v2p5VSjet)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoTau_idDeepTau2018v2p5VSe", genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_idDeepTau2018v2p5VSe)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoTau_idDeepTau2018v2p5VSmu", genElectron_matchedRecoElectron_nearestGenMatchedRecoTau_idDeepTau2018v2p5VSmu)


        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_deltaR", genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_deltaR)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoBoostedTauIdx", genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTauIdx)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_pt", genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_pt)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_eta", genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_eta)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_phi", genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_phi)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_decayMode", genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_decayMode)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_idAntiEle2018", genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_idAntiEle2018)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_idAntiMu", genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_idAntiMu)
        self.out.fillBranch("GenElectronFromHiggsTau_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_rawBoostedDeepTauRunIIv2p0VSjet",genElectron_matchedRecoElectron_nearestGenMatchedRecoBoostedTau_rawBoostedDeepTauRunIIv2p0VSjet)

        self.out.fillBranch("nGenMuonFromHiggsTau", len(genMuon_pt))
        self.out.fillBranch("GenMuonFromHiggsTau_pt", genMuon_pt)
        self.out.fillBranch("GenMuonFromHiggsTau_eta", genMuon_eta)
        self.out.fillBranch("GenMuonFromHiggsTau_phi", genMuon_phi)
        self.out.fillBranch("GenMuonFromHiggsTau_mass", genMuon_mass)
        self.out.fillBranch("GenMuonFromHiggsTau_charge", genMuon_charge)
        self.out.fillBranch("GenMuonFromHiggsTau_genPartIdx", genMuon_genPartIdx)
        self.out.fillBranch("GenMuonFromHiggsTau_tauAncestorIdx", genMuon_tauAncestorIdx)
        self.out.fillBranch("GenMuonFromHiggsTau_higgsAncestorIdx", genMuon_higgsAncestorIdx)
        self.out.fillBranch("GenMuonFromHiggsTau_resonanceAncestorIdx", genMuon_resonanceAncestorIdx)
        self.out.fillBranch("nGenMatchedRecoMuonFromHiggsTau",int(nGenMatchedRecoMuonFromHiggsTau))
        
        self.out.fillBranch("GenMuonFromHiggsTau_matchedRecoMuonIdx", genMuon_matchedRecoMuonIdx)
        self.out.fillBranch("GenMuonFromHiggsTau_hasMatchedRecoMuon", genMuon_hasMatchedRecoMuon)
        self.out.fillBranch("GenMuonFromHiggsTau_matchedRecoMuon_pt", genMuon_matchedRecoMuon_pt)
        self.out.fillBranch("GenMuonFromHiggsTau_matchedRecoMuon_eta", genMuon_matchedRecoMuon_eta)
        self.out.fillBranch("GenMuonFromHiggsTau_matchedRecoMuon_phi", genMuon_matchedRecoMuon_phi)
        self.out.fillBranch("GenMuonFromHiggsTau_matchedRecoMuon_looseId", genMuon_matchedRecoMuon_looseId)
        self.out.fillBranch("GenMuonFromHiggsTau_matchedRecoMuon_mediumId", genMuon_matchedRecoMuon_mediumId)
        self.out.fillBranch("GenMuonFromHiggsTau_matchedRecoMuon_tightId", genMuon_matchedRecoMuon_tightId)
        self.out.fillBranch("GenMuonFromHiggsTau_matchedRecoMuon_genPartFlav", genMuon_matchedRecoMuon_genPartFlav)

        return True


def run_one_file(args):
    input_file, output_dir = args

    base = os.path.splitext(os.path.basename(input_file))[0]
    postfix = ""

    p = PostProcessor(
        output_dir,
        inputFiles=[input_file],      
        cut=None,
        modules=[
            TruthMatchLeptonEfficiencyProducer(
                isData=False,
                debug=False,
                higgs_pdgid=25,
                resonance_pdgid=35,
                max_dr_ele=0.3,
                max_dr_muon=0.3,
                max_dr_hps_tau=0.1,
                max_dr_boosted_tau=0.3,
                json_path=f"{base}_Gen.json",
            )
        ],
        provenance=True,
        fwkJobReport=False,
        noOut=False,
        postfix="",
        haddFileName=None,
        outputbranchsel="Datadrop.txt",
        # maxEntries=10,
        # jsonInput=runsAndLumis(),
    )
    p.run()
    return input_file

if __name__ == "__main__":
    outputDir = "/nfs_scratch/mithakor/ObjectReco_ID_Efficiency/AddedHEEPID"

    inputFiles = [
        "/hdfs/store/user/mithakor/2024_Signal_pythiafixed_original_merged/GluGlutoRadiontoHHto2B2Tau_M-1000_narrow_TuneCP5_13p6TeV_madgraph-pythia8.root",
        "/hdfs/store/user/mithakor/2024_Signal_pythiafixed_original_merged/GluGlutoRadiontoHHto2B2Tau_M-1500_narrow_TuneCP5_13p6TeV_madgraph-pythia8.root",
        "/hdfs/store/user/mithakor/2024_Signal_pythiafixed_original_merged/GluGlutoRadiontoHHto2B2Tau_M-2000_narrow_TuneCP5_13p6TeV_madgraph-pythia8.root",
        "/hdfs/store/user/mithakor/2024_Signal_pythiafixed_original_merged/GluGlutoRadiontoHHto2B2Tau_M-2500_narrow_TuneCP5_13p6TeV_madgraph-pythia8.root",
        "/hdfs/store/user/mithakor/2024_Signal_pythiafixed_original_merged/GluGlutoRadiontoHHto2B2Tau_M-3000_narrow_TuneCP5_13p6TeV_madgraph-pythia8.root",
        "/hdfs/store/user/mithakor/2024_Signal_pythiafixed_original_merged/GluGlutoRadiontoHHto2B2Tau_M-4000_narrow_TuneCP5_13p6TeV_madgraph-pythia8.root",
        "/hdfs/store/user/mithakor/2024_Signal_pythiafixed_original_merged/GluGlutoRadiontoHHto2B2Tau_M-4500_narrow_TuneCP5_13p6TeV_madgraph-pythia8.root" 
    ]

    n_workers = min(40, len(inputFiles)) 

    with mp.get_context("spawn").Pool(processes=n_workers) as pool:
        for finished in pool.imap_unordered(run_one_file, [(f, outputDir) for f in inputFiles]):
            print(f"Finished: {finished}")

    print("DONE")