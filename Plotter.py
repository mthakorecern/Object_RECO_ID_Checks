#!/usr/bin/env python3

import os
import re
import json
import argparse
import multiprocessing as mp
from array import array

import ROOT

ROOT.gROOT.SetBatch(True)


def parse_bins(s):
    return [float(x) for x in s.split(",")]


def sample_label(path):
    base = os.path.basename(path)
    match = re.search(r"M-(\d+)", base)

    if match:
        return "M-{}".format(match.group(1))

    return base.replace(".root", "")


# ----------------------------------------------------------------------
# Electron bitmap helper
# ----------------------------------------------------------------------
def bitmap_passes(bitmap, wp_value, skip_cuts):

    for cut_idx in range(10):

        if cut_idx in skip_cuts:
            continue

        cut_value = (int(bitmap) >> (3 * cut_idx)) & 0x7

        if cut_value < wp_value:
            return False

    return True


# ----------------------------------------------------------------------
# Acceptance helpers
# ----------------------------------------------------------------------
def passes_ele_acceptance(pt, sc_eta, pt_min, eta_max):# exclude_gap=False):

    if pt < pt_min:
        return False

    abs_eta = abs(sc_eta)

    if abs_eta > eta_max:
        return False

    # if exclude_gap:
    #     if 1.479 < abs_eta < 1.566:
    #         return False

    return True



def passes_mu_acceptance(pt, eta, pt_min, eta_max):

    if pt < pt_min:
        return False

    if abs(eta) > eta_max:
        return False

    return True



def passes_tau_acceptance(pt, eta, pt_min, eta_max):

    if pt < pt_min:
        return False

    if abs(eta) > eta_max:
        return False

    return True

def parse_allowed_decay_modes(s):
    if s.strip() == "":
        return None

    return set(int(x) for x in s.split(","))


def passes_tau_decay_mode(decay_mode, allowed_modes):
    if allowed_modes is None:
        return True

    return int(decay_mode) in allowed_modes


# ----------------------------------------------------------------------
# Histograms
# ----------------------------------------------------------------------
def make_hist(name, bins):
    return ROOT.TH1D(name, name, len(bins) - 1, array("d", bins))



def make_efficiency_graph(name, h_num, h_den):
    """
    Build efficiency graph safely.
    Skips bins with zero denominator.
    Keeps graph independent of temporary ROOT objects.
    """

    g = ROOT.TGraphAsymmErrors()
    g.SetName(name)

    # This warning is expected if h_den has empty bins.
    # ROOT will skip those bins.
    g.BayesDivide(h_num, h_den)

    return g


# ----------------------------------------------------------------------
# Tau ID helpers
# ----------------------------------------------------------------------
def hps_tau_wp_passes_vsjet(id_value, wp):
    """
    For HPS Tau ordinal working points.

    Example:
      wp = 0, 1, 2, 3, 4, 5

    Pass condition:
      id_value >= wp
    """
    return int(id_value) >= int(wp)

def hps_tau_wp_passes_vse(id_value, wp):
    return int(id_value) >= int(wp)

def hps_tau_wp_passes_vsmu(id_value, wp):

    return int(id_value) >= int(wp)

def boosted_tau_raw_passes(raw_score, threshold):
    return float(raw_score) >= float(threshold)

HEEP_CUTS = [
    "MinPtCut",                              # bit 0
    "GsfEleSCEtaMultiRangeCut",             # bit 1
    "GsfEleEBEECut_1",                      # bit 2
    "GsfEleEBEECut_2",                      # bit 3
    "GsfEleFull5x5SigmaIEtaIEtaWithSatCut", # bit 4
    "GsfEleFull5x5E2x5OverE5x5WithSatCut",  # bit 5
    "GsfEleHadronicOverEMLinearCut",        # bit 6
    "GsfEleTrkPtIsoCut",                    # bit 7
    "GsfEleEmHadD1IsoRhoCut",               # bit 8
    "GsfEleDxyCut",                         # bit 9
    "GsfEleMissingHitsCut",                 # bit 10
    "GsfEleEcalDrivenCut",                  # bit 11
]

def heep_bitmap_passes(bitmap, skip_cuts=None, ncuts=12):
    """
    HEEP compressed bitmap:
      1 bit per cut
      bit = 1 means pass
      bit = 0 means fail
    """

    if skip_cuts is None:
        skip_cuts = set()

    bitmap = int(bitmap)

    for cut_idx in range(ncuts):

        if cut_idx in skip_cuts:
            continue

        passed = (bitmap >> cut_idx) & 0x1

        if passed == 0:
            return False

    return True

def process_file(job):

    (
        input_file,
        bins,
        ele_wp_value,
        ele_skip_cuts,
        use_ele_cutbased,
        hps_tau_wp_vsjet,
        hps_tau_wp_vse,
        hps_tau_wp_vsmu,
        boosted_tau_raw_threshold,
        allowed_tau_decay_modes,
        allowed_boosted_tau_decay_modes,
        require_tau_newdm_id,
    ) = job

    import ROOT

    ROOT.gROOT.SetBatch(True)

    label = sample_label(input_file)

    f = ROOT.TFile.Open(input_file)

    if not f or f.IsZombie():
        raise RuntimeError("Could not open {}".format(input_file))

    t = f.Get("Events")

    if not t:
        raise RuntimeError("Could not find Events tree")

    # ------------------------------------------------------------------
    # Electron
    # ------------------------------------------------------------------
    h_ele_den = make_hist("h_ele_den_" + label, bins)
    h_ele_reco = make_hist("h_ele_reco_" + label, bins)
    h_ele_recoid = make_hist("h_ele_recoid_" + label, bins)
    h_ele_recoheep = make_hist("h_ele_recoheep_" + label, bins)

    # ------------------------------------------------------------------
    # Muon
    # ------------------------------------------------------------------
    h_mu_den = make_hist("h_mu_den_" + label, bins)
    h_mu_reco = make_hist("h_mu_reco_" + label, bins)
    h_mu_recoid = make_hist("h_mu_recoid_" + label, bins)

    # ------------------------------------------------------------------
    # HPS Tau
    # ------------------------------------------------------------------
    h_tau_den = make_hist("h_tau_den_" + label, bins)
    h_tau_reco = make_hist("h_tau_reco_" + label, bins)
    h_tau_recoid = make_hist("h_tau_recoid_" + label, bins)

    # ------------------------------------------------------------------
    # boostedTau
    # ------------------------------------------------------------------
    h_btau_den = make_hist("h_btau_den_" + label, bins)
    h_btau_reco = make_hist("h_btau_reco_" + label, bins)
    h_btau_recoid = make_hist("h_btau_recoid_" + label, bins)

    nentries = t.GetEntries()

    for iev in range(nentries):

        t.GetEntry(iev)

        # ==============================================================
        # ELECTRONS
        # ==============================================================
        for i in range(int(t.nGenElectronFromHiggsTau)):
            gen_pt = float(t.GenElectronFromHiggsTau_pt[i])
            gen_eta = float(t.GenElectronFromHiggsTau_eta[i])

            if gen_pt < 10:
                continue

            if abs(gen_eta) > 2.5:
                continue

            h_ele_den.Fill(gen_pt)

            if int(t.GenElectronFromHiggsTau_hasMatchedRecoElectron[i]) != 1:
                continue

            reco_pt = float(t.GenElectronFromHiggsTau_matchedRecoElectron_pt[i])
            # reco_sc_eta = float(t.GenElectronFromHiggsTau_matchedRecoElectron_superclusterEta[i])
            reco_eta = float(t.GenElectronFromHiggsTau_matchedRecoElectron_eta[i])


            if not passes_ele_acceptance(reco_pt, reco_eta, 10, 2.5):
                continue

            h_ele_reco.Fill(gen_pt)

            if use_ele_cutbased:

                cutbased = int(
                    t.GenElectronFromHiggsTau_matchedRecoElectron_cutBased[i]
                )

                ele_pass = cutbased >= ele_wp_value

            else:

                bitmap = int(
                    t.GenElectronFromHiggsTau_matchedRecoElectron_vidNestedWPBitmap[i]
                )

                ele_pass = bitmap_passes(
                    bitmap,
                    ele_wp_value,
                    ele_skip_cuts,
                )

            if ele_pass:
                h_ele_recoid.Fill(gen_pt)
                heep_bitmap = int(t.GenElectronFromHiggsTau_matchedRecoElectron_vidNestedWPBitmapHEEP[i])
                heep_pass = heep_bitmap_passes(heep_bitmap, skip_cuts={7, 8}, ncuts=12)
                if heep_pass:
                    h_ele_recoheep.Fill(gen_pt)

        for i in range(int(t.nGenMuonFromHiggsTau)):

            gen_pt = float(t.GenMuonFromHiggsTau_pt[i])
            gen_eta = float(t.GenMuonFromHiggsTau_eta[i])
            if gen_pt < 15:
                continue

            if abs(gen_eta) > 2.4:
                continue

            h_mu_den.Fill(gen_pt)

            if int(t.GenMuonFromHiggsTau_hasMatchedRecoMuon[i]) != 1:
                continue

            reco_pt = float(t.GenMuonFromHiggsTau_matchedRecoMuon_pt[i])
            reco_eta = float(t.GenMuonFromHiggsTau_matchedRecoMuon_eta[i])

            if not passes_mu_acceptance(
                reco_pt,
                reco_eta,
                15,
                2.4,
            ):
                continue

            h_mu_reco.Fill(gen_pt)

            loose_id = int(
                t.GenMuonFromHiggsTau_matchedRecoMuon_looseId[i]
            )

            if loose_id == 1:
                h_mu_recoid.Fill(gen_pt)

        # ==============================================================
        # HPS TAUS
        # ==============================================================
        for i in range(int(t.nGenVisTauFromHiggsTau)):

            gen_pt = float(t.GenVisTauFromHiggsTau_pt[i])
            gen_eta = float(t.GenVisTauFromHiggsTau_eta[i])

            # GEN fiducial denominator
            if gen_pt < 20:
                continue

            if abs(gen_eta) > 2.5:
                continue

            h_tau_den.Fill(gen_pt)

            # Matched reco HPS tau
            if int(t.GenVisTauFromHiggsTau_hasMatchedRecoTau[i]) != 1:
                continue

            reco_pt = float(
                t.GenVisTauFromHiggsTau_matchedRecoTau_pt[i]
            )

            reco_eta = float(
                t.GenVisTauFromHiggsTau_matchedRecoTau_eta[i]
            )

            reco_dz = float(
                t.GenVisTauFromHiggsTau_matchedRecoTau_dz[i]
            )

            if abs(reco_dz) > 0.2 :
                continue 


            # Reco pT/eta acceptance
            if not passes_tau_acceptance(reco_pt, reco_eta, 20, 2.5):
                continue

            # Reco decay mode requirement
            reco_decay_mode = int(
                t.GenVisTauFromHiggsTau_matchedRecoTau_decayMode[i]
            )

            if not passes_tau_decay_mode(
                reco_decay_mode,
                allowed_tau_decay_modes,
            ):
                continue

            # Optional HPS Tau newDM decay-mode ID
            if require_tau_newdm_id:
                newdm_id = int(
                    t.GenVisTauFromHiggsTau_matchedRecoTau_idDecayModeNewDMs[i]
                )

                if newdm_id != 1:
                    continue

            # Reconstruction numerator:
            # matched reco tau + reco pT/eta + decay-mode requirement
            h_tau_reco.Fill(gen_pt)

            # Reco + ID numerator:
            # same as reco numerator + DeepTau VSjet WP
            hps_tau_vsjet_id = int(
                t.GenVisTauFromHiggsTau_matchedRecoTau_idDeepTau2018v2p5VSjet[i]
            )
            hps_tau_vse_id = int(
                t.GenVisTauFromHiggsTau_matchedRecoTau_idDeepTau2018v2p5VSe[i]
            )
            hps_tau_vsmu_id = int(
                t.GenVisTauFromHiggsTau_matchedRecoTau_idDeepTau2018v2p5VSmu[i]
            )

            if hps_tau_wp_passes_vsjet(hps_tau_vsjet_id, hps_tau_wp_vsjet) and hps_tau_wp_passes_vse(hps_tau_vse_id, hps_tau_wp_vse) and hps_tau_wp_passes_vsmu(hps_tau_vsmu_id, hps_tau_wp_vsmu)  :
                h_tau_recoid.Fill(gen_pt)

        # ==============================================================
        # BOOSTED TAUS
        # ==============================================================
        for i in range(int(t.nGenVisTauFromHiggsTau)):

            gen_pt = float(t.GenVisTauFromHiggsTau_pt[i])
            gen_eta = float(t.GenVisTauFromHiggsTau_eta[i])

            # GEN fiducial denominator
            if gen_pt < 25:
                continue

            if abs(gen_eta) > 2.5:
                continue

            h_btau_den.Fill(gen_pt)

            # Matched reco boostedTau
            if int(
                t.GenVisTauFromHiggsTau_hasMatchedRecoBoostedTau[i]
            ) != 1:
                continue

            reco_pt = float(
                t.GenVisTauFromHiggsTau_matchedRecoBoostedTau_pt[i]
            )

            reco_eta = float(
                t.GenVisTauFromHiggsTau_matchedRecoBoostedTau_eta[i]
            )

            # Reco pT/eta acceptance
            if not passes_tau_acceptance(reco_pt, reco_eta, 25, 2.5):
                continue

            # Reco decay mode requirement
            reco_decay_mode = int(
                t.GenVisTauFromHiggsTau_matchedRecoBoostedTau_decayMode[i]
            )

            if not passes_tau_decay_mode(
                reco_decay_mode,
                allowed_boosted_tau_decay_modes,
            ):
                continue

            # Reconstruction numerator:
            # matched boostedTau + reco pT/eta + decay-mode requirement
            h_btau_reco.Fill(gen_pt)

            anti_ele = int(
                t.GenVisTauFromHiggsTau_matchedRecoBoostedTau_idAntiEle2018[i]
            )

            anti_mu = int(
                t.GenVisTauFromHiggsTau_matchedRecoBoostedTau_idAntiMu[i]
            )

            raw_boosted_deeptau = float(
                t.GenVisTauFromHiggsTau_matchedRecoBoostedTau_rawBoostedDeepTauRunIIv2p0VSjet[i]
            )

            # pass_anti_ele = ((anti_ele & 2) == 2)
            # pass_anti_mu = ((anti_mu & 1) == 1)

            pass_raw_boosted_deeptau = boosted_tau_raw_passes(
                raw_boosted_deeptau,
                boosted_tau_raw_threshold,
            )

            if pass_raw_boosted_deeptau:    #pass_anti_ele and pass_anti_mu and
                h_btau_recoid.Fill(gen_pt)

    # ------------------------------------------------------------------
    # Summary
    # ------------------------------------------------------------------
    summary = {
        "label": label,

        "electron": {
            "den": int(h_ele_den.Integral()),
            "reco": int(h_ele_reco.Integral()),
            "recoid": int(h_ele_recoid.Integral()),
            "recoheep": int(h_ele_recoheep.Integral()),

        },

        "muon": {
            "den": int(h_mu_den.Integral()),
            "reco": int(h_mu_reco.Integral()),
            "recoid": int(h_mu_recoid.Integral()),
        },

        "tau": {
            "den": int(h_tau_den.Integral()),
            "reco": int(h_tau_reco.Integral()),
            "recoid": int(h_tau_recoid.Integral()),
        },

        "boostedTau": {
            "den": int(h_btau_den.Integral()),
            "reco": int(h_btau_reco.Integral()),
            "recoid": int(h_btau_recoid.Integral()),
        },
    }

    # detach from file
    for h in [
        h_ele_den,
        h_ele_reco,
        h_ele_recoid,
        h_ele_recoheep,
        h_mu_den,
        h_mu_reco,
        h_mu_recoid,
        h_tau_den,
        h_tau_reco,
        h_tau_recoid,
        h_btau_den,
        h_btau_reco,
        h_btau_recoid,
    ]:
        h.SetDirectory(0)

    f.Close()

    return {
        "label": label,
        "summary": summary,

        "h_ele_den": h_ele_den,
        "h_ele_reco": h_ele_reco,
        "h_ele_recoid": h_ele_recoid,
        "h_ele_recoheep": h_ele_recoheep,


        "h_mu_den": h_mu_den,
        "h_mu_reco": h_mu_reco,
        "h_mu_recoid": h_mu_recoid,

        "h_tau_den": h_tau_den,
        "h_tau_reco": h_tau_reco,
        "h_tau_recoid": h_tau_recoid,

        "h_btau_den": h_btau_den,
        "h_btau_reco": h_btau_reco,
        "h_btau_recoid": h_btau_recoid,
    }



def draw_graphs(results, outdir):

    colors = [
        ROOT.kBlack,
        ROOT.kRed + 1,
        ROOT.kBlue + 1,
        ROOT.kGreen + 2,
        ROOT.kMagenta + 1,
        ROOT.kOrange + 7,
        ROOT.kYellow + 2,

    ]

    plot_configs = [
        ("electron", "Reco", "h_ele_reco", "h_ele_den"),
        ("electron", "Reco+ID", "h_ele_recoid", "h_ele_den"),
        ("electron", "Reco+HEEP", "h_ele_recoheep", "h_ele_den"),


        ("muon", "Reco", "h_mu_reco", "h_mu_den"),
        ("muon", "Reco+ID", "h_mu_recoid", "h_mu_den"),

        ("tau", "Reco", "h_tau_reco", "h_tau_den"),
        ("tau", "Reco+ID", "h_tau_recoid", "h_tau_den"),

        ("boostedTau", "Reco", "h_btau_reco", "h_btau_den"),
        ("boostedTau", "Reco+ID", "h_btau_recoid", "h_btau_den"),
    ]

    # Important: keep references alive until the function exits.
    all_graphs = []
    all_canvases = []

    for obj, mode, num_key, den_key in plot_configs:

        c = ROOT.TCanvas(
            "c_{}_{}".format(obj, mode),
            "c_{}_{}".format(obj, mode),
            900,
            700,
        )
        c.SetGrid()

        leg = ROOT.TLegend(0.50, 0.15, 0.90, 0.45)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.04)

        first = True
        graphs_this_canvas = []

        for idx, res in enumerate(results):

            h_num = res[num_key]
            h_den = res[den_key]

            # Skip samples with completely empty denominator.
            if h_den.Integral() <= 0:
                print(
                    "[draw_graphs] Skipping {} {} for {}: denominator is empty".format(
                        obj,
                        mode,
                        res["label"],
                    )
                )
                continue

            g = make_efficiency_graph(
                "g_{}_{}_{}".format(obj, mode, idx),
                h_num,
                h_den,
            )

            # If all denominator-filled bins were skipped, do not draw.
            if g.GetN() == 0:
                print(
                    "[draw_graphs] Skipping {} {} for {}: graph has zero points".format(
                        obj,
                        mode,
                        res["label"],
                    )
                )
                continue

            color = colors[idx % len(colors)]

            g.SetLineColor(color)
            g.SetMarkerColor(color)
            g.SetMarkerStyle(20 + idx)
            g.SetLineWidth(2)

            graphs_this_canvas.append(g)
            all_graphs.append(g)

            if first:
                g.SetTitle("{} {} efficiency".format(obj, mode))
                g.Draw("AP")

                g.GetXaxis().SetTitle("GEN p_{T} [GeV]")
                g.GetYaxis().SetTitle("Efficiency")
                g.GetYaxis().SetRangeUser(0.0, 1.05)

                first = False
            else:
                g.Draw("P SAME")

            leg.AddEntry(g, res["label"], "lep")

        # If no valid graph was drawn, skip SaveAs.
        if first:
            print(
                "[draw_graphs] No valid graphs for {} {}. Skipping canvas.".format(
                    obj,
                    mode,
                )
            )
            c.Close()
            continue

        leg.Draw()
        c.Update()

        # Keep objects attached to canvas so PyROOT does not garbage collect.
        c._graphs = graphs_this_canvas
        c._legend = leg

        out_png = os.path.join(
            outdir,
            "{}_{}_efficiency.png".format(obj, mode),
        )


        c.SaveAs(out_png)

        all_canvases.append(c)

    # Keep references alive until all SaveAs calls are done.
    return all_graphs, all_canvases


# ----------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------
def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--inputs",
        nargs="+",
        required=True,
    )

    parser.add_argument(
        "--output-dir",
        default="Efficiency_plots",
    )

    parser.add_argument(
        "--bins",
        default="10, 20, 40, 80, 120, 160, 200, 250, 300, 350, 400, 500, 600, 700, 800, 900, 1000",
    )

    parser.add_argument(
        "--use-ele-cutbased",
        action="store_true",
    )

    parser.add_argument(
        "--ele-skip-cuts",
        default="",
    )
    parser.add_argument(
    "--tau-decay-modes",
    default="0,1,10,11",
    help=(
        "Allowed HPS Tau decay modes. "
        "Use empty string '' to disable the decay-mode cut."
    ),
    )

    parser.add_argument(
        "--boosted-tau-decay-modes",
        default="0,1,10,11",
        help=(
            "Allowed boostedTau decay modes. "
            "Use empty string '' to disable the decay-mode cut."
        ),
    )

    parser.add_argument(
        "--require-tau-newdm-id",
        action="store_true",
        help="Require matched HPS Tau_idDecayModeNewDMs == 1.",
    )
    parser.add_argument(
        "--boosted-tau-raw-threshold",
        type=float,
        default=0.85,
        help=(
            "Raw boostedTau DeepTau threshold. "
            "Pass condition is rawBoostedDeepTauRunIIv2p0VSjet >= threshold."
        ),
    )

    parser.add_argument(
        "--nproc",
        type=int,
        default=4,
    )

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    bins = parse_bins(args.bins)

    ele_skip_cuts = set()

    if args.ele_skip_cuts != "":
        ele_skip_cuts = set(
            int(x)
            for x in args.ele_skip_cuts.split(",")
        )
    allowed_tau_decay_modes = parse_allowed_decay_modes(args.tau_decay_modes)

    allowed_boosted_tau_decay_modes = parse_allowed_decay_modes(args.boosted_tau_decay_modes)

    ele_wp_value = 2

    # DeepTau working points
    hps_tau_wp_vsjet = 4
    hps_tau_wp_vse = 2
    hps_tau_wp_vsmu = 1
    boosted_tau_raw_threshold = args.boosted_tau_raw_threshold

    jobs = []

    for infile in args.inputs:

        jobs.append((
            infile,
            bins,
            ele_wp_value,
            ele_skip_cuts,
            args.use_ele_cutbased,
            hps_tau_wp_vsjet,
            hps_tau_wp_vse,
            hps_tau_wp_vsmu,
            boosted_tau_raw_threshold,
            allowed_tau_decay_modes,
            allowed_boosted_tau_decay_modes,
            args.require_tau_newdm_id,
        ))

    ctx = mp.get_context("spawn")

    if args.nproc > 1:

        with ctx.Pool(processes=args.nproc) as pool:
            results = pool.map(process_file, jobs)

    else:
        results = [process_file(job) for job in jobs]

    graphs, canvases = draw_graphs(results, args.output_dir)
    
    with open(
        os.path.join(args.output_dir, "summary.json"),
        "w",
    ) as f:

        json.dump(
            [r["summary"] for r in results],
            f,
            indent=4,
        )

    print("Done")


if __name__ == "__main__":
    main()