// ── Virtual-screening benchmark metrics ─────────────────────────────────────
//
// Given a set of docked compounds with a score and an active/decoy label, these
// functions quantify how well the docking ranks actives ahead of decoys — the
// numbers a DUD-E / CASF-style retrospective validation reports.
//
// Convention: **lower score = better binder** (docking energies are negative),
// so actives are expected near the TOP of an ascending sort. `label = true`
// marks an active.
//
// This is the evaluation harness only. The compound library + labels come from
// the user (e.g. a DUD-E target downloaded from the web); feed each ligand's
// `dock(...)` score and its known label into these functions.

/// Area under the ROC curve: probability a random active scores *better*
/// (lower) than a random decoy. 1.0 = perfect ranking, 0.5 = random.
///
/// Uses the Mann-Whitney rank formula with average ranks for ties; O(n log n).
pub fn roc_auc(scores: &[f64], labels: &[bool]) -> Option<f64> {
    let n = scores.len();
    if n == 0 || n != labels.len() {
        return None;
    }
    let n_pos = labels.iter().filter(|&&l| l).count();
    let n_neg = n - n_pos;
    if n_pos == 0 || n_neg == 0 {
        return None; // AUC undefined without both classes
    }

    // Rank by ascending score (rank 1 = lowest = best). Ties share the average
    // rank. Sum the ranks of the positives.
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&a, &b| scores[a].partial_cmp(&scores[b]).unwrap());

    let mut ranks = vec![0.0_f64; n];
    let mut i = 0;
    while i < n {
        let mut j = i + 1;
        while j < n && scores[idx[j]] == scores[idx[i]] {
            j += 1;
        }
        // average rank for the tie group [i, j)  (1-based)
        let avg = ((i + 1 + j) as f64) / 2.0;
        for k in i..j {
            ranks[idx[k]] = avg;
        }
        i = j;
    }

    let r_pos: f64 = (0..n).filter(|&k| labels[k]).map(|k| ranks[k]).sum();
    // Mann-Whitney U for positives ranked low (= better). With ascending ranks,
    // smaller rank = better, so the probability active < decoy is:
    let u = r_pos - (n_pos * (n_pos + 1)) as f64 / 2.0; // # pairs decoy ranked above active
    let auc_active_worse = u / (n_pos as f64 * n_neg as f64);
    // auc_active_worse = P(active ranked higher/worse than decoy). We want the
    // complement: P(active scores lower/better than decoy).
    Some(1.0 - auc_active_worse)
}

/// Enrichment factor at the top `fraction` (0<χ≤1) of the ranked list:
/// (actives in top χ / size of top χ) / (total actives / N).
/// EF = 1 is random; the theoretical max is 1/χ (or N/A, whichever is smaller).
pub fn enrichment_factor(scores: &[f64], labels: &[bool], fraction: f64) -> Option<f64> {
    let n = scores.len();
    if n == 0 || n != labels.len() || !(0.0..=1.0).contains(&fraction) || fraction == 0.0 {
        return None;
    }
    let total_active = labels.iter().filter(|&&l| l).count();
    if total_active == 0 {
        return None;
    }

    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&a, &b| scores[a].partial_cmp(&scores[b]).unwrap());

    // Top χ fraction (at least one compound).
    let top = ((fraction * n as f64).round() as usize).max(1).min(n);
    let actives_in_top = idx[..top].iter().filter(|&&k| labels[k]).count();

    let observed = actives_in_top as f64 / top as f64;
    let baseline = total_active as f64 / n as f64;
    Some(observed / baseline)
}

/// BEDROC (Boltzmann-Enhanced Discrimination of ROC), Truchon & Bayly 2007.
/// Early-recognition metric with weighting parameter `alpha` (common: 20 ≈ top
/// 8 %, 80 ≈ top 2 %). 1.0 = all actives ranked first, 0.0 = ranked last.
pub fn bedroc(scores: &[f64], labels: &[bool], alpha: f64) -> Option<f64> {
    let n = scores.len();
    if n == 0 || n != labels.len() || alpha <= 0.0 {
        return None;
    }
    let n_pos = labels.iter().filter(|&&l| l).count();
    if n_pos == 0 || n_pos == n {
        return None;
    }
    let big_n = n as f64;
    let n_a = n_pos as f64;
    let ra = n_a / big_n; // ratio of actives

    // Ranks (1-based) of actives, best = lowest score = rank 1.
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&a, &b| scores[a].partial_cmp(&scores[b]).unwrap());

    let mut sum_exp = 0.0;
    for (rank0, &i) in idx.iter().enumerate() {
        if labels[i] {
            let rank = (rank0 + 1) as f64;
            sum_exp += (-alpha * rank / big_n).exp();
        }
    }

    // RIE = mean(exp(-α r/N)) / expected-under-random
    let denom_rie = (1.0 / big_n) * ((1.0 - (-alpha).exp()) / ((alpha / big_n).exp() - 1.0));
    let rie = (sum_exp / n_a) / denom_rie;

    // BEDROC from RIE (closed form).
    let sinh_half = (alpha / 2.0).sinh();
    let cosh_half = (alpha / 2.0).cosh();
    let factor = ra * sinh_half / (cosh_half - (alpha / 2.0 - alpha * ra).cosh());
    let offset = 1.0 / (1.0 - (alpha * (1.0 - ra)).exp());
    Some(rie * factor + offset)
}

// ─────────────────────────────────────────────────────────────────────────────
// Tests
// ─────────────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn auc_perfect_ranking() {
        // actives all score lower (better) than decoys → AUC = 1.
        let scores = vec![-5.0, -4.0, -1.0, 0.0];
        let labels = vec![true, true, false, false];
        let auc = roc_auc(&scores, &labels).unwrap();
        assert!((auc - 1.0).abs() < 1e-9, "auc={auc}");
    }

    #[test]
    fn auc_worst_ranking() {
        // actives score worst → AUC = 0.
        let scores = vec![0.0, -1.0, -4.0, -5.0];
        let labels = vec![true, true, false, false];
        let auc = roc_auc(&scores, &labels).unwrap();
        assert!(auc.abs() < 1e-9, "auc={auc}");
    }

    #[test]
    fn auc_half_when_actives_split_best_and_worst() {
        // one active best (rank 1), one active worst (rank 4): each active beats
        // exactly one decoy → AUC = 0.5.
        let scores = vec![-5.0, -4.0, -3.0, -2.0];
        let labels = vec![true, false, false, true];
        let auc = roc_auc(&scores, &labels).unwrap();
        assert!((auc - 0.5).abs() < 1e-9, "auc={auc}");
    }

    #[test]
    fn auc_better_than_half_when_actives_lead() {
        // actives at ranks 1 and 3 → 3 of 4 active-decoy pairs favour the active.
        let scores = vec![-5.0, -4.0, -3.0, -2.0];
        let labels = vec![true, false, true, false];
        let auc = roc_auc(&scores, &labels).unwrap();
        assert!((auc - 0.75).abs() < 1e-9, "auc={auc}");
    }

    #[test]
    fn auc_handles_ties() {
        // all tied → AUC = 0.5 (no discrimination)
        let scores = vec![-1.0, -1.0, -1.0, -1.0];
        let labels = vec![true, false, true, false];
        let auc = roc_auc(&scores, &labels).unwrap();
        assert!((auc - 0.5).abs() < 1e-9, "tied auc={auc}");
    }

    #[test]
    fn auc_requires_both_classes() {
        assert!(roc_auc(&[-1.0, -2.0], &[true, true]).is_none());
        assert!(roc_auc(&[-1.0, -2.0], &[false, false]).is_none());
    }

    #[test]
    fn ef_perfect_top_fraction() {
        // 2 actives among 10, both in the top 20% → EF = (2/2)/(2/10) = 5.
        let mut scores = vec![0.0; 10];
        let mut labels = vec![false; 10];
        scores[0] = -10.0; labels[0] = true;
        scores[1] = -9.0; labels[1] = true;
        for i in 2..10 { scores[i] = -(i as f64) * 0.1; }
        let ef = enrichment_factor(&scores, &labels, 0.2).unwrap();
        assert!((ef - 5.0).abs() < 1e-9, "ef={ef}");
    }

    #[test]
    fn ef_random_is_one() {
        // actives spread uniformly → EF ≈ 1 at 50%.
        let scores = vec![-1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -7.0, -8.0];
        let labels = vec![true, false, true, false, true, false, true, false];
        let ef = enrichment_factor(&scores, &labels, 0.5).unwrap();
        assert!((ef - 1.0).abs() < 1e-9, "ef={ef}");
    }

    #[test]
    fn bedroc_perfect_near_one() {
        // all actives at the very top → BEDROC close to 1.
        let mut scores = vec![0.0; 20];
        let mut labels = vec![false; 20];
        for i in 0..4 { scores[i] = -10.0 + i as f64; labels[i] = true; }
        for i in 4..20 { scores[i] = i as f64; }
        let b = bedroc(&scores, &labels, 20.0).unwrap();
        assert!(b > 0.9, "perfect bedroc should be ~1: {b}");
    }

    #[test]
    fn bedroc_worst_near_zero() {
        // all actives at the bottom → BEDROC close to 0.
        let mut scores = vec![0.0; 20];
        let mut labels = vec![false; 20];
        for i in 0..16 { scores[i] = -(i as f64); }
        for i in 16..20 { scores[i] = 100.0 + i as f64; labels[i] = true; }
        let b = bedroc(&scores, &labels, 20.0).unwrap();
        assert!(b < 0.1, "worst bedroc should be ~0: {b}");
    }

    #[test]
    fn bedroc_monotonic_in_quality() {
        // Moving an active earlier should not decrease BEDROC.
        let labels = vec![true, false, false, false, false, false, false, false];
        let early = vec![-8.0, -1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -7.0]; // active first
        let late = vec![-0.5, -1.0, -2.0, -3.0, -4.0, -5.0, -6.0, -7.0];  // active last
        let b_early = bedroc(&early, &labels, 20.0).unwrap();
        let b_late = bedroc(&late, &labels, 20.0).unwrap();
        assert!(b_early > b_late, "earlier active should score higher: {b_early} vs {b_late}");
    }
}
