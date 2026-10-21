//! Bounded classical block designs and restricted randomization.
//!
//! Treatment-design inputs are distinct rows of numeric treatment settings.
//! Split/strip inputs are factor rows [low,high,role0/1], both roles required.
//! Taguchi input is one scalar level count per factor; L4/L8/L9/L18 only.
//! These generate plans and structural diagnostics, not response-model inference.
type Mat = Vec<Vec<f64>>;
const MAX_ROWS: usize = 16384;
const MAX_CELLS: usize = 1_000_000;
pub fn run(op: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    let w = data.first().map_or(0, Vec::len);
    if w == 0
        || w > 64
        || data.len() > 4096
        || data.len() > MAX_CELLS / w
        || data
            .iter()
            .any(|r| r.len() != w || r.iter().any(|x| !x.is_finite()))
        || args.len() > 8
        || args.iter().any(|x| !x.is_finite())
    {
        return Err(
            "data must be finite rectangular nonempty and bounded; arguments finite".into(),
        );
    }
    match op {
        "rcbd_design" => rcbd(data, args),
        "balanced_incomplete_block" => bibd(data, args),
        "resolvable_incomplete_block" => resolvable(data, args),
        "split_plot_design" => restricted(data, args, false),
        "strip_plot_design" => restricted(data, args, true),
        "taguchi_orthogonal_array" => orthogonal(data, args),
        _ => Err("unknown classical-design operation".into()),
    }
}
fn nums(a: &[f64]) -> String {
    format!(
        "[{}]",
        a.iter().map(f64::to_string).collect::<Vec<_>>().join(",")
    )
}
fn ints(a: &[usize]) -> String {
    format!(
        "[{}]",
        a.iter().map(usize::to_string).collect::<Vec<_>>().join(",")
    )
}
fn mat(a: &Mat) -> String {
    format!(
        "[{}]",
        a.iter().map(|r| nums(r)).collect::<Vec<_>>().join(",")
    )
}
fn imat(a: &[Vec<usize>]) -> String {
    format!(
        "[{}]",
        a.iter().map(|r| ints(r)).collect::<Vec<_>>().join(",")
    )
}
fn uint(x: f64, name: &str) -> Result<usize, String> {
    if (0. ..=9_007_199_254_740_991.).contains(&x) && x.fract() == 0. {
        Ok(x as usize)
    } else {
        Err(format!(
            "{name} must be an exact nonnegative integer<=2^53-1"
        ))
    }
}
fn seed(args: &[f64], i: usize) -> Result<u64, String> {
    Ok(uint(args.get(i).copied().unwrap_or(0.), "seed")? as u64)
}
struct Rng(u64);
impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9e3779b97f4a7c15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        z ^ (z >> 31)
    }
    fn shuffle<T>(&mut self, a: &mut [T]) {
        for i in (1..a.len()).rev() {
            let n = (i + 1) as u64;
            let threshold = n.wrapping_neg() % n;
            let mut z = self.next();
            while z < threshold {
                z = self.next();
            }
            a.swap(i, (z % n) as usize);
        }
    }
}
fn treatment_input(data: &[Vec<f64>], max: usize) -> Result<(), String> {
    if !(2..=max).contains(&data.len()) || data[0].len() > 16 {
        return Err(format!(
            "requires2..{max}distinct treatments and<=16setting columns"
        ));
    }
    for (i, r) in data.iter().enumerate() {
        if data[..i].contains(r) {
            return Err("treatment setting rows must be distinct".into());
        }
    }
    Ok(())
}
fn rows_bound(n: usize, p: usize) -> Result<(), String> {
    if n > MAX_ROWS || n > MAX_CELLS / p {
        Err("design exceeds16384runs or1000000output cells".into())
    } else {
        Ok(())
    }
}
fn block_result(
    data: &[Vec<f64>],
    blocks: &[Vec<usize>],
    classes: Option<&[usize]>,
    s: u64,
    method: &str,
    search: usize,
) -> Result<String, String> {
    let v = data.len();
    let b = blocks.len();
    let k = blocks[0].len();
    rows_bound(b * k, data[0].len())?;
    let mut d = Vec::new();
    let mut tid = Vec::new();
    let mut bid = Vec::new();
    let mut cl = Vec::new();
    let mut inc = vec![vec![0; b]; v];
    for (j, block) in blocks.iter().enumerate() {
        for &t in block {
            d.push(data[t].clone());
            tid.push(t);
            bid.push(j);
            inc[t][j] += 1;
            if let Some(c) = classes {
                cl.push(c[j]);
            }
        }
    }
    let mut concurrence = vec![vec![0; v]; v];
    for (i, row) in concurrence.iter_mut().enumerate() {
        for (j, x) in row.iter_mut().enumerate() {
            *x = inc[i].iter().zip(&inc[j]).map(|(a, b)| a * b).sum();
        }
    }
    let replication: Vec<_> = inc.iter().map(|r| r.iter().sum()).collect();
    let mut info = vec![vec![0.; v]; v];
    for i in 0..v {
        for j in 0..v {
            info[i][j] = (if i == j { replication[i] as f64 } else { 0. })
                - concurrence[i][j] as f64 / k as f64;
        }
    }
    let extra = if classes.is_some() {
        format!(",\"parallel_class_ids\":{}", ints(&cl))
    } else {
        String::new()
    };
    Ok(format!("{{\"version\":1,\"method\":\"{method}\",\"seed\":{s},\"rng\":\"SplitMix64 rejection-sampled Fisher-Yates\",\"design\":{},\"treatment_ids\":{},\"block_ids\":{},\"blocks\":{},\"n_blocks\":{b},\"block_size\":{k},\"replication\":{},\"pair_concurrence\":{},\"treatment_information\":{},\"treatment_information_rank\":{},\"run_order\":{},\"search_complete\":true,\"search_nodes\":{search},\"inference\":\"plan and block-adjusted treatment information only;no response analysis\"{extra}}}",mat(&d),ints(&tid),ints(&bid),imat(blocks),ints(&replication),imat(&concurrence),mat(&info),v-1,ints(&(1..=d.len()).collect::<Vec<_>>())))
}
fn rcbd(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    treatment_input(data, 256)?;
    if args.is_empty() || args.len() > 2 {
        return Err("RCBD requires[blocks,seed=0]".into());
    }
    let b = uint(args[0], "blocks")?;
    if !(1..=1024).contains(&b) {
        return Err("blocks must be1..1024".into());
    }
    rows_bound(b * data.len(), data[0].len())?;
    let s = seed(args, 1)?;
    let mut rng = Rng(s);
    let mut blocks = vec![(0..data.len()).collect::<Vec<_>>(); b];
    for block in &mut blocks {
        rng.shuffle(block);
    }
    rng.shuffle(&mut blocks);
    block_result(data, &blocks, None, s, "randomized_complete_block", 0)
}
struct Search {
    v: usize,
    target: usize,
    candidates: Vec<Vec<usize>>,
    remaining: Vec<Vec<usize>>,
    rep: Vec<usize>,
    chosen: Vec<Vec<usize>>,
    nodes: usize,
    work: usize,
}
impl Search {
    fn eligible(&mut self, idx: usize) -> Result<bool, String> {
        self.work += 1;
        if self.work > 5_000_000 {
            return Err("search_limit: BIBD candidate feasibility budget5000000 exhausted; existence not determined".into());
        }
        let c = &self.candidates[idx];
        Ok(c.iter().all(|&t| self.rep[t] > 0)
            && c.iter()
                .enumerate()
                .all(|(i, &a)| c[..i].iter().all(|&b| self.remaining[a][b] > 0)))
    }
    fn alter(&mut self, idx: usize, add: bool) {
        let c = &self.candidates[idx];
        for (i, &a) in c.iter().enumerate() {
            if add {
                self.rep[a] -= 1;
            } else {
                self.rep[a] += 1;
            }
            for &b in &c[..i] {
                if add {
                    self.remaining[a][b] -= 1;
                    self.remaining[b][a] -= 1;
                } else {
                    self.remaining[a][b] += 1;
                    self.remaining[b][a] += 1;
                }
            }
        }
    }
    fn solve(&mut self) -> Result<bool, String> {
        self.nodes += 1;
        if self.nodes > 100_000 {
            return Err(
                "search_limit: BIBD node budget100000 exhausted; existence not determined".into(),
            );
        }
        if self.chosen.len() == self.target {
            return Ok(self.rep.iter().all(|x| *x == 0)
                && self.remaining.iter().flatten().all(|x| *x == 0));
        }
        let mut best: Option<Vec<usize>> = None;
        for a in 0..self.v {
            for b in 0..a {
                if self.remaining[a][b] == 0 {
                    continue;
                }
                let mut list = Vec::new();
                for ci in 0..self.candidates.len() {
                    if self.candidates[ci].contains(&a)
                        && self.candidates[ci].contains(&b)
                        && self.eligible(ci)?
                    {
                        list.push(ci);
                    }
                }
                if list.is_empty() {
                    return Ok(false);
                }
                if best.as_ref().is_none_or(|old| list.len() < old.len()) {
                    best = Some(list);
                }
            }
        }
        let Some(list) = best else {
            return Ok(false);
        };
        for ci in list {
            self.alter(ci, true);
            self.chosen.push(self.candidates[ci].clone());
            if self.solve()? {
                return Ok(true);
            }
            self.chosen.pop();
            self.alter(ci, false);
        }
        Ok(false)
    }
}
fn combinations(v: usize, k: usize, start: usize, c: &mut Vec<usize>, out: &mut Vec<Vec<usize>>) {
    if c.len() == k {
        out.push(c.clone());
        return;
    }
    for i in start..=v - (k - c.len()) {
        c.push(i);
        combinations(v, k, i + 1, c, out);
        c.pop();
    }
}
fn bibd(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    treatment_input(data, 9)?;
    if args.is_empty() || args.len() > 3 {
        return Err("BIBD requires[block_size,lambda=1,seed=0]".into());
    }
    let v = data.len();
    let k = uint(args[0], "block_size")?;
    let lam = uint(args.get(1).copied().unwrap_or(1.), "lambda")?;
    if k < 2 || k >= v || k > 4 || !(1..=3).contains(&lam) {
        return Err("unsupported: exact BIBD search supports2<=k<=4,k<v,lambda1..3,v<=9".into());
    }
    if !(lam * (v - 1)).is_multiple_of(k - 1) || !(lam * v * (v - 1)).is_multiple_of(k * (k - 1)) {
        return Err("no_solution: BIBD integer replication/block equations are violated".into());
    }
    let r = lam * (v - 1) / (k - 1);
    let b = v * r / k;
    if b > 128 {
        return Err("unsupported: BIBD requires<=128blocks".into());
    }
    let s = seed(args, 2)?;
    // Certified algebraic constructions avoid exponential search for these
    // small saturated families; neither is a full factorial treatment plan.
    let direct = if k == 4 && lam == 3 && v == 8 {
        Some(
            (1_usize..8)
                .flat_map(|normal| {
                    (0..2).map(move |side| {
                        (0_usize..8)
                            .filter(|point| ((point & normal).count_ones() as usize) % 2 == side)
                            .collect::<Vec<_>>()
                    })
                })
                .collect::<Vec<_>>(),
        )
    } else if k == 4 && lam == 3 && v == 9 {
        // Ordered differences of the two base blocks sum to3 for every
        // nonzero element of Z9; their translation orbits form18blocks.
        Some(
            [[0, 1, 2, 4], [0, 1, 4, 6]]
                .into_iter()
                .flat_map(|base| {
                    (0..9)
                        .map(move |shift| base.iter().map(|x| (x + shift) % 9).collect::<Vec<_>>())
                })
                .collect::<Vec<_>>(),
        )
    } else {
        None
    };
    if let Some(mut blocks) = direct {
        let mut rng = Rng(s);
        let mut labels: Vec<_> = (0..v).collect();
        rng.shuffle(&mut labels);
        for block in &mut blocks {
            for t in block.iter_mut() {
                *t = labels[*t];
            }
            rng.shuffle(block);
        }
        rng.shuffle(&mut blocks);
        return block_result(
            data,
            &blocks,
            None,
            s,
            if v == 8 {
                "balanced_incomplete_block_binary_affine"
            } else {
                "balanced_incomplete_block_cyclic_difference_family"
            },
            0,
        );
    }
    let mut candidates = Vec::new();
    combinations(v, k, 0, &mut Vec::new(), &mut candidates);
    let mut remaining = vec![vec![lam; v]; v];
    for (i, row) in remaining.iter_mut().enumerate() {
        row[i] = 0;
    }
    let mut search = Search {
        v,
        target: b,
        candidates,
        remaining,
        rep: vec![r; v],
        chosen: Vec::new(),
        nodes: 0,
        work: 0,
    };
    // Every possible BIBD can be relabeled so that its first block is0..k-1.
    // Fixing that block preserves existence while removing treatment-label symmetry.
    search.alter(0, true);
    search.chosen.push(search.candidates[0].clone());
    if !search.solve()? {
        return Err(
            "no_solution: exhaustive bounded BIBD search completed without a design".into(),
        );
    }
    let mut rng = Rng(s);
    let mut relabel: Vec<_> = (0..v).collect();
    rng.shuffle(&mut relabel);
    for c in &mut search.chosen {
        for t in c.iter_mut() {
            *t = relabel[*t];
        }
        rng.shuffle(c);
    }
    rng.shuffle(&mut search.chosen);
    block_result(
        data,
        &search.chosen,
        None,
        s,
        "balanced_incomplete_block_exact_search",
        search.nodes,
    )
}
fn resolvable(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    treatment_input(data, 49)?;
    if args.len() > 1 {
        return Err("resolvable design requires[seed=0]".into());
    }
    let v = data.len();
    let q = [2, 3, 5, 7]
        .into_iter()
        .find(|q| q * q == v)
        .ok_or("unsupported: affine-plane resolvable design requiresq^2treatments,q=2,3,5,7")?;
    let s = seed(args, 0)?;
    let mut rng = Rng(s);
    let mut labels: Vec<_> = (0..v).collect();
    rng.shuffle(&mut labels);
    let mut class_order: Vec<_> = (0..=q).collect();
    rng.shuffle(&mut class_order);
    let mut blocks = Vec::new();
    let mut classes = Vec::new();
    for c in class_order {
        let mut order: Vec<_> = (0..q).collect();
        rng.shuffle(&mut order);
        for b in order {
            let mut block: Vec<_> = if c == q {
                (0..q).map(|y| labels[b * q + y]).collect()
            } else {
                (0..q).map(|x| labels[x * q + (c * x + b) % q]).collect()
            };
            rng.shuffle(&mut block);
            blocks.push(block);
            classes.push(c);
        }
    }
    block_result(
        data,
        &blocks,
        Some(&classes),
        s,
        "affine_plane_resolvable_incomplete_block",
        0,
    )
}
fn restricted(data: &[Vec<f64>], args: &[f64], strip: bool) -> Result<String, String> {
    if data[0].len() != 3 || data.len() > 12 || args.len() > 2 {
        return Err(
            "restricted design requires<=12factor rows[low,high,role0/1] and[replicates=1,seed=0]"
                .into(),
        );
    }
    let mut counts = [0; 2];
    for r in data {
        if r[0] >= r[1] || ![0., 1.].contains(&r[2]) {
            return Err("factor low must be<high and role must be0or1".into());
        }
        counts[r[2] as usize] += 1;
    }
    if counts.contains(&0) {
        return Err("both randomization roles must contain at least one factor".into());
    }
    let reps = uint(args.first().copied().unwrap_or(1.), "replicates")?;
    if !(1..=256).contains(&reps) {
        return Err("replicates must be1..256".into());
    }
    let na = 1 << counts[0];
    let nb = 1 << counts[1];
    let n = reps * na * nb;
    rows_bound(n, data.len())?;
    let s = seed(args, 1)?;
    let mut rng = Rng(s);
    let mut reporder: Vec<_> = (0..reps).collect();
    rng.shuffle(&mut reporder);
    let mut d = Vec::new();
    let mut repids = Vec::new();
    let mut aid = Vec::new();
    let mut bid = Vec::new();
    let mut stage2 = Vec::new();
    let mut ats = Vec::new();
    let mut bts = Vec::new();
    for rep in reporder {
        let mut ao: Vec<_> = (0..na).collect();
        let mut bo: Vec<_> = (0..nb).collect();
        rng.shuffle(&mut ao);
        rng.shuffle(&mut bo);
        let mut cells = vec![vec![0; nb]; na];
        for a in ao {
            let mut subplot = bo.clone();
            if !strip {
                rng.shuffle(&mut subplot);
            }
            for b in subplot {
                let mut at = 0;
                let mut bt = 0;
                let mut row = Vec::new();
                for r in data {
                    let high = if r[2] == 0. {
                        let x = (a >> at) & 1;
                        at += 1;
                        x
                    } else {
                        let x = (b >> bt) & 1;
                        bt += 1;
                        x
                    };
                    row.push(r[high]);
                }
                cells[a][b] = d.len();
                d.push(row);
                repids.push(rep);
                aid.push(rep * na + a);
                bid.push(rep * nb + b);
                ats.push(a);
                bts.push(b);
            }
        }
        if strip {
            rng.shuffle(&mut bo);
            for b in bo {
                let mut aa: Vec<_> = (0..na).collect();
                rng.shuffle(&mut aa);
                for a in aa {
                    stage2.push(cells[a][b]);
                }
            }
        }
    }
    let df = if strip {
        vec![
            (reps - 1) * (na - 1),
            (reps - 1) * (nb - 1),
            (reps - 1) * (na - 1) * (nb - 1),
        ]
    } else {
        vec![(reps - 1) * (na - 1), na * (reps - 1) * (nb - 1)]
    };
    let extra = if strip {
        format!(
            ",\"row_strip_ids\":{},\"column_strip_ids\":{},\"stage1_order\":{},\"stage2_order\":{}",
            ints(&aid),
            ints(&bid),
            ints(&(0..n).collect::<Vec<_>>()),
            ints(&stage2)
        )
    } else {
        format!(
            ",\"whole_plot_ids\":{},\"whole_treatment_ids\":{},\"subplot_treatment_ids\":{}",
            ints(&aid),
            ints(&ats),
            ints(&bts)
        )
    };
    Ok(format!("{{\"version\":1,\"method\":\"{}\",\"seed\":{s},\"design\":{},\"replicate_ids\":{},\"run_order\":{},\"replication_error_df\":{},\"inference_ready\":{},\"diagnostic\":\"complete two-level factorial treatment structure with restricted randomization;response inference is not performed;replicates1has no replication error degrees of freedom\"{extra}}}",if strip{"classical_crossed_strip_plot"}else{"classical_nested_split_plot"},mat(&d),ints(&repids),ints(&(1..=n).collect::<Vec<_>>()),ints(&df),reps>1))
}
const L18: [[usize; 8]; 18] = [
    [1, 1, 1, 1, 1, 1, 1, 1],
    [1, 1, 2, 2, 2, 2, 2, 2],
    [1, 1, 3, 3, 3, 3, 3, 3],
    [1, 2, 1, 1, 2, 2, 3, 3],
    [1, 2, 2, 2, 3, 3, 1, 1],
    [1, 2, 3, 3, 1, 1, 2, 2],
    [1, 3, 1, 2, 1, 3, 2, 3],
    [1, 3, 2, 3, 2, 1, 3, 1],
    [1, 3, 3, 1, 3, 2, 1, 2],
    [2, 1, 1, 3, 3, 2, 2, 1],
    [2, 1, 2, 1, 1, 3, 3, 2],
    [2, 1, 3, 2, 2, 1, 1, 3],
    [2, 2, 1, 2, 3, 1, 3, 2],
    [2, 2, 2, 3, 1, 2, 1, 3],
    [2, 2, 3, 1, 2, 3, 2, 1],
    [2, 3, 1, 3, 2, 3, 1, 2],
    [2, 3, 2, 1, 3, 1, 2, 3],
    [2, 3, 3, 2, 1, 2, 3, 1],
];
fn contrast(level: usize, c: usize, l: usize) -> f64 {
    let x = if level <= c {
        1.
    } else if level == c + 1 {
        -((c + 1) as f64)
    } else {
        0.
    };
    x * (l as f64 / ((c + 1) * (c + 2)) as f64).sqrt()
}
fn orthogonal(data: &[Vec<f64>], args: &[f64]) -> Result<String, String> {
    if data[0].len() != 1 || data.len() > 8 || args.is_empty() || args.len() > 2 {
        return Err(
            "orthogonal array requires scalar2/3levels per factor and[array4/8/9/18,seed=0]".into(),
        );
    }
    let levels = data
        .iter()
        .map(|r| uint(r[0], "levels"))
        .collect::<Result<Vec<_>, _>>()?;
    if levels.iter().any(|v| !matches!(v, 2 | 3)) {
        return Err("only2/3level factors supported".into());
    }
    let n = uint(args[0], "array")?;
    let available = match n {
        4 => vec![2; 3],
        8 => vec![2; 7],
        9 => vec![3; 4],
        18 => vec![2, 3, 3, 3, 3, 3, 3, 3],
        _ => return Err("unsupported: catalog only L4,L8,L9,L18".into()),
    };
    let mut used = vec![false; available.len()];
    let mut cols = Vec::new();
    for &l in &levels {
        let col = available
            .iter()
            .enumerate()
            .position(|(i, v)| *v == l && !used[i])
            .ok_or("requested factor levels exceed selected array column capacities")?;
        used[col] = true;
        cols.push(col);
    }
    let raw: Vec<Vec<usize>> = match n {
        4 | 8 => (0..n)
            .map(|r| (1..n).map(|c| (r & c).count_ones() as usize % 2).collect())
            .collect(),
        9 => (0..3)
            .flat_map(|a| (0..3).map(move |b| vec![a, b, (a + b) % 3, (2 * a + b) % 3]))
            .collect(),
        18 => L18
            .iter()
            .map(|r| r.iter().map(|x| x - 1).collect())
            .collect(),
        _ => unreachable!(),
    };
    let mut coded: Vec<Vec<usize>> = raw
        .iter()
        .map(|r| cols.iter().map(|&i| r[i]).collect())
        .collect();
    let s = seed(args, 1)?;
    Rng(s).shuffle(&mut coded);
    // Verify strength2 rather than asserting an array's name guarantees it.
    for a in 0..levels.len() {
        for b in 0..a {
            let mut counts = vec![0; levels[a] * levels[b]];
            for r in &coded {
                counts[r[a] * levels[b] + r[b]] += 1;
            }
            if counts.iter().any(|x| *x != n / counts.len()) {
                return Err("internal orthogonal-array pair balance failure".into());
            }
        }
    }
    let terms: Vec<_> = levels
        .iter()
        .enumerate()
        .flat_map(|(i, l)| (0..l - 1).map(move |c| (i, c)))
        .collect();
    let x: Mat = coded
        .iter()
        .map(|r| {
            std::iter::once(1.)
                .chain(terms.iter().map(|&(f, c)| contrast(r[f], c, levels[f])))
                .collect()
        })
        .collect();
    let k = terms.len() + 1;
    let gram: Mat = (0..k)
        .map(|i| {
            (0..k)
                .map(|j| x.iter().map(|r| r[i] * r[j]).sum::<f64>() / n as f64)
                .collect()
        })
        .collect();
    let mut interactions = Vec::new();
    for (i, &(a, c)) in terms.iter().enumerate() {
        for (j, &(b, d)) in terms.iter().enumerate().take(i) {
            if a != b {
                interactions.push((i + 1, j + 1, vec![a, c, b, d]));
            }
        }
    }
    let alias: Mat = (1..k)
        .map(|i| {
            interactions
                .iter()
                .map(|(a, b, _)| x.iter().map(|r| r[i] * r[*a] * r[*b]).sum::<f64>() / n as f64)
                .collect()
        })
        .collect();
    let max_alias = alias
        .iter()
        .flatten()
        .map(|x| x.abs())
        .fold(0_f64, f64::max);
    let d: Mat = coded
        .iter()
        .map(|r| r.iter().map(|x| *x as f64).collect())
        .collect();
    let interaction_terms: Vec<_> = interactions.into_iter().map(|(_, _, t)| t).collect();
    Ok(format!("{{\"version\":1,\"method\":\"Taguchi_L{n}_orthogonal_array\",\"seed\":{s},\"design\":{},\"level_counts\":{},\"column_ids\":{},\"run_order\":{},\"strength2\":true,\"main_effect_gram\":{},\"main_effect_rank\":{k},\"main_effect_terms\":{},\"interaction_terms\":{},\"main_interaction_alias\":{},\"max_main_interaction_alias\":{max_alias},\"scope\":\"catalog array only;0based categorical levels;main-effect orthogonality does not imply interaction estimability;no S/N analysis or inner/outer crossing\"}}",mat(&d),ints(&levels),ints(&cols),ints(&(1..=n).collect::<Vec<_>>()),mat(&gram),imat(&terms.iter().map(|(a,b)|vec![*a,*b]).collect::<Vec<_>>()),imat(&interaction_terms),mat(&alias)))
}
