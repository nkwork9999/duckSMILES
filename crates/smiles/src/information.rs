//! Graph information from exact integer characteristic-polynomial coefficients.
//! Mathematical IPC/AvgIPC definitions; pinned RDKit/NumPy can be numerically
//! unstable on large graphs. Entropy alone uses f64. No3D or runtime dependency.
//! Independent tests use determinants/interpolation, closed forms and Decimal85.

mod factor;

pub(crate) struct Values {
    pub(crate) weighted: f64,
    pub(crate) average: f64,
}

#[derive(Clone, Default)]
struct Integer {
    negative: bool,
    // Little-endian base10^9 limbs; empty means canonical unsigned zero.
    digits: Vec<u32>,
}

impl Integer {
    const BASE: u64 = 1_000_000_000;

    fn multiply(&self, other: &Self) -> Self {
        if self.digits.is_empty() || other.digits.is_empty() {
            return Self::default();
        }
        let mut digits = vec![0u32; self.digits.len() + other.digits.len()];
        for (i, &left) in self.digits.iter().enumerate() {
            let mut carry = 0u64;
            for (j, &right) in other.digits.iter().enumerate() {
                let value = u64::from(left) * u64::from(right) + u64::from(digits[i + j]) + carry;
                digits[i + j] = (value % Self::BASE) as u32;
                carry = value / Self::BASE;
            }
            // This position has not been touched by a previous outer row.
            digits[i + other.digits.len()] = carry as u32;
        }
        let mut result = Self {
            negative: self.negative != other.negative,
            digits,
        };
        result.normalize();
        result
    }

    fn clear(&mut self) {
        self.digits.clear();
        self.negative = false;
    }

    fn assign(&mut self, other: &Self) {
        self.digits.clone_from(&other.digits);
        self.negative = other.negative;
    }

    fn normalize(&mut self) {
        while self.digits.last() == Some(&0) {
            self.digits.pop();
        }
        if self.digits.is_empty() {
            self.negative = false;
        }
    }

    fn add_assign(&mut self, other: &Self) {
        if other.digits.is_empty() {
            return;
        }
        if self.digits.is_empty() {
            self.assign(other);
            return;
        }
        if self.negative == other.negative {
            let length = self.digits.len().max(other.digits.len());
            self.digits.resize(length, 0);
            let mut carry = 0;
            for i in 0..length {
                let sum = u64::from(self.digits[i])
                    + u64::from(other.digits.get(i).copied().unwrap_or(0))
                    + carry;
                self.digits[i] = (sum % Self::BASE) as u32;
                carry = sum / Self::BASE;
            }
            if carry != 0 {
                self.digits.push(carry as u32);
            }
        } else {
            let ordering = self
                .digits
                .len()
                .cmp(&other.digits.len())
                .then_with(|| self.digits.iter().rev().cmp(other.digits.iter().rev()));
            if ordering.is_eq() {
                self.clear();
                return;
            }
            let other_larger = ordering.is_lt();
            let old_length = self.digits.len();
            let length = old_length.max(other.digits.len());
            self.digits.resize(length, 0);
            let mut borrow = 0i64;
            for i in 0..length {
                let left = i64::from(self.digits[i]);
                let right = i64::from(other.digits.get(i).copied().unwrap_or(0));
                let mut difference = if other_larger {
                    right - left
                } else {
                    left - right
                } - borrow;
                borrow = i64::from(difference < 0);
                if difference < 0 {
                    difference += Self::BASE as i64;
                }
                self.digits[i] = difference as u32;
            }
            debug_assert_eq!(borrow, 0, "magnitude subtraction invariant");
            if other_larger {
                self.negative = other.negative;
            }
            self.normalize();
        }
    }

    fn divide_exact(&mut self, divisor: u32) -> bool {
        let mut remainder = 0u64;
        for digit in self.digits.iter_mut().rev() {
            let value = remainder * Self::BASE + u64::from(*digit);
            *digit = (value / u64::from(divisor)) as u32;
            remainder = value % u64::from(divisor);
        }
        self.normalize();
        remainder == 0
    }

    fn negate(&mut self) {
        if !self.digits.is_empty() {
            self.negative = !self.negative;
        }
    }

    fn decimal(&self) -> String {
        use std::fmt::Write;
        let Some(last) = self.digits.last() else {
            return "0".into();
        };
        let mut result = if self.negative {
            "-".to_owned()
        } else {
            String::new()
        };
        write!(result, "{last}").unwrap();
        for digit in self.digits.iter().rev().skip(1) {
            write!(result, "{digit:09}").unwrap();
        }
        result
    }
}

fn coefficients(adjacency: &[Vec<usize>]) -> Option<Vec<Integer>> {
    if adjacency.len() <= 32 {
        return dense_coefficients(adjacency);
    }
    // Validate before graph traversal. Factoring changes the integer algorithm,
    // never the domain, signed coefficients or subsequent entropy reduction.
    validated_edges(adjacency)?;
    factor::coefficients(adjacency, 32).or_else(|| dense_coefficients(adjacency))
}

fn validated_edges(adjacency: &[Vec<usize>]) -> Option<Vec<bool>> {
    let n = adjacency.len();
    if n == 0 || n > crate::workbench::MAX_ATOMS {
        return None;
    }
    let mut edges = vec![false; n * n];
    for (i, row) in adjacency.iter().enumerate() {
        for &j in row {
            if j >= n || i == j || edges[i * n + j] {
                return None;
            }
            edges[i * n + j] = true;
        }
    }
    if edges.iter().filter(|&&e| e).count() > 2 * crate::workbench::MAX_BONDS {
        return None;
    }
    for i in 0..n {
        for j in 0..i {
            if edges[i * n + j] != edges[j * n + i] {
                return None;
            }
        }
    }
    Some(edges)
}

fn dense_coefficients(adjacency: &[Vec<usize>]) -> Option<Vec<Integer>> {
    let n = adjacency.len();
    let edges = validated_edges(adjacency)?;
    let mut current = vec![Integer::default(); n * n];
    let mut next = current.clone();
    for (i, &edge) in edges.iter().enumerate() {
        if edge {
            current[i].digits.push(1);
        }
    }
    let mut coefficients = Vec::with_capacity(n + 1);
    coefficients.push(Integer {
        negative: false,
        digits: vec![1],
    });
    // Faddeev: C=A B_(k-1), c_k=-tr(C)/k, B_k=C+c_k I.
    // Every division is exact over integers. All matrices are polynomials in
    // symmetric A, so compute one triangle and copy the exact integer values.
    for k in 1..=n {
        let mut coefficient = Integer::default();
        for i in 0..n {
            coefficient.add_assign(&current[i * n + i]);
        }
        if !coefficient.divide_exact(k as u32) {
            return None;
        }
        coefficient.negate();
        coefficients.push(coefficient.clone());
        for i in 0..n {
            current[i * n + i].add_assign(&coefficient);
        }
        if k == n {
            if current.iter().any(|c| !c.digits.is_empty()) {
                return None;
            }
            break;
        }
        for i in 0..n {
            for j in 0..=i {
                let index = i * n + j;
                let sum = &mut next[index];
                sum.clear();
                for &neighbor in &adjacency[i] {
                    sum.add_assign(&current[neighbor * n + j]);
                }
                if i != j {
                    let (upper, lower) = next.split_at_mut(index);
                    upper[j * n + i].assign(&lower[0]);
                }
            }
        }
        std::mem::swap(&mut current, &mut next);
    }
    Some(coefficients)
}

fn entropy(coefficients: &[Integer]) -> Option<Values> {
    let values: Vec<f64> = coefficients
        .iter()
        .map(|c| c.decimal().parse::<f64>().map(f64::abs))
        .collect::<Result<_, _>>()
        .ok()?;
    let total: f64 = values.iter().sum();
    let average = -values
        .iter()
        .filter(|&&v| v != 0.)
        .map(|&v| {
            let p = v / total;
            p * p.ln()
        })
        .sum::<f64>()
        / std::f64::consts::LN_2;
    let weighted = total * average;
    if !weighted.is_finite() || !average.is_finite() {
        return None;
    }
    Some(Values { weighted, average })
}

pub(crate) fn calculate(mol: &crate::parser::Molecule) -> Option<Values> {
    let adjacency: Vec<Vec<usize>> = mol
        .adjacency()
        .into_iter()
        .map(|row| row.into_iter().map(|(neighbor, _)| neighbor).collect())
        .collect();
    entropy(&coefficients(&adjacency)?)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn from_decimal(text: &str) -> Integer {
        let magnitude = text.strip_prefix('-').unwrap_or(text);
        let mut value = Integer {
            negative: text.starts_with('-'),
            digits: magnitude
                .as_bytes()
                .rchunks(9)
                .map(|c| std::str::from_utf8(c).unwrap().parse().unwrap())
                .collect(),
        };
        value.normalize();
        value
    }

    #[test]
    fn signed_addition_and_division_agree_with_independent_i128_arithmetic() {
        let mut state = 0x96a5_b34c_17de_823fu64;
        let next = |state: &mut u64| {
            *state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            *state
        };
        for _ in 0..20_000 {
            let mut number = || {
                let magnitude = (((next(&mut state) as u128) << 64) | next(&mut state) as u128)
                    & ((1u128 << 120) - 1);
                if (next(&mut state) >> 32) & 1 == 0 {
                    magnitude as i128
                } else {
                    -(magnitude as i128)
                }
            };
            let (a, b) = (number(), number());
            let expected = a + b;
            let mut actual = from_decimal(&a.to_string());
            actual.add_assign(&from_decimal(&b.to_string()));
            assert_eq!(actual.decimal(), expected.to_string());
            let divisor = (next(&mut state) % 256 + 1) as u32;
            assert_eq!(
                actual.divide_exact(divisor),
                expected % i128::from(divisor) == 0
            );
            let quotient = expected / i128::from(divisor);
            assert_eq!(actual.decimal(), quotient.to_string());
            actual.negate();
            assert_eq!(actual.decimal(), (-quotient).to_string());
        }
    }

    #[test]
    fn factor_integer_products_match_signed_i128_and_large_closed_form() {
        let mut state = 0xabcd_6378_2ee5_19b4u64;
        let mut number = || {
            state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (state as i64 >> 1) as i128
        };
        for _ in 0..5_000 {
            let (a, b) = (number(), number());
            let actual = from_decimal(&a.to_string()).multiply(&from_decimal(&b.to_string()));
            assert_eq!(actual.decimal(), (a * b).to_string());
        }
        let nines = "9".repeat(270);
        let squared = from_decimal(&nines).multiply(&from_decimal(&nines));
        assert_eq!(
            squared.decimal(),
            "9".repeat(269) + "8" + &"0".repeat(269) + "1"
        );
        let zero = from_decimal("0").multiply(&from_decimal("-999999999999999999999"));
        assert_eq!(zero.decimal(), "0");
        assert!(!zero.negative);
    }

    #[test]
    fn factor_cut_identity_matches_all_58_independent_coefficient_vectors() {
        let mut count = 0;
        for line in include_str!("../tests/fixtures/descriptor_ipc_exact_graphs.tsv")
            .lines()
            .filter(|l| !l.starts_with('#'))
        {
            let fields: Vec<_> = line.split('\t').collect();
            let n: usize = fields[1].parse().unwrap();
            let mut adjacency = vec![Vec::new(); n];
            for edge in fields[2].split(';').filter(|s| !s.is_empty()) {
                let pair: Vec<usize> = edge.split(',').map(|s| s.parse().unwrap()).collect();
                adjacency[pair[0]].push(pair[1]);
                adjacency[pair[1]].push(pair[0]);
            }
            // Force tiny leaves to exercise the algebra. This independent
            // coefficient test has its own larger work allowance; production
            // uses cutoff32 and the strict 4n+64 limit with dense fallback.
            let actual = super::factor::coefficients_limited(&adjacency, 3, 65_536)
                .expect("factorization must return coefficients");
            assert_eq!(
                actual
                    .iter()
                    .map(Integer::decimal)
                    .collect::<Vec<_>>()
                    .join(","),
                fields[3],
                "{}",
                fields[0]
            );
            count += 1;
        }
        assert_eq!(count, 58);
    }

    #[test]
    fn factor_random_block_graphs_match_the_independent_dense_recurrence() {
        fn next(state: &mut u64) -> usize {
            *state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            (*state >> 32) as usize
        }
        let mut state = 0xc636_48f2_136e_5a97u64;
        for case in 0..128 {
            let mut graph = vec![Vec::new()];
            for _ in 0..3 + next(&mut state) % 6 {
                let anchor = next(&mut state) % graph.len();
                let size = 2 + next(&mut state) % 4;
                let mut block = vec![anchor];
                for _ in 1..size {
                    block.push(graph.len());
                    graph.push(Vec::new());
                }
                for i in 0..block.len() {
                    for j in i + 1..block.len() {
                        graph[block[i]].push(block[j]);
                        graph[block[j]].push(block[i]);
                    }
                }
            }
            if case % 3 == 0 {
                graph.push(Vec::new());
            }
            let direct = dense_coefficients(&graph).unwrap();
            let factored = super::factor::coefficients(&graph, 3).unwrap();
            assert_eq!(
                factored.iter().map(Integer::decimal).collect::<Vec<_>>(),
                direct.iter().map(Integer::decimal).collect::<Vec<_>>(),
                "case{case}"
            );
        }
    }

    #[test]
    fn factor_largest_original_input_matches_every_independent_signed_coefficient() {
        let line = include_str!("../tests/fixtures/information_largest_coefficients.tsv")
            .lines()
            .find(|line| !line.starts_with('#'))
            .unwrap();
        let (source, expected) = line.split_once('\t').unwrap();
        let (mol, _) = crate::workbench::prepare(source).unwrap();
        let adjacency: Vec<Vec<usize>> = mol
            .adjacency()
            .iter()
            .map(|row| row.iter().map(|&(j, _)| j).collect())
            .collect();
        assert_eq!(adjacency.len(), 598);
        let actual = super::factor::coefficients(&adjacency, 32).unwrap();
        assert_eq!(
            actual
                .iter()
                .map(Integer::decimal)
                .collect::<Vec<_>>()
                .join(","),
            expected
        );
    }

    #[test]
    fn carries_borrows_and_zero_remain_correct_beyond_machine_integer_width() {
        let nines = "9".repeat(270);
        let mut a = from_decimal(&nines);
        a.add_assign(&from_decimal("1"));
        assert_eq!(a.decimal(), "1".to_owned() + &"0".repeat(270));
        assert!(a.divide_exact(2));
        let half = "5".to_owned() + &"0".repeat(269);
        assert_eq!(a.decimal(), half);
        a.add_assign(&from_decimal(&("-".to_owned() + &half)));
        assert_eq!(a.decimal(), "0");
        a.negate();
        assert_eq!(a.decimal(), "0");
        let mut b = from_decimal(&nines);
        assert!(!b.divide_exact(2));
        assert_eq!(b.decimal(), "4".to_owned() + &"9".repeat(269));
        let mut c = from_decimal(&nines);
        assert!(c.divide_exact(9));
        assert_eq!(c.decimal(), "1".repeat(270));
    }

    #[test]
    fn integer_coefficients_and_entropy_match_independent_exact_oracles() {
        let mut count = 0;
        for line in include_str!("../tests/fixtures/descriptor_ipc_exact_graphs.tsv")
            .lines()
            .filter(|l| !l.starts_with('#'))
        {
            let fields: Vec<_> = line.split('\t').collect();
            let n: usize = fields[1].parse().unwrap();
            let mut adjacency = vec![Vec::new(); n];
            for edge in fields[2].split(';').filter(|s| !s.is_empty()) {
                let pair: Vec<usize> = edge.split(',').map(|s| s.parse().unwrap()).collect();
                adjacency[pair[0]].push(pair[1]);
                adjacency[pair[1]].push(pair[0]);
            }
            let signed = coefficients(&adjacency)
                .unwrap_or_else(|| panic!("no exact result for {}", fields[0]));
            assert_eq!(
                signed
                    .iter()
                    .map(Integer::decimal)
                    .collect::<Vec<_>>()
                    .join(","),
                fields[3],
                "{} signed coefficients",
                fields[0]
            );
            let actual = entropy(&signed).unwrap();
            let ipc: f64 = fields[4].parse().unwrap();
            let average: f64 = fields[5].parse().unwrap();
            assert!(
                (actual.weighted - ipc).abs() <= 1e-12 * ipc.abs().max(1.),
                "{} IPC {} vs {ipc}",
                fields[0],
                actual.weighted
            );
            assert!(
                (actual.average - average).abs() <= 1e-12,
                "{} AvgIPC {} vs {average}",
                fields[0],
                actual.average
            );
            count += 1;
        }
        assert_eq!(count, 58);
    }

    #[test]
    fn invalid_adjacencies_have_no_exact_result() {
        for adjacency in [
            vec![],
            vec![vec![0]],
            vec![vec![1], vec![]],
            vec![vec![1, 1], vec![0, 0]],
            vec![vec![2], vec![0]],
        ] {
            assert!(coefficients(&adjacency).is_none());
        }
        let mut large_invalid = vec![Vec::new(); 33];
        large_invalid[0] = vec![33];
        assert!(coefficients(&large_invalid).is_none());
        large_invalid[0] = vec![1];
        assert!(coefficients(&large_invalid).is_none());
        large_invalid[1] = vec![0, 0];
        assert!(coefficients(&large_invalid).is_none());
        assert!(coefficients(&vec![Vec::new(); 641]).is_none());
    }
}
