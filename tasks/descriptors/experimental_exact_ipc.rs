//! Experimental mathematical IPC; not included by the shipped crate/catalog.
//! Independent expected coefficients use determinants/interpolation or formulas.

struct Exact {
    coefficients: Vec<String>,
    ipc: f64,
    average: f64,
}

#[derive(Clone, Default)]
struct Integer {
    negative: bool,
    // Little-endian base10^9 limbs; empty means canonical unsigned zero.
    digits: Vec<u32>,
}

impl Integer {
    const BASE: u64 = 1_000_000_000;

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
            assert_eq!(borrow, 0, "magnitude subtraction invariant");
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

fn compute(adjacency: &[Vec<usize>]) -> Option<Exact> {
    let n = adjacency.len();
    if n == 0 || n > 256 {
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
    if edges.iter().filter(|&&e| e).count() > 1024 {
        return None;
    }
    for i in 0..n {
        for j in 0..i {
            if edges[i * n + j] != edges[j * n + i] {
                return None;
            }
        }
    }
    let mut current = vec![Integer::default(); n * n];
    let mut next = current.clone();
    for (i, &edge) in edges.iter().enumerate() {
        if edge {
            current[i].digits.push(1);
        }
    }
    let mut coefficients = vec!["1".to_owned()];
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
        coefficients.push(coefficient.decimal());
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
    let values: Vec<f64> = coefficients
        .iter()
        .map(|c| c.parse::<f64>().map(f64::abs))
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
    let ipc = total * average;
    if !ipc.is_finite() || !average.is_finite() {
        return None;
    }
    Some(Exact {
        coefficients,
        ipc,
        average,
    })
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
        for line in include_str!("exact-ipc-reference.tsv")
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
            let actual =
                compute(&adjacency).unwrap_or_else(|| panic!("no exact result for {}", fields[0]));
            assert_eq!(
                actual.coefficients.join(","),
                fields[3],
                "{} signed coefficients",
                fields[0]
            );
            let ipc: f64 = fields[4].parse().unwrap();
            let average: f64 = fields[5].parse().unwrap();
            assert!(
                (actual.ipc - ipc).abs() <= 1e-12 * ipc.abs().max(1.),
                "{} IPC {} vs {ipc}",
                fields[0],
                actual.ipc
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
            assert!(compute(&adjacency).is_none());
        }
    }
}
