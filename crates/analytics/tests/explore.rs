use ducksmiles_analytics::explore;

// This tiny JSON reader verifies the public serialized result independently of
// production serialization, without adding a dependency to the extension.
#[derive(Debug)]
enum Json {
    Number(f64),
    Bool(bool),
    Array(Vec<Json>),
    Object(Vec<(String, Json)>),
}

impl Json {
    fn field(&self, name: &str) -> &Json {
        match self {
            Self::Object(fields) => &fields.iter().find(|(key, _)| key == name).unwrap().1,
            _ => panic!("expected object"),
        }
    }
    fn number(&self) -> f64 {
        match self {
            Self::Number(value) => *value,
            _ => panic!("expected number"),
        }
    }
    fn array(&self) -> &[Json] {
        match self {
            Self::Array(values) => values,
            _ => panic!("expected array"),
        }
    }
    fn vector(&self) -> Vec<f64> {
        self.array().iter().map(Json::number).collect()
    }
    fn matrix(&self) -> Vec<Vec<f64>> {
        self.array().iter().map(Json::vector).collect()
    }
    fn boolean(&self) -> bool {
        match self {
            Self::Bool(value) => *value,
            _ => panic!("expected bool"),
        }
    }
}

fn parse(text: &str) -> Json {
    fn value(bytes: &[u8], pos: &mut usize) -> Json {
        match bytes[*pos] {
            b'[' => {
                *pos += 1;
                let mut values = Vec::new();
                while bytes[*pos] != b']' {
                    if !values.is_empty() {
                        assert_eq!(bytes[*pos], b',');
                        *pos += 1;
                    }
                    values.push(value(bytes, pos));
                }
                *pos += 1;
                Json::Array(values)
            }
            b'{' => {
                *pos += 1;
                let mut fields = Vec::new();
                while bytes[*pos] != b'}' {
                    if !fields.is_empty() {
                        assert_eq!(bytes[*pos], b',');
                        *pos += 1;
                    }
                    assert_eq!(bytes[*pos], b'"');
                    *pos += 1;
                    let start = *pos;
                    while bytes[*pos] != b'"' {
                        *pos += 1;
                    }
                    let key = String::from_utf8(bytes[start..*pos].to_vec()).unwrap();
                    *pos += 1;
                    assert_eq!(bytes[*pos], b':');
                    *pos += 1;
                    fields.push((key, value(bytes, pos)));
                }
                *pos += 1;
                Json::Object(fields)
            }
            b't' => {
                assert_eq!(&bytes[*pos..*pos + 4], b"true");
                *pos += 4;
                Json::Bool(true)
            }
            b'f' => {
                assert_eq!(&bytes[*pos..*pos + 5], b"false");
                *pos += 5;
                Json::Bool(false)
            }
            _ => {
                let start = *pos;
                while *pos < bytes.len()
                    && matches!(bytes[*pos], b'0'..=b'9' | b'-' | b'+' | b'.' | b'e' | b'E')
                {
                    *pos += 1;
                }
                let number = std::str::from_utf8(&bytes[start..*pos])
                    .unwrap()
                    .parse::<f64>()
                    .unwrap();
                assert!(number.is_finite(), "JSON numbers must be finite");
                Json::Number(number)
            }
        }
    }
    let mut pos = 0;
    let result = value(text.as_bytes(), &mut pos);
    assert_eq!(pos, text.len());
    result
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-9 * (1.0 + expected.abs()),
        "actual {actual}, expected {expected}"
    );
}

// Catches incorrect covariance normalization or eigenvalue ordering. The
// centered covariance is [[10/3,8/3],[8/3,10/3]], with eigenvalues 6, 2/3.
#[test]
fn pca_recovers_hand_computed_covariance_spectrum() {
    let data = vec![
        vec![-2.0, -1.0],
        vec![-1.0, -2.0],
        vec![1.0, 2.0],
        vec![2.0, 1.0],
    ];
    let output = parse(&explore::run("pca", &data, &[2.0, 0.0]).unwrap());
    close(output.field("version").number(), 1.0);
    assert_eq!(output.field("mean").vector(), vec![0.0, 0.0]);
    assert_eq!(output.field("scale").vector(), vec![1.0, 1.0]);
    let variance = output.field("explained_variance").vector();
    close(variance[0], 6.0);
    close(variance[1], 2.0 / 3.0);
    let ratio = output.field("explained_variance_ratio").vector();
    close(ratio[0], 0.9);
    close(ratio[1], 0.1);
    let first = &output.field("components").matrix()[0];
    close(first[0], 1.0 / 2.0_f64.sqrt());
    close(first[1], 1.0 / 2.0_f64.sqrt());
}

// Catches transposed loadings, score projection mistakes, and loss of
// orthogonality on a three-dimensional covariance with unequal variances.
#[test]
fn pca_loadings_are_orthonormal_and_reconstruct_centered_rows() {
    let data = vec![
        vec![2.0, 1.0, 4.0],
        vec![3.0, 4.0, 0.0],
        vec![4.0, 2.0, 2.0],
        vec![8.0, 5.0, 1.0],
    ];
    let output = parse(&explore::run("pca", &data, &[3.0, 0.0]).unwrap());
    let components = output.field("components").matrix();
    let scores = output.field("scores").matrix();
    let mean = output.field("mean").vector();
    for (i, a) in components.iter().enumerate() {
        for (j, b) in components.iter().enumerate() {
            close(
                a.iter().zip(b).map(|(x, y)| x * y).sum(),
                if i == j { 1.0 } else { 0.0 },
            );
        }
    }
    for (row, score) in data.iter().zip(scores) {
        for col in 0..3 {
            let reconstructed: f64 = score.iter().zip(&components).map(|(s, c)| s * c[col]).sum();
            close(reconstructed + mean[col], row[col]);
        }
    }
}

// Catches division by zero on constants and use of population rather than
// sample standard deviation. The two variable columns are perfectly correlated.
#[test]
fn standardized_pca_handles_constant_descriptors() {
    let data = vec![
        vec![1.0, 10.0, 7.0],
        vec![2.0, 20.0, 7.0],
        vec![3.0, 30.0, 7.0],
    ];
    let output = parse(&explore::run("pca", &data, &[2.0, 1.0]).unwrap());
    assert_eq!(output.field("mean").vector(), vec![2.0, 20.0, 7.0]);
    assert_eq!(output.field("scale").vector(), vec![1.0, 10.0, 1.0]);
    let variance = output.field("explained_variance").vector();
    close(variance[0], 2.0);
    close(variance[1], 0.0);
    close(output.field("explained_variance_ratio").vector()[0], 1.0);
    let constant = parse(&explore::run("pca", &[vec![4.0], vec![4.0]], &[]).unwrap());
    assert_eq!(
        constant.field("explained_variance_ratio").vector(),
        vec![0.0]
    );
}

// Catches summation rounding a constant decimal column's mean away from its
// input value, which would otherwise invent variance during standardization.
#[test]
fn pca_does_not_invent_variance_for_decimal_constants() {
    let output = parse(&explore::run("pca", &vec![vec![0.1]; 3], &[1.0, 1.0]).unwrap());
    assert_eq!(output.field("mean").vector(), vec![0.1]);
    assert_eq!(output.field("scale").vector(), vec![1.0]);
    assert_eq!(output.field("explained_variance").vector(), vec![0.0]);
    assert_eq!(output.field("scores").matrix(), vec![vec![0.0]; 3]);
}

// Catches a rounded, unrepresentable feature mean silently inflating covariance
// when the feature offset is large relative to its observed variation.
#[test]
fn pca_rejects_features_that_cannot_be_reliably_centered_in_original_units() {
    let data = vec![
        vec![1e16],
        vec![1e16 + 2.0],
        vec![1e16 + 4.0],
        vec![1e16 + 6.0],
    ];
    assert!(explore::run("pca", &data, &[1.0, 0.0]).is_err());
    assert!(explore::run("pca", &data, &[1.0, 1.0]).is_err());
}

// Catches an absolute eigensolver threshold that erases valid small variances.
#[test]
fn pca_preserves_small_scale_covariance() {
    let data = vec![
        vec![-2e-12, -1e-12],
        vec![-1e-12, -2e-12],
        vec![1e-12, 2e-12],
        vec![2e-12, 1e-12],
    ];
    let output = parse(&explore::run("pca", &data, &[2.0]).unwrap());
    close(output.field("explained_variance").vector()[0] / 1e-24, 6.0);
    close(output.field("explained_variance_ratio").vector()[0], 0.9);
}

// Catches assigning by the wrong distance or failing to recompute centers.
#[test]
fn kmeans_separates_groups_with_correct_centers_and_inertia() {
    let data = vec![
        vec![0.0, 0.0],
        vec![0.0, 2.0],
        vec![10.0, 10.0],
        vec![10.0, 12.0],
        vec![-10.0, 10.0],
        vec![-10.0, 12.0],
    ];
    let output = parse(&explore::run("kmeans", &data, &[3.0, 100.0, 42.0]).unwrap());
    let labels = output.field("labels").vector();
    assert_eq!(labels[0], labels[1]);
    assert_eq!(labels[2], labels[3]);
    assert_eq!(labels[4], labels[5]);
    assert_ne!(labels[0], labels[2]);
    assert_ne!(labels[0], labels[4]);
    assert_ne!(labels[2], labels[4]);
    let centers = output.field("centers").matrix();
    for (label, expected) in [
        (labels[0], vec![0.0, 1.0]),
        (labels[2], vec![10.0, 11.0]),
        (labels[4], vec![-10.0, 11.0]),
    ] {
        assert_eq!(centers[label as usize], expected);
    }
    close(output.field("inertia").number(), 6.0);
    assert!(output.field("converged").boolean());
    assert!(output.field("iterations").number() <= 100.0);
}

// Catches nondeterministic initialization and empty-cluster NaNs on duplicates.
#[test]
fn kmeans_is_deterministic_and_retains_nonempty_duplicate_clusters() {
    let data = vec![vec![1.0], vec![1.0], vec![1.0], vec![1.0]];
    let first = explore::run("kmeans", &data, &[3.0, 10.0, 8.0]).unwrap();
    assert_eq!(
        first,
        explore::run("kmeans", &data, &[3.0, 10.0, 8.0]).unwrap()
    );
    let output = parse(&first);
    let labels = output.field("labels").vector();
    for cluster in 0..3 {
        assert!(labels.contains(&(cluster as f64)));
    }
    assert_eq!(output.field("centers").matrix(), vec![vec![1.0]; 3]);
    close(output.field("inertia").number(), 0.0);
    assert!(output.field("converged").boolean());
}

// Catches roundoff in constant-cluster means causing endless assignment and
// repair cycles for descriptor values that are not exact binary fractions.
#[test]
fn kmeans_converges_on_inexact_constant_values() {
    let output = parse(&explore::run("kmeans", &vec![vec![0.1]; 8], &[3.0, 10.0]).unwrap());
    assert!(output.field("converged").boolean());
    assert_eq!(output.field("centers").matrix(), vec![vec![0.1]; 3]);
    close(output.field("inertia").number(), 0.0);
}

// Catches invalid numeric/shape inputs escaping into JSON, and accidental casts
// that silently round integer parameters or accept unsupported options.
#[test]
fn exploration_rejects_invalid_shapes_values_and_parameters() {
    let valid = vec![vec![1.0, 2.0], vec![3.0, 4.0], vec![5.0, 6.0]];
    for bad in [
        vec![],
        vec![vec![]],
        vec![vec![1.0], vec![1.0, 2.0]],
        vec![vec![f64::NAN]],
        vec![vec![f64::INFINITY]],
    ] {
        assert!(explore::run("pca", &bad, &[]).is_err());
        assert!(explore::run("kmeans", &bad, &[]).is_err());
    }
    assert!(explore::run("pca", &[vec![1.0]], &[]).is_err());
    for params in [
        vec![0.0],
        vec![1.5],
        vec![3.0],
        vec![1.0, 2.0],
        vec![1.0, 0.0, 5.0],
        vec![f64::NAN],
    ] {
        assert!(
            explore::run("pca", &valid, &params).is_err(),
            "pca {params:?}"
        );
    }
    for params in [
        vec![0.0],
        vec![4.0],
        vec![1.5],
        vec![2.0, 0.0],
        vec![2.0, 1.2],
        vec![2.0, 10.0, -1.0],
        vec![2.0, 10.0, 1.2],
        vec![2.0, 10.0, 0.0, 1.0],
    ] {
        assert!(
            explore::run("kmeans", &valid, &params).is_err(),
            "kmeans {params:?}"
        );
    }
    assert!(explore::run("unknown", &valid, &[]).is_err());
    assert!(explore::run("pca", &[vec![f64::MAX], vec![-f64::MAX]], &[]).is_err());
    assert!(explore::run("kmeans", &[vec![f64::MAX], vec![-f64::MAX]], &[1.0]).is_err());
}

// Catches unbounded feature/cell allocations before expensive matrix work.
#[test]
fn exploration_rejects_excessive_dimensions_and_work() {
    let too_wide = vec![vec![0.0; 257]; 2];
    assert!(explore::run("pca", &too_wide, &[]).is_err());
    assert!(explore::run("kmeans", &too_wide, &[]).is_err());
    let too_many_cells = vec![vec![0.0; 256]; 4097];
    assert!(explore::run("pca", &too_many_cells, &[]).is_err());
    assert!(explore::run("kmeans", &too_many_cells, &[]).is_err());
    assert!(explore::run("kmeans", &[vec![1.0], vec![2.0]], &[1.0, 10001.0]).is_err());
}
