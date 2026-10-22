//! Bounded static SVG reports. Integer frequencies are virtual replications;
//! hidden rows remain in summaries. No analytic weights or automatic bins.
use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write;

const MAX_ROWS: usize = 100_000;
const MAX_GROUPS: usize = 16;
const MAX_N: u64 = 1_000_000;
const MAX_MARKS: usize = 20_000;
const MAX_SVG: usize = 2_097_152;
const MAX_ID: f64 = 9_007_199_254_740_991.;
type Result<T> = std::result::Result<T, String>;

#[derive(Clone)]
struct Row {
    value: f64,
    freq: u64,
    hidden: bool,
    id: u64,
}
#[derive(Default)]
struct Group {
    id: i64,
    source: usize,
    excluded: usize,
    missing: usize,
    zero: usize,
    hidden: usize,
    n: u64,
    rows: Vec<Row>,
}
#[derive(Clone, Copy)]
struct Panel {
    left: f64,
    top: f64,
    width: f64,
    height: f64,
}
impl Panel {
    fn json(self) -> String {
        format!(
            "{{\"left\":{},\"top\":{},\"width\":{},\"height\":{}}}",
            self.left, self.top, self.width, self.height
        )
    }
}
fn checked(x: f64) -> Result<f64> {
    if x.is_finite() {
        Ok(if x == 0. { 0. } else { x })
    } else {
        Err("unrepresentable report metric or coordinate".into())
    }
}
fn integer(x: f64, lo: f64, hi: f64) -> Result<f64> {
    if x.is_finite() && x.fract() == 0. && x >= lo && x <= hi {
        Ok(x)
    } else {
        Err("invalid report integer or flag".into())
    }
}
fn arr(v: &[f64]) -> String {
    format!(
        "[{}]",
        v.iter()
            .map(|x| x.to_string())
            .collect::<Vec<_>>()
            .join(",")
    )
}
fn escaped(s: &str) -> String {
    let mut out = String::new();
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\r' => out.push_str("\\r"),
            '\t' => out.push_str("\\t"),
            _ => out.push(c),
        }
    }
    out
}
fn label(x: f64) -> String {
    if x == 0. {
        "0".into()
    } else if x.abs() >= 1e6 || x.abs() < 0.001 {
        format!("{x:.5e}")
    } else {
        let s = format!("{x:.5}");
        s.trim_end_matches('0').trim_end_matches('.').to_string()
    }
}
fn parse(data: &[Vec<f64>]) -> Result<Vec<Group>> {
    if data.is_empty() || data.len() > MAX_ROWS {
        return Err("report requires 1..100000 source rows".into());
    }
    let mut groups: BTreeMap<i64, Group> = BTreeMap::new();
    let mut ids = BTreeSet::new();
    let mut total = 0u64;
    for r in data {
        if r.len() != 7 || r.iter().any(|x| !x.is_finite()) {
            return Err("report matrix requires exactly seven finite columns".into());
        }
        let id = integer(r[6], 0., MAX_ID)? as u64;
        if !ids.insert(id) {
            return Err("duplicate report row_id".into());
        }
        let group = integer(r[1], -MAX_ID, MAX_ID)? as i64;
        let freq = integer(r[2], 0., MAX_N as f64)? as u64;
        let excluded = integer(r[3], 0., 1.)? == 1.;
        let hidden = integer(r[4], 0., 1.)? == 1.;
        let observed = integer(r[5], 0., 1.)? == 1.;
        total = total.checked_add(freq).ok_or("frequency total overflow")?;
        if total > MAX_N {
            return Err("report total frequency exceeds 1000000".into());
        }
        let g = groups.entry(group).or_insert_with(|| Group {
            id: group,
            ..Group::default()
        });
        g.source += 1;
        g.excluded += usize::from(excluded);
        g.missing += usize::from(!observed);
        g.zero += usize::from(freq == 0);
        g.hidden += usize::from(hidden);
        if observed && !excluded && freq > 0 {
            g.n += freq;
            g.rows.push(Row {
                value: r[0],
                freq,
                hidden,
                id,
            });
        }
        if groups.len() > MAX_GROUPS {
            return Err("report exceeds 16 source groups".into());
        }
    }
    let mut groups: Vec<_> = groups.into_values().collect();
    if groups.iter().all(|g| g.n == 0) {
        return Err("report requires at least one selected observation".into());
    }
    for g in &mut groups {
        g.rows
            .sort_by(|a, b| a.value.total_cmp(&b.value).then(a.id.cmp(&b.id)));
    }
    Ok(groups)
}
fn group_prefix(g: &Group) -> String {
    format!("{{\"group_id\":{},\"status\":\"{}\",\"source_rows\":{},\"selected_rows\":{},\"N\":{},\"excluded_rows\":{},\"missing_rows\":{},\"zero_frequency_rows\":{},\"hidden_rows\":{}",g.id,if g.n==0 {"empty"} else {"ok"},g.source,g.rows.len(),g.n,g.excluded,g.missing,g.zero,g.hidden)
}
fn order(g: &Group, k: u64) -> f64 {
    let mut sum = 0;
    for r in &g.rows {
        sum += r.freq;
        if sum >= k {
            return r.value;
        }
    }
    g.rows.last().unwrap().value
}
fn quantile(g: &Group, p: f64) -> Result<f64> {
    let rank = (g.n + 1) as f64 * p;
    if rank <= 1. {
        return Ok(g.rows[0].value);
    }
    if rank >= g.n as f64 {
        return Ok(g.rows.last().unwrap().value);
    }
    let i = rank.floor() as u64;
    let f = rank - i as f64;
    let a = order(g, i);
    let b = order(g, i + 1);
    // Convex interpolation avoids overflowing b-a when opposite signs are large.
    checked((1. - f) * a + f * b)
}
fn domain(lo: f64, hi: f64) -> Result<[f64; 2]> {
    if lo < hi {
        let range = checked(hi - lo)?;
        if range <= 0. {
            return Err("unrepresentable value domain".into());
        }
        return Ok([lo, hi]);
    }
    let pad = (lo.abs() * 0.05).max(1.);
    let a = checked(lo - pad)?;
    let b = checked(hi + pad)?;
    if !(a < lo && b > hi) {
        return Err("cannot represent constant-value padding".into());
    }
    checked(b - a)?;
    Ok([a, b])
}
fn transform(p: Panel, d: [f64; 2], vertical: bool) -> Result<(f64, f64)> {
    let span = checked(d[1] - d[0])?;
    let scale = checked(if vertical {
        -p.height / span
    } else {
        p.width / span
    })?;
    let offset = checked(if vertical {
        p.top + p.height - scale * d[0]
    } else {
        p.left - scale * d[0]
    })?;
    if scale == 0. {
        return Err("value transform is not representable".into());
    }
    Ok((scale, offset))
}
fn map(p: Panel, d: [f64; 2], v: f64, vertical: bool) -> Result<f64> {
    let fraction = checked((v - d[0]) / (d[1] - d[0]))?;
    checked(if vertical {
        p.top + p.height * (1. - fraction)
    } else {
        p.left + p.width * fraction
    })
}
fn text(svg: &mut String, x: f64, y: f64, anchor: &str, content: &str) {
    write!(
        svg,
        "<text x=\"{x}\" y=\"{y}\" text-anchor=\"{anchor}\">{content}</text>"
    )
    .unwrap();
}
fn line(svg: &mut String, x1: f64, y1: f64, x2: f64, y2: f64, attributes: &str) {
    write!(
        svg,
        "<line x1=\"{x1}\" y1=\"{y1}\" x2=\"{x2}\" y2=\"{y2}\" {attributes}/>"
    )
    .unwrap();
}
fn start(width: f64, height: f64, title: &str, desc: &str) -> String {
    format!("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{width}\" height=\"{height}\" viewBox=\"0 0 {width} {height}\" role=\"img\" aria-labelledby=\"report-title report-description\"><title id=\"report-title\">{title}</title><desc id=\"report-description\">{desc}</desc><style>text{{font-family:sans-serif;font-size:12px;fill:#17212b}}line{{stroke:#17212b;stroke-width:1.5}}.axis{{stroke:#17212b;stroke-width:1}}.box{{fill:white;stroke:#17212b;stroke-width:2}}.bar{{fill:#2877a8;stroke:#17212b;stroke-width:0.6}}.median{{stroke:#111;stroke-width:3}}</style><rect width=\"{width}\" height=\"{height}\" fill=\"white\"/>")
}
fn precise_label(value: f64) -> String {
    let plain = value.to_string();
    let scientific = format!("{value:e}");
    if plain.len() <= scientific.len() {
        plain
    } else {
        scientific
    }
}
fn tick_labels(d: [f64; 2], ticks: usize) -> Result<Vec<(f64, String)>> {
    let mut result = Vec::new();
    for i in 0..=ticks {
        let f = i as f64 / ticks as f64;
        let value = checked(if i == 0 {
            d[0]
        } else if i == ticks {
            d[1]
        } else {
            (1. - f) * d[0] + f * d[1]
        })?;
        result.push((value, label(value)));
    }
    if result.windows(2).any(|pair| pair[0].0 == pair[1].0) {
        return Err("axis has fewer distinct representable ticks".into());
    }
    if result.windows(2).any(|pair| pair[0].1 == pair[1].1) {
        for (value, text) in &mut result {
            *text = precise_label(*value);
        }
    }
    Ok(result)
}
fn axis_label_width(d: [f64; 2]) -> Result<f64> {
    let labels = tick_labels(d, 4)
        .or_else(|_| tick_labels(d, 2))
        .or_else(|_| tick_labels(d, 1))?;
    Ok(labels.iter().map(|(_, s)| s.len()).max().unwrap() as f64 * 7.)
}
fn axes(svg: &mut String, p: Panel, d: [f64; 2], vertical: bool, axis_name: &str) -> Result<()> {
    let length = if vertical { p.height } else { p.width };
    let mut chosen = None;
    for ticks in [4, 2, 1] {
        if let Ok(labels) = tick_labels(d, ticks) {
            let text_width = labels.iter().map(|(_, s)| s.len()).max().unwrap() as f64 * 7.;
            if length / ticks as f64 >= if vertical { 20. } else { text_width + 10. } {
                chosen = Some(labels);
                break;
            }
        }
    }
    let labels = chosen.ok_or("axis labels cannot fit readably")?;
    let role = if axis_name == "Value" {
        "value-tick-label"
    } else {
        "response-tick-label"
    };
    if vertical {
        line(
            svg,
            p.left,
            p.top,
            p.left,
            p.top + p.height,
            "class=\"axis\"",
        );
    } else {
        line(
            svg,
            p.left,
            p.top + p.height,
            p.left + p.width,
            p.top + p.height,
            "class=\"axis\"",
        );
    }
    for (value, label) in labels {
        let v = map(p, d, value, vertical)?;
        let (x, y, anchor) = if vertical {
            line(svg, p.left - 4., v, p.left, v, "class=\"axis\"");
            (p.left - 7., v + 4., "end")
        } else {
            line(
                svg,
                v,
                p.top + p.height,
                v,
                p.top + p.height + 4.,
                "class=\"axis\"",
            );
            (v, p.top + p.height + 18., "middle")
        };
        write!(svg,"<text data-role=\"{role}\" data-value=\"{value}\" x=\"{x}\" y=\"{y}\" text-anchor=\"{anchor}\">{label}</text>").unwrap();
    }
    if vertical {
        text(svg, p.left, p.top - 12., "start", axis_name);
    } else {
        text(
            svg,
            p.left + p.width / 2.,
            p.top + p.height + 36.,
            "middle",
            axis_name,
        );
    }
    Ok(())
}

fn final_json(
    op: &str,
    groups: &[Group],
    canvas: [f64; 2],
    policy: &str,
    details: &str,
    svg: String,
    options: (&str, &str),
) -> Result<String> {
    let [width, height] = canvas;
    let (implemented, remaining) = options;
    if svg.len() > MAX_SVG {
        return Err("report SVG exceeds 2097152 bytes".into());
    }
    let source: usize = groups.iter().map(|g| g.source).sum();
    let selected: usize = groups.iter().map(|g| g.rows.len()).sum();
    let n: u64 = groups.iter().map(|g| g.n).sum();
    Ok(format!("{{\"version\":1,\"schema\":\"ducksmiles-cpu-plot-v1\",\"operation\":\"{op}\",\"media_type\":\"image/svg+xml\",\"width\":{width},\"height\":{height},\"source_rows\":{source},\"selected_rows\":{selected},\"N\":{n},\"selection_policy\":\"observed=1 and excluded=0 and frequency>0; hidden retains summaries and suppresses individual marks; row-state counts overlap\",\"policy\":\"{policy}\",\"limits\":{{\"source_rows\":100000,\"groups\":16,\"sum_frequency_per_call\":1000000,\"histogram_bins\":512,\"rendered_individual_markers\":20000,\"svg_utf8_bytes\":2097152}},\"implemented_options\":{implemented},\"remaining_options\":{remaining},{details},\"svg\":\"{}\"}}",escaped(&svg)))
}
fn histogram(groups: &[Group], a: &[f64], w: f64, h: f64) -> Result<String> {
    if a.len() < 8 {
        return Err(
            "histogram requires version,width,height,response,orientation,edge_count,edges".into(),
        );
    }
    let response = integer(a[3], 0., 2.)? as usize;
    let orientation = integer(a[4], 0., 1.)? as usize;
    let count = integer(a[5], 2., 513.)? as usize;
    if a.len() != 6 + count {
        return Err("histogram edge_count does not match arguments".into());
    }
    let edges = &a[6..];
    for pair in edges.windows(2) {
        if !pair[0].is_finite()
            || !pair[1].is_finite()
            || pair[0] >= pair[1]
            || checked(pair[1] - pair[0])? <= 0.
        {
            return Err("histogram edges must be finite, increasing and representable".into());
        }
    }
    let value_domain = domain(edges[0], edges[count - 1])?;
    let range = value_domain[1] - value_domain[0];
    let mut summaries = Vec::new();
    let mut max_response = 0f64;
    for g in groups {
        let mut bins = vec![0u64; count - 1];
        for r in &g.rows {
            if r.value < edges[0] || r.value > edges[count - 1] {
                return Err("selected value outside histogram edges".into());
            }
            let index = if r.value == edges[count - 1] {
                count - 2
            } else {
                edges.partition_point(|e| *e <= r.value) - 1
            };
            bins[index] += r.freq;
        }
        let mut probs = Vec::new();
        let mut percents = Vec::new();
        let mut densities = Vec::new();
        let mut responses = Vec::new();
        if g.n > 0 {
            for i in 0..bins.len() {
                let probability = bins[i] as f64 / g.n as f64;
                let percent = checked(100. * probability)?;
                let density = checked(probability / (edges[i + 1] - edges[i]))?;
                if probability > 0. && density == 0. {
                    return Err("positive density is not representable".into());
                }
                let y = [bins[i] as f64, percent, density][response];
                max_response = max_response.max(y);
                probs.push(probability);
                percents.push(percent);
                densities.push(density);
                responses.push(y);
            }
        }
        summaries.push((bins, probs, percents, densities, responses));
    }
    if max_response <= 0. {
        return Err("histogram response scale must be positive".into());
    }
    let active = groups.iter().filter(|g| g.n > 0).count();
    let columns = if active <= 2 {
        active
    } else {
        (active as f64).sqrt().ceil() as usize
    };
    let rows = active.div_ceil(columns);
    let cell_w = (w - 32.) / columns as f64;
    let cell_h = (h - 64.) / rows as f64;
    let vertical = orientation == 1;
    let value_labels = axis_label_width(value_domain)?;
    let response_labels = axis_label_width([0., max_response])?;
    let left_margin = (if vertical {
        value_labels + 16.
    } else {
        response_labels + 16.
    })
    .max(70.);
    let right_margin = (if vertical {
        response_labels / 2. + 8.
    } else {
        value_labels / 2. + 8.
    })
    .max(28.);
    let panels: Vec<Panel> = (0..active)
        .map(|i| Panel {
            left: 16. + (i % columns) as f64 * cell_w + left_margin,
            top: 42. + (i / columns) as f64 * cell_h + 36.,
            width: cell_w - left_margin - right_margin,
            height: cell_h - 100.,
        })
        .collect();
    if panels.iter().any(|p| p.width < 100. || p.height < 80.) {
        return Err("histogram facets need plot area at least 100 by 80 pixels".into());
    }
    let mut svg=start(w,h,"Histogram","Explicit bin edges; left-closed right-open bins, final right edge included. Counts and probabilities include hidden observations.");
    text(&mut svg, w / 2., 22., "middle", "Histogram");
    let mut output = Vec::new();
    let mut pi = 0;
    for (gi, g) in groups.iter().enumerate() {
        let prefix = group_prefix(g);
        if g.n == 0 {
            output.push(format!("{prefix},\"counts\":null,\"probabilities\":null,\"percents\":null,\"densities\":null,\"responses\":null,\"bins\":null,\"panel\":null}}"));
            continue;
        }
        let p = panels[pi];
        if format!("Group {} · N={}", g.id, g.n).len() as f64 * 7. > cell_w - 20. {
            return Err("histogram facet title cannot fit readably".into());
        }
        pi += 1;
        let (bins, probs, percents, densities, responses) = &summaries[gi];
        let vertical = orientation == 1;
        let (vs, vo) = transform(p, value_domain, vertical)?;
        let value_origin = value_domain[0];
        let svg_origin = if vertical { p.top + p.height } else { p.left };
        let rs = checked(if vertical {
            p.width / max_response
        } else {
            -p.height / max_response
        })?;
        let ro = if vertical { p.left } else { p.top + p.height };
        if rs == 0. {
            return Err("response transform is not representable".into());
        }
        let value_len = if vertical { p.height } else { p.width };
        if edges
            .windows(2)
            .any(|e| ((e[1] - e[0]) / range) * value_len < 2.)
        {
            return Err("histogram bins must occupy at least 2 pixels on value axis".into());
        }
        write!(svg,"<defs><clipPath id=\"panel-clip-{gi}\"><rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\"/></clipPath></defs>",p.left,p.top,p.width,p.height).unwrap();
        text(
            &mut svg,
            p.left + p.width / 2.,
            p.top - 27.,
            "middle",
            &format!("Group {} · N={}", g.id, g.n),
        );
        axes(&mut svg, p, value_domain, vertical, "Value")?;
        axes(
            &mut svg,
            p,
            [0., max_response],
            !vertical,
            ["Count", "Percent", "Density"][response],
        )?;
        let mut bj = Vec::new();
        write!(svg, "<g clip-path=\"url(#panel-clip-{gi})\">").unwrap();
        for i in 0..bins.len() {
            let v1 = map(p, value_domain, edges[i], vertical)?;
            let v2 = map(p, value_domain, edges[i + 1], vertical)?;
            let length =
                checked(responses[i] / max_response * if vertical { p.width } else { p.height })?;
            let (x, y, bw, bh) = if vertical {
                (p.left, v2, length, v1 - v2)
            } else {
                (v1, p.top + p.height - length, v2 - v1, length)
            };
            write!(svg,"<rect id=\"hist-{gi}-{i}\" class=\"bar\" data-role=\"histogram-bin\" data-group-id=\"{}\" data-bin-index=\"{i}\" data-count=\"{}\" x=\"{x}\" y=\"{y}\" width=\"{bw}\" height=\"{bh}\"><title>Group {}, bin [{},{}{}: count {}, probability {}, density {}</title></rect>",g.id,bins[i],g.id,edges[i],edges[i+1],if i+1==bins.len() {']'} else {')'},bins[i],probs[i],densities[i]).unwrap();
            bj.push(format!("{{\"left\":{},\"right\":{},\"count\":{},\"probability\":{},\"percent\":{},\"density\":{},\"response\":{}}}",edges[i],edges[i+1],bins[i],probs[i],percents[i],densities[i],responses[i]));
        }
        svg.push_str("</g>");
        // Machine-readable, accessible data is also available in each bar title and JSON.
        output.push(format!("{prefix},\"counts\":{},\"probabilities\":{},\"percents\":{},\"densities\":{},\"responses\":{},\"bins\":[{}],\"panel\":{},\"value_domain\":{},\"response_domain\":{},\"affine\":{{\"value_scale\":{vs},\"value_offset\":{vo},\"value_origin\":{value_origin},\"svg_origin\":{svg_origin},\"response_scale\":{rs},\"response_offset\":{ro}}}}}",arr(&bins.iter().map(|c|*c as f64).collect::<Vec<_>>()),arr(probs),arr(percents),arr(densities),arr(responses),bj.join(","),p.json(),arr(&value_domain),arr(&[0.,max_response])));
    }
    svg.push_str("</svg>");
    final_json("report_histogram_svg",groups,[w,h],"explicit_edges_native_left_closed_right_open_final_closed",&format!("\"orientation\":{orientation},\"response_scale\":{response},\"edges\":{},\"groups\":[{}]",arr(edges),output.join(",")),svg,("[\"explicit edges\",\"count percent density\",\"integer frequency\",\"group facets\",\"both orientations\"]","[\"JMP automatic bins and exact tie membership\",\"analytic Weight and fractional Freq\",\"KDE polygon shadowgram fill and overlays\",\"full styling and annotations\"]"))
}
struct BoxSummary {
    five: [f64; 5],
    iqr: f64,
    fences: [f64; 2],
    whiskers: [f64; 2],
    ticks: [f64; 4],
    outliers: Vec<Row>,
}
fn box_plot(groups: &[Group], a: &[f64], w: f64, h: f64) -> Result<String> {
    if a.len() != 7 {
        return Err("box plot requires seven parameters".into());
    }
    let kind = integer(a[3], 0., 1.)? as usize;
    let orientation = integer(a[4], 0., 1.)? as usize;
    let show = integer(a[5], 0., 1.)? == 1.;
    let labels = integer(a[6], 0., 1.)? == 1.;
    let vertical = orientation == 0;
    let mut summaries = Vec::new();
    let mut lo = f64::INFINITY;
    let mut hi = f64::NEG_INFINITY;
    let mut marks = 0;
    for g in groups {
        if g.n == 0 {
            summaries.push(None);
            continue;
        }
        let five = [
            g.rows[0].value,
            quantile(g, 0.25)?,
            quantile(g, 0.5)?,
            quantile(g, 0.75)?,
            g.rows.last().unwrap().value,
        ];
        let iqr = checked(five[3] - five[1])?;
        let extent = checked(1.5 * iqr)?;
        let fences = [checked(five[1] - extent)?, checked(five[3] + extent)?];
        let inside: Vec<&Row> = g
            .rows
            .iter()
            .filter(|r| r.value >= fences[0] && r.value <= fences[1])
            .collect();
        let whiskers = if kind == 1 {
            [five[0], five[4]]
        } else {
            [
                inside.first().ok_or("empty Tukey interior")?.value,
                inside.last().unwrap().value,
            ]
        };
        let outliers: Vec<_> = g
            .rows
            .iter()
            .filter(|r| r.value < fences[0] || r.value > fences[1])
            .cloned()
            .collect();
        if show && kind == 0 {
            marks += outliers.iter().filter(|r| !r.hidden).count();
        }
        let ticks = [
            quantile(g, 0.025)?,
            quantile(g, 0.10)?,
            quantile(g, 0.90)?,
            quantile(g, 0.975)?,
        ];
        lo = lo.min(five[0]);
        hi = hi.max(five[4]);
        summaries.push(Some(BoxSummary {
            five,
            iqr,
            fences,
            whiskers,
            ticks,
            outliers,
        }));
    }
    if marks > MAX_MARKS {
        return Err("report exceeds 20000 rendered individual markers".into());
    }
    let d = domain(lo, hi)?;
    let active = groups.iter().filter(|g| g.n > 0).count();
    let table_h = if labels {
        32. + active as f64 * 18.
    } else {
        0.
    };
    if labels && w < 650. {
        return Err("five-number table requires at least 650 pixels width".into());
    }
    let value_labels = axis_label_width(d)?;
    let group_labels = groups
        .iter()
        .filter(|g| g.n > 0)
        .map(|g| format!("Group {}", g.id).len())
        .max()
        .unwrap() as f64
        * 7.;
    let left_margin = if vertical {
        (value_labels + 16.).max(70.)
    } else {
        (group_labels + 16.).max(value_labels / 2. + 8.)
    };
    let right_margin = if vertical {
        40.
    } else {
        (value_labels / 2. + 8.).max(28.)
    };
    let p = Panel {
        left: left_margin,
        top: 76.,
        width: w - left_margin - right_margin,
        height: h - 146. - table_h,
    };
    if p.width < 100. || p.height < 80. {
        return Err("box plot needs plot area at least 100 by 80 pixels".into());
    }
    let band = if vertical {
        p.width / active as f64
    } else {
        p.height / active as f64
    };
    if band
        < if vertical {
            (group_labels + 8.).max(50.)
        } else {
            50.
        }
    {
        return Err(
            "box group labels and marks cannot fit readably (minimum 50 pixel band)".into(),
        );
    }
    let half = (band * 0.25).min(28.);
    let (vs, vo) = transform(p, d, vertical)?;
    let value_origin = d[0];
    let svg_origin = if vertical { p.top + p.height } else { p.left };
    let mut svg=start(w,h,"Box plot","Type 6 quantiles with integer-frequency replication. Hollow boxes show quartiles, thick lines medians, capped stems observed-data whiskers, circles visible Tukey outlier source rows. Hidden rows remain in summaries.");
    text(&mut svg, w / 2., 22., "middle", "Box plot");
    axes(&mut svg, p, d, vertical, "Value")?;
    write!(svg,"<defs><clipPath id=\"box-clip\"><rect x=\"{}\" y=\"{}\" width=\"{}\" height=\"{}\"/></clipPath></defs>",p.left-4.,p.top-4.,p.width+8.,p.height+8.).unwrap();
    if labels {
        let y = p.top + p.height + 60.;
        text(&mut svg, 16., y, "start", "Group");
        for (i, s) in ["Minimum", "Q1", "Median", "Q3", "Maximum"]
            .iter()
            .enumerate()
        {
            text(&mut svg, 150. + i as f64 * (w - 190.) / 5., y, "start", s);
        }
    }
    let mut output = Vec::new();
    let mut index = 0;
    for (gi, g) in groups.iter().enumerate() {
        let prefix = group_prefix(g);
        let Some(s) = &summaries[gi] else {
            output.push(format!("{prefix},\"five_number\":null,\"iqr\":null,\"fences\":null,\"whiskers\":null,\"quantile_ticks\":null,\"outliers\":[],\"panel\":null}}"));
            continue;
        };
        let center = if vertical {
            p.left + (index as f64 + 0.5) * band
        } else {
            p.top + (index as f64 + 0.5) * band
        };
        let q1 = map(p, d, s.five[1], vertical)?;
        let median = map(p, d, s.five[2], vertical)?;
        let q3 = map(p, d, s.five[3], vertical)?;
        let wl = map(p, d, s.whiskers[0], vertical)?;
        let wu = map(p, d, s.whiskers[1], vertical)?;
        if vertical {
            text(
                &mut svg,
                center,
                p.top + p.height + 18.,
                "middle",
                &format!("Group {}", g.id),
            );
            text(
                &mut svg,
                center,
                p.top - 28.,
                "middle",
                &format!("N={}", g.n),
            );
        } else {
            text(
                &mut svg,
                p.left - 10.,
                center + 4.,
                "end",
                &format!("Group {}", g.id),
            );
        }
        write!(
            svg,
            "<g clip-path=\"url(#box-clip)\" data-group-id=\"{}\">",
            g.id
        )
        .unwrap();
        write!(svg,"<title>Group {}; N={}; Minimum {}; Q1 {}; median {}; Q3 {}; maximum {}; whiskers {} to {}</title>",g.id,g.n,s.five[0],s.five[1],s.five[2],s.five[3],s.five[4],s.whiskers[0],s.whiskers[1]).unwrap();
        let common = format!("data-group-id=\"{}\"", g.id);
        if vertical {
            line(
                &mut svg,
                center,
                wl,
                center,
                wu,
                &format!("id=\"stem-{gi}\" data-role=\"whisker-stem\" {common}"),
            );
            for (end, y) in [("lower", wl), ("upper", wu)] {
                line(&mut svg,center-half*0.6,y,center+half*0.6,y,&format!("id=\"cap-{gi}-{end}\" data-role=\"whisker-cap\" data-end=\"{end}\" {common}"));
            }
            write!(svg,"<rect id=\"box-{gi}\" class=\"box\" data-role=\"box\" {common} x=\"{}\" y=\"{q3}\" width=\"{}\" height=\"{}\"/>",center-half,half*2.,q1-q3).unwrap();
            line(
                &mut svg,
                center - half,
                median,
                center + half,
                median,
                &format!("id=\"median-{gi}\" class=\"median\" data-role=\"median\" {common}"),
            );
        } else {
            line(
                &mut svg,
                wl,
                center,
                wu,
                center,
                &format!("id=\"stem-{gi}\" data-role=\"whisker-stem\" {common}"),
            );
            for (end, x) in [("lower", wl), ("upper", wu)] {
                line(&mut svg,x,center-half*0.6,x,center+half*0.6,&format!("id=\"cap-{gi}-{end}\" data-role=\"whisker-cap\" data-end=\"{end}\" {common}"));
            }
            write!(svg,"<rect id=\"box-{gi}\" class=\"box\" data-role=\"box\" {common} x=\"{q1}\" y=\"{}\" width=\"{}\" height=\"{}\"/>",center-half,q3-q1,half*2.).unwrap();
            line(
                &mut svg,
                median,
                center - half,
                median,
                center + half,
                &format!("id=\"median-{gi}\" class=\"median\" data-role=\"median\" {common}"),
            );
        }
        if kind == 1 {
            for (i, value) in s.ticks.iter().enumerate() {
                let v = map(p, d, *value, vertical)?;
                let attrs=format!("id=\"tick-{gi}-{i}\" data-role=\"quantile-tick\" data-probability=\"{}\" {common}",[0.025,0.1,0.9,0.975][i]);
                if vertical {
                    line(&mut svg, center - half, v, center + half, v, &attrs)
                } else {
                    line(&mut svg, v, center - half, v, center + half, &attrs)
                }
            }
        }
        if kind == 0 && show {
            for r in s.outliers.iter().filter(|r| !r.hidden) {
                let v = map(p, d, r.value, vertical)?;
                let (x, y) = if vertical { (center, v) } else { (v, center) };
                write!(svg,"<circle id=\"outlier-{}\" data-role=\"outlier\" {common} data-row-id=\"{}\" data-frequency=\"{}\" data-value=\"{}\" cx=\"{x}\" cy=\"{y}\" r=\"3\" fill=\"white\" stroke=\"#17212b\" stroke-width=\"1.5\"><title>Row {}, value {}, frequency {}</title></circle>",r.id,r.id,r.freq,r.value,r.id,r.value,r.freq).unwrap();
            }
        }
        svg.push_str("</g>");
        if labels {
            let y = p.top + p.height + 78. + index as f64 * 18.;
            text(&mut svg, 16., y, "start", &g.id.to_string());
            let mut cells = s.five.map(label);
            if (0..5).any(|i| (i + 1..5).any(|j| s.five[i] != s.five[j] && cells[i] == cells[j])) {
                cells = s.five.map(precise_label);
            }
            let column_width = (w - 190.) / 5.;
            if cells
                .iter()
                .any(|s| s.len() as f64 * 7. + 8. > column_width)
            {
                return Err(
                    "five-number table labels cannot fit readably; increase canvas width".into(),
                );
            }
            for (i, content) in cells.iter().enumerate() {
                let x = 150. + i as f64 * column_width;
                let value = s.five[i];
                let statistic = ["minimum", "q1", "median", "q3", "maximum"][i];
                write!(svg,"<text data-role=\"five-number-value\" data-group-id=\"{}\" data-statistic=\"{statistic}\" data-value=\"{value}\" x=\"{x}\" y=\"{y}\" text-anchor=\"start\">{content}</text>",g.id).unwrap();
            }
        }
        let out = s
            .outliers
            .iter()
            .map(|r| {
                format!(
                    "{{\"row_id\":{},\"value\":{},\"frequency\":{},\"hidden\":{}}}",
                    r.id, r.value, r.freq, r.hidden
                )
            })
            .collect::<Vec<_>>()
            .join(",");
        output.push(format!("{prefix},\"five_number\":{},\"iqr\":{},\"fences\":{},\"whiskers\":{},\"quantile_ticks\":{},\"outliers\":[{out}],\"panel\":{},\"value_domain\":{},\"affine\":{{\"value_scale\":{vs},\"value_offset\":{vo},\"value_origin\":{value_origin},\"svg_origin\":{svg_origin}}},\"group_coordinate\":{center},\"box_half_width\":{half}}}",arr(&s.five),s.iqr,arr(&s.fences),arr(&s.whiskers),arr(&s.ticks),p.json(),arr(&d)));
        index += 1;
    }
    svg.push_str("</svg>");
    final_json("report_box_plot_svg",groups,[w,h],"Hyndman_Fan_type_6_integer_replication_closed_Tukey_fences",&format!("\"outlier_classification\":\"closed_Tukey_fences_independent_of_box_type\",\"outlier_visibility\":\"{}\",\"box_type\":{kind},\"orientation\":{orientation},\"show_outliers\":{show},\"show_five_number_labels\":{labels},\"groups\":[{}]",if kind==1 {"not_applicable_quantile"} else if show {"shown"} else {"suppressed"},output.join(",")),svg,("[\"type6 quartiles\",\"closed Tukey fences and actual-data whiskers\",\"quantile ticks and full-range whiskers\",\"integer frequency\",\"hidden outlier provenance\",\"five-number table\",\"both orientations\"]","[\"analytic Weight and fractional Freq\",\"mean confidence diamond notches shortest half\",\"styles placement jitter\",\"exact JMP hidden and fence-boundary reference fixtures\"]"))
}
/// Render numerical summaries and deterministic self-contained SVG without dependencies.
pub fn run(operation: &str, data: &[Vec<f64>], args: &[f64]) -> Result<String> {
    if args.len() < 3 || args.iter().any(|x| !x.is_finite()) || args[0] != 1. {
        return Err("report requires finite version-1 parameters".into());
    }
    let width = integer(args[1], 320., 4096.)?;
    let height = integer(args[2], 240., 4096.)?;
    let groups = parse(data)?;
    match operation {
        "report_histogram_svg" => histogram(&groups, args, width, height),
        "report_box_plot_svg" => box_plot(&groups, args, width, height),
        _ => Err("unknown static report operation".into()),
    }
}
