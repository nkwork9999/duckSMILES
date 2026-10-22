#!/usr/bin/env python3
"""Independent CPU report oracle, public ABI/SQL checks, and real SVG rendering.

Only Python's standard library is required. Small frequency samples are expanded
and passed to statistics.quantiles(method='exclusive'), with JMP endpoint clamps;
the production cumulative-frequency implementation is never imported. Literal
hand fixtures separately distinguish type6 from type7 and data whiskers from
fences. Actual generated SVG is parsed and rasterized, not searched as source.

New report operations are a bounded subset, not full JMP graphic equivalence.
"""
import argparse
import copy
import ctypes
import datetime
from decimal import Decimal, localcontext
import hashlib
import json
import math
from pathlib import Path
import platform
import random
import re
import shutil
import statistics
import struct
import subprocess
import sys
import traceback
import xml.etree.ElementTree as ET


ROOT = Path(__file__).resolve().parents[1]
SVG_NS = "http://www.w3.org/2000/svg"
CHECKS = 0
PRIMARY = [
    "https://www.jmp.com/support/help/en/19.2/jmp/statistical-details-for-quantiles.shtml",
    "https://www.jmp.com/support/help/en/19.2/jmp/outlier-box-plot.shtml",
    "https://www.jmp.com/support/help/en/19.2/jmp/quantile-box-plot.shtml",
    "https://www.jmp.com/support/help/en/19.2/jmp/histogram.shtml",
    "https://www.jmp.com/support/help/en/19.2/jmp/launch-the-distribution-platform.shtml",
]


def require(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise AssertionError(message)


def near(actual, wanted, label, *, absolute=3e-12, relative=3e-12):
    require(isinstance(actual, (int, float)) and not isinstance(actual, bool),
            (label, "not numeric", actual))
    require(math.isfinite(actual), (label, "nonfinite", actual))
    tolerance = max(absolute, relative * abs(wanted), 4 * math.ulp(float(wanted)))
    require(abs(actual - wanted) <= tolerance, (label, actual, wanted, tolerance))


def same(actual, wanted, label):
    if isinstance(wanted, (list, tuple)):
        require(isinstance(actual, (list, tuple)) and len(actual) == len(wanted),
                (label, actual, wanted))
        for i, (a, b) in enumerate(zip(actual, wanted)):
            same(a, b, f"{label}[{i}]")
    elif isinstance(wanted, (int, float)) and not isinstance(wanted, bool):
        near(actual, wanted, label)
    else:
        require(actual == wanted, (label, actual, wanted))


def strict_json(raw):
    def reject(token):
        raise AssertionError(f"non-JSON numeric constant {token}")
    def unique(pairs):
        out = {}
        for key, value in pairs:
            require(key not in out, ("duplicate JSON key", key))
            out[key] = value
        return out
    result = json.loads(raw, parse_constant=reject, object_pairs_hook=unique)
    def finite(value):
        if isinstance(value, (int, float)) and not isinstance(value, bool):
            require(math.isfinite(value), ("nonfinite JSON number", value))
        elif isinstance(value, list):
            for v in value:
                finite(v)
        elif isinstance(value, dict):
            for v in value.values():
                finite(v)
    finite(result)
    return result


def source_rows(values, frequency=None, group=0, first_id=0):
    frequency = frequency if frequency is not None else [1] * len(values)
    return [[v, group, f, 0, 0, 1, first_id + i]
            for i, (v, f) in enumerate(zip(values, frequency))]


def selected_rows(rows, group):
    return [r for r in rows if r[1] == group and r[5] == 1 and r[3] == 0 and r[2] > 0]


def sample_quantiles(rows):
    # This reference physically expands only tiny independent test fixtures.
    # statistics.quantiles uses exclusive quantiles but extrapolates at its
    # endpoints, so clamp those to the sample range per the JMP formula.
    sample = sorted(v for r in rows for v in [r[0]] * int(r[2]))
    require(0 < len(sample) <= 4096, ("small independent sample bound", len(sample)))
    if len(sample) == 1:
        quantiles = [sample[0]] * 39
    else:
        # Center the independent library call before its weighted arithmetic:
        # statistics.quantiles on raw 1e12-scale values can move Q1 by one ULP.
        # The exact Decimal audit below documents that reference-only pitfall.
        origin = sample[0]
        quantiles = statistics.quantiles([v - origin for v in sample], n=40, method="exclusive")
        quantiles = [min(sample[-1], max(sample[0], origin + q)) for q in quantiles]
    return sample, {p: quantiles[index] for p, index in
                    [(0.025, 0), (0.10, 3), (0.25, 9), (0.5, 19),
                     (0.75, 29), (0.90, 35), (0.975, 38)]}


def box_oracle(rows):
    sample, quantiles = sample_quantiles(rows)
    q1, median, q3 = [quantiles[p] for p in [0.25, 0.5, 0.75]]
    iqr = q3 - q1
    fences = [q1 - 1.5 * iqr, q3 + 1.5 * iqr]
    inliers = [v for v in sample if fences[0] <= v <= fences[1]]
    outliers = [r for r in rows if r[0] < fences[0] or r[0] > fences[1]]
    return {"N": len(sample), "five_number": [sample[0], q1, median, q3, sample[-1]],
            "quartiles": [q1, median, q3], "IQR": iqr, "fences": fences,
            "whiskers": [min(inliers), max(inliers)], "quantiles": quantiles,
            "outliers": sorted(outliers, key=lambda r: r[6])}


def histogram_oracle(rows, edges):
    # Independently count each interval, not the production edge-search helper.
    n = sum(r[2] for r in rows)
    counts = [sum(r[2] for r in rows if lo <= r[0] < hi or
                  (i == len(edges) - 2 and r[0] == hi))
              for i, (lo, hi) in enumerate(zip(edges, edges[1:]))]
    require(sum(counts) == n, ("fixture outside histogram domain", rows, edges))
    return {"N": n, "counts": counts,
            "probabilities": [c / n for c in counts],
            "percents": [100 * c / n for c in counts],
            "densities": [c / n / (hi - lo) for c, lo, hi in
                          zip(counts, edges, edges[1:])]}


HISTOGRAM_LITERALS = [
    ("negative_edge_ties", [-2, -1, -.5, 0, .5, 1, 2], None,
     [-2, -1, 0, 1, 2], [1, 2, 2, 2]),
    ("unequal_bin_density", [0, .5, 1, 2, 3], None,
     [0, 1, 3], [2, 3]),
    ("replication_frequency", [0, .5, 1, 2], [2, 1, 3, 1],
     [0, 1, 2], [3, 4]),
    ("singleton_histogram", [42], None, [41, 43], [1]),
]
BOX_LITERALS = [
    ("type6_not_type7", [1, 2, 3, 4], None,
     [1.25, 2.5, 3.75], [-2.5, 7.5], [1, 4], []),
    ("actual_data_whiskers", [0, 1, 2, 3, 4, 5, 100], None,
     [1, 3, 5], [-5, 11], [0, 5], [100]),
    ("fence_equalities", [-5, 1, 2, 3, 4, 5, 11], None,
     [1, 3, 5], [-5, 11], [-5, 11], []),
    ("zero_iqr_outlier", [0, 0, 0, 0, 0, 0, 10], None,
     [0, 0, 0], [0, 0], [0, 0], [10]),
    ("singleton_box", [42], None,
     [42, 42, 42], [42, 42], [42, 42], []),
    ("two_observations", [-1, 1], None,
     [-1, 0, 1], [-4, 4], [-1, 1], []),
    ("frequency_box", [1, 2, 10], [1, 2, 1],
     [1.25, 2, 8], [-8.875, 18.125], [1, 10], []),
]


def literal_reference_checks():
    for name, values, frequency, edges, counts in HISTOGRAM_LITERALS:
        oracle = histogram_oracle(source_rows(values, frequency), edges)
        same(oracle["counts"], counts, name + ".counts")
    for name, values, frequency, quartiles, fences, whiskers, outliers in BOX_LITERALS:
        oracle = box_oracle(source_rows(values, frequency))
        for key, expected in [("quartiles", quartiles), ("fences", fences),
                              ("whiskers", whiskers)]:
            same(oracle[key], expected, name + "." + key)
        same([r[0] for r in oracle["outliers"]], outliers, name + ".outliers")


def large_offset_reference_audit(all_cases):
    case = next(c for c in all_cases if c["name"] == "large_offset_small_range_labels")
    values = sorted(r[0] for r in case["data"])
    exact = []
    with localcontext() as context:
        context.prec = 90
        for probability in [Decimal(".25"), Decimal(".5"), Decimal(".75")]:
            rank = (len(values) + 1) * probability
            index = int(rank)
            fraction = rank - index
            value = ((1 - fraction) * Decimal.from_float(values[index - 1])
                     + fraction * Decimal.from_float(values[index]))
            exact.append({"probability": str(probability), "exact_real_quantile": str(value),
                          "rounded_f64": float(value), "rounded_f64_hex": float(value).hex()})
    wanted = [q["rounded_f64"] for q in exact]
    oracle = box_oracle(case["data"])
    # Here exact equality is intentional: these are rounded binary64 values
    # from a 90-digit reference, not generic tolerance-based statistical checks.
    require(oracle["quartiles"] == wanted, ("centered reference versus Decimal90", oracle, exact))
    require(oracle["IQR"] == wanted[2] - wanted[0], "IQR must use returned rounded quartiles")
    return {"method": "Decimal90, exact Decimal.from_float inputs, official type6 ranks; round each output quantile once to f64",
            "input_hex": [v.hex() for v in values], "quantiles": exact,
            "rounded_quantile_IQR": wanted[2] - wanted[0],
            "uncentered_statistics_quartiles": statistics.quantiles(values, n=4, method="exclusive"),
            "lesson": "The uncentered Python reference Q1 is one ULP wrong here; production matched the Decimal90-rounded quartiles. No tolerance was loosened to accept the difference."}


class Abi:
    def __init__(self, path):
        self.path = Path(path).resolve()
        self.library = ctypes.CDLL(str(self.path))
        double_p = ctypes.POINTER(ctypes.c_double)
        self.run = self.library.ds_analytics_run
        self.run.argtypes = [ctypes.c_char_p, ctypes.c_size_t, double_p,
                             ctypes.c_size_t, ctypes.c_size_t, double_p,
                             ctypes.c_size_t, ctypes.POINTER(ctypes.c_size_t)]
        self.run.restype = ctypes.c_void_p
        self.free = self.library.ds_analytics_free
        self.free.argtypes = [ctypes.c_void_p, ctypes.c_size_t]
        self.free.restype = None
        self.calls = 0

    def nullable(self, operation, data, parameters):
        columns = len(data[0]) if data else 0
        require(all(len(row) == columns for row in data), "C ABI fixture must be rectangular")
        flat = [float(v) for row in data for v in row]
        numbers = (ctypes.c_double * len(flat))(*flat)
        args = (ctypes.c_double * len(parameters))(*parameters)
        name = operation.encode()
        length = ctypes.c_size_t(12345)
        pointer = self.run(name, len(name), numbers, len(data), columns,
                           args, len(parameters), ctypes.byref(length))
        self.calls += 1
        if not pointer:
            require(length.value == 0, ("NULL ABI buffer with nonzero length", length.value))
            return None
        try:
            return strict_json(ctypes.string_at(pointer, length.value))
        finally:
            self.free(pointer, length.value)


def sql_value(value):
    if value is None:
        return "NULL"
    if isinstance(value, (list, tuple)):
        return "[" + ",".join(sql_value(v) for v in value) + "]"
    if isinstance(value, bool):
        return "1.0" if value else "0.0"
    if math.isnan(float(value)):
        return "'NaN'::DOUBLE"
    if math.isinf(float(value)):
        return ("'Infinity'" if value > 0 else "'-Infinity'") + "::DOUBLE"
    return repr(float(value))


def expression(operation, data, parameters):
    return ("chem_analyze('" + operation.replace("'", "''") + "',"
            + sql_value(data) + "::DOUBLE[][]," + sql_value(parameters) + "::DOUBLE[])")


class Sql:
    def __init__(self, cli, extension):
        self.cli = str(Path(cli).resolve())
        self.extension = Path(extension).resolve()
        self.queries = 0
        self.calls = 0

    def query(self, sql):
        extension = str(self.extension).replace("'", "''")
        process = subprocess.run([self.cli, "-unsigned", "-batch", "-bail", "-json", ":memory:"],
                                 input=f"LOAD '{extension}';\n{sql}\n", text=True,
                                 capture_output=True, timeout=120)
        self.queries += 1
        require(process.returncode == 0, ("DuckDB SQL failure", process.returncode, process.stderr))
        return strict_json(process.stdout)

    def many(self, requests):
        sql = " UNION ALL ".join(f"SELECT {i} AS id,{expression(*r)} AS result"
                                  for i, r in enumerate(requests))
        rows = self.query(sql + " ORDER BY id;")
        require(len(rows) == len(requests), ("SQL result length", len(rows), len(requests)))
        output = []
        for i, row in enumerate(rows):
            require(row["id"] == i, ("SQL selected row identity", row, i))
            output.append(None if row["result"] is None else strict_json(row["result"]))
        self.calls += len(output)
        return output

    def nullable(self, operation, data, parameters):
        return self.many([(operation, data, parameters)])[0]


def input_sha(path):
    path = Path(path).resolve()
    return {"path": str(path), "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}


def parse_svg(svg, width, height):
    require(isinstance(svg, str) and 0 < len(svg.encode()) <= 2097152,
            ("SVG UTF-8 size", type(svg)))
    require("<!DOCTYPE" not in svg.upper() and "<!ENTITY" not in svg.upper(),
            "SVG must not contain a DTD/entity declaration")
    root = ET.fromstring(svg)
    require(root.tag == "{" + SVG_NS + "}svg", ("SVG namespace", root.tag))
    near(float(root.attrib["width"]), width, "SVG canvas width")
    near(float(root.attrib["height"]), height, "SVG canvas height")
    viewbox = [float(x) for x in root.attrib["viewBox"].replace(",", " ").split()]
    same(viewbox, [0, 0, width, height], "SVG viewBox")
    allowed = {"svg", "g", "title", "desc", "rect", "line", "circle", "text", "style",
               "defs", "clipPath", "polyline", "polygon"}
    numeric = {"x", "y", "x1", "y1", "x2", "y2", "cx", "cy", "r", "width",
               "height", "stroke-width", "font-size", "opacity", "fill-opacity"}
    ids = set()
    for node in root.iter():
        require(node.tag.startswith("{" + SVG_NS + "}"), ("foreign SVG namespace", node.tag))
        local = node.tag.rsplit("}", 1)[-1]
        require(local in allowed, ("unsafe/unexpected SVG element", local))
        if local == "style":
            css = node.text or ""
            require(not any(v in css.lower() for v in ["@import", "javascript:", "expression("]),
                    "External/executable inline SVG CSS")
            for target in re.findall(r"url\((.*?)\)", css, re.IGNORECASE):
                require(re.fullmatch(r"#[A-Za-z0-9_-]+", target) is not None,
                        ("External SVG CSS URL", target))
            fonts = re.findall(r"font-size\s*:\s*([0-9.]+)px", css, re.IGNORECASE)
            require(fonts and all(float(size) >= 12 for size in fonts),
                    ("Unreadably small/unverified inline CSS font size", fonts))
        for key, value in node.attrib.items():
            require(not key.lower().startswith("on") and "href" not in key.lower(),
                    ("external/event SVG attribute", key))
            require("javascript:" not in value.lower(), ("unsafe SVG value", value))
            if "url(" in value.lower():
                require(re.fullmatch(r"url\(#[A-Za-z0-9_-]+\)", value) is not None,
                        ("external SVG URL", key, value))
            if key in numeric:
                number = float(value)
                require(math.isfinite(number), ("nonfinite SVG geometry", key, value))
                if key in {"width", "height", "r", "stroke-width", "font-size"}:
                    require(number >= 0, ("negative SVG extent", key, number))
                if key == "font-size":
                    require(number >= 12, ("unreadably small SVG text", number))
            if key == "points":
                points = [float(x) for x in value.replace(",", " ").split()]
                require(len(points) % 2 == 0 and all(math.isfinite(x) for x in points),
                        ("malformed SVG points", value))
            if key == "id":
                require(value not in ids, ("duplicate SVG id", value))
                ids.add(value)
        for key in ["x", "x1", "x2", "cx"]:
            if key in node.attrib:
                require(-.001 <= float(node.attrib[key]) <= width + .001,
                        ("SVG x outside canvas", key, node.attrib[key]))
        for key in ["y", "y1", "y2", "cy"]:
            if key in node.attrib:
                require(-.001 <= float(node.attrib[key]) <= height + .001,
                        ("SVG y outside canvas", key, node.attrib[key]))
        if local == "rect":
            require(float(node.get("x", 0)) + float(node.get("width", 0)) <= width + .001,
                    ("SVG rectangle exceeds width", node.attrib))
            require(float(node.get("y", 0)) + float(node.get("height", 0)) <= height + .001,
                    ("SVG rectangle exceeds height", node.attrib))
    require(any(n.tag == "{" + SVG_NS + "}title" and (n.text or "").strip() for n in root),
            "SVG needs a nonempty title")
    return root


def roles(root, role, group=None):
    return [n for n in root.iter() if n.get("data-role") == role
            and (group is None or float(n.get("data-group-id", "nan")) == group)]


def coordinate(actual, expected, label):
    near(float(actual), expected, label, absolute=.0011, relative=1e-10)


def assert_rejected_svg(svg, width, height, mutation):
    root = ET.fromstring(svg)
    mutation(root)
    try:
        parse_svg(ET.tostring(root, encoding="unicode"), width, height)
    except (AssertionError, ValueError, ET.ParseError):
        return
    raise AssertionError("Independent SVG validator accepted a malicious/malformed mutation")


def rasterize(svg_path, png_path, magick, width, height, anchors, renderer):
    command = ([magick, str(svg_path), str(png_path)] if renderer == "magick" else
               ["/usr/bin/sips", "-s", "format", "png", str(svg_path), "--out", str(png_path)])
    process = subprocess.run(command,
                             capture_output=True, text=True, timeout=30)
    require(process.returncode == 0, ("actual SVG renderer failed", process.stderr))
    require(not process.stderr.strip(), ("actual SVG renderer warning", process.stderr))
    data = png_path.read_bytes()
    require(data[:8] == b"\x89PNG\r\n\x1a\n", "rendered artifact is not PNG")
    same(list(struct.unpack(">II", data[16:24])), [width, height], "PNG dimensions")
    # Independently decode the produced PNG through the installed rasterizer.
    # This checks actual raster content, not an SVG string or byte snapshot.
    decoded = subprocess.run([magick, str(png_path), "-depth", "8", "rgba:-"],
                             capture_output=True, timeout=30)
    require(decoded.returncode == 0, ("PNG decode failed", decoded.stderr))
    rgba = decoded.stdout
    require(len(rgba) == width * height * 4, ("PNG RGBA extent", len(rgba)))
    nonwhite = sum(rgba[i + 3] > 0 and min(rgba[i:i + 3]) < 240
                   for i in range(0, len(rgba), 4))
    require(nonwhite > 100, ("blank actual rendered report", nonwhite))
    border_ink = 0
    for y in range(height):
        for x in ([0, 1, width - 2, width - 1] if 2 <= y < height - 2 else range(width)):
            i = (y * width + x) * 4
            border_ink += rgba[i + 3] > 0 and min(rgba[i:i + 3]) < 240
    require(border_ink == 0, ("accepted report clips visible content at canvas border", border_ink))
    for label, x, y in anchors:
        dark = False
        for py in range(max(0, round(y) - 2), min(height, round(y) + 3)):
            for px in range(max(0, round(x) - 2), min(width, round(x) + 3)):
                i = (py * width + px) * 4
                dark |= rgba[i + 3] > 0 and min(rgba[i:i + 3]) < 245
        require(dark, ("semantic SVG mark absent from actual raster", label, x, y))
    return {"svg": input_sha(svg_path), "png": input_sha(png_path),
            "width": width, "height": height, "nonwhite_pixels": nonwhite,
            "raster_mark_anchors_checked": len(anchors), "svg_rasterizer": renderer,
            "png_decoder": "ImageMagick_RGBA"}


def state_counts(rows):
    selected = [r for r in rows if r[5] == 1 and r[3] == 0 and r[2] > 0]
    return {"source_rows": len(rows), "selected_rows": len(selected),
            "N": sum(r[2] for r in selected),
            "excluded_rows": sum(r[3] == 1 for r in rows),
            "missing_rows": sum(r[5] == 0 for r in rows),
            "zero_frequency_rows": sum(r[2] == 0 for r in rows),
            "hidden_rows": sum(r[4] == 1 for r in rows)}


def panel_check(panel, width, height, name):
    for key in ["left", "top", "width", "height"]:
        require(isinstance(panel[key], (int, float)) and math.isfinite(panel[key]),
                (name, "panel field", key, panel))
    require(panel["left"] >= 0 and panel["top"] >= 0
            and panel["width"] > 0 and panel["height"] > 0,
            (name, "panel has no positive drawable area", panel))
    require(panel["left"] + panel["width"] <= width
            and panel["top"] + panel["height"] <= height,
            (name, "panel outside canvas", panel))


def affine_value(group, value, axis):
    """Reconstruct coordinates from domain and panel, never returned affine."""
    lower, upper = group["value_domain"]
    panel = group["panel"]
    fraction = (value - lower) / (upper - lower)
    return (panel["left"] + panel["width"] * fraction if axis == "x" else
            panel["top"] + panel["height"] * (1 - fraction))


def check_line(line, expected, label):
    for key, value in zip(["x1", "y1", "x2", "y2"], expected):
        coordinate(line.attrib[key], value, label + "." + key)


def check_histogram(result, case, root):
    edges = case["parameters"][6:]
    scale = int(case["parameters"][3])
    orientation = int(case["parameters"][4])
    rows = case["data"]
    expected = {}
    for group in result["groups"]:
        selected = selected_rows(rows, group["group_id"])
        if selected:
            expected[group["group_id"]] = histogram_oracle(selected, edges)
    response_key = ["counts", "percents", "densities"][scale]
    global_max = max(v for oracle in expected.values() for v in oracle[response_key])
    anchors = []
    for group in result["groups"]:
        gid = group["group_id"]
        if gid not in expected:
            require(group["counts"] is None and group["panel"] is None,
                    (case["name"], "empty histogram geometry", group))
            require(not roles(root, "histogram-bin", gid), "Empty group rendered histogram bins")
            continue
        oracle = expected[gid]
        for actual_key, reference_key in [("counts", "counts"), ("probabilities", "probabilities"),
                                          ("percents", "percents"), ("densities", "densities"),
                                          ("responses", response_key)]:
            same(group[actual_key], oracle[reference_key], case["name"] + "." + actual_key)
        same(group["value_domain"], [edges[0], edges[-1]], "histogram edge domain")
        same(group["response_domain"], [0, global_max], "shared histogram response domain")
        p = group["panel"]
        bins = group["bins"]
        require(len(bins) == len(edges) - 1, ("histogram JSON bin count", bins))
        marks = roles(root, "histogram-bin", gid)
        require(len(marks) == len(bins), ("SVG histogram bin count", len(marks), len(bins)))
        by_index = {int(n.attrib["data-bin-index"]): n for n in marks}
        require(len(by_index) == len(bins) and sorted(by_index) == list(range(len(bins))),
                ("histogram bin provenance", by_index.keys()))
        for j, (bin_result, lo, hi, count, response) in enumerate(zip(
                bins, edges, edges[1:], oracle["counts"], oracle[response_key])):
            for field, value in [("left", lo), ("right", hi), ("count", count),
                                 ("probability", oracle["probabilities"][j]),
                                 ("percent", oracle["percents"][j]),
                                 ("density", oracle["densities"][j]), ("response", response)]:
                near(bin_result[field], value, f"{case['name']}.bins[{j}].{field}")
            mark = by_index[j]
            require(mark.tag == "{" + SVG_NS + "}rect", "Histogram semantic mark must be a rectangle")
            near(float(mark.attrib["data-count"]), count, "SVG bin count provenance")
            if orientation == 0:
                x = affine_value(group, lo, "x")
                w = affine_value(group, hi, "x") - x
                h = p["height"] * response / global_max
                y = p["top"] + p["height"] - h
            else:
                y = affine_value(group, hi, "y")
                h = affine_value(group, lo, "y") - y
                w = p["width"] * response / global_max
                x = p["left"]
            for key, value in [("x", x), ("y", y), ("width", w), ("height", h)]:
                coordinate(mark.attrib[key], value, f"histogram {gid}/{j}.{key}")
            if count > 0 and w > 5 and h > 5:
                anchors.append((f"bar-{gid}-{j}", x + w / 2, y + h / 2))
        value_axis = "x" if orientation == 0 else "y"
        expected_value_scale = (p["width"] if value_axis == "x" else -p["height"]) / (edges[-1] - edges[0])
        expected_value_offset = affine_value(group, edges[0], value_axis) - expected_value_scale * edges[0]
        expected_response_scale = (-p["height"] if orientation == 0 else p["width"]) / global_max
        expected_response_offset = p["top"] + p["height"] if orientation == 0 else p["left"]
        for key, wanted in [("value_scale", expected_value_scale), ("value_offset", expected_value_offset),
                            ("response_scale", expected_response_scale), ("response_offset", expected_response_offset)]:
            near(group["affine"][key], wanted, "independent histogram affine." + key)
        same(group["affine"]["value_origin"], edges[0], "stable histogram value origin")
        same(group["affine"]["svg_origin"], affine_value(group, edges[0], value_axis), "stable histogram SVG origin")
    require(len(roles(root, "histogram-bin")) == sum(len(v["counts"]) for v in expected.values()),
            "Unexpected histogram group marks")
    return anchors


def check_box(result, case, root):
    box_type, orientation, show_outliers = map(int, case["parameters"][3:6])
    rows = case["data"]
    all_selected = [r for r in rows if r[5] == 1 and r[3] == 0 and r[2] > 0]
    observed_min = min(r[0] for r in all_selected)
    observed_max = max(r[0] for r in all_selected)
    expected_visible_outliers = 0
    anchors = []
    for group in result["groups"]:
        gid = group["group_id"]
        selected = selected_rows(rows, gid)
        if not selected:
            require(group["five_number"] is None and group["panel"] is None,
                    ("Empty box group has geometry", group))
            require(not roles(root, "box", gid), "Empty group rendered a box")
            continue
        oracle = box_oracle(selected)
        for key, wanted in [("five_number", oracle["five_number"]), ("iqr", oracle["IQR"]),
                            ("fences", oracle["fences"]),
                            ("quantile_ticks", [oracle["quantiles"][q] for q in [.025, .1, .9, .975]])]:
            same(group[key], wanted, case["name"] + "." + key)
        whiskers = oracle["whiskers"] if box_type == 0 else [oracle["five_number"][0], oracle["five_number"][-1]]
        same(group["whiskers"], whiskers, case["name"] + ".whiskers")
        actual_outliers = sorted(group["outliers"], key=lambda o: o["row_id"])
        require(len(actual_outliers) == len(oracle["outliers"]), "Box outlier metadata count")
        for actual, wanted in zip(actual_outliers, oracle["outliers"]):
            for key, value in [("row_id", wanted[6]), ("value", wanted[0]), ("frequency", wanted[2])]:
                near(actual[key], value, "outlier provenance." + key)
            require(isinstance(actual["hidden"], bool) and actual["hidden"] == bool(wanted[4]),
                    ("outlier hidden provenance", actual))
        domain = group["value_domain"]
        require(len(domain) == 2 and domain[0] < domain[1], ("nonpositive box axis domain", domain))
        require(domain[0] <= observed_min <= observed_max <= domain[1],
                ("Box domain clips observed values", domain, observed_min, observed_max))
        if observed_min < observed_max:
            same(domain, [observed_min, observed_max], "shared observed box value domain")
        axis = "y" if orientation == 0 else "x"
        p = group["panel"]
        gc, half = group["group_coordinate"], group["box_half_width"]
        require(half > 0, ("nonpositive box width", half))
        require((p["left"] <= gc - half <= gc + half <= p["left"] + p["width"]) if axis == "y" else
                (p["top"] <= gc - half <= gc + half <= p["top"] + p["height"]),
                ("box category band outside panel", gc, half, p))
        q1, median, q3 = oracle["quartiles"]
        mapped = {k: affine_value(group, v, axis) for k, v in
                  [("q1", q1), ("q3", q3), ("median", median),
                   ("lower", whiskers[0]), ("upper", whiskers[1])]}
        boxes = roles(root, "box", gid)
        require(len(boxes) == 1 and boxes[0].tag == "{" + SVG_NS + "}rect", "One box rectangle per group")
        if orientation == 0:
            rect = [gc - half, mapped["q3"], 2 * half, mapped["q1"] - mapped["q3"]]
            med = [gc - half, mapped["median"], gc + half, mapped["median"]]
            stem = [gc, mapped["lower"], gc, mapped["upper"]]
        else:
            rect = [mapped["q1"], gc - half, mapped["q3"] - mapped["q1"], 2 * half]
            med = [mapped["median"], gc - half, mapped["median"], gc + half]
            stem = [mapped["lower"], gc, mapped["upper"], gc]
        for key, wanted in zip(["x", "y", "width", "height"], rect):
            coordinate(boxes[0].attrib[key], wanted, "SVG box." + key)
        medians = roles(root, "median", gid)
        require(len(medians) == 1, "One visible median per box")
        check_line(medians[0], med, "SVG median")
        anchors.append(("median-" + str(gid), (med[0] + med[2]) / 2, (med[1] + med[3]) / 2))
        stems = roles(root, "whisker-stem", gid)
        require(len(stems) == 1, "One full whisker stem per box")
        check_line(stems[0], stem, "SVG whisker stem")
        caps = roles(root, "whisker-cap", gid)
        require(len(caps) == 2, "Two whisker caps per box")
        cap_by_end = {n.attrib["data-end"]: n for n in caps}
        require(set(cap_by_end) == {"lower", "upper"}, "Whisker cap endpoint provenance")
        for end in ["lower", "upper"]:
            cap = ([gc - .6 * half, mapped[end], gc + .6 * half, mapped[end]] if orientation == 0 else
                   [mapped[end], gc - .6 * half, mapped[end], gc + .6 * half])
            check_line(cap_by_end[end], cap, "SVG " + end + " whisker cap")
        ticks = roles(root, "quantile-tick", gid)
        require(len(ticks) == (4 if box_type == 1 else 0), "Quantile tick style/count")
        tick_by_probability = {float(n.attrib["data-probability"]): n for n in ticks}
        if box_type == 1:
            require(set(tick_by_probability) == {.025, .1, .9, .975}, "Quantile tick provenance")
            for prob, tick in tick_by_probability.items():
                value = affine_value(group, oracle["quantiles"][prob], axis)
                if orientation == 0:
                    coordinate(tick.attrib["y1"], value, "quantile tick y1")
                    coordinate(tick.attrib["y2"], value, "quantile tick y2")
                else:
                    coordinate(tick.attrib["x1"], value, "quantile tick x1")
                    coordinate(tick.attrib["x2"], value, "quantile tick x2")
        visible = [r for r in oracle["outliers"] if box_type == 0 and show_outliers and not r[4]]
        markers = roles(root, "outlier", gid)
        require(len(markers) == len(visible), ("outlier marker visibility", gid, len(markers), len(visible)))
        expected_visible_outliers += len(visible)
        by_id = {int(n.attrib["data-row-id"]): n for n in markers}
        require(set(by_id) == {r[6] for r in visible}, "Outlier row-ID preservation")
        for row in visible:
            mark = by_id[row[6]]
            require(mark.tag == "{" + SVG_NS + "}circle", "Outliers are source-row circles")
            near(float(mark.attrib["data-frequency"]), row[2], "outlier frequency metadata")
            near(float(mark.attrib["data-value"]), row[0], "outlier value metadata")
            position = affine_value(group, row[0], axis)
            cx, cy = (gc, position) if orientation == 0 else (position, gc)
            coordinate(mark.attrib["cx"], cx, "outlier circle cx")
            coordinate(mark.attrib["cy"], cy, "outlier circle cy")
            anchors.append(("outlier-" + str(row[6]), cx - float(mark.attrib["r"]), cy))
        expected_scale = (p["width"] if axis == "x" else -p["height"]) / (domain[1] - domain[0])
        expected_offset = affine_value(group, domain[0], axis) - expected_scale * domain[0]
        near(group["affine"]["value_scale"], expected_scale, "independent box affine scale")
        near(group["affine"]["value_offset"], expected_offset, "independent box affine offset")
        same(group["affine"]["value_origin"], domain[0], "stable box value origin")
        same(group["affine"]["svg_origin"], affine_value(group, domain[0], axis), "stable box SVG origin")
    require(len(roles(root, "outlier")) == expected_visible_outliers, "Unexpected outlier group")
    return anchors


def check_large_offset_five_number_text(case, result, root):
    """Catch the observed visible-table collapse, even before semantic attrs exist."""
    wanted = box_oracle(selected_rows(case["data"], 0))["five_number"]
    names = ["Minimum", "Q1", "Median", "Q3", "Maximum"]
    statistics_names = ["minimum", "q1", "median", "q3", "maximum"]
    texts = [n for n in root.iter() if n.tag == "{" + SVG_NS + "}text"]
    labelled = roles(root, "five-number-value", 0)
    actual = []
    if labelled:
        by_statistic = {n.attrib["data-statistic"]: n for n in labelled}
        require(set(by_statistic) == set(statistics_names), "Five-number visible cell provenance")
        for statistic, value in zip(statistics_names, wanted):
            cell = by_statistic[statistic]
            same(float(cell.attrib["data-value"]), value, "Five-number exact numerical metadata")
            actual.append("".join(cell.itertext()).strip())
    else:
        # Existing frozen ABI has no per-cell semantic data attributes. Locate
        # its actual user-visible cells under their header columns. This makes
        # the RED prove the rounding defect, not merely missing new metadata.
        for name in names:
            headers = [n for n in texts if "".join(n.itertext()).strip() == name]
            require(len(headers) == 1, ("Visible five-number column heading", name))
            header = headers[0]
            cells = [n for n in texts
                     if abs(float(n.attrib["x"]) - float(header.attrib["x"])) < .001
                     and float(n.attrib["y"]) > float(header.attrib["y"])]
            require(cells, ("Missing visible cell under five-number header", name))
            cell = min(cells, key=lambda n: float(n.attrib["y"]))
            actual.append("".join(cell.itertext()).strip())
    require(len(set(actual)) == len(set(wanted)),
            ("Visible five-number table collapses distinct large-offset values", actual, wanted))
    for name, text, value in zip(names, actual, wanted):
        require(float(text.replace("−", "-")) == value,
                ("Visible five-number value must retain roundtrip precision for this narrow fixture",
                 name, text, value))


def verify_result(case, result):
    require(isinstance(result, dict), (case["name"], "public result is NULL/not object", result))
    require(result["version"] == 1 and result["schema"] == "ducksmiles-cpu-plot-v1", "Versioned report JSON schema")
    same(result["operation"], case["operation"], "public operation provenance")
    same(result["media_type"], "image/svg+xml", "report media type")
    width, height = map(int, case["parameters"][1:3])
    same([result["width"], result["height"]], [width, height], "canvas dimensions")
    totals = state_counts(case["data"])
    for key in ["source_rows", "selected_rows", "N"]:
        same(result[key], totals[key], "top-level denominator." + key)
    group_ids = sorted({r[1] for r in case["data"]})
    same([g["group_id"] for g in result["groups"]], group_ids, "ascending group order including empty groups")
    for group in result["groups"]:
        group_rows = [r for r in case["data"] if r[1] == group["group_id"]]
        counts = state_counts(group_rows)
        for key, wanted in counts.items():
            same(group[key], wanted, "source-state denominator." + key)
        same(group["status"], "ok" if counts["N"] else "empty", "group status")
        if counts["N"]:
            panel_check(group["panel"], width, height, case["name"])
    root = parse_svg(result["svg"], width, height)
    tick_labels = roles(root, "value-tick-label")
    require(len(tick_labels) >= 2, "Accepted report needs distinguishable numeric value-axis labels")
    text_values = {}
    for tick in tick_labels:
        text = "".join(tick.itertext()).strip()
        value = float(tick.attrib["data-value"])
        require(text, "Empty value-axis tick label")
        parsed = float(text.replace("−", "-"))
        domains = [g["value_domain"] for g in result["groups"] if g["status"] == "ok"]
        span = max(hi - lo for lo, hi in domains)
        near(parsed, value, "actual tick text represents its numerical value",
             absolute=max(1e-14, span * 1e-6), relative=0)
        text_values.setdefault(text, set()).add(value)
    require(all(len(values) == 1 for values in text_values.values()),
            ("different axis values collapse to identical displayed text", text_values))
    require(len(text_values) >= 2, "All value-axis tick text collapsed to one string")
    if case["operation"] == "report_histogram_svg":
        anchors = check_histogram(result, case, root)
    else:
        same(result["outlier_classification"], "closed_Tukey_fences_independent_of_box_type", "outlier classification scope")
        visibility = ("not_applicable_quantile" if case["parameters"][3] == 1 else
                      "shown" if case["parameters"][5] == 1 else "suppressed")
        same(result["outlier_visibility"], visibility, "outlier display scope")
        anchors = check_box(result, case, root)
        if case["name"] == "large_offset_small_range_labels":
            check_large_offset_five_number_text(case, result, root)
    return root, anchors


def cases():
    output = []
    def hist(name, rows, edges, scale=0, orientation=0, render=True, size=(800, 500)):
        output.append({"name": name, "operation": "report_histogram_svg", "data": rows,
                       "parameters": [1, *size, scale, orientation, len(edges), *edges],
                       "render": render})
    def box(name, rows, kind=0, orientation=0, show=1, render=True, size=(800, 500)):
        output.append({"name": name, "operation": "report_box_plot_svg", "data": rows,
                       "parameters": [1, *size, kind, orientation, show, 1], "render": render})
    for name, values, frequency, edges, _ in HISTOGRAM_LITERALS:
        hist(name, source_rows(values, frequency), edges)
    hist("unequal_density_response", source_rows([0, .5, 1, 2, 3]), [0, 1, 3], scale=2)
    hist("frequency_percent_response", source_rows([0, .5, 1, 2], [2, 1, 3, 1]), [0, 1, 2], scale=1)
    for name, values, frequency, *_ in BOX_LITERALS:
        box(name, source_rows(values, frequency))
    box("quantile_full_range_and_four_ticks", source_rows([0, 1, 2, 3, 4, 5, 100]), kind=1)
    box("interior_quantile_ticks", source_rows(list(range(41))), kind=1)
    box("suppressed_outlier_keeps_statistics", source_rows([0, 1, 2, 3, 4, 5, 100]), show=0)
    # Independent flags deliberately overlap: hidden/excluded/missing are not
    # a single priority-coded status. Empty group9 remains in the denominator.
    states = source_rows([0, 1, 2, 3, 4, 5, 100], first_id=10)
    states[-1][4] = 1
    states += [[999, 0, 1, 1, 0, 1, 1001], [888, 0, 1, 0, 0, 0, 1002],
               [777, 0, 0, 0, 0, 1, 1003], [666, 0, 2, 1, 1, 0, 1004]]
    states += source_rows([-3, -1, 0, 5], [1, 2, 1, 1], group=7, first_id=1100)
    states += [[12345, 9, 1, 1, 0, 1, 2000], [12345, 9, 1, 0, 0, 0, 2001]]
    hist("grouped_source_states", states, [-10, 0, 2, 10, 101], size=(1000, 800))
    hist("grouped_horizontal_density", states, [-10, 0, 2, 10, 101], scale=2,
         orientation=1, size=(1000, 800))
    for kind, orientation in [(0, 0), (0, 1), (1, 0), (1, 1)]:
        box(f"grouped_states_kind{kind}_orientation{orientation}", states,
            kind=kind, orientation=orientation, size=(1000, 800))
    adjacent = [-2, math.nextafter(-1, -math.inf), -1, math.nextafter(-1, math.inf),
                math.nextafter(0, -math.inf), 0, math.nextafter(0, math.inf),
                math.nextafter(1, -math.inf), 1, math.nextafter(1, math.inf), 2]
    hist("adjacent_floating_edge_membership", source_rows(adjacent), [-2, -1, 0, 1, 2], render=False)
    hist("narrow_distinct_axis_labels", source_rows([1, 1 + 1e-8, 1 + 2e-8]),
         [1, 1 + 1e-8, 1 + 2e-8], size=(1000, 600))
    box("large_offset_small_range_labels", source_rows([1e12 + x for x in [0, .005, .01, .015, .02]]),
        size=(1000, 600))
    box("large_offset_ulp_spaced_labels", source_rows([1e16 + x for x in [0, 2, 4, 6, 8]]),
        orientation=1, size=(1200, 600))
    long_ids = source_rows([1, 2, 3, 4], group=9007199254740000)
    long_ids += source_rows([2, 3, 4, 6], group=9007199254740002, first_id=20)
    hist("long_group_identifiers", long_ids, [0, 2, 4, 6], size=(1200, 800))
    box("long_horizontal_group_identifiers", long_ids, orientation=1, size=(1200, 600))
    hist("horizontal_histogram", source_rows([-2, -1, -.5, 0, .5, 1, 2]), [-2, -1, 0, 1, 2],
         orientation=1)
    box("horizontal_box_with_outlier", source_rows([0, 1, 2, 3, 4, 5, 100]), orientation=1)
    hist("resized_histogram", source_rows([0, .5, 1, 2, 3]), [0, 1, 3],
         size=(1200, 700), render=False)
    box("resized_box", source_rows([1, 2, 3, 4]), size=(1200, 700), render=False)
    transformed = [[4 * r[0] + 7, *r[1:]] for r in source_rows([0, .5, 1, 2, 3])]
    hist("affine_histogram_density", transformed, [7, 11, 19], scale=2, render=False)
    transformed = [[4 * r[0] + 7, *r[1:]] for r in source_rows([0, 1, 2, 3, 4, 5, 100])]
    box("affine_box_summary_and_fences", transformed, render=False)
    shuffled = copy.deepcopy(states)
    random.Random(20261022).shuffle(shuffled)
    hist("permuted_histogram_rows", shuffled, [-10, 0, 2, 10, 101], render=False, size=(1000, 800))
    box("permuted_box_rows", shuffled, render=False, size=(1000, 800))
    return output


def malformed_cases():
    valid = source_rows([0, 1, 2])
    hist = [1, 800, 500, 0, 0, 3, 0, 1, 2]
    box = [1, 800, 500, 0, 0, 1, 1]
    output = []
    def add(name, rows=valid, params=hist, op="report_histogram_svg", sql_only=False):
        output.append({"name": name, "operation": op, "data": rows,
                       "parameters": params, "sql_only": sql_only})
    for name, column, value in [("negative_frequency", 2, -1), ("fractional_frequency", 2, .5),
                                 ("nonboolean_excluded", 3, 2), ("nonboolean_hidden", 4, -.5),
                                 ("nonboolean_observed", 5, 2), ("fractional_group", 1, .5),
                                 ("negative_row_id", 6, -1), ("fractional_row_id", 6, .5),
                                 ("frequency_budget", 2, 1000001),
                                 ("nonfinite_value", 0, float("nan"))]:
        rows = copy.deepcopy(valid)
        rows[0][column] = value
        add(name, rows)
    duplicate = copy.deepcopy(valid)
    duplicate[1][6] = duplicate[0][6]
    add("duplicate_row_id", duplicate)
    add("empty_matrix", [])
    add("wrong_column_count", [r[:-1] for r in valid])
    add("extra_column", [r + [0] for r in valid])
    add("no_observations", [[0, 0, 1, 0, 0, 0, 1]])
    add("all_excluded", [[0, 0, 1, 1, 0, 1, 1]])
    add("all_zero_frequency", [[0, 0, 0, 0, 0, 1, 1]])
    add("too_many_groups", [r for g in range(17) for r in source_rows([0, 1], group=g, first_id=g * 10)])
    for name, params in [("wrong_version", [2, *hist[1:]]),
                         ("small_canvas", [1, 50, 50, *hist[3:]]),
                         ("unsupported_scale", [1, 800, 500, 3, 0, 3, 0, 1, 2]),
                         ("bad_orientation", [1, 800, 500, 0, 2, 3, 0, 1, 2]),
                         ("duplicate_edges", [1, 800, 500, 0, 0, 3, 0, 1, 1]),
                         ("decreasing_edges", [1, 800, 500, 0, 0, 3, 2, 1, 0]),
                         ("truncated_edges", hist[:-1]),
                         ("selected_outside_domain", [1, 800, 500, 0, 0, 3, 0, 1, 1.5]),
                         ("nonfinite_parameter", [1, 800, 500, float("inf"), *hist[4:]]),
                         ("unreadably_many_bins", [1, 320, 240, 0, 0, 513, *[i / 256 for i in range(513)]])]:
        add(name, params=params)
    add("invalid_box_type", params=[1, 800, 500, 2, 0, 1, 1], op="report_box_plot_svg")
    add("invalid_box_visibility", params=[1, 800, 500, 0, 0, 2, 1], op="report_box_plot_svg")
    add("truncated_box_arguments", params=box[:-1], op="report_box_plot_svg")
    add("nested_SQL_NULL_cell", [[None, *valid[0][1:]]], sql_only=True)
    add("nested_SQL_NULL_row", [None], sql_only=True)
    add("top_SQL_NULL_matrix", None, sql_only=True)
    add("top_SQL_NULL_parameters", params=None, sql_only=True)
    add("ragged_SQL_matrix", [valid[0], valid[1][:-1]], sql_only=True)
    add("unknown_operation", op="report_not_a_real_operation")
    return output


def save_artifact(case, layer, result, artifacts, magick, renderer):
    root, anchors = verify_result(case, result)
    destination = artifacts / layer
    destination.mkdir(parents=True, exist_ok=True)
    path = destination / case["name"]
    svg_path = path.with_suffix(".svg")
    json_path = path.with_suffix(".json")
    svg_path.write_text(result["svg"], encoding="utf-8")
    json_path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    record = {"name": case["name"], "operation": case["operation"], "boundary": layer,
              "input": case["data"], "parameters": case["parameters"],
              "result_json": input_sha(json_path), "svg": input_sha(svg_path),
              "xml_verified": True, "rasterized": False}
    if case["render"]:
        record["render"] = rasterize(svg_path, path.with_suffix(".png"), magick,
                                      int(result["width"]), int(result["height"]), anchors, renderer)
        record["rasterized"] = True
    # These mutations prove that the independent artifact parser rejects actual
    # unsafe/nonfinite payloads, without treating a handcrafted SVG as output.
    if case["name"] == "type6_not_type7":
        assert_rejected_svg(result["svg"], result["width"], result["height"],
                            lambda n: ET.SubElement(n, "{" + SVG_NS + "}script"))
        assert_rejected_svg(result["svg"], result["width"], result["height"],
                            lambda n: n.set("onload", "alert(1)"))
        assert_rejected_svg(result["svg"], result["width"], result["height"],
                            lambda n: ET.SubElement(n, "{" + SVG_NS + "}rect", x="NaN", y="0", width="1", height="1"))
        record["malicious_svg_mutations_rejected"] = 3
    return record


def sql_vectors(runner):
    records = []
    for operation in ["report_histogram_svg", "report_box_plot_svg"]:
        if operation == "report_histogram_svg":
            matrix = "[[i::DOUBLE,0.,1.,0.,0.,1.,i::DOUBLE],[i+1.,0.,1.,0.,0.,1.,i+100.]]"
            parameters = "[1.,800.,500.,0.,0.,3.,i-1.,i+.5,i+2.]"
        else:
            matrix = "[[i-1.,0.,1.,0.,0.,1.,i::DOUBLE],[i+1.,0.,1.,0.,0.,1.,i+100.]]"
            parameters = "[1.,800.,500.,0.,0.,1.,1.]"
        rows = runner.query(
            "SELECT i,chem_analyze('" + operation + "',CASE WHEN i%3=0 THEN NULL::DOUBLE[][] ELSE "
            + matrix + " END," + parameters + ") AS result FROM range(8) t(i) WHERE i<>4 ORDER BY i;")
        runner.calls += len(rows)
        same([r["i"] for r in rows], [0, 1, 2, 3, 5, 6, 7], "actual SQL selected-vector rows")
        require(sum(r["result"] is None for r in rows) == 3, "Actual SQL selected-vector NULL bitmap")
        for row in rows:
            i = row["i"]
            if i % 3 == 0:
                require(row["result"] is None, "Selected SQL NULL row became a report")
                continue
            values = [i, i + 1] if operation == "report_histogram_svg" else [i - 1, i + 1]
            data = source_rows(values)
            data[0][6], data[1][6] = i, i + 100
            params = ([1, 800, 500, 0, 0, 3, i - 1, i + .5, i + 2]
                      if operation == "report_histogram_svg" else [1, 800, 500, 0, 0, 1, 1])
            case = {"name": operation + "_selected_" + str(i), "operation": operation,
                    "data": data, "parameters": params, "render": False}
            result = strict_json(row["result"])
            verify_result(case, result)
            records.append({"operation": operation, "selected_row": i, "report_N": result["N"],
                            "source_row_ids": [i, i + 100], "xml_verified": True})
    return records


def alias_checks(runner, all_cases):
    names = ["chem_report_box_plot_svg", "chem_report_histogram_svg"]
    catalog = runner.query(
        "SELECT function_name,description,parameters,parameter_types,examples FROM duckdb_functions() "
        "WHERE function_name IN ('chem_report_histogram_svg','chem_report_box_plot_svg') ORDER BY function_name;")
    same([row["function_name"] for row in catalog], names, "dedicated SQL aliases registered")
    for row in catalog:
        same(row["parameters"], ["data", "parameters"], "alias named parameters")
        same(row["parameter_types"], ["DOUBLE[][]", "DOUBLE[]"], "alias parameter types")
        require(row["description"] and row["examples"], "alias documentation/examples")
    for operation in ["report_histogram_svg", "report_box_plot_svg"]:
        case = next(c for c in all_cases if c["operation"] == operation)
        controlled = "chem_" + operation + "(" + sql_value(case["data"]) + "::DOUBLE[][]," + sql_value(case["parameters"]) + "::DOUBLE[])"
        rows = runner.query("SELECT " + controlled + " AS result;")
        runner.calls += 1
        result = strict_json(rows[0]["result"])
        verify_result(case, result)
    return {"catalog": catalog, "controlled_alias_calls": 2}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--library", type=Path)
    parser.add_argument("--duckdb", type=Path)
    parser.add_argument("--extension", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--magick", default=shutil.which("magick"))
    parser.add_argument("--rasterizer", choices=["magick", "sips"], default="magick",
                        help="Renderer for the original SVG; ImageMagick independently decodes either PNG")
    parser.add_argument("--require-aliases", action="store_true")
    parser.add_argument("--reference-only", action="store_true")
    parser.add_argument("--case", action="append", default=[],
                        help="Run a named focused actual-output regression; default runs every positive fixture")
    args = parser.parse_args()
    require(bool(args.duckdb) == bool(args.extension), "--duckdb and --extension must be supplied together")
    require(args.reference_only or args.library or args.duckdb,
            "Supply a real public boundary or explicitly choose --reference-only")
    require(not args.require_aliases or args.duckdb,
            "--require-aliases needs an actual --duckdb/--extension boundary")
    require(not args.reference_only or not (args.library or args.duckdb),
            "Reference-only mode cannot be mixed with actual boundary flags")
    args.output = args.output.resolve()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    output = {"schema": "jmp-report-independent-verification-v1", "status": "RUNNING",
              "checked_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
              "python": sys.version, "platform": platform.platform(), "primary_sources": PRIMARY,
              "verifier": input_sha(__file__), "reference_only": args.reference_only,
              "reference_method": "11 hand-derived literal cases; origin-centered standard Python statistics.quantiles exclusive with JMP endpoint clamps; tiny physically expanded integer frequency samples; independently counted intervals; Decimal90 audit of large-offset rounding",
              "verification_boundaries": [], "cases": [], "nullable_cases": [],
              "call_count_method": "C ABI counts native invocations. SQL counts requested scalar evaluations/selected rows, not unobservable optimizer callbacks.",
              "visual_qa": {"status": "pending_human_or_agent_image_inspection",
                            "note": "XML/pixel checks do not substitute for viewing the produced artifacts."},
              "whole_JMP_feature_equivalence_claim": False}
    current_case = None
    current_result = None
    try:
        literal_reference_checks()
        output["literal_reference_fixtures"] = 11
        all_cases = cases()
        output["large_offset_reference_audit"] = large_offset_reference_audit(all_cases)
        require(len({c["name"] for c in all_cases}) == len(all_cases), "Unique fixture names")
        output["independent_expected_summaries"] = []
        for case in all_cases:
            groups = []
            for gid in sorted({r[1] for r in case["data"]}):
                selected = selected_rows(case["data"], gid)
                expected = {"group_id": gid, **state_counts([r for r in case["data"] if r[1] == gid])}
                if selected:
                    expected.update(histogram_oracle(selected, case["parameters"][6:])
                                    if case["operation"] == "report_histogram_svg" else box_oracle(selected))
                groups.append(expected)
            output["independent_expected_summaries"].append({"name": case["name"], "groups": groups})
        if args.case:
            require(set(args.case) <= {c["name"] for c in all_cases}, "Unknown requested fixture name")
            all_cases = [c for c in all_cases if c["name"] in args.case]
        output["full_positive_suite"] = not bool(args.case)
        if not args.reference_only:
            require(args.magick is not None, "ImageMagick required for actual CPU artifact verification")
            require(args.rasterizer != "sips" or Path("/usr/bin/sips").is_file(),
                    "Explicit sips rasterizer requires the existing macOS command")
            version = subprocess.run([args.magick, "-version"], capture_output=True, text=True, check=True)
            output["renderer"] = {"path": args.magick, "version": version.stdout.splitlines()[0],
                                  "svg_rasterizer": args.rasterizer,
                                  "original_svg_used_without_rewriting": True,
                                  "png_decoder": "ImageMagick_RGBA"}
            boundaries = {}
            if args.library:
                boundaries["public_c_abi"] = Abi(args.library)
                output["library"] = input_sha(args.library)
            if args.duckdb:
                boundaries["native_sql"] = Sql(args.duckdb, args.extension)
                output["duckdb"] = input_sha(args.duckdb)
                output["extension"] = input_sha(args.extension)
            output["verification_boundaries"] = list(boundaries)
            actual = {}
            for layer, boundary in boundaries.items():
                if isinstance(boundary, Sql):
                    results = boundary.many([(c["operation"], c["data"], c["parameters"]) for c in all_cases])
                else:
                    results = [boundary.nullable(c["operation"], c["data"], c["parameters"]) for c in all_cases]
                actual[layer] = results
                for case, result in zip(all_cases, results):
                    current_case, current_result = case, result
                    output["cases"].append(save_artifact(case, layer, result,
                                                         args.output.parent / "artifacts", args.magick,
                                                         args.rasterizer))
                invalid = [c for c in malformed_cases() if not c["sql_only"] or isinstance(boundary, Sql)]
                if isinstance(boundary, Sql):
                    invalid_results = boundary.many([(c["operation"], c["data"], c["parameters"]) for c in invalid])
                else:
                    invalid_results = [boundary.nullable(c["operation"], c["data"], c["parameters"]) for c in invalid]
                for case, result in zip(invalid, invalid_results):
                    require(result is None, (layer, "malformed public input accepted", case["name"], result))
                    output["nullable_cases"].append({"name": case["name"], "boundary": layer, "returned_NULL": True})
                repeat = boundary.nullable(all_cases[0]["operation"], all_cases[0]["data"], all_cases[0]["parameters"])
                require(repeat == results[0], (layer, "Repeated actual report call is nondeterministic"))
                if isinstance(boundary, Sql):
                    output["selected_SQL_vectors"] = sql_vectors(boundary)
                    if args.require_aliases:
                        output["aliases"] = alias_checks(boundary, all_cases)
                    output["SQL_queries"] = boundary.queries
                output[layer + "_calls"] = boundary.calls
            if set(actual) == {"public_c_abi", "native_sql"}:
                require(actual["public_c_abi"] == actual["native_sql"],
                        "Actual SQL and public C ABI produced different report JSON/SVG")
                output["C_ABI_SQL_full_output_equal"] = True
            output["status"] = "PASS"
        else:
            output["status"] = "REFERENCE_ONLY"
        output["positive_fixtures"] = len(all_cases)
        output["assertions"] = CHECKS
        output["full_public_boundaries_verified"] = set(output["verification_boundaries"]) == {"public_c_abi", "native_sql"}
        output["rasterized_artifacts"] = sum(c["rasterized"] for c in output["cases"])
    except Exception:
        output["status"] = "FAIL"
        output["assertions_before_failure"] = CHECKS
        output["failure"] = traceback.format_exc()
        if current_case is not None:
            output["failure_case"] = current_case
        if isinstance(current_result, dict):
            json_path = args.output.with_name(args.output.stem + "-failure-result.json")
            json_path.write_text(json.dumps(current_result, indent=2, allow_nan=False) + "\n")
            output["failure_result"] = input_sha(json_path)
            if isinstance(current_result.get("svg"), str):
                svg_path = args.output.with_name(args.output.stem + "-failure.svg")
                svg_path.write_text(current_result["svg"])
                output["failure_SVG"] = input_sha(svg_path)
        args.output.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
        raise
    args.output.write_text(json.dumps(output, indent=2, allow_nan=False) + "\n")
    print(json.dumps({"status": output["status"], "positive_fixtures": output["positive_fixtures"],
                      "boundaries": output["verification_boundaries"], "assertions": CHECKS,
                      "rasterized_artifacts": output["rasterized_artifacts"], "output": str(args.output)}))


if __name__ == "__main__":
    main()
