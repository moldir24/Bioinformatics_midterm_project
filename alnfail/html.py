"""Turn the report folder into one page a person can read: `report.html`.

`summary.md` is exact but tiring: wide tables of numbers, and the figures are
separate files. This module converts it into a single self-contained HTML file
(the figures are embedded, so the one file can be downloaded, e-mailed or
printed to PDF). Nothing is recomputed: every number on the page is a number
from `summary.md`. What is added is presentation:

    * the figures appear inside the section they belong to;
    * each "wrong per 1,000 (interval)" value is also drawn on a log-scale
      ruler next to the rate the aligner claimed, with a three-way verdict
      (above the claim / within the claim / too few reads to tell);
    * the failure-map tables are tinted by value, failed checks are marked;
    * every section starts with two sentences on how to read it.

    python -m alnfail html --report results/report
    python -m alnfail.html results/report            # the same, without the main CLI

The Markdown reader below understands only what `report.py` writes (headings,
paragraphs, pipe tables, bullet lists, `code`, **bold**, _italic_). It is not
a general Markdown converter and does not try to be.
"""
from __future__ import annotations

import base64
import datetime
import html
import math
import os
import re
import sys

TITLE = "Where do aligners actually fail?"

# the same colours as the figures (alnfail/plots.py), so that a name in a table
# and a line in a plot are visibly the same thing
ALIGNER_COLOUR = {"minimap2": "#2a78d6", "bwamem2": "#eb6834", "bowtie2": "#1baf7a", "strobealign": "#eda100",
                  "winnowmap": "#e87ba4", "bwa": "#008300", "ngmlr": "#4a3aa7"}

# which figures belong to which numbered section of summary.md (by file-name prefix)
FIGURE_SECTION = (("fig_accuracy", 1), ("fig_divergence", 1), ("fig_mapq", 2), ("fig_strata", 3),
                  ("fig_calibration", 4), ("fig_transfer", 5), ("fig_resources", 7))
FIGURES_FIRST = {3, 4, 5, 7}  # sections where the picture is easier to take in than the table

HOW_TO_READ = {
    1: "One row is one aligner on one sequencing platform; correct, misplaced and unmapped add up to 100 %. "
       "The last columns look only at reads the aligner called confident (MAPQ 30 or more), where it promises "
       "at most 1 wrong alignment per 1,000.",
    2: "The same question band by band. A band of reported MAPQ claims an error rate no higher than that of its "
       "lower edge: MAPQ 20 claims 10 per 1,000, MAPQ 40 claims 0.1 per 1,000. An honest aligner stays at or "
       "left of the line in every band.",
    3: "Reads are grouped by where in the genome they truly came from. Darker cells mean more reads placed in "
       "the wrong position. Read along a row to compare aligners in the same context, down a column to see "
       "what one aligner finds hard.",
    4: "The simulator is only useful if its reads look like real ones. Each row is a property measured the "
       "same way on real and on simulated reads; rows marked as outside tolerance are properties the "
       "simulation does not reproduce and belong in the limitations.",
    5: "Real reads have no known origin, so accuracy cannot be measured on them directly. These tables compare "
       "what can be measured without truth, and whether the genomic contexts that are hard in simulation are "
       "also the ones where aligners disagree on real reads.",
    6: "Every filter applied to the real reads before they were used, the share of reads or bases it removed, "
       "and the reason the filter exists.",
    7: "Wall-clock time and peak memory of the index and mapping jobs, as measured by the workflow.",
}
RULER_LEGEND = ("On the ruler, the dot is the measured rate, the bar is its 95 % interval and the vertical line is "
                "the claim. The scale is logarithmic, from 0.001 to 1,000 wrong per 1,000. A bar with no dot means "
                "no wrong alignment was seen: the bar then shows how high the true rate could still be. Undecided means the "
                "interval crosses the line, so this many reads cannot settle the question.")

NICE_HEADER = {
    "stratum_set": "context type", "stratum": "context", "agree_pct": "agree %",
    "n_compared": "reads compared", "confident_disagree_per_1000": "confident disagreements per 1,000",
    "strata_compared": "contexts compared", "proxy_validity_spearman": "proxy validity (Spearman)",
    "transfer_spearman": "transfer (Spearman)", "median_ratio_real_over_sim": "median ratio, real ÷ simulated",
    "simulated_value": "simulated value", "real_value": "real value", "median_seconds": "median seconds",
    "total_cpu_seconds": "CPU seconds, total", "peak_rss_mb": "peak memory, MB",
}
WORD_COLUMNS = {"stratum_set", "stratum", "metric", "step", "kind"}  # identifiers shown as words
VERDICT = {"above": ("✕", "above the claim"), "within": ("✓", "within the claim"), "open": ("?", "undecided")}
# the ruler and its verdict say the same thing as this column (the claim there is 1, so the ratio equals the rate)
REDUNDANT_COLUMNS = {"× the promised 1 per 1,000"}
TALL_TABLE = 14  # tables longer than this scroll inside a box with a fixed header

LOG_MIN, LOG_MAX = 0.001, 1000.0  # ruler domain, wrong alignments per 1,000
RULER_W, RULER_H, RULER_PAD = 150, 22, 5


# ---------------------------------------------------------------- reading summary.md
def _split_row(line: str) -> list[str]:
    line = line.strip()
    if line.startswith("|"):
        line = line[1:]
    if line.endswith("|"):
        line = line[:-1]
    return [cell.strip() for cell in line.split("|")]


def parse_markdown(text: str) -> list[tuple]:
    """Blocks of the Markdown subset: ('h1', s), ('h2', s), ('p', s), ('ul', [s]), ('table', header, rows)."""
    blocks: list[tuple] = []
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        line = lines[i]
        if not line.strip():
            i += 1
        elif line.startswith("## "):
            blocks.append(("h2", line[3:].strip()))
            i += 1
        elif line.startswith("# "):
            blocks.append(("h1", line[2:].strip()))
            i += 1
        elif line.startswith("|"):
            rows = []
            while i < len(lines) and lines[i].startswith("|"):
                rows.append(_split_row(lines[i]))
                i += 1
            is_rule = len(rows) > 1 and all(set(cell) <= set("-: ") and cell for cell in rows[1])
            blocks.append(("table", rows[0], rows[2:] if is_rule else rows[1:]))
        elif line.startswith("- "):
            items = []
            while i < len(lines) and lines[i].startswith("- "):
                items.append(lines[i][2:].strip())
                i += 1
            blocks.append(("ul", items))
        else:
            words = []
            while i < len(lines) and lines[i].strip() and not lines[i].startswith(("#", "|", "- ")):
                words.append(lines[i].strip())
                i += 1
            blocks.append(("p", " ".join(words)))
    return blocks


def inline(text: str) -> str:
    """Escape, then `code`, **bold** and _italic_ (underscores inside identifiers are left alone)."""
    out = html.escape(text, quote=False)
    out = re.sub(r"`([^`]+)`", r"<code>\1</code>", out)
    out = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", out)
    return re.sub(r"(?<![\w])_([^_]+)_(?![\w])", r"<em>\1</em>", out)


def _slug(text: str) -> str:
    return re.sub(r"[^a-z0-9]+", "-", text.lower()).strip("-") or "section"


# ---------------------------------------------------------------- numbers
_NUMBER = re.compile(r"^-?[\d,]*\.?\d+(?:e[+-]?\d+)?%?$", re.IGNORECASE)
_RATE = re.compile(r"^\s*([\d.]+)\s*\(([\d.]+)[–-]([\d.]+)\)\s*$")
_SUPERSCRIPT = str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹")


def is_number(cell: str) -> bool:
    return bool(_NUMBER.match(cell.strip()))


def to_float(cell: str) -> float | None:
    try:
        return float(cell.strip().rstrip("%").replace(",", ""))
    except ValueError:
        return None


def pretty_number(cell: str) -> str:
    """'1e+03' -> '1,000' and '1.707e-05' -> '1.71 × 10⁻⁵'; everything else is returned as written."""
    match = re.match(r"^(-?[\d.]+)e([+-]?\d+)$", cell.strip(), re.IGNORECASE)
    if not match:
        return cell
    value = float(cell)
    if abs(value) >= 1:
        return f"{value:,.0f}"
    mantissa, exponent = f"{value:.2e}".split("e")
    return f"{float(mantissa):g} × 10{str(int(exponent)).translate(_SUPERSCRIPT)}"


def parse_rate(cell: str) -> tuple[float, float, float] | None:
    """'0.62 (0.43–0.89)' -> (0.62, 0.43, 0.89); None for anything else ('no reads')."""
    match = _RATE.match(cell)
    return tuple(float(g) for g in match.groups()) if match else None


def verdict(low: float, high: float, claim: float) -> str:
    """Compare a 95 % interval with a claimed upper limit.

    above    the whole interval is above the claim: the claim is broken
    within   the whole interval is at or below it: the claim holds
    open     the interval straddles the claim: this many reads cannot decide
    """
    if low > claim:
        return "above"
    if high <= claim:
        return "within"
    return "open"


# ---------------------------------------------------------------- the ruler
def _x(value: float) -> float:
    value = min(max(value, LOG_MIN), LOG_MAX)
    span = math.log10(LOG_MAX) - math.log10(LOG_MIN)
    return RULER_PAD + (math.log10(value) - math.log10(LOG_MIN)) / span * (RULER_W - 2 * RULER_PAD)


def ruler(rate: tuple[float, float, float], claim: float) -> str:
    """Inline SVG: interval bar and point estimate on a log axis, with the claimed rate as a vertical line."""
    estimate, low, high = rate
    state = verdict(low, high, claim)
    mid = RULER_H / 2
    # an interval that runs past either end of the axis is drawn to the very edge of the picture
    x1 = 0 if low <= LOG_MIN else _x(low)
    x2 = RULER_W if high >= LOG_MAX else _x(high)
    label = (f"{estimate:g} wrong per 1,000, 95 % interval {low:g} to {high:g}; "
             f"claimed at most {claim:g}: {VERDICT[state][1]}")
    ticks = "".join(f'<line x1="{_x(10.0 ** e):.1f}" x2="{_x(10.0 ** e):.1f}" y1="{RULER_H - 4}" y2="{RULER_H}" class="tick"/>'
                    for e in range(int(math.log10(LOG_MIN)), int(math.log10(LOG_MAX)) + 1))
    dot = f'<circle cx="{_x(estimate):.1f}" cy="{mid}" r="4" class="dot"/>' if estimate > 0 else ""
    return (f'<svg class="ruler {state}" width="{RULER_W}" height="{RULER_H}" viewBox="0 0 {RULER_W} {RULER_H}" '
            f'role="img" aria-label="{html.escape(label)}"><title>{html.escape(label)}</title>{ticks}'
            f'<line x1="{x1:.1f}" x2="{max(x2, x1 + 2):.1f}" y1="{mid}" y2="{mid}" class="bar"/>'
            f'<line x1="{_x(claim):.1f}" x2="{_x(claim):.1f}" y1="1" y2="{RULER_H - 1}" class="claim"/>{dot}</svg>')


def ruler_axis() -> str:
    """The scale shown once, in the column header."""
    labels = "".join(f'<text x="{_x(v):.1f}" y="9" text-anchor="{anchor}">{text}</text>'
                     for v, text, anchor in ((LOG_MIN, "0.001", "start"), (1.0, "1", "middle"), (LOG_MAX, "1,000", "end")))
    return (f'<svg class="axis" width="{RULER_W}" height="12" viewBox="0 0 {RULER_W} 12" aria-hidden="true">'
            f"{labels}</svg>")


def verdict_badge(state: str) -> str:
    icon, text = VERDICT[state]
    return f'<span class="verdict {state}"><span class="icon" aria-hidden="true">{icon}</span> {text}</span>'


# ---------------------------------------------------------------- tables
def _swatch(name: str) -> str:
    colour = ALIGNER_COLOUR.get(name.strip())
    return f'<span class="swatch" style="background:{colour}"></span>' if colour else ""


def _header_text(name: str) -> str:
    return NICE_HEADER.get(name, name.replace("_", " "))


def _words(cell: str) -> str:
    return cell.replace(">=", "≥ ").replace("<", "< ").replace("_", " ").replace("  ", " ").strip()


def is_heat_table(header: list[str]) -> bool:
    """The failure-map tables: context type, context, then one column of rates per aligner."""
    return header[:2] == ["stratum_set", "stratum"]


def heat_values(rows: list[list[str]]) -> list[float]:
    return [to_float(cell) or 0.0 for row in rows for cell in row[2:]]


def render_table(header: list[str], rows: list[list[str]], heat_max: float = 0.0) -> str:
    """One table of summary.md as HTML, with the extras described at the top of this file.

    heat_max is the largest rate over all failure-map tables of the page, so that
    the same tint means the same rate in every one of them.
    """
    rows = [row + [""] * (len(header) - len(row)) for row in rows]
    keep = [i for i, name in enumerate(header) if name not in REDUNDANT_COLUMNS]
    header, rows = [header[i] for i in keep], [[row[i] for i in keep] for row in rows]
    n = len(header)
    column = {name: i for i, name in enumerate(header)}
    numeric = [bool(rows) and all(is_number(r[i]) or not r[i] for r in rows) and any(r[i] for r in rows)
               for i in range(n)]

    # which column holds "wrong per 1,000 (interval)", and what each row claims
    rate_col = next((i for i, name in enumerate(header) if name.startswith("wrong per 1,000")), None)
    claim_col = column.get("claimed at most per 1,000")
    heat = is_heat_table(header)
    heat_max = heat_max or max(heat_values(rows) + [0.0])
    # checks that failed are the reason to look at such a table: show them first
    for flag in ("ok", "transfers"):
        if flag in column:
            rows = sorted(rows, key=lambda r: r[column[flag]] != "False")

    head = []
    for i, name in enumerate(header):
        cls = ' class="num"' if numeric[i] or (heat and i >= 2) else ""
        head.append(f"<th{cls}>{_swatch(name) if heat else ''}{inline(_header_text(name))}</th>")
        if i == rate_col:
            head.append(f'<th class="rulerhead">against the claim{ruler_axis()}</th><th>verdict</th>')

    body = []
    for row in rows:
        cells, flagged = [], False
        claim = 1.0 if claim_col is None else (to_float(row[claim_col]) or 0.0)
        for i, cell in enumerate(row):
            name = header[i]
            if name in ("ok", "transfers") and cell in ("True", "False"):
                good = cell == "True"
                flagged = flagged or not good
                text = "✓ yes" if good else "✕ no"
                cells.append(f'<td class="{"yes" if good else "no"}">{text}</td>')
            elif heat and i >= 2:
                value = to_float(cell)
                alpha = 0.6 * math.sqrt(value / heat_max) if value and heat_max > 0 else 0.0
                style = f' style="background:rgba(42,120,214,{alpha:.3f})"' if alpha else ""
                cells.append(f'<td class="num"{style}>{html.escape(cell)}</td>')
            elif name == "aligner":
                cells.append(f'<td class="name">{_swatch(cell)}{html.escape(cell)}</td>')
            elif name in WORD_COLUMNS:
                cells.append(f"<td>{html.escape(_words(cell))}</td>")
            elif numeric[i]:
                cells.append(f'<td class="num">{html.escape(pretty_number(cell))}</td>')
            elif i == rate_col:
                cells.append(f'<td class="rate">{html.escape(cell)}</td>')
            else:
                cells.append(f"<td>{inline(cell)}</td>")
            if i == rate_col:
                rate = parse_rate(cell)
                if rate and claim > 0:
                    state = verdict(rate[1], rate[2], claim)
                    flagged = flagged or state == "above"
                    cells.append(f'<td class="rulercell">{ruler(rate, claim)}</td><td>{verdict_badge(state)}</td>')
                else:
                    cells.append("<td></td><td></td>")
        row_class = ' class="flagged"' if flagged else ""
        body.append(f"<tr{row_class}>{''.join(cells)}</tr>")

    tall = len(rows) > TALL_TABLE
    note = f'<p class="rows">{len(rows)} rows; scroll inside the table to see them all.</p>' if tall else ""
    return (f'<div class="tablebox{" tall" if tall else ""}"><table><thead><tr>{"".join(head)}</tr></thead>'
            f'<tbody>{"".join(body)}</tbody></table></div>{note}')


# ---------------------------------------------------------------- figures
def _figure(report_dir: str, name: str, caption: str) -> str:
    with open(os.path.join(report_dir, name), "rb") as fh:
        data = base64.b64encode(fh.read()).decode("ascii")
    return (f'<figure><img src="data:image/png;base64,{data}" alt="{html.escape(caption)}">'
            f"<figcaption>{inline(caption)}</figcaption></figure>")


def _figure_section(name: str) -> int | None:
    return next((number for prefix, number in FIGURE_SECTION if name.startswith(prefix)), None)


# ---------------------------------------------------------------- the page
CSS = """
:root { --surface:#fcfcfb; --ink:#0b0b0b; --ink2:#52514e; --muted:#898781; --rule:#e1e0d9; --wash:#f3f2ee;
        --blue:#1c5cab; --alarm:#b3261e; --prose:42rem;
        --serif: Charter, "Bitstream Charter", "Iowan Old Style", "Sitka Text", Cambria, Georgia, serif;
        --sans: "Segoe UI", "Helvetica Neue", Helvetica, Arial, "Noto Sans", sans-serif; }
* { box-sizing:border-box; }
html { color-scheme:light; scroll-behavior:smooth; }
body { margin:0; background:var(--surface); color:var(--ink); font:17px/1.6 var(--serif); }
.page { display:grid; grid-template-columns:13rem minmax(0,1fr); gap:2.5rem; max-width:88rem; margin:0 auto; padding:0 2rem 6rem; }
nav { position:sticky; top:0; align-self:start; max-height:100vh; overflow:auto; padding:3.5rem 0 2rem; font:14px/1.4 var(--sans); }
nav p { margin:0 0 .75rem; color:var(--muted); }
nav ol { list-style:none; margin:0; padding:0; border-left:1px solid var(--rule); }
nav a { display:block; padding:.35rem 0 .35rem .9rem; margin-left:-1px; border-left:1px solid transparent; color:var(--ink2); text-decoration:none; }
nav a:hover, nav a:focus-visible { color:var(--ink); border-left-color:var(--ink); }
main { min-width:0; padding-top:3.5rem; }
header h1 { font:600 2.6rem/1.1 var(--sans); letter-spacing:-.02em; margin:0 0 1rem; }
header p { max-width:var(--prose); margin:.5rem 0; }
header .built { color:var(--ink2); font:14px/1.5 var(--sans); }
section { margin-top:4rem; }
h2 { font:600 1.5rem/1.25 var(--sans); letter-spacing:-.01em; margin:0 0 .75rem; padding-top:1.25rem; border-top:2px solid var(--ink); }
section > p, section > ul { max-width:var(--prose); }
p.how { color:var(--ink2); }
p.subhead { margin:2.25rem 0 .5rem; max-width:none; font:15px/1.4 var(--sans); }
code { font:.86em ui-monospace, "SF Mono", Menlo, Consolas, monospace; background:var(--wash); padding:.08em .3em; border-radius:3px; }
a { color:var(--blue); }
:focus-visible { outline:2px solid var(--blue); outline-offset:2px; }
.tablebox { overflow-x:auto; margin:1rem 0 1.5rem; border-top:1px solid var(--ink); border-bottom:1px solid var(--rule); }
.tablebox.tall { max-height:32rem; overflow-y:auto; }
table { border-collapse:collapse; width:100%; font:14px/1.35 var(--sans); font-variant-numeric:tabular-nums; }
th, td { padding:.42rem .6rem; text-align:left; vertical-align:middle; }
th { font-weight:600; color:var(--ink2); background:var(--surface); border-bottom:1px solid var(--ink); vertical-align:bottom; position:sticky; top:0; z-index:1; }
td { border-bottom:1px solid var(--rule); }
tbody tr:last-child td { border-bottom:0; }
p.rows { margin:-1rem 0 1.5rem; color:var(--muted); font:13px/1.4 var(--sans); }
.num { text-align:right; white-space:nowrap; }
th.num { white-space:normal; }
td.name, td.rate { white-space:nowrap; }
td.yes { color:var(--ink2); white-space:nowrap; }
td.no { color:var(--alarm); font-weight:600; white-space:nowrap; }
tr.flagged td { background:rgba(179,38,30,.055); }
.swatch { display:inline-block; width:.62em; height:.62em; border-radius:50%; margin-right:.45em; }
.rulerhead { white-space:nowrap; }
.rulerhead svg { display:block; margin-top:.2rem; }
.axis text { font:10px var(--sans); fill:var(--muted); }
.rulercell { padding-top:.2rem; padding-bottom:.2rem; }
.ruler { display:block; overflow:visible; }
.ruler .tick { stroke:var(--rule); stroke-width:1; }
.ruler .claim { stroke:var(--ink); stroke-width:1.5; }
.ruler .bar { stroke-width:4; stroke-linecap:butt; }
.ruler .dot { stroke:var(--surface); stroke-width:1.5; }
.ruler.above .bar { stroke:var(--alarm); } .ruler.above .dot { fill:var(--alarm); }
.ruler.within .bar { stroke:var(--blue); } .ruler.within .dot { fill:var(--blue); }
.ruler.open .bar { stroke:var(--muted); } .ruler.open .dot { fill:var(--muted); }
.verdict { white-space:nowrap; }
.verdict .icon { display:inline-block; width:1.1em; font-weight:700; }
.verdict.above { color:var(--alarm); font-weight:600; }
.verdict.within { color:var(--blue); }
.verdict.open { color:var(--ink2); }
figure { margin:1.5rem 0 2.5rem; }
figure img { display:block; max-width:100%; height:auto; }
figcaption { max-width:var(--prose); margin-top:.5rem; color:var(--ink2); font:14px/1.45 var(--sans); }
@media (max-width:60rem) {
  .page { display:block; padding:0 1.1rem 4rem; }
  nav { position:static; max-height:none; padding:1.5rem 0 0; }
  main { padding-top:1.5rem; } header h1 { font-size:2rem; }
}
@media (prefers-reduced-motion:reduce) { html { scroll-behavior:auto; } }
@media print {
  .page { display:block; max-width:none; padding:0; } nav { display:none; }
  .tablebox.tall { max-height:none; overflow:visible; } th { position:static; } p.rows { display:none; }
  section { margin-top:2rem; } figure, tr { break-inside:avoid; }
}
"""


def build_html(report_dir: str, out_path: str | None = None) -> str:
    """Write report.html from summary.md and the PNG figures in `report_dir`; returns the path written."""
    summary = os.path.join(report_dir, "summary.md")
    if not os.path.exists(summary):
        raise FileNotFoundError(f"{summary} not found: run the workflow (or `python -m alnfail report`) first")
    with open(summary, encoding="utf-8") as fh:
        blocks = parse_markdown(fh.read())

    # split into the introduction and the sections; the "Figures" list only supplies captions
    intro: list[tuple] = []
    sections: list[dict] = []
    captions: dict[str, str] = {}
    for block in blocks:
        if block[0] == "h1":
            continue
        if block[0] == "h2":
            match = re.match(r"(\d+)\.\s*(.*)", block[1])
            sections.append({"title": block[1], "number": int(match.group(1)) if match else None, "blocks": []})
        elif not sections:
            intro.append(block)
        elif sections[-1]["title"].lower() == "figures" and block[0] == "ul":
            for item in block[1]:
                match = re.match(r"`([^`]+)`:\s*(.*)", item)
                if match:
                    captions[match.group(1)] = match.group(2)
        else:
            sections[-1]["blocks"].append(block)
    sections = [s for s in sections if s["title"].lower() != "figures"]

    # every PNG in the folder is shown once: in its section if that section exists, otherwise at the end
    on_disk = sorted(f for f in os.listdir(report_dir) if f.lower().endswith(".png"))
    ordered = [f for f in captions if f in on_disk] + [f for f in on_disk if f not in captions]
    numbers = {s["number"] for s in sections}
    placed: dict[int, list[str]] = {}
    leftover: list[str] = []
    for name in ordered:
        number = _figure_section(name)
        (placed.setdefault(number, []) if number in numbers else leftover).append(name)

    def figures(names: list[str]) -> str:
        return "".join(_figure(report_dir, n, captions.get(n, n)) for n in names)

    heat_max = max([v for s in sections for b in s["blocks"] if b[0] == "table" and is_heat_table(b[1])
                    for v in heat_values(b[2])] + [0.0])

    def render(block: tuple) -> str:
        if block[0] == "table":
            return render_table(block[1], block[2], heat_max)
        if block[0] == "ul":
            return "<ul>" + "".join(f"<li>{inline(item)}</li>" for item in block[1]) + "</ul>"
        if block[1].startswith("**"):  # a bold lead-in names the table that follows
            return f'<p class="subhead">{inline(block[1])}</p>'
        return f"<p>{inline(block[1])}</p>"

    parts, toc = [], []
    for section in sections:
        anchor = _slug(section["title"])
        toc.append(f'<li><a href="#{anchor}">{html.escape(section["title"])}</a></li>')
        content = "".join(render(b) for b in section["blocks"])
        has_ruler = any(b[0] == "table" and any(h.startswith("wrong per 1,000") for h in b[1]) for b in section["blocks"])
        how = HOW_TO_READ.get(section["number"], "")
        how_html = f'<p class="how">{html.escape(how)}{" " + html.escape(RULER_LEGEND) if has_ruler else ""}</p>' if how else ""
        pictures = figures(placed.get(section["number"], []))
        ordered_content = pictures + content if section["number"] in FIGURES_FIRST else content + pictures
        parts.append(f'<section id="{anchor}"><h2>{html.escape(section["title"])}</h2>{how_html}{ordered_content}</section>')
    if leftover:
        toc.append('<li><a href="#more-figures">More figures</a></li>')
        parts.append(f'<section id="more-figures"><h2>More figures</h2>{figures(leftover)}</section>')

    built = datetime.datetime.now().strftime("%d %B %Y, %H:%M").lstrip("0")
    source = html.escape(os.path.normpath(report_dir))
    page = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{html.escape(TITLE)}</title>
<style>{CSS}</style>
</head>
<body>
<div class="page">
<nav aria-label="Sections"><p>On this page</p><ol>{"".join(toc)}</ol></nav>
<main>
<header>
<h1>{html.escape(TITLE)}</h1>
{"".join(render(b) for b in intro)}
<p class="built">Built on {built} from <code>{source}</code>. This page restates <code>summary.md</code>; it adds no numbers of its own.</p>
</header>
{"".join(parts)}
</main>
</div>
</body>
</html>
"""
    out_path = out_path or os.path.join(report_dir, "report.html")
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(page)
    return out_path


def main(argv: list[str] | None = None) -> None:
    argv = sys.argv[1:] if argv is None else argv
    if not 1 <= len(argv) <= 2:
        sys.exit("usage: python -m alnfail.html <report folder> [output.html]")
    print(f"written: {build_html(*argv)}")


if __name__ == "__main__":
    main()
