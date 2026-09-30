#!/usr/bin/env python3
"""Update the self-test result table in README.md from fingerprint.txt."""

import os
import re
import sys
import unicodedata
from pathlib import Path

ROOT = Path(__file__).resolve().parent
FINGERPRINT = ROOT / "fingerprint.txt"
README = ROOT / "README.md"

HEADING = "### 4.2 自测结果"
HEADERS = ["k", "n", "用时", "峰值内存"]
ALIGNS = ["left", "right", "right", "right"]


def display_width(text):
    return sum(2 if unicodedata.east_asian_width(ch) in ("W", "F") else 1 for ch in text)


def pad(text, width, align):
    gap = max(0, width - display_width(text))
    if align == "right":
        return " " * gap + text
    return text + " " * gap


def parse_fingerprint(path):
    entries = []
    current = None
    for raw in path.read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        match = re.match(r"k=(\d+),\s*n=(\d+)", line)
        if match:
            current = {"k": int(match.group(1)), "n": int(match.group(2))}
            entries.append(current)
            continue
        if current is None:
            continue
        match = re.match(r"Time:\s*([0-9]+(?:\.[0-9]+)?)s", line)
        if match:
            current["time"] = float(match.group(1))
            continue
        match = re.match(r"Memory:\s*([0-9]+(?:\.[0-9]+)?)\s*MB", line)
        if match:
            current["memory"] = float(match.group(1))
    complete = [entry for entry in entries if "time" in entry and "memory" in entry]
    if not complete:
        raise ValueError(f"no complete benchmark entries found in {path.name}")
    return complete


def build_table(entries):
    rows = [
        [str(entry["k"]), f"{entry['n']:,}", f"{entry['time']:.2f} s", f"{entry['memory']:.2f} MB"]
        for entry in entries
    ]
    widths = [
        max(3, display_width(HEADERS[i]), *(display_width(row[i]) for row in rows))
        for i in range(len(HEADERS))
    ]

    header = "| " + " | ".join(
        pad(HEADERS[i], widths[i], ALIGNS[i]) for i in range(len(HEADERS))
    ) + " |"
    separator = "| " + " | ".join(
        (":" + "-" * widths[i]) if ALIGNS[i] == "left" else ("-" * widths[i] + ":")
        for i in range(len(HEADERS))
    ) + " |"
    body = [
        "| " + " | ".join(pad(row[i], widths[i], ALIGNS[i]) for i in range(len(HEADERS))) + " |"
        for row in rows
    ]
    return [header, separator, *body]


def replace_table(lines, table):
    try:
        start = next(i for i, line in enumerate(lines) if line.strip() == HEADING)
    except StopIteration:
        raise ValueError(f"heading {HEADING!r} not found in README.md")

    header = None
    for i in range(start + 1, len(lines)):
        if lines[i].lstrip().startswith("|") and "峰值内存" in lines[i]:
            header = i
            break
    if header is None:
        raise ValueError("result table not found after the heading")

    end = header
    while end < len(lines) and lines[end].lstrip().startswith("|"):
        end += 1
    return lines[:header] + table + lines[end:]


def update_peak_rss(text, value):
    pattern = re.compile(r"(最大测试规模的峰值 RSS 为：\s*```text\s*\n)([0-9.]+)( MB\s*\n```)")
    updated, count = pattern.subn(
        lambda match: match.group(1) + f"{value:.2f}" + match.group(3), text
    )
    if count == 0:
        raise ValueError("peak RSS anchor not found in README.md")
    return updated


def main():
    entries = parse_fingerprint(FINGERPRINT)
    original = README.read_bytes().decode("utf-8")
    newline = "\r\n" if "\r\n" in original else "\n"

    lines = replace_table(original.splitlines(), build_table(entries))
    updated = newline.join(lines)
    if original.endswith(("\n", "\r")):
        updated += newline
    updated = update_peak_rss(updated, max(entry["memory"] for entry in entries))

    if updated == original:
        print("README.md already up to date")
        return 0

    tmp = README.with_name(README.name + ".tmp")
    tmp.write_bytes(updated.encode("utf-8"))
    os.replace(tmp, README)
    print(f"Updated README.md with {len(entries)} benchmark result(s)")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as exc:
        print(f"warning: {exc}", file=sys.stderr)
        sys.exit(1)
