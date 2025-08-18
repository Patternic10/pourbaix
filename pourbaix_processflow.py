#!/usr/bin/env python3
"""
Generate a simple flowchart PNG for the temperature-dependent Pourbaix pipeline.

Output: pourbaix_simple_flowchart.png (in current working directory)

It tries Graphviz first (pip install graphviz; also install system Graphviz),
and falls back to Matplotlib if Graphviz is unavailable.
"""

import os
import sys

OUT_PNG = "pourbaix_simple_flowchart.png"

def try_graphviz(path_png: str) -> bool:
    try:
        from graphviz import Digraph
    except Exception as e:
        return False

    dot = Digraph(comment="Simple Pourbaix Flowchart", format="png")
    dot.attr(rankdir="TB", size="6")

    # Nodes
    dot.node("A", "Start", shape="oval")
    dot.node("B", "Inputs: elements, T, concentration", shape="box")
    dot.node("C", "Fetch Pourbaix entries\n(MPRester)", shape="box")
    dot.node("D", "Is T ≈ 298 K?", shape="diamond")
    dot.node("E", "Use entries as-is", shape="box")
    dot.node("F", "Shift energies for T", shape="box")
    dot.node("G", "Adjust solids ΔG(T)", shape="box")
    dot.node("H", "Adjust ions ΔG(T)", shape="box")
    dot.node("I", "Build Pourbaix diagram", shape="box")
    dot.node("J", "Plot & Save diagram", shape="box")
    dot.node("K", "End", shape="oval")

    # Edges
    dot.edge("A", "B")
    dot.edge("B", "C")
    dot.edge("C", "D")
    dot.edge("D", "E", label="Yes")
    dot.edge("D", "F", label="No")
    dot.edge("F", "G")
    dot.edge("F", "H")
    dot.edge("E", "I")
    dot.edge("G", "I")
    dot.edge("H", "I")
    dot.edge("I", "J")
    dot.edge("J", "K")

    # Render directly to the requested filename (without extension)
    stem, _ = os.path.splitext(os.path.abspath(path_png))
    try:
        dot.render(stem, format="png", cleanup=True)
        return True
    except Exception:
        # If render fails (e.g., system Graphviz not installed), report failure
        return False


def fallback_matplotlib(path_png: str) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

    def add_box(ax, x, y, text, w=4.8, h=0.9):
        rect = FancyBboxPatch((x - w/2, y - h/2), w, h,
                              boxstyle="round,pad=0.02,rounding_size=0.03",
                              fc="white", ec="black", lw=1.2)
        ax.add_patch(rect)
        ax.text(x, y, text, ha="center", va="center", fontsize=10)
        return (x, y, w, h)

    def add_arrow(ax, src, dst, label=None):
        x1, y1, w1, h1 = src
        x2, y2, w2, h2 = dst
        # simple vertical routing
        start = (x1, y1 - h1/2) if y2 < y1 else (x1, y1 + h1/2)
        end   = (x2, y2 + h2/2) if y2 < y1 else (x2, y2 - h2/2)
        arr = FancyArrowPatch(start, end, arrowstyle="-|>", mutation_scale=12, lw=1.1, color="black")
        ax.add_patch(arr)
        if label:
            ax.text((start[0]+end[0])/2, (start[1]+end[1])/2 + 0.2, label, ha="center", va="bottom", fontsize=9)

    fig, ax = plt.subplots(figsize=(8, 10))
    ax.set_axis_off()
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 18)

    # y positions (top to bottom)
    Y = {
        "A": 16.5, "B": 15.0, "C": 13.5, "D": 12.0,
        "E": 10.5, "F": 10.5, "G": 9.0, "H": 9.0,
        "I": 7.0, "J": 5.5, "K": 4.0
    }

    bA = add_box(ax, 5,  Y["A"], "Start", w=3.0)
    bB = add_box(ax, 5,  Y["B"], "Inputs: elements, T, concentration")
    bC = add_box(ax, 5,  Y["C"], "Fetch Pourbaix entries\n(MPRester)")
    bD = add_box(ax, 5,  Y["D"], "Is T ≈ 298 K?")

    # Branches
    bE = add_box(ax, 3,  Y["E"], "Use entries as-is")
    bF = add_box(ax, 7,  Y["F"], "Shift energies for T")
    bG = add_box(ax, 7,  Y["G"], "Adjust solids ΔG(T)")
    bH = add_box(ax, 7,  Y["H"], "Adjust ions ΔG(T)")
    bI = add_box(ax, 5,  Y["I"], "Build Pourbaix diagram")
    bJ = add_box(ax, 5,  Y["J"], "Plot & Save diagram")
    bK = add_box(ax, 5,  Y["K"], "End", w=3.0)

    # Arrows
    add_arrow(ax, bA, bB)
    add_arrow(ax, bB, bC)
    add_arrow(ax, bC, bD)

    # From decision
    # Left branch (Yes → use as-is)
    add_arrow(ax, bD, bE)
    ax.text(4.0, (Y["D"]+Y["E"])/2 + 0.2, "Yes", fontsize=9, ha="center")

    # Right branch (No → shift, then solids & ions)
    add_arrow(ax, bD, bF)
    ax.text(6.0, (Y["D"]+Y["F"])/2 + 0.2, "No", fontsize=9, ha="center")
    add_arrow(ax, bF, bG)
    add_arrow(ax, bF, bH)

    # Merge down to build
    add_arrow(ax, bE, bI)
    add_arrow(ax, bG, bI)
    add_arrow(ax, bH, bI)

    add_arrow(ax, bI, bJ)
    add_arrow(ax, bJ, bK)

    fig.tight_layout()
    fig.savefig(path_png, dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    # Try Graphviz first; fall back to Matplotlib if needed
    ok = try_graphviz(OUT_PNG)
    if not ok:
        try:
            fallback_matplotlib(OUT_PNG)
        except Exception as e:
            sys.stderr.write(
                "Failed to generate with Graphviz and Matplotlib.\n"
                "If using Graphviz, please:\n"
                "  pip install graphviz\n"
                "  and install the system package (e.g., apt-get install graphviz / brew install graphviz).\n"
            )
            raise

    print("Saved:", os.path.abspath(OUT_PNG))
