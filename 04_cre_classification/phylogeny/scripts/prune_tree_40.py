#!/usr/bin/env python3

from pathlib import Path
from Bio import Phylo

ROOT = Path(__file__).resolve().parent.parent

DATA_DIR = ROOT / "data"
RESULTS_DIR = ROOT / "results"

TREE_IN = DATA_DIR / "301Fly_HOG_UCLDtree.nw"
SPECIES_FILE = DATA_DIR / "species_40_tree_names.txt"

TREE_OUT = RESULTS_DIR / "301Fly_HOG_UCLDtree_40species.nw"
ORDER_OUT = RESULTS_DIR / "species_order_40_tree_names.txt"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

wanted = {
    x.strip()
    for x in SPECIES_FILE.read_text().splitlines()
    if x.strip()
}

tree = Phylo.read(TREE_IN, "newick")

tips_before = {
    terminal.name
    for terminal in tree.get_terminals()
}

missing = sorted(wanted - tips_before)

print(f"Tips in original tree: {len(tips_before)}")
print(f"Requested species:     {len(wanted)}")

if missing:
    print("\nERROR: requested species missing from tree:")
    for name in missing:
        print("  ", name)
    raise SystemExit(1)

for terminal in list(tree.get_terminals()):
    if terminal.name not in wanted:
        tree.prune(terminal)

tips_after = tree.get_terminals()

print(f"Tips after pruning:     {len(tips_after)}")

if len(tips_after) != 40:
    raise SystemExit(
        f"ERROR: expected 40 tips, got {len(tips_after)}"
    )

order = [terminal.name for terminal in tips_after]

Phylo.write(tree, TREE_OUT, "newick")

ORDER_OUT.write_text(
    "\n".join(order) + "\n"
)

print(f"\nWrote tree:  {TREE_OUT}")
print(f"Wrote order: {ORDER_OUT}")

print("\nPhylogenetic tip order:")
for i, name in enumerate(order, start=1):
    print(f"{i:2d}\t{name}")
