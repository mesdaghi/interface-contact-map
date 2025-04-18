# interface-contact-map
This script extracts and visualizes the **Cα–Cα contact map** at the **interface between two chains** in a PDB file using [ConKit](https://github.com/rigdenlab/conkit) and Biopython.


##  What it does

- Parses a multi-chain PDB file.
- Selects two specified chains.
- Calculates residue–residue contacts (based on Cα distance cutoff).
- Constructs a contact map of the inter-chain interface.
- Saves:
  - A `.map` file in CASP contact format.
  - A `.png` image of the contact map.

##  Requirements

- Python 3.x
- [Biopython](https://biopython.org/)
- [ConKit](https://github.com/rigdenlab/conkit)
- Matplotlib (for visualization)

Install them via pip:

```bash
pip install biopython conkit matplotlib
