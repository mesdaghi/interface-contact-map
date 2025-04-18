from conkit.io import pdb
from conkit.core.contactmap import ContactMap
from conkit.core.sequence import Sequence
from Bio.PDB import PDBParser
from conkit.core.contactmap import ContactMap
from conkit.core.sequence import Sequence
from conkit.core.contact import Contact
from conkit.io.casp import CaspParser
from conkit.plot import ContactMapFigure

parser = PDBParser(QUIET=True)
structure = parser.get_structure('protein', 'oca2_wo_x10_homo.pdb')

# Create a chain-wise residue list
chain_residues = {}
for model in structure:
    for chain in model:
        chain_id = chain.id
        chain_residues[chain_id] = []
        for residue in chain:
            if 'CA' in residue:
                chain_residues[chain_id].append(residue)

# Choose two chains
chainA_id, chainB_id = 'B', 'C'
residues_A = chain_residues[chainA_id]
residues_B = chain_residues[chainB_id]

# Create contact map and sequence
contact_map = ContactMap("interface_map")
sequence_str = "A" * (len(residues_A) + len(residues_B))
contact_map.sequence = Sequence("interface_seq", sequence_str)

# Distance cutoff for Cα–Cα contacts
cutoff = 8.0

# Add contacts between chain A and B
for res_a in residues_A:
    for res_b in residues_B:
        dist = res_a['CA'] - res_b['CA']
        if dist <= cutoff:
            resnum_a = res_a.get_id()[1]
            resnum_b = res_b.get_id()[1]
            contact = Contact(resnum_a, resnum_b, 1.0)
            contact_map.add(contact)

# Sort contacts by residue IDs
contact_map.sort("id")

parser = CaspParser()
with open("interface_contacts.map", "w") as f:
    parser.write(f, contact_map)

fig = ContactMapFigure(contact_map)
fig.draw()
fig.savefig("interface_contact_map_oca2_wo_x10_homo.png")
