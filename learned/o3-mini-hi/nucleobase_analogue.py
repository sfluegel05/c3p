"""
Classifies: CHEBI:67142 nucleobase analogue
"""
"""
Classifies: Nucleobase analogue
A nucleobase analogue is defined as a molecule that can substitute for a normal nucleobase in nucleic acids.
This program uses several heuristics:
  1. The molecule must be a valid RDKit molecule.
  2. Its molecular weight should be between 80 and 350 Da.
  3. It must contain at least 2 nitrogen atoms.
  4. It must be small (≤35 heavy atoms).
  5. It must be nearly rigid (≤3 rotatable bonds) so that only molecules with limited conformational flexibility pass.
  6. It should have 1 or 2 rings (nucleobases are typically mono- or bicyclic).
  7. It must contain at least one appropriate aromatic heterocyclic ring. In our case we scan for an aromatic ring of size 5 or 6 that has at least 2 nitrogen atoms.
  8. The ratio of carbon atoms to nitrogen atoms must be ≤4.0.
  9. Finally, since natural nucleobases consist almost entirely of an aromatic ring, we require that the largest qualifying aromatic ring represents at least 75% of the total heavy atoms.
  
These filters are intended to eliminate false positives (molecules that just “look like” nucleobase analogues due to a small aromatic bit) and false negatives (molecules with too many flexible substituents).
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    
    Heuristics used:
      1. Valid molecule that RDKit can handle.
      2. Molecular weight between 80 and 350 Da.
      3. Contains ≥2 nitrogen atoms.
      4. Contains ≤35 heavy atoms.
      5. Contains ≤3 rotatable bonds (i.e. is rigid).
      6. Contains 1 or 2 rings (typical for nucleobases).
      7. Contains at least one aromatic heterocyclic ring (of size 5 or 6) that has at least 2 nitrogen atoms.
      8. The ratio of carbon atoms to nitrogen atoms is ≤4.0.
      9. The largest aromatic nucleobase ring must account for at least 75% of the molecule’s heavy atoms.
      
    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule qualifies as a nucleobase analogue, False otherwise.
        str: A reason explaining the classification decision.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization failed: {str(e)}"
        
    # Check molecular weight.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 80 or mol_wt > 350:
        return False, f"Molecular weight {mol_wt:.1f} Da is outside expected range (80–350 Da)"
    
    # Count nitrogen atoms (require at least 2).
    n_nitrogen = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_nitrogen < 2:
        return False, "Too few nitrogen atoms; a nucleobase analogue should contain at least 2 nitrogen atoms"
    
    # Count heavy atoms (should be ≤35).
    heavy_atom_count = mol.GetNumHeavyAtoms()
    if heavy_atom_count > 35:
        return False, f"Too many heavy atoms ({heavy_atom_count}); nucleobase analogues are typically small (≤35 heavy atoms)"
    
    # Count rotatable bonds (should be ≤3 for rigidity).
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable > 3:
        return False, f"Too many rotatable bonds ({n_rotatable}); nucleobase analogues are expected to be rigid (≤3 rotatable bonds)"
    
    # Retrieve ring information.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    n_rings = len(atom_rings)
    if n_rings < 1:
        return False, "No ring structures detected; nucleobase analogues are heterocyclic"
    if n_rings > 2:
        return False, f"Too many rings ({n_rings}); natural nucleobases are usually mono- or bicyclic"
    
    # Look for at least one appropriate aromatic heterocyclic ring.
    # We require the ring to be of size 5 or 6 (typical for pyrimidine or the 6-membered part of purines)
    # and it must contain at least 2 nitrogen atoms.
    candidate_ring_sizes = []
    for ring in atom_rings:
        ring_size = len(ring)
        if ring_size not in (5, 6):
            continue
        # Check that every atom in the ring is aromatic.
        if not all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            continue
        # Count the number of nitrogen atoms in the ring.
        ring_nitrogens = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() == 7)
        if ring_nitrogens >= 2:
            candidate_ring_sizes.append(ring_size)
    
    if not candidate_ring_sizes:
        return False, "No appropriate aromatic heterocyclic ring (size 5 or 6 with at least 2 nitrogen atoms) found"
    
    # For the purpose of the next check, consider the largest candidate ring.
    largest_ring_size = max(candidate_ring_sizes)
    
    # Determine what fraction of the molecule is made up of the aromatic core.
    aromatic_core_fraction = largest_ring_size / heavy_atom_count
    if aromatic_core_fraction < 0.75:
        return False, ("The molecule has too many non‐aromatic substituents (largest aromatic ring contains "
                       f"{largest_ring_size} of {heavy_atom_count} heavy atoms, fraction {aromatic_core_fraction:.2f}); "
                       "nucleobase analogues are expected to be predominantly an aromatic core.")
    
    # Check the ratio of carbon to nitrogen atoms.
    n_carbon = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if n_nitrogen > 0 and (n_carbon / n_nitrogen) > 4.0:
        return False, (f"High carbon-to-nitrogen ratio ({n_carbon}/{n_nitrogen} = {n_carbon/n_nitrogen:.2f}); "
                       "nucleobase analogues are generally nitrogen rich")
    
    return True, ("Molecule exhibits features of a nucleobase analogue: appropriate molecular weight, nitrogen content, "
                  "small size, rigidity, limited ring system with an aromatic heterocyclic core (5–6 membered ring with ≥2 nitrogens) "
                  "that comprises most of the structure, and a balanced C/N ratio.")

# Example testing section (uncomment below to test several cases)
# if __name__ == "__main__":
#     test_smiles = [
#         "C=12C(=NC=NC1NCC(=O)N)NC=N2",      # N(6)-carbamoylmethyladenine (True)
#         "Nc1nc(=O)[nH]cc1CO",                # 5-(hydroxymethyl)cytosine (True)
#         "Nc1ncc2nc[nH]c2n1",                # 2-aminopurine (True)
#         "CCC[NH2+][C@H]1CCc2[nH+]c(N)sc2C1",  # pramipexole(2+) (False)
#         "[C@H](CNC1=C(C(NC(N1)=O)=O)/N=C/C(=O)C)([C@@H](CCO)O)O",  # 6-[(1,4-dideoxy-D-ribityl)amino]-... (False)
#         "[H]C(=O)Nc1cncnc1"                 # 5-formamidopyrimidine (False)
#     ]
#     for sm in test_smiles:
#         res, reason = is_nucleobase_analogue(sm)
#         print(f"SMILES: {sm}\nResult: {res}\nReason: {reason}\n")