"""
Classifies: CHEBI:67142 nucleobase analogue
"""
"""
Classifies: Nucleobase analogue
A nucleobase analogue is defined as a molecule that can substitute for a normal nucleobase in nucleic acids.
This improved program uses several heuristics:
 - The molecule must be a valid RDKit molecule.
 - Its molecular weight should be within roughly 80–350 Da.
 - It should contain at least 2 nitrogen atoms.
 - It should be small (≤35 heavy atoms).
 - It should be nearly rigid (≤3 rotatable bonds).
 - The molecule must have at most 3 rings (usually 1 or 2 in natural nucleobases).
 - There must be at least one aromatic, heterocyclic ring (of size 5–7 members with at least one non-carbon atom).
These thresholds are a compromise in order to reduce false positives (large, polyaromatic molecules) and 
avoid misclassifying flexible nucleosides (with large ribityl substituents) as nucleobase analogues.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    
    The algorithm uses these heuristics:
      1. The SMILES string must be valid and successfully sanitized.
      2. Molecular weight should be within a narrow range (80–350 Da).
      3. The molecule should contain at least 2 nitrogen atoms.
      4. The molecule should be small in terms of heavy atoms (≤35 heavy atoms).
      5. The molecule should be relatively rigid (≤3 rotatable bonds).
      6. The molecule should contain at most 3 rings overall.
      7. There must be at least one aromatic heterocyclic ring (of size 5–7 and containing ≥1 heteroatom).
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a nucleobase analogue, False otherwise.
        str: A reason explaining the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    try:
        # Ensure the molecule is sanitized
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization failed: {str(e)}"
    
    # Check molecular weight (80-350 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 80 or mol_wt > 350:
        return False, f"Molecular weight {mol_wt:.1f} Da is outside expected range (80–350 Da)"
    
    # Count nitrogen atoms: require at least 2 (nucleobases naturally have 2 or more).
    n_nitrogen = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_nitrogen < 2:
        return False, "Too few nitrogen atoms; a nucleobase analogue should contain at least 2 nitrogen atoms"
    
    # Check heavy atom count to ensure the molecule is small.
    heavy_atom_count = mol.GetNumHeavyAtoms()
    if heavy_atom_count > 35:
        return False, f"Too many heavy atoms ({heavy_atom_count}); nucleobase analogues are typically small (≤35 heavy atoms)"
    
    # Check rotatable bonds: nucleobases are rigid. Allow up to 3.
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable > 3:
        return False, f"Too many rotatable bonds ({n_rotatable}); nucleobase analogues are expected to be rigid (≤3 rotatable bonds)"
    
    # Retrieve ring information: natural nucleobases typically have 1–2 rings.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    n_rings = len(atom_rings)
    if n_rings == 0:
        return False, "No ring structures found; nucleobase analogues are heterocyclic"
    if n_rings > 3:
        return False, f"Too many rings ({n_rings}); nucleobase analogues typically have 1 or 2 aromatic rings"
    
    # Look for at least one appropriate aromatic heterocyclic ring.
    heterocycle_found = False
    for ring in atom_rings:
        ring_size = len(ring)
        if ring_size in [5, 6, 7]:
            # Check that every atom in the ring is aromatic.
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                # Count heteroatoms (non-carbon atoms) in the ring.
                hetero_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6)
                if hetero_count >= 1:
                    heterocycle_found = True
                    break

    if not heterocycle_found:
        return False, "No appropriate aromatic heterocyclic ring (size 5–7 with at least one heteroatom) found"
    
    return True, ("Molecule exhibits features of a nucleobase analogue: appropriate molecular weight, nitrogen count, small size, rigidity, "
                  "limited ring system, and the presence of an aromatic heterocyclic ring (5–7 members).")

# Example test calls (uncomment to test)
# if __name__ == "__main__":
#     test_smiles = [
#         "Nc1nc(=O)[nH]cc1CO",    # 5-(hydroxymethyl)cytosine (should be True)
#         "Nc1ncc2nc[nH]c2n1",      # 2-aminopurine (should be True)
#         "Oc1ccc2c3ccnc4ccc(=O)n(c2c1)c34"  # 9-hydroxycanthin-6-one (should be False)
#     ]
#     for sm in test_smiles:
#         res, reason = is_nucleobase_analogue(sm)
#         print(f"SMILES: {sm}\nResult: {res}\nReason: {reason}\n")