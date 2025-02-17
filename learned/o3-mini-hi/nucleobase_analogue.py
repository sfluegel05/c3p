"""
Classifies: CHEBI:67142 nucleobase analogue
"""
"""
Classifies: Nucleobase analogue
A nucleobase analogue is defined as a molecule that can substitute for a normal nucleobase in nucleic acids.
This improved program uses several heuristics:
 - The molecule must be valid.
 - The molecular weight should be within roughly 80–500 Da.
 - It should contain at least 2 nitrogen atoms.
 - It should be small in terms of heavy atoms (<= 30 heavy atoms).
 - It should be nearly rigid (<= 2 rotatable bonds).
 - There must be at least one aromatic, heterocyclic ring of common nucleobase ring sizes (5–7 atoms).
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    
    The algorithm uses these heuristics:
      1. The SMILES string must be valid and successfully sanitized.
      2. Molecular weight should be within a reasonable range (80–500 Da).
      3. The molecule should contain at least 2 nitrogen atoms.
      4. The molecule should have a low number of heavy atoms (<= 30) 
         to ensure it is small.
      5. The molecule should be relatively rigid (<= 2 rotatable bonds).
      6. There must be at least one aromatic heterocyclic ring (size 5–7) containing at least one non-carbon atom.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule qualifies as a nucleobase analogue, False otherwise.
        str: A reason explaining the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    try:
        # Ensure molecule is sanitized
        Chem.SanitizeMol(mol)
    except Exception as e:
        return False, f"Sanitization failed: {str(e)}"
    
    # Check molecular weight (80-500 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 80 or mol_wt > 500:
        return False, f"Molecular weight {mol_wt:.1f} Da is outside expected nucleobase range (80–500 Da)"
    
    # Count nitrogen atoms; natural nucleobases have at least 2.
    n_nitrogen = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_nitrogen < 2:
        return False, "Too few nitrogen atoms; nucleobase analogues should contain at least 2 nitrogen atoms"
    
    # Check heavy atom count: nucleobase analogues are small.
    heavy_atom_count = mol.GetNumHeavyAtoms()
    if heavy_atom_count > 30:
        return False, f"Too many heavy atoms ({heavy_atom_count}); nucleobase analogues are typically small (<=30 heavy atoms)"
    
    # Check rotatable bonds: nucleobases are rigid.
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable > 2:
        return False, f"Too many rotatable bonds ({n_rotatable}); nucleobase analogues are expected to be rigid"
    
    # Retrieve ring information and check for aromatic heterocyclic rings.
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if not atom_rings:
        return False, "No ring structures found; nucleobase analogues are heterocyclic"
    
    heterocycle_found = False
    for ring in atom_rings:
        # Consider rings of common nucleobase sizes (5, 6, or 7 members)
        if len(ring) in [5, 6, 7]:
            # Verify that all atoms in the ring are aromatic.
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                # Count atoms that are not carbon.
                hetero_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6)
                if hetero_count >= 1:
                    heterocycle_found = True
                    break

    if not heterocycle_found:
        return False, "No appropriate aromatic heterocyclic ring (size 5–7 with at least one heteroatom) found"
    
    # Passed all criteria: likely a nucleobase analogue.
    return True, "Molecule exhibits features of a nucleobase analogue: appropriate molecular weight, nitrogen content, small size, rigidity, and aromatic heterocyclic ring(s)."

# Example test calls (you can uncomment these lines to test)
# if __name__ == "__main__":
#     # Example: hypoxanthine
#     result, reason = is_nucleobase_analogue("O=c1[nH]cnc2nc[nH]c12")
#     print("Hypoxanthine:", result, reason)
#     
#     # Example: an invalid SMILES
#     result, reason = is_nucleobase_analogue("invalid_smiles")
#     print("Invalid SMILES:", result, reason)