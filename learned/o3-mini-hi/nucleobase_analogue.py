"""
Classifies: CHEBI:67142 nucleobase analogue
"""
"""
Classifies: Nucleobase analogue
A nucleobase analogue is defined as a molecule that can substitute for a normal nucleobase in nucleic acids.
This program uses several heuristics:
 - The molecule should be small (molecular weight roughly between 80 and 500 Da).
 - It should contain at least two nitrogen atoms.
 - It should contain at least one aromatic, heterocyclic ring of common size (5–7 atoms).
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleobase_analogue(smiles: str):
    """
    Determines if a molecule is a nucleobase analogue based on its SMILES string.
    
    The algorithm uses these heuristics:
      1. The SMILES string must be valid.
      2. Molecular weight should be within a reasonable range (80–500 Da).
      3. The molecule should contain at least two nitrogen atoms.
      4. There must be at least one aromatic ring of common nucleobase ring sizes (5–7 atoms)
         that is heterocyclic (contains at least one non-carbon atom).
    
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
    
    # Check the molecular weight; nucleobases usually fall in a low-molecular weight range.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 80 or mol_wt > 500:
        return False, f"Molecular weight {mol_wt:.1f} Da is outside expected nucleobase range (80–500 Da)"
    
    # Count nitrogen atoms (atomic number 7); natural nucleobases have at least 2.
    n_nitrogen = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
    if n_nitrogen < 2:
        return False, "Too few nitrogen atoms; nucleobase analogues should contain at least 2 nitrogen atoms"
    
    # Retrieve ring information and check for aromatic heterocyclic rings
    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    if not atom_rings:
        return False, "No ring structures found; nucleobase analogues are heterocyclic"
    
    heterocycle_found = False
    for ring in atom_rings:
        # Consider rings common in nucleobases (sizes 5, 6, or 7)
        if len(ring) in [5, 6, 7]:
            # Verify that all atoms in the ring are aromatic
            if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
                # Count atoms that are not carbon (i.e. heteroatoms)
                hetero_count = sum(1 for idx in ring if mol.GetAtomWithIdx(idx).GetAtomicNum() != 6)
                if hetero_count >= 1:
                    heterocycle_found = True
                    break

    if not heterocycle_found:
        return False, "No aromatic heterocyclic ring (of size 5–7) with heteroatoms found"

    # Passed all criteria: likely a nucleobase analogue.
    return True, "Molecule exhibits features of a nucleobase analogue: appropriate molecular weight, sufficient nitrogen content, and at least one aromatic heterocyclic ring."

# Example test calls (you can uncomment these lines to test)
# if __name__ == "__main__":
#     # Example: hypoxanthine
#     result, reason = is_nucleobase_analogue("O=c1[nH]cnc2nc[nH]c12")
#     print("Hypoxanthine:", result, reason)
#     
#     # Example: an invalid SMILES
#     result, reason = is_nucleobase_analogue("invalid_smiles")
#     print("Invalid SMILES:", result, reason)