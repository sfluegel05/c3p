"""
Classifies: CHEBI:16219 cucurbitacin
"""
"""
Classifies: Cucurbitacin (tetracyclic triterpenoids from the cucurbitane family)
Based on heuristics: 
  – must be a valid molecule,
  – should contain at least 4 rings (the tetracyclic core),
  – should have at least 25 carbon atoms,
  – and should show at least one carbonyl (C=O) group.
Note: Because cucurbitacins can feature additional sugar moieties or oxygenations, 
      this classifier uses a few rough criteria and may miss or mis‐classify edge cases.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_cucurbitacin(smiles: str):
    """
    Determines if a molecule is a cucurbitacin based on its SMILES string.
    Cucurbitacins are a class of tetracyclic triterpenoids (derived from cucurbitane).
    This heuristic check includes:
      - Verification of a valid molecule.
      - At least 4 rings present.
      - At least 25 carbon atoms (consistent with a triterpenoid skeleton).
      - Presence of one or more carbonyl groups (C=O), as many cucurbitacins contain ketones.
      
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule meets the heuristic criteria for a cucurbitacin; False otherwise.
        str: A reason message explaining the classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Count rings using the ring Info from RDKit
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    if num_rings < 4:
        return False, f"Detected only {num_rings} rings; expected at least 4 for a tetracyclic core"
    
    # Count the number of carbon atoms
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 25:
        return False, f"Too few carbon atoms ({carbon_count}); expected at least 25 as in a triterpene"
    
    # Search for carbonyl groups (C=O) using a SMARTS pattern.
    # Note: Many cucurbitacins contain one or more ketone groups.
    carbonyl_pattern = Chem.MolFromSmarts("[CX3]=O")
    carbonyl_matches = mol.GetSubstructMatches(carbonyl_pattern)
    if len(carbonyl_matches) < 1:
        return False, "No carbonyl group (C=O) found; cucurbitacins typically have ketone functionality"

    # (Optional) Check molecular weight. Cucurbitacins aglycones are often in the 400–600 Da range.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) compared with typical cucurbitacins"
    
    return True, ("The molecule has at least four rings, a sufficient number of carbon atoms "
                  "and at least one carbonyl group, consistent with a cucurbitacin (tetracyclic triterpenoid) scaffold.")

# Example usage:
if __name__ == "__main__":
    # Example: cucurbitacin D (one of the provided examples)
    example_smiles = "[H][C@@]12C[C@H](O)C(=O)C(C)(C)C1=CC[C@@]1([H])[C@]3(C)C[C@@H](O)[C@]([H])([C@@](C)(O)C(=O)\\C=C\\C(C)(C)O)[C@@]3(C)CC(=O)[C@@]21C"
    result, reason = is_cucurbitacin(example_smiles)
    print("Cucurbitacin check:", result)
    print("Reason:", reason)