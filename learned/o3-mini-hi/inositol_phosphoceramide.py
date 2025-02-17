"""
Classifies: CHEBI:60245 inositol phosphoceramide
"""
"""
Classifies: Inositol Phosphoceramide
Definition: A phosphosphingolipid in which an inositol residue and the ceramide moiety are linked via a phosphodiester bridge.
The ceramide moiety contains substituents (typically an amide bond and long aliphatic chains).
In this implementation:
  - A more relaxed SMARTS pattern is used for the inositol ring: a cyclohexane with one hydroxyl on each carbon.
  - We verify that at least one phosphorus atom connects via an oxygen to an inositol carbon.
  - We require an amide bond indicative of the ceramide portion.
  - Finally, we check the overall molecular size.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_inositol_phosphoceramide(smiles: str):
    """
    Determines if a molecule is an inositol phosphoceramide based on its SMILES string.
    
    The method checks:
      1) A proper inositol ring substructure: a cyclohexane bearing a hydroxyl group on each carbon.
         (We use a relaxed SMARTS pattern to account for variable chirality annotation.)
      2) A phosphodiester linkage: the presence of a phosphate (P) wherein at least one oxygen
         connects to an atom that is part of the inositol ring.
      3) A ceramide-like moiety: evidenced by the presence of an amide bond (C(=O)N).
      4) Sufficient molecular size (by weight and carbon count) to reflect the long aliphatic chains.
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an inositol phosphoceramide, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Use a relaxed SMARTS pattern for inositol.
    # Myo-inositol is a cyclohexane with a hydroxyl group on each carbon.
    # The pattern below ignores chirality and matches a hexagon where each carbon is substituted with an OH.
    inositol_smarts = "OC1C(O)C(O)C(O)C(O)C1O"
    inositol_query = Chem.MolFromSmarts(inositol_smarts)
    inositol_matches = mol.GetSubstructMatches(inositol_query)
    if not inositol_matches:
        return False, "Missing inositol ring substructure"
        
    # Collect indices of all atoms that are part of any inositol match.
    inositol_atom_indices = set()
    for match in inositol_matches:
        inositol_atom_indices.update(match)
    
    # Check for phosphodiester linkage:
    # Look for a phosphorus (P, atomic number 15) atom that is bonded to oxygen(s),
    # where at least one such oxygen is connected (through another bond) to an atom in the inositol ring.
    phosphate_found = False
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 15:
            continue
        for neighbor in atom.GetNeighbors():
            # Only consider oxygen neighbors (atomic number 8)
            if neighbor.GetAtomicNum() != 8:
                continue
            # For this oxygen, check its other connected atoms
            other_neighbors = [nbr for nbr in neighbor.GetNeighbors() if nbr.GetIdx() != atom.GetIdx()]
            for o_nbr in other_neighbors:
                if o_nbr.GetIdx() in inositol_atom_indices:
                    phosphate_found = True
                    break
            if phosphate_found:
                break
        if phosphate_found:
            break
    if not phosphate_found:
        return False, "Phosphodiester linkage (connecting phosphate to inositol) not detected"
    
    # Check for an amide bond (indicative of the ceramide portion).
    amide_smarts = "C(=O)N"
    amide_query = Chem.MolFromSmarts(amide_smarts)
    if not mol.HasSubstructMatch(amide_query):
        return False, "Missing ceramide amide bond"
    
    # Check overall molecular size: weight and carbon count
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for inositol phosphoceramide"
    
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 25:
        return False, f"Too few carbons ({carbon_count}) to be an inositol phosphoceramide"
    
    # If all conditions are met, classify the molecule accordingly.
    return True, "Molecule contains an inositol ring, a phosphodiester bridge to phosphate, and a ceramide-like amide bond with sufficient molecular size"

# Example usage (for testing purposes):
if __name__ == "__main__":
    test_smiles = "CCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)[C@H](O)CCCCCCCCCCCCCC"
    result, reason = is_inositol_phosphoceramide(test_smiles)
    print("Classification:", result)
    print("Reason:", reason)