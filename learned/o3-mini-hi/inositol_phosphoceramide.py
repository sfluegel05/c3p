"""
Classifies: CHEBI:60245 inositol phosphoceramide
"""
"""
Classifies: Inositol Phosphoceramide
Definition: A phosphosphingolipid in which an inositol residue and the ceramide moiety are linked via a phosphodiester bridge.
The ceramide moiety contains an amide and long aliphatic chains.
In this implementation:
  - A tight inositol SMARTS pattern is used (a chiral six‐membered ring with six hydroxyl groups).
  - We verify that at least one phosphorus atom (“P”) connects via oxygen to an inositol carbon.
  - We require an amide bond indicative of the ceramide portion.
  - Finally, we check overall molecular size.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_inositol_phosphoceramide(smiles: str):
    """
    Determines if a molecule is an inositol phosphoceramide based on its SMILES string.
    
    The method checks:
      1) A proper inositol ring: a chiral six‐membered ring with one –OH on each carbon.
      2) A phosphodiester linkage: a phosphate (P) where at least one oxygen (O)
         is directly attached (via its second bond) to an inositol carbon.
      3) A ceramide-like moiety: indicated by an amide bond.
      4) Sufficient molecular size (molecular weight and number of carbons).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an inositol phosphoceramide, False otherwise.
        str: Explanation of the classification decision.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define a more specific SMARTS for inositol:
    # This pattern implies a six-membered ring (with chiral centers) with hydroxyl substitutions.
    inositol_smarts = "[C@H]1([OH])C([OH])C([OH])C([OH])C([OH])C1([OH])"
    inositol_query = Chem.MolFromSmarts(inositol_smarts)
    inositol_matches = mol.GetSubstructMatches(inositol_query)
    if not inositol_matches:
        return False, "Missing inositol ring substructure"
    
    # Gather all atom indices found in any inositol match for later verification.
    inositol_atom_indices = set()
    for match in inositol_matches:
        inositol_atom_indices.update(match)
    
    # Check for a phosphodiester linkage.
    # We look for at least one phosphorus (atomic number 15) that is bonded to oxygens,
    # and at least one of these oxygens is linked (via its other bond) to an atom in the inositol ring.
    phosphate_found = False
    # Iterate over all atoms; find P atoms.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 15:
            continue  # Not phosphorus
        # For phosphorus, check its neighbors (typically oxygen atoms) 
        # and see if any oxygen is connected via its other bond to an inositol atom.
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() != 8:
                continue  # We want oxygen
            # Get atoms bonded to this oxygen besides the phosphorus.
            oxy_neighbors = [nbr for nbr in neighbor.GetNeighbors() if nbr.GetIdx() != atom.GetIdx()]
            for o_nbr in oxy_neighbors:
                if o_nbr.GetIdx() in inositol_atom_indices:
                    phosphate_found = True
                    break
            if phosphate_found:
                break
        if phosphate_found:
            break
    if not phosphate_found:
        return False, "Phosphodiester linkage (connecting phosphate to inositol) not detected"
    
    # Check for an amide bond (indicative of the ceramide acyl linkage).
    amide_smarts = "C(=O)N"
    amide_query = Chem.MolFromSmarts(amide_smarts)
    if not mol.HasSubstructMatch(amide_query):
        return False, "Missing ceramide amide bond"
    
    # Check overall molecular size: weight and carbon count (long aliphatic chains present)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for inositol phosphoceramide"
    
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 25:
        return False, f"Too few carbons ({carbon_count}) to be an inositol phosphoceramide"
    
    # If all checks pass, classify as an inositol phosphoceramide.
    return True, "Molecule contains a correctly linked inositol ring, a phosphodiester bridge, and a ceramide-like amide bond with sufficient molecular size"

# Example usage (for testing purposes):
if __name__ == "__main__":
    test_smiles = "CCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)[C@H](O)CCCCCCCCCCCCCC"
    result, reason = is_inositol_phosphoceramide(test_smiles)
    print("Classification:", result)
    print("Reason:", reason)