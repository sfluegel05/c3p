"""
Classifies: CHEBI:60245 inositol phosphoceramide
"""
"""
Classifies: Inositol Phosphoceramide
Definition: A phosphosphingolipid in which an inositol residue and the ceramide moiety are linked via a phosphodiester bridge.
The ceramide moiety contains an amide and long aliphatic chains.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_inositol_phosphoceramide(smiles: str):
    """
    Determines if a molecule is an inositol phosphoceramide based on its SMILES string.
    
    The method searches for:
      1) An inositol ring: a 6-membered ring with hydroxyl substitutions.
      2) A phosphodiester group: using a phosphate pattern (OP(=O)(O)O).
      3) A ceramide-like moiety: at least one amide bond as found in ceramides.
    Additionally, the molecule must be of sufficient molecular size (molecular weight and carbon count).
    
    Args:
        smiles (str): SMILES string of the molecule.
    
    Returns:
        bool: True if the molecule is classified as an inositol phosphoceramide, False otherwise.
        str: Explanation of the classification decision.
    """
    
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # SMARTS for a generic inositol ring.
    # This pattern uses a 6-membered ring with one hydroxyl on each carbon.
    inositol_smarts = "OC1C(O)C(O)C(O)C(O)C1O"
    inositol_query = Chem.MolFromSmarts(inositol_smarts)
    if not mol.HasSubstructMatch(inositol_query):
        return False, "Missing inositol ring substructure"
    
    # SMARTS for a phosphate group attached via oxygen.
    # This pattern finds an -O-P(=O)(O)O group.
    phosphate_smarts = "[O]P(=O)(O)O"
    phosphate_query = Chem.MolFromSmarts(phosphate_smarts)
    if not mol.HasSubstructMatch(phosphate_query):
        return False, "Missing phosphate (phosphodiester) group"
    
    # SMARTS for an amide bond (indicative of the ceramide acyl linkage).
    amide_smarts = "C(=O)N"
    amide_query = Chem.MolFromSmarts(amide_smarts)
    if not mol.HasSubstructMatch(amide_query):
        return False, "Missing ceramide amide bond"
    
    # Check that the molecule is large enough. Inositol phosphoceramides are lipids and quite heavy.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a phosphoceramide"
    
    # Count carbon atoms: ceramides contain long aliphatic chains.
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 25:
        return False, f"Too few carbons ({carbon_count}) to be an inositol phosphoceramide"
    
    # If all the checks pass, we classify the molecule as an inositol phosphoceramide.
    return True, "Molecule contains an inositol ring, a phosphodiester linkage, and a ceramide-like amide bond with sufficient molecular size"
    
# Example usage (for testing purposes):
if __name__ == "__main__":
    test_smiles = "CCCCCCCCCCCCCCCCCCCCCCC(=O)N[C@@H](COP(O)(=O)O[C@H]1[C@H](O)[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O)[C@H](O)CCCCCCCCCCCCCC"
    result, reason = is_inositol_phosphoceramide(test_smiles)
    print("Classification:", result)
    print("Reason:", reason)