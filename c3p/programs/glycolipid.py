"""
Classifies: CHEBI:33563 glycolipid
"""
"""
Classifies: glycolipid
Definition: Any member of the class of 1,2-di-O-acylglycerols joined at oxygen 3 by a glycosidic linkage to a carbohydrate part (usually a mono-, di- or tri-saccharide). Some substances (e.g. bacterial glycolipids) may have the sugar part acylated by one or more fatty acids and the glycerol part may be absent.
This program uses heuristic substructure matching to decide if a given SMILES looks like a glycolipid.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycolipid(smiles: str):
    """
    Determines if a molecule is a glycolipid based on its SMILES string.
    The algorithm uses a couple of heuristic substructure matches:
      1. It checks for the presence of a carbohydrate moiety (using a hexopyranose-like ring pattern).
      2. It checks for an acyl group attachment (either as an ester or an amide) which is common to glycolipids.
      3. It checks that the molecular weight is moderately high (e.g. >500 Da) so that the molecule contains long chain(s).
      
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is classified as a glycolipid, False otherwise.
        str: Reason for classification.
    """
    # Parse SMILES with RDKit.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a SMARTS pattern for a common carbohydrate ring (hexopyranose pattern).
    # This is a heuristic pattern that will match many sugars.
    sugar_pattern = Chem.MolFromSmarts("C1OC(O)C(O)C(O)C1O")
    if sugar_pattern is None:
        return False, "Error in sugar SMARTS pattern"
    
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) == 0:
        return False, "No carbohydrate (sugar ring) moiety found"
    
    # Define a SMARTS pattern for an ester group: acyl attached via an oxygen (common in 1,2-di-O-acylglycerols).
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=O)[#6]")
    if ester_pattern is None:
        return False, "Error in ester SMARTS pattern"
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    
    # Also define a pattern for an amide group (common in ceramide-based glycolipids).
    amide_pattern = Chem.MolFromSmarts("NC(=O)[#6]")
    if amide_pattern is None:
        return False, "Error in amide SMARTS pattern"
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    
    # Require at least one acyl linkage: either ester or amide.
    if (len(ester_matches) + len(amide_matches)) < 1:
        return False, "No acyl linkage (ester or amide) found that could anchor a fatty acid chain"
    
    # Check molecular weight: glycolipids often have a moderately high weight (>500 Da).
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for a glycolipid"
    
    # Count rotatable bonds and long chain fragments as extra support (optional heuristic)
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 3:
        return False, "Not enough rotatable bonds to support presence of one or more long fatty acyl chains"
    
    # If everything passes, we classify the molecule as a glycolipid.
    return True, "Molecule contains a carbohydrate moiety and acyl linkage(s) consistent with a glycolipid structure"

# For testing (uncomment the following lines and run the script):
# if __name__ == "__main__":
#     test_smiles = "OC[C@H]1O[C@@H](O[C@@H]2[C@@H](O)[C@H](O)[C@@H](CO)O[C@H]2OCCCCCCCC\\C=C/CCCCCCCC(O)=O"
#     result, reason = is_glycolipid(test_smiles)
#     print("Is glycolipid?", result)
#     print("Reason:", reason)