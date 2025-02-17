"""
Classifies: CHEBI:16038 phosphatidylethanolamine
"""
"""
Classifies: Phosphatidylethanolamine
Definition: A class of glycerophospholipids in which a phosphatidyl group is esterified 
to the hydroxy group of ethanolamine.
Heuristics used:
  - Check for a phosphate (P) that is bonded to an oxygen chain pattern "OCCN" (ethanolamine head group).
  - Verify the presence of at least two ester bonds (for the fatty acyl chains).
  - Confirm a phosphorus atom is present.
  - Optional molecular weight threshold (e.g. >500 Da) to filter out very small molecules.
Note: This is a heuristic approach.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylethanolamine(smiles: str):
    """
    Determines if a molecule is a phosphatidylethanolamine based on its SMILES string.
    A phosphatidylethanolamine is characterized by a phosphatidyl group (with a phosphate)
    bearing an ethanolamine substituent (identified by the -OCCN- pattern) and two fatty
    acyl chains (ester bonds) esterified to a glycerol backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as phosphatidylethanolamine, False otherwise.
        str: Explanation for the classification result.
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Check for ethanolamine head group via a phosphate-ethanolamine fragment.
    # We look for a P atom directly connected to an oxygen that leads to a 2-carbon chain ending in N.
    # This SMARTS "P-OCCN" is a simplification and should pick up typical phosphoethanolamine head groups,
    # including cases where the nitrogen might be substituted.
    headgroup_pattern = Chem.MolFromSmarts("P-OCCN")
    if not mol.HasSubstructMatch(headgroup_pattern):
        return False, "Missing phosphoethanolamine head group pattern (P-OCCN not found)"
    
    # Check for at least two ester bonds (fatty acyl chains).
    # This pattern detects a carbonyl (C=O) directly bonded to an oxygen.
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Insufficient ester groups for fatty acid chains (found {len(ester_matches)}, require at least 2)"
    
    # Check that a phosphorus atom is present
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not p_atoms:
        return False, "No phosphorus atom found in the structure"
    
    # Check molecular weight; typical phospholipids are > 500 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for a phosphatidylethanolamine (>500 Da expected)"
    
    return True, "Structure contains a phosphate bonded to an ethanolamine head group and at least 2 ester groups typical of phosphatidylethanolamine"

# Example usage (you can remove or comment out these lines when integrating as a module)
if __name__ == "__main__":
    # A test SMILES string from one of the provided examples:
    test_smiles = "P(OCC(OC(=O)CCCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCCCCCC/C=C\\CCCCCCCC)(OCCN(C)C)(O)=O"
    result, reason = is_phosphatidylethanolamine(test_smiles)
    print(result, reason)