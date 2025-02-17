"""
Classifies: CHEBI:16038 phosphatidylethanolamine
"""
"""
Classifies: Phosphatidylethanolamine
Definition: A class of glycerophospholipids in which a phosphatidyl group is esterified 
to the hydroxy group of ethanolamine.
Heuristics used (revised):
  - Look for a phosphate head group directly attached via an oxygen to an ethanolamine-like fragment.
    We use the SMARTS "P(=O)(O)(OCC[NX3])" which should catch both the unmodified and N–substituted
    ethanolamine head groups.
  - Require at least two ester bonds (each detected as [CX3](=O)[OX2]).
  - Confirm that a phosphorus atom is present.
  - Enforce a molecular weight threshold of 400 Da.
Note: This heuristic approach may not capture all edge cases.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylethanolamine(smiles: str):
    """
    Determines if a molecule is a phosphatidylethanolamine based on its SMILES string.
    A phosphatidylethanolamine (PE) is characterized by:
      1. A phosphate group esterified to an ethanolamine-like head group.
      2. At least two ester bonds (indicating two fatty acyl chains).
      3. The presence of a phosphorus atom.
      4. A molecular weight typically high enough (set here as ≥ 400 Da).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is classified as a phosphatidylethanolamine, False otherwise.
        str: Explanation for the classification result.
    """
    # Parse the input SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for the phosphoethanolamine head group.
    # We use a simplified but precise SMARTS pattern:
    # "P(=O)(O)(OCC[NX3])" looks for:
    # - A phosphorus (P) double bonded to an oxygen (=O)
    # - Two oxygen substituents, one of which is directly connected to "CC[NX3]"
    #   where [NX3] is a trivalent nitrogen (this covers both NH2 and N-alkylated amines).
    headgroup_smarts = "P(=O)(O)(OCC[NX3])"
    headgroup_pattern = Chem.MolFromSmarts(headgroup_smarts)
    if headgroup_pattern is None:
        return False, "Error parsing headgroup SMARTS pattern"
    if not mol.HasSubstructMatch(headgroup_pattern):
        return False, "Missing or incorrect phosphoethanolamine head group pattern"
    
    # 2. Check for at least two ester bonds.
    # An ester bond here is defined as an acyl group: a carbonyl (C=O) bonded to an oxygen.
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    if ester_pattern is None:
        return False, "Error parsing ester SMARTS pattern"
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Insufficient ester groups for fatty acyl chains (found {len(ester_matches)}, require at least 2)"
    
    # 3. Verify the presence of at least one phosphorus atom.
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not p_atoms:
        return False, "No phosphorus atom found in the structure"
    
    # 4. Check molecular weight – common phosphatidylethanolamines are heavier molecules.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Molecular weight too low for a phosphatidylethanolamine ({mol_wt:.1f} Da; expected >= 400 Da)"
    
    return True, "Structure contains a phosphoethanolamine head group, at least 2 ester bonds, phosphorus, and acceptable molecular weight"

# For testing purposes (can be removed or commented out when used as a module)
if __name__ == "__main__":
    # Example SMILES string from one of the provided examples:
    test_smiles = "P(OCC(OC(=O)CCCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCCCCCC/C=C\\CCCCCCCC)(OCCN(C)C)(O)=O"
    result, reason = is_phosphatidylethanolamine(test_smiles)
    print(result, reason)