"""
Classifies: CHEBI:16038 phosphatidylethanolamine
"""
"""
Classifies: Phosphatidylethanolamine
Definition: A class of glycerophospholipids in which a phosphatidyl group is esterified 
to the hydroxy group of ethanolamine.
Heuristics used (revised):
  - Look for a phosphate group that is bonded to an oxygen–ethyl chain ending in a non–quaternary amine.
    The SMARTS "P(=O)(O)(OCC[N;H;!$([N]C(=O)])])" aims to restrict matches to true phosphoethanolamine head groups.
  - Require at least two ester bonds (each is defined as a carbonyl “C(=O)” bonded to an oxygen).
  - Confirm that a phosphorus atom is present.
  - Use a molecular weight threshold of 400 Da (lower than before) to allow for short-chain examples.
Note: This is a heuristic approach and may miss edge cases.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylethanolamine(smiles: str):
    """
    Determines if a molecule is a phosphatidylethanolamine based on its SMILES string.
    A phosphatidylethanolamine (PE) is characterized by:
      1. A phosphate group that is esterified with an ethanolamine head group.
         (We require a substructure of a phosphate connected via an oxygen to an ethyl chain ending in an amine,
         making sure the amine is not acylated as in phosphatidylserine and is not quaternary as in phosphatidylcholine.)
      2. Two fatty acyl chains (indicated by at least two ester bonds).
      3. Presence of a phosphorus atom.
      4. A typical molecular weight (we set a lower bound of 400 Da).
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as a phosphatidylethanolamine, False otherwise.
        str: Explanation for the classification result.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # 1. Check for the phosphoethanolamine head group.
    # The SMARTS here looks for a phosphate (P with =O and O substituents) bonded to an oxygen
    # that is connected to a two-carbon chain terminating in an amine.
    # The substructure "OCC[N;H;!$([N]C(=O)])]" forces the nitrogen to have at least one hydrogen
    # and to NOT be connected to a carbonyl (which would be expected in a serine head group).
    headgroup_smarts = "P(=O)(O)(OCC[N;H;!$([N]C(=O)])])"
    headgroup_pattern = Chem.MolFromSmarts(headgroup_smarts)
    if not mol.HasSubstructMatch(headgroup_pattern):
        return False, "Missing or incorrect phosphoethanolamine head group pattern"
    
    # 2. Check for at least two ester bonds.
    # An ester bond is identified through the pattern: carbonyl (C=O) bonded to an oxygen.
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 2:
        return False, f"Insufficient ester groups for fatty acyl chains (found {len(ester_matches)}, require at least 2)"
    
    # 3. Verify the presence of phosphorus.
    p_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15]
    if not p_atoms:
        return False, "No phosphorus atom found in the structure"
    
    # 4. Check that the molecular weight is high enough.
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400:
        return False, f"Molecular weight too low for a phosphatidylethanolamine ({mol_wt:.1f} Da; expected >= 400 Da)"
    
    return True, "Structure contains a phosphate bonded to an ethanolamine head group and at least 2 ester bonds typical of phosphatidylethanolamine"

# Example usage (can be removed or commented out when integrating as a module)
if __name__ == "__main__":
    # Example SMILES string from one of the provided examples:
    test_smiles = "P(OCC(OC(=O)CCCCCC/C=C\\C/C=C\\C/C=C\\C/C=C\\CC)COC(=O)CCCCCCC/C=C\\CCCCCCCC)(OCCN(C)C)(O)=O"
    result, reason = is_phosphatidylethanolamine(test_smiles)
    print(result, reason)