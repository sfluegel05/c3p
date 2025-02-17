"""
Classifies: CHEBI:24654 hydroxy fatty acid
"""
"""
Classifies: Hydroxy Fatty Acid
Definition: Any fatty acid carrying one or more hydroxy substituents.
A molecule qualifies if it has a terminal carboxylic acid group – that is, the acid group’s carbon is attached to exactly one other carbon (typical of the fatty acid endgroup) – as well as at least one additional hydroxyl group outside the acid.
We also require that the molecule is aliphatic (a high fraction of its heavy atoms are carbons) so that peptides or polyaromatic compounds are ruled out.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acid based on its SMILES string.
    
    Criteria:
      - Must contain a carboxylic acid group in the form C(=O)[OX2H]
      - The acid group must be terminal: the carbonyl carbon must be attached to exactly one carbon (typical for fatty acids).
      - Must contain at least one –OH group that is not part of the acid group.
      - The molecule should be mainly aliphatic. We require that the fraction of heavy (non-H) atoms that are carbon is at least 0.70.
      - Optionally, we enforce a minimal number of carbon atoms (e.g. 4) to avoid very small acids.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is classified as a hydroxy fatty acid, False otherwise.
        str: Reason for the classification decision.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens to help detect hydroxyl groups.
    mol = Chem.AddHs(mol)
    
    # Define SMARTS for the carboxylic acid group.
    acid_pattern = Chem.MolFromSmarts("C(=O)[OX2H]")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if not acid_matches:
        return False, "No carboxylic acid group found (not a fatty acid)"
    
    # We now filter these acid matches to select only those where the acid group is terminal.
    valid_acid_indices = []
    for match in acid_matches:
        # In the acid SMARTS "C(=O)[OX2H]", match[0] is the carboxyl carbon and match[1] is the acid hydroxyl oxygen.
        acid_carbon = mol.GetAtomWithIdx(match[0])
        # Count neighboring carbons (ignoring oxygens or other atoms).
        carbon_neighbors = [nbr for nbr in acid_carbon.GetNeighbors() if nbr.GetAtomicNum() == 6]
        if len(carbon_neighbors) == 1:
            # Accept as a valid fatty acid acid group.
            valid_acid_indices.append(match)
    
    if not valid_acid_indices:
        return False, "No terminal carboxylic acid group found (acid carbon not attached to exactly one carbon)"
    
    # For our classification, pick the first valid acid match.
    acid_match = valid_acid_indices[0]
    acid_oh_idx = acid_match[1]  # the oxygen atom that is part of the acid group
    
    # Check that the molecule has a minimal number of carbons.
    carbons = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6]
    if len(carbons) < 4:
        return False, "Too few carbons to be considered a fatty acid"
    
    # Check aliphatic content: the fraction of heavy (non-hydrogen) atoms that are carbon.
    heavy_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() != 1]
    n_carbon = sum(1 for atom in heavy_atoms if atom.GetAtomicNum() == 6)
    if heavy_atoms and (n_carbon / len(heavy_atoms)) < 0.70:
        return False, "Molecule is not sufficiently aliphatic to be a fatty acid"
     
    # Define a SMARTS pattern for hydroxyl groups.
    # [OX2H] matches an oxygen with two connections and one hydrogen (typical for –OH).
    oh_pattern = Chem.MolFromSmarts("[OX2H]")
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    
    # Filter out the hydroxyl that is part of the acid group.
    non_acid_oh = []
    for match in oh_matches:
        # match is a tuple with one element: the index of an oxygen atom
        if match[0] != acid_oh_idx:
            non_acid_oh.append(match[0])
    
    if not non_acid_oh:
        return False, "No hydroxy substituent found outside the carboxylic acid group"
    
    # If all conditions are met, classify as a hydroxy fatty acid.
    return True, "Molecule contains a terminal carboxylic acid group, sufficient aliphatic chain, and hydroxy substituent(s), classifying it as a hydroxy fatty acid"

# Example usage:
if __name__ == '__main__':
    # Test a few example SMILES strings.
    test_smiles = [
        "C(C(O)=O)C[C@@H](/C=C/C=C\\C/C=C\\C/C=C\\C=C\\[C@H](C/C=C\\CC)O)O",  # resolvin D6 (should be True)
        "OC(=CC)C(O)=O",  # 2-hydroxy-2-butenoic acid (should be True)
        "OC(=O)CCC(CCCC)CC",  # 4-Ethyloctanoic acid (should be False)
    ]
    for s in test_smiles:
        result, reason = is_hydroxy_fatty_acid(s)
        print(f"SMILES: {s}")
        print(f"Result: {result}, Reason: {reason}\n")