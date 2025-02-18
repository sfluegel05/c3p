"""
Classifies: CHEBI:59644 oxo fatty acid
"""
"""
Classifies: Oxo fatty acid
Definition: Any fatty acid containing at least one aldehydic or ketonic group
            in addition to the carboxylic acid group.
Additional criteria (heuristic):
  - The molecule must contain only typical fatty acid elements (C, O, and H).
  - It must be acyclic (no rings).
  - There must be exactly one carboxylic acid group (terminal –C(=O)[O]H).
  - The molecule must contain at least a minimal number of carbons (>=5) for meaningful fatty acid structure.
  - In addition, an extra oxo function (ketone or aldehyde) must be present
    (its carbonyl is not part of the acid).
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_oxo_fatty_acid(smiles: str):
    """
    Determines if a molecule is an oxo fatty acid based on its SMILES string.
    
    The molecule must contain:
      - A single carboxylic acid group (as in fatty acids) that should appear at a terminus.
      - At least one additional oxo function (ketone or aldehyde) on a different carbon.
    In addition, the molecule must be acyclic and contain only C, O and H atoms.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if the molecule is an oxo fatty acid, False otherwise.
        str: A reason for the classification.
    """
    # Parse input SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string."
    
    # Criterion: The molecule should contain only C, H, and O.
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in (1, 6, 8):
            return False, "Molecule contains elements other than C, H, and O."
    
    # Criterion: Must be acyclic (linear fatty acids are expected to have no rings)
    if mol.GetRingInfo().NumRings() > 0:
        return False, "Molecule contains rings and does not appear to be a linear fatty acid."
    
    # Criterion: The molecule should have a minimal number of carbons (at least 5)
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_count < 5:
        return False, f"Molecule contains too few carbons ({carbon_count}) to be a fatty acid."
    
    # Criterion: Check for exactly one carboxylic acid group.
    # SMARTS pattern for carboxylic acid: a carbonyl carbon (sp2) with an -OH group.
    acid_smarts = "[CX3](=O)[OX2H]"
    acid_pattern = Chem.MolFromSmarts(acid_smarts)
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    if len(acid_matches) == 0:
        return False, "Molecule does not contain a carboxylic acid group."
    elif len(acid_matches) > 1:
        return False, "Molecule contains more than one carboxylic acid group, not a typical fatty acid."
    
    # Record the carbonyl atom index from the acid group so that we do not count it as an extra oxo.
    acid_carbonyl_indices = {match[0] for match in acid_matches}  # assuming first atom is carbonyl
    
    # Criterion: Look for an extra carbonyl group (either ketone or aldehyde) outside the acid.
    extra_oxo_found = False
    
    # Check for a ketone using SMARTS: ketone should be a carbonyl carbon (C=O) with two carbon neighbors.
    ketone_smarts = "[#6][CX3](=O)[#6]"
    ketone_pattern = Chem.MolFromSmarts(ketone_smarts)
    for match in mol.GetSubstructMatches(ketone_pattern):
        # match is a tuple (R, C, R') where the middle is the carbonyl carbon.
        carbonyl_index = match[1]
        if carbonyl_index in acid_carbonyl_indices:
            continue  # skip the carboxylic acid carbonyl
        extra_oxo_found = True
        break
    
    # If no ketone was found, look for an aldehyde.
    if not extra_oxo_found:
        # SMARTS for aldehyde: a carbonyl carbon bonded to one carbon and one hydrogen.
        aldehyde_smarts = "[#6][CX3H1](=O)"
        aldehyde_pattern = Chem.MolFromSmarts(aldehyde_smarts)
        for match in mol.GetSubstructMatches(aldehyde_pattern):
            carbonyl_index = match[1]
            if carbonyl_index in acid_carbonyl_indices:
                continue
            extra_oxo_found = True
            break

    if not extra_oxo_found:
        return False, "Molecule does not contain an additional aldehydic or ketonic group outside the carboxylic acid."
    
    return True, "Molecule is an oxo fatty acid: it is acyclic, contains only C, H and O, has a single carboxylic acid group and an extra oxo function."

# (Optionally, one could include a __main__ block for testing.)
if __name__ == '__main__':
    test_smiles = [
        "CC(C)C(=O)CC(O)=O",        # 4-methyl-3-oxopentanoic acid (true)
        "O=C(CCCC)CCC(O)=O",         # 4-keto-n-caprylic acid (expected false based on extra criteria)
        "N[C@@H](CCCCNCC(=O)[C@@H](O)[C@H](O)[C@H](O)CO)C(O)=O",  # fructosyllysine (false, contains N)
        "CCCC(=O)C(O)=O",           # 2-oxopentanoic acid (true)
    ]
    for smi in test_smiles:
        result, reason = is_oxo_fatty_acid(smi)
        print(f"SMILES: {smi}\nResult: {result}\nReason: {reason}\n")