"""
Classifies: CHEBI:17478 aldehyde
"""
#!/usr/bin/env python
"""
Classifies: Aldehyde – a compound containing the functional group RC(=O)H
where a carbonyl carbon is bonded to a hydrogen and to one R group.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_aldehyde(smiles: str):
    """
    Determines if a molecule is an aldehyde based on its SMILES string.
    An aldehyde is defined as a compound containing the functional group –CHO 
    (i.e. a carbonyl group with one hydrogen attached).

    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule contains at least one aldehyde group, False otherwise.
        str: Reason for classification.
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens so we can detect the hydrogen bound to the carbonyl carbon.
    mol = Chem.AddHs(mol)

    # Define the aldehyde SMARTS pattern.
    # [CX3H1](=O) represents a trigonal carbon (sp2 hybridized) with exactly one hydrogen,
    # double bonded to an oxygen atom.
    aldehyde_pattern = Chem.MolFromSmarts("[CX3H1](=O)")
    if aldehyde_pattern is None:
        return False, "Failed to create aldehyde pattern"
    
    # Search for aldehyde substructures within the molecule.
    matches = mol.GetSubstructMatches(aldehyde_pattern)
    if not matches:
        return False, "No aldehyde group ([CX3H1](=O)) found in the molecule"
    
    # If multiple aldehyde groups are found, we still classify the compound as aldehyde.
    count = len(matches)
    reason = f"Found {count} aldehyde group{'s' if count > 1 else ''} matching the pattern [CX3H1](=O)"
    
    return True, reason

# (Optional) Testing code; remove or comment out when using in a larger project.
if __name__ == "__main__":
    test_smiles = [
        "O=CC(CCC=C(C)C)C",  # 5-Heptenal, 2,6-dimethyl-
        "Oc1c(C=O)ccc2ccccc12",  # 1-hydroxy-2-naphthaldehyde
        "CCCCCCCCCCCCCCCCCC=O",  # octadecanal
        "O=CCCCCCCCCCCCCCCCCCCCCCCC",  # tetracosanal
        "CCCC/C=C/C=O",  # (E)-hept-2-enal
        "CC(C)CC=O",  # 3-methylbutanal
        "C(C#C)=O",  # prop-2-ynal
    ]
    
    for sm in test_smiles:
        result, explanation = is_aldehyde(sm)
        print(f"SMILES: {sm}\nResult: {result}\nExplanation: {explanation}\n")