"""
Classifies: CHEBI:32877 primary amine
"""
from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.
    A primary amine is characterized by the presence of the -NH2 group attached to a hydrocarbyl group.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a primary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Use a smarter and more inclusive SMARTS pattern to recognize primary amines, allowing for aromatic rings and complex scaffolds
    # [NX3H2]: nitrogen with exactly two hydrogens
    # [CX4][!#8;!#7]: carbon attached, ensuring that it’s not directly oxygen or another nitrogen forming complications
    primary_amine_patterns = [
        Chem.MolFromSmarts("[NX3H2][CX4]"),
        Chem.MolFromSmarts("[NX3H2][CX3](=O)"), # example for carbonyl linkages
        Chem.MolFromSmarts("[NX3H2][c]")  # aromatic connections
    ]
    
    for pattern in primary_amine_patterns:
        matches = mol.GetSubstructMatches(pattern)
        
        if matches:
            # Verify that each matched nitrogen is attached to a hydrocarbyl group
            for match in matches:
                nitrogen_atom = mol.GetAtomWithIdx(match[0])
                # Ensure that this nitrogen is primarily bonded to acceptable hydrocarbyl group environments
                associated_hydrocarbons = any(
                    [b.GetOtherAtom(nitrogen_atom).GetAtomicNum() == 6 for b in nitrogen_atom.GetBonds()]
                )
                
                if associated_hydrocarbons:
                    return True, f"Primary amine group (-NH2) identified with pattern {pattern}"

    return False, "Primary amine group not found"