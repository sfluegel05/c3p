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

    # Use a SMARTS pattern to recognize primary amines, attaching NH2 to any hydrocarbyl group
    primary_amine_pattern = Chem.MolFromSmarts("[NX3H2][CH0-3]")  
    # NX3H2: nitrogen with exactly two hydrogens
    # CH0-3: carbon that can be part of aliphatic or aromatic hydrocarbyl
    
    matches = mol.GetSubstructMatches(primary_amine_pattern)
    
    if matches:
        # Verify that each matched nitrogen is attached to a hydrocarbyl group
        for match in matches:
            nitrogen_atom = mol.GetAtomWithIdx(match[0])
            carbon_atom = mol.GetAtomWithIdx(match[1])

            # Ensure that this nitrogen is only bonded to hydrogen and carbon atoms
            is_primary_amine = all(
                [b.GetOtherAtom(nitrogen_atom).GetAtomicNum() in [1, 6] for b in nitrogen_atom.GetBonds()]
            )
            
            if is_primary_amine:
                return True, "Contains primary amine group (-NH2) attached to a hydrocarbyl group"
    
    return False, "Primary amine group not found"