from rdkit import Chem
from rdkit.Chem import AllChem

def is_monohydroxyflavone(smiles: str):
    """
    Determines if a molecule is a monohydroxyflavone (flavone with a single hydroxy substituent).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a monohydroxyflavone, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for flavone core structure 
    # Pattern captures the essential flavone skeleton including the ketone and fused rings
    flavone_pattern = Chem.MolFromSmarts('O=C1CC(c2ccccc2)Oc2ccccc12')
    if not mol.HasSubstructMatch(flavone_pattern):
        return False, "Not a flavone structure"

    # Count number of hydroxy groups directly attached to aromatic carbons
    oh_pattern = Chem.MolFromSmarts('Oc1:c')
    oh_matches = mol.GetSubstructMatches(oh_pattern)
    
    # Count number of phenolic hydroxy groups
    if len(oh_matches) == 0:
        return False, "No hydroxy groups found"
    elif len(oh_matches) > 1:
        return False, f"Contains {len(oh_matches)} hydroxy groups, more than one"
    
    # Check if OH is attached to flavone core
    flavone_atoms = set(mol.GetSubstructMatch(flavone_pattern))
    oh_atom = mol.GetAtomWithIdx(oh_matches[0][0])
    oh_neighbor = oh_atom.GetNeighbors()[0]
    
    if oh_neighbor.GetIdx() not in flavone_atoms:
        return False, "Hydroxy group not attached to flavone core"
        
    # Check for non-phenolic OH groups
    other_oh = Chem.MolFromSmarts('[OH]')
    total_oh = len(mol.GetSubstructMatches(other_oh))
    if total_oh > len(oh_matches):
        return False, "Contains additional non-phenolic OH groups"
        
    return True, "Valid monohydroxyflavone"
# Pr=None
# Recall=0.0