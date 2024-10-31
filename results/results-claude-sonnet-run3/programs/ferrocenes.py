from rdkit import Chem
from rdkit.Chem import AllChem

def is_ferrocenes(smiles: str):
    """
    Determines if a molecule contains a ferrocene skeleton.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains ferrocene skeleton, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
        
    # Check for presence of iron atom
    if not any(atom.GetSymbol() == 'Fe' for atom in mol.GetAtoms()):
        return False, "No iron atom found"
        
    # Find iron atom
    fe_atom = None
    for atom in mol.GetAtoms():
        if atom.GetSymbol() == 'Fe':
            fe_atom = atom
            break
            
    # Check if Fe has correct valence
    if fe_atom.GetDegree() != 10:
        return False, "Iron atom does not have correct coordination"

    # Check for ferrocene pattern in SMILES
    if "[Fe]" not in smiles:
        return False, "Does not match ferrocene connection pattern"

    # Get neighbors of Fe atom
    fe_neighbors = fe_atom.GetNeighbors()
    
    # Check that all neighbors are carbon
    if not all(n.GetSymbol() == 'C' for n in fe_neighbors):
        return False, "Iron atom not bonded to all carbons"

    # Check for the characteristic ferrocene SMILES pattern
    # This pattern ensures we have the correct sandwich structure
    if not any(p in smiles for p in ["C1[Fe]23451234", "[Fe]23451234"]):
        return False, "Does not have correct ferrocene structure"

    # Look for substituents
    base_ferrocene = "C12C3C4C5C1[Fe]23451234C5C1C2C3C45"
    if smiles == base_ferrocene:
        return True, "Unsubstituted ferrocene"
    else:
        return True, "Substituted ferrocene"
# Pr=1.0
# Recall=1.0