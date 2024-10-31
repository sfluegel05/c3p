from rdkit import Chem
from rdkit.Chem import AllChem

def is_phenolate_anion(smiles: str):
    """
    Determines if a molecule is a phenolate anion (deprotonated phenol).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a phenolate anion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find all phenol oxygen atoms with negative charge
    phenolate_oxygens = []
    for atom in mol.GetAtoms():
        if (atom.GetSymbol() == 'O' and 
            atom.GetFormalCharge() == -1 and
            len(atom.GetNeighbors()) == 1 and
            atom.GetNeighbors()[0].GetIsAromatic()):
            phenolate_oxygens.append(atom)
            
    if not phenolate_oxygens:
        return False, "No phenolate anion groups found"

    # Check that each phenolate oxygen is attached to an aromatic ring
    for oxygen in phenolate_oxygens:
        carbon = oxygen.GetNeighbors()[0]
        if not carbon.IsInRing() or not carbon.GetIsAromatic():
            return False, "Phenolate oxygen not attached to aromatic ring"
            
    num_phenolates = len(phenolate_oxygens)
    if num_phenolates == 1:
        return True, "Molecule contains 1 phenolate anion group"
    else:
        return True, f"Molecule contains {num_phenolates} phenolate anion groups"
# Pr=None
# Recall=None