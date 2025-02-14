"""
Classifies: CHEBI:51689 enone
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_enone(smiles: str):
    """
    Determines if a molecule is an enone based on its SMILES string.
    An enone is an alpha,beta-unsaturated ketone with the general formula R1R2C=CR3-C(=O)R4 (R4 != H)

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an enone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for the enone motif
    enone_pattern = Chem.MolFromSmarts("[C]=[C]-[C]=O")
    
    # Find the enone substructure
    matches = mol.GetSubstructMatches(enone_pattern)
    if not matches:
      return False, "No alpha,beta-unsaturated ketone (C=C-C=O) motif found."
    
    # Check R4 is not H
    for match in matches:
        carbonyl_carbon_index = match[2]
        carbonyl_carbon = mol.GetAtomWithIdx(carbonyl_carbon_index)
        has_hydrogen = False
        for neighbor in carbonyl_carbon.GetNeighbors():
            if neighbor.GetSymbol() == 'H':
                has_hydrogen=True
                break
        if not has_hydrogen:
           return True, "Contains an alpha,beta-unsaturated ketone (C=C-C=O) motif with R4 != H"
        
    return False, "C=O carbon has a hydrogen substituent"