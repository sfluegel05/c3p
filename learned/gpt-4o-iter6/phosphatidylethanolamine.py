"""
Classifies: CHEBI:16038 phosphatidylethanolamine
"""
from rdkit import Chem

def is_phosphatidylethanolamine(smiles: str):
    """
    Determines if a molecule is a phosphatidylethanolamine based on its SMILES string.
    A phosphatidylethanolamine is a glycerophospholipid where a phosphatidyl group is esterified to
    the hydroxy group of ethanolamine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylethanolamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern with carbon atoms marked
    glycerol_pattern = Chem.MolFromSmarts("O[C@H]([O])CO")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found with correct stereochemistry"
        
    # Look for phosphate group (-P(=O)(O)-)
    phosphate_pattern = Chem.MolFromSmarts("P(=O)(O)")
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Phosphate group not found"
    
    # Look for ethanolamine group (-OCCN)
    ethanolamine_pattern = Chem.MolFromSmarts("OCCN")
    if not mol.HasSubstructMatch(ethanolamine_pattern):
        return False, "Ethanolamine group not found"

    return True, "Structure matches phosphatidylethanolamine class"

# Example usage
smiles_example = "P(OC[C@H](OC(=O)CCCCCCCCCCCCCCC)COC(=O)CCCCCCC/C=C\\CCCCCCC)(OCCN)(O)=O"
result, reason = is_phosphatidylethanolamine(smiles_example)
print(result, reason)