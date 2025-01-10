"""
Classifies: CHEBI:18303 phosphatidyl-L-serine
"""
from rdkit import Chem

def is_phosphatidyl_L_serine(smiles: str):
    """
    Classifies a molecule as a phosphatidyl-L-serine based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is a phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """
    
    # Convert the SMILES string to an RDKit molecule object.
    mol = Chem.MolFromSmiles(smiles)

    # Return false with a reason if the molecule object cannot be created.
    if mol is None:
        return False, "Invalid SMILES string"

    # Define the SMARTS pattern for the phosphatidyl group linked to serine.
    phosphatidyl_serine_pattern = Chem.MolFromSmarts(
        "O[P](=O)(OCC[C@H](N)C(=O)O)OC[C@H](COC(=O)[C,C](C)C)COC(=O)[C,C](C)C"
    )

    # Check if the molecule matches the phosphatidyl-L-serine core structure.
    if not mol.HasSubstructMatch(phosphatidyl_serine_pattern):
        return False, "No phosphatidyl-L-serine structure found"

    return True, "Molecule meets the criteria for phosphatidyl-L-serine"