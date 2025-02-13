"""
Classifies: CHEBI:23086 chalcones
"""
from rdkit import Chem

def is_chalcones(smiles: str):
    """
    Determines if a molecule is a chalcone based on its SMILES string.
    Chalcones are defined as 1,3-diphenylpropenone or benzylideneacetophenone derivatives.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a chalcone, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Enhanced and more specific SMARTS patterns for the chalcone core ArCH=CHC(=O)Ar
    chalcone_patterns = [
        Chem.MolFromSmarts("c1ccccc1C=CC(=O)c2ccccc2"),  # Basic chalcone pattern
        Chem.MolFromSmarts("c1ccccc1C=CC(=O)Ar"),  # General pattern with Ar as any aromatic
        Chem.MolFromSmarts("c1cccnc1C=CC(=O)c2cccnc2"),  # Nitrogen in aromatic rings
        Chem.MolFromSmarts("C=CC(=O)c1ccc([OH])cc1"),  # Hydroxy substituted aromatic rings
        Chem.MolFromSmarts("c1ccc(C=CC(=O)c2ccccc2)cc1"),  # Another common substitution pattern
        Chem.MolFromSmarts("C=CC(=O)C2=C(OC)C=CC=C2"),  # Consider oxygen substitutions on aromatic rings
        Chem.MolFromSmarts("c1ccc2c(c1)c(C=CC(=O)c3ccc4c(c3)OCO4)cc2"),  # Benzodioxole and related systems
    ]

    for pattern in chalcone_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Molecule matches chalcone core structure ArCH=CHC(=O)Ar"
        
    return False, "Does not match chalcone core structure"