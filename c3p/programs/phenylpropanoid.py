"""
Classifies: CHEBI:26004 phenylpropanoid
"""
from rdkit import Chem

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    Phenylpropanoids have an aromatic structure based on a phenylpropane skeleton 
    and may include flavonoids, coumarins, etc.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phenylpropanoid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Refined patterns for key aromatic systems in phenylpropanoids
    aromatic_system_patterns = [
        Chem.MolFromSmarts("c1ccccc1C=C"),   # Basic phenylpropane structure
        Chem.MolFromSmarts("c1ccccc1CO"),    # Phenyl with hydroxymethyl linkage (common in chromen-based compounds)
        Chem.MolFromSmarts("c1ccc2c(c1)cc(=O)oc2"),  # Core flavonoid/coumarin
    ]

    # Check for aromatic backbone
    backbone_match = any(mol.HasSubstructMatch(pattern) for pattern in aromatic_system_patterns)
    if not backbone_match:
        return False, "No suitable phenylpropane skeleton found"

    # Specific functional group and subclass patterns
    subclass_patterns = [
        Chem.MolFromSmarts("c1c(oc2ccccc2)c(cc1)"),  # Flavone/core chromone scaffold
        Chem.MolFromSmarts("c1c(oc2ccc(cc2)o1)"),    # Coumarin or 2H-chromen-2-one flag
        Chem.MolFromSmarts("Oc1ccc(cc1)C=C"),        # Phenolic compound with propenyl linkage
    ]

    # Check if any subclass pattern matches
    subclass_match = any(mol.HasSubstructMatch(pattern) for pattern in subclass_patterns)
    
    if not subclass_match:
        return False, "No specific phenylpropanoid subclass structures found"

    return True, "Structure contains a phenylpropane skeleton with functional groups typical of phenylpropanoids."