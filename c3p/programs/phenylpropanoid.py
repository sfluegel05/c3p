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

    # Key aromatic systems common in many phenylpropanoids
    aromatic_system_patterns = [
        Chem.MolFromSmarts("c1ccccc1"),  # Phenyl
        Chem.MolFromSmarts("c1ccccc1C=C"),  # Phenyl with propene linkage
        Chem.MolFromSmarts("c1ccc2c(c1)cccc2"),  # Naphthalene-like
        Chem.MolFromSmarts("c1c2c(ccc1)c(=O)oc2"),  # Flavonoid core
        Chem.MolFromSmarts("c1c2c(ccc1)oc(=O)c2"),  # Coumarin
        Chem.MolFromSmarts("c1cc(c2cocc2)ccc1"),  # Coumarin
        Chem.MolFromSmarts("c1ccc(cc1)-C=C-C(=O)"),  # Cinnamic acid-like
        Chem.MolFromSmarts("c1c(oc2ccccc2)c(cc1)"),  # Chromone
    ]
    
    # Check for aromatic systems
    if not any(mol.HasSubstructMatch(pattern) for pattern in aromatic_system_patterns):
        return False, "No suitable aromatic system found"

    # Common functional groups in phenylpropanoids
    functional_group_patterns = [
        Chem.MolFromSmarts("[OX2H]"),    # Hydroxyl group (phenolic OH)
        Chem.MolFromSmarts("C(=O)C"),    # Carbonyl groups
        Chem.MolFromSmarts("C(=O)OC"),   # Ester linkage
        Chem.MolFromSmarts("Oc1ccccc1"), # Phenolic hydroxyl directly on aromatic rings
        Chem.MolFromSmarts("C=C[OX2]"),  # Enolic or diradical oxygens
        Chem.MolFromSmarts("C(C)=C"),    # Alkenyl (from side propane chain)
        Chem.MolFromSmarts("OCc1ccc(cc1)"),  # Ether linkage to aromatic
    ]

    # Check for functional groups
    if not any(mol.HasSubstructMatch(pattern) for pattern in functional_group_patterns):
        return False, "No common phenylpropanoid functional groups found"
    
    return True, "Structure contains an aromatic system and functional groups characteristic of a phenylpropanoid."