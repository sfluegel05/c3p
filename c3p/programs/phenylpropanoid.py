"""
Classifies: CHEBI:26004 phenylpropanoid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phenylpropanoid(smiles: str):
    """
    Determines if a molecule is a phenylpropanoid based on its SMILES string.
    Phenylpropanoids contain an aromatic ring connected to a propane chain or similar structures like flavonoids.
    
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

    # Look for aromatic ring
    aromatic_pattern = Chem.MolFromSmarts("c1ccccc1")  # Simple benzene ring pattern
    if not mol.HasSubstructMatch(aromatic_pattern):
        return False, "No aromatic benzene ring found"

    # Look for aliphatic three-carbon chain attached to the aromatic ring
    # The pattern may look for a chain like -C-C-C- attached to an aromatic carbon
    propane_chain_pattern = Chem.MolFromSmarts("c-C-C-C")
    if not mol.HasSubstructMatch(propane_chain_pattern):
        return False, "No phenylpropane backbone found"

    # Check for common functional groups in phenylpropanoids
    # Hydroxy (-OH), Methoxy (-OCH3), Carbonyl (=O) are common.
    functional_group_patterns = [
        Chem.MolFromSmarts("c-[OH]"),       # Phenolic OH group
        Chem.MolFromSmarts("c-CO"),         # Carbonyl
        Chem.MolFromSmarts("c-COC")         # Methoxy
    ]
    function_group_count = sum(1 for pattern in functional_group_patterns 
                               if mol.HasSubstructMatch(pattern))

    if function_group_count < 1:
        return False, "No common phenylpropanoid functional groups found"

    return True, "Structure contains aromatic ring with a phenylpropane backbone and functional groups common to phenylpropanoids."