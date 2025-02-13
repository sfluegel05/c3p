"""
Classifies: CHEBI:24654 hydroxy fatty acid
"""
"""
Classifies: CHEBI:36804 hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a hydroxy fatty acid based on its SMILES string.
    A hydroxy fatty acid is defined as a fatty acid with one or more hydroxy groups
    attached to the aliphatic carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for carboxylic acid group
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"

    # Look for aliphatic carbon chain
    chain_pattern = Chem.MolFromSmarts("[CH3;!$(C(-!C))]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    if len(chain_matches) == 0:
        return False, "No aliphatic carbon chain found"

    # Look for hydroxy group(s) attached to the carbon chain
    hydroxy_pattern = Chem.MolFromSmarts("[$([#8&H1]);!$(*=[#8])]")
    hydroxy_matches = mol.GetSubstructMatches(hydroxy_pattern)
    if len(hydroxy_matches) == 0:
        return False, "No hydroxy groups found"

    # Check if hydroxy groups are attached to the carbon chain
    chain_atoms = set(chain_matches[0])
    for match in hydroxy_matches:
        hydroxy_atom = mol.GetAtoms()[match][0]
        if not any(bond.GetOtherAtomIdx(hydroxy_atom.GetIdx()) in chain_atoms for bond in hydroxy_atom.GetBonds()):
            return False, "Hydroxy group(s) not attached to aliphatic chain"

    # Additional checks
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Chain too short to be a fatty acid"

    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200:
        return False, "Molecular weight too low for a fatty acid"

    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count < 10:
        return False, "Too few carbons for a fatty acid"

    return True, "Contains carboxylic acid group, aliphatic carbon chain, and hydroxy group(s) attached to the chain"