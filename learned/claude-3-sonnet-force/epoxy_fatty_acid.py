"""
Classifies: CHEBI:61498 epoxy fatty acid
"""
"""
Classifies: CHEBI:46787 epoxy fatty acid
A heterocyclic fatty acid containing an epoxide ring as part of its structure.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, AllChem

def is_epoxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is an epoxy fatty acid based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an epoxy fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for epoxide pattern (O1C[C@]2CCCC[C@@]2)
    epoxide_pattern = Chem.MolFromSmarts("O1C[C@]2CCCC[C@@]2")
    if not mol.HasSubstructMatch(epoxide_pattern):
        return False, "No epoxide ring found"

    # Look for carboxylic acid pattern (-C(=O)O)
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(acid_pattern):
        return False, "No carboxylic acid group found"

    # Check if epoxide and acid are connected
    epoxide_atoms = mol.GetSubstructMatches(epoxide_pattern)
    acid_atoms = mol.GetSubstructMatches(acid_pattern)
    for epoxide_atom in epoxide_atoms:
        for acid_atom in acid_atoms:
            path = Chem.FindPathBetweenFragments(mol, epoxide_atom, acid_atom)
            if path:
                break
        else:
            continue
        break
    else:
        return False, "Epoxide ring and carboxylic acid group not connected"

    # Look for carbon chain (at least 8 carbons)
    carbon_chain_pattern = Chem.MolFromSmarts("[C;H3,H2,H1,H0]~[C;H3,H2,H1,H0]~[C;H3,H2,H1,H0]~[C;H3,H2,H1,H0]~[C;H3,H2,H1,H0]~[C;H3,H2,H1,H0]~[C;H3,H2,H1,H0]~[C;H3,H2,H1,H0]")
    if not mol.HasSubstructMatch(carbon_chain_pattern):
        return False, "Carbon chain too short for fatty acid"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)

    if c_count < 8:
        return False, "Too few carbons for fatty acid"
    if o_count < 2:
        return False, "Too few oxygens for epoxy fatty acid"

    return True, "Contains an epoxide ring, a carboxylic acid group, and a long carbon chain, all connected"