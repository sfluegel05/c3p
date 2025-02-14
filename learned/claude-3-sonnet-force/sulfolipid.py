"""
Classifies: CHEBI:61384 sulfolipid
"""
"""
Classifies: CHEBI:36975 sulfolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid based on its SMILES string.
    A sulfolipid is a compound containing a sulfonic acid residue joined by a carbon-sulfur bond to a lipid.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfolipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for sulfonic acid group (-SO3H)
    sulfonic_acid_pattern = Chem.MolFromSmarts("S(=O)(=O)O")
    if not mol.HasSubstructMatch(sulfonic_acid_pattern):
        return False, "No sulfonic acid group found"
    
    # Look for lipid chains (long carbon chains attached to sulfur)
    lipid_chain_pattern = Chem.MolFromSmarts("[SX4]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    if not mol.HasSubstructMatch(lipid_chain_pattern):
        return False, "No lipid chain found attached to sulfur"
    
    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short to be lipid chains"
    
    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 10:
        return False, "Too few carbons for sulfolipid"
    if o_count < 4:
        return False, "Too few oxygens for sulfolipid"
    
    return True, "Contains a sulfonic acid group joined by a carbon-sulfur bond to a lipid chain"