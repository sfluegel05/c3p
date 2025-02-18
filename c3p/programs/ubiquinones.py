"""
Classifies: CHEBI:16389 ubiquinones
"""
"""
Classifies: CHEBI:27102 ubiquinones
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_ubiquinones(smiles: str):
    """
    Determines if a molecule is a ubiquinone based on its SMILES string.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a ubiquinone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for benzoquinone core with proper substitution pattern
    # Allows for both methoxy (-OC) and hydroxy (-O) groups
    # Also matches charged forms with [O-]
    benzoquinone_pattern = Chem.MolFromSmarts("O=C1C([O,OC])=C([O,OC])C(=O)C=C1")
    if not mol.HasSubstructMatch(benzoquinone_pattern):
        return False, "No properly substituted benzoquinone core found"

    # At least one position must be methoxy (-OC)
    methoxy_pattern = Chem.MolFromSmarts("O=C1[$(C([O,OC])=C([OC])),$([$(C(OC)=C([O,OC])))]C(=O)C=C1")
    if not mol.HasSubstructMatch(methoxy_pattern):
        return False, "Missing at least one methoxy group"

    # Look for methyl group on the quinone ring
    methyl_pattern = Chem.MolFromSmarts("O=C1C([O,OC])=C([O,OC])C(=O)C(C)=C1")
    if not mol.HasSubstructMatch(methyl_pattern):
        return False, "Missing methyl group on quinone ring"

    # Check for proper side chain
    # Either an isoprenoid unit or a long chain that could be derived from it
    isoprenoid_pattern = Chem.MolFromSmarts("[$(CC=C(C)C),$(CC=C(C)CC),$(CCC=C(C)C)]")
    long_chain_pattern = Chem.MolFromSmarts("CCCCCCC") # At least 7 carbons in chain
    
    if not (mol.HasSubstructMatch(isoprenoid_pattern) or mol.HasSubstructMatch(long_chain_pattern)):
        return False, "Missing required side chain"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if o_count < 4:
        return False, "Insufficient oxygen atoms for ubiquinone"
    
    # Check for fused ring systems that would indicate wrong class
    fused_rings_pattern = Chem.MolFromSmarts("[r6][r5]")
    if mol.HasSubstructMatch(fused_rings_pattern):
        return False, "Contains fused ring system not typical of ubiquinones"

    # Basic structure confirmed
    return True, "Contains benzoquinone core with proper substitution pattern and characteristic side chain"