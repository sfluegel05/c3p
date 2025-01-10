"""
Classifies: CHEBI:26255 prenylquinone
"""
"""
Classifies: CHEBI:XXXXX prenylquinone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_prenylquinone(smiles: str):
    """
    Determines if a molecule is a prenylquinone based on its SMILES string.
    A prenylquinone is a quinone substituted by a polyprenyl-derived side-chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a prenylquinone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for quinone core (1,4-benzoquinone or 1,4-naphthoquinone)
    # More flexible patterns that allow substitutions
    benzoquinone_pattern = Chem.MolFromSmarts("[#6]1([#6]=O)[#6]=[#6][#6]([#6]=O)[#6]=[#6]1")
    naphthoquinone_pattern = Chem.MolFromSmarts("[#6]1([#6]=O)[#6]=[#6][#6]([#6]=O)[#6]=[#6]2[#6]=[#6][#6]=[#6][#6]=[#6]12")
    
    if not mol.HasSubstructMatch(benzoquinone_pattern) and not mol.HasSubstructMatch(naphthoquinone_pattern):
        return False, "No quinone core found"

    # Look for polyprenyl side-chain (isoprenoid chain)
    # More flexible pattern for isoprene units
    isoprene_pattern = Chem.MolFromSmarts("[C;H2]=C(-[C;H3])-[C;H2,H1]")
    isoprene_matches = mol.GetSubstructMatches(isoprene_pattern)
    
    # Check for at least one isoprene unit
    if len(isoprene_matches) < 1:
        return False, f"Found {len(isoprene_matches)} isoprene units, need at least 1"

    # Check if isoprene units are connected in a chain
    # Get all atoms in isoprene units
    isoprene_atoms = set(atom for match in isoprene_matches for atom in match)
    
    # Check connectivity
    connected = False
    for atom_idx in isoprene_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        neighbors = [n.GetIdx() for n in atom.GetNeighbors() if n.GetIdx() in isoprene_atoms]
        if len(neighbors) > 1:
            connected = True
            break

    if not connected:
        return False, "Isoprene units not connected in a chain"

    # Check molecular weight - prenylquinones typically >300 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for prenylquinone"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 15:
        return False, "Too few carbons for prenylquinone"
    if o_count < 2:
        return False, "Must have at least 2 oxygens (quinone core)"

    return True, "Contains quinone core with polyprenyl side-chain"