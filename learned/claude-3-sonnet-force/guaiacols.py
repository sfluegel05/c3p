"""
Classifies: CHEBI:134251 guaiacols
"""
"""
Classifies: CHEBI:26056 guaiacol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_guaiacol(smiles: str):
    """
    Determines if a molecule is a guaiacol based on its SMILES string.
    A guaiacol is any phenol carrying an additional methoxy substituent at the ortho-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a guaiacol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find phenol substructure
    phenol_pattern = Chem.MolFromSmarts("c1ccc(O)cc1")
    phenol_matches = mol.GetSubstructMatches(phenol_pattern)
    
    if not phenol_matches:
        return False, "No phenol substructure found"
    
    # Find methoxy substituent
    methoxy_pattern = Chem.MolFromSmarts("COc")
    methoxy_matches = mol.GetSubstructMatches(methoxy_pattern)
    
    if not methoxy_matches:
        return False, "No methoxy substituent found"
    
    # Check if methoxy is in ortho position to hydroxyl group
    for phenol_match in phenol_matches:
        hydroxy_atom = phenol_match[Chem.MolFromSmarts("c1ccc(O)cc1").GetAtomMapNumber(5) - 1]
        for methoxy_match in methoxy_matches:
            methoxy_atom = methoxy_match[0]
            if mol.GetBondBetweenAtoms(hydroxy_atom, methoxy_atom) is not None:
                # Check if methoxy and hydroxyl are ortho to each other
                bonds = mol.GetBondBetweenAtoms(hydroxy_atom, methoxy_atom).GetBeginAtomIdx(), mol.GetBondBetweenAtoms(hydroxy_atom, methoxy_atom).GetEndAtomIdx()
                ring_info = mol.GetRingInfo()
                for ring in ring_info.AtomRings():
                    if hydroxy_atom in ring and methoxy_atom in ring and len(ring) == 6 and bonds[0] in ring and bonds[1] in ring and abs(ring.index(bonds[0]) - ring.index(bonds[1])) == 1:
                        return True, "Contains a phenol ring with a methoxy substituent at the ortho position"
    
    return False, "No methoxy substituent found at the ortho position of a phenol ring"