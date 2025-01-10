"""
Classifies: CHEBI:36141 quinone
"""
"""
Classifies: quinone
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_quinone(smiles: str):
    """
    Determines if a molecule is a quinone based on its SMILES string.
    A quinone is defined as a compound with a fully conjugated cyclic dione structure
    derived from an aromatic compound by conversion of an even number of -CH= groups
    into -C(=O)- groups with any necessary rearrangement of double bonds.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a quinone, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Kekulize the molecule to get explicit double bonds
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Chem.KekulizeException:
        return False, "Failed to kekulize molecule"
    
    # Define a SMARTS pattern for a conjugated cyclic dione (quinone core)
    quinone_pattern = Chem.MolFromSmarts('C1=CC(=O)C=CC(=O)C=C1')  # 1,4-benzoquinone core
    # Extend the pattern to include fused ring systems (naphthoquinones, anthraquinones)
    quinone_patterns = [
        Chem.MolFromSmarts('C1=CC(=O)C=CC(=O)C=C1'),  # para-benzoquinone
        Chem.MolFromSmarts('C1=CC(=O)C=CC(=C1)C=O'),  # ortho-benzoquinone
        Chem.MolFromSmarts('C1=CC=CC2=C1C(=O)C=CC2=O'),  # naphthoquinone
        Chem.MolFromSmarts('C1=CC=CC2=C1C=CC3=C2C(=O)C=CC3=O'),  # anthraquinone
    ]
    
    # Check if the molecule matches any of the quinone patterns
    for pattern in quinone_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains a quinone core structure"
    
    # Alternatively, search for aromatic rings with two ketone substitutions
    aromatic_rings = []
    for ring in mol.GetRingInfo().AtomRings():
        if all(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_rings.append(ring)
    
    if not aromatic_rings:
        return False, "No aromatic rings found"
    
    for ring in aromatic_rings:
        ketone_count = 0
        for atom_idx in ring:
            atom = mol.GetAtomWithIdx(atom_idx)
            # Check for carbonyl carbons within the ring
            if atom.GetAtomicNum() == 6 and atom.GetDegree() == 3:
                # Look for C=O double bond
                for nbr in atom.GetNeighbors():
                    bond = mol.GetBondBetweenAtoms(atom_idx, nbr.GetIdx())
                    if nbr.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        ketone_count += 1
                        break
        if ketone_count == 2:
            return True, "Aromatic ring with two ketone groups found"
    
    return False, "Does not contain a quinone core structure"
    

__metadata__ = {
    'chemical_class': {
        'name': 'quinone',
        'definition': 'Compounds having a fully conjugated cyclic dione structure, such as that of benzoquinones, derived from aromatic compounds by conversion of an even number of -CH= groups into -C(=O)- groups with any necessary rearrangement of double bonds (polycyclic and heterocyclic analogues are included).',
    },
    'examples': [
        'O=C1C=CC(=O)C=CC1=O',  # Benzoquinone
        'O=C1C=CC2=CC=CC=C2C1=O',  # Naphthoquinone
        'O=C1C=CC2=C1C=CC(=O)C=C2',  # Anthraquinone
    ],
    'version': '1.0',
    'author': 'Assistant',
}