from rdkit import Chem
from rdkit.Chem import AllChem

def is_tocotrienol(smiles: str):
    """
    Determines if a molecule is a tocotrienol (tocol with 3 double bonds in hydrocarbon chain).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a tocotrienol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check for chromanol core structure (more specific pattern)
    chromanol_pattern = Chem.MolFromSmarts('O1[CH2][CH2][CH2]c2c1cccc2')
    if not mol.HasSubstructMatch(chromanol_pattern):
        return False, "Missing chromanol core structure"
    
    # Check for phenol group (hydroxyl attached to aromatic ring)
    phenol_pattern = Chem.MolFromSmarts('Oc1:c:c:c:c:c1')
    if not mol.HasSubstructMatch(phenol_pattern):
        return False, "Missing phenol group"

    # Check for the presence of three double bonds in an isoprenoid chain
    isoprenoid_chain = Chem.MolFromSmarts('[CH3]-C=C-[CH2]-[CH2]-C=C-[CH2]-[CH2]-C=C-[CH2]-[CH2]')
    if not mol.HasSubstructMatch(isoprenoid_chain):
        # Try alternative pattern with different bond configurations
        isoprenoid_chain2 = Chem.MolFromSmarts('C-C=C-C-C-C=C-C-C-C=C-C-C')
        if not mol.HasSubstructMatch(isoprenoid_chain2):
            return False, "Missing required isoprenoid chain with three double bonds"

    # Count double bonds in aliphatic chain
    double_bonds = mol.GetSubstructMatches(Chem.MolFromSmarts('C=C'))
    
    # Filter out aromatic bonds
    aromatic_atoms = set()
    for atom in mol.GetAromaticAtoms():
        aromatic_atoms.add(atom.GetIdx())
        
    aliphatic_double_bonds = []
    for bond in double_bonds:
        if bond[0] not in aromatic_atoms and bond[1] not in aromatic_atoms:
            aliphatic_double_bonds.append(bond)

    if len(aliphatic_double_bonds) < 3:
        return False, f"Found {len(aliphatic_double_bonds)} aliphatic double bonds, expected at least 3"

    # Check for chain attachment at position 2
    chain_attachment = Chem.MolFromSmarts('O1[C]([CH2][CH2]c2c1cccc2)(CC=C)')
    if not mol.HasSubstructMatch(chain_attachment):
        return False, "Missing required side chain attachment at position 2"

    return True, "Valid tocotrienol structure with chromanol core, phenol group, and isoprenoid chain with three double bonds"
# Pr=None
# Recall=0.0