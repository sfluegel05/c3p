"""
Classifies: CHEBI:35346 11beta-hydroxy steroid
"""
"""
Classifies: 11beta-hydroxy steroid
"""
from rdkit import Chem
from rdkit.Chem import rdMolTransforms

def is_11beta_hydroxy_steroid(smiles: str):
    """
    Determines if a molecule is an 11beta-hydroxy steroid based on its SMILES string.
    An 11beta-hydroxy steroid is any steroid that has a hydroxy group at position 11 with beta-configuration.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an 11beta-hydroxy steroid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define steroid nucleus SMARTS pattern (rings A/B/C/D)
    steroid_pattern = Chem.MolFromSmarts('''
    C1CCC2C(C1)CCC3C2CCC4C3CCCC4
    ''')
    if not mol.HasSubstructMatch(steroid_pattern):
        return False, "Does not contain steroid nucleus"

    # Template steroid molecule with atom numbering
    template_smi = '''
    [H][C@]12CC[C@H]3C[C@@H](CC4=CC(=O)CC[C@]34C)[C@@H]1CC[C@H]2O
    '''
    template = Chem.MolFromSmiles(template_smi)
    # Add hydrogens to template and molecule
    template = Chem.AddHs(template)
    mol = Chem.AddHs(mol)

    # Generate 3D coordinates for alignment
    Chem.EmbedMolecule(template)
    Chem.EmbedMolecule(mol)

    # Align mol to template
    match = mol.GetSubstructMatch(template)
    if not match:
        return False, "Unable to align steroid nucleus"

    # Map atom indices from template to mol
    atom_mapping = {template_idx: mol_idx for template_idx, mol_idx in enumerate(match)}

    # Atom 10 in template is position 11 (0-based indexing)
    atom_11_idx = atom_mapping.get(10)
    if atom_11_idx is None:
        return False, "Position 11 not found"

    atom_11 = mol.GetAtomWithIdx(atom_11_idx)

    # Check for hydroxy group at position 11
    has_oh = False
    for neighbor in atom_11.GetNeighbors():
        if neighbor.GetAtomicNum() == 8:
            # Check if oxygen is connected to hydrogen (OH group)
            for n in neighbor.GetNeighbors():
                if n.GetAtomicNum() == 1:
                    has_oh = True
                    oxygen_idx = neighbor.GetIdx()
                    break
            if has_oh:
                break
    if not has_oh:
        return False, "No hydroxy group at position 11"

    # Check stereochemistry at position 11
    chiral_tag = atom_11.GetChiralTag()
    if chiral_tag == Chem.CHI_UNSPECIFIED:
        return False, "No stereochemistry specified at position 11"

    # Determine if configuration is beta (usually R configuration in steroids)
    stereo = Chem.FindMolChiralCenters(mol, includeUnassigned=False, includeCIP=True)
    atom11_stereo = None
    for idx, config in stereo:
        if idx == atom_11_idx:
            atom11_stereo = config
            break
    if atom11_stereo is None:
        return False, "No chiral center at position 11"

    if atom11_stereo != 'R':
        return False, f"Wrong stereochemistry at position 11: {atom11_stereo}, expected beta (R)"

    return True, "Contains hydroxy group at position 11 with beta-configuration"