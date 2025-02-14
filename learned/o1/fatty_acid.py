"""
Classifies: CHEBI:35366 fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid(smiles: str):
    """
    Determines if a molecule is a fatty acid based on its SMILES string.
    A fatty acid is any aliphatic monocarboxylic acid derived from or contained
    in esterified form in an animal or vegetable fat, oil or wax.
    Natural fatty acids commonly have a chain of 4 to 28 carbons (usually unbranched
    and even-numbered), which may be saturated or unsaturated.
    By extension, the term is sometimes used to embrace all acyclic aliphatic carboxylic acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carboxylic acid group(s)
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[O;H1,-]')
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    num_carboxylic_acids = len(carboxylic_acid_matches)
    if num_carboxylic_acids == 0:
        return False, "No carboxylic acid group found"

    # Allow molecules with one or more carboxylic acid groups (mono- or diacids)
    # Exclude molecules with functional groups not typical in fatty acids
    common_atoms = {1, 6, 7, 8, 9, 15, 16, 17, 35, 53}  # H, C, N, O, F, P, S, Cl, Br, I
    uncommon_atoms = []
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() not in common_atoms:
            uncommon_atoms.append(atom.GetSymbol())
    if uncommon_atoms:
        return False, f"Contains uncommon atoms ({', '.join(set(uncommon_atoms))}) not typical in fatty acids"

    # Check if the molecule is primarily aliphatic
    # Calculate the number of aromatic rings
    aromatic_rings = mol.GetRingInfo().NumAromaticRings()
    if aromatic_rings > 0:
        return False, "Contains aromatic rings, not typical for fatty acids"

    # Check the carbon chain length (longest aliphatic chain)
    aliphatic_carbons = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and not atom.IsInRing()]
    if not aliphatic_carbons:
        return False, "No aliphatic carbon chain found"

    # Find the longest carbon chain
    paths = Chem.GetSymmSSSR(mol)
    max_chain_length = 0
    for bond_path in Chem.FindAllPathsOfLengthN(mol, N=30, useBonds=True):
        chain_length = 0
        for bond_idx in bond_path:
            bond = mol.GetBondWithIdx(bond_idx)
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetAtomicNum() == 6 and end_atom.GetAtomicNum() == 6:
                chain_length += 1
        if chain_length > max_chain_length:
            max_chain_length = chain_length

    # If unable to find chain length using bonds, fall back to counting aliphatic carbons
    if max_chain_length == 0:
        max_chain_length = len(aliphatic_carbons)

    if max_chain_length < 4:
        return False, f"Aliphatic chain too short ({max_chain_length} carbons)"
    if max_chain_length > 28:
        return False, f"Aliphatic chain too long ({max_chain_length} carbons)"

    return True, "Molecule is classified as a fatty acid"