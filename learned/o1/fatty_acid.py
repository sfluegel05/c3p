"""
Classifies: CHEBI:35366 fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_fatty_acid(smiles: str):
    """
    Determines if a molecule is a fatty acid based on its SMILES string.
    A fatty acid is any aliphatic monocarboxylic acid derived from or contained
    in esterified form in an animal or vegetable fat, oil, or wax.
    Natural fatty acids commonly have a chain of 4 to 28 carbons (usually unbranched
    and even-numbered), which may be saturated or unsaturated.
    The term can also include all acyclic aliphatic carboxylic acids.
    
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

    # Check for carboxylic acid group
    carboxylic_acid_pattern = Chem.MolFromSmarts('C(=O)[O,H]')  # Matches both protonated and deprotonated
    carboxylic_acid_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    if len(carboxylic_acid_matches) == 0:
        return False, "No carboxylic acid group found"

    # Ensure only one carboxylic acid group (monocarboxylic acid)
    if len(carboxylic_acid_matches) > 1:
        return False, "More than one carboxylic acid group found"

    # Identify the carboxylic acid carbon atom
    ca_atoms = []
    for match in carboxylic_acid_matches:
        ca_atoms.append(match[0])  # First atom is the carbonyl carbon

    # Exclude molecules with amide bonds or peptide bonds
    amide_pattern = Chem.MolFromSmarts('C(=O)N')
    if mol.HasSubstructMatch(amide_pattern):
        return False, "Contains amide bond, not typical for fatty acids"

    # Exclude molecules with ester groups (excluding the possibility of esterified fatty acids)
    ester_pattern = Chem.MolFromSmarts('C(=O)O[C]')
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) > 0:
        return False, "Contains ester group, not a free fatty acid"

    # Find the longest aliphatic carbon chain connected to the carboxylic acid carbon
    max_chain_length = 0
    for ca_atom in ca_atoms:
        paths = Chem.FindAllPathsOfLengthN(mol, 100, useBonds=False, useHS=False)
        for path in paths:
            if ca_atom in path:
                # Check if path is aliphatic
                is_aliphatic = True
                for idx in path:
                    atom = mol.GetAtomWithIdx(idx)
                    if atom.GetAtomicNum() != 6:  # Not carbon
                        is_aliphatic = False
                        break
                    if atom.IsAromatic():
                        is_aliphatic = False
                        break
                if is_aliphatic:
                    # Update max chain length
                    chain_length = len(path)
                    if chain_length > max_chain_length:
                        max_chain_length = chain_length

    if max_chain_length < 4:
        return False, f"Aliphatic chain too short ({max_chain_length} carbons)"

    # Check for presence of heteroatoms (excluding the carboxyl group)
    heteroatoms = set()
    for atom in mol.GetAtoms():
        atomic_num = atom.GetAtomicNum()
        if atomic_num not in [1, 6, 8]:  # Allowing only H, C, O
            heteroatoms.add(atom.GetSymbol())
    if heteroatoms:
        allowed_heteroatoms = {'O', 'F', 'Cl', 'Br', 'I', 'S', 'N'}
        if not heteroatoms.issubset(allowed_heteroatoms):
            return False, f"Contains uncommon heteroatoms: {', '.join(heteroatoms)}"

    # Allow molecules with rings if they are small and not aromatic
    ring_info = mol.GetRingInfo()
    aromatic_ring_count = 0
    for ring in ring_info.AtomRings():
        if any(mol.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring):
            aromatic_ring_count += 1

    if aromatic_ring_count > 1:
        return False, "Contains multiple aromatic rings, not typical for fatty acids"

    # If it has rings, ensure that the ring does not disrupt the main aliphatic chain
    # For simplicity, we can allow rings if the main chain length is sufficient
    if max_chain_length < 4 and ring_info.NumRings() > 0:
        return False, "Main aliphatic chain is too short and contains rings"

    return True, "Molecule is classified as a fatty acid"