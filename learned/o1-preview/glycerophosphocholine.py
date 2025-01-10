"""
Classifies: CHEBI:36313 glycerophosphocholine
"""
"""
Classifies: glycerophosphocholine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycerophosphocholine(smiles: str):
    """
    Determines if a molecule is a glycerophosphocholine based on its SMILES string.
    A glycerophosphocholine has a glycerol backbone with fatty acid chains attached
    at positions 1 and 2, and a phosphocholine group at position 3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophosphocholine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern: C1-C2-C3 with appropriate attachments
    glycerol_pattern = Chem.MolFromSmarts("C[C@H](O*)C")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"

    # Check for ester or ether bonds at positions 1 and 2
    # Position 1 and 2 O attached to carbonyl (ester) or alkyl (ether)
    ester_or_ether_pattern = Chem.MolFromSmarts("[C;R0][O;R0][C,S][!#1]")
    ester_or_ether_matches = mol.GetSubstructMatches(ester_or_ether_pattern)
    if len(ester_or_ether_matches) < 2:
        return False, "Less than two ester or ether bonds attached to glycerol backbone"

    # Look for phosphocholine group
    phosphocholine_pattern = Chem.MolFromSmarts("P(=O)([O-])OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine group found"

    # Ensure phosphocholine is attached to the glycerol backbone
    phosphocholine_matches = mol.GetSubstructMatch(phosphocholine_pattern)
    glycerol_matches = mol.GetSubstructMatch(glycerol_pattern)
    if not phosphocholine_matches or not glycerol_matches:
        return False, "Phosphocholine or glycerol backbone not properly matched"

    # Verify connectivity between glycerol backbone and phosphocholine group
    glycerol_atom = glycerol_matches[1]  # Central carbon in glycerol
    phospho_oxygen = phosphocholine_matches[2]  # Oxygen connecting to glycerol
    if not mol.GetBondBetweenAtoms(glycerol_atom, phospho_oxygen):
        return False, "Phosphocholine group not connected to glycerol backbone"

    # Check for at least two long carbon chains (fatty acids or alkyl chains)
    # This can be approximated by counting the number of carbons connected via ester or ether bonds
    chain_lengths = []
    for match in ester_or_ether_matches:
        oxygen_atom = match[1]
        attached_fragments = Chem.FragmentOnBonds(mol, [mol.GetBondBetweenAtoms(match[1], match[2]).GetIdx()])
        frags = Chem.GetMolFrags(attached_fragments, asMols=True)
        for frag in frags:
            if frag.HasSubstructMatch(Chem.MolFromSmarts("[C][C][C]")):
                num_carbons = sum(1 for atom in frag.GetAtoms() if atom.GetAtomicNum() == 6)
                chain_lengths.append(num_carbons)
                break

    if len([l for l in chain_lengths if l >= 4]) < 2:
        return False, "Less than two fatty acid chains attached"

    return True, "Molecule is a glycerophosphocholine with glycerol backbone, fatty acid chains, and phosphocholine group"

__metadata__ = {
    'chemical_class': {
        'id': None,
        'name': 'glycerophosphocholine',
        'definition': 'The glycerol phosphate ester of a phosphocholine. A nutrient with many different roles in human health.',
        'parents': []
    }
}