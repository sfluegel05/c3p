"""
Classifies: CHEBI:64583 sphingomyelin
"""
"""
Classifies: CHEBI:15837 sphingomyelin
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_sphingomyelin(smiles: str):
    """
    Determines if a molecule is a sphingomyelin based on its SMILES string.
    A sphingomyelin is a sphingolipid with a sphingoid base backbone that has:
    - An amide linkage between the amino group and a fatty acid
    - A phosphorylcholine group esterified to the terminal hydroxyl group

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingomyelin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for sphingoid base backbone
    # Pattern for sphingoid base backbone with amino and hydroxyl groups
    sphingoid_pattern = Chem.MolFromSmarts(
        "[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#6]-[#7]-[#6]-[#6]-[#8]"
    )
    if not mol.HasSubstructMatch(sphingoid_pattern):
        return False, "No sphingoid base backbone found"

    # Check for amide linkage with fatty acid
    # Pattern for amide bond attached to the amino group of sphingoid base
    amide_pattern = Chem.MolFromSmarts("[$([NX3][C;R0]=[O;R0])]")
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if len(amide_matches) == 0:
        return False, "No amide linkage with fatty acid found"

    # Check for phosphorylcholine group
    # Pattern for phosphorylcholine group attached via ester linkage
    phosphorylcholine_pattern = Chem.MolFromSmarts(
        "[O][P](=O)([O-])[O][C][C][N+](C)(C)C"
    )
    if not mol.HasSubstructMatch(phosphorylcholine_pattern):
        return False, "No phosphorylcholine group found"

    # Optional: Verify the overall structure
    # Ensure that the phosphorylcholine group is connected to the terminal hydroxyl group
    # and the fatty acid is connected via amide bond to the amino group

    # Check that the amide and phosphorylcholine groups are attached to the sphingoid base
    sphingoid_match = mol.GetSubstructMatch(sphingoid_pattern)
    sphingoid_atoms = set(sphingoid_match)

    # Find amide nitrogen atom
    amide_nitrogens = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7 and any(bond.GetBondType() == Chem.rdchem.BondType.DOUBLE and bond.GetOtherAtom(atom).GetAtomicNum() == 6 for bond in atom.GetBonds())]
    if not amide_nitrogens:
        return False, "Amide nitrogen not found"

    # Verify amide nitrogen is in sphingoid backbone
    if amide_nitrogens[0].GetIdx() not in sphingoid_atoms:
        return False, "Amide nitrogen not part of sphingoid backbone"

    # Find phosphorylcholine oxygen atom
    phosphorylcholine_matches = mol.GetSubstructMatches(phosphorylcholine_pattern)
    if not phosphorylcholine_matches:
        return False, "Phosphorylcholine group not properly attached"
    phosphorylcholine_atoms = phosphorylcholine_matches[0]

    # Get the ester oxygen atom connecting to sphingoid base
    ester_oxygen_idx = phosphorylcholine_atoms[0]
    ester_oxygen = mol.GetAtomWithIdx(ester_oxygen_idx)

    # Verify ester oxygen is connected to terminal carbon of sphingoid base
    ester_bonds = ester_oxygen.GetBonds()
    connected_atoms = [bond.GetOtherAtom(ester_oxygen) for bond in ester_bonds]
    sphingoid_terminal_oxygen = [atom for atom in connected_atoms if atom.GetIdx() in sphingoid_atoms]
    if not sphingoid_terminal_oxygen:
        return False, "Phosphorylcholine group not connected to terminal hydroxyl group"

    return True, "Molecule is a sphingomyelin with correct structural features"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:15837',
        'name': 'sphingomyelin',
        'definition': 'Any of a class of phospholipids in which the amino group of a sphingoid base is in amide linkage with one of several fatty acids, while the terminal hydroxy group of the sphingoid base is esterified to phosphorylcholine.',
        'parents': ['CHEBI:64713', 'CHEBI:65409']
    }
}