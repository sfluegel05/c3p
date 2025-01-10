"""
Classifies: CHEBI:17517 phosphatidylglycerol
"""
"""
Classifies: phosphatidylglycerol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_phosphatidylglycerol(smiles: str):
    """
    Determines if a molecule is a phosphatidylglycerol based on its SMILES string.
    Phosphatidylglycerols consist of a glycerol backbone with fatty acyl chains at sn-1 and sn-2 positions,
    and a phosphatidyl group attached to another glycerol at the sn-3 position.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for phosphatidylglycerol headgroup
    # Pattern for phosphate connected to glycerol
    pg_headgroup = Chem.MolFromSmarts("COP(=O)(OCC(O)CO)O")
    if not mol.HasSubstructMatch(pg_headgroup):
        return False, "No phosphatidylglycerol headgroup found"

    # Check for glycerol backbone with two ester-linked fatty acyl chains
    # Pattern for glycerol backbone esterified at sn-1 and sn-2 positions
    esterified_glycerol = Chem.MolFromSmarts("C(COC(=O)[#6])[C@@H](O)COP(=O)(OCC(O)CO)O")
    if not mol.HasSubstructMatch(esterified_glycerol):
        return False, "No glycerol backbone with two ester-linked fatty acids found"

    # Find all ester bonds in the molecule
    ester_bond_pattern = Chem.MolFromSmarts("C(=O)O[C;!$(*C(=O))]")
    ester_bond_matches = mol.GetSubstructMatches(ester_bond_pattern)
    if len(ester_bond_matches) < 2:
        return False, f"Less than two ester bonds found ({len(ester_bond_matches)} found)"

    # Analyze fatty acyl chains attached via ester bonds
    fatty_acid_lengths = []
    for match in ester_bond_matches:
        ester_oxygen_idx = match[2]  # Index of the oxygen in C(=O)O
        fatty_acid_length = 0
        visited_atoms = set()
        atoms_to_visit = [mol.GetAtomWithIdx(ester_oxygen_idx)]
        while atoms_to_visit:
            atom = atoms_to_visit.pop()
            if atom.GetIdx() in visited_atoms:
                continue
            visited_atoms.add(atom.GetIdx())
            if atom.GetAtomicNum() == 6:  # Carbon atom
                fatty_acid_length += 1
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIdx() not in visited_atoms and neighbor.GetAtomicNum() in [6, 1]:
                        atoms_to_visit.append(neighbor)
            elif atom.GetAtomicNum() == 1:  # Hydrogen atom
                continue
            else:
                for neighbor in atom.GetNeighbors():
                    if neighbor.GetIdx() not in visited_atoms:
                        atoms_to_visit.append(neighbor)
        fatty_acid_lengths.append(fatty_acid_length)

    if len(fatty_acid_lengths) < 2:
        return False, "Less than two fatty acyl chains found"
    if min(fatty_acid_lengths) < 8:
        return False, f"Fatty acyl chains are too short (lengths: {fatty_acid_lengths})"

    # Passed all checks
    return True, "Molecule matches phosphatidylglycerol structure"

__metadata__ = {   
    'chemical_class': {   
        'name': 'phosphatidylglycerol',
        'definition': 'A glycerophosphoglycerol that is glycerol in which the hydrogen of one of the primary hydroxy groups has been replaced by a phosphatidyl group.',
    },
}