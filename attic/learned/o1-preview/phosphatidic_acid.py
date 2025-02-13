"""
Classifies: CHEBI:16337 phosphatidic acid
"""
"""
Classifies: phosphatidic acid
"""
from rdkit import Chem

def is_phosphatidic_acid(smiles: str):
    """
    Determines if a molecule is a phosphatidic acid based on its SMILES string.
    A phosphatidic acid has a glycerol backbone where one hydroxyl group is esterified with phosphoric acid
    and the other two are esterified with fatty acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a phosphatidic acid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns
    # Phosphatidic acid core pattern: glycerol backbone esterified with two fatty acids and one phosphate group
    phosphatidic_acid_pattern = Chem.MolFromSmarts("""
    [C;X4](O[C;X4](CO[P](=O)(O)O)C(=O)[C;X4])[C;X4](OC(=O)[C;X4])[O]
    """)
    # Alternatively, since the above pattern may be too restrictive, we'll piece together checks

    # Step 1: Find glycerol backbone (three connected carbons)
    backbone_pattern = Chem.MolFromSmarts("C-C-C")
    backbone_matches = mol.GetSubstructMatches(backbone_pattern)
    if not backbone_matches:
        return False, "No glycerol backbone found"

    # Step 2: For each backbone, check attachments
    for match in backbone_matches:
        c1_idx, c2_idx, c3_idx = match
        c1 = mol.GetAtomWithIdx(c1_idx)
        c2 = mol.GetAtomWithIdx(c2_idx)
        c3 = mol.GetAtomWithIdx(c3_idx)

        # Check that each carbon has an oxygen attached
        oxygens_c1 = [nbr for nbr in c1.GetNeighbors() if nbr.GetAtomicNum() == 8]
        oxygens_c2 = [nbr for nbr in c2.GetNeighbors() if nbr.GetAtomicNum() == 8]
        oxygens_c3 = [nbr for nbr in c3.GetNeighbors() if nbr.GetAtomicNum() == 8]

        if not oxygens_c1 or not oxygens_c2 or not oxygens_c3:
            continue  # Not all carbons have oxygen attached

        # Step 3: Check attachments of oxygens
        # Initialize counts
        phosphate_attached = False
        ester_count = 0

        # Check oxygen attached to c1
        o_c1 = oxygens_c1[0]
        if is_ester_linkage(mol, o_c1, c1_idx):
            ester_count += 1
        elif is_phosphate_linkage(mol, o_c1):
            phosphate_attached = True

        # Check oxygen attached to c2
        o_c2 = oxygens_c2[0]
        if is_ester_linkage(mol, o_c2, c2_idx):
            ester_count += 1
        elif is_phosphate_linkage(mol, o_c2):
            phosphate_attached = True

        # Check oxygen attached to c3
        o_c3 = oxygens_c3[0]
        if is_ester_linkage(mol, o_c3, c3_idx):
            ester_count += 1
        elif is_phosphate_linkage(mol, o_c3):
            phosphate_attached = True

        if phosphate_attached and ester_count >=2:
            return True, "Contains glycerol backbone with phosphate group and two fatty acid chains attached via ester bonds"

    return False, "Does not match phosphatidic acid structure"

def is_ester_linkage(mol, oxygen_atom, carbon_idx):
    """
    Checks if the oxygen atom is part of an ester linkage connected to the specified carbon atom.

    Args:
        mol: RDKit molecule object
        oxygen_atom: RDKit atom object (oxygen)
        carbon_idx: Index of the carbon atom in the ester linkage

    Returns:
        bool: True if oxygen is part of an ester linkage, False otherwise
    """
    for neighbor in oxygen_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 6 and neighbor.GetIdx() != carbon_idx:
            # Check if this carbon neighbor has a double bond to an oxygen (carbonyl group)
            for bond in neighbor.GetBonds():
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    other_atom = bond.GetOtherAtom(neighbor)
                    if other_atom.GetAtomicNum() == 8:
                        return True
    return False

def is_phosphate_linkage(mol, oxygen_atom):
    """
    Checks if the oxygen atom is connected to a phosphate group.

    Args:
        mol: RDKit molecule object
        oxygen_atom: RDKit atom object (oxygen)

    Returns:
        bool: True if oxygen is connected to a phosphate group, False otherwise
    """
    for neighbor in oxygen_atom.GetNeighbors():
        if neighbor.GetAtomicNum() == 15:
            # Confirm it's a phosphate group by checking for P=O bonds
            double_bonded_oxygens = 0
            single_bonded_oxygens = 0
            for bond in neighbor.GetBonds():
                other_atom = bond.GetOtherAtom(neighbor)
                if other_atom.GetAtomicNum() == 8:
                    if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        double_bonded_oxygens += 1
                    elif bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                        single_bonded_oxygens += 1
            if double_bonded_oxygens >= 1 and single_bonded_oxygens >= 2:
                return True
    return False

__metadata__ = {   'chemical_class': {   'id': 'CHEBI:49105',
                              'name': 'phosphatidic acid',
                              'definition': 'A derivative of glycerol in which one hydroxy group, commonly but not necessarily primary, is esterified with phosphoric acid and the other two are esterified with fatty acids.',
                              'parents': []},
        'config': {   'llm_model_name': 'lbl/claude-sonnet',
                      'f1_threshold': 0.8,
                      'max_attempts': 5},
        'message': None,
        'attempt': 2,
        'success': True,
        'best': True,
        'error': '',
        'stdout': None}