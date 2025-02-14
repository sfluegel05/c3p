"""
Classifies: CHEBI:28963 amino sugar
"""
"""
Classifies: CHEBI:17754 amino sugar
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_amino_sugar(smiles: str):
    """
    Determines if a molecule is an amino sugar based on its SMILES string.
    An amino sugar is any sugar having one or more alcoholic hydroxy groups
    replaced by substituted or unsubstituted amino groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an amino sugar, False otherwise
        str: Reason for classification
    """
    # Parse SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Normalize molecule (add hydrogens)
    mol = Chem.AddHs(mol)

    # Define SMARTS pattern for pyranose and furanose rings
    pyranose_pattern = Chem.MolFromSmarts("OC1[C@@H]([O,H])[C@@H]([O,H])[C@H]([O,H])[C@@H]([O,H])O1")
    furanose_pattern = Chem.MolFromSmarts("OC1[C@H]([O,H])[C@@H]([O,H])[C@H]([O,H])O1")

    # Define pattern for amino sugar (hydroxy group replaced by amino group)
    amino_substitution_pattern = Chem.MolFromSmarts("[C;R][N;!H0]")  # Ring carbon bonded to nitrogen with at least one hydrogen or substituent

    # Flag to track if amino sugar is found
    is_amino_sugar_found = False

    # Find sugar rings (pyranose and furanose)
    sugar_rings = []
    pyranose_matches = mol.GetSubstructMatches(pyranose_pattern)
    furanose_matches = mol.GetSubstructMatches(furanose_pattern)
    sugar_rings.extend(pyranose_matches)
    sugar_rings.extend(furanose_matches)

    # Check for cyclic amino sugars
    if sugar_rings:
        for ring in sugar_rings:
            ring_atoms = set(ring)
            amino_replaced = False
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 6:  # Carbon atom
                    # Check for hydroxy group
                    has_hydroxy = False
                    has_amino = False
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetAtomicNum() == 8:
                            # Oxygen atom could be hydroxy group
                            if mol.GetBondBetweenAtoms(atom.GetIdx(), neighbor.GetIdx()).GetBondType() == Chem.rdchem.BondType.SINGLE:
                                has_hydroxy = True
                        elif neighbor.GetAtomicNum() == 7:
                            has_amino = True
                    # If amino group is present where hydroxy group is expected
                    if has_amino and not has_hydroxy:
                        amino_replaced = True
            if amino_replaced:
                return True, "Contains sugar ring with amino group(s) replacing hydroxy group(s)"

    # Check for acyclic amino sugars
    # Define acyclic sugar backbone pattern (open-chain monosaccharides)
    acyclic_sugar_pattern = Chem.MolFromSmarts("[H][C@@H](O)[C@@H](O)[C@H](O)[C@H](O)[C@H](O)C=O")  # Aldoses
    acyclic_sugar_ketone_pattern = Chem.MolFromSmarts("[H][C@@H](O)[C@@H](O)[C@H](O)[C@@H](O)C(=O)C")  # Ketoses

    acyclic_matches = mol.GetSubstructMatches(acyclic_sugar_pattern) + mol.GetSubstructMatches(acyclic_sugar_ketone_pattern)

    if acyclic_matches:
        for match in acyclic_matches:
            amino_replaced = False
            for idx in match:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 6 and atom.GetDegree() > 1:
                    # Check for hydroxy group
                    has_hydroxy = False
                    has_amino = False
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetAtomicNum() == 8:
                            has_hydroxy = True
                        elif neighbor.GetAtomicNum() == 7:
                            has_amino = True
                    if has_amino and not has_hydroxy:
                        amino_replaced = True
            if amino_replaced:
                return True, "Contains acyclic amino sugar with amino group(s) replacing hydroxy group(s)"

    # Additional check: any sugar ring with amino substitution
    # Define a generic sugar ring pattern
    generic_sugar_pattern = Chem.MolFromSmarts("C1[C,O][C,O][C,O][C,O][C,O]1")
    generic_sugar_matches = mol.GetSubstructMatches(generic_sugar_pattern)
    if generic_sugar_matches:
        for ring in generic_sugar_matches:
            ring_atoms = set(ring)
            amino_replaced = False
            for idx in ring:
                atom = mol.GetAtomWithIdx(idx)
                if atom.GetAtomicNum() == 6:
                    has_amino = False
                    has_hydroxy = False
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetAtomicNum() == 8:
                            has_hydroxy = True
                        elif neighbor.GetAtomicNum() == 7:
                            has_amino = True
                    if has_amino and not has_hydroxy:
                        amino_replaced = True
            if amino_replaced:
                return True, "Contains sugar ring with amino group(s) replacing hydroxy group(s)"

    return False, "Does not contain amino sugar ring or acyclic amino sugar"