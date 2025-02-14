"""
Classifies: CHEBI:25608 nucleoside phosphate
"""
"""
Classifies: nucleoside phosphate
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_nucleoside_phosphate(smiles: str):
    """
    Determines if a molecule is a nucleoside phosphate based on its SMILES string.
    A nucleoside phosphate is a nucleobase-containing molecule where one or more of
    the sugar hydroxy groups has been converted into a mono- or poly-phosphate.
    This includes both nucleotides and non-nucleotide nucleoside phosphates.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule is a nucleoside phosphate, False otherwise.
        str: Reason for classification.
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define a generalized nucleobase pattern: aromatic heterocycle with nitrogen atoms
    nucleobase_smarts = '[nR][cR]1[cR,nR][cR,nR][cR,nR][cR,nR][cR,nR]1'  # Simplified pattern for purine/pyrimidine rings
    nucleobase_pattern = Chem.MolFromSmarts(nucleobase_smarts)

    # Check for nucleobase
    if not mol.HasSubstructMatch(nucleobase_pattern):
        return False, "No nucleobase found"

    # Define a generalized furanose sugar ring (five-membered ring containing oxygen)
    sugar_smarts = '[#6R1]-[#8R1]-[#6R1]-[#6R1]-[#6R1]'  # Five-membered ring with one oxygen
    sugar_pattern = Chem.MolFromSmarts(sugar_smarts)

    # Check for sugar moiety
    if not mol.HasSubstructMatch(sugar_pattern):
        return False, "No sugar moiety found"

    # Check for glycosidic bond between sugar and nucleobase
    # Find all sugar rings
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    nucleobase_matches = mol.GetSubstructMatches(nucleobase_pattern)

    glycosidic_bond_found = False
    for sugar_match in sugar_matches:
        sugar_atom_indices = set(sugar_match)
        for nucleobase_match in nucleobase_matches:
            nucleobase_atom_indices = set(nucleobase_match)
            # Check for bonds between sugar and nucleobase atoms
            for sugar_atom_idx in sugar_atom_indices:
                sugar_atom = mol.GetAtomWithIdx(sugar_atom_idx)
                for bond in sugar_atom.GetBonds():
                    neighbor = bond.GetOtherAtom(sugar_atom)
                    neighbor_idx = neighbor.GetIdx()
                    if neighbor_idx in nucleobase_atom_indices:
                        glycosidic_bond_found = True
                        break
                if glycosidic_bond_found:
                    break
            if glycosidic_bond_found:
                break
        if glycosidic_bond_found:
            break

    if not glycosidic_bond_found:
        return False, "No glycosidic bond between nucleobase and sugar found"

    # Define phosphate group pattern
    phosphate_smarts = 'P(=O)(O)[O-0]'  # Phosphate group with at least two oxygens
    phosphate_pattern = Chem.MolFromSmarts(phosphate_smarts)

    # Check for phosphate group attached to sugar
    phosphate_found = False
    for sugar_match in sugar_matches:
        for sugar_atom_idx in sugar_match:
            sugar_atom = mol.GetAtomWithIdx(sugar_atom_idx)
            if sugar_atom.GetAtomicNum() == 8:  # Oxygen atom in sugar ring
                for bond in sugar_atom.GetBonds():
                    neighbor = bond.GetOtherAtom(sugar_atom)
                    if neighbor.GetAtomicNum() == 15:  # Phosphorus atom
                        # Check if neighbor matches phosphate pattern
                        phosphate_group = Chem.MolFromSmarts('[O][P](=O)([O])[O]')
                        match = mol.GetSubstructMatch(phosphate_group, maxMatches=1, useChirality=False)
                        if match:
                            phosphate_found = True
                            break
                if phosphate_found:
                    break
        if phosphate_found:
            break

    if not phosphate_found:
        return False, "No phosphate group attached to sugar"

    return True, "Molecule is a nucleoside phosphate with nucleobase, sugar, and phosphate group(s)"