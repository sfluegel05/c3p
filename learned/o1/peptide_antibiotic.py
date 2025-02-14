"""
Classifies: CHEBI:25903 peptide antibiotic
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is a peptide antibiotic based on its SMILES string.
    A peptide antibiotic is a peptide that exhibits antimicrobial properties.

    This function checks for features common in peptide antibiotics:
    - Presence of multiple peptide bonds
    - The molecule is cyclic (cyclic peptide)
    - Presence of thioether bridges or unusual amino acid modifications

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is likely a peptide antibiotic, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define peptide bond pattern (amide bond between C=O and N)
    peptide_bond_pattern = Chem.MolFromSmarts("N[C;!H0]=[O,N]")  # Tweaked pattern for peptide bond
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)
    num_peptide_bonds = len(peptide_bonds)

    if num_peptide_bonds < 5:
        return False, f"Contains {num_peptide_bonds} peptide bonds; too few for peptide antibiotic"

    # Check for cyclic peptide (peptide bonds forming a ring)
    ring_info = mol.GetRingInfo()
    is_cyclic_peptide = False
    for bond in mol.GetBonds():
        if bond.IsInRing():
            atom1 = bond.GetBeginAtom()
            atom2 = bond.GetEndAtom()
            # Check if the bond is a peptide bond in a ring
            if ((atom1.GetSymbol() == 'C' and atom2.GetSymbol() == 'N') or
                (atom1.GetSymbol() == 'N' and atom2.GetSymbol() == 'C')):
                if bond.GetBondType() == Chem.rdchem.BondType.SINGLE or bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    is_cyclic_peptide = True
                    break

    if not is_cyclic_peptide:
        return False, "Not a cyclic peptide"

    # Define thioether bridge pattern (e.g., lanthionine bridge)
    thioether_pattern = Chem.MolFromSmarts("C-S-C")
    has_thioether_bridge = mol.HasSubstructMatch(thioether_pattern)

    # Check for unusual amino acids or modifications (e.g., dehydroalanine)
    dehydroalanine_pattern = Chem.MolFromSmarts("C=CN")
    has_unusual_amino_acids = mol.HasSubstructMatch(dehydroalanine_pattern)

    # If either thioether bridges or unusual amino acids are present, it strengthens the case
    if has_thioether_bridge or has_unusual_amino_acids:
        return True, "Likely a peptide antibiotic (cyclic peptide with characteristic modifications)"

    return False, "Does not contain characteristic modifications of peptide antibiotics"