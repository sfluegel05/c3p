"""
Classifies: CHEBI:32877 primary amine
"""
"""
Classifies: CHEBI:36975 primary amine
A compound formally derived from ammonia by replacing one hydrogen atom by a hydrocarbyl group.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary amine, False otherwise
        str: Reason for classification
    """
    try:
        # Parse SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False, "Invalid SMILES string"

        # Get list of nitrogen atoms
        n_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]

        # Check if at least one nitrogen atom meets primary amine criteria
        for n_atom in n_atoms:
            # Check for double/triple bonds to nitrogen or other unsaturated bonds
            if sum(bond.GetBondTypeAsDouble() for bond in mol.GetBondEdges(n_atom.GetIdx())) > 1 \
                or any(bond.GetBondType() == Chem.BondType.DOUBLE or bond.GetBondType() == Chem.BondType.TRIPLE for bond in mol.GetBonds()):
                continue

            # Check for exactly one hydrogen atom attached to nitrogen
            n_hydrogens = sum(1 for bond in mol.GetBondEdges(n_atom.GetIdx()) if mol.GetAtomWithIdx(bond[1]).GetAtomicNum() == 1)
            if n_hydrogens != 1:
                continue

            # Check for alkyl or aryl substituent
            alkyl_pattern = Chem.MolFromSmarts("[CH3][CH2]*[NH]")
            aryl_pattern = Chem.MolFromSmarts("a[NH]")
            if not mol.HasSubstructMatch(alkyl_pattern) and not mol.HasSubstructMatch(aryl_pattern):
                continue

            # If all checks pass, this is a primary amine
            return True, "Molecule has characteristic structural features of a primary amine"

        # If no nitrogen atoms meet criteria, return false
        return False, "No nitrogen atoms meet the criteria for a primary amine"

    except Exception as e:
        return False, f"Error: {str(e)}"