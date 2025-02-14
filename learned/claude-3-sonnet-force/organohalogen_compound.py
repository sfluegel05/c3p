"""
Classifies: CHEBI:17792 organohalogen compound
"""
"""
Classifies: CHEBI:50486 organohalogen compound
An organohalogen compound is defined as a compound containing at least one carbon-halogen bond (where X is a halogen atom).
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdchem

def is_organohalogen_compound(smiles: str):
    """
    Determines if a molecule is an organohalogen compound based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organohalogen compound, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for carbon-halogen bonds
    has_c_x_bond = False
    for bond in mol.GetBonds():
        atom1 = bond.GetBeginAtom()
        atom2 = bond.GetEndAtom()
        if atom1.GetAtomicNum() == 6 and atom2.GetAtomicNum() in [9, 17, 35, 53]:  # Carbon and halogen (F, Cl, Br, I)
            has_c_x_bond = True
            break
        elif atom2.GetAtomicNum() == 6 and atom1.GetAtomicNum() in [9, 17, 35, 53]:
            has_c_x_bond = True
            break

    if has_c_x_bond:
        # Exclude compounds with trifluoromethyl (-CF3) or trifluoroacetamido (-NHCOCF3) groups
        cf3_pattern = Chem.MolFromSmarts("C(F)(F)F")
        nhcocf3_pattern = Chem.MolFromSmarts("C(=O)N(C(F)(F)F)")
        if mol.HasSubstructMatch(cf3_pattern) or mol.HasSubstructMatch(nhcocf3_pattern):
            return False, "Contains halogen but no carbon-halogen bond"

        # Handle halogenated ring systems or aromatic compounds
        aromatic_rings = mol.GetAromaticRings()
        has_halogenated_ring = any(ring.HasAtomWithQuery(lambda atom: atom.GetAtomicNum() in [9, 17, 35, 53]) for ring in aromatic_rings)
        if has_halogenated_ring:
            return True, "Contains halogenated aromatic ring system"

        # Check for perfluorinated alkanes
        perfluorinated = all(atom.GetAtomicNum() in [6, 9] for atom in mol.GetAtoms())
        if perfluorinated:
            return False, "Perfluorinated alkane"

        return True, "Contains at least one carbon-halogen bond"
    else:
        return False, "No carbon-halogen bonds found"

    # Additional checks for specific functional groups or substructures
    # ...

    # Exclude non-carbon-halogen bonds (e.g., halogen-nitrogen, halogen-oxygen)
    # ...

    # Utilize RDKit functionalities for advanced structure analysis
    # ...

    # Handle exceptional cases or incorporate domain knowledge
    # ...