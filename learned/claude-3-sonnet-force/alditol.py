"""
Classifies: CHEBI:17522 alditol
"""
"""
Classifies: CHEBI:17937 alditol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_alditol(smiles: str):
    """
    Determines if a molecule is an alditol based on its SMILES string.
    An alditol is a carbohydrate that is an acyclic polyol with the general formula HOCH2[CH(OH)]nCH2OH.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alditol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for alditol backbone pattern
    alditol_pattern = Chem.MolFromSmarts("[CH2O][CH](O)[CH2](O)[CH2](O)~[CH2](O)~[CH2](O)")
    if not mol.HasSubstructMatch(alditol_pattern):
        return False, "Missing alditol backbone pattern"
    
    # Count hydroxy groups
    num_hydroxy = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "O" and atom.GetHybridization() == Chem.HybridizationType.SP3)
    
    # Count carbon atoms
    num_carbon = mol.GetNumAtoms(onlyObsMask=Chem.AtomProp.Query(["#6"]))
    
    # Check if number of hydroxy groups matches the expected formula
    if num_hydroxy != num_carbon + 2:
        return False, f"Incorrect number of hydroxy groups for alditol (expected {num_carbon + 2}, got {num_hydroxy})"
    
    # Check for absence of other heteroatoms
    het_atom_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() not in ["C", "H", "O"])
    if het_atom_count > 0:
        return False, "Molecule contains heteroatoms other than C, H, O"

    return True, "Acyclic polyol matching the formula HOCH2[CH(OH)]nCH2OH"

__metadata__ = {
    'chemical_class': {
        'id': 'CHEBI:17937',
        'name': 'alditol',
        'definition': 'A carbohydrate that is an acyclic polyol having the general formula HOCH2[CH(OH)]nCH2OH (formally derivable from an aldose by reduction of the carbonyl group).',
        'parents': ['CHEBI:36973', 'CHEBI:16646']
    },
    # ... (rest of metadata omitted for brevity)
}