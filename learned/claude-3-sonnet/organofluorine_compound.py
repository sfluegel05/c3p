"""
Classifies: CHEBI:37143 organofluorine compound
"""
"""
Classifies: CHEBI:24467 - Organofluorine compound
An organofluorine compound is a compound containing at least one carbon-fluorine bond.
"""

from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_organofluorine_compound(smiles: str):
    """
    Determines if a molecule is an organofluorine compound based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an organofluorine compound, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count fluorine atoms and carbon-fluorine bonds
    n_fluorine = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 9)
    n_c_f_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBeginAtom().GetAtomicNum() == 6 and
                      bond.GetEndAtom().GetAtomicNum() == 9)
    
    # Check for presence of carbon-fluorine bonds
    if n_c_f_bonds > 0:
        return True, f"Contains {n_c_f_bonds} carbon-fluorine bond(s)"
    
    # Check for fluorinated alkyl groups
    fluorinated_alkyl_pattern = Chem.MolFromSmarts("[C]C(F)(F)(F),[C]C(F)(F)")
    if mol.HasSubstructMatch(fluorinated_alkyl_pattern):
        return True, "Contains fluorinated alkyl group(s)"
    
    # Check for fluorinated aromatic rings
    fluorinated_aromatic_pattern = Chem.MolFromSmarts("c1ccccc1F,c1ccccc1C(F)(F)F")
    if mol.HasSubstructMatch(fluorinated_aromatic_pattern):
        return True, "Contains fluorinated aromatic ring(s)"
    
    # Check for fluorinated heterocycles
    fluorinated_heterocycle_pattern = Chem.MolFromSmarts("*1ccccc1F,*1ccccc1C(F)(F)F")
    if mol.HasSubstructMatch(fluorinated_heterocycle_pattern):
        return True, "Contains fluorinated heterocyclic ring(s)"
    
    # Filter out inorganic compounds and polymers
    if mol.GetNumAtoms() <= 10 or rdMolDescriptors.CalcMolWt(mol) > 1000:
        return False, "Inorganic compound or polymer"
    
    # If no carbon-fluorine bonds or fluorinated groups found
    if n_fluorine == 0:
        return False, "No fluorine atoms present"
    else:
        return False, "Fluorine present but no carbon-fluorine bonds found"