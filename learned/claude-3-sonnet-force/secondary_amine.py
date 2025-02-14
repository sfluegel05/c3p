"""
Classifies: CHEBI:32863 secondary amine
"""
"""
Classifies: CHEBI:35550 secondary amine
A compound formally derived from ammonia by replacing two hydrogen atoms by hydrocarbyl groups.
"""

from rdkit import Chem
from rdkit.Chem import AllChem

def is_secondary_amine(smiles: str):
    """
    Determines if a molecule is a secondary amine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Count nitrogen atoms
    n_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if len(n_atoms) != 1:
        return False, "Must have exactly one nitrogen atom"
    
    # Get the nitrogen atom
    n_atom = n_atoms[0]
    
    # Check if nitrogen is aromatic
    if n_atom.GetIsAromatic():
        return False, "Nitrogen atom is aromatic, should be aliphatic"
    
    # Get the hydrocarbyl groups attached to the nitrogen atom
    hydrocarbyl_groups = [bond.GetOtherAtom(n_atom) for bond in n_atom.GetBonds()]
    hydrocarbyl_groups = [atom for atom in hydrocarbyl_groups if atom.GetAtomicNum() in [6, 7]]
    
    # Check if there are exactly two hydrocarbyl groups
    if len(hydrocarbyl_groups) != 2:
        return False, "Nitrogen atom does not have exactly two hydrocarbyl groups attached"
    
    # Check if the hydrocarbyl groups are alkyl or aryl
    for atom in hydrocarbyl_groups:
        if atom.GetAtomicNum() == 6:  # Carbon
            if not any(bond.GetIsAromatic() for bond in atom.GetBonds()):
                continue  # Alkyl group
            else:
                aromatic_rings = Chem.GetSymmSSSR(mol)
                if not any(atom.IsInRingOfSize(r.NumAtoms()) for r in aromatic_rings):
                    return False, "Hydrocarbyl group is not alkyl or aryl"
        else:  # Nitrogen
            return False, "Nitrogen atom attached to another nitrogen atom"
    
    return True, "Contains a secondary aliphatic amine group"