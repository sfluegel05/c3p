from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
from rdkit.Chem import rdMolDescriptors

def is_alkatetraene(smiles: str):
    """
    Determines if a molecule is an alkatetraene (acyclic hydrocarbon with 4 C=C double bonds)
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is an alkatetraene, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
        
    # Check if molecule is cyclic
    if rdMolDescriptors.CalcNumRings(mol) > 0:
        return False, "Molecule contains rings"
        
    # Check if molecule contains only C and H
    atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
    if not all(symbol == 'C' for symbol in atoms):
        return False, "Molecule contains non-carbon atoms"
        
    # Count number of double bonds
    double_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            double_bonds += 1
            
    if double_bonds != 4:
        return False, f"Molecule has {double_bonds} double bonds, not 4"
        
    # Check if all double bonds are C=C
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
            begin_atom = bond.GetBeginAtom()
            end_atom = bond.GetEndAtom()
            if begin_atom.GetSymbol() != 'C' or end_atom.GetSymbol() != 'C':
                return False, "Double bonds are not all C=C"
                
    return True, "Acyclic hydrocarbon with 4 C=C double bonds"
# Pr=1.0
# Recall=1.0