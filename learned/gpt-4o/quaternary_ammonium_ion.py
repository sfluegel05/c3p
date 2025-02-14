"""
Classifies: CHEBI:35267 quaternary ammonium ion
"""
from rdkit import Chem

def is_quaternary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a quaternary ammonium ion (QAI) based on its SMILES string.
    A QAI is typically represented by a positively charged nitrogen bonded typically to four univalent groups, usually organic groups.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a quaternary ammonium ion, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for quaternary ammonium characteristic: [N+] with four single bonds
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:  # Check if the atom is Nitrogen
           
            positive_charge = atom.GetFormalCharge() == 1
            if positive_charge or not all(map(lambda bond: bond.IsInRing(), atom.GetBonds())):
                bonded_atoms = [bond.GetOtherAtom(atom) for bond in atom.GetBonds()]
                
                # Relax previous checks to account for any aliphatic univalent attachment
                if len(bonded_atoms) == 4:  # Ensure four bonds
                    if any(bonded_atom.GetAtomicNum() == 6 for bonded_atom in bonded_atoms):  # Ensure at least one carbon
                        # Avoid limiting only to carbon due to need for diversity in univalency e.g. [F-], [O-]
                        return True, "Contains a quaternary ammonium nitrogen with positive charge"
    
    return False, "Does not contain a quaternary ammonium nitrogen"