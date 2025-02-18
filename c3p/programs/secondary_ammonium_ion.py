"""
Classifies: CHEBI:137419 secondary ammonium ion
"""
"""
Classifies: CHEBI:39143 secondary ammonium ion
"""
from rdkit import Chem
from rdkit.Chem import Mol

def is_secondary_ammonium_ion(smiles: str):
    """
    Determines if a molecule is a secondary ammonium ion based on its SMILES string.
    A secondary ammonium ion is a cation derived from protonation of a secondary amine (NH with two R groups).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a secondary ammonium ion, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check if molecule is a cation (has positive formal charge)
    if Chem.GetFormalCharge(mol) <= 0:
        return False, "Molecule is not a cation"

    # Iterate through all nitrogen atoms
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7 and atom.GetFormalCharge() == 1:
            # Count carbon neighbors and total hydrogens (explicit + implicit)
            carbon_count = sum(1 for neighbor in atom.GetNeighbors() 
                             if neighbor.GetAtomicNum() == 6)
            total_h = atom.GetTotalNumHs()
            
            # Check valence (should be 4 for sp3 hybridized N+)
            # Secondary ammonium: 2 carbons + 2 hydrogens = 4 valence
            if (carbon_count == 2 and 
                total_h == 2 and 
                atom.GetTotalValence() == 4 and 
                atom.GetHybridization() == Chem.HybridizationType.SP3):
                
                # Check for electron-withdrawing groups:
                # Direct neighbors (N/O) or carbons in carbonyl groups
                has_ewg = False
                for neighbor in atom.GetNeighbors():
                    # Direct N/O neighbors
                    if neighbor.GetAtomicNum() in [7, 8]:
                        has_ewg = True
                        break
                    # Check for carbonyl groups (C=O) attached to neighbor carbons
                    if neighbor.GetAtomicNum() == 6:
                        for bond in neighbor.GetBonds():
                            if bond.GetBondType() == Chem.BondType.DOUBLE:
                                other_atom = bond.GetOtherAtom(neighbor)
                                if other_atom.GetAtomicNum() == 8:
                                    has_ewg = True
                                    break
                        if has_ewg:
                            break
                
                if not has_ewg:
                    return True, "Positively charged nitrogen with two carbon substituents and two hydrogens"
    
    return False, "No secondary ammonium group (N+ with two carbons and two hydrogens) found"