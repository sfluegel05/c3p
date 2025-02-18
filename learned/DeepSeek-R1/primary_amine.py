"""
Classifies: CHEBI:32877 primary amine
"""
"""
Classifies: CHEBI:15734 primary amine
"""
from rdkit import Chem

def is_primary_amine(smiles: str):
    """
    Determines if a molecule is a primary amine based on its SMILES string.
    A primary amine has exactly one hydrocarbyl group attached to an NH2/NH3+ group,
    and is not part of amide/sulfonamide groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a primary amine, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    valid_amines = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:
            # Count bonds to carbon
            carbon_bonds = sum(1 for bond in atom.GetBonds() 
                             if bond.GetOtherAtom(atom).GetAtomicNum() == 6)
            
            # Must have exactly one carbon bond
            if carbon_bonds != 1:
                continue
                
            # Check hydrogen count and charge state
            h_count = atom.GetTotalNumHs()
            charge = atom.GetFormalCharge()
            
            # Allow NH2 (neutral) or NH3+ (protonated)
            valid_charge = (h_count == 2 and charge == 0) or (h_count == 3 and charge == 1)
            if not valid_charge:
                continue
                
            # Get the attached carbon atom
            attached_carbon = next(bond.GetOtherAtom(atom) for bond in atom.GetBonds()
                                  if bond.GetOtherAtom(atom).GetAtomicNum() == 6)
            
            # Check if carbon is part of carbonyl group (amide exclusion)
            has_carbonyl = any(bond.GetBondType() == Chem.BondType.DOUBLE and 
                              bond.GetOtherAtom(attached_carbon).GetAtomicNum() == 8 
                              for bond in attached_carbon.GetBonds())
            
            # Check if carbon is part of sulfonamide/sulfonyl group
            has_sulfonyl = any(bond.GetBondType() == Chem.BondType.DOUBLE and 
                              bond.GetOtherAtom(attached_carbon).GetAtomicNum() == 16 
                              for bond in attached_carbon.GetBonds())
            
            if not has_carbonyl and not has_sulfonyl:
                valid_amines += 1

    if valid_amines > 0:
        return True, f"Found {valid_amines} primary amine group(s)"
    return False, "No primary amine group (NH2 with single carbon bond) found"