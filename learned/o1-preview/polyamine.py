"""
Classifies: CHEBI:88061 polyamine
"""
"""
Classifies: polyamine
"""
from rdkit import Chem

def is_polyamine(smiles: str):
    """
    Determines if a molecule is a polyamine based on its SMILES string.
    A polyamine is defined as any organic compound containing two or more amino groups.
    
    An amino group is a nitrogen atom (excluding those in amide, nitro, nitrile, or aromatic groups)
    that is sp3 hybridized and bonded only to carbons or hydrogens.
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if molecule is a polyamine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Initialize amino group count
    amino_group_count = 0
    
    # Iterate over all nitrogen atoms in the molecule
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() != 7:
            continue  # Skip if not nitrogen
        
        # Exclude nitrogen atoms in aromatic rings
        if atom.GetIsAromatic():
            continue
        
        # Exclude nitrogen atoms not sp3 hybridized
        if atom.GetHybridization() != Chem.rdchem.HybridizationType.SP3:
            continue
        
        # Exclude nitrogen atoms bonded to heteroatoms other than nitrogen
        bonded_to_heteroatom = False
        for neighbor in atom.GetNeighbors():
            atomic_num = neighbor.GetAtomicNum()
            if atomic_num not in [1, 6, 7]:  # Hydrogen, Carbon, Nitrogen
                bonded_to_heteroatom = True
                break
        if bonded_to_heteroatom:
            continue
        
        # Exclude nitrogen atoms in amides, nitro groups, nitriles
        is_amide = False
        is_nitro = False
        is_nitrile = False
        
        for bond in atom.GetBonds():
            neighbor = bond.GetOtherAtom(atom)
            # Check for amide N-C(=O)
            if neighbor.GetAtomicNum() == 6 and bond.GetBondType() == Chem.rdchem.BondType.SINGLE:
                # Check if the carbon is double-bonded to oxygen
                for nbr_bond in neighbor.GetBonds():
                    nbr_neighbor = nbr_bond.GetOtherAtom(neighbor)
                    if nbr_neighbor.GetAtomicNum() == 8 and nbr_bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                        is_amide = True
                        break
            # Check for nitro group N(=O)=O
            elif neighbor.GetAtomicNum() == 8 and bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                is_nitro = True
            # Check for nitrile N#C
            elif neighbor.GetAtomicNum() == 6 and bond.GetBondType() == Chem.rdchem.BondType.TRIPLE:
                is_nitrile = True
        
        if is_amide or is_nitro or is_nitrile:
            continue  # Exclude these functional groups
        
        # If passed all checks, count as amino group
        amino_group_count += 1
    
    if amino_group_count >= 2:
        return True, f"Contains {amino_group_count} amino groups"
    else:
        return False, f"Contains {amino_group_count} amino group(s), which is less than 2"