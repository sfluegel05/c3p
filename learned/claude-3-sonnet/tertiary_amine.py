"""
Classifies: CHEBI:32876 tertiary amine
"""
"""
Classifies: tertiary amine compounds
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_tertiary_amine(smiles: str):
    """
    Determines if a molecule contains a tertiary amine group.
    A tertiary amine has a nitrogen atom with three single bonds to carbon atoms.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains tertiary amine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Exclude patterns - functional groups that aren't tertiary amines
    exclude_patterns = [
        '[NX3](=[OX1])',  # N-oxide
        '[NX3](=[SX1])',  # N-thiooxide
        '[NX3][CX3](=[OX1])',  # Amide
        '[NX3][CX3](=[SX1])',  # Thioamide
        '[NX3][SX4](=[OX1])(=[OX1])',  # Sulfonamide
        '[n]',  # Aromatic nitrogen
        '[N+]',  # Quaternary nitrogen
        '[N-]',  # Negatively charged nitrogen
        '[NX3]=[NX2]',  # Diazo
        '[NX3]=[CX3]',  # Imine
        '[NX3]#[CX2]',  # Nitrile
        '[NX3]=[NX3]',  # Azoxy
        '[NX3][NX3]',  # Hydrazine
        '[NX3][OX2]',  # N-O bond
        '[NX3][SX2]',  # N-S bond
        '[NX3]([CX3](=O))[CX3](=O)',  # N with two carbonyls
        '[NX3]C=[NX2]',  # Amidine
        '[NX3]C(=[NX2])N',  # Guanidine
        '[N]1[C](=O)[C]1',  # Aziridine
        '[N]1[C](=O)[CH2][C]1',  # β-lactam
        '[N+](=[O])[O-]',  # Nitro
    ]

    # Convert patterns to mol objects
    exclude_mols = [Chem.MolFromSmarts(pattern) for pattern in exclude_patterns]

    # Find tertiary amine pattern
    tertiary_amine_pattern = Chem.MolFromSmarts('[NX3;H0;!$(N[*]=*);!$(N#*);!$(N=*);!$([N][!#6]);!$(N[*]=[!#6]);!$(N[*]#[!#6]);!$(N@[*]@[*]@N)]([#6])[#6][#6]')
    
    if not mol.HasSubstructMatch(tertiary_amine_pattern):
        return False, "No tertiary amine group found"

    # Check for excluding patterns
    for pattern in exclude_mols:
        if pattern is not None and mol.HasSubstructMatch(pattern):
            return False, "Contains excluded functional group"

    # Additional checks for each matching nitrogen
    matches = mol.GetSubstructMatches(tertiary_amine_pattern)
    for match in matches:
        n_idx = match[0]
        n_atom = mol.GetAtomWithIdx(n_idx)
        
        # Check hybridization
        if n_atom.GetHybridization() != Chem.HybridizationType.SP3:
            continue
            
        # Verify all bonds are single
        all_single = True
        carbon_count = 0
        for bond in n_atom.GetBonds():
            if bond.GetBondType() != Chem.BondType.SINGLE:
                all_single = False
                break
            # Count carbons connected by single bonds
            other_atom = bond.GetOtherAtom(n_atom)
            if other_atom.GetAtomicNum() == 6:
                # Check the carbon isn't part of a carbonyl
                if not any(b.GetBondType() != Chem.BondType.SINGLE for b in other_atom.GetBonds()):
                    carbon_count += 1
                    
        if all_single and carbon_count == 3:
            return True, "Contains a nitrogen atom with three single bonds to carbon atoms"

    return False, "No valid tertiary amine group found"