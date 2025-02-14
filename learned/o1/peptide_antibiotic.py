"""
Classifies: CHEBI:25903 peptide antibiotic
"""
from rdkit import Chem
from rdkit.Chem import rdchem

def is_peptide_antibiotic(smiles: str):
    """
    Determines if a molecule is a peptide antibiotic based on its SMILES string.
    A peptide antibiotic is a peptide that exhibits antimicrobial properties.
    
    This function checks for features common in peptide antibiotics:
    - Presence of peptide bonds (amide bonds)
    - Presence of cyclic peptides (rings involving peptide bonds)
    - Contains unusual amino acids or modifications (e.g., D-amino acids, thiazole rings)
    
    Args:
        smiles (str): SMILES string of the molecule
    
    Returns:
        bool: True if the molecule is likely a peptide antibiotic, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Detect peptide bonds (amide bonds between C=O and N)
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    peptide_bonds = mol.GetSubstructMatches(peptide_bond_pattern)
    num_peptide_bonds = len(peptide_bonds)
    
    # Detect cyclic peptides (rings involving peptide bonds)
    is_cyclic_peptide = False
    sssr = mol.GetRingInfo().AtomRings()
    for ring in sssr:
        ring_atoms = [mol.GetAtomWithIdx(idx) for idx in ring]
        # Check if ring contains a peptide bond
        contains_peptide_bond = False
        for bond in mol.GetBonds():
            if bond.GetBeginAtomIdx() in ring and bond.GetEndAtomIdx() in ring:
                if bond.GetBondType() == rdchem.BondType.SINGLE:
                    begin_atom = bond.GetBeginAtom()
                    end_atom = bond.GetEndAtom()
                    if ((begin_atom.GetSymbol() == 'C' and begin_atom.GetDoubleProp('_TriposPartialCharge', 0) > 0 and
                         end_atom.GetSymbol() == 'N') or
                        (begin_atom.GetSymbol() == 'N' and
                         end_atom.GetSymbol() == 'C' and end_atom.GetDoubleProp('_TriposPartialCharge', 0) > 0)):
                        contains_peptide_bond = True
                        break
        if contains_peptide_bond:
            is_cyclic_peptide = True
            break

    # Detect unusual amino acids or modifications
    # Examples include thiazole, oxazole rings commonly found in peptide antibiotics
    unusual_motifs = [
        Chem.MolFromSmarts("c1ncsc1"),  # Thiazole ring
        Chem.MolFromSmarts("c1nco1"),   # Oxazole ring
        Chem.MolFromSmarts("C=CSC=N"),  # Dehydroalanine
        Chem.MolFromSmarts("C1=CSC=N1") # Thiazoline ring
    ]
    contains_unusual_motifs = False
    for motif in unusual_motifs:
        if mol.HasSubstructMatch(motif):
            contains_unusual_motifs = True
            break
    
    # Apply classification rules
    if num_peptide_bonds == 0:
        return False, f"No peptide bonds detected"
    elif is_cyclic_peptide and contains_unusual_motifs:
        return True, "Cyclic peptide with unusual amino acids; likely a peptide antibiotic"
    elif is_cyclic_peptide:
        return True, "Cyclic peptide; likely a peptide antibiotic"
    elif contains_unusual_motifs:
        return True, "Linear peptide with unusual amino acids; may be a peptide antibiotic"
    elif num_peptide_bonds >= 5:
        return True, f"Contains {num_peptide_bonds} peptide bonds; possible peptide antibiotic"
    else:
        return False, f"Contains {num_peptide_bonds} peptide bonds; insufficient evidence for peptide antibiotic"